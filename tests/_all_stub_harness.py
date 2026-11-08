"""Run the ``mlmm all`` CLI end to end on a two-atom fixture with stubbed children.

Every child command is replaced at ``all._run_cli_main`` (and the DFT child at
its subprocess call) by a stub that writes the files the real child writes, and
the parent ML/MM calculator is an energy-only stub. Input preparation, layer
handling, the reading of child results, summary assembly and publication run
unchanged; no force field, MLIP or DFT runs. Energies are fixed functions of the
C-atom position, so every exit gives the same numbers on every run.

Used by the root ``summary.json`` golden tests and by the behaviour tests that
replaced the source-text checks of the ``all`` callback.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence

import numpy as np
import yaml
from ase import Atoms
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.core import utils
from mlmm.core.result_commit import MLMM_RUN_ID_ENV, with_current_run_id
from mlmm.workflows import all as all_workflow
from mlmm.workflows import dft

BOHR2ANG = 0.529177210903


def coords(frame: float):
    """C and O positions (Å) of reaction-coordinate frame 0 (R) .. 2 (P)."""
    return [(frame * 0.125, 0.25, -0.5), (1.25 + frame * 0.25, 0.5, 0.75)]


def frame_of(x_carbon_ang: float) -> float:
    return float(x_carbon_ang) / 0.125


def mlmm_energy(frame: float) -> float:
    """ML/MM energy (Eh): R = -1.0, TS = -0.975, P = -0.99."""
    return -1.0 + 0.02 * frame * (2.0 - frame) + 0.005 * frame


def dft_energy(frame: float) -> float:
    return -100.0 + 0.03 * frame * (2.0 - frame) + 0.004 * frame


def xyz_text(frame: float, energy: Optional[float] = None, *, values=None) -> str:
    rows = values if values is not None else coords(frame)
    energy = mlmm_energy(frame) if energy is None else energy
    atoms = "".join(
        f"{element} {x:.6f} {y:.6f} {z:.6f}\n"
        for element, (x, y, z) in zip(("C", "O"), rows)
    )
    return f"2\nE={energy:.12f} unit=hartree frame{frame:g}\n{atoms}"


def pdb_text(frame: float, *, values=None, bfactors=(0.0, 0.0)) -> str:
    rows = values if values is not None else coords(frame)
    return "".join(
        f"HETATM{serial:5d} {name:<4s} MOL A   1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{b:6.2f}          {element:>2s}\n"
        for serial, (name, element, (x, y, z), b) in enumerate(
            zip(("C1", "O1"), ("C", "O"), rows, bfactors), 1
        )
    ) + "END\n"


def write(path: Path, text: str) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.next")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)
    return path


def write_json(path: Path, payload: Mapping[str, Any]) -> Path:
    return write(path, json.dumps(with_current_run_id(dict(payload))))


def read_frame(path: Path) -> float:
    """Frame of the first structure in an XYZ or PDB file (from the C atom)."""
    text = Path(path).read_text(encoding="utf-8").splitlines()
    if Path(path).suffix.lower() == ".pdb":
        line = next(row for row in text if row.startswith(("ATOM", "HETATM")))
        return frame_of(float(line[30:38]))
    return frame_of(float(text[2].split()[1]))


def read_values(path: Path):
    """C and O positions of the first structure in an XYZ file."""
    rows = Path(path).read_text(encoding="utf-8").splitlines()[2:4]
    return [tuple(float(v) for v in row.split()[1:4]) for row in rows]


def option(args: Sequence[str], flag: str) -> Optional[str]:
    args = list(args)
    return args[args.index(flag) + 1] if flag in args else None


TS_CONVERGED = {
    "optimization_status": "converged", "converged": True,
    "execution_status": "completed", "scientific_status": "success",
    "hessian_status": "completed", "saddle_validation": "first_order",
    "n_imaginary_modes": 1, "n_negative_modes": 1,
    "imaginary_frequencies_cm": [-512.5], "frequency_zero_cutoff_cm": 10.0,
    "imaginary_mode_criterion": "frequency", "imaginary_frequency_threshold_cm": -20.0,
    "reaction_mode_index": 0, "reaction_mode_frequency_cm": -512.5,
    "n_opt_cycles": 7, "max_cycles": 100, "stop_reason": None,
}

TS_NOT_CONVERGED = {
    "optimization_status": "not_converged", "converged": False,
    "execution_status": "completed", "scientific_status": "failed",
    "hessian_status": "skipped", "n_imaginary_modes": None,
    "n_opt_cycles": 100, "max_cycles": 100, "stop_reason": "max_cycles",
}


@dataclass
class Behaviour:
    """What each stubbed child reports.

    ``opt`` maps the endpoint role (``R``/``P``, the child out-dir name) to
    ``converged``, ``not_converged``, ``error`` (the child raises) or
    ``nonfinite`` (the child writes NaN coordinates).
    """

    tsopt: Dict[str, Any] = field(default_factory=lambda: dict(TS_CONVERGED))
    opt: Dict[str, str] = field(default_factory=dict)
    path_opt_converged: Sequence[bool] = ()
    path_search_segments: Optional[List[Dict[str, Any]]] = None
    opt_cycles: Mapping[str, int] = field(default_factory=lambda: {"R": 3, "P": 4})


class EnergyCalculator:
    """Parent-side calculator stub: an energy from the C-atom position."""

    def get_energy(self, _atoms, coords3):
        values = np.asarray(coords3, dtype=float).reshape(-1)
        return {"energy": mlmm_energy(frame_of(values[0] * BOHR2ANG))}

    def save_scf_checkpoint(self, destination, _geometry):
        destination = Path(destination)
        write(destination, "checkpoint fixture\n")
        write(destination.with_suffix(destination.suffix + ".json"), "{}\n")

    def close(self) -> None:
        pass


@dataclass
class AllRun:
    out: Path
    result: Any
    calls: List[tuple]

    @property
    def summary(self) -> Dict[str, Any]:
        return json.loads((self.out / "summary.json").read_text(encoding="utf-8"))

    def child_names(self) -> List[str]:
        return [name for name, _args in self.calls]

    def child_args(self, name: str) -> List[List[str]]:
        return [args for child, args in self.calls if child == name]


def _thermo_payload(frame: float, role: str) -> Dict[str, Any]:
    energy = mlmm_energy(frame)
    correction = 0.010 + 0.002 * frame
    return {
        "temperature_K": 298.15, "pressure_atm": 1.0,
        "point_group": "Cinfv", "point_group_source": "detected",
        "symmetry_number": 1, "symmetry_number_source": "point_group",
        "num_imag_freq": 1 if role == "TS" else 0,
        "n_negative_modes": 1 if role == "TS" else 0,
        "electronic_energy_ha": energy,
        "zpe_correction_ha": 0.005,
        "thermal_correction_energy_ha": 0.007,
        "thermal_correction_enthalpy_ha": 0.008,
        "thermal_correction_free_energy_ha": correction,
        "sum_EE_and_thermal_free_energy_ha": energy + correction,
        "imaginary_frequencies_cm": [-512.5] if role == "TS" else [],
    }


def run_all(
    tmp_path: Path,
    monkeypatch,
    *,
    frames: Sequence[float] = (0.0, 2.0),
    extra: Sequence[str] = (),
    behaviour: Optional[Behaviour] = None,
    xyz_inputs: bool = False,
    wrap: Optional[Callable[[], None]] = None,
    parm: bool = True,
    model_pdb: bool = True,
    bfactors: Sequence[float] = (0.0, 0.0),
) -> AllRun:
    """Run ``mlmm all`` and return the output dir, CLI result and child calls.

    With ``parm=False`` no ``--parm`` is passed, so ``all`` builds the parm7 with
    ``_build_mm_parm7`` (the caller replaces it). With ``model_pdb=False`` no
    ``--model-pdb`` is passed. ``bfactors`` are the B-factors (the layer encoding)
    of the C and O atoms of the input PDBs.
    """

    behaviour = behaviour or Behaviour()
    monkeypatch.delenv(MLMM_RUN_ID_ENV, raising=False)
    monkeypatch.setattr(utils, "_CONVERT_FILES_ENABLED", True)
    pdbs = [
        write(tmp_path / f"input{i}.pdb", pdb_text(f, bfactors=bfactors))
        for i, f in enumerate(frames)
    ]
    if xyz_inputs:
        inputs = [write(tmp_path / f"input{i}.xyz", xyz_text(f)) for i, f in enumerate(frames)]
    else:
        inputs = pdbs
    # Only the ML/MM workspace construction is replaced; this parm7 fixture is
    # not a force field and never reaches a calculator or topology reader.
    parm_file = write(tmp_path / "fixture.parm7", "dispatch fixture; not a force field\n")
    model = Atoms("CO", positions=coords(frames[0]))
    workspace = SimpleNamespace(
        atoms_model=model, atoms_model_lh=model.copy(), model_pdb=pdbs[0],
        link_pairs=[], cleanup=lambda: None,
    )
    monkeypatch.setattr(dft, "_prepare_ml_region_workspace", lambda **kwargs: workspace)
    calls: List[tuple] = []

    def child(name, _cli, args, **_kwargs):
        args = [str(arg) for arg in args]
        calls.append((name, args))
        out = Path(option(args, "--out-dir"))
        source = Path(option(args, "-i")) if option(args, "-i") else None
        if name == "path_opt":
            pair = sum(1 for call in calls if call[0] == "path_opt") - 1
            left = read_frame(source)
            right = read_frame(Path(args[args.index("-i") + 2]))
            middle = (left + right) / 2.0
            write(out / "final_geometries_trj.xyz",
                  "".join(xyz_text(f) for f in (left, middle, right)))
            write(out / "hei.xyz", xyz_text(middle))
            write(out / "hei.pdb", pdb_text(middle))
            converged = (
                behaviour.path_opt_converged[pair]
                if pair < len(behaviour.path_opt_converged) else True
            )
            write_json(out / "result.json", {
                "optimization_status": "converged" if converged else "not_converged",
                "converged": converged,
                "preopt_requested": False, "preopt_converged": None,
                "path_optimizers": [],
                "stage_outcomes": [{"stage": "path-opt", "item_id": "gsm_mep",
                                    "executed": True, "converged": converged,
                                    "usable": True, "required": True}],
            })
        elif name == "path_search":
            _path_search_child(out, behaviour, args)
        elif name == "tsopt":
            ts_frame = read_frame(source)
            write(out / "final_geometry.xyz", xyz_text(ts_frame))
            write(out / "final_geometry.pdb", pdb_text(ts_frame))
            payload = dict(behaviour.tsopt)
            payload.setdefault("energy_hartree", mlmm_energy(ts_frame))
            payload["files"] = {"final_geometry": str(out / "final_geometry.xyz")}
            write_json(out / "result.json", payload)
        elif name == "irc":
            irc_frames = (0.0, 1.0, 2.0)
            trajectory = write(
                out / "finished_irc_trj.xyz", "".join(xyz_text(f) for f in irc_frames),
            )
            write(out / "forward_irc_trj.xyz", "".join(xyz_text(f) for f in irc_frames[1:]))
            write(out / "backward_irc_trj.xyz", "".join(xyz_text(f) for f in irc_frames[:2]))
            write_json(out / "result.json", {
                "files": {"finished_irc": str(trajectory)},
                "forward_requested": True, "backward_requested": True,
                "forward_integration_converged": True,
                "backward_integration_converged": True,
                "forward_integration_stop_reason": "energy_minimum",
                "backward_integration_stop_reason": "energy_minimum",
                "n_frames_forward": 2, "n_frames_backward": 2,
            })
        elif name == "opt":
            role = out.name
            mode = behaviour.opt.get(role, "converged")
            values = read_values(source) if source.suffix == ".xyz" else coords(read_frame(source))
            write(out / "iterate.xyz", xyz_text(0, -1.0, values=[(x + 10.0, y, z) for x, y, z in values]))
            if mode == "error":
                raise RuntimeError(f"endpoint numerical failure ({role})")
            if mode == "nonfinite":
                values = [(float("nan"), y, z) for x, y, z in values]
            frame = frame_of(values[0][0]) if np.isfinite(values[0][0]) else 0.0
            text = xyz_text(frame, values=values)
            write(out / "final_geometry.xyz", text)
            write(out / "final_geometry.pdb", pdb_text(frame, values=values))
            status = "not_converged" if mode == "not_converged" else "converged"
            write_json(out / "result.json", {
                "optimization_status": status,
                "converged": status == "converged",
                "n_opt_cycles": behaviour.opt_cycles.get(role, 3),
                "max_cycles": 3000,
            })
        elif name == "freq":
            role = out.name
            frame = read_frame(source)
            payload = _thermo_payload(frame, role)
            if "--dump" in args:
                write(out / "thermoanalysis.yaml", yaml.safe_dump(payload, sort_keys=False))
            write_json(out / "result.json", {"num_imag_freq": payload["num_imag_freq"]})
        else:
            raise AssertionError(f"Unexpected child command: {name}")
        return 0

    real_run = subprocess.run

    def dft_subprocess(cmd, *args, **kwargs):
        cmd = [str(part) for part in cmd]
        if cmd[1:4] != ["-m", "mlmm", "dft"]:
            return real_run(cmd, *args, **kwargs)
        calls.append(("dft", cmd[4:]))
        out = Path(option(cmd, "--out-dir"))
        frame = read_frame(Path(option(cmd, "-i")))
        energy = dft_energy(frame)
        write(out / "result.yaml", yaml.safe_dump({
            "energy": {"hartree": energy, "converged": True},
            "mlmm_energy": {"E_total_ml_dft_mm_hartree": energy - 0.5},
        }, sort_keys=False))
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    def plot(_source, outputs, **_kwargs):
        write(Path(outputs[0]), "plot fixture\n")

    def image(_figure, destination, **_kwargs):
        write(Path(destination), "image fixture\n")
        return Path(destination)

    monkeypatch.setattr(all_workflow, "_run_cli_main", child)
    monkeypatch.setattr(subprocess, "run", dft_subprocess)
    monkeypatch.setattr(all_workflow, "_mlmm_calc", lambda **_cfg: EnergyCalculator())
    monkeypatch.setattr(all_workflow, "run_trj2fig", plot)
    monkeypatch.setattr(all_workflow, "close_matplotlib_figures", lambda: None)
    monkeypatch.setattr(all_workflow, "write_plotly_image", image)
    if wrap is not None:
        wrap()

    out = tmp_path / "out"
    args = ["all", *[a for path in inputs for a in ("-i", str(path))]]
    if xyz_inputs:
        args += ["--ref-pdb", str(pdbs[0])]
    args += ["-q", "0", "-m", "1", *(["--parm", str(parm_file)] if parm else []),
             *(["--model-pdb", str(pdbs[0])] if model_pdb else []),
             "--out-dir", str(out), "--no-preopt", *extra]
    result = CliRunner().invoke(root_cli, args)
    return AllRun(out=out, result=result, calls=calls)


def _path_search_child(out: Path, behaviour: Behaviour, args: Sequence[str]) -> None:
    """Write the files that ``path-search`` leaves for ``all``.

    The child summary is completed by the real path-search contract helper, so
    the root summary that ``all`` republishes from it has the real child keys.
    """

    from mlmm.core.defaults import MLMM_CALC_KW
    from mlmm.io.summary import method_references
    from mlmm.workflows import path_search as path_search_workflow

    segments = behaviour.path_search_segments
    if segments is None:
        segments = [{
            "index": 1, "tag": "seg_001", "kind": "seg",
            "barrier_kcal": (mlmm_energy(1.0) - mlmm_energy(0.0)) * 627.509474,
            "delta_kcal": (mlmm_energy(2.0) - mlmm_energy(0.0)) * 627.509474,
            "bond_changes": "C1-O1 : formed", "converged": True,
        }]
    frames = (0.0, 1.0, 2.0)

    def models(selected) -> str:
        return "".join(
            f"MODEL     {m + 1:4d}\n" + pdb_text(f).replace("END\n", "ENDMDL\n")
            for m, f in enumerate(selected)
        ) + "END\n"

    write(out / "mep_trj.xyz", "".join(xyz_text(f) for f in frames))
    write(out / "mep_trj.pdb", models(frames))
    write(out / "mep_plot.png", "plot fixture\n")
    write(out / "energy_diagram_MEP.png", "image fixture\n")
    for segment in segments:
        if segment.get("kind") == "bridge":
            continue
        index = int(segment["index"])
        write(out / f"mep_seg_{index:02d}_trj.xyz", "".join(xyz_text(f) for f in frames))
        write(out / f"mep_seg_{index:02d}.pdb", models((0.0, 2.0)))
        write(out / f"hei_seg_{index:02d}.xyz", xyz_text(1.0))
        write(out / f"hei_seg_{index:02d}.pdb", pdb_text(1.0))
    energies = [mlmm_energy(f) for f in frames]
    summary: Dict[str, Any] = {
        "out_dir": str(out),
        "n_images": len(frames),
        "n_segments": len(segments),
        "segments": [dict(segment) for segment in segments],
    }
    if any(segment.get("kind", "seg") == "seg" for segment in segments):
        summary["energy_diagrams"] = [{
            "name": "energy_diagram_MEP", "labels": ["R", "TS1", "P"],
            "energies_kcal": [(e - energies[0]) * 627.509474 for e in energies],
            "ylabel": "ΔE (kcal/mol)", "energies_au": energies,
            "image": str(out / "energy_diagram_MEP.png"),
        }]
    reports = [
        path_search_workflow.SegmentReport(
            tag=segment["tag"], barrier_kcal=segment["barrier_kcal"],
            delta_kcal=segment["delta_kcal"], summary=segment.get("bond_changes", ""),
            kind=segment.get("kind", "seg"), seg_index=int(segment["index"]),
            converged=segment.get("converged"),
        )
        for segment in segments
    ]
    path_search_workflow._enrich_path_summary_contract(
        summary, segments=reports, out_dir=out,
        calc_cfg={**MLMM_CALC_KW, "model_charge": int(option(args, "-q")),
                  "model_mult": int(option(args, "-m"))},
        command="mlmm path-search " + " ".join(args),
    )
    summary["search_max_depth"] = 10
    summary["path_optimizers"] = []
    summary["references"] = method_references({
        "pipeline_mode": "path-search", "mep_mode": "gsm", "path_opt_mode": "grad",
        "preopt": False, "path_optimizers": [],
        "mlip_backend": summary.get("mlip_backend"),
        "mlip_model": summary.get("mlip_model"), "mlip_task": summary.get("mlip_task"),
    })
    write_json(out / "summary.json", summary)
    write(out / "summary.log", "path-search summary fixture\n")


__all__ = [
    "AllRun", "Behaviour", "TS_CONVERGED", "TS_NOT_CONVERGED", "coords",
    "dft_energy", "mlmm_energy", "option", "pdb_text", "read_frame", "read_values", "run_all",
    "write", "write_json", "xyz_text",
]
