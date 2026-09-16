"""Parent CLI publication integration, not an ML/MM calculation integration.

MM workspace preparation and scan/path-opt calculations are stubbed. Input
preparation, MEP concatenation/conversion, public copies and manifest stay real.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest
from ase import Atoms
from click.testing import CliRunner

from mlmm.core import utils
from mlmm.core.result_commit import MLMM_RUN_ID_ENV, with_current_run_id
from mlmm.workflows import all as all_workflow
from mlmm.workflows import dft


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.next")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)
    return path


def _coords(index: int):
    return [(index * 0.125, 0.25, -0.5), (1.25 + index * 0.25, 0.5, 0.75)]


def _frame(index: int) -> str:
    atoms = "".join(
        f"{element} {x:.3f} {y:.3f} {z:.3f}\n"
        for element, (x, y, z) in zip(("C", "O"), _coords(index))
    )
    return f"2\nE={index * 0.01:.3f} unit=hartree frame{index}\n{atoms}"


def _pdb(residue: str) -> str:
    return "".join(
        f"HETATM{serial:5d} {name:<4s} {residue:3s} A   1    "
        f"{x:8.3f}{8.0:8.3f}{9.0:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
        for serial, name, element, x in ((1, "C1", "C", 8.0), (2, "O1", "O", 9.0))
    ) + "END\n"


def _assert_pdb(path: Path, frames, residue: str) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    assert sum(line.startswith("MODEL") for line in lines) == len(frames)
    assert sum(line.startswith("ENDMDL") for line in lines) == len(frames)
    atoms = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    assert len(atoms) == 2 * len(frames)
    assert [line[12:16].strip() for line in atoms] == ["C1", "O1"] * len(frames)
    assert [line[76:78].strip() for line in atoms] == ["C", "O"] * len(frames)
    assert {line[17:20] for line in atoms} == {residue}
    actual = [tuple(float(line[start:start + 8]) for start in (30, 38, 46)) for line in atoms]
    assert actual == [xyz for index in frames for xyz in _coords(index)]


@pytest.mark.parametrize("plot_writes", [False, True])
@pytest.mark.parametrize("case", ["scan_preopt", "direct_pdb", "no_convert"])
def test_all_path_opt_mep_pdb_publication(tmp_path: Path, monkeypatch, case: str, plot_writes: bool) -> None:
    monkeypatch.delenv(MLMM_RUN_ID_ENV, raising=False)
    # Isolate inherited process state. Keep conversion enabled even in the
    # negative case, so the parent's explicit --no-convert-files guard is tested.
    monkeypatch.setattr(utils, "_CONVERT_FILES_ENABLED", True)
    out = tmp_path / "out"
    stale_diagram = _write(out / "_work/path_opt/energy_diagram_MEP.png", "old image")
    _write(out / "_work/path_opt/mep_plot.png", "old plot")
    convert = case != "no_convert"
    inputs = [_write(tmp_path / f"input{i}.pdb", _pdb(label)) for i, label in enumerate(("RAW", "MID", "END"))]
    # Explicit model/parm avoid AmberTools; this fake parm is never parsed.
    parm = _write(tmp_path / "unused.parm7", "publication fixture; not a force field\n")
    model = Atoms("CO", positions=_coords(0))
    workspace = SimpleNamespace(
        atoms_model=model, atoms_model_lh=model.copy(), model_pdb=inputs[0],
        link_pairs=[], cleanup=lambda: None,
    )
    monkeypatch.setattr(dft, "_prepare_ml_region_workspace", lambda **_k: workspace)
    references = []
    child_calls = []

    def fake_child(name, _cli, args, **_kwargs):
        child_calls.append(name)
        child_out = Path(args[args.index("--out-dir") + 1])
        if name == "scan":
            _write(child_out / "preopt/result.xyz", _frame(0))
            _write(child_out / "preopt/result.pdb", _pdb("PRE"))
            for index in (1, 2):
                _write(child_out / f"stage_{index:02d}/result.xyz", _frame(index))
                _write(child_out / f"stage_{index:02d}/result.pdb", _pdb(f"S{index:02d}"))
            payload = {"scientific_status": "success", "preopt_converged": True, "stages": []}
        elif name == "path_opt":
            pair = len(references)
            refs = [Path(args[i + 1]) for i, token in enumerate(args) if token == "--ref-pdb"]
            references.append(refs)
            left = Path(args[args.index("-i") + 1])
            assert left.suffix == (".pdb" if case == "direct_pdb" else ".xyz")
            _write(child_out / "final_geometries_trj.xyz", _frame(pair) + _frame(pair + 1))
            payload = {"stage_outcomes": [{"item_id": "gsm_mep", "converged": True}], "preopt_converged": True}
        else:
            raise AssertionError(f"Unexpected computational child: {name}")
        _write(child_out / "result.json", json.dumps(with_current_run_id(payload)))

    def no_calculator(*_args, **_kwargs):
        raise AssertionError("This publication test must not construct a calculator")

    monkeypatch.setattr(all_workflow, "_run_cli_main", fake_child)
    monkeypatch.setattr(all_workflow, "_mlmm_calc", no_calculator)
    def plot(_source, outputs, **kwargs):
        if not plot_writes:
            raise RuntimeError("plot renderer unavailable")
        Path(outputs[0]).write_text("current plot")
    monkeypatch.setattr(all_workflow, "run_trj2fig", plot)
    monkeypatch.setattr(all_workflow, "close_matplotlib_figures", lambda: None)
    monkeypatch.setattr(all_workflow, "_write_segment_energy_diagram", lambda *_a, **_k: None)
    monkeypatch.setattr(all_workflow._path_search, "_has_bond_change", lambda *_a, **_k: (False, ""))
    selected = inputs if case == "direct_pdb" else inputs[:1]
    args = [arg for path in selected for arg in ("-i", str(path))]
    args += ["-q", "0", "-m", "1", "--parm", str(parm), "--model-pdb", str(inputs[0]),
             "--out-dir", str(out), "--convert-files" if convert else "--no-convert-files"]
    args += (["--no-preopt"] if case == "direct_pdb" else
             ["--preopt", "--scan-lists", "[(1,2,1.5)]", "--scan-lists", "[(1,2,1.75)]"])
    result = CliRunner().invoke(all_workflow.cli, args)
    assert result.exit_code == 0, f"{result.output}\n{result.exception!r}"
    assert child_calls == (["scan"] if case != "direct_pdb" else []) + ["path_opt", "path_opt"]
    summary = json.loads((out / "summary.json").read_text())
    manifest = json.loads((out / "_work/_run_manifest.json").read_text())
    assert summary["run_id"] == manifest["run_id"]
    assert stale_diagram.read_text() == "old image"
    assert not (out / "energy_diagram_MEP.png").exists()
    assert (out / "mep_plot.png").exists() is plot_writes
    if plot_writes:
        assert (out / "mep_plot.png").read_text() == "current plot"
    assert ("mep_plot.png" in summary["key_output_files"]) is plot_writes
    assert "energy_diagram_MEP.png" not in summary["key_output_files"]
    assert summary["n_images"] == 3 and summary["n_segments"] == 2
    assert (out / "mep_trj.xyz").read_text() == "".join(_frame(i) for i in range(3))
    assert not (out / "mep_trj.pdb").exists()
    assert (out / "mep.pdb").exists() is convert
    assert ("mep.pdb" in summary["key_output_files"]) is convert
    assert ("output.public.mep.pdb" in manifest["produced"]) is convert
    path_dir = out / "_work/path_opt"
    labels = ["RAW", "MID"] if case == "direct_pdb" else ["RAW", "S01"]
    if case == "direct_pdb":
        assert references == [[], []]  # Native PDB inputs supply their own templates.
    else:
        # ML/MM intentionally keeps the original layered template for preopt.
        assert [len(refs) for refs in references] == [2, 2]
        assert [[next(line[17:20] for line in ref.read_text().splitlines() if line.startswith("HETATM"))
                 for ref in refs] for refs in references] == [["RAW", "S01"], ["S01", "S02"]]
    for index in (1, 2):
        segment_pdb = path_dir / f"mep_seg_{index:02d}.pdb"
        assert segment_pdb.exists() is convert
        if convert:
            _assert_pdb(segment_pdb, [index - 1, index], labels[index - 1])
    if convert:
        _assert_pdb(out / "mep.pdb", [0, 1, 2], labels[0])
