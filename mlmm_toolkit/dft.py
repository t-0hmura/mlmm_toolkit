# mlmm_toolkit/dft.py

"""ML/MM-aware single-point DFT for the ML region with energy recombination.

Example:
    mlmm dft -i enzyme.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0

For detailed documentation, see: docs/dft.md
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import shutil
import sys
import tempfile
import textwrap
import time
import traceback

import click
import numpy as np
import yaml

from ase import Atoms
from ase.io import read

from pysisyphus.constants import AU2EV

from .mlmm_calc import hessianffCalculator
from .opt import (
    GEOM_KW as OPT_GEOM_KW,
    CALC_KW as OPT_CALC_KW,
    _parse_freeze_atoms as _parse_freeze_atoms_opt,
    _normalize_geom_freeze as _normalize_geom_freeze_opt,
)
from .utils import (
    apply_layer_freeze_constraints,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    format_freeze_atoms_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    prepare_input_structure,
    resolve_charge_spin_or_raise,
    parse_indices_string,
    build_model_pdb_from_bfactors,
    build_model_pdb_from_indices,
)

from functools import reduce

HARTREE_TO_KCALMOL = 627.5094740631
EV2AU = 1.0 / AU2EV

# -----------------------------------------------
# Defaults (override via CLI / YAML)
# -----------------------------------------------
DFT_KW: Dict[str, Any] = {
    "conv_tol": 1e-9,
    "max_cycle": 100,
    "grid_level": 3,
    "verbose": 4,
    "out_dir": "./result_dft/",
}


# -----------------------------------------------
# Helper classes & utilities
# -----------------------------------------------
@dataclass
class MLRegionWorkspace:
    tmpdir: tempfile.TemporaryDirectory
    input_pdb: Path
    real_parm7: Path
    real_rst7: Path
    model_pdb: Path
    model_parm7: Path
    model_rst7: Path
    selection_indices: List[int]
    link_pairs: List[Tuple[int, int]]  # 1-based REAL indices (ml_idx, mm_idx)
    atoms_real: Atoms
    atoms_model: Atoms
    atoms_model_lh: Atoms

    def cleanup(self) -> None:
        self.tmpdir.cleanup()


def _parse_func_basis(s: str) -> Tuple[str, str]:
    if not s or "/" not in s:
        raise click.BadParameter("Expected 'FUNC/BASIS' (e.g., 'wb97m-v/def2-tzvpd').")
    func, basis = s.split("/", 1)
    func = func.strip()
    basis = basis.strip()
    if not func or not basis:
        raise click.BadParameter("Functional or basis is empty. Example: --func-basis 'wb97m-v/6-31g**'")
    return func, basis


def _atoms_to_xyz_string(atoms: Atoms, comment: str) -> str:
    lines = [str(len(atoms)), comment]
    for sym, (x, y, z) in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        lines.append(f"{sym:<2s} {x:15.8f} {y:15.8f} {z:15.8f}")
    return "\n".join(lines) + "\n"


def _atoms_to_pyscf_atoms(atoms: Atoms) -> List[Tuple[str, Tuple[float, float, float]]]:
    entries: List[Tuple[str, Tuple[float, float, float]]] = []
    for sym, coord in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        entries.append((sym, (float(coord[0]), float(coord[1]), float(coord[2]))))
    return entries


def _load_model_region_ids(model_pdb: Path) -> set[str]:
    ids: set[str] = set()
    with model_pdb.open() as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                ids.add(f"{line[12:16].strip()} {line[17:20].strip()} {line[22:26].strip()}")
    if not ids:
        raise ValueError("No atoms found in model_pdb to define the ML region.")
    return ids


def _load_input_atoms(input_pdb: Path) -> List[Dict[str, Any]]:
    atoms: List[Dict[str, Any]] = []
    with input_pdb.open() as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            elem = line[76:78].strip()
            if not elem:
                elem = line[12:16].strip()[0]
            atoms.append(
                {
                    "idx": int(line[6:11]),
                    "id": f"{line[12:16].strip()} {line[17:20].strip()} {line[22:26].strip()}",
                    "elem": elem,
                    "coord": np.array(
                        [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                        dtype=float,
                    ),
                }
            )
    if not atoms:
        raise ValueError("No ATOM/HETATM records found in the input PDB.")
    return atoms


def _detect_link_pairs(
    leap_atoms: Sequence[Dict[str, Any]],
    ml_region_ids: set[str],
    manual_links: Optional[Sequence[Sequence[str]]],
) -> List[Tuple[int, int]]:
    if manual_links:
        processed = [(" ".join(q.split()[:3]), " ".join(m.split()[:3])) for q, m in manual_links]
        ml_indices: List[int] = []
        mm_indices: List[int] = []
        for atom in leap_atoms:
            for qnm, mnm in processed:
                if atom["id"] == qnm:
                    ml_indices.append(atom["idx"])
                elif atom["id"] == mnm:
                    mm_indices.append(atom["idx"])
        if len(set(ml_indices)) != len(ml_indices) or len(set(mm_indices)) != len(mm_indices):
            raise ValueError("Duplicated ML or MM indices detected in link_mlmm specification.")
        return list(zip(ml_indices, mm_indices))

    threshold = 1.7
    ml_set = {atom["idx"] for atom in leap_atoms if atom["id"] in ml_region_ids}
    coords = {atom["idx"]: atom["coord"] for atom in leap_atoms}
    elems = {atom["idx"]: atom["elem"] for atom in leap_atoms}
    ml_indices: List[int] = []
    mm_indices: List[int] = []
    for qidx in ml_set:
        for atom in leap_atoms:
            midx = atom["idx"]
            if midx in ml_set:
                continue
            if np.linalg.norm(coords[midx] - coords[qidx]) < threshold and (
                (elems[midx] == "C" and elems[qidx] == "C")
                or (elems[midx] == "N" and elems[qidx] == "C")
                or (elems[midx] == "C" and elems[qidx] == "N")
            ):
                ml_indices.append(qidx)
                mm_indices.append(midx)
    if len(set(ml_indices)) != len(ml_indices) or len(set(mm_indices)) != len(mm_indices):
        raise ValueError("Automatic link detection produced duplicate ML/MM indices; specify link_mlmm explicitly.")
    return list(zip(ml_indices, mm_indices))


def _append_link_hydrogens(atoms_model: Atoms, atoms_real: Atoms, link_pairs: Sequence[Tuple[int, int]]) -> Atoms:
    atoms_with_link = atoms_model.copy()
    for ml_idx1, mm_idx1 in link_pairs:
        ml_idx = ml_idx1 - 1
        mm_idx = mm_idx1 - 1
        ml_elem = atoms_real[ml_idx].symbol.strip().upper()
        if ml_elem == "C":
            dist = 1.09
        elif ml_elem == "N":
            dist = 1.01
        else:
            raise ValueError(
                f"Unsupported link-atom parent element '{ml_elem}' (only C or N are allowed)."
            )
        vec = atoms_real[mm_idx].position - atoms_real[ml_idx].position
        R = np.linalg.norm(vec)
        if R < 1e-8:
            continue
        pos = atoms_real[ml_idx].position + (vec / R) * dist
        atoms_with_link += Atoms("H", positions=[pos])
    return atoms_with_link


def _prepare_ml_region_workspace(
    *,
    input_pdb: Path,
    real_parm7: Path,
    model_pdb: Path,
    link_mlmm: Optional[Sequence[Sequence[str]]],
) -> MLRegionWorkspace:
    tmpdir = tempfile.TemporaryDirectory()
    tmp = Path(tmpdir.name)

    input_copy = tmp / "input.pdb"
    real_copy = tmp / "real.parm7"
    model_copy = tmp / "model.pdb"
    shutil.copyfile(input_pdb, input_copy)
    shutil.copyfile(real_parm7, real_copy)
    shutil.copyfile(model_pdb, model_copy)

    real_top = None
    try:
        import parmed as pmd

        real_top = pmd.load_file(str(real_copy))
        start_struct = pmd.load_file(str(input_copy))
        real_top.coordinates = start_struct.coordinates
        real_top.box = None
        real_top.save(str(real_copy), overwrite=True)
        real_rst7 = tmp / "real.rst7"
        real_top.save(str(real_rst7), overwrite=True)
    except Exception as exc:  # pragma: no cover - requires parmed
        tmp.cleanup()
        raise RuntimeError(f"Failed to prepare sanitized Amber inputs: {exc}") from exc

    ml_region_ids = _load_model_region_ids(model_copy)
    leap_atoms = _load_input_atoms(input_copy)
    ml_ids = [atom["idx"] for atom in leap_atoms if atom["id"] in ml_region_ids]
    if not ml_ids:
        tmp.cleanup()
        raise ValueError("No overlap between model_pdb atoms and the input PDB was found.")

    link_pairs = _detect_link_pairs(leap_atoms, ml_region_ids, link_mlmm)
    selection_indices = [idx - 1 for idx in ml_ids]

    model_parm7 = tmp / "model.parm7"
    model_rst7 = tmp / "model.rst7"
    selection = selection_indices
    if len(selection) == len(real_top.atoms):
        shutil.copyfile(real_copy, model_parm7)
        shutil.copyfile(real_rst7, model_rst7)
    else:
        model = real_top[selection]
        model.box = None
        model.save(str(model_parm7), overwrite=True)
        model.save(str(model_rst7), overwrite=True)

    atoms_real = read(str(input_copy))
    atoms_model = read(str(model_copy))
    if len(atoms_model) != len(selection_indices):
        tmp.cleanup()
        raise ValueError(
            "model_pdb atom count does not match the detected ML-region selection from the input PDB."
        )
    for i, ridx in enumerate(selection_indices):
        atoms_model[i].position = atoms_real[ridx].position

    atoms_model_lh = _append_link_hydrogens(atoms_model, atoms_real, link_pairs)

    return MLRegionWorkspace(
        tmpdir=tmpdir,
        input_pdb=input_copy,
        real_parm7=real_copy,
        real_rst7=real_rst7,
        model_pdb=model_copy,
        model_parm7=model_parm7,
        model_rst7=model_rst7,
        selection_indices=selection_indices,
        link_pairs=link_pairs,
        atoms_real=atoms_real,
        atoms_model=atoms_model,
        atoms_model_lh=atoms_model_lh,
    )


def _hartree_to_kcalmol(Eh: float) -> float:
    return float(Eh * HARTREE_TO_KCALMOL)


class FlowList(list):
    pass


def _flow_seq_representer(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)


yaml.SafeDumper.add_representer(FlowList, _flow_seq_representer)


def _format_row_for_echo(row: List[Any]) -> str:
    def _fmt(x: Any) -> str:
        if x is None:
            return "null"
        if isinstance(x, float):
            return f"{x:.10g}"
        return str(x)

    return "[" + ", ".join(_fmt(v) for v in row) + "]"


# ---- PySCF helpers copied from the legacy DFT script ----
def fast_iao_mullikan_spin_pop(mol, dm, iaos, verbose=None):
    import numpy
    from pyscf.lib import logger as pyscf_logger
    from pyscf.lo.iao import reference_mol
    from pyscf.scf import uhf as scf_uhf

    if verbose is None:
        verbose = pyscf_logger.DEBUG

    pmol = reference_mol(mol)
    if getattr(mol, 'pbc_intor', None):
        ovlpS = mol.pbc_intor('int1e_ovlp')
    else:
        ovlpS = mol.intor_symmetric('int1e_ovlp')

    cs = numpy.dot(iaos.T.conj(), ovlpS)
    s_iao = numpy.dot(cs, iaos)
    iao_inv = numpy.linalg.solve(s_iao, cs)

    if isinstance(dm, numpy.ndarray) and dm.ndim == 2:
        spin_pop_ao = numpy.zeros(s_iao.shape[0], dtype=float)
        Ms = numpy.zeros(pmol.natm, dtype=float)
        return spin_pop_ao, Ms

    dm_a = reduce(numpy.dot, (iao_inv, dm[0], iao_inv.conj().T))
    dm_b = reduce(numpy.dot, (iao_inv, dm[1], iao_inv.conj().T))
    return scf_uhf.mulliken_spin_pop(pmol, [dm_a, dm_b], s_iao, verbose)


def _compute_atomic_charges(mol, mf) -> Dict[str, Optional[List[float]]]:
    from pyscf.scf import hf as scf_hf
    from pyscf.lo import iao as lo_iao

    dm = mf.make_rdm1()
    S = mf.get_ovlp()
    dm_tot = dm[0] + dm[1] if (isinstance(dm, np.ndarray) and dm.ndim == 3) else dm

    try:
        _, mull_chg = scf_hf.mulliken_pop(mol, dm_tot, s=S, verbose=0)
        mull_q = np.asarray(mull_chg, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[Mulliken] WARNING: Failed to compute Mulliken charges: {e}", err=True)
        mull_q = None

    try:
        _, low_chg = scf_hf.mulliken_pop_meta_lowdin_ao(mol, dm_tot, verbose=0, s=S)
        low_q = np.asarray(low_chg, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[Löwdin] WARNING: Failed to compute meta-Löwdin charges: {e}", err=True)
        low_q = None

    iao_q: Optional[List[float]] = None
    try:
        mo = mf.mo_coeff
        mo_occ = mf.mo_occ
        if isinstance(mo, np.ndarray) and mo.ndim == 2:
            occ_idx = np.asarray(mo_occ) > 0
            orbocc = mo[:, occ_idx]
        else:
            occ_idx = np.asarray(mo_occ[0]) > 0
            orbocc = mo[0][:, occ_idx]
        iaos = lo_iao.iao(mol, orbocc, minao="minao")
        _, iao_chg = lo_iao.fast_iao_mullikan_pop(mol, dm, iaos, verbose=0)
        iao_q = np.asarray(iao_chg, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[IAO] WARNING: Failed to compute IAO charges: {e}", err=True)
        iao_q = None

    return {
        "mulliken": mull_q,
        "lowdin": low_q,
        "iao": iao_q,
    }


def _compute_atomic_spin_densities(mol, mf) -> Dict[str, Optional[List[float]]]:
    from pyscf.scf import uhf as scf_uhf
    from pyscf.lo import iao as lo_iao

    dm = mf.make_rdm1()
    S = mf.get_ovlp()
    nat = mol.natm

    if not (isinstance(dm, np.ndarray) and dm.ndim == 3):
        zeros = [0.0] * nat
        return {"mulliken": zeros, "lowdin": zeros, "iao": zeros}

    try:
        _, Ms_mull = scf_uhf.mulliken_spin_pop(mol, dm, s=S, verbose=0)
        mull = np.asarray(Ms_mull, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[Spin Mulliken] WARNING: Failed to compute Mulliken spin densities: {e}", err=True)
        mull = None

    try:
        _, Ms_low = scf_uhf.mulliken_spin_pop_meta_lowdin_ao(mol, dm, verbose=0, s=S)
        low = np.asarray(Ms_low, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[Spin Löwdin] WARNING: Failed to compute meta-Löwdin spin densities: {e}", err=True)
        low = None

    iao_ms: Optional[List[float]] = None
    try:
        mo = mf.mo_coeff
        mo_occ = mf.mo_occ
        if isinstance(mo, np.ndarray) and mo.ndim == 2:
            occ_idx = np.asarray(mo_occ) > 0
            orbocc = mo[:, occ_idx]
        else:
            occ_idx = np.asarray(mo_occ[0]) > 0
            orbocc = mo[0][:, occ_idx]
        iaos = lo_iao.iao(mol, orbocc, minao="minao")
        _, Ms_iao = fast_iao_mullikan_spin_pop(mol, dm, iaos, verbose=0)
        iao_ms = np.asarray(Ms_iao, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[Spin IAO] WARNING: Failed to compute IAO spin densities: {e}", err=True)
        iao_ms = None

    return {"mulliken": mull, "lowdin": low, "iao": iao_ms}
# -----------------------------------------------
# CLI
# -----------------------------------------------


@click.command(
    help="Single-point ML-region DFT with ML(dft)/MM energy recombination.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i",
    "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Full enzyme PDB used for ML/MM (must be .pdb).",
)
@click.option(
    "--real-parm7",
    "real_parm7",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Amber parm7 topology for the full system.",
)
@click.option(
    "--model-pdb",
    "model_pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help="PDB defining the ML region (atom IDs must match the enzyme PDB). "
         "Optional when --detect-layer is enabled.",
)
@click.option(
    "--model-indices",
    "model_indices_str",
    type=str,
    default=None,
    show_default=False,
    help="Comma-separated atom indices for the ML region (ranges allowed like 1-5). "
         "Used when --model-pdb is omitted.",
)
@click.option(
    "--model-indices-one-based/--model-indices-zero-based",
    "model_indices_one_based",
    default=True,
    show_default=True,
    help="Interpret --model-indices as 1-based (default) or 0-based.",
)
@click.option(
    "--detect-layer/--no-detect-layer",
    "detect_layer",
    default=True,
    show_default=True,
    help="Detect ML/MM layers from input PDB B-factors (B=0/10/20). "
         "If disabled, you must provide --model-pdb or --model-indices.",
)
@click.option("-q", "--charge", type=int, required=True, help="Charge of the ML region.")
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=1,
    show_default=True,
    help="Spin multiplicity (2S+1) for the ML region.",
)
@click.option(
    "--freeze-atoms",
    "freeze_atoms_text",
    type=str,
    default=None,
    help="Comma-separated 1-based indices to freeze (e.g., '1,3,5').",
)
@click.option(
    "--func-basis",
    "func_basis",
    type=str,
    default="wb97m-v/6-31g**",
    show_default=True,
    help='Exchange-correlation functional and basis set as "FUNC/BASIS".',
)
@click.option("--max-cycle", type=int, default=DFT_KW["max_cycle"], show_default=True)
@click.option("--conv-tol", type=float, default=DFT_KW["conv_tol"], show_default=True)
@click.option("--grid-level", type=int, default=DFT_KW["grid_level"], show_default=True)
@click.option(
    "--out-dir",
    type=click.Path(path_type=Path, dir_okay=True, file_okay=False),
    default=Path(DFT_KW["out_dir"]),
    show_default=True,
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Optional YAML overrides under geom/calc/mlmm/dft.",
)
def cli(
    input_path: Path,
    real_parm7: Path,
    model_pdb: Optional[Path],
    model_indices_str: Optional[str],
    model_indices_one_based: bool,
    detect_layer: bool,
    charge: int,
    spin: int,
    freeze_atoms_text: Optional[str],
    func_basis: str,
    max_cycle: int,
    conv_tol: float,
    grid_level: int,
    out_dir: Path,
    args_yaml: Optional[Path],
) -> None:
    if input_path.suffix.lower() != ".pdb":
        raise click.BadParameter("Input structure must be a PDB file for ML/MM DFT.")

    prepared_input = prepare_input_structure(input_path)
    workspace: Optional[MLRegionWorkspace] = None

    model_indices: Optional[List[int]] = None
    if model_indices_str:
        try:
            model_indices = parse_indices_string(model_indices_str, one_based=model_indices_one_based)
        except click.BadParameter as e:
            raise click.ClickException(str(e))

    try:
        yaml_cfg = load_yaml_dict(args_yaml)

        geom_kw = deepcopy(OPT_GEOM_KW)
        calc_kw = deepcopy(OPT_CALC_KW)
        dft_kw = dict(DFT_KW)

        apply_yaml_overrides(
            yaml_cfg,
            [
                (geom_kw, (("geom",),)),
                (calc_kw, (("calc",), ("mlmm",))),
                (dft_kw, (("dft",),)),
            ],
        )

        geom_kw["coord_type"] = "cart"
        geom_kw["freeze_atoms"] = _normalize_geom_freeze_opt(geom_kw.get("freeze_atoms"))
        freeze_atoms_cli = _parse_freeze_atoms_opt(freeze_atoms_text)
        merged_freeze = merge_freeze_atom_indices(geom_kw, freeze_atoms_cli)
        calc_kw["freeze_atoms"] = merged_freeze

        out_dir_path = Path(out_dir).resolve()

        layer_source_pdb = input_path.resolve()
        model_pdb_cfg = model_pdb if model_pdb is not None else calc_kw.get("model_pdb")
        model_pdb_path: Optional[Path] = None

        if detect_layer:
            try:
                model_pdb_path, layer_info = build_model_pdb_from_bfactors(layer_source_pdb, out_dir_path)
                calc_kw["use_bfactor_layers"] = True
                click.echo(
                    f"[layer] Detected B-factor layers: ML={len(layer_info.get('ml_indices', []))}, "
                    f"MovableMM={len(layer_info.get('movable_mm_indices', []))}, "
                    f"FrozenMM={len(layer_info.get('frozen_indices', []))}"
                )
            except Exception as e:
                if model_pdb_cfg is None and not model_indices:
                    raise click.ClickException(str(e))
                click.echo(f"[layer] WARNING: {e} Falling back to explicit ML region.", err=True)
                detect_layer = False

        if not detect_layer:
            if model_pdb_cfg is None and not model_indices:
                raise click.ClickException("Provide --model-pdb or --model-indices when --no-detect-layer.")
            if model_pdb_cfg is not None:
                model_pdb_path = Path(model_pdb_cfg)
            else:
                try:
                    model_pdb_path = build_model_pdb_from_indices(layer_source_pdb, out_dir_path, model_indices or [])
                except Exception as e:
                    raise click.ClickException(str(e))
            calc_kw["use_bfactor_layers"] = False

        if model_pdb_path is None:
            raise click.ClickException("Failed to resolve model PDB for the ML region.")

        calc_kw["input_pdb"] = str(input_path.resolve())
        calc_kw["real_parm7"] = str(real_parm7.resolve())
        calc_kw["model_pdb"] = str(model_pdb_path.resolve())
        merged_freeze = apply_layer_freeze_constraints(
            geom_kw,
            calc_kw,
            layer_info if detect_layer else None,
            echo_fn=click.echo,
        )
        calc_kw["model_charge"] = int(charge)
        calc_kw["model_mult"] = int(spin)

        dft_kw["conv_tol"] = float(conv_tol)
        dft_kw["max_cycle"] = int(max_cycle)
        dft_kw["grid_level"] = int(grid_level)
        dft_kw["out_dir"] = str(out_dir)

        xc, basis = _parse_func_basis(func_basis)

        click.echo(pretty_block("geom", format_freeze_atoms_for_echo(geom_kw, key="freeze_atoms")))
        click.echo(pretty_block("calc", {k: calc_kw[k] for k in sorted(calc_kw.keys()) if k not in {"freeze_atoms"}}))
        click.echo(pretty_block("dft", {
            "charge": charge,
            "multiplicity": spin,
            "xc": xc,
            "basis": basis,
            "conv_tol": dft_kw["conv_tol"],
            "max_cycle": dft_kw["max_cycle"],
            "grid_level": dft_kw["grid_level"],
            "out_dir": dft_kw["out_dir"],
        }))

        time_start = time.perf_counter()
        out_dir_path = Path(dft_kw["out_dir"]).resolve()
        out_dir_path.mkdir(parents=True, exist_ok=True)

        workspace = _prepare_ml_region_workspace(
            input_pdb=Path(calc_kw["input_pdb"]),
            real_parm7=Path(calc_kw["real_parm7"]),
            model_pdb=Path(calc_kw["model_pdb"]),
            link_mlmm=calc_kw.get("link_mlmm"),
        )

        xyz_path = out_dir_path / "ml_region_with_linkH.xyz"
        xyz_path.write_text(_atoms_to_xyz_string(workspace.atoms_model_lh, "ML region + link-H"))
        click.echo(f"[write] Wrote '{xyz_path}'.")

        try:
            from pyscf import gto
        except Exception as exc:
            raise click.ClickException(f"PySCF import failed: {exc}") from exc

        mol = gto.Mole()
        mol.verbose = int(dft_kw.get("verbose", 4))
        mol.build(
            atom=_atoms_to_pyscf_atoms(workspace.atoms_model_lh),
            unit="Angstrom",
            charge=int(charge),
            spin=int(spin - 1),
            basis=basis,
        )

        using_gpu = False
        engine_label = "pyscf(cpu)"
        try:
            import gpu4pyscf

            gpu4pyscf.activate()
            from gpu4pyscf import dft as gdf

            mf = gdf.RKS(mol) if spin == 1 else gdf.UKS(mol)
            using_gpu = True
            engine_label = "gpu4pyscf"
        except Exception:
            from pyscf import dft as pdft

            mf = pdft.RKS(mol) if spin == 1 else pdft.UKS(mol)

        mf.xc = xc
        mf.max_cycle = int(dft_kw["max_cycle"])
        mf.conv_tol = float(dft_kw["conv_tol"])
        try:
            mf.grids.level = int(dft_kw["grid_level"])
        except Exception as exc:
            click.echo(f"[grids] WARNING: Could not set grids.level={dft_kw['grid_level']}: {exc}", err=True)
        try:
            mf.chkfile = None
        except Exception:
            pass
        if xc.lower().endswith("-v") or "vv10" in xc.lower():
            mf.nlc = "vv10"

        click.echo("\n=== ML-region DFT single-point started ===\n")
        tic_scf = time.time()
        e_tot = mf.kernel()
        toc_scf = time.time()
        click.echo("\n=== ML-region DFT single-point finished ===\n")

        converged = bool(getattr(mf, "converged", False))
        if e_tot is None:
            e_tot = float(getattr(mf, "e_tot", np.nan))
        e_h = float(e_tot)
        e_kcal = _hartree_to_kcalmol(e_h)

        charges = _compute_atomic_charges(mol, mf)
        spins = _compute_atomic_spin_densities(mol, mf)

        def _round(xs: Optional[List[float]]) -> Optional[List[float]]:
            if xs is None:
                return None
            return [0.0 if (x == x) and abs(x) < 1e-10 else float(x) for x in xs]

        charges = {k: _round(v) for k, v in charges.items()}
        spins = {k: _round(v) for k, v in spins.items()}

        charges_table: List[List[Any]] = []
        spins_table: List[List[Any]] = []
        for i in range(mol.natm):
            elem = mol.atom_symbol(i)
            charges_table.append([
                i,
                elem,
                None if charges["mulliken"] is None else charges["mulliken"][i],
                None if charges["lowdin"] is None else charges["lowdin"][i],
                None if charges["iao"] is None else charges["iao"][i],
            ])
            spins_table.append([
                i,
                elem,
                None if spins["mulliken"] is None else spins["mulliken"][i],
                None if spins["lowdin"] is None else spins["lowdin"][i],
                None if spins["iao"] is None else spins["iao"][i],
            ])

        click.echo("\ncharges [index, element, mulliken, lowdin, iao]:")
        for row in charges_table:
            click.echo(f"- {_format_row_for_echo(row)}")

        click.echo("\nspin_densities [index, element, mulliken, lowdin, iao]:")
        for row in spins_table:
            click.echo(f"- {_format_row_for_echo(row)}")

        mm_device = calc_kw.get("mm_device", "cpu")
        mm_cuda_idx = int(calc_kw.get("mm_cuda_idx", 0))
        mm_threads = int(calc_kw.get("mm_threads", 16))

        atoms_real = workspace.atoms_real.copy()
        atoms_model = workspace.atoms_model.copy()

        calc_real = hessianffCalculator(
            parm7=str(workspace.real_parm7),
            rst7=str(workspace.real_rst7),
            device=mm_device,
            cuda_idx=mm_cuda_idx,
            threads=mm_threads,
        )
        atoms_real.calc = calc_real
        e_real_low = atoms_real.get_potential_energy()

        calc_model = hessianffCalculator(
            parm7=str(workspace.model_parm7),
            rst7=str(workspace.model_rst7),
            device=mm_device,
            cuda_idx=mm_cuda_idx,
            threads=mm_threads,
        )
        atoms_model.calc = calc_model
        e_model_low = atoms_model.get_potential_energy()

        e_real_low_au = e_real_low * EV2AU
        e_model_low_au = e_model_low * EV2AU
        e_total_au = e_real_low_au + e_h - e_model_low_au
        e_total_kcal = _hartree_to_kcalmol(e_total_au)

        result_yaml = {
            "input": {
                "charge": int(charge),
                "multiplicity": int(spin),
                "xc": xc,
                "basis": basis,
                "conv_tol": dft_kw["conv_tol"],
                "max_cycle": dft_kw["max_cycle"],
                "grid_level": dft_kw["grid_level"],
                "out_dir": str(out_dir_path),
            },
            "energy": {
                "hartree": e_h,
                "kcal_per_mol": e_kcal,
                "converged": converged,
                "scf_time_sec": round(toc_scf - tic_scf, 3),
                "engine": engine_label,
                "used_gpu": bool(using_gpu),
            },
            "mlmm_energy": {
                "E_real_low_eV": e_real_low,
                "E_model_low_eV": e_model_low,
                "E_real_low_hartree": e_real_low_au,
                "E_model_low_hartree": e_model_low_au,
                "E_total_ml_dft_mm_hartree": e_total_au,
                "E_total_ml_dft_mm_kcal_per_mol": e_total_kcal,
            },
            "charges [index, element, mulliken, lowdin, iao]": [FlowList(r) for r in charges_table],
            "spin_densities [index, element, mulliken, lowdin, iao]": [FlowList(r) for r in spins_table],
        }

        result_file = out_dir_path / "result.yaml"
        result_file.write_text(yaml.safe_dump(result_yaml, sort_keys=False, allow_unicode=True))
        click.echo(f"[write] Wrote '{result_file}'.")

        click.echo(f"\nE_DFT (Hartree): {e_h:.12f}")
        click.echo(f"E_DFT (kcal/mol): {e_kcal:.6f}")
        click.echo(f"E_total ML(dft)/MM (Hartree): {e_total_au:.12f}")
        click.echo(f"E_total ML(dft)/MM (kcal/mol): {e_total_kcal:.6f}")

        if not converged:
            click.echo("WARNING: SCF did not converge.", err=True)
            sys.exit(3)

        click.echo(format_elapsed("[time] Elapsed Time for DFT", time_start))

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except click.ClickException:
        raise
    except Exception as exc:
        tb = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
        click.echo("Unhandled error during ML/MM DFT:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)
    finally:
        prepared_input.cleanup()
        if workspace is not None:
            workspace.cleanup()


if __name__ == "__main__":
    cli()
