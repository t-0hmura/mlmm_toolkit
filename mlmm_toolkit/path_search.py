# mlmm_toolkit/path_search.py

"""
ML/MM recursive GSM segmentation for multistep minimum-energy paths.

Example:
    mlmm path-search -i R.pdb P.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0

For detailed documentation, see: docs/path_search.md
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import sys
import traceback
import textwrap
import tempfile
import os
import shutil
import time  # timing
import re    # used in _segment_base_id

import click
import numpy as np
import yaml
from click.core import ParameterSource

from pysisyphus.helpers import geom_loader
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import AU2KCALPERMOL, BOHR2ANG

# Biopython for PDB I/O and parsing
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO

from .mlmm_calc import mlmm, MLMMCore, MLMMASECalculator
from .path_opt import GS_KW as _PATH_GS_KW, STOPT_KW as _PATH_STOPT_KW, DMF_KW as _PATH_DMF_KW, _run_dmf_mep, _select_hei_index
from .opt import (
    GEOM_KW as _OPT_GEOM_KW,
    CALC_KW as _OPT_CALC_KW,
    LBFGS_KW as _OPT_LBFGS_KW,
    _parse_freeze_atoms as _parse_freeze_atoms_opt,
    _normalize_geom_freeze as _normalize_geom_freeze_opt,
)
from .utils import (
    apply_layer_freeze_constraints,
    convert_xyz_to_pdb,
    load_yaml_dict,
    deep_update,
    apply_yaml_overrides,
    pretty_block,
    strip_inherited_keys,
    format_freeze_atoms_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    build_energy_diagram,
    prepare_input_structure,
    resolve_charge_spin_or_raise,
    PreparedInputStructure,
    parse_indices_string,
    build_model_pdb_from_bfactors,
    build_model_pdb_from_indices,
    read_bfactors_from_pdb,
    has_valid_layer_bfactors,
    parse_layer_indices_from_bfactors,
    collect_single_option_values,
)
from .cli_utils import resolve_yaml_sources, load_merged_yaml_cfg, link_or_copy_file
from .preflight import validate_existing_files
from .trj2fig import run_trj2fig  # auto-generate an energy plot when a _trj.xyz is produced
from .summary_log import write_summary_log
from .bond_changes import compare_structures, summarize_changes
from .align_freeze_atoms import align_and_refine_sequence_inplace, kabsch_R_t

# -----------------------------------------------
# Configuration defaults
# -----------------------------------------------

# Geometry (input handling) — reuse opt.py defaults
GEOM_KW: Dict[str, Any] = deepcopy(_OPT_GEOM_KW)

# ML/MM calculator settings — reuse opt.py defaults
CALC_KW: Dict[str, Any] = deepcopy(_OPT_CALC_KW)

# GrowingString (path representation)
GS_KW: Dict[str, Any] = deepcopy(_PATH_GS_KW)

# StringOptimizer (GSM optimization control)
STOPT_KW: Dict[str, Any] = deepcopy(_PATH_STOPT_KW)
STOPT_KW.update({
    "out_dir": "./result_path_search/",
})

# LBFGS settings
LBFGS_KW: Dict[str, Any] = deepcopy(_OPT_LBFGS_KW)
LBFGS_KW.update({
    "out_dir": "./result_path_search/",
})

# Covalent-bond change detection
BOND_KW: Dict[str, Any] = {
    "device": "cuda",            # str, device for UMA graph analysis during bond detection
    "bond_factor": 1.20,         # float, scaling of covalent radii for bond cutoff
    "margin_fraction": 0.05,     # float, fractional margin to tolerate small deviations
    "delta_fraction": 0.05,      # float, change threshold to flag bond formation/breaking
}

# DMF (Direct Max Flux) defaults
DMF_KW: Dict[str, Any] = deepcopy(_PATH_DMF_KW)

# Global search control
SEARCH_KW: Dict[str, Any] = {
    "max_depth": 10,               # int, recursion depth cap for segment refinement
    "stitch_rmsd_thresh": 1.0e-4,  # float, ≤ this → treat endpoints as duplicates on stitching
    "bridge_rmsd_thresh": 1.0e-4,  # float, > this → insert a bridge GSM
    "rmsd_align": True,            # bool, retained for compatibility (ignored internally)
    "max_nodes_segment": 10,       # int, max_nodes override for recursive segments
    "max_nodes_bridge": 5,         # int, max_nodes override for bridge GSMs
    "kink_max_nodes": 3,           # int, nodes for linear interpolation when skipping GSM at a kink
    "refine_mode": None,           # Optional[str], HEI seed policy (peak|minima)
}

# Multi-structure loader
def _load_structures(
    inputs: Sequence[PreparedInputStructure],
    coord_type: str,
    base_freeze: Sequence[int],
) -> List[Any]:
    """
    Load multiple geometries and assign `freeze_atoms`; return a list of geometries.
    """
    geoms: List[Any] = []
    for prepared in inputs:
        geom_path = prepared.geom_path
        g = geom_loader(geom_path, coord_type=coord_type)
        cfg: Dict[str, Any] = {"freeze_atoms": list(base_freeze)}
        freeze = merge_freeze_atom_indices(cfg)
        g.freeze_atoms = np.array(freeze, dtype=int)
        geoms.append(g)
    return geoms


# Helpers shared with opt.py for freeze parsing/normalization
_parse_freeze_atoms = _parse_freeze_atoms_opt
_normalize_geom_freeze = _normalize_geom_freeze_opt


def _write_xyz_trj_with_energy(images: Sequence, energies: Sequence[float], path: Path) -> None:
    """
    Write an XYZ `_trj.xyz` with the energy on line 2 of each block.
    """
    blocks: List[str] = []
    E = np.array(energies, dtype=float)
    for geom, e in zip(images, E):
        s = geom.as_xyz()
        lines = s.splitlines()
        if len(lines) >= 2 and lines[0].strip().isdigit():
            lines[1] = f"{e:.12f}"
        s_mod = "\n".join(lines)
        if not s_mod.endswith("\n"):
            s_mod += "\n"
        blocks.append(s_mod)
    with open(path, "w") as f:
        f.write("".join(blocks))


def _maybe_convert_to_pdb(in_path: Path, ref_pdb_path: Optional[Path], out_path: Optional[Path] = None) -> Optional[Path]:
    """
    If any input is PDB, convert the given `.xyz/_trj.xyz` to PDB using `ref_pdb_path`.
    Return the output path on success, else None.
    """
    try:
        if ref_pdb_path is None or (not in_path.exists()) or in_path.suffix.lower() not in (".xyz", "_trj.xyz"):
            return None
        out_pdb = out_path if out_path is not None else in_path.with_suffix(".pdb")
        convert_xyz_to_pdb(in_path, ref_pdb_path, out_pdb)
        click.echo(f"[convert] Wrote '{out_pdb}'.")
        return out_pdb
    except Exception as e:
        click.echo(f"[convert] WARNING: Failed to convert '{in_path.name}' to PDB: {e}", err=True)
        return None


def _kabsch_rmsd(A: np.ndarray, B: np.ndarray, align: bool = True, indices: Optional[Sequence[int]] = None) -> float:
    """
    RMSD between A and B (no rigid alignment; `align` is ignored). Optional subset selection via `indices`.
    """
    assert A.shape == B.shape and A.shape[1] == 3
    if indices is not None and len(indices) > 0:
        idx = np.array(sorted({int(i) for i in indices if 0 <= int(i) < A.shape[0]}), dtype=int)
        if idx.size == 0:
            idx = np.arange(A.shape[0], dtype=int)
        A = A[idx]
        B = B[idx]
    diff = A - B
    return float(np.sqrt((diff * diff).sum() / A.shape[0]))




def _has_bond_change(x, y, bond_cfg: Dict[str, Any]) -> Tuple[bool, str]:
    """
    Determine whether covalent bonds are forming or breaking between `x` and `y`.
    """
    res = compare_structures(
        x, y,
        device=bond_cfg.get("device", "cuda"),
        bond_factor=float(bond_cfg.get("bond_factor", 1.20)),
        margin_fraction=float(bond_cfg.get("margin_fraction", 0.05)),
        delta_fraction=float(bond_cfg.get("delta_fraction", 0.05)),
    )
    formed = len(res.formed_covalent) > 0
    broken = len(res.broken_covalent) > 0
    summary = summarize_changes(x, res, one_based=True)
    return (formed or broken), summary


# ---------- Minimal GS configuration helper ----------


# -----------------------------------------------
# Kink detection & interpolation helpers
# -----------------------------------------------

def _max_displacement_between(ga, gb, align: bool = True, indices: Optional[Sequence[int]] = None) -> float:
    """
    Maximum per-atom displacement (Å) between two structures (no alignment).
    """
    A = np.asarray(ga.coords3d, dtype=float)
    B = np.asarray(gb.coords3d, dtype=float)
    if A.shape != B.shape or A.shape[1] != 3:
        raise ValueError("Geometries must have the same number of atoms for displacement.")
    disp = np.linalg.norm(A - B, axis=1)
    return float(np.max(disp))


def _new_geom_from_coords(atoms: Sequence[str], coords: np.ndarray, coord_type: str, freeze_atoms: Sequence[int]) -> Any:
    """
    Create a pysisyphus Geometry from Bohr coords via temporary XYZ; attach `freeze_atoms`.
    """
    lines = [str(len(atoms)), ""]
    coords_ang = np.asarray(coords, dtype=float) * BOHR2ANG
    for sym, (x, y, z) in zip(atoms, coords_ang):
        lines.append(f"{sym} {x:.15f} {y:.15f} {z:.15f}")
    s = "\n".join(lines) + "\n"
    tmp = tempfile.NamedTemporaryFile("w+", suffix=".xyz", delete=False)
    try:
        tmp.write(s)
        tmp.flush()
        tmp.close()
        g = geom_loader(Path(tmp.name), coord_type=coord_type)
        g.freeze_atoms = np.array(sorted(set(map(int, freeze_atoms))), dtype=int)
        return g
    finally:
        try:
            os.unlink(tmp.name)
        except Exception:
            pass


def _make_linear_interpolations(gL, gR, n_internal: int) -> List[Any]:
    """
    Return `n_internal` linearly interpolated structures between gL → gR (excluding endpoints).
    Atom order follows `gL`.
    """
    A = np.asarray(gL.coords3d, dtype=float)
    B = np.asarray(gR.coords3d, dtype=float)
    assert A.shape == B.shape and A.shape[1] == 3, "Atom counts must match for interpolation."
    atoms = [a for a in gL.atoms]
    coord_type = gL.coord_type
    faL = getattr(gL, "freeze_atoms", np.array([], dtype=int))
    faR = getattr(gR, "freeze_atoms", np.array([], dtype=int))
    freeze_union = sorted(set(map(int, faL)) | set(map(int, faR)))
    interps: List[Any] = []
    for k in range(1, n_internal + 1):
        t = k / (n_internal + 1.0)
        C = (1.0 - t) * A + t * B
        interps.append(_new_geom_from_coords(atoms, C, coord_type, freeze_union))
    return interps


# ---- Segment/bridge tagging helpers ----

def _tag_images(images: Sequence[Any], **attrs: Any) -> None:
    """
    Attach arbitrary attributes to Geometry images.
    """
    for im in images:
        for k, v in attrs.items():
            try:
                setattr(im, k, v)
            except Exception:
                pass


def _segment_base_id(tag: str) -> str:
    """
    Extract base id 'seg_XXX' from a tag like 'seg_000_refine'; fallback to `tag` or 'seg'.
    """
    m = re.search(r"(seg_\d{3})", tag or "")
    return m.group(1) if m else (tag or "seg")


def _is_local_minimum(idx: int, energies: Sequence[float]) -> bool:
    if idx < 0 or idx >= len(energies):
        return False
    if idx == 0:
        return len(energies) > 1 and energies[1] > energies[0]
    if idx == len(energies) - 1:
        return energies[-2] > energies[-1]
    return energies[idx - 1] > energies[idx] and energies[idx + 1] > energies[idx]


def _find_nearest_local_minimum(
    hei_idx: int,
    direction: int,
    energies: Sequence[float],
) -> Optional[int]:
    i = hei_idx + direction
    while 0 <= i < len(energies):
        if _is_local_minimum(i, energies):
            return i
        i += direction
    return None


@dataclass
class GSMResult:
    images: List[Any]
    energies: List[float]
    hei_idx: int


# ---- Per-segment summary for the console report ----
@dataclass
class SegmentReport:
    tag: str
    barrier_kcal: float
    delta_kcal: float
    summary: str  # summarize_changes string (empty for bridges)
    kind: str = "seg"          # "seg" or "bridge"
    seg_index: int = 0         # 1-based index along final MEP (assigned later)


def _run_gsm_between(
    gA,
    gB,
    shared_calc,
    gs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],  # reference PDB for conversion
) -> GSMResult:
    """
    Run GSM between `gA`–`gB`, save segment outputs, and return images/energies/HEI index.
    """
    # Attach calculator to endpoints
    for g in (gA, gB):
        g.set_calculator(shared_calc)

    gs = GrowingString(
        images=[gA, gB],
        calc_getter=(lambda: shared_calc),
        **gs_cfg,
    )

    _opt_args = dict(opt_cfg)
    seg_dir = out_dir / f"{tag}_mep"
    seg_dir.mkdir(parents=True, exist_ok=True)
    _opt_args["out_dir"] = str(seg_dir)

    optimizer = StringOptimizer(
        geometry=gs,
        **{k: v for k, v in _opt_args.items() if k != "type"}
    )

    click.echo(f"\n=== [{tag}] GSM started ===\n")
    optimizer.run()
    click.echo(f"\n=== [{tag}] GSM finished ===\n")

    energies = list(map(float, np.array(gs.energy, dtype=float)))
    images = list(gs.images)

    # Choose HEI: prefer internal local maxima; fallback to highest internal node
    E = np.array(energies, dtype=float)
    nE = len(E)
    local_max_candidates = [i for i in range(1, nE - 1) if (E[i] > E[i - 1] and E[i] > E[i + 1])]
    if local_max_candidates:
        hei_idx = int(max(local_max_candidates, key=lambda i: E[i]))
    else:
        hei_idx = int(np.argmax(E[1:-1])) + 1 if nE >= 3 else int(np.argmax(E))

    # Write trajectory
    final_trj = seg_dir / "final_geometries_trj.xyz"
    wrote_with_energy = True
    try:
        _write_xyz_trj_with_energy(images, energies, final_trj)
        click.echo(f"[{tag}] Wrote '{final_trj}'.")
    except Exception:
        wrote_with_energy = False
        with open(final_trj, "w") as f:
            f.write(gs.as_xyz())
        click.echo(f"[{tag}] Wrote '{final_trj}'.")

    # Energy plot for the segment
    try:
        if wrote_with_energy:
            run_trj2fig(final_trj, [seg_dir / "mep_plot.png"], unit="kcal", reference="init", reverse_x=False)
            click.echo(f"[{tag}] Saved energy plot → '{seg_dir / 'mep_plot.png'}'")
        else:
            click.echo(f"[{tag}] WARNING: Energies missing; skipping plot.", err=True)
    except Exception as e:
        click.echo(f"[{tag}] WARNING: Failed to plot energy: {e}", err=True)

    # If PDB input exists, convert intermediate _trj.xyz to PDB
    _maybe_convert_to_pdb(final_trj, ref_pdb_path, seg_dir / "final_geometries.pdb")

    # Write HEI structure (XYZ with energy in line 2)
    try:
        hei_geom = images[hei_idx]
        hei_E = float(E[hei_idx])
        hei_xyz = seg_dir / "hei.xyz"
        s = hei_geom.as_xyz()
        lines = s.splitlines()
        if len(lines) >= 2 and lines[0].strip().isdigit():
            lines[1] = f"{hei_E:.12f}"
            s_out = "\n".join(lines)
            if not s_out.endswith("\n"):
                s_out += "\n"
        else:
            s_out = s if s.endswith("\n") else (s + "\n")
        with open(hei_xyz, "w") as f:
            f.write(s_out)
        click.echo(f"[{tag}] Wrote '{hei_xyz}'.")
        _maybe_convert_to_pdb(hei_xyz, ref_pdb_path, seg_dir / "hei.pdb")
    except Exception as e:
        click.echo(f"[{tag}] WARNING: Failed to write HEI structure: {e}", err=True)

    return GSMResult(images=images, energies=energies, hei_idx=hei_idx)


def _run_dmf_between(
    gA,
    gB,
    shared_calc,
    calc_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],
    max_nodes: int,
    dmf_cfg: Dict[str, Any],
) -> GSMResult:
    """Run DMF for a segment and convert outputs to pysisyphus Geometries."""
    from pysisyphus.constants import ANG2BOHR
    from ase.io import read as ase_read, write as ase_write
    from io import StringIO

    seg_dir = out_dir / f"{tag}_mep"
    seg_dir.mkdir(parents=True, exist_ok=True)

    fix_atoms: List[int] = []
    try:
        fix_atoms = sorted(
            {int(i) for g in [gA, gB] for i in getattr(g, "freeze_atoms", [])}
        )
    except Exception:
        pass

    # Convert pysisyphus geometries to ASE Atoms for DMF
    def _geom_to_ase(g):
        return ase_read(StringIO(g.as_xyz()), format="xyz")

    geoms_for_dmf = [gA, gB]

    try:
        from ase.calculators.mixing import SumCalculator
        from dmf import DirectMaxFlux, interpolate_fbenm
    except Exception as e:
        raise RuntimeError(f"DMF mode requires pydmf and cyipopt: {e}") from e

    from .harmonic_constraints import HarmonicFixAtoms
    from .utils import deep_update, convert_xyz_to_pdb

    ref_images = [_geom_to_ase(g) for g in geoms_for_dmf]
    charge = int(calc_cfg.get("model_charge", 0))
    spin = int(calc_cfg.get("model_mult", 1))
    for img in ref_images:
        img.info["charge"] = charge
        img.info["spin"] = spin

    # Build ASE calculator from the shared PySisyphus calculator
    ase_calc = MLMMASECalculator(core=shared_calc.core)

    dmf_cfg_local = deep_update(dict(DMF_KW), dmf_cfg)
    fbenm_opts = dict(dmf_cfg_local.get("fbenm_options", {}))
    cfbenm_opts = dict(dmf_cfg_local.get("cfbenm_options", {}))
    dmf_opts = dict(dmf_cfg_local.get("dmf_options", {}))
    update_teval = bool(dmf_opts.pop("update_teval", False))
    k_fix = float(dmf_cfg_local.get("k_fix", 300.0))

    mxflx_fbenm = interpolate_fbenm(
        ref_images,
        nmove=max(1, int(max_nodes)),
        fbenm_only_endpoints=bool(dmf_cfg_local.get("fbenm_only_endpoints", False)),
        correlated=bool(dmf_cfg_local.get("correlated", False)),
        sequential=bool(dmf_cfg_local.get("sequential", False)),
        output_file=str(seg_dir / "dmf_fbenm_ipopt.out"),
        fbenm_options=fbenm_opts,
        cfbenm_options=cfbenm_opts,
        dmf_options=dmf_opts,
    )
    coefs = mxflx_fbenm.coefs.copy()

    mxflx = DirectMaxFlux(
        ref_images,
        coefs=coefs,
        nmove=max(1, int(max_nodes)),
        update_teval=update_teval,
        remove_rotation_and_translation=bool(dmf_opts.get("remove_rotation_and_translation", False)),
        mass_weighted=bool(dmf_opts.get("mass_weighted", False)),
        parallel=bool(dmf_opts.get("parallel", False)),
        eps_vel=float(dmf_opts.get("eps_vel", 0.01)),
        eps_rot=float(dmf_opts.get("eps_rot", 0.01)),
        beta=float(dmf_opts.get("beta", 10.0)),
    )

    for image in mxflx.images:
        if "charge" not in image.info:
            image.info["charge"] = charge
        if "spin" not in image.info:
            image.info["spin"] = spin
        if fix_atoms:
            ref_positions = image.get_positions()[fix_atoms]
            harmonic_calc = HarmonicFixAtoms(indices=fix_atoms, ref_positions=ref_positions, k_fix=k_fix)
            image.calc = SumCalculator([ase_calc, harmonic_calc])
        else:
            image.calc = ase_calc

    mxflx.add_ipopt_options({"output_file": str(seg_dir / "dmf_ipopt.out")})
    max_cycles_dmf = dmf_cfg_local.get("max_cycles")
    if max_cycles_dmf is not None:
        try:
            max_iter = int(max_cycles_dmf)
            if max_iter > 0:
                mxflx.add_ipopt_options({"max_iter": max_iter})
        except Exception:
            pass

    click.echo(f"\n=== [{tag}] DMF started ===\n")
    mxflx.solve(tol="tight")
    click.echo(f"\n=== [{tag}] DMF finished ===\n")

    # Evaluate energies using PySisyphus calculator
    energies = []
    for image in mxflx.images:
        elems = image.get_chemical_symbols()
        coords_bohr = np.asarray(image.get_positions(), dtype=float).reshape(-1, 3) * ANG2BOHR
        energies.append(float(shared_calc.get_energy(elems, coords_bohr)["energy"]))
    hei_idx = _select_hei_index(energies)

    # Write trajectory
    final_trj = seg_dir / "final_geometries_trj.xyz"
    _write_xyz_trj_with_energy_from_ase(mxflx.images, energies, final_trj)
    click.echo(f"[{tag}] Wrote '{final_trj}'.")
    _maybe_convert_to_pdb(final_trj, ref_pdb_path, seg_dir / "final_geometries.pdb")

    try:
        run_trj2fig(final_trj, [seg_dir / "mep_plot.png"], unit="kcal", reference="init", reverse_x=False)
        click.echo(f"[{tag}] Saved energy plot → '{seg_dir / 'mep_plot.png'}'")
    except Exception as e:
        click.echo(f"[{tag}] WARNING: Failed to plot energy: {e}", err=True)

    # Convert ASE images back to pysisyphus Geometries
    from pysisyphus.helpers import geom_loader as gl
    imgs = []
    for atoms in mxflx.images:
        buf = StringIO()
        ase_write(buf, atoms, format="xyz")
        buf.seek(0)
        # Write temp xyz and load as geom
        tmp_xyz = seg_dir / f"_tmp_dmf_{len(imgs)}.xyz"
        with open(tmp_xyz, "w") as f:
            f.write(buf.getvalue())
        g = gl(tmp_xyz, coord_type=gA.coord_type)
        try:
            g.freeze_atoms = np.array(getattr(gA, "freeze_atoms", []), dtype=int)
        except Exception:
            pass
        g.set_calculator(shared_calc)
        imgs.append(g)
        tmp_xyz.unlink(missing_ok=True)

    return GSMResult(images=imgs, energies=energies, hei_idx=hei_idx)


def _write_xyz_trj_with_energy_from_ase(images, energies, path: Path) -> None:
    """Write an ASE Atoms list with energies as an XYZ trajectory."""
    from ase.io import write as ase_write
    from io import StringIO
    blocks = []
    for atoms, E in zip(images, energies):
        buf = StringIO()
        ase_write(buf, atoms, format="xyz")
        s = buf.getvalue()
        lines = s.splitlines()
        if len(lines) >= 2 and lines[0].strip().isdigit():
            lines[1] = f"{E:.12f}"
        blocks.append("\n".join(lines) + "\n")
    with open(path, "w") as f:
        f.write("".join(blocks))


def _run_mep_between(
    gA,
    gB,
    shared_calc,
    gs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],
    mep_mode_kind: str = "gsm",
    calc_cfg: Optional[Dict[str, Any]] = None,
    max_nodes: int = 10,
    dmf_cfg: Optional[Dict[str, Any]] = None,
) -> GSMResult:
    """Dispatcher: run GSM or DMF between two geometries."""
    if mep_mode_kind == "dmf":
        return _run_dmf_between(
            gA, gB, shared_calc,
            calc_cfg=calc_cfg or {},
            out_dir=out_dir, tag=tag,
            ref_pdb_path=ref_pdb_path,
            max_nodes=max_nodes,
            dmf_cfg=dmf_cfg or dict(DMF_KW),
        )
    return _run_gsm_between(gA, gB, shared_calc, gs_cfg, opt_cfg, out_dir, tag=tag, ref_pdb_path=ref_pdb_path)


def _optimize_single(
    g,
    shared_calc,
    lbfgs_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],  # for PDB conversion
):
    """
    Run single-structure optimization (LBFGS) and return the final Geometry.
    """
    g.set_calculator(shared_calc)

    seg_dir = out_dir / f"{tag}_lbfgs_opt"
    seg_dir.mkdir(parents=True, exist_ok=True)
    args = dict(lbfgs_cfg)
    args["out_dir"] = str(seg_dir)

    opt = LBFGS(g, **args)

    click.echo(f"\n=== [{tag}] Single-structure LBFGS started ===\n")
    opt.run()
    click.echo(f"\n=== [{tag}] Single-structure LBFGS finished ===\n")

    try:
        final_xyz = Path(opt.final_fn) if isinstance(opt.final_fn, (str, Path)) else Path(opt.final_fn)
        _maybe_convert_to_pdb(final_xyz, ref_pdb_path)
        g_final = geom_loader(final_xyz, coord_type=g.coord_type)
        try:
            g_final.freeze_atoms = np.array(getattr(g, "freeze_atoms", []), dtype=int)
        except Exception:
            pass
        g_final.set_calculator(shared_calc)
        return g_final
    except Exception:
        return g


def _refine_between(
    gL,
    gR,
    shared_calc,
    gs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],  # for PDB conversion
    mep_mode_kind: str = "gsm",
    calc_cfg: Optional[Dict[str, Any]] = None,
    max_nodes: int = 10,
    dmf_cfg: Optional[Dict[str, Any]] = None,
) -> GSMResult:
    """
    Refine End1–End2 via GSM or DMF (force climb=True for GSM).
    """
    gs_refine_cfg = {**gs_cfg, "climb": True, "climb_lanczos": True}
    return _run_mep_between(
        gL, gR, shared_calc, gs_refine_cfg, opt_cfg, out_dir, tag=f"{tag}_refine",
        ref_pdb_path=ref_pdb_path, mep_mode_kind=mep_mode_kind,
        calc_cfg=calc_cfg, max_nodes=max_nodes, dmf_cfg=dmf_cfg,
    )


def _maybe_bridge_segments(
    tail_g,
    head_g,
    shared_calc,
    gs_cfg: Dict[str, Any],  # bridge-specific GS config
    opt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    rmsd_thresh: float,
    ref_pdb_path: Optional[Path],  # for PDB conversion
    mep_mode_kind: str = "gsm",
    calc_cfg: Optional[Dict[str, Any]] = None,
    max_nodes: int = 5,
    dmf_cfg: Optional[Dict[str, Any]] = None,
) -> Optional[GSMResult]:
    """
    Run a bridge GSM/DMF if two segment endpoints are farther than the threshold.
    """
    rmsd = _kabsch_rmsd(np.array(tail_g.coords3d), np.array(head_g.coords3d), align=False)
    if rmsd <= rmsd_thresh:
        return None
    click.echo(f"[{tag}] Gap detected between segments (RMSD={rmsd:.4e} Å) — bridging via {mep_mode_kind.upper()}.")
    return _run_mep_between(
        tail_g, head_g, shared_calc, gs_cfg, opt_cfg, out_dir, tag=f"{tag}_bridge",
        ref_pdb_path=ref_pdb_path, mep_mode_kind=mep_mode_kind,
        calc_cfg=calc_cfg, max_nodes=max_nodes, dmf_cfg=dmf_cfg,
    )


def _stitch_paths(
    parts: List[Tuple[List[Any], List[float]]],
    stitch_rmsd_thresh: float,
    bridge_rmsd_thresh: float,
    shared_calc,
    gs_cfg,   # GS config for bridges (climb=False, max_nodes=search.max_nodes_bridge)
    opt_cfg,
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],  # for PDB conversion
    bond_cfg: Optional[Dict[str, Any]] = None,  # detect bond changes between adjacent parts
    segment_builder: Optional[Callable[[Any, Any, str], "CombinedPath"]] = None,  # builds a recursive segment
    segments_out: Optional[List["SegmentReport"]] = None,  # append inserted segment summaries in order
    bridge_pair_index: Optional[int] = None,   # pair index to tag bridge frames across pairs
    mep_mode_kind: str = "gsm",
    calc_cfg: Optional[Dict[str, Any]] = None,
    dmf_cfg: Optional[Dict[str, Any]] = None,
) -> Tuple[List[Any], List[float]]:
    """
    Concatenate path parts (images, energies). Insert bridge GSMs when needed.
    If covalent changes are detected across an interface, build and insert a *new* recursive segment
    using `segment_builder` instead of bridging. Update `segments_out` accordingly.
    """
    all_imgs: List[Any] = []
    all_E: List[float] = []

    def _last_known_seg_tag_from_images(imgs: List[Any]) -> Optional[str]:
        for im in reversed(imgs):
            t = getattr(im, "mep_seg_tag", None)
            if t:
                return t
        return None

    def _first_known_seg_tag_from_images(imgs: List[Any]) -> Optional[str]:
        for im in imgs:
            t = getattr(im, "mep_seg_tag", None)
            if t:
                return t
        return None

    def append_part(imgs: List[Any], Es: List[float]) -> None:
        nonlocal all_imgs, all_E
        if not imgs:
            return
        if not all_imgs:
            all_imgs.extend(imgs)
            all_E.extend(Es)
            return
        tail = all_imgs[-1]
        head = imgs[0]

        adj_changed, adj_summary = False, ""
        if segment_builder is not None and bond_cfg is not None:
            try:
                adj_changed, adj_summary = _has_bond_change(tail, head, bond_cfg)
            except Exception:
                adj_changed, adj_summary = False, ""

        if adj_changed and segment_builder is not None:
            click.echo(f"[{tag}] Covalent changes detected at interface — inserting a new recursive segment.")
            if adj_summary:
                click.echo(textwrap.indent(adj_summary, prefix="  "))
            sub = segment_builder(tail, head, f"{tag}_mid")
            seg_imgs, seg_E = sub.images, sub.energies
            if segments_out is not None and getattr(sub, "segments", None):
                segments_out.extend(sub.segments)
            if seg_imgs:
                if _kabsch_rmsd(np.array(all_imgs[-1].coords3d), np.array(seg_imgs[0].coords3d), align=False) <= stitch_rmsd_thresh:
                    seg_imgs = seg_imgs[1:]
                    seg_E = seg_E[1:]
                all_imgs.extend(seg_imgs)
                all_E.extend(seg_E)
            if _kabsch_rmsd(np.array(all_imgs[-1].coords3d), np.array(imgs[0].coords3d), align=False) <= stitch_rmsd_thresh:
                imgs = imgs[1:]
                Es = Es[1:]
            all_imgs.extend(imgs)
            all_E.extend(Es)
            return

        rmsd = _kabsch_rmsd(np.array(tail.coords3d), np.array(head.coords3d), align=False)
        if rmsd <= stitch_rmsd_thresh:
            all_imgs.extend(imgs[1:])
            all_E.extend(Es[1:])
        elif rmsd > bridge_rmsd_thresh:
            left_tag_recent = _last_known_seg_tag_from_images(all_imgs) or "segL"
            right_tag_upcoming = _first_known_seg_tag_from_images(imgs) or "segR"
            left_base = _segment_base_id(left_tag_recent)
            right_base = _segment_base_id(right_tag_upcoming)
            bridge_name_base = f"{left_base}_{right_base}"

            br = _maybe_bridge_segments(
                tail, head, shared_calc, gs_cfg, opt_cfg, out_dir, tag=bridge_name_base,
                rmsd_thresh=bridge_rmsd_thresh, ref_pdb_path=ref_pdb_path,
                mep_mode_kind=mep_mode_kind, calc_cfg=calc_cfg, dmf_cfg=dmf_cfg,
            )
            if br is not None:
                _tag_images(br.images, mep_seg_tag=f"{bridge_name_base}_bridge", mep_seg_kind="bridge",
                            mep_has_bond_changes=False, pair_index=bridge_pair_index)
                b_imgs, b_E = br.images, br.energies
                if _kabsch_rmsd(np.array(all_imgs[-1].coords3d), np.array(b_imgs[0].coords3d), align=False) <= stitch_rmsd_thresh:
                    b_imgs = b_imgs[1:]
                    b_E = b_E[1:]
                if b_imgs:
                    all_imgs.extend(b_imgs)
                    all_E.extend(b_E)

                if segments_out is not None:
                    try:
                        barrier_kcal = (max(br.energies) - br.energies[0]) * AU2KCALPERMOL
                        delta_kcal = (br.energies[-1] - br.energies[0]) * AU2KCALPERMOL
                    except Exception:
                        barrier_kcal = float("nan")
                        delta_kcal = float("nan")
                    bridge_report = SegmentReport(
                        tag=f"{bridge_name_base}_bridge",
                        barrier_kcal=float(barrier_kcal),
                        delta_kcal=float(delta_kcal),
                        summary="",
                        kind="bridge"
                    )
                    insert_pos: Optional[int] = None
                    try:
                        for j, sr in enumerate(segments_out):
                            if sr.tag == right_tag_upcoming:
                                insert_pos = j
                                break
                    except Exception:
                        insert_pos = None
                    if insert_pos is None:
                        segments_out.append(bridge_report)
                    else:
                        segments_out.insert(insert_pos, bridge_report)

            if _kabsch_rmsd(np.array(all_imgs[-1].coords3d), np.array(imgs[0].coords3d), align=False) <= stitch_rmsd_thresh:
                imgs = imgs[1:]
                Es = Es[1:]
            all_imgs.extend(imgs)
            all_E.extend(Es)
        else:
            all_imgs.extend(imgs)
            all_E.extend(Es)

    for (imgs, Es) in parts:
        append_part(imgs, Es)

    return all_imgs, all_E


# -----------------------------------------------
# Recursive search (core)
# -----------------------------------------------

@dataclass
class CombinedPath:
    images: List[Any]
    energies: List[float]
    segments: List[SegmentReport]  # segment summaries in final output order


def _build_multistep_path(
    gA,
    gB,
    shared_calc,
    geom_cfg: Dict[str, Any],
    gs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    single_opt_cfg: Dict[str, Any],
    bond_cfg: Dict[str, Any],
    search_cfg: Dict[str, Any],
    refine_mode_kind: str,
    out_dir: Path,
    ref_pdb_path: Optional[Path],
    depth: int,
    seg_counter: List[int],
    branch_tag: str,
    pair_index: Optional[int] = None,
    mep_mode_kind: str = "gsm",
    calc_cfg: Optional[Dict[str, Any]] = None,
    dmf_cfg: Optional[Dict[str, Any]] = None,
) -> CombinedPath:
    """
    Recursively construct a multistep MEP from A–B and return it (A→B order).
    """
    seg_max_nodes = int(search_cfg.get("max_nodes_segment", gs_cfg.get("max_nodes", 10)))
    gs_seg_cfg = {**gs_cfg, "max_nodes": seg_max_nodes}

    if depth > int(search_cfg.get("max_depth", 10)):
        click.echo(f"[{branch_tag}] Reached maximum recursion depth. Returning current endpoints only.")
        gsm = _run_mep_between(
            gA, gB, shared_calc, gs_seg_cfg, opt_cfg, out_dir, tag=f"seg_{seg_counter[0]:03d}_maxdepth",
            ref_pdb_path=ref_pdb_path, mep_mode_kind=mep_mode_kind,
            calc_cfg=calc_cfg, max_nodes=seg_max_nodes, dmf_cfg=dmf_cfg,
        )
        seg_counter[0] += 1
        _tag_images(gsm.images, pair_index=pair_index)
        return CombinedPath(images=gsm.images, energies=gsm.energies, segments=[])

    seg_id = seg_counter[0]
    seg_counter[0] += 1
    tag0 = f"seg_{seg_id:03d}"

    gs_seg_cfg_first = {**gs_seg_cfg, "climb": True, "climb_lanczos": True}
    gsm0 = _run_mep_between(
        gA, gB, shared_calc, gs_seg_cfg_first, opt_cfg, out_dir, tag=tag0,
        ref_pdb_path=ref_pdb_path, mep_mode_kind=mep_mode_kind,
        calc_cfg=calc_cfg, max_nodes=seg_max_nodes, dmf_cfg=dmf_cfg,
    )

    hei = int(gsm0.hei_idx)
    if not (1 <= hei <= len(gsm0.images) - 2):
        click.echo(f"[{tag0}] WARNING: HEI is at an endpoint (idx={hei}). Returning the raw GSM path.")
        _tag_images(gsm0.images, pair_index=pair_index)
        return CombinedPath(images=gsm0.images, energies=gsm0.energies, segments=[])

    if refine_mode_kind == "minima":
        left_idx = _find_nearest_local_minimum(hei_idx=hei, direction=-1, energies=gsm0.energies)
        right_idx = _find_nearest_local_minimum(hei_idx=hei, direction=1, energies=gsm0.energies)
        if left_idx is None:
            left_idx = hei - 1
        if right_idx is None:
            right_idx = hei + 1
        click.echo(f"[{tag0}] Using nearest local minima around HEI (left idx={left_idx}, right idx={right_idx}).")
        left_img = gsm0.images[left_idx]
        right_img = gsm0.images[right_idx]
    else:
        left_img = gsm0.images[hei - 1]
        right_img = gsm0.images[hei + 1]
        click.echo(f"[{tag0}] Refining HEI±1 (peak mode).")

    left_end = _optimize_single(left_img, shared_calc, single_opt_cfg, out_dir, tag=f"{tag0}_left", ref_pdb_path=ref_pdb_path)
    right_end = _optimize_single(right_img, shared_calc, single_opt_cfg, out_dir, tag=f"{tag0}_right", ref_pdb_path=ref_pdb_path)

    try:
        lr_changed, lr_summary = _has_bond_change(left_end, right_end, bond_cfg)
    except Exception as e:
        click.echo(f"[{tag0}] WARNING: Failed to evaluate bond changes for kink detection: {e}", err=True)
        lr_changed, lr_summary = True, ""
    use_kink = (not lr_changed)

    if use_kink:
        n_inter = int(search_cfg.get("kink_max_nodes", 3))
        click.echo(f"[{tag0}] Kink detected (no covalent changes between End1 and End2). "
                   f"Using {n_inter} linear interpolation nodes + single-structure optimizations instead of GSM.")
        inter_geoms = _make_linear_interpolations(left_end, right_end, n_inter)
        opt_inters: List[Any] = []
        for i, g_int in enumerate(inter_geoms, 1):
            g_int.set_calculator(shared_calc)
            g_opt = _optimize_single(g_int, shared_calc, single_opt_cfg, out_dir, tag=f"{tag0}_kink_int{i}", ref_pdb_path=ref_pdb_path)
            opt_inters.append(g_opt)
        step_imgs = [left_end] + opt_inters + [right_end]
        step_E = [float(img.energy) for img in step_imgs]
        ref1 = GSMResult(images=step_imgs, energies=step_E, hei_idx=int(np.argmax(step_E)))
        step_tag_for_report = f"{tag0}_kink"
    else:
        click.echo(f"[{tag0}] Kink not detected (covalent changes present between End1 and End2).")
        if lr_summary:
            click.echo(textwrap.indent(lr_summary, prefix="  "))
        ref1 = _refine_between(left_end, right_end, shared_calc, gs_seg_cfg, opt_cfg, out_dir, tag=tag0,
                               ref_pdb_path=ref_pdb_path, mep_mode_kind=mep_mode_kind,
                               calc_cfg=calc_cfg, max_nodes=seg_max_nodes, dmf_cfg=dmf_cfg)
        step_tag_for_report = f"{tag0}_refine"

    step_imgs, step_E = ref1.images, ref1.energies

    _changed, step_summary = _has_bond_change(step_imgs[0], step_imgs[-1], bond_cfg)
    _tag_images(step_imgs, mep_seg_tag=step_tag_for_report, mep_seg_kind="seg",
                mep_has_bond_changes=bool(_changed), pair_index=pair_index)

    left_changed, left_summary = _has_bond_change(gA, left_end, bond_cfg)
    right_changed, right_summary = _has_bond_change(right_end, gB, bond_cfg)

    click.echo(f"[{tag0}] Covalent changes (A vs left_end): {'Yes' if left_changed else 'No'}")
    if left_changed:
        click.echo(textwrap.indent(left_summary, prefix="  "))
    click.echo(f"[{tag0}] Covalent changes (right_end vs B): {'Yes' if right_changed else 'No'}")
    if right_changed:
        click.echo(textwrap.indent(right_summary, prefix="  "))

    try:
        barrier_kcal = (max(step_E) - step_E[0]) * AU2KCALPERMOL
        delta_kcal = (step_E[-1] - step_E[0]) * AU2KCALPERMOL
    except Exception:
        barrier_kcal = float("nan")
        delta_kcal = float("nan")

    seg_report = SegmentReport(
        tag=step_tag_for_report,
        barrier_kcal=float(barrier_kcal),
        delta_kcal=float(delta_kcal),
        summary=step_summary if _changed else "(no covalent changes detected)",
        kind="seg"
    )

    parts: List[Tuple[List[Any], List[float]]] = []
    seg_reports: List[SegmentReport] = []

    if left_changed:
        subL = _build_multistep_path(
            gA, left_end, shared_calc, geom_cfg, gs_cfg, opt_cfg,
            single_opt_cfg, bond_cfg, search_cfg, refine_mode_kind,
            out_dir, ref_pdb_path, depth + 1, seg_counter, branch_tag=f"{branch_tag}L",
            pair_index=pair_index,
            mep_mode_kind=mep_mode_kind, calc_cfg=calc_cfg, dmf_cfg=dmf_cfg,
        )
        _tag_images(subL.images, pair_index=pair_index)
        parts.append((subL.images, subL.energies))
        seg_reports.extend(subL.segments)

    parts.append((step_imgs, step_E))
    seg_reports.append(seg_report)

    if right_changed:
        subR = _build_multistep_path(
            right_end, gB, shared_calc, geom_cfg, gs_cfg, opt_cfg,
            single_opt_cfg, bond_cfg, search_cfg, refine_mode_kind,
            out_dir, ref_pdb_path, depth + 1, seg_counter, branch_tag=f"{branch_tag}R",
            pair_index=pair_index,
            mep_mode_kind=mep_mode_kind, calc_cfg=calc_cfg, dmf_cfg=dmf_cfg,
        )
        _tag_images(subR.images, pair_index=pair_index)
        parts.append((subR.images, subR.energies))
        seg_reports.extend(subR.segments)

    bridge_max_nodes = int(search_cfg.get("max_nodes_bridge", 5))
    gs_bridge_cfg = {**gs_cfg, "max_nodes": bridge_max_nodes, "climb": False, "climb_lanczos": False}

    def _segment_builder(tail_g, head_g, _tag: str) -> CombinedPath:
        sub = _build_multistep_path(
            tail_g, head_g,
            shared_calc,
            geom_cfg, gs_cfg, opt_cfg,
            single_opt_cfg,
            bond_cfg, search_cfg, refine_mode_kind,
            out_dir=out_dir,
            ref_pdb_path=ref_pdb_path,
            depth=depth + 1,
            seg_counter=seg_counter,
            branch_tag=f"{branch_tag}B",
            pair_index=pair_index,
            mep_mode_kind=mep_mode_kind, calc_cfg=calc_cfg, dmf_cfg=dmf_cfg,
        )
        _tag_images(sub.images, pair_index=pair_index)
        return sub

    stitched_imgs, stitched_E = _stitch_paths(
        parts,
        stitch_rmsd_thresh=float(search_cfg["stitch_rmsd_thresh"]),
        bridge_rmsd_thresh=float(search_cfg["bridge_rmsd_thresh"]),
        shared_calc=shared_calc,
        gs_cfg=gs_bridge_cfg,
        opt_cfg=opt_cfg,
        out_dir=out_dir,
        tag=tag0,
        ref_pdb_path=ref_pdb_path,
        bond_cfg=bond_cfg,
        segment_builder=_segment_builder,
        segments_out=seg_reports,
        bridge_pair_index=pair_index,
        mep_mode_kind=mep_mode_kind, calc_cfg=calc_cfg, dmf_cfg=dmf_cfg,
    )

    _tag_images(stitched_imgs, pair_index=pair_index)

    return CombinedPath(images=stitched_imgs, energies=stitched_E, segments=seg_reports)


# -----------------------------------------------
# Full-system merge helpers (Biopython)
# -----------------------------------------------

def _atom_key_from_res_atom(res: PDB.Residue.Residue, atom: PDB.Atom.Atom) -> Tuple[str, str, str, str, str]:
    """
    Build a key for atom identity:
    (RESNAME, RESSEQ, ICODE, CHAIN, ATOMNAME) — uppercase where applicable.
    - RESSEQ is numeric (without insertion code).
    - ICODE is '' when blank (or ' ' in PDB).
    """
    resname = (res.get_resname() or "").strip().upper()
    het, resseq, icode = res.id
    icode_txt = "" if (icode == " " or icode is None) else str(icode).strip().upper()
    resseq_txt = str(int(resseq))
    chain_id = (res.get_parent().id or "").strip().upper()
    atname = atom.get_name().strip().upper()
    return (resname, resseq_txt, icode_txt, chain_id, atname)


def _structure_to_arrays(struct: PDB.Structure.Structure) -> Tuple[np.ndarray, List[PDB.Atom.Atom], List[Tuple[str, str, str, str, str]], Dict[Tuple[str,str,str,str,str], int]]:
    """
    Extract: coordinates (Å), atom list, key list, and key→index map from a Biopython Structure.
    Keys are (RESNAME, RESSEQ, ICODE, CHAIN, ATOMNAME).
    """
    atoms: List[PDB.Atom.Atom] = [a for a in struct.get_atoms()]
    coords = np.array([a.get_coord() for a in atoms], dtype=float)
    keys: List[Tuple[str, str, str, str, str]] = []
    key2idx: Dict[Tuple[str,str,str,str,str], int] = {}
    for idx, a in enumerate(atoms):
        res = a.get_parent()
        k = _atom_key_from_res_atom(res, a)
        keys.append(k)
        if k not in key2idx:
            key2idx[k] = idx
    return coords, atoms, keys, key2idx


def _load_structures_and_chain_align(ref_paths: Sequence[Path]) -> Tuple[List[PDB.Structure.Structure], List[np.ndarray], List[List[PDB.Atom.Atom]], List[Dict[Tuple[str,str,str,str,str], int]]]:
    """
    Load all full templates and rigidly chain-align them into the coordinate frame of the first template.
    """
    parser = PDBParser(QUIET=True)
    structs: List[PDB.Structure.Structure] = [parser.get_structure(f"ref{i:02d}", str(p)) for i, p in enumerate(ref_paths)]
    N_expected = None
    coords_list: List[np.ndarray] = []
    atoms_list: List[List[PDB.Atom.Atom]] = []
    keymaps: List[Dict[Tuple[str,str,str,str,str], int]] = []

    for s in structs:
        coords, atoms, _keys, key2idx = _structure_to_arrays(s)
        if N_expected is None:
            N_expected = coords.shape[0]
        else:
            if coords.shape[0] != N_expected:
                raise click.BadParameter(f"[merge] Atom count mismatch among --ref-pdb templates: {N_expected} vs {coords.shape[0]}")
        coords_list.append(coords)
        atoms_list.append(atoms)
        keymaps.append(key2idx)

    aligned_coords: List[np.ndarray] = []
    aligned_coords.append(coords_list[0].copy())
    for j in range(1, len(coords_list)):
        P = aligned_coords[j-1]
        Q = coords_list[j]
        R, t = kabsch_R_t(P, Q)
        Qa = (Q @ R) + t
        aligned_coords.append(Qa)

    return structs, aligned_coords, atoms_list, keymaps


def _pocket_keys_from_pdb(pocket_pdb: Path) -> List[Tuple[str, str, str, str, str]]:
    """
    Return atom identity keys for a pocket PDB file.
    """
    parser = PDBParser(QUIET=True)
    st = parser.get_structure("pocket", str(pocket_pdb))
    keys: List[Tuple[str, str, str, str, str]] = []
    for a in st.get_atoms():
        res = a.get_parent()
        k = _atom_key_from_res_atom(res, a)
        keys.append(k)
    return keys


def _write_model_block(structure: PDB.Structure.Structure,
                       remark_lines: List[str]) -> str:
    """
    Render a single MODEL block (without 'MODEL/ENDMDL') with provided REMARK lines.
    """
    io = PDBIO()
    io.set_structure(structure)
    from io import StringIO
    buf = StringIO()
    io.save(buf)
    body = "\n".join([ln for ln in buf.getvalue().splitlines() if ln.strip() != "END"])
    rem = ""
    for line in remark_lines:
        rem += f"REMARK   1 {line}\n"
    return rem + body + ("\n" if not body.endswith("\n") else "")


def _chunk_remark_indices(indices: List[int], width: int = 60) -> List[str]:
    """
    Wrap pocket atom indices into REMARK lines with limited width.
    """
    s = ",".join(map(str, indices))
    out: List[str] = []
    cur = ""
    for tok in s.split(","):
        add = (tok if not cur else "," + tok)
        if len(cur) + len(add) > width:
            out.append(f"POCKET_ATOM_INDICES {cur}")
            cur = tok
        else:
            cur += tok if not cur else "," + tok
    if cur:
        out.append(f"POCKET_ATOM_INDICES {cur}")
    return out


def _merge_pair_to_full(pair_images: List[Any],
                        pocket_ref_pdb: Path,
                        structA: PDB.Structure.Structure,
                        structB: PDB.Structure.Structure,
                        coordsA_aligned: np.ndarray,
                        coordsB_aligned: np.ndarray,
                        keymapA: Dict[Tuple[str,str,str,str,str], int],
                        keymapB: Dict[Tuple[str,str,str,str,str], int],
                        out_path: Optional[Path],
                        drop_first: bool = False,
                        seg_indices_for_frames: Optional[List[int]] = None,
                        seg_report_lookup: Optional[Dict[int, SegmentReport]] = None,
                        include_pocket_indices_for_first_model: bool = False) -> Tuple[List[str], List[int]]:
    """
    Merge a pocket-only trajectory for a *pair* into the corresponding full templates (A,B),
    generating MODEL blocks and (optionally) writing a PDB. Returns (blocks, 1-based active indices).
    """
    pocket_keys = _pocket_keys_from_pdb(pocket_ref_pdb)

    match_tpl_idx: List[int] = []
    for k in pocket_keys:
        ia = keymapA.get(k, None)
        ib = keymapB.get(k, None)
        if ia is None or ib is None:
            match_tpl_idx.append(-1)
        else:
            match_tpl_idx.append(int(ia))

    active_full_idx = sorted({i for i in match_tpl_idx if i >= 0})
    active_one_based = [i+1 for i in active_full_idx]

    Nfull = coordsA_aligned.shape[0]
    if Nfull != coordsB_aligned.shape[0]:
        raise click.BadParameter("[merge] Template A/B atom count mismatch.")

    atomsA: List[PDB.Atom.Atom] = [a for a in structA.get_atoms()]

    start_k = 1 if drop_first and len(pair_images) > 0 else 0

    model_blocks: List[str] = []

    M = len(pair_images)
    if M == 0:
        return model_blocks, active_one_based

    seg_idx_seq: List[int] = []
    if seg_indices_for_frames is not None and len(seg_indices_for_frames) == M:
        seg_idx_seq = seg_indices_for_frames[start_k:]
    else:
        seg_idx_seq = [0] * (M - start_k)

    for kk, k in enumerate(range(M)):
        if k < start_k:
            continue
        tfrac = 0.0 if M == 1 else (k / (M - 1.0))

        C = (1.0 - tfrac) * coordsA_aligned + tfrac * coordsB_aligned

        idx_sel = [j for j in range(len(match_tpl_idx)) if match_tpl_idx[j] >= 0]
        if len(idx_sel) >= 3:
            Y = np.array([C[match_tpl_idx[j]] for j in idx_sel], dtype=float)
            P_bohr = np.array(pair_images[k].coords3d, dtype=float)
            P = P_bohr * BOHR2ANG
            P_sel = np.array([P[j] for j in idx_sel], dtype=float)
            R, t = kabsch_R_t(Y, P_sel)
            Paligned = (P @ R) + t
            for jj, pidx in enumerate(idx_sel):
                full_i = match_tpl_idx[pidx]
                if 0 <= full_i < Nfull:
                    C[full_i] = Paligned[pidx]
        else:
            P = np.array(pair_images[k].coords3d, dtype=float) * BOHR2ANG
            for j, full_i in enumerate(match_tpl_idx):
                if full_i >= 0 and 0 <= full_i < Nfull and j < P.shape[0]:
                    C[full_i] = P[j]

        for i, a in enumerate(atomsA):
            a.set_coord(C[i])
            a.set_bfactor(100.0 if i in active_full_idx else 0.0)

        remark_lines: List[str] = []
        remark_lines.append(f"PAIR_MERGE FRAC {tfrac:.6f}")

        if include_pocket_indices_for_first_model and kk == 0:
            remark_lines.extend(_chunk_remark_indices([i for i in active_one_based], width=60))

        seg_idx = seg_idx_seq[kk] if kk < len(seg_idx_seq) else 0
        if seg_idx and seg_report_lookup is not None:
            rep = seg_report_lookup.get(seg_idx, None)
            if rep is not None:
                remark_lines.append(
                    f"MEP_SEG_INDEX {int(seg_idx):02d} TAG {rep.tag} KIND {rep.kind} "
                    f"DELTAE_BARRIER_KCAL {rep.barrier_kcal:.6f} DELTAE_KCAL {rep.delta_kcal:.6f}"
                )
                if rep.kind != "bridge" and rep.summary and rep.summary.strip() and rep.summary.strip() != "(no covalent changes detected)":
                    for ln in rep.summary.strip().splitlines():
                        remark_lines.append(f"SEG_BONDS {ln.strip()}")
            else:
                remark_lines.append(f"MEP_SEG_INDEX {int(seg_idx):02d}")

        model_blocks.append(_write_model_block(structA, remark_lines))

    if out_path is not None:
        with open(out_path, "w") as f:
            for m, blk in enumerate(model_blocks, start=1):
                f.write(f"MODEL     {m}\n")
                f.write(blk)
                f.write("ENDMDL\n")
            f.write("END\n")
        click.echo(f"[merge] Wrote pair-merged PDB → '{out_path}'")

    return model_blocks, active_one_based


def _merge_final_and_write(final_images: List[Any],
                           pocket_inputs: Sequence[Path],
                           ref_pdbs: Sequence[Path],
                           segments: List[SegmentReport],
                           out_dir: Path) -> None:
    """
    Merge the entire pocket MEP into full templates (for all pairs) and write outputs.
    """
    # NOTE: In --align True mode, caller may pass a single ref PDB replicated per pair.
    if len(ref_pdbs) != len(pocket_inputs):
        raise click.BadParameter("--ref-pdb must match the number of --input after preprocessing (caller should replicate the first ref for all pairs when --align True).")

    structs, aligned_coords, atoms_list, keymaps = _load_structures_and_chain_align(ref_pdbs)

    seg_lookup: Dict[int, SegmentReport] = {int(s.seg_index): s for s in segments if int(s.seg_index) > 0}

    n_pairs = len(pocket_inputs) - 1
    groups: List[Tuple[int, List[Any]]] = []
    cur_idx = None
    cur_list: List[Any] = []
    for im in final_images:
        pi = getattr(im, "pair_index", None)
        if pi is None:
            pi = 0
        if cur_idx is None:
            cur_idx = int(pi)
            cur_list = [im]
        elif int(pi) == int(cur_idx):
            cur_list.append(im)
        else:
            groups.append((int(cur_idx), cur_list))
            cur_idx = int(pi)
            cur_list = [im]
    if cur_list:
        groups.append((int(cur_idx), cur_list))

    for (pi, _) in groups:
        if not (0 <= pi < n_pairs):
            raise click.BadParameter(f"[merge] Illegal pair_index {pi} (n_pairs={n_pairs}).")

    final_blocks: List[str] = []
    wrote_indices = False

    for gi, (pi, imgs) in enumerate(groups):
        pocket_ref = Path(pocket_inputs[pi])
        structA = structs[pi]
        structB = structs[pi+1]
        coordsA = aligned_coords[pi]
        coordsB = aligned_coords[pi+1]
        keymapA = keymaps[pi]
        keymapB = keymaps[pi+1]
        seg_indices_for_frames = [int(getattr(im, "mep_seg_index", 0) or 0) for im in imgs]

        blocks, active_one_based = _merge_pair_to_full(
            pair_images=imgs,
            pocket_ref_pdb=pocket_ref,
            structA=structA,
            structB=structB,
            coordsA_aligned=coordsA,
            coordsB_aligned=coordsB,
            keymapA=keymapA,
            keymapB=keymapB,
            out_path=None,
            drop_first=(gi > 0),
            seg_indices_for_frames=seg_indices_for_frames,
            seg_report_lookup=seg_lookup,
            include_pocket_indices_for_first_model=(not wrote_indices)
        )
        if blocks:
            wrote_indices = True
            final_blocks.extend(blocks)

    final_path = out_dir / "mep_w_ref.pdb"
    with open(final_path, "w") as f:
        for m, blk in enumerate(final_blocks, start=1):
            f.write(f"MODEL     {m}\n")
            f.write(blk)
            f.write("ENDMDL\n")
        f.write("END\n")
    click.echo(f"[merge] Wrote concatenated full-system trajectory → '{final_path}'")

    # Per-segment merged MEPs (bond-change segments only) + HEI merged only for bond-change segments
    for s in segments:
        seg_idx = int(s.seg_index)
        seg_frames: List[Any] = [im for im in final_images if int(getattr(im, "mep_seg_index", 0) or 0) == seg_idx]
        if not seg_frames:
            continue

        # Determine pair index for this segment (assume consistent within the segment)
        pi_vals = sorted({int(getattr(im, "pair_index", 0)) for im in seg_frames})
        pi = pi_vals[0]
        pocket_ref = Path(pocket_inputs[pi])
        structA = structs[pi]
        structB = structs[pi+1]
        coordsA = aligned_coords[pi]
        coordsB = aligned_coords[pi+1]
        keymapA = keymaps[pi]
        keymapB = keymaps[pi+1]

        # Per-segment merged MEP only when covalent changes are present
        if s.kind != "bridge" and s.summary and s.summary.strip() != "(no covalent changes detected)":
            seg_indices_for_frames = [seg_idx] * len(seg_frames)
            blocks, _ = _merge_pair_to_full(
                pair_images=seg_frames,
                pocket_ref_pdb=pocket_ref,
                structA=structA,
                structB=structB,
                coordsA_aligned=coordsA,
                coordsB_aligned=coordsB,
                keymapA=keymapA,
                keymapB=keymapB,
                out_path=None,
                drop_first=False,
                seg_indices_for_frames=seg_indices_for_frames,
                seg_report_lookup=seg_lookup,
                include_pocket_indices_for_first_model=True,
            )
            out_seg = out_dir / f"mep_w_ref_seg_{seg_idx:02d}.pdb"
            with open(out_seg, "w") as f:
                for m, blk in enumerate(blocks, start=1):
                    f.write(f"MODEL     {m}\n")
                    f.write(blk)
                    f.write("ENDMDL\n")
                f.write("END\n")
            click.echo(f"[merge] Wrote per-segment merged trajectory → '{out_seg}'")

        # Per-segment HEI merged to reference (only for bond-change segments)
        if s.kind != "bridge" and s.summary and s.summary.strip() != "(no covalent changes detected)":
            try:
                energies_eh = [float(getattr(im, "energy")) for im in seg_frames]
                imax = int(np.argmax(np.array(energies_eh, dtype=float)))
                hei_frame = seg_frames[imax]
                blocks_hei, _ = _merge_pair_to_full(
                    pair_images=[hei_frame],
                    pocket_ref_pdb=pocket_ref,
                    structA=structA,
                    structB=structB,
                    coordsA_aligned=coordsA,
                    coordsB_aligned=coordsB,
                    keymapA=keymapA,
                    keymapB=keymapB,
                    out_path=None,
                    drop_first=False,
                    seg_indices_for_frames=[seg_idx],
                    seg_report_lookup=seg_lookup,
                    include_pocket_indices_for_first_model=True,
                )
                out_hei = out_dir / f"hei_w_ref_seg_{seg_idx:02d}.pdb"
                with open(out_hei, "w") as f:
                    for m, blk in enumerate(blocks_hei, start=1):
                        f.write(f"MODEL     {m}\n")
                        f.write(blk)
                        f.write("ENDMDL\n")
                    f.write("END\n")
                click.echo(f"[merge] Wrote merged HEI for segment → '{out_hei}'")
            except Exception as e:
                click.echo(f"[merge] WARNING: Failed to write merged HEI for segment {seg_idx:02d}: {e}", err=True)


# -----------------------------------------------
# CLI
# -----------------------------------------------

@click.command(
    help="Multistep MEP search via recursive GSM segmentation.",
    context_settings={
        "help_option_names": ["-h", "--help"],
        # Allow a single '-i' followed by multiple paths (as extra args)
        "ignore_unknown_options": True,
        "allow_extra_args": True,
    },
)
@click.option(
    "-i", "--input",
    "input_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True,   # allow: -i A -i B -i C   or   -i A B C
    required=True,
    help=("Two or more structures in reaction order. "
          "Either repeat '-i' (e.g., '-i A -i B -i C') or use a single '-i' "
          "followed by multiple space-separated paths (e.g., '-i A B C').")
)
@click.option(
    "--parm",
    "real_parm7",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Amber parm7 topology covering the full enzyme complex.",
)
@click.option(
    "--model-pdb",
    "model_pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help="PDB describing atoms that belong to the ML (high-level) region. "
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
@click.option(
    "-q",
    "--charge",
    type=int,
    required=True,
    help="Total system charge.",
)
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=None,
    show_default=False,
    help="Spin multiplicity (2S+1). Defaults to 1 when omitted.",
)
@click.option(
    "--mep-mode",
    "mep_mode",
    type=click.Choice(["gsm", "dmf"], case_sensitive=False),
    default="gsm",
    show_default=True,
    help="MEP method: gsm (GrowingString) or dmf (Direct Max Flux).",
)
@click.option(
    "--refine-mode",
    type=click.Choice(["peak", "minima"], case_sensitive=False),
    default=None,
    show_default=True,
    help=(
        "Refinement seed around the highest-energy image: "
        "'peak' uses HEI±1, 'minima' uses nearest local minima. "
        "Defaults to peak for gsm and minima for dmf."
    ),
)
@click.option(
    "--freeze-atoms",
    "freeze_atoms_text",
    type=str,
    default=None,
    show_default=False,
    help="Comma-separated 1-based atom indices to freeze (e.g., '1,3,5').",
)
@click.option(
    "--hess-cutoff",
    "hess_cutoff",
    type=float,
    default=None,
    show_default=False,
    help="Distance cutoff (Å) from ML region for MM atoms to include in Hessian calculation. "
         "Applied to movable MM atoms and can be combined with --detect-layer.",
)
@click.option(
    "--movable-cutoff",
    "movable_cutoff",
    type=float,
    default=None,
    show_default=False,
    help="Distance cutoff (Å) from ML region for movable MM atoms. MM atoms beyond this are frozen. "
         "Providing --movable-cutoff disables --detect-layer.",
)
@click.option("--max-nodes", type=int, default=10, show_default=True,
              help=("Number of internal nodes (string has max_nodes+2 images including endpoints). "
                    "Used for *segment* GSM unless overridden by YAML search.max_nodes_segment."))
@click.option("--max-cycles", type=int, default=300, show_default=True, help="Maximum GSM optimization cycles.")
@click.option(
    "--climb/--no-climb",
    default=True,
    show_default=True,
    help="Enable transition-state search after path growth.",
)
@click.option(
    "--dump/--no-dump",
    default=False,
    show_default=True,
    help="Dump GSM/single-optimization trajectories during the run.",
)
@click.option(
    "--opt-mode",
    "opt_mode",
    type=click.Choice(["light", "heavy"], case_sensitive=False),
    default="light",
    show_default=True,
    help="Single-structure optimizer: light (=LBFGS) or heavy (=RFO).",
)
@click.option("--out-dir", "out_dir", type=str, default="./result_path_search/", show_default=True, help="Output directory.")
@click.option(
    "--thresh",
    type=str,
    default=None,
    help="Convergence preset for GSM/StringOptimizer and single LBFGS runs (gau_loose|gau|gau_tight|gau_vtight|baker|never).",
)
@click.option(
    "--config",
    "config_yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Base YAML configuration file applied before explicit CLI options.",
)
@click.option(
    "--show-config/--no-show-config",
    "show_config",
    default=False,
    show_default=True,
    help="Print resolved configuration and continue execution.",
)
@click.option(
    "--dry-run/--no-dry-run",
    "dry_run",
    default=False,
    show_default=True,
    help="Validate options and print the execution plan without running path search.",
)
@click.option(
    "--preopt/--no-preopt",
    "pre_opt",
    default=True,
    show_default=True,
    help="If False, skip initial single-structure optimizations of inputs."
)
# Input alignment switch (default True)
@click.option(
    "--align/--no-align",
    "align",
    default=True,
    show_default=True,
    help=("After pre-optimization, align all inputs to the *first* input and match freeze_atoms "
          "using the align_freeze_atoms API. When --align is True and --ref-pdb is provided, "
          "the first reference PDB will be used for all pairs in the final merge.")
)
# Full template PDBs for merge
@click.option(
    "--ref-pdb",
    "ref_pdb_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True,
    default=None,
    help=("Full-size template PDBs in the same reaction order as --input. "
          "With --align True, only the *first* provided reference PDB is used for all pairs "
          "in the final merge (you may pass just one).")
)
@click.pass_context
def cli(
    ctx: click.Context,
    input_paths: Sequence[Path],
    real_parm7: Path,
    model_pdb: Optional[Path],
    model_indices_str: Optional[str],
    model_indices_one_based: bool,
    detect_layer: bool,
    charge: Optional[int],
    spin: Optional[int],
    mep_mode: str,
    refine_mode: Optional[str],
    freeze_atoms_text: Optional[str],
    hess_cutoff: Optional[float],
    movable_cutoff: Optional[float],
    max_nodes: int,
    max_cycles: int,
    climb: bool,
    dump: bool,
    opt_mode: str,
    out_dir: str,
    thresh: Optional[str],
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    pre_opt: bool,
    align: bool,
    ref_pdb_paths: Optional[Sequence[Path]],
) -> None:
    prepared_inputs: List[PreparedInputStructure] = []
    # --- Robustly accept both styles for -i/--input and --ref-pdb ---
    argv_all = sys.argv[1:]  # drop program name
    i_vals = collect_single_option_values(argv_all, ("-i", "--input"), label="-i/--input")
    if i_vals:
        i_parsed = validate_existing_files(
            i_vals,
            option_name="-i/--input",
            hint="When using '-i', list only existing file paths (multiple paths may follow a single '-i').",
        )
        input_paths = tuple(i_parsed)

    ref_vals = collect_single_option_values(argv_all, ("--ref-pdb",), label="--ref-pdb")
    if ref_vals:
        ref_parsed = validate_existing_files(
            ref_vals,
            option_name="--ref-pdb",
            hint="When using '--ref-pdb', multiple files may follow a single option.",
        )
        ref_pdb_paths = tuple(ref_parsed)
    # --- end of robust parsing fix ---

    def _is_param_explicit(name: str) -> bool:
        try:
            source = ctx.get_parameter_source(name)
            return source not in (None, ParameterSource.DEFAULT)
        except Exception:
            return False

    config_yaml, override_yaml, used_legacy_yaml = resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )
    merged_yaml_cfg, _, _ = load_merged_yaml_cfg(
        config_yaml=config_yaml,
        override_yaml=None,
    )

    time_start = time.perf_counter()  # start timing
    command_str = "mlmm path-search " + " ".join(sys.argv[1:])
    try:
        # --------------------------
        # 0) Input validation (multi-structure)
        # --------------------------
        if len(input_paths) < 2:
            raise click.BadParameter("Provide at least two structures for --input in reaction order (reactant [intermediates ...] product).")

        do_merge = bool(ref_pdb_paths) and len(ref_pdb_paths) > 0
        if do_merge:
            if align:
                # With --align True we will use only the *first* ref PDB for all pairs (replicated later).
                pass
            else:
                if len(ref_pdb_paths) != len(input_paths):
                    raise click.BadParameter("--ref-pdb must be given for each --input (same count and order). "
                                             "Alternatively, use --align to allow using only the first reference PDB for all pairs.")

        p_list = [Path(p) for p in input_paths]
        for p in p_list:
            if p.suffix.lower() != ".pdb":
                raise click.BadParameter(
                    f"'{p}' is not a PDB file. ML/MM path-search currently requires PDB inputs for all states."
                )
        prepared_inputs = [prepare_input_structure(p) for p in p_list]
        # --------------------------
        # 1) Resolve settings (defaults < config < CLI(explicit) < override)
        # --------------------------
        config_layer_cfg = load_yaml_dict(config_yaml)
        override_layer_cfg = load_yaml_dict(override_yaml)

        mep_mode_kind = mep_mode.lower().strip()
        refine_mode_kind = refine_mode.strip().lower() if refine_mode else None

        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        gs_cfg = dict(GS_KW)
        opt_cfg = dict(STOPT_KW)
        lbfgs_cfg = dict(LBFGS_KW)
        bond_cfg = dict(BOND_KW)
        search_cfg = dict(SEARCH_KW)
        dmf_cfg = dict(DMF_KW)

        apply_yaml_overrides(
            config_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",), ("mlmm",))),
                (gs_cfg, (("gs",),)),
                (opt_cfg, (("opt",),)),
                (lbfgs_cfg, (("sopt", "lbfgs"), ("lbfgs",))),
                (bond_cfg, (("bond",),)),
                (search_cfg, (("search",),)),
                (dmf_cfg, (("dmf",),)),
            ],
        )

        try:
            geom_freeze = _normalize_geom_freeze(geom_cfg.get("freeze_atoms"))
        except click.BadParameter as e:
            click.echo(f"ERROR: {e}", err=True)
            sys.exit(1)
        geom_cfg["freeze_atoms"] = geom_freeze

        try:
            cli_freeze = _parse_freeze_atoms(freeze_atoms_text)
        except click.BadParameter as e:
            click.echo(f"ERROR: {e}", err=True)
            sys.exit(1)

        model_indices: Optional[List[int]] = None
        if model_indices_str:
            try:
                model_indices = parse_indices_string(model_indices_str, one_based=model_indices_one_based)
            except click.BadParameter as e:
                click.echo(f"ERROR: {e}", err=True)
                sys.exit(1)
        if cli_freeze:
            merge_freeze_atom_indices(geom_cfg, cli_freeze)

        freeze_atoms_final = list(geom_cfg.get("freeze_atoms") or [])
        calc_cfg["freeze_atoms"] = freeze_atoms_final

        resolved_charge = charge
        resolved_spin = spin
        for prepared in prepared_inputs:
            resolved_charge, resolved_spin = resolve_charge_spin_or_raise(
                prepared,
                resolved_charge,
                resolved_spin,
            )
        charge_value = calc_cfg.get("model_charge", resolved_charge)
        if charge_value is None:
            charge_value = resolved_charge
        calc_cfg["model_charge"] = int(charge_value)
        if _is_param_explicit("charge"):
            calc_cfg["model_charge"] = int(resolved_charge)

        spin_value = calc_cfg.get("model_mult", resolved_spin)
        if spin_value is None:
            spin_value = resolved_spin
        calc_cfg["model_mult"] = int(spin_value)
        if _is_param_explicit("spin"):
            calc_cfg["model_mult"] = int(resolved_spin)

        first_input = p_list[0]
        calc_cfg["input_pdb"] = str(first_input)
        calc_cfg["real_parm7"] = str(real_parm7)

        detect_layer_effective = bool(calc_cfg.get("use_bfactor_layers", detect_layer))
        if _is_param_explicit("detect_layer"):
            detect_layer_effective = bool(detect_layer)

        if _is_param_explicit("max_nodes"):
            gs_cfg["max_nodes"] = int(max_nodes)
            search_cfg["max_nodes_segment"] = int(max_nodes)
        if _is_param_explicit("max_cycles"):
            opt_cfg["max_cycles"] = int(max_cycles)
            opt_cfg["stop_in_when_full"] = int(max_cycles)
            dmf_cfg["max_cycles"] = int(max_cycles)
        if _is_param_explicit("climb"):
            gs_cfg["climb"] = bool(climb)
            gs_cfg["climb_lanczos"] = bool(climb)
        if _is_param_explicit("dump"):
            opt_cfg["dump"] = bool(dump)
            lbfgs_cfg["dump"] = bool(dump)
        if _is_param_explicit("out_dir"):
            opt_cfg["out_dir"] = out_dir
            lbfgs_cfg["out_dir"] = out_dir
        if _is_param_explicit("thresh") and thresh is not None:
            opt_cfg["thresh"] = str(thresh)
            lbfgs_cfg["thresh"] = str(thresh)
        if _is_param_explicit("hess_cutoff") and hess_cutoff is not None:
            calc_cfg["hess_cutoff"] = float(hess_cutoff)
        if _is_param_explicit("movable_cutoff") and movable_cutoff is not None:
            calc_cfg["movable_cutoff"] = float(movable_cutoff)
            detect_layer_effective = False
        if _is_param_explicit("refine_mode"):
            search_cfg["refine_mode"] = refine_mode_kind

        apply_yaml_overrides(
            override_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",), ("mlmm",))),
                (gs_cfg, (("gs",),)),
                (opt_cfg, (("opt",),)),
                (lbfgs_cfg, (("sopt", "lbfgs"), ("lbfgs",))),
                (bond_cfg, (("bond",),)),
                (search_cfg, (("search",),)),
                (dmf_cfg, (("dmf",),)),
            ],
        )

        refine_mode_kind = search_cfg.get("refine_mode")
        if refine_mode_kind is None:
            refine_mode_kind = "peak" if mep_mode_kind == "gsm" else "minima"
        else:
            refine_mode_kind = str(refine_mode_kind).strip().lower()
            if refine_mode_kind not in {"peak", "minima"}:
                raise click.BadParameter(f"Unknown --refine-mode '{refine_mode_kind}'.")
        search_cfg["refine_mode"] = refine_mode_kind

        out_dir_path = Path(opt_cfg.get("out_dir", out_dir)).resolve()
        detect_layer_effective = bool(calc_cfg.get("use_bfactor_layers", detect_layer_effective))

        model_pdb_effective: Optional[Path] = None
        if _is_param_explicit("model_pdb") and model_pdb is not None:
            model_pdb_effective = Path(model_pdb)
        else:
            model_pdb_cfg = calc_cfg.get("model_pdb")
            if isinstance(model_pdb_cfg, (str, Path)) and str(model_pdb_cfg).strip():
                model_pdb_effective = Path(model_pdb_cfg)

        hess_cutoff_effective = calc_cfg.get("hess_cutoff")
        movable_cutoff_effective = calc_cfg.get("movable_cutoff")
        if movable_cutoff_effective is not None:
            if detect_layer_effective:
                click.echo("[layer] movable_cutoff is set; disabling detect-layer.", err=True)
            detect_layer_effective = False

        layer_source_pdb = first_input
        if detect_layer_effective and layer_source_pdb.suffix.lower() != ".pdb":
            click.echo("ERROR: --detect-layer requires a PDB input.", err=True)
            sys.exit(1)

        if dry_run:
            layer_info_preview: Optional[Dict[str, List[int]]] = None
            model_region_source = "bfactor"

            if detect_layer_effective:
                try:
                    bfactors = read_bfactors_from_pdb(layer_source_pdb)
                    if not bfactors:
                        raise ValueError(f"No ATOM/HETATM records found in {layer_source_pdb}.")
                    if not has_valid_layer_bfactors(bfactors):
                        raise ValueError(
                            "Invalid or missing layer B-factors (expected ~0/10/20). "
                            "Provide --no-detect-layer with --model-pdb/--model-indices."
                        )
                    layer_info_preview = parse_layer_indices_from_bfactors(bfactors)
                    if not layer_info_preview.get("ml_indices"):
                        raise ValueError("No ML atoms detected from B-factors (value ~0).")
                except Exception as e:
                    if model_pdb_effective is None and not model_indices:
                        click.echo(f"ERROR: {e}", err=True)
                        sys.exit(1)
                    click.echo(f"[layer] WARNING: {e} Falling back to explicit ML region.", err=True)
                    detect_layer_effective = False

            if not detect_layer_effective:
                if model_pdb_effective is not None:
                    model_region_source = "model_pdb"
                elif model_indices:
                    model_region_source = "model_indices"
                    if layer_source_pdb.suffix.lower() != ".pdb":
                        click.echo("ERROR: --model-indices requires a PDB input.", err=True)
                        sys.exit(1)
                    n_atoms = 0
                    with layer_source_pdb.open("r", encoding="utf-8", errors="ignore") as fh:
                        for line in fh:
                            if line.startswith(("ATOM  ", "HETATM")):
                                n_atoms += 1
                    bad_idx = [i for i in model_indices if i < 0 or i >= n_atoms]
                    if bad_idx:
                        click.echo(
                            f"ERROR: model index out of range: {bad_idx[0]} (valid: 0 <= idx < {n_atoms})",
                            err=True,
                        )
                        sys.exit(1)
                else:
                    click.echo("ERROR: Provide --model-pdb or --model-indices when --no-detect-layer.", err=True)
                    sys.exit(1)

            if show_config:
                click.echo(
                    pretty_block(
                        "yaml_layers",
                        {
                            "config": None if config_yaml is None else str(config_yaml),
                            "override": None if override_yaml is None else str(override_yaml),
                            "merged_keys": sorted(merged_yaml_cfg.keys()),
                        },
                    )
                )

            dry_payload: Dict[str, Any] = {
                "input_count": len(p_list),
                "input_first": str(p_list[0]) if p_list else None,
                "input_last": str(p_list[-1]) if p_list else None,
                "output_dir": str(out_dir_path),
                "mep_mode": mep_mode_kind,
                "refine_mode": refine_mode_kind,
                "opt_mode": str(opt_mode),
                "detect_layer": bool(detect_layer_effective),
                "model_region_source": model_region_source,
                "model_indices_count": 0 if not model_indices else len(model_indices),
                "pre_opt": bool(pre_opt),
                "align": bool(align),
                "max_depth": int(search_cfg.get("max_depth", SEARCH_KW["max_depth"])),
                "max_nodes_segment": int(search_cfg.get("max_nodes_segment", gs_cfg.get("max_nodes", 0))),
                "will_run_path_search": True,
                "will_write_summary": True,
            }
            if layer_info_preview is not None:
                dry_payload["layer_counts"] = {
                    "ml": len(layer_info_preview.get("ml_indices", [])),
                    "movable_mm": len(layer_info_preview.get("movable_mm_indices", [])),
                    "frozen": len(layer_info_preview.get("frozen_indices", [])),
                    "unassigned": len(layer_info_preview.get("unassigned_indices", [])),
                }

            click.echo(pretty_block("dry_run_plan", dry_payload))
            click.echo("[dry-run] Validation complete. Path search execution was skipped.")
            return

        model_pdb_path: Optional[Path] = None
        layer_info: Optional[Dict[str, List[int]]] = None

        if detect_layer_effective:
            try:
                model_pdb_path, layer_info = build_model_pdb_from_bfactors(layer_source_pdb, out_dir_path)
                calc_cfg["use_bfactor_layers"] = True
                click.echo(
                    f"[layer] Detected B-factor layers: ML={len(layer_info.get('ml_indices', []))}, "
                    f"MovableMM={len(layer_info.get('movable_mm_indices', []))}, "
                    f"FrozenMM={len(layer_info.get('frozen_indices', []))}"
                )
            except Exception as e:
                if model_pdb_effective is None and not model_indices:
                    click.echo(f"ERROR: {e}", err=True)
                    sys.exit(1)
                click.echo(f"[layer] WARNING: {e} Falling back to explicit ML region.", err=True)
                detect_layer_effective = False

        if not detect_layer_effective:
            if model_pdb_effective is None and not model_indices:
                click.echo("ERROR: Provide --model-pdb or --model-indices when --no-detect-layer.", err=True)
                sys.exit(1)
            if model_pdb_effective is not None:
                model_pdb_path = Path(model_pdb_effective)
            else:
                if layer_source_pdb.suffix.lower() != ".pdb":
                    click.echo("ERROR: --model-indices requires a PDB input.", err=True)
                    sys.exit(1)
                try:
                    model_pdb_path = build_model_pdb_from_indices(layer_source_pdb, out_dir_path, model_indices or [])
                except Exception as e:
                    click.echo(f"ERROR: {e}", err=True)
                    sys.exit(1)
            calc_cfg["use_bfactor_layers"] = False

        if model_pdb_path is None:
            click.echo("ERROR: Failed to resolve model PDB for the ML region.", err=True)
            sys.exit(1)

        calc_cfg["model_pdb"] = str(model_pdb_path)
        freeze_atoms_final = apply_layer_freeze_constraints(
            geom_cfg,
            calc_cfg,
            layer_info,
            echo_fn=click.echo,
        )

        # Distance-based overrides for Hessian-target and movable MM selection.
        if hess_cutoff_effective is not None:
            calc_cfg["hess_cutoff"] = float(hess_cutoff_effective)
        if movable_cutoff_effective is not None:
            calc_cfg["movable_cutoff"] = float(movable_cutoff_effective)
            calc_cfg["use_bfactor_layers"] = False

        for key in ("input_pdb", "real_parm7", "model_pdb", "mm_fd_dir"):
            val = calc_cfg.get(key)
            if isinstance(val, (str, Path)):
                calc_cfg[key] = str(Path(val).expanduser().resolve())

        opt_cfg["stop_in_when_full"] = int(opt_cfg.get("max_cycles", STOPT_KW["max_cycles"]))
        out_dir_path = Path(opt_cfg.get("out_dir", out_dir)).resolve()
        echo_geom = format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")
        echo_calc = strip_inherited_keys(calc_cfg, CALC_KW, mode="same")
        echo_calc = format_freeze_atoms_for_echo(echo_calc, key="freeze_atoms")
        echo_gs   = strip_inherited_keys(gs_cfg, GS_KW, mode="same")
        echo_opt  = strip_inherited_keys({**opt_cfg, "out_dir": str(out_dir_path)}, STOPT_KW, mode="same")
        echo_lbfgs = strip_inherited_keys(lbfgs_cfg, LBFGS_KW, mode="same")
        echo_bond = strip_inherited_keys(bond_cfg, BOND_KW, mode="same")
        echo_search = strip_inherited_keys(search_cfg, SEARCH_KW, mode="same")

        click.echo(pretty_block("geom", echo_geom))
        click.echo(pretty_block("calc", echo_calc))
        click.echo(pretty_block("gs",   echo_gs))
        click.echo(pretty_block("opt",  echo_opt))
        click.echo(pretty_block("lbfgs", echo_lbfgs))
        click.echo(pretty_block("bond", echo_bond))
        click.echo(pretty_block("search", echo_search))
        # Echo pre-optimization and alignment flags
        click.echo(
            pretty_block(
                "run_flags",
                {
                    "pre_opt": bool(pre_opt),
                    "align": bool(align),
                    "mep_mode": mep_mode_kind,
                    "refine_mode": refine_mode_kind,
                    "opt_mode": str(opt_mode),
                },
            )
        )

        if show_config:
            click.echo(
                pretty_block(
                    "yaml_layers",
                    {
                        "config": None if config_yaml is None else str(config_yaml),
                        "override": None if override_yaml is None else str(override_yaml),
                        "merged_keys": sorted(merged_yaml_cfg.keys()),
                    },
                )
            )

        if int(opt_cfg.get("max_cycles", 0)) <= 0:
            click.echo("[INFO] max_cycles <= 0: skipping path search.")
            return

        # --------------------------
        # 2) Prepare inputs
        # --------------------------
        out_dir_path.mkdir(parents=True, exist_ok=True)

        geoms = _load_structures(
            inputs=prepared_inputs,
            coord_type=geom_cfg.get("coord_type", "cart"),
            base_freeze=geom_cfg.get("freeze_atoms", []),
        )

        shared_calc = mlmm(**calc_cfg)
        for g in geoms:
            g.set_calculator(shared_calc)

        # If any input is PDB, treat as "PDB input" for final output handling.
        ref_pdb_for_segments: Optional[Path] = None
        for p in p_list:
            if p.suffix.lower() == ".pdb":
                ref_pdb_for_segments = p.resolve()
                break

        if pre_opt:
            new_geoms: List[Any] = []
            for i, g in enumerate(geoms):
                tag = f"init{i:02d}"
                g_opt = _optimize_single(g, shared_calc, lbfgs_cfg, out_dir_path, tag=tag, ref_pdb_path=ref_pdb_for_segments)
                new_geoms.append(g_opt)
            geoms = new_geoms
        else:
            click.echo("[init] Skipping endpoint pre-optimization as requested by --no-preopt.")

        # Align all inputs to the first structure, guided by freeze constraints, when requested
        align_thresh = str(opt_cfg.get("thresh", "gau"))
        if align:
            try:
                click.echo("\n=== Aligning all inputs to the first structure (freeze-guided scan + relaxation) ===\n")
                _ = align_and_refine_sequence_inplace(
                    geoms,
                    thresh=align_thresh,
                    shared_calc=shared_calc,
                    out_dir=out_dir_path / "align_refine",
                    verbose=True,
                )
                click.echo("[align] Completed input alignment.")
            except Exception as e:
                click.echo(f"[align] WARNING: Alignment failed; continuing without alignment: {e}", err=True)
        else:
            click.echo("[align] Skipping input alignment as requested by --no-align.")

        # --------------------------
        # 3) Run recursive search for each adjacent pair and stitch
        # --------------------------
        click.echo("\n=== Multistep MEP search (multi-structure) started ===\n")
        seg_counter = [0]

        bridge_max_nodes = int(search_cfg.get("max_nodes_bridge", 5))
        gs_bridge_cfg = {**gs_cfg, "max_nodes": bridge_max_nodes, "climb": False, "climb_lanczos": False}

        combined_imgs: List[Any] = []
        combined_Es: List[float] = []
        seg_reports_all: List[SegmentReport] = []

        def _segment_builder_for_pairs(tail_g, head_g, _tag: str) -> CombinedPath:
            sub = _build_multistep_path(
                tail_g, head_g,
                shared_calc,
                geom_cfg, gs_cfg, opt_cfg,
                lbfgs_cfg,
                bond_cfg, search_cfg, refine_mode_kind,
                out_dir=out_dir_path,
                ref_pdb_path=ref_pdb_for_segments,
                depth=0,
                seg_counter=seg_counter,
                branch_tag="B",
                pair_index=None,
                mep_mode_kind=mep_mode_kind, calc_cfg=calc_cfg, dmf_cfg=dmf_cfg,
            )
            return sub

        for i in range(len(geoms) - 1):
            gA, gB = geoms[i], geoms[i + 1]
            pair_tag = f"pair_{i:02d}"
            click.echo(f"\n--- Processing pair {i:02d}: image {i} → {i+1} ---")
            pair_path = _build_multistep_path(
                gA, gB,
                shared_calc,
                geom_cfg, gs_cfg, opt_cfg,
                lbfgs_cfg,
                bond_cfg, search_cfg, refine_mode_kind,
                out_dir=out_dir_path,
                ref_pdb_path=ref_pdb_for_segments,
                depth=0,
                seg_counter=seg_counter,
                branch_tag=pair_tag,
                pair_index=i,
                mep_mode_kind=mep_mode_kind, calc_cfg=calc_cfg, dmf_cfg=dmf_cfg,
            )

            if i == 0:
                combined_imgs = list(pair_path.images)
                combined_Es = list(pair_path.energies)
                seg_reports_all.extend(pair_path.segments)
            else:
                parts = [(combined_imgs, combined_Es), (pair_path.images, pair_path.energies)]
                combined_imgs, combined_Es = _stitch_paths(
                    parts=parts,
                    stitch_rmsd_thresh=float(search_cfg["stitch_rmsd_thresh"]),
                    bridge_rmsd_thresh=float(search_cfg["bridge_rmsd_thresh"]),
                    shared_calc=shared_calc,
                    gs_cfg=gs_bridge_cfg,
                    opt_cfg=opt_cfg,
                    out_dir=out_dir_path,
                    tag=pair_tag,
                    ref_pdb_path=ref_pdb_for_segments,
                    bond_cfg=bond_cfg,
                    segment_builder=_segment_builder_for_pairs,
                    segments_out=seg_reports_all,
                    bridge_pair_index=i,
                    mep_mode_kind=mep_mode_kind, calc_cfg=calc_cfg, dmf_cfg=dmf_cfg,
                )
                seg_reports_all.extend(pair_path.segments)

        click.echo("\n=== Multistep MEP search (multi-structure) finished ===\n")

        combined_all = CombinedPath(images=combined_imgs, energies=combined_Es, segments=seg_reports_all)

        # --------------------------
        # 4) Outputs
        # --------------------------
        for idx, srep in enumerate(combined_all.segments, 1):
            srep.seg_index = idx
        tag_to_index = {s.tag: int(s.seg_index) for s in combined_all.segments}
        for im in combined_all.images:
            tag = getattr(im, "mep_seg_tag", None)
            if tag and tag in tag_to_index:
                try:
                    setattr(im, "mep_seg_index", int(tag_to_index[tag]))
                except Exception:
                    pass

        # Always write mep_trj.xyz for downstream compatibility; convert to PDB when possible.
        pdb_input = ref_pdb_for_segments is not None
        final_trj = out_dir_path / "mep_trj.xyz"
        _write_xyz_trj_with_energy(combined_all.images, combined_all.energies, final_trj)
        click.echo(f"[write] Wrote '{final_trj}'.")
        try:
            run_trj2fig(final_trj, [out_dir_path / "mep_plot.png"], unit="kcal", reference="init", reverse_x=False)
            click.echo(f"[plot] Saved energy plot → '{out_dir_path / 'mep_plot.png'}'")
        except Exception as e:
            click.echo(f"[plot] WARNING: Failed to plot final energy: {e}", err=True)

        if pdb_input:
            try:
                final_pdb = out_dir_path / "mep.pdb"
                convert_xyz_to_pdb(final_trj, ref_pdb_for_segments, final_pdb)
                click.echo(f"[convert] Wrote '{final_pdb}'.")
            except Exception as e:
                click.echo(f"[convert] WARNING: Failed to convert final MEP to PDB: {e}", err=True)

        # ---- Pocket-only per-segment trajectories & HEIs ----
        try:
            # Map frames → segment indices
            frame_seg_indices: List[int] = [int(getattr(im, "mep_seg_index", 0) or 0) for im in combined_all.images]
            seg_to_frames: Dict[int, List[int]] = {}
            for ii, sidx in enumerate(frame_seg_indices):
                if sidx <= 0:
                    continue
                seg_to_frames.setdefault(int(sidx), []).append(ii)

            for s in combined_all.segments:
                seg_idx = int(s.seg_index)
                idxs = seg_to_frames.get(seg_idx, [])
                if not idxs:
                    continue

                # (A) Only for bond-change segments: pocket-only per-segment path
                if s.kind != "bridge" and s.summary and s.summary.strip() != "(no covalent changes detected)":
                    seg_imgs = [combined_all.images[j] for j in idxs]
                    seg_Es = [combined_all.energies[j] for j in idxs]
                    seg_trj = out_dir_path / f"mep_seg_{seg_idx:02d}_trj.xyz"
                    _write_xyz_trj_with_energy(seg_imgs, seg_Es, seg_trj)
                    click.echo(f"[write] Wrote per-segment pocket trajectory → '{seg_trj}'")
                    if ref_pdb_for_segments is not None:
                        _maybe_convert_to_pdb(seg_trj, ref_pdb_for_segments, out_path=out_dir_path / f"mep_seg_{seg_idx:02d}.pdb")

                # (B) HEI pocket files only for bond-change segments
                if s.kind != "bridge" and s.summary and s.summary.strip() != "(no covalent changes detected)":
                    energies_seg = [combined_all.energies[j] for j in idxs]
                    imax_rel = int(np.argmax(np.array(energies_seg, dtype=float)))
                    imax_abs = idxs[imax_rel]
                    hei_img = combined_all.images[imax_abs]
                    hei_E = [combined_all.energies[imax_abs]]
                    hei_trj = out_dir_path / f"hei_seg_{seg_idx:02d}.xyz"
                    _write_xyz_trj_with_energy([hei_img], hei_E, hei_trj)
                    click.echo(f"[write] Wrote segment HEI (pocket) → '{hei_trj}'")
                    if ref_pdb_for_segments is not None:
                        _maybe_convert_to_pdb(hei_trj, ref_pdb_for_segments, out_path=out_dir_path / f"hei_seg_{seg_idx:02d}.pdb")
        except Exception as e:
            click.echo(f"[write] WARNING: Failed to emit per-segment pocket outputs: {e}", err=True)
        # ---- END ----

        if do_merge:
            click.echo("\n=== Full-system merge (pocket → templates) started ===\n")
            # With --align True, use only the first reference PDB for all pairs (replicate it).
            if align:
                if not ref_pdb_paths or len(ref_pdb_paths) < 1:
                    raise click.BadParameter("--ref-pdb must provide at least one file when performing final merge with --align True.")
                first_ref = Path(ref_pdb_paths[0])
                ref_list_for_merge = [first_ref for _ in input_paths]  # replicate for all (length == n_inputs)
            else:
                ref_list_for_merge = [Path(p) for p in ref_pdb_paths]

            _merge_final_and_write(
                final_images=list(combined_all.images),
                pocket_inputs=[Path(p) for p in input_paths],
                ref_pdbs=ref_list_for_merge,
                segments=combined_all.segments,
                out_dir=out_dir_path
            )
            click.echo("\n=== Full-system merge finished ===\n")

        summary = {
            "out_dir": str(out_dir_path),
            "n_images": len(combined_all.images),
            "n_segments": len(combined_all.segments),
            "segments": [
                {
                    "index": int(s.seg_index),
                    "tag": s.tag,
                    "kind": s.kind,
                    "barrier_kcal": float(s.barrier_kcal),
                    "delta_kcal": float(s.delta_kcal),
                    "bond_changes": (s.summary if (s.kind != "bridge") else "")
                } for s in combined_all.segments
            ],
        }

        # --------------------------
        # 5) Console summary
        # --------------------------
        try:
            overall_changed, overall_summary = _has_bond_change(combined_all.images[0], combined_all.images[-1], bond_cfg)
        except Exception:
            overall_changed, overall_summary = False, ""

        click.echo("\n=== MEP Summary ===\n")

        click.echo("\n[overall] Covalent-bond changes between first and last image:")
        if overall_changed and overall_summary.strip():
            click.echo(textwrap.indent(overall_summary.strip(), prefix="  "))
        else:
            click.echo("  (no covalent changes detected)")

        if combined_all.segments:
            click.echo("\n[segments] Along the final MEP order (ΔE‡, ΔE). Bridges are shown between connected segments:")
            for i, seg in enumerate(combined_all.segments, 1):
                kind_label = "BRIDGE" if seg.kind == "bridge" else "SEG"
                click.echo(f"  [{i:02d}] ({kind_label}) {seg.tag}  |  ΔE‡ = {seg.barrier_kcal:.2f} kcal/mol,  ΔE = {seg.delta_kcal:.2f} kcal/mol")
                if seg.kind != "bridge" and seg.summary.strip():
                    click.echo(textwrap.indent(seg.summary.strip(), prefix="      "))
        else:
            click.echo("\n[segments] (no segment reports)")

        # --------------------------
        # 6) Energy diagram from bond-change segments (state labeling; compressed)
        # --------------------------
        diagram_payload: Optional[Dict[str, Any]] = None
        try:
            # Map each segment index → list of frame indices
            frame_seg_indices: List[int] = [int(getattr(im, "mep_seg_index", 0) or 0) for im in combined_all.images]
            seg_to_frames: Dict[int, List[int]] = {}
            for ii, sidx in enumerate(frame_seg_indices):
                if sidx <= 0:
                    continue
                seg_to_frames.setdefault(int(sidx), []).append(ii)

            # Build TS groups (each bond-change segment starts a group)
            ts_groups: List[Dict[str, Any]] = []
            ts_count = 0
            current: Optional[Dict[str, Any]] = None

            for s in combined_all.segments:
                idxs = seg_to_frames.get(int(s.seg_index), [])
                if not idxs:
                    continue

                if s.kind == "seg" and s.summary and s.summary.strip() != "(no covalent changes detected)":
                    # New TS group
                    ts_count += 1
                    imax = max(idxs, key=lambda j: combined_all.energies[j])
                    ts_e = float(combined_all.energies[imax])
                    first_im_e = float(combined_all.energies[idxs[-1]])
                    current = {
                        "ts_label": f"TS{ts_count}",
                        "ts_energy": ts_e,
                        "first_im_energy": first_im_e,
                        "tail_im_energy": first_im_e,
                        "has_extra": False,
                        "index": ts_count,
                    }
                    ts_groups.append(current)
                else:
                    # Kink/bridge: fold into current group as "extra" and update tail energy
                    if current is not None:
                        current["tail_im_energy"] = float(combined_all.energies[idxs[-1]])
                        current["has_extra"] = True
                    else:
                        # pre-TS region without bond change → ignore
                        pass

            # Clip endpoints to first/last bond-change segment edges
            start_idx_for_diag = 0
            end_idx_for_diag = len(combined_all.energies) - 1
            bc_segments_in_order: List[SegmentReport] = [
                s for s in combined_all.segments
                if (s.kind == "seg" and s.summary and s.summary.strip() != "(no covalent changes detected)")
            ]
            if bc_segments_in_order:
                first_bc = bc_segments_in_order[0]
                last_bc = bc_segments_in_order[-1]
                idxs_first_bc = seg_to_frames.get(int(first_bc.seg_index), [])
                idxs_last_bc = seg_to_frames.get(int(last_bc.seg_index), [])
                if idxs_first_bc:
                    start_idx_for_diag = int(idxs_first_bc[0])
                if idxs_last_bc:
                    end_idx_for_diag = int(idxs_last_bc[-1])

            # Compose compressed labels/energies & human-readable chain
            labels: List[str] = ["R"]
            energies_eh: List[float] = [float(combined_all.energies[start_idx_for_diag])]
            chain_tokens: List[str] = ["R"]

            for i, g in enumerate(ts_groups, start=1):
                last_group = (i == len(ts_groups))

                # TS
                labels.append(g["ts_label"])
                energies_eh.append(g["ts_energy"])
                chain_tokens.extend(["-->", g["ts_label"]])

                # For the last TS group: compress directly to P (no IMs)
                if last_group:
                    continue

                # IM1 (always keep)
                labels.append(f"IM{i}_1")
                energies_eh.append(g["first_im_energy"])
                chain_tokens.extend(["-->", f"IM{i}_1"])

                # IM2 (represent all extra kink/bridge before next TS)
                if g["has_extra"]:
                    labels.append(f"IM{i}_2")
                    energies_eh.append(g["tail_im_energy"])
                    chain_tokens.extend(["-|-->", f"IM{i}_2"])

            # Product
            labels.append("P")
            energies_eh.append(float(combined_all.energies[end_idx_for_diag]))
            chain_tokens.extend(["-->", "P"])

            # Convert to kcal/mol relative to R
            e0 = energies_eh[0]
            energies_kcal = [(e - e0) * AU2KCALPERMOL for e in energies_eh]
            energies_au = list(energies_eh)
            diagram_payload = {
                "name": "energy_diagram_MEP",
                "labels": list(labels),
                "energies_kcal": energies_kcal,
                "ylabel": "ΔE (kcal/mol)",
                "energies_au": energies_au,
                "image": str(out_dir_path / "energy_diagram_MEP.png"),
            }

            # Log exact inputs to build_energy_diagram, and the human-readable chain
            labels_repr = "[" + ", ".join(f'"{lab}"' for lab in labels) + "]"
            energies_repr = "[" + ", ".join(f"{val:.6f}" for val in energies_kcal) + "]"
            click.echo(f"[diagram] build_energy_diagram.labels = {labels_repr}")
            click.echo(f"[diagram] build_energy_diagram.energies_kcal = {energies_repr}")

            fig = build_energy_diagram(
                energies=energies_kcal,
                labels=labels,
                ylabel="ΔE (kcal/mol)",
                baseline=True,
                showgrid=False,
            )

            try:
                png_path = out_dir_path / "energy_diagram_MEP.png"
                fig.write_image(str(png_path), scale=2)
                click.echo(f"[diagram] Wrote energy diagram (PNG) → '{png_path}'")
            except Exception as e:
                click.echo(f"[diagram] NOTE: PNG export skipped (install 'kaleido' to enable): {e}", err=True)

            chain_text = " ".join(chain_tokens)
            click.echo(f"[diagram] State label sequence: {chain_text}")

        except Exception as e:
            click.echo(f"[diagram] WARNING: Failed to build energy diagram: {e}", err=True)

        # --------------------------
        # 7) Summary (YAML + log)
        # --------------------------
        if diagram_payload is not None:
            summary["energy_diagrams"] = [diagram_payload]

        with open(out_dir_path / "summary.yaml", "w") as f:
            yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)
        click.echo(f"[write] Wrote '{out_dir_path / 'summary.yaml'}'.")

        try:
            freeze_atoms_for_log: List[int] = []
            try:
                freeze_atoms_for_log = sorted(
                    {
                        int(i)
                        for g in getattr(combined_all, "images", [])
                        for i in getattr(g, "freeze_atoms", [])
                    }
                )
            except Exception:
                freeze_atoms_for_log = []

            diag_for_log: Dict[str, Any] = diagram_payload or {}
            mep_info = {
                "n_images": len(combined_all.images),
                "n_segments": len(combined_all.segments),
                "traj_pdb": str(out_dir_path / "mep.pdb") if (out_dir_path / "mep.pdb").exists() else None,
                "mep_plot": str(out_dir_path / "mep_plot.png") if (out_dir_path / "mep_plot.png").exists() else None,
                "diagram": diag_for_log,
            }
            summary_payload = {
                "root_out_dir": str(out_dir_path),
                "path_dir": str(out_dir_path),
                "path_module_dir": "path_search",
                "pipeline_mode": "path-search",
                "refine_path": True,
                "tsopt": False,
                "thermo": False,
                "dft": False,
                "opt_mode": opt_mode,
                "mep_mode": "path-search",
                "uma_model": calc_cfg.get("uma_model"),
                "command": command_str,
                "charge": calc_cfg.get("model_charge"),
                "spin": calc_cfg.get("model_mult"),
                "freeze_atoms": freeze_atoms_for_log,
                "mep": mep_info,
                "segments": summary.get("segments", []),
                "energy_diagrams": summary.get("energy_diagrams", []),
                "key_files": {},
            }
            write_summary_log(out_dir_path / "summary.log", summary_payload)
            click.echo(f"[write] Wrote '{out_dir_path / 'summary.log'}'.")
        except Exception as e:
            click.echo(f"[write] WARNING: Failed to write summary.log: {e}", err=True)

        # summary.md and key_* outputs are disabled.
        # --------------------------
        # 8) Elapsed time
        # --------------------------
        click.echo(format_elapsed("[time] Elapsed for Path Search", time_start))

    except ZeroStepLength:
        click.echo("ERROR: Proposed step length dropped below the minimum allowed (ZeroStepLength).", err=True)
        sys.exit(2)
    except OptimizationError as e:
        click.echo(f"ERROR: Path search failed — {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled error during path search:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)
    finally:
        for prepared in prepared_inputs:
            prepared.cleanup()
