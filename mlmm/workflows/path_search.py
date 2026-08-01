"""Recursive ML/MM GSM/DMF paths with multistep segmentation."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import gc
import logging
import sys
import textwrap
import tempfile
import os

logger = logging.getLogger(__name__)
import time  # timing
import re    # used in _segment_base_id

import click
from mlmm.core.output import emit
import numpy as np
import torch

from pysisyphus.helpers import geom_loader
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import AU2KCALPERMOL, BOHR2ANG


from mlmm.backends.mlmm_calc import mlmm, MLMMASECalculator
from mlmm.core.defaults import (
    BOND_KW as _BOND_KW_DEFAULT,
    fresh_dmf_config,
    OUT_DIR_PATH_SEARCH,
    SEARCH_KW as _SEARCH_KW_DEFAULT,
    THRESH_CHOICES,
)
from mlmm.workflows.path_opt import (
    GS_KW as _PATH_GS_KW,
    STOPT_KW as _PATH_STOPT_KW,
    DMF_KW as _PATH_DMF_KW,
    _select_hei_index,
    _shared_frozen_reference,
    resolve_dmf_solve_tol,
)
from mlmm.workflows.opt import (
    GEOM_KW as _OPT_GEOM_KW,
    CALC_KW as _OPT_CALC_KW,
    LBFGS_KW as _OPT_LBFGS_KW,
    _parse_freeze_atoms as _parse_freeze_atoms_opt,
    _normalize_geom_freeze as _normalize_geom_freeze_opt,
)
from mlmm.workflows.opt import _convert_yaml_layer_atoms_1to0
from mlmm.workflows._outcomes import optimizer_converged_bit
from mlmm.workflows.charge_prep import resolve_charge_spin_or_raise
from mlmm.core.utils import (
    apply_layer_freeze_constraints,
    apply_ref_pdb_override,
    convert_xyz_to_pdb,
    is_convert_file_enabled,
    set_convert_file_enabled,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    strip_inherited_keys,
    filter_calc_for_echo,
    format_freeze_atoms_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    build_energy_diagram,
    prepare_input_structure,
    PreparedInputStructure,
    validate_endpoint_atom_identities,
    parse_indices_string,
    resolve_ml_layer_assignment,
)
from mlmm.core.result_commit import commit_json_exact, with_current_run_id
from mlmm.cli.common_options import add_ml_layer_detection_options, add_precision_option, add_workers_options, add_backend_model_option, add_calc_file_option, add_deterministic_option, add_allow_charge_mult_mismatch_option
from mlmm.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, make_is_param_explicit, _write_error_json, render_cli_exception
from mlmm.cli.preflight import validate_existing_files
from mlmm.io.trj2fig import run_trj2fig  # auto-generate an energy plot when a _trj.xyz is produced
from mlmm.io.summary import emit_method_citations, method_references, write_summary_log
from mlmm.domain.bond_changes import compare_structures, summarize_changes
from mlmm.workflows.align_freeze import (
    align_and_refine_sequence_inplace,
    alignment_failed_pair_indices,
)


# Geometry (input handling) — reuse opt.py defaults
GEOM_KW: Dict[str, Any] = deepcopy(_OPT_GEOM_KW)

# ML/MM calculator settings — reuse opt.py defaults
CALC_KW: Dict[str, Any] = deepcopy(_OPT_CALC_KW)

# GrowingString (path representation)
GS_KW: Dict[str, Any] = deepcopy(_PATH_GS_KW)

# StringOptimizer (GSM optimization control)
STOPT_KW: Dict[str, Any] = deepcopy(_PATH_STOPT_KW)
STOPT_KW.update({
    "out_dir": OUT_DIR_PATH_SEARCH,
})

# LBFGS settings
LBFGS_KW: Dict[str, Any] = deepcopy(_OPT_LBFGS_KW)
LBFGS_KW.update({
    "out_dir": OUT_DIR_PATH_SEARCH,
})

# Covalent-bond change detection
BOND_KW: Dict[str, Any] = deepcopy(_BOND_KW_DEFAULT)

# DMF (Direct Max Flux) defaults
DMF_KW: Dict[str, Any] = deepcopy(_PATH_DMF_KW)

# Global search control
SEARCH_KW: Dict[str, Any] = deepcopy(_SEARCH_KW_DEFAULT)


def _normalized_path_tangent(
    coords: Sequence[np.ndarray],
    index: int,
    energies: Optional[Sequence[float]] = None,
) -> Optional[np.ndarray]:
    """Return a normalized Cartesian tangent for a path image.

    Finite image energies enable the improved upwind tangent used by the
    chain-of-states optimizer.  Trajectories without energies fall back to a
    spacing-independent secant bisector; endpoints use their only secant.
    """
    if len(coords) < 2 or not 0 <= int(index) < len(coords):
        return None

    arrays = [np.asarray(item, dtype=float).reshape(-1) for item in coords]

    def _unit(vector: np.ndarray) -> Optional[np.ndarray]:
        norm = float(np.linalg.norm(vector))
        if not np.isfinite(norm) or norm <= 0.0:
            return None
        return vector / norm

    if index == 0:
        return _unit(arrays[1] - arrays[0])
    if index == len(arrays) - 1:
        return _unit(arrays[-1] - arrays[-2])

    incoming_raw = arrays[index] - arrays[index - 1]
    outgoing_raw = arrays[index + 1] - arrays[index]
    incoming = _unit(incoming_raw)
    outgoing = _unit(outgoing_raw)
    if incoming is None:
        return outgoing
    if outgoing is None:
        return incoming

    if energies is not None and len(energies) == len(arrays):
        energy_values = np.asarray(energies, dtype=float)
        previous = float(energy_values[index - 1])
        current = float(energy_values[index])
        following = float(energy_values[index + 1])
        if np.all(np.isfinite((previous, current, following))):
            if following > current > previous:
                tangent_raw = outgoing_raw
            elif following < current < previous:
                tangent_raw = incoming_raw
            else:
                next_delta = abs(following - current)
                previous_delta = abs(previous - current)
                delta_max = max(next_delta, previous_delta)
                delta_min = min(next_delta, previous_delta)
                if following >= previous:
                    tangent_raw = outgoing_raw * delta_max + incoming_raw * delta_min
                else:
                    tangent_raw = outgoing_raw * delta_min + incoming_raw * delta_max
            tangent = _unit(tangent_raw)
            if tangent is not None:
                return tangent

    tangent = _unit(incoming + outgoing)
    if tangent is not None:
        return tangent
    return _unit(arrays[index + 1] - arrays[index - 1])


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
        cfg: Dict[str, Any] = {"freeze_atoms": list(base_freeze)}
        freeze = merge_freeze_atom_indices(cfg)
        g = geom_loader(geom_path, coord_type=coord_type, freeze_atoms=freeze)
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
        # path-search set the flag but never read it, so --no-convert-files was a no-op here.
        if not is_convert_file_enabled():
            return None
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
    if A.shape != B.shape or A.ndim != 2 or A.shape[1] != 3:
        raise ValueError(
            f"RMSD coordinate shapes must match (N, 3); got {A.shape} and {B.shape}."
        )
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
        g = geom_loader(Path(tmp.name), coord_type=coord_type, freeze_atoms=freeze_atoms)
        g.freeze_atoms = np.array(sorted(set(map(int, freeze_atoms))), dtype=int)
        return g
    finally:
        try:
            os.unlink(tmp.name)
        except Exception:
            logger.debug("Failed to unlink temp file %s", tmp.name, exc_info=True)


def _make_linear_interpolations(gL, gR, n_internal: int) -> List[Any]:
    """
    Return `n_internal` linearly interpolated structures between gL → gR (excluding endpoints).
    Atom order follows `gL`.
    """
    A = np.asarray(gL.coords3d, dtype=float)
    B = np.asarray(gR.coords3d, dtype=float)
    if A.shape != B.shape or A.ndim != 2 or A.shape[1] != 3:
        raise ValueError(
            "Interpolation coordinate shapes must match (N, 3); "
            f"got {A.shape} and {B.shape}."
        )
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
                logger.debug("Failed to set attribute %s on image", k, exc_info=True)


def _frame_ranges_by_segment(images: Sequence[Any]) -> Dict[int, Dict[str, Any]]:
    """Return additive half-open frame ranges for each tagged MEP segment."""

    result: Dict[int, Dict[str, Any]] = {}
    run_index: Optional[int] = None
    run_kind = ""
    run_start = 0

    def close_run(stop: int) -> None:
        if run_index is None or run_index <= 0:
            return
        entry = result.setdefault(
            run_index,
            {"kind": run_kind or "seg", "frame_ranges": []},
        )
        entry["frame_ranges"].append([run_start, stop])

    for frame_index, image in enumerate(images):
        segment_index = int(getattr(image, "mep_seg_index", 0) or 0)
        segment_kind = str(getattr(image, "mep_seg_kind", "") or "seg")
        if (segment_index, segment_kind) != (run_index, run_kind):
            close_run(frame_index)
            run_index = segment_index
            run_kind = segment_kind
            run_start = frame_index
    close_run(len(images))

    for entry in result.values():
        ranges = entry["frame_ranges"]
        if len(ranges) == 1:
            entry["frame_start"], entry["frame_stop"] = ranges[0]
    return result


def _segment_base_id(tag: str) -> str:
    """
    Extract base id 'seg_XXX' from a tag like 'seg_000_refine'; fallback to `tag` or 'seg'.
    """
    m = re.search(r"(seg_\d{3})", tag or "")
    return m.group(1) if m else (tag or "seg")


def _select_hei_index(energies: Sequence[float]) -> int:
    """Return the global highest-energy-image index."""
    E = np.asarray(energies, dtype=float)
    if E.size == 0:
        raise ValueError("Cannot select an HEI from an empty energy profile.")
    if not np.all(np.isfinite(E)):
        raise ValueError("Cannot select an HEI from non-finite energies.")
    return int(np.argmax(E))


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
    # reported convergence of the string/DMF optimizer that produced this
    # MEP. ``None`` means no readable convergence signal (fail-closed: never
    # promoted to a usable segment by artifact existence alone).
    is_converged: Optional[bool] = None


# ---- Per-segment summary for the console report ----
@dataclass
class SegmentReport:
    tag: str
    barrier_kcal: float
    delta_kcal: float
    summary: str  # summarize_changes string (empty for bridges)
    kind: str = "seg"          # "seg" or "bridge"
    seg_index: int = 0         # 1-based index along final MEP (assigned later)
    # the segment's optimizer convergence, threaded from GSMResult. A
    # reactive segment whose optimizer did not explicitly converge is unusable
    # and cannot make the path aggregate a scientific success.
    converged: Optional[bool] = None


def _run_gsm_between(
    gA,
    gB,
    shared_calc,
    gs_cfg: Dict[str, Any],
    stopt_cfg: Dict[str, Any],
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

    _opt_args = dict(stopt_cfg)
    seg_dir = out_dir / f"{tag}_mep"
    seg_dir.mkdir(parents=True, exist_ok=True)
    _opt_args["out_dir"] = str(seg_dir)

    optimizer = StringOptimizer(
        geometry=gs,
        **{k: v for k, v in _opt_args.items() if k != "type"}
    )

    optimizer.run()
    # a normal (non-raising) run is NOT convergence — capture the
    # StringOptimizer's explicit bit so a max-cycle segment cannot be promoted.
    _gsm_converged = optimizer_converged_bit(optimizer)

    energies = list(map(float, np.array(gs.energy, dtype=float)))
    images = list(gs.images)

    try:
        hei_idx = _select_hei_index(energies)
    except ValueError as exc:
        raise click.ClickException(f"{tag}: {exc}") from exc

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
            emit(f"[{tag}] Saved energy plot → '{seg_dir / 'mep_plot.png'}'", detail=True)
        else:
            click.echo(f"[{tag}] WARNING: Energies missing; skipping plot.", err=True)
    except Exception as e:
        click.echo(f"[{tag}] WARNING: Failed to plot energy: {e}", err=True)

    # If PDB input exists, convert intermediate _trj.xyz to PDB
    _maybe_convert_to_pdb(final_trj, ref_pdb_path, seg_dir / "final_geometries.pdb")

    # Write HEI structure (XYZ with energy in line 2)
    try:
        hei_geom = images[hei_idx]
        hei_E = float(energies[hei_idx])
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

    return GSMResult(images=images, energies=energies, hei_idx=hei_idx, is_converged=_gsm_converged)


def _run_dmf_between(
    gA,
    gB,
    shared_calc,
    calc_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],
    max_nodes: int,
    dmf_cfg: Optional[Dict[str, Any]],
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
        logger.debug("Failed to extract freeze_atoms from endpoints", exc_info=True)

    # Convert pysisyphus geometries to ASE Atoms for DMF
    def _geom_to_ase(g):
        return ase_read(StringIO(g.as_xyz()), format="xyz")

    geoms_for_dmf = [gA, gB]

    dmf_backend = str((dmf_cfg or {}).get("backend", "gpu")).strip().lower()
    try:
        from ase.calculators.mixing import SumCalculator
        if dmf_backend == "cpu":
            from dmf import DirectMaxFlux, interpolate_fbenm
        else:
            from dmf.torch import DirectMaxFlux, interpolate_fbenm
    except Exception as e:
        raise RuntimeError(
            "DMF mode requires cyipopt and pydmf>=1.2 "
            "(`conda install -c conda-forge cyipopt -y`, then `pip install "
            "'pydmf[torch]>=1.2'` for GPU or `pip install 'pydmf>=1.2'` for CPU): "
            f"{e}"
        ) from e

    from mlmm.workflows.restraints import HarmonicFixAtoms
    from mlmm.core.utils import is_verbose

    ref_images = [_geom_to_ase(g) for g in geoms_for_dmf]
    fix_ref_positions = _shared_frozen_reference(ref_images, fix_atoms)
    charge = int(calc_cfg.get("model_charge", 0))
    spin = int(calc_cfg.get("model_mult", 1))
    for img in ref_images:
        img.info["charge"] = charge
        img.info["spin"] = spin

    # Build ASE calculator from the shared PySisyphus calculator
    ase_calc = MLMMASECalculator(core=shared_calc.core)

    dmf_cfg_local = fresh_dmf_config(dmf_cfg)
    fbenm_opts = dict(dmf_cfg_local.get("fbenm_options", {}))
    cfbenm_opts = dict(dmf_cfg_local.get("cfbenm_options", {}))
    dmf_opts = dict(dmf_cfg_local.get("dmf_options", {}))
    ipopt_opts = dict(dmf_cfg_local.get("ipopt_options", {}))
    if "print_level" not in ipopt_opts and not is_verbose():
        ipopt_opts["print_level"] = 0
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
        ipopt_options=ipopt_opts,
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
            harmonic_calc = HarmonicFixAtoms(
                indices=fix_atoms,
                ref_positions=fix_ref_positions,
                k_fix=k_fix,
            )
            image.calc = SumCalculator([ase_calc, harmonic_calc])
        else:
            image.calc = ase_calc

    dmf_ipopt_opts = dict(ipopt_opts)
    dmf_ipopt_opts["output_file"] = str(seg_dir / "dmf_ipopt.out")
    mxflx.add_ipopt_options(dmf_ipopt_opts)
    max_cycles_dmf = dmf_cfg_local.get("max_cycles")
    if max_cycles_dmf is not None:
        try:
            max_iter = int(max_cycles_dmf)
            if max_iter > 0:
                mxflx.add_ipopt_options({"max_iter": max_iter})
        except Exception:
            logger.debug("Failed to set ipopt max_iter option", exc_info=True)

    # IPOPT status 0/1 = converged; any other (max-iter, infeasible) is
    # not. A missing status fails closed to unknown so DMF artifact existence
    # never promotes a nonconverged solve.
    from mlmm.workflows._outcomes import ipopt_status_to_converged
    _dmf_solve_ret = mxflx.solve(
        tol=resolve_dmf_solve_tol(dmf_cfg_local, prefix="[path-search]")
    )
    _dmf_info = (
        _dmf_solve_ret[1]
        if isinstance(_dmf_solve_ret, (tuple, list)) and len(_dmf_solve_ret) >= 2
        and isinstance(_dmf_solve_ret[1], dict)
        else {}
    )
    _dmf_status = _dmf_info.get("status")
    _dmf_converged, _ = ipopt_status_to_converged(_dmf_status)

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
        emit(f"[{tag}] Saved energy plot → '{seg_dir / 'mep_plot.png'}'", detail=True)
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
        try:
            g = gl(tmp_xyz, coord_type=gA.coord_type)
            try:
                g.freeze_atoms = np.array(getattr(gA, "freeze_atoms", []), dtype=int)
            except Exception:
                logger.debug(
                    "Failed to set freeze_atoms on interpolated image",
                    exc_info=True,
                )
            g.set_calculator(shared_calc)
            imgs.append(g)
        finally:
            tmp_xyz.unlink(missing_ok=True)

    return GSMResult(images=imgs, energies=energies, hei_idx=hei_idx, is_converged=_dmf_converged)


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
    stopt_cfg: Dict[str, Any],
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
            dmf_cfg=dmf_cfg,
        )
    return _run_gsm_between(gA, gB, shared_calc, gs_cfg, stopt_cfg, out_dir, tag=tag, ref_pdb_path=ref_pdb_path)


def _optimize_single(
    g,
    shared_calc,
    lbfgs_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],  # for PDB conversion
):
    """
    Run single-structure optimization (LBFGS) and return ``(geometry, converged)``.

    ``converged`` is the optimizer's fail-closed tri-state convergence bit
    (:func:`optimizer_converged_bit`): a non-raising ``run()`` is NOT convergence
    so the caller must gate on this bit rather than assume a completed
    optimization converged (e.g. so a nonconverged kink-segment single-structure
    opt cannot silently become a usable reactive leaf).
    """
    emit(f"\n====== [{tag}] Single-structure LBFGS ======\n", narrative=True)
    g.set_calculator(shared_calc)

    seg_dir = out_dir / f"{tag}_lbfgs_opt"
    seg_dir.mkdir(parents=True, exist_ok=True)
    args = dict(lbfgs_cfg)
    args["out_dir"] = str(seg_dir)

    opt = LBFGS(g, **args)

    opt.run()
    converged = optimizer_converged_bit(opt)

    try:
        final_xyz = Path(opt.final_fn) if isinstance(opt.final_fn, (str, Path)) else Path(opt.final_fn)
        _maybe_convert_to_pdb(final_xyz, ref_pdb_path)
        g_final = geom_loader(
            final_xyz,
            coord_type=g.coord_type,
            freeze_atoms=getattr(g, "freeze_atoms", []),
        )
        try:
            g_final.freeze_atoms = np.array(getattr(g, "freeze_atoms", []), dtype=int)
        except Exception:
            logger.debug("Failed to set freeze_atoms on final geometry", exc_info=True)
        g_final.set_calculator(shared_calc)
        return g_final, converged
    except Exception:
        return g, converged


def _refine_between(
    gL,
    gR,
    shared_calc,
    gs_cfg: Dict[str, Any],
    stopt_cfg: Dict[str, Any],
    out_dir: Path,
    tag: str,
    ref_pdb_path: Optional[Path],  # for PDB conversion
    mep_mode_kind: str = "gsm",
    calc_cfg: Optional[Dict[str, Any]] = None,
    max_nodes: int = 10,
    dmf_cfg: Optional[Dict[str, Any]] = None,
) -> GSMResult:
    """
    Refine End1–End2 via GSM or DMF using the resolved climbing policy.
    """
    return _run_mep_between(
        gL, gR, shared_calc, gs_cfg, stopt_cfg, out_dir, tag=f"{tag}_refine",
        ref_pdb_path=ref_pdb_path, mep_mode_kind=mep_mode_kind,
        calc_cfg=calc_cfg, max_nodes=max_nodes, dmf_cfg=dmf_cfg,
    )


def _maybe_bridge_segments(
    tail_g,
    head_g,
    shared_calc,
    gs_cfg: Dict[str, Any],  # bridge-specific GS config
    stopt_cfg: Dict[str, Any],
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
    emit(f"[{tag}] Gap detected between segments (RMSD={rmsd:.4e} bohr) — bridging via {mep_mode_kind.upper()}.", narrative=True)
    return _run_mep_between(
        tail_g, head_g, shared_calc, gs_cfg, stopt_cfg, out_dir, tag=f"{tag}_bridge",
        ref_pdb_path=ref_pdb_path, mep_mode_kind=mep_mode_kind,
        calc_cfg=calc_cfg, max_nodes=max_nodes, dmf_cfg=dmf_cfg,
    )


def _stitch_paths(
    parts: List[Tuple[List[Any], List[float]]],
    stitch_rmsd_thresh: float,
    bridge_rmsd_thresh: float,
    shared_calc,
    gs_cfg,   # GS config for bridges (climb=False, max_nodes=search.max_nodes_bridge)
    stopt_cfg,
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
            except Exception as _bond_exc:
                click.echo(
                    "[path-search] WARNING: the interface bond-change check failed "
                    f"({_bond_exc}); treating the interface as unchanged, so a reaction step "
                    "may be missing from the path.",
                    err=True,
                )
                adj_changed, adj_summary = False, ""

        if adj_changed and segment_builder is not None:
            emit(f"[{tag}] Covalent changes detected at interface — inserting a new recursive segment.", narrative=True)
            if adj_summary:
                click.echo(textwrap.indent(adj_summary, prefix="  "))
            sub = segment_builder(tail, head, f"{tag}_mid")
            seg_imgs, seg_E = sub.images, sub.energies
            if segments_out is not None and getattr(sub, "segments", None):
                right_tag = _first_known_seg_tag_from_images(imgs)
                insert_pos = next(
                    (
                        index
                        for index, report in enumerate(segments_out)
                        if report.tag == right_tag
                    ),
                    len(segments_out),
                )
                segments_out[insert_pos:insert_pos] = list(sub.segments)
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
                tail, head, shared_calc, gs_cfg, stopt_cfg, out_dir, tag=bridge_name_base,
                rmsd_thresh=bridge_rmsd_thresh, ref_pdb_path=ref_pdb_path,
                mep_mode_kind=mep_mode_kind, calc_cfg=calc_cfg,
                max_nodes=int(gs_cfg.get("max_nodes", 5)), dmf_cfg=dmf_cfg,
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
                        kind="bridge",
                        converged=getattr(br, "is_converged", None),
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


# Recursive search (core)

@dataclass
class CombinedPath:
    images: List[Any]
    energies: List[float]
    segments: List[SegmentReport]  # segment summaries in final output order


def _path_leaves_and_expected(
    segments: Sequence[SegmentReport],
    *,
    raw_artifacts: Sequence[str] = (),
    engine_converged: Optional[bool] = True,
):
    """Build path :class:`LeafOutcome` list and expected reactive-segment IDs.

    Reactive segments (``kind != "bridge"``) are required leaves; bridges are
    optional connectors.  When there is no reactive segment at all — the
    endpoint-HEI branch returns ``segments=[]`` even though an R/P energy diagram
    can still be drawn — an unusable ``raw_path`` leaf is emitted so the aggregate
    mapper cannot promote the diagnostic diagram to success.  The raw
    trajectory/diagram remain reportable as artifacts.
    """

    from mlmm.workflows._outcomes import LeafOutcome, make_leaf

    leaves: List[Any] = []
    reactive = [s for s in segments if getattr(s, "kind", "seg") != "bridge"]
    for s in segments:
        is_reactive = getattr(s, "kind", "seg") != "bridge"
        # a reactive segment is usable only when its optimizer explicitly
        # converged. A nonconverged (max-cycle) StringOptimizer segment retains
        # its trajectory artifact but must not count toward completeness.
        _seg_conv = getattr(s, "converged", None)
        leaves.append(
            make_leaf(
                "path",
                f"segment_{int(s.seg_index)}",
                required=is_reactive,
                executed=True,
                converged=_seg_conv,
            )
        )
    expected = [f"segment_{int(s.seg_index)}" for s in reactive]
    if not reactive:
        reason = "endpoint_hei"
        if engine_converged is False:
            reason = "endpoint_hei;engine_nonconverged"
        leaves.append(
            LeafOutcome(
                stage="path",
                item_id="raw_path",
                required=True,
                executed=True,
                converged=engine_converged if isinstance(engine_converged, bool) else None,
                usable=False,
                reason=reason,
                artifacts=tuple(str(a) for a in raw_artifacts),
            )
        )
    return leaves, expected


def _enrich_path_summary_contract(
    summary: Dict[str, Any],
    *,
    segments: Sequence[SegmentReport],
    out_dir: Path,
    calc_cfg: Dict[str, Any],
    command: str,
) -> Dict[str, Any]:
    """Attach the fail-closed machine contract for standalone path-search."""

    try:
        from mlmm._version import __version__
    except Exception:
        __version__ = "unknown"
    from mlmm.core.utils import (
        RESULT_JSON_SCHEMA_VERSION,
        calculator_provenance,
    )
    from mlmm.workflows._outcomes import (
        aggregate_workflow_truth,
        attach_outcomes,
    )

    summary["mlmm_toolkit_version"] = __version__
    summary["schema_version"] = RESULT_JSON_SCHEMA_VERSION
    summary["pipeline_mode"] = "path-search"
    raw_artifacts = [
        name
        for name in (
            "mep.pdb",
            "mep.cif",
            "mep_plot.png",
            "energy_diagram_MEP.png",
        )
        if (out_dir / name).exists()
    ]
    path_leaves, path_expected = _path_leaves_and_expected(
        segments,
        raw_artifacts=raw_artifacts,
    )
    path_truth = aggregate_workflow_truth(path_leaves, path_expected)

    reactive = [
        segment
        for segment in segments
        if getattr(segment, "kind", "seg") != "bridge"
    ]
    legacy_status = "success" if summary.get("energy_diagrams") else "partial"
    if legacy_status == "success" and not reactive:
        legacy_status = "partial"
    summary["status"] = legacy_status
    attach_outcomes(summary, truth=path_truth, stage_outcomes=path_leaves)
    summary.update(calculator_provenance(calc_cfg))
    summary["charge"] = calc_cfg.get("model_charge")
    summary["spin"] = calc_cfg.get("model_mult")
    summary["command"] = command
    try:
        from mlmm.core.utils import _collect_environment_info

        summary.setdefault("environment", _collect_environment_info())
    except Exception:
        pass
    return summary


def _summary_log_provenance(summary: Dict[str, Any]) -> Dict[str, Any]:
    """Copy provenance and outcome truth from the enriched run summary."""

    return {
        "mlip_backend": summary.get("mlip_backend"),
        "mlip_model": summary.get("mlip_model"),
        "mlip_precision": summary.get("mlip_precision"),
        "status": summary.get("status"),
        "status_reasons": summary.get("status_reasons", []),
        "execution_status": summary.get("execution_status"),
        "scientific_status": summary.get("scientific_status"),
        "scientific_status_reasons": summary.get(
            "scientific_status_reasons", []
        ),
    }


def _trailing_kink_count(segments: Sequence[SegmentReport]) -> int:
    """Return the number of consecutive kink segments at the end of *segments*."""
    count = 0
    for seg in reversed(segments):
        if seg.tag and "kink" in seg.tag:
            count += 1
        else:
            break
    return count


def _build_multistep_path(
    gA,
    gB,
    shared_calc,
    geom_cfg: Dict[str, Any],
    gs_cfg: Dict[str, Any],
    stopt_cfg: Dict[str, Any],
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
    kink_seq_count: int = 0,
) -> CombinedPath:
    """
    Recursively construct a multistep MEP from A–B and return it (A→B order).
    """
    seg_max_nodes = int(
        search_cfg.get(
            "max_nodes_segment",
            gs_cfg.get("max_nodes", GS_KW["max_nodes"]),
        )
    )
    gs_seg_cfg = {**gs_cfg, "max_nodes": seg_max_nodes}
    max_seq_kink = int(search_cfg.get("max_seq_kink", 2))

    if depth > int(search_cfg.get("max_depth", 10)):
        click.echo(f"[{branch_tag}] Reached maximum recursion depth. Returning current endpoints only.")
        gsm = _run_mep_between(
            gA, gB, shared_calc, gs_seg_cfg, stopt_cfg, out_dir, tag=f"seg_{seg_counter[0]:03d}_maxdepth",
            ref_pdb_path=ref_pdb_path, mep_mode_kind=mep_mode_kind,
            calc_cfg=calc_cfg, max_nodes=seg_max_nodes, dmf_cfg=dmf_cfg,
        )
        seg_counter[0] += 1
        _tag_images(gsm.images, pair_index=pair_index)
        return CombinedPath(images=gsm.images, energies=gsm.energies, segments=[])

    seg_id = seg_counter[0]
    seg_counter[0] += 1
    tag0 = f"seg_{seg_id:03d}"

    gsm0 = _run_mep_between(
        gA, gB, shared_calc, gs_seg_cfg, stopt_cfg, out_dir, tag=tag0,
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
        emit(f"[{tag0}] Refining HEI±1 (peak mode).", narrative=True)

    left_end, left_conv = _optimize_single(left_img, shared_calc, single_opt_cfg, out_dir, tag=f"{tag0}_left", ref_pdb_path=ref_pdb_path)
    right_end, right_conv = _optimize_single(right_img, shared_calc, single_opt_cfg, out_dir, tag=f"{tag0}_right", ref_pdb_path=ref_pdb_path)

    try:
        lr_changed, _ = _has_bond_change(left_end, right_end, bond_cfg)
    except Exception as e:
        click.echo(f"[{tag0}] WARNING: Failed to evaluate bond changes for kink detection: {e}", err=True)
        lr_changed, _ = True, ""
    use_kink = (not lr_changed)

    if use_kink:
        n_inter = int(search_cfg.get("kink_max_nodes", 3))
        emit(f"[{tag0}] Kink detected (no covalent changes between End1 and End2). "
                   f"Using {n_inter} linear interpolation nodes + single-structure optimizations instead of GSM.",
                   narrative=True)
        inter_geoms = _make_linear_interpolations(left_end, right_end, n_inter)
        opt_inters: List[Any] = []
        inter_convs: List[Optional[bool]] = []
        for i, g_int in enumerate(inter_geoms, 1):
            g_int.set_calculator(shared_calc)
            g_opt, _inter_conv = _optimize_single(g_int, shared_calc, single_opt_cfg, out_dir, tag=f"{tag0}_kink_int{i}", ref_pdb_path=ref_pdb_path)
            opt_inters.append(g_opt)
            inter_convs.append(_inter_conv)
        step_imgs = [left_end] + opt_inters + [right_end]
        step_E = [float(img.energy) for img in step_imgs]
        _kink_hei = int(np.argmax(step_E[1:-1])) + 1 if len(step_E) > 2 else int(np.argmax(step_E))
        # a kink segment is assembled from single-structure optimizations
        # (not a StringOptimizer). It is usable only when EVERY endpoint/inter
        # optimization explicitly converged; fold their convergence rather than
        # hardcode True, so a nonconverged single-structure opt cannot silently
        # become a usable reactive leaf (fail-closed).
        from mlmm.workflows._outcomes import combine_step_convergence
        _kink_converged = combine_step_convergence(
            [left_conv] + inter_convs + [right_conv]
        )
        ref1 = GSMResult(images=step_imgs, energies=step_E, hei_idx=_kink_hei,
                         is_converged=_kink_converged)
        step_tag_for_report = f"{tag0}_kink"
    else:
        ref1 = _refine_between(left_end, right_end, shared_calc, gs_seg_cfg, stopt_cfg, out_dir, tag=tag0,
                               ref_pdb_path=ref_pdb_path, mep_mode_kind=mep_mode_kind,
                               calc_cfg=calc_cfg, max_nodes=seg_max_nodes, dmf_cfg=dmf_cfg)
        step_tag_for_report = f"{tag0}_refine"

    step_imgs, step_E = ref1.images, ref1.energies

    _changed, step_summary = _has_bond_change(step_imgs[0], step_imgs[-1], bond_cfg)
    _tag_images(step_imgs, mep_seg_tag=step_tag_for_report, mep_seg_kind="seg",
                mep_has_bond_changes=bool(_changed), pair_index=pair_index)

    left_changed, left_summary = _has_bond_change(gA, left_end, bond_cfg)
    right_changed, right_summary = _has_bond_change(right_end, gB, bond_cfg)

    emit(f"[{tag0}] Covalent changes (A vs left_end): {'Yes' if left_changed else 'No'}", narrative=True)
    if left_changed:
        click.echo(textwrap.indent(left_summary, prefix="  "))
    emit(f"[{tag0}] Covalent changes (right_end vs B): {'Yes' if right_changed else 'No'}", narrative=True)
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
        kind="seg",
        converged=getattr(ref1, "is_converged", None),
    )

    parts: List[Tuple[List[Any], List[float]]] = []
    seg_reports: List[SegmentReport] = []

    trailing_kink_run = kink_seq_count
    if left_changed:
        subL = _build_multistep_path(
            gA, left_end, shared_calc, geom_cfg, gs_cfg, stopt_cfg,
            single_opt_cfg, bond_cfg, search_cfg, refine_mode_kind,
            out_dir, ref_pdb_path, depth + 1, seg_counter, branch_tag=f"{branch_tag}L",
            pair_index=pair_index,
            mep_mode_kind=mep_mode_kind, calc_cfg=calc_cfg, dmf_cfg=dmf_cfg,
            kink_seq_count=kink_seq_count,
        )
        _tag_images(subL.images, pair_index=pair_index)
        parts.append((subL.images, subL.energies))
        seg_reports.extend(subL.segments)
        trailing_kink_run = _trailing_kink_count(seg_reports)

    current_kink_run = trailing_kink_run + 1 if use_kink else 0
    if use_kink and current_kink_run >= max_seq_kink:
        warning_msg = (
            f"[{tag0}] Consecutive kink segments were detected. Something seems wrong. "
            "Please check the initial structure and the generated intermediate structures. "
            "Alternatively, try switching the mep-mode. If that still fails, try including intermediate structures in the inputs."
        )
        click.echo(warning_msg)
        gsm = _run_mep_between(
            gA, gB, shared_calc, gs_seg_cfg, stopt_cfg, out_dir, tag=f"seg_{seg_counter[0]:03d}_kinklimit",
            ref_pdb_path=ref_pdb_path, mep_mode_kind=mep_mode_kind,
            calc_cfg=calc_cfg, max_nodes=seg_max_nodes, dmf_cfg=dmf_cfg,
        )
        seg_counter[0] += 1
        _tag_images(gsm.images, pair_index=pair_index)
        return CombinedPath(images=gsm.images, energies=gsm.energies, segments=[])

    parts.append((step_imgs, step_E))
    seg_reports.append(seg_report)

    if right_changed:
        subR = _build_multistep_path(
            right_end, gB, shared_calc, geom_cfg, gs_cfg, stopt_cfg,
            single_opt_cfg, bond_cfg, search_cfg, refine_mode_kind,
            out_dir, ref_pdb_path, depth + 1, seg_counter, branch_tag=f"{branch_tag}R",
            pair_index=pair_index,
            mep_mode_kind=mep_mode_kind, calc_cfg=calc_cfg, dmf_cfg=dmf_cfg,
            kink_seq_count=current_kink_run,
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
            geom_cfg, gs_cfg, stopt_cfg,
            single_opt_cfg,
            bond_cfg, search_cfg, refine_mode_kind,
            out_dir=out_dir,
            ref_pdb_path=ref_pdb_path,
            depth=depth + 1,
            seg_counter=seg_counter,
            branch_tag=f"{branch_tag}B",
            pair_index=pair_index,
            mep_mode_kind=mep_mode_kind, calc_cfg=calc_cfg, dmf_cfg=dmf_cfg,
            kink_seq_count=_trailing_kink_count(seg_reports),
        )
        _tag_images(sub.images, pair_index=pair_index)
        return sub

    stitched_imgs, stitched_E = _stitch_paths(
        parts,
        stitch_rmsd_thresh=float(search_cfg["stitch_rmsd_thresh"]),
        bridge_rmsd_thresh=float(search_cfg["bridge_rmsd_thresh"]),
        shared_calc=shared_calc,
        gs_cfg=gs_bridge_cfg,
        stopt_cfg=stopt_cfg,
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



@click.command(
    help="Multistep MEP search via recursive GSM segmentation.",
    context_settings={
        "help_option_names": ["-h", "--help"],
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
    help=("Two or more PDB/mmCIF structures, or XYZ files with corresponding "
          "--ref-pdb entries, in reaction order. "
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
    help="ML-only, link-H-free PDB subset; atom identity/order must match the "
         "full PDB/parm7. When provided, it defines ML membership; "
         "--detect-layer still reads valid movable/frozen MM B-factors.",
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
    "-q",
    "--charge",
    type=int,
    required=False,
    help="ML region charge. Required unless --ligand-charge is provided.",
)
@click.option("-l", "--ligand-charge", type=str, default=None, show_default=False,
              help="Total charge for unknown ligand residues or a per-resname mapping "
                   "(e.g., GPP:-3,SAM:1), used to derive the ML-region charge when -q "
                   "is omitted (requires PDB input or --ref-pdb).")
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
    help="MEP method: gsm (Growing String) or dmf (Direct Max Flux).",
)
@click.option(
    "--dmf-backend",
    type=click.Choice(["cpu", "gpu"], case_sensitive=False),
    default="gpu",
    show_default=True,
    help="DMF compute backend (--mep-mode dmf only): gpu (dmf.torch / CUDA) or cpu (dmf / NumPy). "
    "On a GPU out-of-memory error, retry with cpu.",
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
@click.option(
    "--max-nodes",
    type=int,
    default=GS_KW["max_nodes"],
    show_default=True,
    help=(
        "Number of movable internal images per GSM or DMF segment "
        "(total images = max_nodes + 2 endpoints); recursive segments may "
        "override it with YAML search.max_nodes_segment."
    ),
)
@click.option(
    "--max-cycles",
    type=int,
    default=300,
    show_default=True,
    help="Maximum MEP optimization cycles.",
)
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
    # Only "grad" (LBFGS) is supported, so YAML/CLI `opt_mode: hess`
    # raises a clean Click error instead of silently running LBFGS.
    type=click.Choice(["grad"], case_sensitive=False),
    default="grad",
    show_default=True,
    help="Single-structure optimizer: grad (=L-BFGS). RFO (hess) not yet wired.",
)
@click.option("-o", "--out-dir", "out_dir", type=str, default=OUT_DIR_PATH_SEARCH, show_default=True, help="Output directory.")
@click.option(
    "--thresh",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default=None,
    help=(
        "Convergence preset for single L-BFGS runs only. "
        "The MEP itself keeps --thresh-gsm / --thresh-dmf."
    ),
)
@click.option(
    "--thresh-gsm",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default=None,
    show_default=False,
    help=(
        "Convergence preset for the GSM string optimizer "
        "(gau_loose|gau|gau_tight|gau_vtight|baker|never). "
        "Defaults to 'gau_loose' when not provided."
    ),
)
@click.option(
    "--thresh-dmf",
    type=str,
    default=None,
    show_default=False,
    help=(
        "IPOPT dual-infeasibility tolerance for the DMF path optimizer: "
        "tight (0.04) | middle (0.10) | loose (0.20) or a positive float. "
        "This is not a Gaussian preset. Defaults to 'tight' when not provided."
    ),
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
    # Default True matches GS_KW.fix_first/fix_last semantics and the
    # mlmm-all.py forwarding: endpoints are typically pre-relaxed before
    # string growth to avoid GSM step inflation.
    default=True,
    show_default=True,
    help="If True, run initial single-structure optimizations of inputs."
)
# Input alignment switch (default True)
@click.option(
    "--align/--no-align",
    "align",
    default=True,
    show_default=True,
    help=("After pre-optimization, align all inputs to the *first* input and match freeze_atoms "
          "using the align_freeze_atoms API.")
)
# Full template PDBs for XYZ→PDB conversion and topology reference
@click.option(
    "--ref-pdb",
    "ref_pdb_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True,
    default=None,
    help=("Full-size template PDBs in the same reaction order as --input. "
          "Required when using XYZ inputs to provide topology and B-factor information.")
)
@click.option(
    "--convert-files/--no-convert-files",
    "convert_files",
    default=True,
    show_default=True,
    help="Convert XYZ/TRJ outputs into PDB companions based on the input format.",
)
@click.option(
    "-b", "--backend",
    type=click.Choice(["uma", "orb", "mace", "aimnet2"], case_sensitive=False),
    default=None,
    show_default=False,
    help="ML backend for the ONIOM high-level region (default: uma).",
)
@click.option(
    "--embedcharge/--no-embedcharge",
    "embedcharge",
    default=False,
    show_default=True,
    help="Unavailable in v0.3.3; retained so older commands fail with an actionable diagnostic.",
)
@click.option(
    "--embedcharge-cutoff",
    "embedcharge_cutoff",
    type=float,
    default=None,
    show_default=False,
    help="Unavailable in v0.3.3 together with the retired electronic-embedding path.",
)
@click.option(
    "--link-atom-method",
    "link_atom_method",
    type=click.Choice(["scaled", "fixed"], case_sensitive=False),
    default=None,
    show_default=False,
    help="Link-atom position mode: scaled (g-factor, default) or fixed (legacy 1.09/1.01 Å).",
)
@click.option(
    "--mm-backend",
    "mm_backend",
    type=click.Choice(["hessian_ff", "openmm"], case_sensitive=False),
    default=None,
    show_default=False,
    help="MM backend (default: hessian_ff). MM Hessians use finite differences by default; set calc.mm_fd: false for the hessian_ff analytical path.",
)
@click.option(
    "--cmap/--no-cmap",
    "use_cmap",
    default=None,
    show_default=False,
    help="Preserve CMAP terms in both real and model MM layers. Default: enabled when present in parm7.",
)
@add_ml_layer_detection_options()
@add_precision_option()
@add_workers_options()
@add_backend_model_option()
@add_calc_file_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
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
    ligand_charge: Optional[str],
    spin: Optional[int],
    mep_mode: str,
    dmf_backend: str,
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
    thresh_gsm: Optional[str],
    thresh_dmf: Optional[str],
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    pre_opt: bool,
    align: bool,
    ref_pdb_paths: Optional[Sequence[Path]],
    convert_files: bool,
    backend: Optional[str],
    embedcharge: bool,
    embedcharge_cutoff: Optional[float],
    link_atom_method: Optional[str],
    mm_backend: Optional[str],
    use_cmap: Optional[bool],
    precision: Optional[str],
    workers: Optional[int],
    workers_per_node: Optional[int],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: Optional[str],
) -> None:
    from mlmm.core.utils import (
        collect_option_values,
        current_cli_args,
        reject_option_like_extra_args,
    )

    argv_all = current_cli_args(ctx)
    _claimed_values = collect_option_values(
        argv_all, ("-i", "--input", "--ref-pdb")
    )
    reject_option_like_extra_args(
        ctx.args,
        allowed_values=_claimed_values,
        consumed_values=[*input_paths, *(ref_pdb_paths or ())],
    )
    set_convert_file_enabled(convert_files)
    prepared_inputs: List[PreparedInputStructure] = []
    # --- Robustly accept both styles for -i/--input and --ref-pdb ---
    i_vals = collect_option_values(argv_all, ("-i", "--input"))
    if i_vals:
        i_parsed = validate_existing_files(
            i_vals,
            option_name="-i/--input",
            hint="When using '-i', list only existing file paths (multiple paths may follow a single '-i').",
        )
        input_paths = tuple(i_parsed)

    ref_vals = collect_option_values(argv_all, ("--ref-pdb",))
    if ref_vals:
        ref_parsed = validate_existing_files(
            ref_vals,
            option_name="--ref-pdb",
            hint="When using '--ref-pdb', multiple files may follow a single option.",
        )
        ref_pdb_paths = tuple(ref_parsed)
    # --- end of robust parsing fix ---

    _is_param_explicit = make_is_param_explicit(ctx)

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
    command_str = "mlmm " + " ".join(argv_all)
    try:
        if len(input_paths) < 2:
            raise click.BadParameter("Provide at least two structures for --input in reaction order (reactant [intermediates ...] product).")

        p_list = [Path(p) for p in input_paths]
        ref_list = list(ref_pdb_paths) if ref_pdb_paths else []
        prepared_inputs = []
        for i, p in enumerate(p_list):
            pi = prepare_input_structure(p)
            if p.suffix.lower() == ".xyz":
                if i < len(ref_list):
                    apply_ref_pdb_override(pi, ref_list[i])
                else:
                    raise click.BadParameter(
                        f"XYZ input '{p}' requires a corresponding --ref-pdb for topology/B-factor info."
                    )
            elif p.suffix.lower() not in {".pdb", ".cif", ".mmcif"}:
                raise click.BadParameter(
                    f"'{p}': unsupported format. Use .pdb/.cif/.mmcif or .xyz (with --ref-pdb)."
                )
            prepared_inputs.append(pi)
        config_layer_cfg = load_yaml_dict(config_yaml)
        override_layer_cfg = load_yaml_dict(override_yaml)

        mep_mode_kind = mep_mode.lower().strip()
        refine_mode_kind = refine_mode.strip().lower() if refine_mode else None

        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        gs_cfg = dict(GS_KW)
        stopt_cfg = dict(STOPT_KW)
        lbfgs_cfg = dict(LBFGS_KW)
        bond_cfg = dict(BOND_KW)
        search_cfg = dict(SEARCH_KW)
        dmf_cfg = fresh_dmf_config()

        apply_yaml_overrides(
            config_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",), ("mlmm",))),
                (gs_cfg, (("gs",),)),
                (stopt_cfg, (("stopt",), ("opt",))),
                (lbfgs_cfg, (("stopt", "lbfgs"), ("lbfgs",))),
                (bond_cfg, (("bond",),)),
                (search_cfg, (("search",),)),
                (dmf_cfg, (("dmf",),)),
            ],
        )

        # CLI explicit overrides (after config YAML, before override YAML)
        if backend is not None:
            calc_cfg["backend"] = str(backend).lower()
        from mlmm.backends import apply_precision_to_calc_cfg
        # Always run so a YAML-set workers>1 also gets the analytical-Hessian guard.
        from mlmm.backends import apply_workers_to_calc_cfg
        apply_workers_to_calc_cfg(calc_cfg, workers, workers_per_node)
        from mlmm.backends import apply_backend_model_to_calc_cfg
        # Unconditional: also pops a raw backend_model token from a --config YAML.
        apply_backend_model_to_calc_cfg(calc_cfg, backend_model)
        # --calc-file overrides --backend with a user ASE Calculator (custom backend).
        from mlmm.backends import apply_calc_file_to_calc_cfg
        apply_calc_file_to_calc_cfg(calc_cfg, calc_file, calc_factory)
        # Unconditional: also dispatches a --config YAML calc.precision
        # after the final backend has been selected.
        apply_precision_to_calc_cfg(calc_cfg, precision)
        geom_cfg["coord_type"] = "cart"  # no microiteration: DLC over MM atoms is meaningless; fixed to cart
        if _is_param_explicit("embedcharge"):
            calc_cfg["embedcharge"] = bool(embedcharge)
        if _is_param_explicit("embedcharge_cutoff"):
            calc_cfg["embedcharge_cutoff"] = embedcharge_cutoff
        if link_atom_method is not None:
            calc_cfg["link_atom_method"] = str(link_atom_method).lower()
        if mm_backend is not None:
            calc_cfg["mm_backend"] = str(mm_backend).lower()
        if use_cmap is not None:
            calc_cfg["use_cmap"] = use_cmap

        try:
            geom_freeze = _normalize_geom_freeze(geom_cfg.get("freeze_atoms"))
        except click.BadParameter as e:
            click.echo(f"ERROR: {e}", err=True)
            sys.exit(1)
        geom_cfg["freeze_atoms"] = geom_freeze
        _convert_yaml_layer_atoms_1to0(calc_cfg)

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
                ligand_charge=ligand_charge,
                prefix="[path-search]",
                model_pdb=model_pdb,
                model_indices_spec=model_indices_str,
                detect_layer=detect_layer,
                yaml_cfg=merged_yaml_cfg,
            )
        # CLI-resolved charge/spin (from -q / -l derivation, or -m / spin_default)
        # always wins over the CALC_KW default carried in calc_cfg.
        calc_cfg["model_charge"] = int(resolved_charge)
        calc_cfg["model_mult"] = int(resolved_spin)

        # The calculator consumes the normalized internal PDB topology.
        calc_cfg["input_pdb"] = str(prepared_inputs[0].source_path)
        calc_cfg["real_parm7"] = str(real_parm7)

        detect_layer_effective = bool(calc_cfg.get("use_bfactor_layers", detect_layer))
        if _is_param_explicit("detect_layer"):
            detect_layer_effective = bool(detect_layer)
            calc_cfg["use_bfactor_layers"] = detect_layer_effective

        if _is_param_explicit("max_nodes"):
            gs_cfg["max_nodes"] = int(max_nodes)
            search_cfg["max_nodes_segment"] = int(max_nodes)
        if _is_param_explicit("max_cycles"):
            stopt_cfg["max_cycles"] = int(max_cycles)
            stopt_cfg["stop_in_when_full"] = int(max_cycles)
            dmf_cfg["max_cycles"] = int(max_cycles)
        if _is_param_explicit("dmf_backend"):
            dmf_cfg["backend"] = str(dmf_backend).lower()
        if _is_param_explicit("climb"):
            gs_cfg["climb"] = bool(climb)
            gs_cfg["climb_lanczos"] = bool(climb)
        if _is_param_explicit("dump"):
            stopt_cfg["dump"] = bool(dump)
            lbfgs_cfg["dump"] = bool(dump)
        if _is_param_explicit("out_dir"):
            stopt_cfg["out_dir"] = out_dir
            lbfgs_cfg["out_dir"] = out_dir
        if _is_param_explicit("thresh") and thresh is not None:
            lbfgs_cfg["thresh"] = str(thresh)
        if _is_param_explicit("thresh_gsm") and thresh_gsm is not None:
            stopt_cfg["thresh"] = str(thresh_gsm)
        if _is_param_explicit("thresh_dmf") and thresh_dmf is not None:
            dmf_cfg["tol"] = str(thresh_dmf)
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
                (stopt_cfg, (("stopt",), ("opt",))),
                (lbfgs_cfg, (("stopt", "lbfgs"), ("lbfgs",))),
                (bond_cfg, (("bond",),)),
                (search_cfg, (("search",),)),
                (dmf_cfg, (("dmf",),)),
            ],
        )
        # The final layer may replace strict method enums or the workers count.
        # Revalidate the fully resolved calculator mapping before dry-run can
        # report success (the constructor repeats this for normal execution).
        apply_workers_to_calc_cfg(calc_cfg, None, None)

        # A dormant YAML DMF section does not affect GSM. An explicit CLI
        # tolerance is still validated as user input, regardless of MEP mode.
        if mep_mode_kind == "dmf" or _is_param_explicit("thresh_dmf"):
            resolve_dmf_solve_tol(dmf_cfg, prefix="[path-search]")

        refine_mode_kind = search_cfg.get("refine_mode")
        if refine_mode_kind is None:
            refine_mode_kind = "peak" if mep_mode_kind == "gsm" else "minima"
        else:
            refine_mode_kind = str(refine_mode_kind).strip().lower()
            if refine_mode_kind not in {"peak", "minima"}:
                raise click.BadParameter(f"Unknown --refine-mode '{refine_mode_kind}'.")
        search_cfg["refine_mode"] = refine_mode_kind

        out_dir_path = Path(stopt_cfg.get("out_dir", out_dir)).resolve()
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

        # For layer detection, prefer --ref-pdb (which carries B-factor layers)
        # over the first input (which may be XYZ).
        if ref_list and ref_list[0]:
            layer_source_pdb = Path(ref_list[0]).resolve()
        else:
            layer_source_pdb = prepared_inputs[0].source_path.resolve()
        if detect_layer_effective and layer_source_pdb.suffix.lower() != ".pdb":
            click.echo("ERROR: --detect-layer requires a PDB input (or --ref-pdb).", err=True)
            sys.exit(1)

        from mlmm.core.embedcharge_policy import reject_retired_embedcharge_cli

        reject_retired_embedcharge_cli(
            calc_cfg,
            cutoff_requested=_is_param_explicit("embedcharge_cutoff"),
        )

        if dry_run:
            if model_pdb_effective is not None:
                model_region_source = "model_pdb"
            elif model_indices:
                model_region_source = "model_indices"
            else:
                model_region_source = "bfactor"

            validation_cfg = dict(calc_cfg)
            try:
                with tempfile.TemporaryDirectory(prefix="mlmm_path_search_validate_") as tmp:
                    _, layer_info_preview = resolve_ml_layer_assignment(
                        source_path=layer_source_pdb,
                        out_dir_path=Path(tmp),
                        model_pdb=model_pdb_effective,
                        model_indices=model_indices,
                        detect_layer=detect_layer_effective,
                        hess_cutoff=hess_cutoff_effective,
                        movable_cutoff=movable_cutoff_effective,
                        calc_cfg=validation_cfg,
                        protected_inputs=(),
                        echo_fn=click.echo,
                    )
            except click.ClickException as exc:
                click.echo(f"ERROR: {exc.message}", err=True)
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
                    force=True)
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
                "backend": calc_cfg.get("backend", "uma"),
                "embedcharge": bool(calc_cfg.get("embedcharge", False)),
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
            emit(
                format_elapsed("[time] Elapsed Time for Path Search", time_start),
                narrative=True,
            )
            return

        try:
            model_pdb_path, layer_info = resolve_ml_layer_assignment(
                source_path=layer_source_pdb,
                out_dir_path=out_dir_path,
                model_pdb=model_pdb_effective,
                model_indices=model_indices,
                detect_layer=detect_layer_effective,
                hess_cutoff=hess_cutoff_effective,
                movable_cutoff=movable_cutoff_effective,
                calc_cfg=calc_cfg,
                protected_inputs=(
                    *p_list,
                    *(prepared.source_path for prepared in prepared_inputs),
                    *ref_list,
                    real_parm7,
                    config_yaml,
                    override_yaml,
                    (
                        Path(calc_cfg["calc_file"])
                        if calc_cfg.get("calc_file")
                        else None
                    ),
                ),
                echo_fn=click.echo,
            )
        except click.ClickException as exc:
            click.echo(f"ERROR: {exc.message}", err=True)
            sys.exit(1)
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

        stopt_cfg["stop_in_when_full"] = int(stopt_cfg.get("max_cycles", STOPT_KW["max_cycles"]))
        out_dir_path = Path(stopt_cfg.get("out_dir", out_dir)).resolve()
        echo_geom = format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")
        echo_calc = format_freeze_atoms_for_echo(filter_calc_for_echo(calc_cfg), key="freeze_atoms")
        echo_gs   = strip_inherited_keys(gs_cfg, GS_KW, mode="same")
        echo_stopt = strip_inherited_keys({**stopt_cfg, "out_dir": str(out_dir_path)}, STOPT_KW, mode="same")
        echo_lbfgs = strip_inherited_keys(lbfgs_cfg, LBFGS_KW, mode="same")
        echo_bond = strip_inherited_keys(bond_cfg, BOND_KW, mode="same")
        echo_search = strip_inherited_keys(search_cfg, SEARCH_KW, mode="same")

        click.echo(pretty_block("geom", echo_geom))
        click.echo(pretty_block("calc", echo_calc))
        click.echo(pretty_block("gs",   echo_gs))
        click.echo(pretty_block("stopt", echo_stopt))
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
                force=True)
            )

        effective_max_cycles = (
            dmf_cfg.get("max_cycles", 0)
            if mep_mode_kind == "dmf"
            else stopt_cfg.get("max_cycles", 0)
        )
        if int(effective_max_cycles) <= 0:
            raise click.BadParameter(
                "--max-cycles must be at least 1.",
                param_hint="--max-cycles",
            )

        validate_endpoint_atom_identities(prepared_inputs)
        out_dir_path.mkdir(parents=True, exist_ok=True)

        geoms = _load_structures(
            inputs=prepared_inputs,
            coord_type=geom_cfg.get("coord_type", "cart"),
            base_freeze=geom_cfg.get("freeze_atoms", []),
        )

        shared_calc = mlmm(**calc_cfg)
        for g in geoms:
            g.set_calculator(shared_calc)

        # Reference PDB for output conversion: prefer --ref-pdb, fall back to input PDBs
        ref_pdb_for_segments: Optional[Path] = None
        if prepared_inputs:
            ref_pdb_for_segments = prepared_inputs[0].source_path.resolve()

        if pre_opt:
            new_geoms: List[Any] = []
            for i, g in enumerate(geoms):
                tag = f"init{i:02d}"
                g_opt, _ = _optimize_single(g, shared_calc, lbfgs_cfg, out_dir_path, tag=tag, ref_pdb_path=ref_pdb_for_segments)
                new_geoms.append(g_opt)
            geoms = new_geoms
        else:
            click.echo("[init] Skipping endpoint pre-optimization as requested by --no-preopt.")

        # Align all inputs to the first structure, guided by freeze constraints, when requested
        align_thresh = str(stopt_cfg.get("thresh", "gau"))
        if align:
            try:
                emit("\n====== Aligning all inputs to the first structure (freeze-guided scan + relaxation) ======\n", narrative=True)
                alignment_results = align_and_refine_sequence_inplace(
                    geoms,
                    thresh=align_thresh,
                    shared_calc=shared_calc,
                    out_dir=out_dir_path / "align_refine",
                    verbose=True,
                )
                failed_pairs = alignment_failed_pair_indices(alignment_results)
                if failed_pairs:
                    raise click.ClickException(
                        "Input alignment did not converge for pair(s): "
                        + ", ".join(str(index) for index in failed_pairs)
                    )
                click.echo("[align] Completed input alignment.")
            except Exception as e:
                raise click.ClickException(f"Input alignment failed: {e}") from e
        else:
            click.echo("[align] Skipping input alignment as requested by --no-align.")

        _mep_search_start = time.perf_counter()
        emit("\n====== Multistep MEP search (multi-structure) started ======\n", narrative=True)
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
                geom_cfg, gs_cfg, stopt_cfg,
                lbfgs_cfg,
                bond_cfg, search_cfg, refine_mode_kind,
                out_dir=out_dir_path,
                ref_pdb_path=ref_pdb_for_segments,
                depth=0,
                seg_counter=seg_counter,
                branch_tag="B",
                pair_index=None,
                mep_mode_kind=mep_mode_kind, calc_cfg=calc_cfg, dmf_cfg=dmf_cfg,
                kink_seq_count=_trailing_kink_count(seg_reports_all),
            )
            return sub

        for i in range(len(geoms) - 1):
            gA, gB = geoms[i], geoms[i + 1]
            pair_tag = f"pair_{i:02d}"
            emit(f"\n--- Processing pair {i:02d}: image {i} → {i+1} ---", narrative=True)
            pair_path = _build_multistep_path(
                gA, gB,
                shared_calc,
                geom_cfg, gs_cfg, stopt_cfg,
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
                    stopt_cfg=stopt_cfg,
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
            emit(
                f"[stage] Pair {i:02d} done: images={len(pair_path.images)}, "
                f"segments={len(pair_path.segments)}",
                detail=True,
            )

        emit(
            "====== Multistep MEP search (multi-structure) finished "
            f"(pairs={max(len(geoms) - 1, 0)}, segments={len(seg_reports_all)}, "
            f"elapsed={time.perf_counter() - _mep_search_start:.1f}s) ======\n",
            narrative=True,
        )

        combined_all = CombinedPath(images=combined_imgs, energies=combined_Es, segments=seg_reports_all)

        for idx, srep in enumerate(combined_all.segments, 1):
            srep.seg_index = idx
        tag_to_index = {s.tag: int(s.seg_index) for s in combined_all.segments}
        for im in combined_all.images:
            tag = getattr(im, "mep_seg_tag", None)
            if tag and tag in tag_to_index:
                try:
                    setattr(im, "mep_seg_index", int(tag_to_index[tag]))
                except Exception:
                    logger.debug("Failed to set mep_seg_index on image", exc_info=True)

        # Always write mep_trj.xyz for downstream compatibility; convert to PDB when possible.
        pdb_input = ref_pdb_for_segments is not None
        final_trj = out_dir_path / "mep_trj.xyz"
        _write_xyz_trj_with_energy(combined_all.images, combined_all.energies, final_trj)
        emit(f"[write] Wrote '{final_trj}'.", detail=True)
        try:
            run_trj2fig(final_trj, [out_dir_path / "mep_plot.png"], unit="kcal", reference="init", reverse_x=False)
            emit(f"[plot] Saved energy plot → '{out_dir_path / 'mep_plot.png'}'", detail=True)
        except Exception as e:
            click.echo(f"[plot] WARNING: Failed to plot final energy: {e}", err=True)

        if pdb_input and is_convert_file_enabled():
            try:
                final_pdb = out_dir_path / "mep.pdb"
                convert_xyz_to_pdb(final_trj, ref_pdb_for_segments, final_pdb)
                emit(f"[convert] Wrote '{final_pdb}'.", detail=True)
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
                    emit(f"[write] Wrote per-segment pocket trajectory → '{seg_trj}'", detail=True)
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
                    emit(f"[write] Wrote segment HEI (pocket) → '{hei_trj}'", detail=True)
                    if ref_pdb_for_segments is not None:
                        _maybe_convert_to_pdb(hei_trj, ref_pdb_for_segments, out_path=out_dir_path / f"hei_seg_{seg_idx:02d}.pdb")
        except Exception as e:
            click.echo(f"[write] WARNING: Failed to emit per-segment pocket outputs: {e}", err=True)
        # ---- END ----

        frame_ranges = _frame_ranges_by_segment(combined_all.images)
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
                    "bond_changes": (s.summary if (s.kind != "bridge") else ""),
                    # the segment's reported optimizer convergence, so the
                    # all-pipeline aggregate can gate on it (a nonconverged segment
                    # keeps its trajectory but cannot make the path a success).
                    "converged": s.converged,
                    **frame_ranges.get(int(s.seg_index), {}),
                } for s in combined_all.segments
            ],
        }

        try:
            overall_changed, overall_summary = _has_bond_change(combined_all.images[0], combined_all.images[-1], bond_cfg)
        except Exception:
            logger.debug(
                "path_search: overall bond-change diff failed; reporting no covalent changes",
                exc_info=True,
            )
            overall_changed, overall_summary = False, ""

        emit("\n====== MEP Summary started ======\n", narrative=True)

        emit("\n[overall] Covalent-bond changes between first and last image:", narrative=True)
        if overall_changed and overall_summary.strip():
            click.echo(textwrap.indent(overall_summary.strip(), prefix="  "))
        else:
            click.echo("  (no covalent changes detected)")

        if combined_all.segments:
            emit("\n[segments] Along the final MEP order (ΔE‡, ΔE). Bridges are shown between connected segments:", narrative=True)
            for i, seg in enumerate(combined_all.segments, 1):
                kind_label = "BRIDGE" if seg.kind == "bridge" else "SEG"
                click.echo(f"  [{i:02d}] ({kind_label}) {seg.tag}  |  ΔE‡ = {seg.barrier_kcal:.2f} kcal/mol,  ΔE = {seg.delta_kcal:.2f} kcal/mol")
                if seg.kind != "bridge" and seg.summary.strip():
                    click.echo(textwrap.indent(seg.summary.strip(), prefix="      "))
        else:
            click.echo("\n[segments] (no segment reports)")

        emit("====== MEP Summary finished ======\n", narrative=True)

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
                emit(f"[diagram] Wrote energy diagram (PNG) → '{png_path}'", detail=True)
            except Exception as e:
                click.echo(f"[diagram] NOTE: PNG export skipped (install 'kaleido' to enable): {e}", err=True)

            chain_text = " ".join(chain_tokens)
            emit(f"[diagram] State label sequence: {chain_text}", detail=True)

        except Exception as e:
            click.echo(f"[diagram] WARNING: Failed to build energy diagram: {e}", err=True)

        if diagram_payload is not None:
            summary["energy_diagrams"] = [diagram_payload]

        _enrich_path_summary_contract(
            summary,
            segments=combined_all.segments,
            out_dir=out_dir_path,
            calc_cfg=calc_cfg,
            command=command_str,
        )
        summary["references"] = method_references(
            {
                "pipeline_mode": "path-search",
                "mep_mode": mep_mode_kind,
                "path_opt_mode": opt_mode,
                "dmf_correlated": bool(dmf_cfg.get("correlated", False)),
            }
        )

        summary = with_current_run_id(summary)
        commit_json_exact(out_dir_path / "summary.json", summary)
        emit(f"[write] Wrote '{out_dir_path / 'summary.json'}'.", detail=True)

        summary_payload_for_citations: Dict[str, Any] = {}
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
                "mep_mode": mep_mode_kind,
                "dmf_correlated": bool(dmf_cfg.get("correlated", False)),
                **_summary_log_provenance(summary),
                "command": command_str,
                "charge": calc_cfg.get("model_charge"),
                "spin": calc_cfg.get("model_mult"),
                "freeze_atoms": freeze_atoms_for_log,
                "mep": mep_info,
                "segments": summary.get("segments", []),
                "energy_diagrams": summary.get("energy_diagrams", []),
                "key_files": {},
            }
            summary_payload_for_citations = summary_payload
            write_summary_log(out_dir_path / "summary.log", summary_payload)
            emit(f"[write] Wrote '{out_dir_path / 'summary.log'}'.", detail=True)
        except Exception as e:
            click.echo(f"[write] WARNING: Failed to write summary.log: {e}", err=True)

        from mlmm.core.utils import is_child_mode

        if summary_payload_for_citations and not is_child_mode():
            emit_method_citations(summary_payload_for_citations)
        emit(
            format_elapsed("[time] Elapsed Time for Path Search", time_start),
            narrative=True,
        )

    except ZeroStepLength as e:
        _write_error_json(Path(out_dir).resolve(), "path-search", e, "ZeroStepLength", time_start)
        click.echo("ERROR: Proposed step length dropped below the minimum allowed (ZeroStepLength).", err=True)
        sys.exit(2)
    except OptimizationError as e:
        _write_error_json(Path(out_dir).resolve(), "path-search", e, "OptimizationError", time_start)
        click.echo(f"ERROR: Path search failed — {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        render_cli_exception(e, label="path search", out_dir=out_dir, command="path-search", time_start=time_start)
    finally:
        for prepared in prepared_inputs:
            prepared.cleanup()
        # Release GPU memory so subsequent pipeline stages don't OOM.
        # `= None` decref's the heavy refs; `del` then removes names from
        # the local frame so torch.nn.Module hooks / closures cannot retain.
        shared_calc = geoms = None
        del shared_calc, geoms
        gc.collect()  # break cyclic refs inside torch.nn.Module
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
