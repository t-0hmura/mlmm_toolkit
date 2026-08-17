"""
ML/MM minimum-energy path optimization via the Growing String Method (GSM) or Direct Max Flux (DMF).

Example:
    mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0
    mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 -q 0 --mep-mode dmf

For detailed documentation, see: docs/path_opt.md
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

import gc
import inspect
import logging
import math
import sys
import traceback
import textwrap

logger = logging.getLogger(__name__)

import click
from mlmm.core.output import emit
import numpy as np
import time
import torch

from pysisyphus.helpers import geom_loader
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError
from pysisyphus.optimizers.LBFGS import LBFGS

from mlmm.backends.mlmm_calc import mlmm, MLMMASECalculator
from mlmm.workflows.opt import (
    GEOM_KW as OPT_GEOM_KW,
    CALC_KW as OPT_CALC_KW,
    LBFGS_KW as OPT_LBFGS_KW,
    _parse_freeze_atoms as _parse_freeze_atoms_opt,
    _normalize_geom_freeze as _normalize_geom_freeze_opt,
)
from mlmm.workflows.opt import _convert_yaml_layer_atoms_1to0
from mlmm.workflows.charge_prep import resolve_charge_spin_or_raise
from mlmm.core.utils import (
    apply_layer_freeze_constraints,
    convert_xyz_to_pdb,
    set_convert_file_enabled,
    is_convert_file_enabled,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    strip_inherited_keys,
    filter_calc_for_echo,
    format_freeze_atoms_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    apply_ref_pdb_override,
    prepare_input_structure,
    PreparedInputStructure,
    validate_endpoint_atom_identities,
    parse_indices_string,
    resolve_ml_layer_assignment,
    echo_resolved_device,
)
from mlmm.cli.common_options import add_ml_layer_detection_options, add_precision_option, add_workers_options, add_backend_model_option, add_calc_file_option, add_deterministic_option, add_allow_charge_mult_mismatch_option
from mlmm.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, make_is_param_explicit, _write_error_json, render_cli_exception
from mlmm.workflows.align_freeze import (
    align_and_refine_sequence_inplace,
    alignment_failed_pair_indices,
)
from mlmm.core.defaults import (
    BFACTOR_FROZEN,
    BFACTOR_ML,
    BFACTOR_MOVABLE_MM,
    DMF_KW as _DMF_KW_DEFAULT,
    fresh_dmf_config,
    GS_KW as _GS_KW_DEFAULT,
    OUT_DIR_PATH_OPT,
    STOPT_KW as _STOPT_KW_DEFAULT,
    THRESH_CHOICES,
)


# Defaults (overridden by YAML/CLI)

# Geometry (input handling) — share defaults with opt.py
GEOM_KW: Dict[str, Any] = deepcopy(OPT_GEOM_KW)

# ML/MM calculator settings — share defaults with opt.py
CALC_KW: Dict[str, Any] = deepcopy(OPT_CALC_KW)

# LBFGS (used for optional endpoint pre-optimization)
LBFGS_KW: Dict[str, Any] = deepcopy(OPT_LBFGS_KW)

# DMF (Direct Max Flux) defaults
DMF_KW: Dict[str, Any] = deepcopy(_DMF_KW_DEFAULT)

# GrowingString (path representation)
GS_KW: Dict[str, Any] = deepcopy(_GS_KW_DEFAULT)

# StringOptimizer (optimization control)
STOPT_KW: Dict[str, Any] = deepcopy(_STOPT_KW_DEFAULT)

def _load_two_endpoints(
    inputs: Sequence[PreparedInputStructure],
    coord_type: str,
    base_freeze: Sequence[int],
) -> Sequence:
    """
    Load the two endpoint structures and set `freeze_atoms` as needed.
    """
    geoms = []
    for prepared in inputs:
        geom_path = prepared.geom_path
        cfg: Dict[str, Any] = {"freeze_atoms": list(base_freeze)}
        freeze = merge_freeze_atom_indices(cfg)
        g = geom_loader(geom_path, coord_type=coord_type, freeze_atoms=freeze)
        g.freeze_atoms = np.array(freeze, dtype=int)
        geoms.append(g)
    return geoms


# Helpers shared with opt.py (imported for consistency)
_parse_freeze_atoms = _parse_freeze_atoms_opt
_normalize_geom_freeze = _normalize_geom_freeze_opt


# B-factor annotation helpers

def _parse_pdb_atoms_for_indexing(pdb_path: Path) -> List[Dict[str, Any]]:
    """Parse PDB ATOM/HETATM records for indexing/matching."""
    atoms: List[Dict[str, Any]] = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                # Ensure line is long enough
                s = line.rstrip("\n")
                s = s + (" " * (80 - len(s))) if len(s) < 80 else s
                serial_str = s[6:11].strip()
                resseq_str = s[22:26].strip()
                try:
                    serial = int(serial_str) if serial_str else None
                except ValueError:
                    serial = None
                try:
                    resseq = int(resseq_str) if resseq_str else None
                except ValueError:
                    resseq = None
                atom = {
                    "line": s,
                    "serial": serial,
                    "name": s[12:16].strip(),
                    "altloc": s[16].strip() if len(s) > 16 else "",
                    "resname": s[17:20].strip(),
                    "chain": s[21].strip() if len(s) > 21 else "",
                    "resseq": resseq,
                    "icode": s[26].strip() if len(s) > 26 else "",
                }
                atoms.append(atom)
    return atoms


def _compute_ml_indices_from_model_and_ref(ref_pdb: Path, model_pdb: Path) -> Set[int]:
    """
    Compute 0-based atom indices (in the order of ATOM/HETATM records of ref_pdb)
    that belong to the ML region defined by model_pdb.

    Matching strategy (robust, fallback-based):
      1) Exact key match on (name, altloc, resname, chain, resseq, icode)
      2) Key match ignoring altloc when model altloc is blank
      3) Serial number match
    """
    ref_atoms = _parse_pdb_atoms_for_indexing(ref_pdb)
    model_atoms = _parse_pdb_atoms_for_indexing(model_pdb)

    # Build maps for the reference structure
    key_to_indices: Dict[Tuple[str, str, str, str, Optional[int], str], List[int]] = {}
    key_wo_alt_to_indices: Dict[Tuple[str, str, str, Optional[int], str], List[int]] = {}
    serial_to_index: Dict[int, int] = {}

    for idx, a in enumerate(ref_atoms):
        key = (a["name"], a["altloc"], a["resname"], a["chain"], a["resseq"], a["icode"])
        key_wo = (a["name"], a["resname"], a["chain"], a["resseq"], a["icode"])
        key_to_indices.setdefault(key, []).append(idx)
        key_wo_alt_to_indices.setdefault(key_wo, []).append(idx)
        if a["serial"] is not None:
            # If duplicated serials exist, keep the first occurrence
            serial_to_index.setdefault(a["serial"], idx)

    ml_indices: Set[int] = set()
    misses = 0

    for ma in model_atoms:
        key = (ma["name"], ma["altloc"], ma["resname"], ma["chain"], ma["resseq"], ma["icode"])
        key_wo = (ma["name"], ma["resname"], ma["chain"], ma["resseq"], ma["icode"])

        idx: Optional[int] = None

        candidates = key_to_indices.get(key)
        if candidates and len(candidates) == 1:
            idx = candidates[0]
        elif candidates and len(candidates) > 1:
            # Multiple; try to disambiguate via serial number if possible
            if ma["serial"] is not None:
                si = serial_to_index.get(ma["serial"])
                if si in candidates:
                    idx = si

        if idx is None and (ma["altloc"] == "" or ma["altloc"] == " "):
            candidates2 = key_wo_alt_to_indices.get(key_wo)
            if candidates2:
                idx = candidates2[0]  # pick first

        if idx is None and ma["serial"] is not None:
            idx = serial_to_index.get(ma["serial"])

        if idx is not None:
            ml_indices.add(idx)
        else:
            misses += 1

    if misses:
        click.echo(f"[annotate] WARNING: {misses} ML atoms from '{model_pdb.name}' could not be mapped to '{ref_pdb.name}'.", err=True)

    return ml_indices


def _apply_bfactor_annotations_inplace(
    pdb_path: Path,
    ml_indices: Set[int],
    freeze_indices: Sequence[int],
    beta_ml: float = BFACTOR_ML,
    beta_freeze: float = BFACTOR_FROZEN,
    beta_both: float = BFACTOR_ML,
) -> None:
    """
    In-place update of B-factors for PDB ATOM/HETATM records.

    Rules (new 3-layer encoding):
      - ML only:       0  (BFACTOR_ML)
      - Freeze only:  20  (BFACTOR_FROZEN)
      - ML ∩ Freeze:   0  (ML takes precedence)
      - Others:       10  (BFACTOR_MOVABLE_MM)

    The index for lookups is the 0-based position among ATOM/HETATM
    records and resets at each MODEL record for multi-model PDBs.
    """
    freeze_set: Set[int] = set(int(i) for i in (freeze_indices or []))
    ml_set: Set[int] = set(int(i) for i in (ml_indices or set()))

    def _format_b(b: float) -> str:
        # PDB tempFactor field is 6.2 width
        return f"{b:6.2f}"

    # Read and process lines
    lines_out: List[str] = []
    atom_idx = 0  # 0-based within a MODEL (or entire file if no MODEL)

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith("MODEL"):
            # Reset index at each model
            atom_idx = 0
            lines_out.append(line)
            continue

        if line.startswith(("ATOM  ", "HETATM")):
            s = line.rstrip("\n")
            # Pad to at least 66 chars so we can safely replace tempFactor (cols 61-66)
            if len(s) < 66:
                s = s + (" " * (66 - len(s)))

            # Decide B-factor for this atom index
            if (atom_idx in ml_set) and (atom_idx in freeze_set):
                b = beta_both
            elif atom_idx in ml_set:
                b = beta_ml
            elif atom_idx in freeze_set:
                b = beta_freeze
            else:
                b = BFACTOR_MOVABLE_MM

            s = s[:60] + _format_b(b) + s[66:]
            # Ensure trailing newline
            s = s if s.endswith("\n") else s + "\n"
            lines_out.append(s)

            atom_idx += 1
        else:
            lines_out.append(line)

    with open(pdb_path, "w") as f:
        f.writelines(lines_out)

    click.echo(
        f"[annotate] Updated B-factors in '{pdb_path}' "
        f"(ML={BFACTOR_ML:.0f}, MovableMM={BFACTOR_MOVABLE_MM:.0f}, "
        f"FrozenMM={BFACTOR_FROZEN:.0f}; {len(ml_set)} ML, {len(freeze_set)} frozen)."
    )



def _select_hei_index(energies: Sequence[float]) -> int:
    """Return the global highest-energy-image index."""
    E = np.array(energies, dtype=float)
    if E.size == 0:
        raise ValueError("Cannot select an HEI from an empty energy profile.")
    if not np.all(np.isfinite(E)):
        raise ValueError("Cannot select an HEI from non-finite energies.")
    return int(np.argmax(E))


@dataclass(frozen=True)
class DMFMepResult:
    """Scientific state returned by one DMF solve."""

    images: Tuple[Any, ...]
    energies: Tuple[float, ...]
    hei_idx: int
    converged: bool
    ipopt_status: Optional[int]
    reason: str


DMF_TOL_PRESETS = ("tight", "middle", "loose")


def resolve_dmf_solve_tol(
    dmf_cfg: Mapping[str, Any], prefix: str = "[path-opt]"
) -> Any:
    """Resolve the ``tol`` argument of the DMF solve from the DMF configuration.

    ``dmf.tol`` accepts the pydmf presets ``tight`` / ``middle`` / ``loose``
    (IPOPT ``dual_inf_tol`` 0.04 / 0.10 / 0.20) or a positive float.  When it is
    unset, an explicitly pinned ``dmf.ipopt_options.dual_inf_tol`` is honoured by
    returning ``None``: pydmf's ``solve`` applies its own ``tol`` after the
    caller's IPOPT options, so a preset passed here would silently replace that
    value.  With neither set, the historical ``tight`` default applies.
    """
    raw = dmf_cfg.get("tol")
    if raw is None:
        pinned = (dmf_cfg.get("ipopt_options") or {}).get("dual_inf_tol")
        return None if pinned is not None else "tight"
    if isinstance(raw, str):
        text = raw.strip().lower()
        if text in DMF_TOL_PRESETS:
            return text
    if isinstance(raw, bool):
        raise click.ClickException(
            f"{prefix} Invalid DMF tolerance '{raw}': expected "
            f"{'|'.join(DMF_TOL_PRESETS)} or a positive float (IPOPT "
            "dual_inf_tol)."
        )
    try:
        value = float(raw)
    except (TypeError, ValueError):
        raise click.ClickException(
            f"{prefix} Invalid DMF tolerance '{raw}': expected "
            f"{'|'.join(DMF_TOL_PRESETS)} or a positive float (IPOPT "
            "dual_inf_tol). Gaussian convergence presets apply to --thresh / "
            "--thresh-gsm, not to the DMF optimizer."
        ) from None
    if not math.isfinite(value) or value <= 0.0:
        raise click.ClickException(
            f"{prefix} Invalid DMF tolerance '{raw}': the IPOPT "
            "dual_inf_tol must be a finite positive number."
        )
    return value


def _dmf_solver_outcome(solve_result: Any) -> Tuple[bool, Optional[int], str]:
    """Normalize cyipopt's ``(x, info)`` result without hiding failures."""

    info: Dict[str, Any] = {}
    if (
        isinstance(solve_result, tuple)
        and len(solve_result) >= 2
        and isinstance(solve_result[1], dict)
    ):
        info = solve_result[1]

    raw_status = info.get("status")
    try:
        status = int(raw_status) if raw_status is not None else None
    except (TypeError, ValueError):
        status = None

    raw_reason = (
        info.get("status_msg")
        or info.get("status_message")
        or info.get("message")
    )
    if isinstance(raw_reason, bytes):
        reason = raw_reason.decode("utf-8", errors="replace")
    elif raw_reason is not None:
        reason = str(raw_reason)
    elif status is None:
        reason = "IPOPT status was not reported."
    else:
        reason = f"IPOPT status {status}."
    return status == 0, status, reason


def _shared_frozen_reference(
    images: Sequence[Any], fix_atoms: Sequence[int]
) -> Optional[np.ndarray]:
    """Copy one validated frozen-coordinate anchor from the first image."""

    indices = tuple(int(index) for index in fix_atoms)
    if not indices:
        return None
    if not images:
        raise ValueError("Frozen-atom restraints require at least one path image.")

    atom_count = len(images[0])
    invalid = sorted({index for index in indices if index < 0 or index >= atom_count})
    if invalid:
        raise ValueError(
            "Frozen atom indices are outside the path image bounds "
            f"[0, {atom_count}): {invalid}"
        )
    if any(len(image) != atom_count for image in images):
        raise ValueError("All DMF path images must contain the same number of atoms.")

    reference = np.asarray(images[0].get_positions(), dtype=float)[list(indices)].copy()
    reference.setflags(write=False)
    return reference


def _build_dmf_result_data(
    result: DMFMepResult,
    calc_cfg: Dict[str, Any],
) -> Dict[str, Any]:
    """Build the JSON payload from structured DMF state only."""

    from mlmm.core.utils import calculator_provenance
    from pysisyphus.constants import AU2KCALPERMOL

    energies = result.energies
    hei_idx = int(result.hei_idx)
    barrier = None
    delta = None
    if energies:
        initial = float(energies[0])
        barrier = (float(energies[hei_idx]) - initial) * AU2KCALPERMOL
        delta = (float(energies[-1]) - initial) * AU2KCALPERMOL

    return {
        "status": "converged" if result.converged else "not_converged",
        "converged": bool(result.converged),
        "mep_mode": "dmf",
        "ipopt_status": result.ipopt_status,
        "reason": result.reason,
        **calculator_provenance(calc_cfg),
        "charge": calc_cfg.get("model_charge"),
        "spin": calc_cfg.get("model_mult"),
        "reactant_energy_hartree": float(energies[0]) if energies else None,
        "product_energy_hartree": float(energies[-1]) if energies else None,
        "image_energies_hartree": [float(energy) for energy in energies],
        "n_images": len(energies),
        "hei_index": hei_idx,
        "hei_energy_hartree": float(energies[hei_idx]) if energies else None,
        "barrier_kcal": round(barrier, 6) if barrier is not None else None,
        "delta_kcal": round(delta, 6) if delta is not None else None,
        "files": {
            "final_geometries_trj_xyz": "final_geometries_trj.xyz",
            "hei_xyz": "hei.xyz",
        },
    }


def _prepare_path_output_dir(path: Path) -> Path:
    """Create the output directory and invalidate prior result envelopes."""

    resolved = Path(path).resolve()
    resolved.mkdir(parents=True, exist_ok=True)
    for name in (
        "result.json",
        "summary.json",
        "final_geometries.pdb",
        "hei.pdb",
        "hei.gjf",
    ):
        (resolved / name).unlink(missing_ok=True)
    return resolved


class _PathOutputCollisionError(click.UsageError):
    """A path-opt output/input collision."""


def _reject_path_output_collisions(
    out_dir: Path,
    protected_inputs: Sequence[Optional[Path]],
) -> None:
    fixed_names = (
        "result.json",
        "summary.json",
        "dmf_initial_trj.xyz",
        "dmf_fbenm_ipopt.out",
        "dmf_ipopt.out",
        "final_geometries_trj.xyz",
        "final_geometries.pdb",
        "hei.xyz",
        "hei.pdb",
        "model_from_bfactor.pdb",
    )
    fixed = {
        (Path(out_dir) / name).resolve(strict=False) for name in fixed_names
    }
    reserved_roots = {
        (Path(out_dir) / name).resolve(strict=False)
        for name in ("preopt", "align_refine")
    }
    for protected in protected_inputs:
        if protected is None:
            continue
        resolved = Path(protected).expanduser().resolve(strict=False)
        if resolved in fixed or any(
            root == resolved or root in resolved.parents
            for root in reserved_roots
        ):
            raise _PathOutputCollisionError(
                f"Input {protected} collides with a reserved path-opt output "
                f"under {out_dir}."
            )


# DMF (Direct Max Flux) MEP optimization

def _release_dmf_interpolation_cache(mxflx_fbenm: Any) -> None:
    """Release an optional PyDMF device cache at the interpolation boundary."""
    release_cache = getattr(mxflx_fbenm, "release_device_cache", None)
    if callable(release_cache):
        release_cache(empty_cache=False)


def _torch_dmf_runtime_kwargs(
    dmf_backend: str,
    dmf_options: Mapping[str, Any],
    fbenm_options: Mapping[str, Any],
    cfbenm_options: Mapping[str, Any],
    *,
    supports_keep_history: bool = True,
) -> Dict[str, Any]:
    """Resolve Torch-only path settings shared by both DMF stages."""
    if dmf_backend != "gpu":
        return {}

    resolved: Dict[str, Any] = {}
    if supports_keep_history:
        resolved["keep_history"] = bool(dmf_options.get("keep_history", False))
    for name in ("device", "dtype"):
        value = dmf_options.get(name)
        value = fbenm_options.get(name, value)
        value = cfbenm_options.get(name, value)
        if value is not None:
            resolved[name] = value
    if "device" not in resolved:
        # The public 'gpu' backend means CUDA. Without an expert device dmf.torch
        # would resolve its own default and silently run on the CPU, so require
        # CUDA here and let the user retry with --dmf-backend cpu.
        if not torch.cuda.is_available():
            raise RuntimeError(
                "--dmf-backend gpu requires a visible CUDA device. Retry with "
                "`--dmf-backend cpu`, or set an explicit expert device in "
                "'dmf.dmf_options.device'."
            )
        resolved["device"] = "cuda"
    return resolved


def _is_cuda_oom(exc: BaseException) -> bool:
    """True if `exc` looks like a CUDA out-of-memory (torch.cuda.OutOfMemoryError or a
    RuntimeError carrying 'out of memory'), so the DMF gpu backend can advise --dmf-backend cpu."""
    if type(exc).__name__ == "OutOfMemoryError":
        return True
    msg = str(exc).lower()
    return "out of memory" in msg or "cuda oom" in msg


def _run_dmf_mep(
    geoms: Sequence,
    shared_calc,
    out_dir_path: Path,
    input_paths: Sequence[Path],
    max_nodes: int,
    fix_atoms: Sequence[int],
    dmf_cfg: Optional[Dict[str, Any]] = None,
    ml_indices_set: Optional[Set[int]] = None,
    freeze_atoms_final: Optional[Sequence[int]] = None,
) -> DMFMepResult:
    """Run Direct Max Flux (DMF) MEP optimization between two endpoints.

    Uses pydmf with harmonic constraints for frozen atoms; the ML/MM ONIOM calculator is
    wrapped as an ASE calculator. The backend is selected by ``dmf_cfg["backend"]``: ``"gpu"``
    (default) imports the PyTorch backend ``dmf.torch`` (CUDA), ``"cpu"`` imports ``dmf`` (NumPy).

    References:
    [1] S.-i. Koda and S. Saito, JCTC, 20, 2798-2811 (2024). doi: 10.1021/acs.jctc.3c01246
    [2] S.-i. Koda and S. Saito, JCTC, 20, 7176-7187 (2024). doi: 10.1021/acs.jctc.4c00792
    [3] S.-i. Koda and S. Saito, JCTC, 21, 3513-3522 (2025). doi: 10.1021/acs.jctc.4c01549
    """
    dmf_backend = str((dmf_cfg or {}).get("backend", "gpu")).strip().lower()
    try:
        from ase.io import read as ase_read, write as ase_write
        from ase.calculators.mixing import SumCalculator
        if dmf_backend == "cpu":
            from dmf import DirectMaxFlux, interpolate_fbenm
        else:
            from dmf.torch import DirectMaxFlux, interpolate_fbenm
    except Exception as e:
        raise RuntimeError(
            "DMF mode (--mep-mode dmf) requires ase, cyipopt, and pydmf>=1.2 "
            "(`conda install -c conda-forge cyipopt -y`, then `pip install "
            "'pydmf[torch]>=1.2'` for GPU or `pip install 'pydmf>=1.2'` for CPU). "
            f"Import error: {e}"
        ) from e

    from mlmm.workflows.restraints import HarmonicFixAtoms

    def _geom_to_ase(g):
        from io import StringIO
        return ase_read(StringIO(g.as_xyz()), format="xyz")

    fix_atoms = list(sorted(set(map(int, fix_atoms))))

    ref_images = [_geom_to_ase(g) for g in geoms]
    fix_ref_positions = _shared_frozen_reference(ref_images, fix_atoms)
    charge = int(shared_calc.core.model_charge)
    spin = int(shared_calc.core.model_mult)
    for img in ref_images:
        img.info["charge"] = charge
        img.info["spin"] = spin

    # Reuse the already-created heavy core for both DMF and final evaluation.
    ase_calc = MLMMASECalculator(core=shared_calc.core)

    dmf_cfg = fresh_dmf_config(dmf_cfg)
    fbenm_opts: Dict[str, Any] = dict(dmf_cfg.get("fbenm_options", {}))
    cfbenm_opts: Dict[str, Any] = dict(dmf_cfg.get("cfbenm_options", {}))
    dmf_opts: Dict[str, Any] = dict(dmf_cfg.get("dmf_options", {}))
    update_teval = bool(dmf_opts.pop("update_teval", False))
    supports_keep_history = (
        dmf_backend == "gpu"
        and "keep_history" in inspect.signature(DirectMaxFlux.__init__).parameters
    )
    torch_dmf_kwargs = _torch_dmf_runtime_kwargs(
        dmf_backend, dmf_opts, fbenm_opts, cfbenm_opts,
        supports_keep_history=supports_keep_history,
    )
    if supports_keep_history:
        dmf_opts.setdefault("keep_history", False)
    k_fix = float(dmf_cfg.get("k_fix", DMF_KW["k_fix"]))

    # Default-mode IPOPT options: print_level=0 silences the per-iteration
    # IPOPT banner + MUMPS license header that would otherwise print 16+
    # times per pipeline (one per align+scan refinement). Under `-v` we
    # let dmf use its own default (print_level=5 = full table). User-
    # supplied `dmf_cfg["ipopt_options"]` always wins.
    from mlmm.core.utils import is_verbose
    ipopt_opts: Dict[str, Any] = dict(dmf_cfg.get("ipopt_options", {}))
    if "print_level" not in ipopt_opts and not is_verbose():
        ipopt_opts["print_level"] = 0

    # Run FB-ENM interpolation
    emit("\n====== DMF: FB-ENM interpolation ======\n", narrative=True)
    mxflx_fbenm = interpolate_fbenm(
        ref_images,
        nmove=max(1, int(max_nodes)),
        fbenm_only_endpoints=bool(dmf_cfg.get("fbenm_only_endpoints", False)),
        correlated=bool(dmf_cfg.get("correlated", False)),
        sequential=bool(dmf_cfg.get("sequential", False)),
        output_file=str(out_dir_path / "dmf_fbenm_ipopt.out"),
        fbenm_options=fbenm_opts,
        cfbenm_options=cfbenm_opts,
        dmf_options=dmf_opts,
        ipopt_options=ipopt_opts,
    )

    initial_trj = out_dir_path / "dmf_initial_trj.xyz"
    ase_write(initial_trj, mxflx_fbenm.images, format="xyz")
    click.echo(f"[write] Wrote '{initial_trj}' ({len(mxflx_fbenm.images)} images).")

    # Convert initial trajectory to PDB if possible
    if input_paths[0].suffix.lower() == ".pdb" and is_convert_file_enabled():
        try:
            initial_pdb = initial_trj.with_suffix(".pdb")
            convert_xyz_to_pdb(initial_trj, input_paths[0].resolve(), initial_pdb)
            click.echo(f"[convert] Wrote '{initial_pdb}'.")
        except Exception as e:
            click.echo(f"[convert] WARNING: {e}", err=True)

    coefs = mxflx_fbenm.coefs.copy()

    # FB-ENM interpolation and the accurate ML/MM solve are separate GPU
    # phases.  Drop the interpolation cache here, after its final use, so one
    # constant upload serves all interpolation callbacks but no stale cache is
    # carried into the accurate stage.
    _release_dmf_interpolation_cache(mxflx_fbenm)
    del mxflx_fbenm
    gc.collect()
    if dmf_backend == "gpu" and torch.cuda.is_available():
        torch.cuda.empty_cache()

    # Create DirectMaxFlux object
    emit("\n====== DMF: Direct Max Flux optimization ======\n", narrative=True)
    mxflx = DirectMaxFlux(
        ref_images,
        coefs=coefs,
        nmove=max(1, int(max_nodes)),
        update_teval=update_teval,
        remove_rotation_and_translation=bool(
            dmf_opts.get("remove_rotation_and_translation", False)
        ),
        mass_weighted=bool(dmf_opts.get("mass_weighted", False)),
        parallel=bool(dmf_opts.get("parallel", False)),
        eps_vel=float(dmf_opts.get("eps_vel", DMF_KW["dmf_options"]["eps_vel"])),
        eps_rot=float(dmf_opts.get("eps_rot", DMF_KW["dmf_options"]["eps_rot"])),
        beta=float(dmf_opts.get("beta", DMF_KW["dmf_options"]["beta"])),
        **torch_dmf_kwargs,
    )

    # Assign calculators to images
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

    accurate_ipopt_opts = dict(ipopt_opts)
    accurate_ipopt_opts["output_file"] = str(out_dir_path / "dmf_ipopt.out")
    max_cycles = dmf_cfg.get("max_cycles") if isinstance(dmf_cfg, dict) else None
    if max_cycles is not None:
        try:
            max_iter = int(max_cycles)
            if max_iter > 0:
                accurate_ipopt_opts["max_iter"] = max_iter
        except Exception:
            logger.debug("Failed to set ipopt max_iter option", exc_info=True)
    mxflx.add_ipopt_options(accurate_ipopt_opts)
    solve_result = mxflx.solve(tol=resolve_dmf_solve_tol(dmf_cfg))
    converged, ipopt_status, reason = _dmf_solver_outcome(solve_result)
    emit("\n====== DMF: optimization finished ======\n", narrative=True)

    # Evaluate final energies using the PySisyphus calculator for consistency
    from pysisyphus.constants import ANG2BOHR
    energies = []
    for image in mxflx.images:
        elems = image.get_chemical_symbols()
        coords_bohr = np.asarray(image.get_positions(), dtype=float).reshape(-1, 3) * ANG2BOHR
        energies.append(float(shared_calc.get_energy(elems, coords_bohr)["energy"]))
    hei_idx = _select_hei_index(energies)

    # Write final trajectory
    final_trj = out_dir_path / "final_geometries_trj.xyz"
    blocks = []
    for _idx, (image, E) in enumerate(zip(mxflx.images, energies)):
        from io import StringIO
        buf = StringIO()
        ase_write(buf, image, format="xyz")
        s = buf.getvalue()
        lines = s.splitlines()
        if len(lines) >= 2 and lines[0].strip().isdigit():
            lines[1] = f"{E:.12f}"
        blocks.append("\n".join(lines) + "\n")
    with open(final_trj, "w") as f:
        f.write("".join(blocks))
    click.echo(f"[write] Wrote '{final_trj}' with energy.")

    # Convert to PDB
    if input_paths[0].suffix.lower() == ".pdb" and is_convert_file_enabled():
        ref_pdb = input_paths[0].resolve()
        try:
            final_pdb = out_dir_path / "final_geometries.pdb"
            convert_xyz_to_pdb(final_trj, ref_pdb, final_pdb)
            click.echo(f"[convert] Wrote '{final_pdb}'.")
            _apply_bfactor_annotations_inplace(
                final_pdb,
                ml_indices=ml_indices_set or set(),
                freeze_indices=freeze_atoms_final or [],
            )
        except Exception as e:
            click.echo(f"[convert] WARNING: {e}", err=True)

    # Write HEI
    hei_geom = mxflx.images[hei_idx]
    hei_E = energies[hei_idx]
    hei_xyz = out_dir_path / "hei.xyz"
    from io import StringIO
    buf = StringIO()
    ase_write(buf, hei_geom, format="xyz")
    s = buf.getvalue()
    lines = s.splitlines()
    if len(lines) >= 2 and lines[0].strip().isdigit():
        lines[1] = f"{hei_E:.12f}"
        s = "\n".join(lines) + "\n"
    with open(hei_xyz, "w") as f:
        f.write(s)
    click.echo(f"[write] Wrote '{hei_xyz}' (HEI index={hei_idx}).")

    if input_paths[0].suffix.lower() == ".pdb" and is_convert_file_enabled():
        ref_pdb = input_paths[0].resolve()
        hei_pdb = out_dir_path / "hei.pdb"
        try:
            convert_xyz_to_pdb(hei_xyz, ref_pdb, hei_pdb)
            click.echo(f"[convert] Wrote '{hei_pdb}'.")
            _apply_bfactor_annotations_inplace(
                hei_pdb,
                ml_indices=ml_indices_set or set(),
                freeze_indices=freeze_atoms_final or [],
            )
        except Exception as e:
            click.echo(f"[convert] WARNING: {e}", err=True)

    images = tuple(mxflx.images)
    result = DMFMepResult(
        images=images,
        energies=tuple(float(energy) for energy in energies),
        hei_idx=int(hei_idx),
        converged=bool(converged),
        ipopt_status=ipopt_status,
        reason=reason,
    )
    # Image calculators own only light wrappers, but clearing them bounds those
    # references while the caller retains the single heavy ``shared_calc`` core.
    for image in images:
        image.calc = None
    del ase_calc, mxflx
    gc.collect()
    if dmf_backend == "gpu" and torch.cuda.is_available():
        torch.cuda.empty_cache()
    return result



@click.command(
    help="MEP optimization via the Growing String method or Direct Max Flux.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    nargs=2,
    required=True,
    help=("Two endpoint structures in PDB/mmCIF, or XYZ with a corresponding "
          "--ref-pdb for each endpoint."),
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
    show_default="1",
    help="Spin multiplicity (2S+1). Defaults to 1 when omitted.",
)
@click.option(
    "--mep-mode",
    type=click.Choice(["gsm", "dmf"], case_sensitive=False),
    default="gsm",
    show_default=True,
    help="MEP optimizer: Growing String Method (gsm) or Direct Max Flux (dmf).",
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
    "--max-nodes",
    type=int,
    default=GS_KW["max_nodes"],
    show_default=True,
    help=(
        "Number of movable internal images for GSM or DMF "
        "(total images = max_nodes + 2 endpoints)."
    ),
)
@click.option("--max-cycles-gsm", type=int, default=None, show_default="300",
              help="Maximum GSM string-optimizer cycles for the MEP stage.")
@click.option("--max-cycles-dmf", type=int, default=None, show_default="300",
              help=("Maximum IPOPT iterations for the DMF MEP stage. This is a solver "
                    "iteration count, not a string-optimizer cycle count."))
@click.option(
    "--climb/--no-climb",
    default=True,
    show_default=True,
    help="Search for a transition state (climbing image) after path growth.",
)
@click.option(
    "--preopt/--no-preopt",
    # Default True matches GS_KW.fix_first/fix_last semantics and the
    # mlmm-all.py forwarding: endpoints are typically pre-relaxed before
    # string growth to avoid GSM step inflation.
    default=True,
    show_default=True,
    help="Pre-optimize the two endpoint structures with L-BFGS before string growth.",
)
@click.option("--preopt-max-cycles", "preopt_max_cycles", type=int, default=10000, show_default=True,
              help="Maximum L-BFGS cycles for endpoint pre-optimization when --preopt is enabled.")
@click.option(
    "--fix-ends/--no-fix-ends",
    default=True,
    show_default=True,
    help="Fix endpoint structures during path growth.",
)
@click.option(
    "--dump/--no-dump",
    default=False,
    show_default=True,
    help="Dump optimizer trajectory/restarts during the run.",
)
@click.option("-o", "--out-dir", "out_dir", type=str, default=OUT_DIR_PATH_OPT, show_default=True,
              help="Output directory.")
@click.option(
    "--thresh",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default=None,
    show_default="gau",
    help=(
        "Convergence preset for endpoint preoptimization only. "
        "The MEP itself keeps --thresh-gsm / --thresh-dmf."
    ),
)
@click.option(
    "--thresh-gsm",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default=None,
    show_default="gau_loose",
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
    show_default="tight",
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
    help="Validate options and print the execution plan without running path optimization.",
)
@click.option(
    "--parm",
    "real_parm7",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Amber parm7 topology for the enzyme complex (MM layers).",
)
@click.option(
    "--model-pdb",
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
    "--freeze-atoms",
    "freeze_atoms_cli",
    type=str,
    default=None,
    help="Comma-separated 1-based indices to freeze (applied to every image).",
)
@click.option(
    "--movable-cutoff",
    "movable_cutoff",
    type=float,
    default=None,
    show_default="use freeze_atoms",
    help="Distance cutoff (Å) from ML region for movable MM atoms. MM atoms beyond this are frozen. "
         "Providing --movable-cutoff disables --detect-layer.",
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
    show_default="uma",
    help="ML backend for the ONIOM high-level region (default: uma).",
)
@click.option(
    "--embedcharge/--no-embedcharge",
    "embedcharge",
    default=False,
    show_default=True,
    help="Enable the experimental, computationally expensive xTB point-charge delta correction for MLIP/MM.",
)
@click.option(
    "--embedcharge-cutoff",
    "embedcharge_cutoff",
    type=float,
    default=None,
    show_default="12.0",
    help="Distance cutoff (Å) from the ML region for MM point charges used by the xTB delta correction.",
)
@click.option(
    "--link-atom-method",
    "link_atom_method",
    type=click.Choice(["scaled", "fixed"], case_sensitive=False),
    default=None,
    show_default="scaled",
    help="Link-atom position mode: scaled (g-factor, default) or fixed (legacy 1.09/1.01 Å).",
)
@click.option(
    "--mm-backend",
    "mm_backend",
    type=click.Choice(["hessian_ff", "openmm"], case_sensitive=False),
    default=None,
    show_default="hessian_ff",
    help="MM backend (default: hessian_ff). MM Hessians use finite differences by default; set calc.mm_fd: false for the hessian_ff analytical path.",
)
@click.option(
    "--cmap/--no-cmap",
    "use_cmap",
    default=None,
    show_default="cmap",
    help="Preserve CMAP terms in both real and model MM layers. Default: enabled when present in parm7.",
)
@click.option(
    "--out-json/--no-out-json",
    "out_json",
    default=False,
    show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
# Full template PDBs for XYZ→PDB conversion and topology reference
@click.option(
    "--ref-pdb",
    "ref_pdb_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True,
    default=None,
    help=("Full-size template PDBs in the same order as --input. "
          "Required when using XYZ inputs to provide topology and B-factor information."),
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
    ref_pdb_paths: Sequence[Path],
    charge: Optional[int],
    ligand_charge: Optional[str],
    spin: Optional[int],
    mep_mode: str,
    dmf_backend: str,
    max_nodes: int,
    max_cycles_gsm: Optional[int],
    max_cycles_dmf: Optional[int],
    climb: bool,
    preopt: bool,
    preopt_max_cycles: int,
    fix_ends: bool,
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    thresh_gsm: Optional[str],
    thresh_dmf: Optional[str],
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    real_parm7: Path,
    model_pdb: Optional[Path],
    model_indices_str: Optional[str],
    model_indices_one_based: bool,
    detect_layer: bool,
    freeze_atoms_cli: Optional[str],
    movable_cutoff: Optional[float],
    convert_files: bool,
    backend: Optional[str],
    embedcharge: bool,
    embedcharge_cutoff: Optional[float],
    link_atom_method: Optional[str],
    mm_backend: Optional[str],
    use_cmap: Optional[bool],
    out_json: bool,
    precision: Optional[str],
    workers: Optional[int],
    workers_per_node: Optional[int],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: Optional[str],
) -> None:
    set_convert_file_enabled(convert_files)
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

    input_paths = tuple(Path(p) for p in input_paths)
    requested_input_paths = input_paths
    prepared_inputs: List[PreparedInputStructure] = []
    time_start = time.perf_counter()
    out_dir_path = Path(out_dir).resolve()
    try:
        if len(input_paths) != 2:
            click.echo("ERROR: Provide exactly two endpoint structures (-i reactant product).", err=True)
            sys.exit(1)

        # Accept XYZ inputs when a matching --ref-pdb is supplied: overlay the XYZ
        # coordinates onto the reference PDB topology so path-opt runs on full
        # ML/MM PDBs (mirrors the --ref-pdb support already in path-search/scan).
        ref_list = list(ref_pdb_paths) if ref_pdb_paths else []
        for i, src in enumerate(input_paths):
            suffix = src.suffix.lower()
            prepared = prepare_input_structure(src)
            if suffix in {".pdb", ".cif", ".mmcif"}:
                pass
            elif suffix == ".xyz":
                if i >= len(ref_list):
                    raise click.UsageError(
                        f"XYZ input '{src.name}' requires a corresponding --ref-pdb "
                        "for topology/B-factor info."
                    )
                apply_ref_pdb_override(prepared, Path(ref_list[i]))
            else:
                click.echo(
                    f"ERROR: '{src.name}': unsupported format. Use .pdb/.cif/.mmcif or .xyz (with --ref-pdb).",
                    err=True,
                )
                sys.exit(1)
            prepared_inputs.append(prepared)
        input_paths = tuple(prep.source_path for prep in prepared_inputs)

        config_layer_cfg = load_yaml_dict(config_yaml)
        override_layer_cfg = load_yaml_dict(override_yaml)

        mep_mode_kind = mep_mode.strip().lower()

        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        gs_cfg = dict(GS_KW)
        stopt_cfg = dict(STOPT_KW)
        lbfgs_cfg = dict(LBFGS_KW)
        dmf_cfg = fresh_dmf_config()

        apply_yaml_overrides(
            config_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",), ("mlmm",))),
                (gs_cfg, (("gs",),)),
                (stopt_cfg, (("stopt",), ("opt",))),
                (lbfgs_cfg, (("opt", "lbfgs"), ("lbfgs",), ("stopt", "lbfgs"))),
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
        if _is_param_explicit("embedcharge"):
            calc_cfg["embedcharge"] = bool(embedcharge)
        if _is_param_explicit("embedcharge_cutoff"):
            calc_cfg["embedcharge_cutoff"] = embedcharge_cutoff
        geom_cfg["coord_type"] = "cart"  # no microiteration: DLC over MM atoms is meaningless; fixed to cart
        if link_atom_method is not None:
            calc_cfg["link_atom_method"] = str(link_atom_method).lower()
        if mm_backend is not None:
            calc_cfg["mm_backend"] = str(mm_backend).lower()
        if use_cmap is not None:
            calc_cfg["use_cmap"] = use_cmap

        if _is_param_explicit("max_nodes"):
            gs_cfg["max_nodes"] = int(max_nodes)
        # The GSM cycle budget also bounds the fully-grown string; DMF's budget
        # is a separate IPOPT iteration count.
        if _is_param_explicit("max_cycles_gsm") and max_cycles_gsm is not None:
            stopt_cfg["max_cycles"] = int(max_cycles_gsm)
            stopt_cfg["stop_in_when_full"] = int(max_cycles_gsm)
        if _is_param_explicit("max_cycles_dmf") and max_cycles_dmf is not None:
            dmf_cfg["max_cycles"] = int(max_cycles_dmf)
        if _is_param_explicit("dmf_backend"):
            dmf_cfg["backend"] = str(dmf_backend).lower()
        if _is_param_explicit("climb"):
            gs_cfg["climb"] = bool(climb)
            gs_cfg["climb_lanczos"] = bool(climb)
        if _is_param_explicit("fix_ends"):
            gs_cfg["fix_first"] = bool(fix_ends)
            gs_cfg["fix_last"] = bool(fix_ends)
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
        if _is_param_explicit("detect_layer"):
            calc_cfg["use_bfactor_layers"] = bool(detect_layer)
        if _is_param_explicit("movable_cutoff") and movable_cutoff is not None:
            calc_cfg["movable_cutoff"] = float(movable_cutoff)
            calc_cfg["use_bfactor_layers"] = False
        if _is_param_explicit("preopt_max_cycles"):
            lbfgs_cfg["max_cycles"] = int(preopt_max_cycles)

        resolved_charge = charge
        resolved_spin = spin
        for prepared in prepared_inputs:
            resolved_charge, resolved_spin = resolve_charge_spin_or_raise(
                prepared,
                resolved_charge,
                resolved_spin,
                ligand_charge=ligand_charge,
                prefix="[path-opt]",
                model_pdb=model_pdb,
                model_indices_spec=model_indices_str,
                detect_layer=detect_layer,
                yaml_cfg=merged_yaml_cfg,
            )
        # CLI-resolved charge/spin (from -q / -l derivation, or -m / spin_default)
        # always wins over the CALC_KW default carried in calc_cfg.
        # Same pattern as opt.py fix (commit f9a2799).
        calc_cfg["model_charge"] = int(resolved_charge)
        calc_cfg["model_mult"] = int(resolved_spin)

        if model_pdb is not None:
            calc_cfg["model_pdb"] = str(model_pdb)
        calc_cfg["input_pdb"] = str(input_paths[0])
        calc_cfg["real_parm7"] = str(real_parm7)

        apply_yaml_overrides(
            override_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",), ("mlmm",))),
                (gs_cfg, (("gs",),)),
                (stopt_cfg, (("stopt",), ("opt",))),
                (lbfgs_cfg, (("opt", "lbfgs"), ("lbfgs",), ("stopt", "lbfgs"))),
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
            resolve_dmf_solve_tol(dmf_cfg)

        try:
            geom_freeze = _normalize_geom_freeze(geom_cfg.get("freeze_atoms"))
        except click.BadParameter as e:
            click.echo(f"ERROR: {e}", err=True)
            sys.exit(1)
        geom_cfg["freeze_atoms"] = geom_freeze
        _convert_yaml_layer_atoms_1to0(calc_cfg)

        try:
            cli_freeze = _parse_freeze_atoms(freeze_atoms_cli)
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

        # Keep optimizer alignment policy deterministic.
        stopt_cfg["align"] = False
        stopt_cfg["stop_in_when_full"] = int(stopt_cfg.get("max_cycles", STOPT_KW["max_cycles"]))

        out_dir_path = Path(stopt_cfg["out_dir"]).resolve()
        preopt_max_cycles_effective = int(lbfgs_cfg.get("max_cycles", preopt_max_cycles))

        # movable_cutoff implies full distance-based layer assignment.
        detect_layer_enabled = bool(calc_cfg.get("use_bfactor_layers", True))
        model_pdb_cfg = calc_cfg.get("model_pdb")
        if calc_cfg.get("movable_cutoff") is not None:
            if detect_layer_enabled:
                click.echo("[layer] movable_cutoff is set; disabling --detect-layer.", err=True)
            detect_layer_enabled = False
            calc_cfg["use_bfactor_layers"] = False

        layer_source_pdb = input_paths[0]
        path_protected_inputs = (
            *requested_input_paths,
            *(prep.source_path for prep in prepared_inputs),
            *ref_list,
            real_parm7,
            (
                Path(calc_cfg["input_pdb"])
                if calc_cfg.get("input_pdb")
                else None
            ),
            (
                Path(calc_cfg["real_parm7"])
                if calc_cfg.get("real_parm7")
                else None
            ),
            model_pdb,
            (
                Path(calc_cfg["model_pdb"])
                if calc_cfg.get("model_pdb")
                else None
            ),
            config_yaml,
            (
                Path(calc_cfg["calc_file"])
                if calc_cfg.get("calc_file")
                else None
            ),
        )
        _reject_path_output_collisions(
            out_dir_path,
            path_protected_inputs,
        )
        if detect_layer_enabled and layer_source_pdb.suffix.lower() != ".pdb":
            click.echo("ERROR: --detect-layer requires a PDB input.", err=True)
            sys.exit(1)

        if show_config:
            click.echo(
                pretty_block(
                    "yaml_layers",
                    {
                        "config": None if config_yaml is None else str(config_yaml),
                        "override_yaml": None if override_yaml is None else str(override_yaml),
                        "merged_keys": sorted(merged_yaml_cfg.keys()),
                    },
                force=True)
            )

        effective_max_cycles = (
            dmf_cfg.get("max_cycles", 0)
            if mep_mode_kind == "dmf"
            else stopt_cfg.get("max_cycles", 0)
        )
        cycles_hint = (
            "--max-cycles-dmf" if mep_mode_kind == "dmf" else "--max-cycles-gsm"
        )
        if int(effective_max_cycles) <= 0:
            raise click.BadParameter(
                f"{cycles_hint} must be at least 1.",
                param_hint=cycles_hint,
            )

        validate_endpoint_atom_identities(prepared_inputs)

        if dry_run:
            if model_pdb_cfg is not None:
                model_region_source = "model_pdb"
            elif model_indices:
                model_region_source = "model_indices"
            elif detect_layer_enabled:
                model_region_source = "bfactor"
            else:
                click.echo("ERROR: Provide --model-pdb or --model-indices when B-factor layer detection is disabled in the configuration.", err=True)
                sys.exit(1)
            if (
                not detect_layer_enabled
                and model_pdb_cfg is None
                and model_indices
                and layer_source_pdb.suffix.lower() != ".pdb"
            ):
                click.echo("ERROR: --model-indices requires a PDB input.", err=True)
                sys.exit(1)
            click.echo(
                pretty_block(
                    "dry_run_plan",
                    {
                        "input_endpoints": [str(p) for p in input_paths],
                        "output_dir": str(out_dir_path),
                        "mep_mode": mep_mode_kind,
                        "fix_ends": bool(gs_cfg.get("fix_first", False) and gs_cfg.get("fix_last", False)),
                        "detect_layer": bool(detect_layer_enabled),
                        "model_region_source": model_region_source,
                        "model_indices_count": 0 if not model_indices else len(model_indices),
                        "preopt": bool(preopt),
                        "preopt_max_cycles": int(preopt_max_cycles_effective),
                        "will_run_path_opt": True,
                        "will_write_summary": True,
                        "backend": calc_cfg.get("backend", "uma"),
                        "embedcharge": bool(calc_cfg.get("embedcharge", False)),
                    },
                )
            )
            click.echo("[dry-run] Validation complete. Path optimization execution was skipped.")
            emit(
                format_elapsed("[time] Elapsed Time for Path Opt", time_start),
                narrative=True,
            )
            return

        try:
            model_pdb_path, layer_info = resolve_ml_layer_assignment(
                source_path=layer_source_pdb,
                out_dir_path=out_dir_path,
                model_pdb=model_pdb_cfg,
                model_indices=model_indices,
                detect_layer=detect_layer_enabled,
                hess_cutoff=calc_cfg.get("hess_cutoff"),
                movable_cutoff=calc_cfg.get("movable_cutoff"),
                calc_cfg=calc_cfg,
                protected_inputs=path_protected_inputs,
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

        for key in ("input_pdb", "real_parm7", "model_pdb", "mm_fd_dir"):
            val = calc_cfg.get(key)
            if val:
                calc_cfg[key] = str(Path(val).expanduser().resolve())

        # For display: resolved configuration (show only non-default values)
        echo_geom = format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")
        echo_calc = format_freeze_atoms_for_echo(filter_calc_for_echo(calc_cfg), key="freeze_atoms")
        echo_gs = strip_inherited_keys(gs_cfg, GS_KW, mode="same")
        echo_stopt = strip_inherited_keys({**stopt_cfg, "out_dir": str(out_dir_path)}, STOPT_KW, mode="same")
        echo_lbfgs = strip_inherited_keys({**lbfgs_cfg, "out_dir": stopt_cfg.get("out_dir")}, LBFGS_KW, mode="same")

        click.echo(pretty_block("geom", echo_geom))
        click.echo(pretty_block("calc", echo_calc))
        if mep_mode_kind == "gsm":
            click.echo(pretty_block("gs", echo_gs))
            click.echo(pretty_block("stopt", echo_stopt))
            click.echo(pretty_block("lbfgs", echo_lbfgs))
        elif mep_mode_kind == "dmf":
            click.echo(pretty_block("dmf", dmf_cfg))
        click.echo(
            pretty_block(
                "run_flags",
                {
                    "mep_mode": mep_mode_kind,
                    "preopt": bool(preopt),
                    "preopt_max_cycles": int(preopt_max_cycles_effective),
                    "fix_ends": bool(gs_cfg.get("fix_first", False) and gs_cfg.get("fix_last", False)),
                },
            )
        )

        out_dir_path = _prepare_path_output_dir(out_dir_path)

        source_paths = [prep.source_path for prep in prepared_inputs]

        # Pre-compute ML-region indices (0-based in ref PDB atom order) for later PDB annotation
        ml_indices_set: Set[int] = set()
        try:
            ref_pdb_for_map = source_paths[0]
            if ref_pdb_for_map.suffix.lower() == ".pdb":
                ml_indices_set = _compute_ml_indices_from_model_and_ref(
                    ref_pdb_for_map.resolve(),
                    Path(calc_cfg["model_pdb"]).resolve(),
                )
                click.echo(f"[annotate] ML-region atoms mapped: {len(ml_indices_set)}")
        except Exception as e:
            click.echo(f"[annotate] WARNING: Failed to pre-compute ML-region indices: {e}", err=True)

        # Load endpoints (if PDB, merge in link-parent freezing)
        geoms = _load_two_endpoints(
            inputs=prepared_inputs,
            coord_type=geom_cfg.get("coord_type", "cart"),
            base_freeze=geom_cfg.get("freeze_atoms", []),
        )

        # Shared ML/MM calculator (reuse the same instance for all images)
        shared_calc = mlmm(**calc_cfg)

        echo_resolved_device()

        # optional endpoint pre-optimization
        if preopt:
            preopt_completed = 0
            preopt_errors: List[str] = []
            try:
                emit("\n====== Pre-optimizing endpoints (LBFGS) ======\n", narrative=True)
                pre_dir_base = out_dir_path / "preopt"
                for i, g in enumerate(geoms):
                    try:
                        g.set_calculator(shared_calc)
                    except Exception:
                        logger.debug("Failed to set calculator on geometry", exc_info=True)
                    subdir = pre_dir_base / f"end{i:02d}"
                    subdir.mkdir(parents=True, exist_ok=True)
                    lbfgs_args = dict(lbfgs_cfg)
                    lbfgs_args.update({
                        "out_dir": str(subdir),
                        "max_cycles": int(preopt_max_cycles_effective),
                    })
                    optimizer = LBFGS(g, **lbfgs_args)
                    optimizer.run()
                    try:
                        final_xyz_path = optimizer.final_fn if isinstance(optimizer.final_fn, Path) else Path(optimizer.final_fn)
                        g_new = geom_loader(
                            final_xyz_path,
                            coord_type=geom_cfg.get("coord_type", "cart"),
                            freeze_atoms=getattr(g, "freeze_atoms", []),
                        )
                        try:
                            g_new.freeze_atoms = np.array(getattr(g, "freeze_atoms", []), dtype=int)
                        except Exception:
                            logger.debug("Failed to set freeze_atoms on new geometry", exc_info=True)
                        geoms[i] = g_new
                        preopt_completed += 1
                    except Exception as e:
                        preopt_errors.append(f"endpoint #{i}: {e}")
                        click.echo(f"[preopt] WARNING: Failed to reload optimized endpoint #{i}: {e}", err=True)
            except Exception as e:
                preopt_errors.append(str(e))
                click.echo(f"[preopt] WARNING: Endpoint pre-optimization stopped: {e}", err=True)
            click.echo(
                f"[preopt] Pre-optimized {preopt_completed}/{len(geoms)} endpoints"
                + (f" ({len(preopt_errors)} error(s))." if preopt_errors else ".")
            )

        # By default, apply external Kabsch alignment (if freeze_atoms exist, use only them)
        align_thresh = str(stopt_cfg.get("thresh", "gau"))
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

        # Collect freeze_atoms for DMF
        # No try/except: swallowing here silently hands DMF an EMPTY frozen set, i.e. runs the
        # segment unconstrained. The expression cannot fail for the Geometry objects built just
        # above (freeze_atoms is an int array, and getattr already covers absence).
        fix_atoms: List[int] = sorted(
            {int(i) for g in geoms for i in getattr(g, "freeze_atoms", [])}
        )

        if mep_mode_kind == "dmf":
            try:
                dmf_res = _run_dmf_mep(
                    geoms,
                    shared_calc,
                    out_dir_path,
                    input_paths,
                    int(gs_cfg["max_nodes"]),
                    fix_atoms,
                    dmf_cfg=dmf_cfg,
                    ml_indices_set=ml_indices_set,
                    freeze_atoms_final=freeze_atoms_final,
                )
            except Exception as e:
                if str(dmf_cfg.get("backend", "gpu")).lower() != "cpu" and _is_cuda_oom(e):
                    click.echo(
                        "[dmf] GPU out of memory. Retry with `--dmf-backend cpu` "
                        "(NumPy backend; slower but not limited by GPU memory).",
                        err=True,
                    )
                else:
                    tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
                    click.echo(f"[dmf] ERROR: DMF optimization failed:\n{textwrap.indent(tb, '  ')}", err=True)
                _write_error_json(
                    out_dir_path,
                    "path-opt",
                    e,
                    "DMFError",
                    time_start,
                )
                sys.exit(3)
            if out_json:
                from mlmm.core.utils import write_result_json

                result_data_dmf = _build_dmf_result_data(dmf_res, calc_cfg)
                for ext in (".pdb", ".gjf"):
                    f = out_dir_path / f"hei{ext}"
                    if f.exists():
                        result_data_dmf["files"][f"hei_{ext[1:]}"] = f.name
                # The DMF path is a required
                # leaf usable only when the IPOPT solve explicitly converged. The
                # scientific_status path routes convergence through the
                # canonical criterion (IPOPT status 0 or 1), matching path_search,
                # so the additive axis is consistent across both DMF producers.
                # The legacy convergence-aware ``status``/``converged`` fields
                # (status==0) in result_data_dmf are intentionally left untouched.
                from mlmm.workflows._outcomes import (
                    aggregate_workflow_truth as _agg_truth,
                    attach_outcomes as _attach,
                    ipopt_status_to_converged,
                    make_leaf as _mk_leaf,
                )
                _dmf_leaf_conv, _dmf_leaf_reason = ipopt_status_to_converged(dmf_res.ipopt_status)
                _dmf_leaf = _mk_leaf(
                    "path-opt",
                    "dmf_mep",
                    executed=True,
                    converged=_dmf_leaf_conv,
                    artifacts=["final_geometries_trj.xyz"],
                    reason=_dmf_leaf_reason or dmf_res.reason or "",
                )
                _attach(
                    result_data_dmf,
                    truth=_agg_truth([_dmf_leaf], ["dmf_mep"]),
                    stage_outcomes=[_dmf_leaf],
                )
                write_result_json(
                    out_dir_path, result_data_dmf,
                    command="path-opt",
                    elapsed_seconds=time.perf_counter() - time_start,
                )

            emit(
                format_elapsed("[time] Elapsed Time for Path Opt", time_start),
                narrative=True,
            )
            return

        for g in geoms:
            g.set_calculator(shared_calc)

        def calc_getter():
            # Used when GrowingString generates new nodes
            return shared_calc

        gs = GrowingString(
            images=geoms,
            calc_getter=calc_getter,
            **gs_cfg,
        )

        # StringOptimizer expects 'out_dir' under the key "out_dir"
        opt_args = dict(stopt_cfg)
        opt_args["out_dir"] = str(out_dir_path)

        optimizer = StringOptimizer(
            geometry=gs,
            **{k: v for k, v in opt_args.items() if k != "type"}  # 'type' is just a tag
        )

        optimizer.run()

        final_trj = out_dir_path / "final_geometries_trj.xyz"
        try:
            try:
                energies = np.array(gs.energy, dtype=float)
                blocks = []
                for _idx, (geom, E) in enumerate(zip(gs.images, energies)):
                    s = geom.as_xyz()
                    lines = s.splitlines()
                    if len(lines) >= 2 and lines[0].strip().isdigit():
                        lines[1] = f"{E:.12f}"
                    s_mod = "\n".join(lines)
                    if not s_mod.endswith("\n"):
                        s_mod += "\n"
                    blocks.append(s_mod)
                annotated = "".join(blocks)
                with open(final_trj, "w") as f:
                    f.write(annotated)
                click.echo(f"[write] Wrote '{final_trj}' with energy.")
            except Exception:
                with open(final_trj, "w") as f:
                    f.write(gs.as_xyz())
                click.echo(f"[write] Wrote '{final_trj}'.")

            if input_paths[0].suffix.lower() == ".pdb" and is_convert_file_enabled():
                ref_pdb = input_paths[0].resolve()

                try:
                    out_pdb = out_dir_path / "final_geometries.pdb"
                    convert_xyz_to_pdb(final_trj, ref_pdb, out_pdb)
                    click.echo(f"[convert] Wrote '{out_pdb}'.")
                    # Annotate B-factors for ML & freeze atoms
                    _apply_bfactor_annotations_inplace(
                        out_pdb,
                        ml_indices=ml_indices_set,
                        freeze_indices=freeze_atoms_final,
                    )
                except Exception as e:
                    click.echo(f"[convert] WARNING: Failed to convert MEP path trajectory to PDB: {e}", err=True)

        except Exception as e:
            click.echo(f"[write] ERROR: Failed to write final trajectory: {e}", err=True)
            sys.exit(4)

        try:
            energies = np.array(gs.energy, dtype=float)
            hei_idx = _select_hei_index(energies)

            hei_geom = gs.images[hei_idx]
            hei_E = float(energies[hei_idx])

            hei_xyz = out_dir_path / "hei.xyz"
            s = hei_geom.as_xyz()
            lines = s.splitlines()
            if len(lines) >= 2 and lines[0].strip().isdigit():
                lines[1] = f"{hei_E:.12f}"
                s = "\n".join(lines) + ("\n" if not s.endswith("\n") else "")
            with open(hei_xyz, "w") as f:
                f.write(s)
            click.echo(f"[write] Wrote '{hei_xyz}'.")

            ref_pdb = None
            if source_paths[0].suffix.lower() == ".pdb":
                ref_pdb = source_paths[0].resolve()
            if ref_pdb is not None and is_convert_file_enabled():
                hei_pdb = out_dir_path / "hei.pdb"
                convert_xyz_to_pdb(hei_xyz, ref_pdb, hei_pdb)
                click.echo(f"[convert] Wrote '{hei_pdb}'.")
                # Annotate B-factors for ML & freeze atoms
                _apply_bfactor_annotations_inplace(
                    hei_pdb,
                    ml_indices=ml_indices_set,
                    freeze_indices=freeze_atoms_final,
                )
            else:
                click.echo("[convert] Skipped 'hei.pdb' (no PDB reference among inputs).")

        except Exception as e:
            click.echo(f"[HEI] ERROR: Failed to dump HEI: {e}", err=True)
            sys.exit(5)

        if out_json:
            from mlmm.core.utils import calculator_provenance, write_result_json
            from pysisyphus.constants import AU2KCALPERMOL as _AU2KCAL
            _gsm_energies = list(map(float, energies))
            _gsm_hei = int(hei_idx)
            _gsm_hei_E = float(_gsm_energies[_gsm_hei])
            _gsm_e0 = float(_gsm_energies[0])
            _gsm_eN = float(_gsm_energies[-1])
            _barrier = (_gsm_hei_E - _gsm_e0) * _AU2KCAL
            _delta = (_gsm_eN - _gsm_e0) * _AU2KCAL
            result_data_gsm: Dict[str, Any] = {
                "mep_mode": "gsm",
                **calculator_provenance(calc_cfg),
                "charge": calc_cfg.get("model_charge"),
                "spin": calc_cfg.get("model_mult"),
                "reactant_energy_hartree": float(_gsm_e0),
                "product_energy_hartree": float(_gsm_eN),
                "image_energies_hartree": [float(e) for e in _gsm_energies],
                "n_images": len(_gsm_energies),
                "hei_index": _gsm_hei,
                "hei_energy_hartree": _gsm_hei_E,
                "barrier_kcal": round(_barrier, 6),
                "delta_kcal": round(_delta, 6),
                "files": {
                    "final_geometries_trj_xyz": "final_geometries_trj.xyz",
                    "hei_xyz": "hei.xyz",
                },
            }
            for ext in (".pdb",):
                f = out_dir_path / f"hei{ext}"
                if f.exists():
                    result_data_gsm["files"][f"hei_{ext[1:]}"] = f.name
            # The GSM path is a required leaf
            # usable only when the StringOptimizer explicitly converged. The legacy
            # convergence-aware ``status``/``converged`` fields are left untouched.
            from mlmm.workflows._outcomes import (
                aggregate_workflow_truth as _agg_truth,
                attach_outcomes as _attach,
                make_leaf as _mk_leaf,
                optimizer_converged_bit as _optimizer_converged_bit,
            )
            _converged = _optimizer_converged_bit(optimizer)
            result_data_gsm["status"] = (
                "converged"
                if _converged is True
                else ("not_converged" if _converged is False else "completed")
            )
            result_data_gsm["converged"] = _converged
            _gsm_leaf = _mk_leaf(
                "path-opt",
                "gsm_mep",
                executed=True,
                converged=_converged,
                artifacts=["final_geometries_trj.xyz"],
            )
            _attach(
                result_data_gsm,
                truth=_agg_truth([_gsm_leaf], ["gsm_mep"]),
                stage_outcomes=[_gsm_leaf],
            )
            write_result_json(
                out_dir_path, result_data_gsm,
                command="path-opt",
                elapsed_seconds=time.perf_counter() - time_start,
            )

        # summary.md and key_* outputs are disabled.
        emit(
            format_elapsed("[time] Elapsed Time for Path Opt", time_start),
            narrative=True,
        )

    except OptimizationError as e:
        _write_error_json(
            out_dir_path, "path-opt", e, "OptimizationError", time_start
        )
        click.echo(f"ERROR: Path optimization failed — {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except _PathOutputCollisionError:
        raise
    except Exception as e:
        render_cli_exception(
            e,
            label="path optimization",
            out_dir=out_dir_path,
            command="path-opt",
            time_start=time_start,
        )
    finally:
        for prepared in prepared_inputs:
            prepared.cleanup()
        # Release GPU memory so subsequent pipeline stages don't OOM.
        # `= None` decref's the heavy refs; `del` then removes names from
        # the local frame so torch.nn.Module hooks / closures cannot retain.
        shared_calc = gs = geoms = None
        del shared_calc, gs, geoms
        gc.collect()  # break cyclic refs inside torch.nn.Module
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


if __name__ == "__main__":
    cli()
