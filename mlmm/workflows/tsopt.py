"""Partial-Hessian transition-state optimization for ML/MM models."""
# DOMAIN_PURE

from __future__ import annotations

import contextlib
import gc
import io
from itertools import count
import logging
import sys
import time
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import click
from mlmm.core.output import emit
import numpy as np
import torch
from ase import Atoms
import ase.units as units
from ase.data import atomic_masses
from ase.io import read as ase_read
from ase.io import write

# ---------------- pysisyphus / mlmm imports ----------------
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.intcoords.exceptions import RebuiltInternalsException
from pysisyphus.constants import BOHR2ANG, AMU2AU, AU2EV
from pysisyphus._array import active_square
from pysisyphus.calculators.Dimer import Dimer  # Dimer calculator (orientation-projected forces)
from pysisyphus.tr_projection import (
    active_tr_basis,
    compact_project_hessian,
    full_cartesian_tr_basis,
    normalize_tr_projection_mode,
    project_hessian_inplace,
)

# RS-I-RFO optimizer for heavy mode
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer
from pysisyphus.tsoptimizers.TRIM import TRIM  # Helgaker trust-region image-min TS opt
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer  # Banerjee P-RFO TS opt
from pysisyphus.TablePrinter import TablePrinter

# local helpers from mlmm
from mlmm.backends.mlmm_calc import mlmm, mlmm_mm_only
from mlmm.core.defaults import (
    OUT_DIR_TSOPT,
    GEOM_KW_DEFAULT,
    MLMM_CALC_KW,
    OPT_BASE_KW,
    LBFGS_KW,
    DIMER_KW,
    HESSIAN_DIMER_KW,
    FREQ_KW,
    RSIRFO_KW,
    MICROITER_KW,
    TSOPT_MODE_ALIASES,
    TS_IMAG_SOFT_WARN_CM,
    BFACTOR_ML,
    BFACTOR_MOVABLE_MM,
    BFACTOR_FROZEN,
    THRESH_CHOICES,
)
from mlmm.io.path_mode_cache import read_reference_mode_candidates
from mlmm.workflows.opt import (
    _parse_freeze_atoms as _parse_freeze_atoms_opt,
    _normalize_geom_freeze as _normalize_geom_freeze_opt,
    _convert_yaml_layer_atoms_1to0,
)
from mlmm.workflows.charge_prep import resolve_charge_spin_or_raise
from mlmm.core.utils import (
    append_xyz_trajectory as _append_xyz_trajectory,
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
    prepare_input_structure,
    apply_ref_pdb_override,
    parse_indices_string,
    resolve_ml_layer_assignment,
    update_pdb_bfactors_from_layers,
    normalize_choice,
    yaml_section_has_key,
    echo_resolved_device,
    emit_optimizer_terminal_status,
    finalize_microiter_macro_convergence,
    optimizer_cycle_count,
    optional_positive_int,
)
from mlmm.workflows._microiteration import (
    MicroiterationOutcome,
    OptimizerOutcome,
    PartitionError,
    build_aggregate,
    micro_reached_force_equilibrium,
    resolve_partition_from_core,
    describe_micro_stop,
    macro_progress_due,
)
from mlmm.cli.common_options import (
    add_ml_layer_detection_options,
    add_precision_option,
    add_backend_model_option,
    add_calc_file_option,
    add_workers_options,
    add_deterministic_option,
    add_coord_type_option,
    add_print_every_option,
    add_allow_charge_mult_mismatch_option,
)
from mlmm.cli.decorators import (
    resolve_yaml_sources,
    load_merged_yaml_cfg,
    make_is_param_explicit,
    _write_error_json,
    render_cli_exception,
)
from mlmm.workflows.freq import (
    _calc_full_hessian_torch as _freq_calc_full_hessian_torch,
    _torch_device,
    _mass_weighted_hessian,
    _align_three_layer_hessian_targets,
    _ordered_hessian_coverage_atoms,
    _reconcile_hessian_analysis_basis,
    _resolve_active_atom_indices,
)
from pysisyphus.normal_modes import (
    DEFAULT_FREQUENCY_ZERO_CUTOFF_CM,
    filter_resolved_modes,
    normalize_frequency_zero_cutoff_cm,
    resolved_frequency_mask,
    resolved_imaginary_mask,
)

logger = logging.getLogger(__name__)


@dataclass
class _OptimizationCycleLedger:
    """Command-level optimization-cycle budget shared by every heavy TS trial."""

    limit: Optional[int]
    spent: int = 0

    def __post_init__(self) -> None:
        self.limit = None if self.limit is None else int(self.limit)
        self.spent = int(self.spent)
        if self.limit is not None and self.limit < 1:
            raise ValueError("Optimization cycle limit must be at least 1.")
        if self.spent < 0 or (
            self.limit is not None and self.spent > self.limit
        ):
            raise ValueError("Initial optimization cycle usage exceeds its limit.")

    @property
    def remaining(self) -> Optional[int]:
        return None if self.limit is None else self.limit - self.spent

    def debit(self, cycles: int) -> None:
        cycles = int(cycles)
        if cycles < 0:
            raise ValueError("Optimization cycle usage cannot be negative.")
        if self.remaining is not None and cycles > self.remaining:
            raise RuntimeError(
                "Optimizer exceeded the command-level --max-cycles budget "
                f"({cycles} requested, {self.remaining} remaining)."
            )
        self.spent += cycles


def _set_cartesian_flatten_coords(geom, cart_coords: np.ndarray) -> None:
    """Install a Cartesian TS trial and accept a completed internals rebuild."""

    try:
        geom.cart_coords = np.asarray(cart_coords, dtype=float).reshape(-1)
    except RebuiltInternalsException:
        # The Cartesian setter installs the requested geometry before it uses
        # this exception to report that its primitive internals were rebuilt.
        geom.clear()


# TS optimizer class map. All three classes inherit from TSHessianOptimizer
# and share the kwargs surface (`_build_rsirfo_kwargs`). Microiter macro
# stays RSIRFO-specific (different image-function / line-search math).
TSOPT_CLASS_MAP = {"rsirfo": RSIRFOptimizer, "trim": TRIM, "rsprfo": RSPRFOptimizer}
PATH_MODE_RESTART_AMPLITUDES_ANG = (-0.10, 0.10, -0.20, 0.20)
FLATTEN_RETRY_HIGHER_ORDER_CHECKS = 3


def _restart_microiteration_carry(
    outcome: Dict[str, Any],
) -> Tuple[Optional[Dict[str, Any]], Optional[int]]:
    """Extract the additive microiteration block + executed micro-cycle total from
    a ``_run_microiter_tsopt`` outcome so a restart/multistart/flatten selection
    can re-anchor the serialized ``microiteration`` object on the run that
    produced the FINAL geometry (never a superseded initial run).

    Returns ``(microiteration_result_object, micro_cycles)`` on the
    microiteration path, or ``(None, None)`` for an ordinary (non-microiter)
    restart whose outcome carries no :class:`MicroiterationOutcome`.
    """

    mi_obj = outcome.get("outcome")
    if mi_obj is None:
        return None, None
    return mi_obj.to_result_object(), int(outcome.get("micro_cycles", 0))


def _optimizer_safeguard_payload(optimizer) -> Dict[str, Any]:
    """Return machine-readable TS-mode selection and rollback diagnostics."""
    return {
        "rejected_mode_loss_trials": int(
            getattr(optimizer, "rejected_mode_loss_steps", 0)
        ),
        "exact_saddle_checks": int(getattr(optimizer, "exact_saddle_checks", 0)),
        "saddle_recovery_steps": int(
            getattr(optimizer, "saddle_recovery_steps", 0)
        ),
        "last_exact_n_imaginary": getattr(
            optimizer, "_last_exact_n_imaginary", None
        ),
        "last_exact_validation": getattr(
            optimizer, "_last_exact_validation", "unavailable"
        ),
        "last_exact_failure_reason": getattr(
            optimizer, "_last_exact_failure_reason", None
        ),
        "initial_reference_root_index": getattr(
            optimizer, "_initial_reference_root_index", None
        ),
        "initial_reference_root_overlap": getattr(
            optimizer, "_initial_reference_root_overlap", None
        ),
        "initial_reference_root_eigenvalue": getattr(
            optimizer, "_initial_reference_root_eigenvalue", None
        ),
        "last_exact_target_mode_index": getattr(
            optimizer, "_last_exact_target_mode_index", None
        ),
        "last_exact_target_mode_overlap": getattr(
            optimizer, "_last_exact_target_mode_overlap", None
        ),
        "last_exact_target_mode_is_negative": getattr(
            optimizer, "_last_exact_target_mode_is_negative", None
        ),
        "last_exact_target_mode_reanchored": bool(
            getattr(optimizer, "_last_exact_target_mode_reanchored", False)
        ),
    }


def _hessian_postprocessing_is_ready(optimizer: Any) -> bool:
    """Whether numerical convergence authorizes terminal PHVA."""
    return bool(
        optimizer is not None
        and getattr(optimizer, "is_converged", False)
        and not getattr(optimizer, "_last_exact_failure_reason", None)
    )


def _saddle_validation_from_count(n_imaginary: Optional[int]) -> str:
    if n_imaginary is None:
        return "unavailable"
    if int(n_imaginary) == 1:
        return "first_order"
    if int(n_imaginary) > 1:
        return "higher_order"
    return "no_imaginary"


def _optimizer_exact_frequency_data(
    optimizer: Any, geometry: Any
) -> Optional[Tuple[np.ndarray, torch.Tensor, Dict[str, Any], Any]]:
    """Reuse terminal exact PHVA owned by the optimizer at this geometry."""
    if optimizer is None:
        return None
    coords = getattr(optimizer, "_last_exact_cart_coords", None)
    freqs = getattr(optimizer, "_last_exact_frequencies_cm", None)
    modes = getattr(optimizer, "_last_exact_modes", None)
    if coords is None or freqs is None or modes is None:
        return None
    current = np.asarray(geometry.cart_coords, dtype=float).reshape(-1)
    checked = np.asarray(coords, dtype=float).reshape(-1)
    if current.shape != checked.shape or not np.allclose(
        current, checked, rtol=0.0, atol=1.0e-12
    ):
        return None
    modes_t = (
        modes.detach().cpu().clone()
        if isinstance(modes, torch.Tensor)
        else torch.as_tensor(np.asarray(modes), dtype=torch.float64).clone()
    )
    projection = dict(getattr(optimizer, "_last_rigid_projection_info", {}) or {})
    projection.update({
        "source": "optimizer_terminal_exact_phva",
        "reused_without_hessian_recalculation": True,
    })
    exact_hessian = getattr(optimizer, "cur_H", None)
    if isinstance(exact_hessian, torch.Tensor):
        exact_hessian = exact_hessian.detach().cpu().clone()
    elif exact_hessian is not None:
        exact_hessian = np.array(exact_hessian, dtype=float, copy=True)
    return (
        np.asarray(freqs, dtype=float).copy(),
        modes_t,
        projection,
        exact_hessian,
    )

def _mirrored_flatten_start(
    saddle_coords: np.ndarray,
    primary_start: np.ndarray,
) -> np.ndarray:
    """Return the opposite signed displacement about a saddle geometry."""
    saddle = np.asarray(saddle_coords, dtype=float)
    primary = np.asarray(primary_start, dtype=float)
    if saddle.shape != primary.shape:
        raise ValueError("saddle and flatten-start coordinates must have equal shapes")
    return 2.0 * saddle - primary


def _flatten_branch_needs_alternate(result: Dict[str, Any]) -> bool:
    optimizer = result["optimizer"]
    target_negative = (
        getattr(optimizer, "_last_exact_target_mode_is_negative", None) is True
    )
    return int(result["n_imag"]) != 1 or not target_negative


def _effective_flatten_iterations(
    configured: int,
    *,
    has_reference_mode: bool,
    n_imag: int,
    target_mode_is_negative: Optional[bool],
) -> Tuple[int, bool]:
    """Return the explicit flatten budget and whether path safety vetoed it."""
    iterations = max(int(configured), 0)
    vetoed = (
        iterations > 0
        and has_reference_mode
        and int(n_imag) > 1
        and target_mode_is_negative is not True
    )
    return (0 if vetoed else iterations), vetoed


def _transported_path_mode_full(
    optimizer,
    geometry,
    fallback: Optional[np.ndarray],
) -> Optional[np.ndarray]:
    """Return the latest path-correlated mode as a full Cartesian unit vector."""
    mode = None
    if geometry.coord_type in ("cart", "cartesian"):
        for attr in (
            "_last_exact_physical_mode",
            "_physical_ts_mode",
            "_saddle_recovery_mode",
        ):
            mode = getattr(optimizer, attr, None)
            if mode is not None:
                break
        if mode is None and len(getattr(optimizer, "ts_modes", ())):
            mode = optimizer.ts_modes[0]
        if isinstance(mode, torch.Tensor):
            mode = mode.detach().cpu().numpy()
        if mode is not None:
            mode = np.asarray(mode, dtype=float).reshape(-1)
            if mode.size != geometry.cart_coords.size:
                try:
                    mode = optimizer.full_from_active(mode)
                except (AttributeError, IndexError, ValueError):
                    mode = None
    if mode is None and fallback is not None:
        mode = np.asarray(fallback, dtype=float).reshape(-1)
    if mode is None or mode.size != geometry.cart_coords.size:
        return None
    norm = float(np.linalg.norm(mode))
    if not np.all(np.isfinite(mode)) or norm <= 0.0:
        return None
    return mode / norm


def _initial_path_root_mode_full(optimizer, geometry) -> Optional[np.ndarray]:
    if geometry.coord_type not in ("cart", "cartesian"):
        return None
    mode = getattr(optimizer, "_initial_reference_root_mode", None)
    if isinstance(mode, torch.Tensor):
        mode = mode.detach().cpu().numpy()
    if mode is None:
        return None
    mode = np.asarray(mode, dtype=float).reshape(-1)
    if mode.size != geometry.cart_coords.size:
        try:
            mode = optimizer.full_from_active(mode)
        except (AttributeError, IndexError, ValueError):
            return None
    norm = float(np.linalg.norm(mode))
    if mode.size != geometry.cart_coords.size or not np.all(np.isfinite(mode)) or norm <= 0.0:
        return None
    return mode / norm


def _path_restart_mode_candidates(
    optimizer,
    geometry,
    reference_modes: Sequence[np.ndarray],
    reference_labels: Optional[Sequence[str]] = None,
) -> List[Tuple[str, np.ndarray]]:
    """Return cached MEP candidates plus a distinct initial soft root."""
    labels = list(reference_labels or ())
    candidates: List[Tuple[str, np.ndarray]] = []
    for index, raw in enumerate(reference_modes):
        mode = np.asarray(raw, dtype=float).reshape(-1)
        norm = float(np.linalg.norm(mode))
        if (
            mode.size != geometry.cart_coords.size
            or not np.all(np.isfinite(mode))
            or not np.isfinite(norm)
            or norm <= 0.0
        ):
            continue
        unit = mode / norm
        if any(abs(float(np.dot(unit, prior))) >= 1.0 - 1.0e-8 for _, prior in candidates):
            continue
        label = labels[index] if index < len(labels) else f"mep-candidate-{index + 1}"
        candidates.append((f"mep-{label}", unit))
    primary = candidates[0][1] if candidates else None
    soft_root = _initial_path_root_mode_full(optimizer, geometry)
    if soft_root is not None and (
        primary is None or abs(float(np.dot(primary, soft_root))) < 0.95
    ):
        if not any(abs(float(np.dot(soft_root, prior))) >= 1.0 - 1.0e-8 for _, prior in candidates):
            candidates.append(("initial-soft-root", soft_root))
    return candidates

def _force_ts_reject_uphill_off(kwargs: Dict[str, Any]) -> Dict[str, Any]:
    """Return TS optimizer kwargs with physical-energy rejection disabled."""
    effective = dict(kwargs)
    effective["reject_uphill"] = False
    return effective


def _resolve_shared_optimizer_value(
    opt_cfg: Dict[str, Any],
    downstream_cfg: Dict[str, Any],
    key: str,
    *,
    opt_explicit: bool,
    downstream_explicit: bool,
    downstream_default: Any,
    downstream_section: str,
) -> None:
    """Resolve one duplicated optimizer setting without silent precedence."""
    if opt_explicit and downstream_explicit and opt_cfg[key] != downstream_cfg[key]:
        raise click.BadParameter(
            f"opt.{key} and {downstream_section}.{key} conflict."
        )
    if opt_explicit:
        value = opt_cfg[key]
    elif downstream_explicit:
        value = downstream_cfg[key]
    else:
        value = downstream_default
    opt_cfg[key] = value
    downstream_cfg[key] = value


def _build_rsirfo_kwargs(
    rsirfo_cfg: Dict[str, Any],
    *,
    max_cycles: Optional[int],
    out_dir: Path,
    macro_thresh: Optional[str] = None,
    mode: str = "rsirfo",
    opt_cfg: Optional[Dict[str, Any]] = None,
    dump: bool = False,
    reference_mode: Optional[Any] = None,
    flatten_enabled: bool = False,
) -> Dict[str, Any]:
    # RSIRFOptimizer rejects RFOptimizer-only DIIS knobs
    # (gediis/gdiis/gdiis_thresh/gediis_thresh/gdiis_test_direction/adapt_step_func);
    # if a user-supplied YAML inherits those from a generic opt block they must be
    # stripped before construction. Centralised here so both macro/micro orchestrators
    # get identical DIIS-quirk handling.
    args = dict(rsirfo_cfg)
    # The shared `opt` block is shared: a setting the user changed there reaches
    # the TS optimizer, exactly as it reaches an ordinary optimizer. Only changed
    # values are passed, so `rsirfo.*` stays authoritative for untouched keys and
    # a default configuration is byte-identical to the pre-merge construction.
    if opt_cfg:
        args.update(strip_inherited_keys(dict(opt_cfg), OPT_BASE_KW, mode="same"))
    args["max_cycles"] = max_cycles
    args["out_dir"] = str(out_dir)
    args["dump"] = bool(dump)
    if macro_thresh is not None:
        args["thresh"] = str(macro_thresh)
    if reference_mode is not None:
        args["reference_mode"] = reference_mode
    args["flatten_enabled"] = bool(flatten_enabled)

    roots = args.get("roots")
    root_single = args.pop("root", None)
    if root_single is not None:
        roots = [int(root_single)]
    if roots is None:
        roots = [0]
    try:
        normalized_roots = [int(root) for root in roots]
    except TypeError as exc:
        raise click.BadParameter(
            "rsirfo.roots must be a list containing exactly one root index."
        ) from exc
    if len(normalized_roots) != 1:
        raise click.BadParameter(
            "rsirfo.roots must contain exactly one root index for a "
            "first-order transition-state search."
        )
    args["roots"] = normalized_roots
    args.pop("rfo_overlaps", None)

    for _diis_kw in ("gediis", "gdiis", "gdiis_thresh", "gediis_thresh", "gdiis_test_direction", "adapt_step_func"):
        args.pop(_diis_kw, None)
    if mode == "rsprfo":
        args.setdefault("min_line_search", False)
        args.setdefault("max_line_search", False)
    else:
        args.pop("min_line_search", None)
        args.pop("max_line_search", None)
    return _force_ts_reject_uphill_off(args)


def _calc_full_hessian_torch(geom, calc_kwargs: Dict[str, Any], device: torch.device) -> torch.Tensor:
    """
    Shared Hessian backend from freq.py; keeps tsopt metadata refresh behavior.
    """
    H, _ = _freq_calc_full_hessian_torch(
        geom,
        calc_kwargs,
        device,
        refresh_geom_meta=True,
    )
    return H


from mlmm.core.calc_eval import calc_energy as _calc_energy  # noqa: E402


def _omega2_to_freqs_cm(omega2: torch.Tensor) -> np.ndarray:
    """Convert eigenvalues (omega^2) to vibrational frequencies in cm^-1."""
    s_new = (units._hbar * 1e10 / np.sqrt(units._e * units._amu) * np.sqrt(AU2EV) / BOHR2ANG)
    hnu = s_new * torch.sqrt(torch.abs(omega2))
    hnu = torch.where(omega2 < 0, -hnu, hnu)
    return (hnu / units.invcm).detach().cpu().numpy()


def _clear_cuda_cache(tensor: Optional[torch.Tensor] = None) -> None:
    """Clear CUDA cache if available and tensor (if provided) is on CUDA."""
    if torch.cuda.is_available():
        if tensor is None or tensor.is_cuda:
            torch.cuda.empty_cache()




def _mw_projected_hessian_inplace(H_t: torch.Tensor,
                                  coords_bohr_t: torch.Tensor,
                                  masses_au_t: torch.Tensor,
                                  freeze_idx: Optional[List[int]] = None,
                                  tr_projection: str = "constrained",
                                  projection_info: Optional[dict] = None,
                                  compact: bool = False):
    """
    Mass-weight H in-place, optionally restrict to active DOF subspace (PHVA) and
    project out TR motions (in that subspace), also in-place.
    With ``compact=True``, return ``(H_compact, lift)`` for the true
    orthogonal complement of the rigid basis.  This keeps root selection away
    from the artificial zero eigenvectors of a same-size ``P H P`` matrix.
    """
    device = H_t.device
    with torch.no_grad():
        N = coords_bohr_t.shape[0]
        if freeze_idx:
            frozen = set(int(i) for i in freeze_idx if 0 <= int(i) < N)
            active_idx = [i for i in range(N) if i not in frozen]
            if len(active_idx) == 0:
                raise RuntimeError("All atoms are frozen; no active DOF left for TR projection.")
            # mass-weight first
            H_t = _mass_weighted_hessian(H_t, masses_au_t)
            # take active DOF submatrix
            mask_dof = torch.ones(3 * N, dtype=torch.bool, device=device)
            for i in frozen:
                mask_dof[3 * i:3 * i + 3] = False
            active_dof = torch.nonzero(mask_dof, as_tuple=False).flatten()
            H_t = active_square(H_t, active_dof)
            del active_dof
            Q, info = active_tr_basis(
                coords_bohr_t,
                masses_au_t,
                active_idx,
                mode=tr_projection,
            )
            if compact:
                H_t, lift = compact_project_hessian(H_t, Q)
            else:
                project_hessian_inplace(H_t, Q)
                lift = None
            del Q, mask_dof, active_idx, frozen
        else:
            # Full DOF: mass-weight + TR projection in-place
            H_t = _mass_weighted_hessian(H_t, masses_au_t)
            Q, info = active_tr_basis(
                coords_bohr_t,
                masses_au_t,
                list(range(int(N))),
                mode=tr_projection,
            )
            if compact:
                H_t, lift = compact_project_hessian(H_t, Q)
            else:
                project_hessian_inplace(H_t, Q)
                lift = None
            del Q
        if projection_info is not None:
            projection_info.clear()
            projection_info.update(info.as_dict())
        _clear_cuda_cache()
        return (H_t, lift) if compact else H_t


def _mode_direction_by_root(H_t: torch.Tensor,
                            coords_bohr_t: torch.Tensor,
                            masses_au_t: torch.Tensor,
                            root: int = 0,
                            freeze_idx: Optional[List[int]] = None,
                            tr_projection: str = "constrained",
                            projection_info: Optional[dict] = None) -> np.ndarray:
    """
    Get the eigenvector (Cartesian space) corresponding to the `root`-th most negative
    eigenvalue (root=0: most negative) of the mass-weighted, TR-projected Hessian.
    PHVA (active-subspace) is applied if freeze_idx is provided: frozen DOFs are zero.
    root==0 prefers torch.lobpcg; fallback to eigh (UPLO='U').
    """
    with torch.no_grad():
        # In-place: mass weight + (active-subspace) TR projection
        Hmw_proj, lift = _mw_projected_hessian_inplace(
            H_t,
            coords_bohr_t,
            masses_au_t,
            freeze_idx=freeze_idx,
            tr_projection=tr_projection,
            projection_info=projection_info,
            compact=True,
        )
        if Hmw_proj.shape[0] == 0:
            raise RuntimeError(
                "No Dimer orientation remains after rigid-null projection."
            )

        # Bounded-peak symmetrization (helper writes both triangles).
        from mlmm.core.utils import symmetrize_inplace
        symmetrize_inplace(Hmw_proj)

        # Solve eigenproblem in the (possibly reduced) space
        if int(root) == 0:
            try:
                w, v_mw_sub = torch.lobpcg(Hmw_proj, k=1, largest=False)
                u_mw_reduced = v_mw_sub[:, 0]
            except Exception:
                evals_f, evecs_f = torch.linalg.eigh(Hmw_proj, UPLO="U")
                u_mw_reduced = evecs_f[:, torch.argmin(evals_f)]
                del evals_f, evecs_f
        else:
            evals, evecs_mw = torch.linalg.eigh(Hmw_proj, UPLO="U")  # ascending
            neg = (evals < 0)
            neg_inds = torch.nonzero(neg, as_tuple=False).view(-1)
            if neg_inds.numel() == 0:
                pick = int(torch.argmin(evals).item())
            else:
                k = max(0, min(int(root), neg_inds.numel() - 1))
                pick = int(neg_inds[k].item())
            u_mw_reduced = evecs_mw[:, pick]
            del evals, evecs_mw

        u_mw_sub = (
            lift.T @ u_mw_reduced if lift is not None else u_mw_reduced
        )

        # Embed back to full 3N (frozen DOF as zeros) if we solved in subspace
        N = coords_bohr_t.shape[0]
        if freeze_idx:
            frozen = set(int(i) for i in freeze_idx if 0 <= int(i) < N)
            mask_dof = torch.ones(3 * N, dtype=torch.bool, device=Hmw_proj.device)
            for i in frozen:
                mask_dof[3 * i:3 * i + 3] = False
            u_mw_full = torch.zeros(3 * N, dtype=Hmw_proj.dtype, device=Hmw_proj.device)
            u_mw_full[mask_dof] = u_mw_sub
            u_mw = u_mw_full
            del mask_dof, frozen
        else:
            u_mw = u_mw_sub

        # Convert mass-weighted → Cartesian & normalize
        masses_amu_t = (masses_au_t / AMU2AU).to(dtype=Hmw_proj.dtype, device=Hmw_proj.device)
        m3 = torch.repeat_interleave(masses_amu_t, 3).clamp(min=1e-10)
        inv_sqrt_m = torch.sqrt(1.0 / m3)
        v = inv_sqrt_m * u_mw
        v = v / torch.linalg.norm(v)
        mode = v.reshape(-1, 3).detach().cpu().numpy()

        del masses_amu_t, m3, inv_sqrt_m, v, u_mw, u_mw_sub, u_mw_reduced
        _clear_cuda_cache()
        return mode


def _calc_gradient(geom, calc_kwargs: Dict[str, Any]) -> np.ndarray:
    """
    Return true Cartesian gradient (shape 3N,) in Hartree/Bohr.
    """
    kw = dict(calc_kwargs or {})
    kw["out_hess_torch"] = False
    calc = mlmm(**kw)
    geom.set_calculator(calc)
    g = np.array(geom.cart_gradient, dtype=float).reshape(-1)
    geom.set_calculator(None)
    del calc
    _clear_cuda_cache()
    return g


def _frequencies_cm_and_modes(H_t: torch.Tensor,
                              atomic_numbers: List[int],
                              coords_bohr: np.ndarray,
                              device: torch.device,
                              freeze_idx: Optional[List[int]] = None,
                              tr_projection: str = "constrained",
                              projection_info: Optional[dict] = None,
                              frequency_zero_cutoff_cm: float = DEFAULT_FREQUENCY_ZERO_CUTOFF_CM) -> Tuple[np.ndarray, torch.Tensor]:
    """
    In-place PHVA/TR projection (active-subspace if freeze_idx) and diagonalization.
    Returns:
      freqs_cm : (nmode,) numpy (negatives are imaginary)
      modes    : (nmode, 3N) torch (mass-weighted eigenvectors embedded to full 3N)
    """
    with torch.no_grad():
        Z = np.array(atomic_numbers, dtype=int)
        N = int(len(Z))
        masses_amu = np.array([atomic_masses[z] for z in Z])  # amu
        masses_au_t = torch.as_tensor(masses_amu * AMU2AU, dtype=H_t.dtype, device=device)
        coords_bohr_t = torch.as_tensor(coords_bohr.reshape(-1, 3), dtype=H_t.dtype, device=device)

        # in-place mass-weight + (active-subspace) TR projection, reduced to the
        # orthogonal complement of the rigid basis so every physical root survives
        Hmw, lift = _mw_projected_hessian_inplace(
            H_t,
            coords_bohr_t,
            masses_au_t,
            freeze_idx=freeze_idx,
            tr_projection=tr_projection,
            projection_info=projection_info,
            compact=True,
        )

        # Bounded-peak symmetrization (helper writes both triangles).
        from mlmm.core.utils import symmetrize_inplace
        symmetrize_inplace(Hmw)
        omega2, Vsub = torch.linalg.eigh(Hmw, UPLO="U")

        if lift is not None:
            Vsub = lift.T @ Vsub  # (3N_act or 3N, nmode)
        del lift

        # embed modes to full 3N
        if freeze_idx:
            frozen = set(int(i) for i in freeze_idx if 0 <= int(i) < N)
            mask_dof = torch.ones(3 * N, dtype=torch.bool, device=Hmw.device)
            for i in frozen:
                mask_dof[3 * i:3 * i + 3] = False
            modes = torch.zeros((Vsub.shape[1], 3 * N), dtype=Hmw.dtype, device=Hmw.device)
            modes[:, mask_dof] = Vsub.T
            del mask_dof, frozen
        else:
            modes = Vsub.T  # (nsel, 3N)

        # convert to cm^-1
        freqs_cm = _omega2_to_freqs_cm(omega2)
        freqs_cm, modes = filter_resolved_modes(
            freqs_cm, modes, frequency_zero_cutoff_cm
        )

        del omega2, Vsub, masses_amu, masses_au_t, coords_bohr_t, Hmw
        _clear_cuda_cache(H_t)
        return freqs_cm, modes


def _write_mode_trj_and_pdb(geom,
                            mode_vec_3N: np.ndarray,
                            out_trj: Path,
                            out_pdb: Path,
                            amplitude_ang: float = 0.25,
                            n_frames: int = 20,
                            comment: str = "imag mode",
                            ref_pdb: Optional[Path] = None) -> None:
    """
    Write a single imaginary mode trajectory both as _trj.xyz (XYZ-like) and .pdb.

    If `ref_pdb` is provided and is a .pdb file, the .pdb is generated by
    converting the _trj.xyz using the input PDB as the template.
    """
    ref_ang = geom.coords3d * BOHR2ANG
    mode = mode_vec_3N.reshape(-1, 3).copy()
    mode /= np.linalg.norm(mode)

    # _trj.xyz (XYZ-like concatenation) — always write
    with out_trj.open("w", encoding="utf-8") as f:
        for i in range(n_frames):
            phase = np.sin(2.0 * np.pi * i / n_frames)
            coords = ref_ang + phase * amplitude_ang * mode
            f.write(f"{len(geom.atoms)}\n{comment} frame={i+1}/{n_frames}\n")
            for sym, (x, y, z) in zip(geom.atoms, coords):
                f.write(f"{sym:2s} {x: .8f} {y: .8f} {z: .8f}\n")

    # .pdb — use ref_pdb template when available
    if ref_pdb is not None and ref_pdb.suffix.lower() == ".pdb" and is_convert_file_enabled():
        try:
            convert_xyz_to_pdb(out_trj, ref_pdb, out_pdb)
            return
        except Exception as exc:
            click.echo(
                "[convert] WARNING: mode PDB fell back to plain ASE output "
                f"without the reference topology: {exc}",
                err=True,
            )

    # Fallback: MODEL/ENDMDL via ASE (no topology)
    atoms0 = Atoms(geom.atoms, positions=ref_ang, pbc=False)
    for i in range(n_frames):
        phase = np.sin(2.0 * np.pi * i / n_frames)
        ai = atoms0.copy()
        ai.set_positions(ref_ang + phase * amplitude_ang * mode)
        write(out_pdb, ai, append=(i != 0))


def _write_all_imag_modes(
    geom,
    freqs_cm: np.ndarray,
    modes: torch.Tensor,
    neg_freq_thresh_cm: float,
    vib_dir: Path,
    *,
    ref_pdb: Optional[Path] = None,
    filename_prefix: str = "imag",
    amplitude_ang: float = 0.8,
    n_frames: int = 20,
) -> int:
    """
    Write all resolved imaginary modes to vib_dir.

    Returns:
        Number of mode trajectories written.
    """
    neg_idx = np.flatnonzero(
        resolved_imaginary_mask(freqs_cm, neg_freq_thresh_cm)
    )
    if len(neg_idx) == 0:
        return 0

    masses_amu = np.array([atomic_masses[int(z)] for z in geom.atomic_numbers], dtype=float)
    sqrt_m3 = np.sqrt(np.repeat(masses_amu, 3))
    order = np.argsort(freqs_cm[neg_idx])  # most negative first
    written = 0

    for rank, rel_i in enumerate(order, start=1):
        mode_idx = int(neg_idx[int(rel_i)])
        freq = float(freqs_cm[mode_idx])
        mode_mw = modes[mode_idx].detach().cpu().numpy().reshape(-1)
        v_cart = mode_mw / sqrt_m3
        norm = float(np.linalg.norm(v_cart))
        if norm <= 0.0:
            del mode_mw, v_cart
            continue
        v_cart = v_cart / norm

        stem = f"{filename_prefix}_{rank:02d}_{freq:+.2f}cm-1"
        out_trj = vib_dir / f"{stem}_trj.xyz"
        out_pdb = vib_dir / f"{stem}.pdb"
        _write_mode_trj_and_pdb(
            geom,
            v_cart,
            out_trj,
            out_pdb,
            amplitude_ang=amplitude_ang,
            n_frames=n_frames,
            comment=f"imag#{rank} mode={mode_idx} {freq:+.2f} cm^-1",
            ref_pdb=ref_pdb,
        )
        del mode_mw, v_cart
        written += 1

    del masses_amu, sqrt_m3, order, neg_idx
    _clear_cuda_cache()
    return written


def _certified_negative_frequencies(
    freqs_cm: np.ndarray,
    neg_freq_thresh_cm: float,
) -> List[float]:
    """Return exact-PHVA negative roots outside the shared zero window."""

    values = np.asarray(freqs_cm, dtype=float)
    return [
        float(value)
        for value in np.sort(
            values[resolved_imaginary_mask(values, neg_freq_thresh_cm)]
        )
    ]


def _certified_saddle_order(
    freqs_cm: np.ndarray,
    neg_freq_thresh_cm: float,
) -> int:
    """Count negative roots using the configured saddle threshold."""

    return len(_certified_negative_frequencies(freqs_cm, neg_freq_thresh_cm))



def _active_indices(N: int, freeze_idx: Optional[List[int]]) -> List[int]:
    if not freeze_idx:
        return list(range(N))
    fz = set(int(i) for i in freeze_idx if 0 <= int(i) < N)
    return [i for i in range(N) if i not in fz]


def _active_mask_dof(N: int, freeze_idx: Optional[List[int]]) -> np.ndarray:
    mask = np.ones(3 * N, dtype=bool)
    if freeze_idx:
        for i in freeze_idx:
            if 0 <= int(i) < N:
                mask[3 * int(i):3 * int(i) + 3] = False
    return mask


def _mask_dof_from_active_idx(N: int, active_idx: List[int]) -> np.ndarray:
    mask = np.zeros(3 * N, dtype=bool)
    for i in active_idx:
        j = int(i)
        if 0 <= j < N:
            mask[3 * j:3 * j + 3] = True
    return mask


def _extract_active_block(H_full: torch.Tensor, mask_dof: np.ndarray) -> torch.Tensor:
    """
    Return the active-DOF block as a torch.Tensor sharing device/dtype.
    """
    idx = np.flatnonzero(np.asarray(mask_dof, dtype=bool))
    return active_square(H_full, idx)


def _mw_tr_project_active_inplace(H_act: torch.Tensor,
                                  coords_full_t: torch.Tensor,
                                  masses_full_au_t: torch.Tensor,
                                  active_idx: List[int],
                                  tr_projection: str = "constrained",
                                  projection_info: Optional[dict] = None,
                                  compact: bool = False):
    """
    Mass-weight & project TR in the *active* subspace (in-place).
    """
    with torch.no_grad():
        # mass-weight
        masses_act_au_t = masses_full_au_t[active_idx]
        masses_amu_t = (masses_act_au_t / AMU2AU).to(dtype=H_act.dtype, device=H_act.device)
        m3 = torch.repeat_interleave(masses_amu_t, 3).clamp(min=1e-10)
        inv_sqrt_m_col = torch.sqrt(1.0 / m3).view(1, -1)
        inv_sqrt_m_row = inv_sqrt_m_col.view(-1, 1)
        H_act.mul_(inv_sqrt_m_row)
        H_act.mul_(inv_sqrt_m_col)
        # TR basis & projection
        Q, info = active_tr_basis(
            coords_full_t,
            masses_full_au_t,
            active_idx,
            mode=tr_projection,
        )
        if compact:
            H_act, lift = compact_project_hessian(H_act, Q)
        else:
            project_hessian_inplace(H_act, Q)
            lift = None
        if projection_info is not None:
            projection_info.clear()
            projection_info.update(info.as_dict())
        del masses_act_au_t, masses_amu_t, m3, inv_sqrt_m_col, inv_sqrt_m_row, Q
        return (H_act, lift) if compact else H_act


def _frequencies_from_Hact(H_act: torch.Tensor,
                           atomic_numbers: List[int],
                           coords_bohr: np.ndarray,
                           active_idx: List[int],
                           device: torch.device,
                           tr_projection: str = "constrained",
                           projection_info: Optional[dict] = None,
                           frequency_zero_cutoff_cm: float = DEFAULT_FREQUENCY_ZERO_CUTOFF_CM) -> np.ndarray:
    """
    Frequencies (cm^-1) computed from active-block Hessian with active-space TR projection.
    """
    with torch.no_grad():
        coords_full = torch.as_tensor(coords_bohr.reshape(-1, 3), dtype=H_act.dtype, device=device)
        masses_full_au = torch.as_tensor(
            [atomic_masses[int(z)] * AMU2AU for z in np.array(atomic_numbers, int)],
            dtype=H_act.dtype,
            device=device,
        )
        Hmw, _lift = _mw_tr_project_active_inplace(
            H_act.clone(),
            coords_full,
            masses_full_au,
            active_idx,
            tr_projection=tr_projection,
            projection_info=projection_info,
            compact=True,
        )
        # Bounded-peak symmetrization (helper writes both triangles).
        from mlmm.core.utils import symmetrize_inplace
        symmetrize_inplace(Hmw)
        omega2 = torch.linalg.eigvalsh(Hmw, UPLO="U")
        freqs_cm = _omega2_to_freqs_cm(omega2)
        freqs_cm = freqs_cm[
            resolved_frequency_mask(freqs_cm, frequency_zero_cutoff_cm)
        ]
        del coords_full, masses_full_au, Hmw, omega2, _lift
        _clear_cuda_cache(H_act)
        return freqs_cm


def _modes_from_Hact_embedded(H_act: torch.Tensor,
                              atomic_numbers: List[int],
                              coords_bohr: np.ndarray,
                              active_idx: List[int],
                              device: torch.device,
                              tr_projection: str = "constrained",
                              projection_info: Optional[dict] = None,
                              frequency_zero_cutoff_cm: float = DEFAULT_FREQUENCY_ZERO_CUTOFF_CM) -> Tuple[np.ndarray, torch.Tensor]:
    """
    Diagonalize active-block Hessian with mass-weight/TR in active space and return:
      freqs_cm : (nmode,)
      modes    : (nmode, 3N) mass-weighted eigenvectors embedded to full 3N (torch)
    """
    with torch.no_grad():
        N = len(atomic_numbers)
        coords_full = torch.as_tensor(coords_bohr.reshape(-1, 3), dtype=H_act.dtype, device=device)
        masses_full_au = torch.as_tensor(
            [atomic_masses[int(z)] * AMU2AU for z in np.array(atomic_numbers, int)],
            dtype=H_act.dtype,
            device=device,
        )
        Hmw, lift = _mw_tr_project_active_inplace(
            H_act.clone(),
            coords_full,
            masses_full_au,
            active_idx,
            tr_projection=tr_projection,
            projection_info=projection_info,
            compact=True,
        )
        # Bounded-peak symmetrization (helper writes both triangles).
        from mlmm.core.utils import symmetrize_inplace
        symmetrize_inplace(Hmw)
        omega2, Vsub = torch.linalg.eigh(Hmw, UPLO="U")
        if lift is not None:
            # Lift the reduced eigenvectors back to the active DOF space.
            Vsub = lift.T @ Vsub  # (3N_act, nmode)
        del lift

        # Embed to full 3N (mass-weighted eigenvectors)
        modes_full = torch.zeros((Vsub.shape[1], 3 * N), dtype=Hmw.dtype, device=device)
        mask_dof = _active_mask_dof(N, list(set(range(N)) - set(active_idx)))  # give frozen list
        mask_t = torch.as_tensor(mask_dof, dtype=torch.bool, device=device)
        modes_full[:, mask_t] = Vsub.T
        # frequencies
        freqs_cm = _omega2_to_freqs_cm(omega2)
        freqs_cm, modes_full = filter_resolved_modes(
            freqs_cm, modes_full, frequency_zero_cutoff_cm
        )

        del coords_full, masses_full_au, Hmw, omega2, Vsub, mask_t
        _clear_cuda_cache(H_act)
        return freqs_cm, modes_full


def _mode_direction_by_root_from_Hact(H_act: torch.Tensor,
                                      coords_bohr: np.ndarray,
                                      atomic_numbers: List[int],
                                      masses_au_t: torch.Tensor,
                                      active_idx: List[int],
                                      device: torch.device,
                                      root: int = 0,
                                      tr_projection: str = "constrained",
                                      projection_info: Optional[dict] = None) -> np.ndarray:
    """
    TS direction from the *active* Hessian block. Mass-weighting/TR are done in the
    active space. Result is embedded back to full 3N in Cartesian space.
    """
    with torch.no_grad():
        N = len(atomic_numbers)
        coords_full = torch.as_tensor(coords_bohr.reshape(-1, 3), dtype=H_act.dtype, device=device)
        masses_act_au = masses_au_t[active_idx].to(device=device, dtype=H_act.dtype)
        # mass-weight + TR in active space
        Hmw, lift = _mw_tr_project_active_inplace(
            H_act.clone(),
            coords_full,
            masses_au_t.to(device=device, dtype=H_act.dtype),
            active_idx,
            tr_projection=tr_projection,
            projection_info=projection_info,
            compact=True,
        )
        if Hmw.shape[0] == 0:
            raise RuntimeError(
                "No Dimer orientation remains after rigid-null projection."
            )
        # Bounded-peak symmetrization (helper writes both triangles).
        from mlmm.core.utils import symmetrize_inplace
        symmetrize_inplace(Hmw)

        # eigenvector for requested root
        if int(root) == 0:
            try:
                w, V = torch.lobpcg(Hmw, k=1, largest=False)
                u_mw_reduced = V[:, 0]
            except Exception:
                vals, vecs = torch.linalg.eigh(Hmw, UPLO="U")
                u_mw_reduced = vecs[:, torch.argmin(vals)]
                del vals, vecs
        else:
            vals, vecs = torch.linalg.eigh(Hmw, UPLO="U")
            neg = (vals < 0)
            neg_inds = torch.nonzero(neg, as_tuple=False).view(-1)
            if neg_inds.numel() == 0:
                pick = int(torch.argmin(vals).item())
            else:
                k = max(0, min(int(root), neg_inds.numel() - 1))
                pick = int(neg_inds[k].item())
            u_mw_reduced = vecs[:, pick]
            del vals, vecs

        u_mw = lift.T @ u_mw_reduced if lift is not None else u_mw_reduced

        # Mass un-weight to Cartesian in the active space, then embed to full 3N
        masses_act_amu = (masses_act_au / AMU2AU).to(dtype=H_act.dtype, device=device)
        m3 = torch.repeat_interleave(masses_act_amu, 3)
        v_cart_act = u_mw / torch.sqrt(m3)
        v_cart_act = v_cart_act / torch.linalg.norm(v_cart_act)

        full = torch.zeros(3 * N, dtype=H_act.dtype, device=device)
        mask_dof = _active_mask_dof(N, list(set(range(N)) - set(active_idx)))
        mask_t = torch.as_tensor(mask_dof, dtype=torch.bool, device=device)
        full[mask_t] = v_cart_act
        mode = full.reshape(-1, 3).detach().cpu().numpy()

        del coords_full, masses_act_au, masses_act_amu, m3, v_cart_act, full, mask_t, Hmw, u_mw
        _clear_cuda_cache(H_act)
        return mode


def _representative_atoms_for_mode(mode: torch.Tensor, flatten_k: int) -> np.ndarray:
    """
    Return indices of the top-k atoms with largest displacement norm in mode.
    """
    vec = mode.reshape(-1, 3)
    norms = torch.linalg.norm(vec, dim=1)
    k = min(int(flatten_k), vec.shape[0])
    if k <= 0:
        return np.zeros(0, dtype=int)
    topk = torch.topk(norms, k=k, largest=True)
    return topk.indices.detach().cpu().numpy()


def _select_flatten_targets_for_geom(
    freqs_cm: np.ndarray,
    modes: torch.Tensor,
    coords_bohr: np.ndarray,
    neg_freq_thresh_cm: float,
    root: int,
    flatten_sep_cutoff: float,
    flatten_k: int,
    primary_idx: Optional[int] = None,
) -> List[int]:
    """
    Select a subset of imaginary modes to flatten for a geometry.
    """
    neg_idx_all = np.where(freqs_cm < -abs(neg_freq_thresh_cm))[0]
    if len(neg_idx_all) <= 1:
        return []

    order = np.argsort(freqs_cm[neg_idx_all])
    sorted_neg = neg_idx_all[order]
    if primary_idx is None or int(primary_idx) not in set(map(int, sorted_neg)):
        root_clamped = max(0, min(int(root), len(order) - 1))
        primary_idx = int(sorted_neg[root_clamped])
    else:
        primary_idx = int(primary_idx)
    candidates = [int(i) for i in sorted_neg if int(i) != int(primary_idx)]
    if not candidates:
        return []

    coords_ang = torch.as_tensor(
        coords_bohr.reshape(-1, 3) * BOHR2ANG,
        dtype=modes.dtype,
        device=modes.device,
    )

    targets: List[int] = []
    reps_list: List[np.ndarray] = []

    for idx in candidates:
        rep = _representative_atoms_for_mode(modes[idx], flatten_k)
        if rep.size == 0:
            continue
        rep_coords = coords_ang[rep]
        if not reps_list:
            targets.append(idx)
            reps_list.append(rep)
            continue

        accept = True
        for prev_rep in reps_list:
            prev_coords = coords_ang[prev_rep]
            dmat = torch.cdist(rep_coords, prev_coords)
            min_dist = float(torch.min(dmat).item())
            if min_dist < float(flatten_sep_cutoff):
                accept = False
                break
        if accept:
            targets.append(idx)
            reps_list.append(rep)

    return targets


def _flatten_once_with_modes_for_geom(
    geom,
    masses_amu: np.ndarray,
    calc_kwargs: dict,
    freqs_cm: np.ndarray,
    modes: torch.Tensor,
    neg_freq_thresh_cm: float,
    flatten_amp_ang: float,
    flatten_sep_cutoff: float,
    flatten_k: int,
    root: int,
    reference_mode: Optional[np.ndarray] = None,
) -> bool:
    """
    Flatten extra imaginary modes for a geometry (single pass).
    """
    neg_idx_all = np.where(freqs_cm < -abs(neg_freq_thresh_cm))[0]
    if len(neg_idx_all) <= 1:
        return False

    primary_idx = None
    if reference_mode is not None:
        reference = np.asarray(reference_mode, dtype=float).reshape(-1)
        reference_norm = float(np.linalg.norm(reference))
        if reference.size == modes.shape[1] and reference_norm > 0.0:
            reference /= reference_norm
            m3_flat = np.repeat(masses_amu, 3)
            overlaps = []
            for idx in neg_idx_all:
                cart_mode = modes[int(idx)].detach().cpu().numpy().reshape(-1)
                cart_mode = cart_mode / np.sqrt(m3_flat)
                cart_norm = float(np.linalg.norm(cart_mode))
                overlaps.append(
                    0.0
                    if cart_norm <= 0.0
                    else abs(float(np.dot(reference, cart_mode / cart_norm)))
                )
            primary_idx = int(neg_idx_all[int(np.argmax(overlaps))])

    targets = _select_flatten_targets_for_geom(
        freqs_cm,
        modes,
        geom.cart_coords,
        neg_freq_thresh_cm,
        root,
        flatten_sep_cutoff,
        flatten_k,
        primary_idx=primary_idx,
    )
    if not targets:
        return False

    amp_bohr = float(flatten_amp_ang) / BOHR2ANG
    energy_ref = _calc_energy(geom, calc_kwargs)

    for idx in targets:
        v_mw = modes[idx].detach().cpu().numpy().reshape(-1, 3)
        m3 = np.repeat(masses_amu, 3).reshape(-1, 3)
        v_cart = v_mw / np.sqrt(m3)
        v_cart /= np.linalg.norm(v_cart)

        # Dividing the mass-weighted eigenvector by sqrt(mass) above already
        # gives its Cartesian direction. A second mass factor would rotate the
        # displacement away from the negative-curvature mode.
        disp = amp_bohr * v_cart
        ref = geom.cart_coords.reshape(-1, 3)

        plus = ref + disp
        minus = ref - disp

        _set_cartesian_flatten_coords(geom, plus)
        E_plus = _calc_energy(geom, calc_kwargs)

        _set_cartesian_flatten_coords(geom, minus)
        E_minus = _calc_energy(geom, calc_kwargs)

        use_plus = E_plus <= E_minus
        _set_cartesian_flatten_coords(geom, plus if use_plus else minus)
        energy_keep = E_plus if use_plus else E_minus
        click.echo(
            f"[Flatten] mode={idx} freq={freqs_cm[idx]:+.2f} cm^-1 "
            f"E_disp={energy_keep:.8f} Ha "
            f"ΔE={energy_keep - energy_ref:+.8f} Ha"
        )

    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return True


def _resolve_validated_hessian_analysis_atoms(
    calc_cfg: Dict[str, Any],
    n_atoms: int,
    active_dof_mode: str,
    freeze_atoms_final: List[int],
    *,
    validate_coverage: bool,
) -> List[int]:
    """Resolve the requested PHVA basis and reject known missing curvature.

    ``partial`` deliberately means ML + every movable MM atom.  A finite
    Hessian cutoff can evaluate a strict subset of that basis, so detect the
    incompatible request before an expensive TS optimization.  The exact
    final-Hessian reconciliation remains the authoritative backstop.
    """
    active_indices, layer_sets = _resolve_active_atom_indices(
        calc_cfg, n_atoms, active_dof_mode
    )
    requested = (
        set(range(int(n_atoms)))
        if active_indices is None
        else set(int(i) for i in active_indices)
    )
    requested -= set(int(i) for i in freeze_atoms_final)
    if not requested:
        raise click.ClickException(
            "Final frequency analysis requires at least one active atom."
        )

    if validate_coverage and any(layer_sets.values()):
        evaluated = set(layer_sets["ml"]) | set(layer_sets["hess_mm"])
        evaluated -= set(int(i) for i in freeze_atoms_final)
        missing = sorted(requested - evaluated)
        if missing:
            preview = ", ".join(str(i + 1) for i in missing[:12])
            suffix = " …" if len(missing) > 12 else ""
            raise click.ClickException(
                "The requested frequency-analysis basis is wider than the "
                "configured Hessian coverage; missing 1-based atom indices: "
                f"{preview}{suffix}. Increase --radius-hessian or choose a "
                "narrower --active-dof-mode."
            )
    return sorted(requested)


# The upper-triangle in-place path avoids allocating a second full-size d⊗d^T
# temporary for the Cartesian Hessian.
# CHEMISTRY-RULE:7 Bofill update advanced-indexing on active Cartesian Hessian block.
def _bofill_update_active(H_act: torch.Tensor,
                          delta_act: np.ndarray,
                          g_new_act: np.ndarray,
                          g_old_act: np.ndarray,
                          eps: float = 1e-12) -> torch.Tensor:
    """
    Memory-efficient Bofill update on the *active* Cartesian Hessian block.
    Apply symmetric rank-1/2 updates directly **in place** using only the **upper triangle**
    index set (and mirror to the lower) to avoid allocating large NxN temporaries.
    Explicit symmetrization is applied at eigendecomposition sites.
    """
    device = H_act.device
    dtype = H_act.dtype

    # as torch vectors
    d = torch.as_tensor(delta_act, dtype=dtype, device=device).reshape(-1)
    g0 = torch.as_tensor(g_old_act, dtype=dtype, device=device).reshape(-1)
    g1 = torch.as_tensor(g_new_act, dtype=dtype, device=device).reshape(-1)
    y = g1 - g0

    # Use current symmetric H_act for matvec (no extra allocation)
    Hd = H_act @ d
    xi = y - Hd

    d_dot_xi = torch.dot(d, xi)
    d_norm2 = torch.dot(d, d)
    xi_norm2 = torch.dot(xi, xi)

    # guards
    if torch.abs(d_dot_xi) > eps:
        denom_ms = d_dot_xi
    else:
        sign = torch.sign(d_dot_xi)
        denom_ms = (sign if sign != 0 else torch.tensor(1.0, device=device)) * eps
    denom_psb_d4 = d_norm2 * d_norm2 if d_norm2 > eps else eps
    denom_psb_d2 = d_norm2 if d_norm2 > eps else eps
    denom_phi = d_norm2 * xi_norm2 if (d_norm2 > eps and xi_norm2 > eps) else (1.0)

    phi = 1.0 - (d_dot_xi * d_dot_xi) / denom_phi
    phi = torch.clamp(phi, 0.0, 1.0)

    # coefficients for rank updates
    alpha = (1.0 - phi) / denom_ms                      # for xi xi^T
    beta  = - phi * (d_dot_xi / denom_psb_d4)           # for d d^T
    gamma = phi / denom_psb_d2                          # for d xi^T + xi d^T

    n = H_act.shape[0]
    iu0, iu1 = torch.triu_indices(n, n, device=device)
    is_diag = (iu0 == iu1)
    off = ~is_diag

    # Diagonal contributions (i == j): alpha*xi_i^2 + beta*d_i^2 + 2*gamma*d_i*xi_i
    if is_diag.any():
        idx = iu0[is_diag]
        diag_inc = (alpha * xi[idx] * xi[idx]
                    + beta * d[idx] * d[idx]
                    + 2.0 * gamma * d[idx] * xi[idx])
        # CHEMISTRY-RULE:7 write back by ASSIGNMENT, never `.add_`. `H_act[idx, idx]` with a
        # tensor index is advanced indexing, which returns a COPY, so an in-place add on it is
        # silently discarded and the Hessian never updates.
        H_act[idx, idx] = H_act[idx, idx] + diag_inc

    # Off-diagonal (i < j): symmetric update
    if off.any():
        i = iu0[off]
        j = iu1[off]
        inc = (alpha * xi[i] * xi[j]
               + beta * d[i] * d[j]
               + gamma * (d[i] * xi[j] + xi[i] * d[j]))
        H_act[i, j] = H_act[i, j] + inc
        H_act[j, i] = H_act[j, i] + inc

    return H_act


#                        HessianDimer (extended)

def _warn_if_leading_imaginary_mode_is_soft(ims: Any) -> None:
    """Warn when the imaginary mode that certifies the saddle is very soft.

    Certification counts imaginary modes (``n_imag == 1``); it does not assess
    their character. This only warns—the status and counting rule are unchanged.
    """
    if ims is None or len(ims) == 0:
        return
    leading = min(float(x) for x in ims)
    if abs(leading) >= TS_IMAG_SOFT_WARN_CM:
        return
    emit(
        f"[tsopt] WARNING: the leading imaginary mode is {leading:.2f} cm^-1, "
        f"below {TS_IMAG_SOFT_WARN_CM:.0f} cm^-1. Visualize the mode and "
        f"confirm IRC connectivity before treating this as a transition state.",
        narrative=True,
    )


def _tsopt_terminal_status(optimizer: Any, *, saddle_verified: bool) -> str:
    """Return numerical optimizer status independently of saddle order."""
    del saddle_verified
    if getattr(optimizer, "is_stalled", False):
        return "stalled"
    if getattr(optimizer, "is_converged", False):
        return "converged"
    return "not_converged"


def _heavy_ts_terminal_status(
    *,
    optimizer_converged: bool,
    n_imag: Optional[int],
    stalled: bool,
) -> str:
    """Return numerical heavy-optimizer status; n_imag is separate metadata."""
    del n_imag
    if stalled:
        return "stalled"
    return "converged" if optimizer_converged else "not_converged"

def _finalize_dimer_saddle_status(
    runner: Any,
    freqs_cm: np.ndarray,
    neg_freq_thresh_cm: float,
) -> np.ndarray:
    """Record the threshold-consistent final exact-Hessian verdict."""

    neg_idx = np.flatnonzero(
        resolved_imaginary_mask(freqs_cm, neg_freq_thresh_cm)
    )
    certified = _certified_negative_frequencies(freqs_cm, neg_freq_thresh_cm)
    runner.n_imaginary_modes = len(certified)
    runner.imaginary_frequencies_cm = certified
    runner.saddle_order_verified = len(certified) == 1
    if len(certified) > 1:
        click.echo(
            _unexpected_saddle_order_message(len(certified)),
            err=True,
        )
    return neg_idx


def _unexpected_saddle_order_message(n_imag: int) -> str:
    """Return the concise recovery hint for a non-first-order result."""

    if n_imag == 0:
        return "[tsopt] No imaginary mode detected. Try all --refine-path."
    if n_imag > 1:
        return (
            f"[tsopt] WARNING: Higher-order stationary point (n_imag={n_imag}). "
            "Try --flatten or all --refine-path."
        )
    raise ValueError("n_imag must differ from 1")


def _dimer_mode_export_message(
    n_written: int,
    n_imag: int,
    threshold_cm: float,
    min_frequency_cm: float,
) -> tuple[str, bool]:
    """Return the final mode-export message and whether it is diagnostic."""

    del threshold_cm, min_frequency_cm
    if n_written:
        return f"[tsopt] Wrote {n_written} final imaginary mode(s).", False
    if n_imag == 0:
        return _unexpected_saddle_order_message(n_imag), True
    return (
        "[tsopt] ERROR: Failed to write imaginary mode trajectory.",
        True,
    )


class HessianDimer:
    """
    Dimer-based TS search with periodic Hessian updates.

    Extensions in this implementation:
      - `root` parameter: choose which imaginary mode to follow (0 = most negative).
      - Pass-through kwargs: `dimer_kwargs` and `lbfgs_kwargs` to tune internals.
      - Optional cap on total LBFGS steps across segments: `max_total_cycles`.
      - PHVA (active DOF subspace) + TR projection for mode picking,
        respecting ``freeze_atoms``, with in-place operations. Root 0
        unconditionally uses LOBPCG with the existing dense fallback.
      - The flatten loop uses a *Bofill*-updated active Hessian block, so the
        expensive exact Hessian is evaluated only once before the flatten loop and
        once at the end for the final frequency analysis.
      - Calculator kwargs accept ``freeze_atoms`` and ``hessian_calc_mode`` and
        default to ``return_partial_hessian=True`` (active-block Hessian when frozen).
    """

    def __init__(self,
                 fn: str,
                 out_dir: str = OUT_DIR_TSOPT,
                 thresh_loose: str = "gau_loose",
                 thresh: str = "baker",
                 update_interval_hessian: int = 500,
                 # neg_freq_thresh_cm selects modes for trajectory output, flattening,
                 # and recovery. Saddle-order certification counts every
                 # negative root of the exact compact PHVA spectrum.
                 neg_freq_thresh_cm: float = 5.0,
                 flatten_amp_ang: float = 0.10,
                 flatten_max_iter: int = 50,
                 mem: int = 100000,
                 use_lobpcg: bool = True,  # compatibility no-op; root 0 always uses LOBPCG
                 calc_kwargs: Optional[dict] = None,
                 device: str = "auto",
                 dump: bool = False,
                 #
                 # New:
                 root: int = 0,
                 dimer_kwargs: Optional[Dict[str, Any]] = None,
                 lbfgs_kwargs: Optional[Dict[str, Any]] = None,
                 max_total_cycles: Optional[int] = None,
                 #
                # Pass geom kwargs so freeze-atoms and YAML geometry overrides apply on the light path
                 geom_kwargs: Optional[Dict[str, Any]] = None,
                 # New: Use partial Hessian for imaginary mode detection in flatten loop
                 partial_hessian_flatten: bool = True,
                 # Spatial separation for flatten mode selection
                 flatten_sep_cutoff: float = 0.0,
                 flatten_k: int = 10,
                 flatten_loop_bofill: bool = False,
                 ml_only_hessian_dimer: bool = False,
                 analysis_active_atoms: Optional[List[int]] = None,
                 source_path: Optional[Path] = None,
                 skip_final_freq: bool = False,
                 ) -> None:

        update_interval_hessian = int(update_interval_hessian)
        if update_interval_hessian < 1:
            raise ValueError("update_interval_hessian must be at least 1")

        self.fn = fn
        self.source_path = Path(source_path) if source_path is not None else None
        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.vib_dir = self.out_dir / "vib"
        self.vib_dir.mkdir(parents=True, exist_ok=True)

        self.thresh_loose = thresh_loose
        self.thresh = thresh
        self.update_interval_hessian = update_interval_hessian
        self.neg_freq_thresh_cm = normalize_frequency_zero_cutoff_cm(
            neg_freq_thresh_cm
        )
        self.flatten_amp_ang = float(flatten_amp_ang)
        self.flatten_max_iter = int(flatten_max_iter)
        self.mem = int(mem)
        # Retain the public attribute for compatibility; mode selection ignores
        # it and preserves the established root-0 LOBPCG-with-eigh-fallback path.
        self.use_lobpcg = bool(use_lobpcg)
        self.root = int(root)
        self.dimer_kwargs = dict(dimer_kwargs or {})
        self.lbfgs_kwargs = dict(lbfgs_kwargs or {})
        self.max_total_cycles = (
            None if max_total_cycles is None else int(max_total_cycles)
        )
        self.partial_hessian_flatten = bool(partial_hessian_flatten)
        # Spatial separation for flatten mode selection
        self.flatten_sep_cutoff = float(flatten_sep_cutoff)
        self.flatten_k = int(flatten_k)
        self.flatten_loop_bofill = bool(flatten_loop_bofill)
        self.ml_only_hessian_dimer = bool(ml_only_hessian_dimer)
        self.skip_final_freq = bool(skip_final_freq)

        # Total cycles across all flatten + opt loops/segments.
        self._cycles_spent = 0

        # Honest convergence state of the last dimer loop that ran: True iff that
        # loop ended because the optimizer reported convergence, False iff it ended
        # on global cycle-budget exhaustion. Threaded out to the result.json status
        # so the dimer (grad) branch reports converged/not_converged like the
        # RS-I-RFO branch instead of a neutral "completed" literal.
        self.is_converged = False

        # Additive stall state propagated from a child LBFGS whose
        # energy plateaued (energy stopped decreasing while its force/step
        # criteria stayed unmet).  A stall stops all later segments/loops and
        # is never reported as a converged TS.
        self.is_stalled = False
        self.stop_reason = ""
        self.flatten_skip_reason: Optional[str] = None
        self.saddle_order_verified = False
        self.n_imaginary_modes: Optional[int] = None
        self.imaginary_frequencies_cm: List[float] = []
        self.hessian_status = "not_run"
        self.hessian_error: Optional[str] = None

        # Hessian caching for 0-step convergence (avoid redundant recalculation)
        self._raw_hessian_cache_cpu: Optional[torch.Tensor] = None
        self._raw_hessian_coords_cpu: Optional[np.ndarray] = None
        self._raw_hessian_identity: Optional[Dict[str, Any]] = None
        self._last_active_idx: Optional[List[int]] = None
        self._last_active_mask_dof: Optional[np.ndarray] = None

        # ML/MM calculator settings
        self.calc_kwargs = dict(calc_kwargs or {})
        self.calc_kwargs.setdefault("out_hess_torch", False)

        # Geometry & masses (use provided geom kwargs so freeze_atoms etc. apply)
        gkw = dict(geom_kwargs or {})
        coord_type = str(gkw.pop("coord_type", "cart")).lower()
        if coord_type != "cart":
            raise ValueError(
                "HessianDimer uses Cartesian 3N Hessian, mode, and Bofill kernels; "
                "coord_type must be 'cart'."
            )
        freeze_geom = list(gkw.get("freeze_atoms", [])) if "freeze_atoms" in gkw else []
        freeze_calc_raw = self.calc_kwargs.get("freeze_atoms") or []
        try:
            freeze_calc = [int(i) for i in freeze_calc_raw]
        except TypeError:
            freeze_calc = [int(freeze_calc_raw)]
        merged_freeze = sorted({int(i) for i in (freeze_geom + freeze_calc)})
        if merged_freeze:
            gkw["freeze_atoms"] = merged_freeze
        elif "freeze_atoms" in gkw:
            gkw["freeze_atoms"] = []
        self.calc_kwargs["freeze_atoms"] = merged_freeze

        self.calc_kwargs_partial = dict(self.calc_kwargs)
        self.calc_kwargs_partial["mm_hessian_mode"] = "none"
        self.calc_kwargs_partial["return_partial_hessian"] = False
        self.calc_kwargs_partial["out_hess_torch"] = True
        self.calc_kwargs_full = dict(self.calc_kwargs)
        self.calc_kwargs_full.setdefault("mm_fd", True)
        self.calc_kwargs_full["return_partial_hessian"] = False
        self.calc_kwargs_full["out_hess_torch"] = True
        # ML-only Hessian kwargs: skip MM Hessian entirely, use ML partial Hessian only
        self.calc_kwargs_ml_only = dict(self.calc_kwargs)
        self.calc_kwargs_ml_only["mm_hessian_mode"] = "none"
        self.calc_kwargs_ml_only["return_partial_hessian"] = True
        self.calc_kwargs_ml_only["out_hess_torch"] = True
        self.calc_kwargs_ml_only["hess_cutoff"] = 0.0  # ML atoms only in Hessian
        self.geom = geom_loader(fn, coord_type=coord_type, **gkw)
        self.tr_projection = self.geom.tr_projection
        self.rigid_projection_info: Dict[str, Any] = {}
        # If partial Hessian is requested (explicitly or via B-factor layers),
        # avoid full 3N Hessian allocations in light TS dimer runs.
        if self.calc_kwargs.get("return_partial_hessian"):
            self.calc_kwargs_partial["return_partial_hessian"] = True
            self.calc_kwargs_full["return_partial_hessian"] = True
        elif self.partial_hessian_flatten and self.calc_kwargs.get("use_bfactor_layers"):
            self.calc_kwargs_partial["return_partial_hessian"] = True
            self.calc_kwargs_full["return_partial_hessian"] = True
        self.masses_amu = np.array([atomic_masses[z] for z in self.geom.atomic_numbers])
        self.masses_au_t = torch.as_tensor(self.masses_amu * AMU2AU, dtype=torch.float32)

        # --- Preserve freeze list (for PHVA) ---
        self.freeze_atoms: List[int] = [int(i) for i in self.geom.freeze_atoms]
        requested = (
            set(range(len(self.geom.atomic_numbers)))
            if analysis_active_atoms is None
            else set(int(i) for i in analysis_active_atoms)
        )
        requested -= set(self.freeze_atoms)
        if not requested:
            raise ValueError("Final frequency analysis requires at least one active atom.")
        if min(requested) < 0 or max(requested) >= len(self.geom.atomic_numbers):
            raise ValueError("analysis_active_atoms contains an out-of-range atom index.")
        self.analysis_active_atoms = sorted(requested)

        # Device
        self.device = _torch_device(device)
        self.masses_au_t = self.masses_au_t.to(self.device)

        # temp file for Dimer orientation (N_raw)
        self.mode_path = self.out_dir / ".dimer_mode.dat"

        self.dump = bool(dump)
        self.optim_all_path = self.out_dir / "optimization_all_trj.xyz"

    @property
    def termination_status(self) -> str:
        """Public terminal outcome: ``stalled`` > ``converged`` > ``not_converged``."""
        if self.is_stalled:
            return "stalled"
        if self.is_converged:
            return "converged"
        return "not_converged"

    # ----- One dimer segment for up to n_steps; returns (steps_done, converged) -----
    def _dimer_segment(self, threshold: str, n_steps: int) -> Tuple[int, bool]:
        # Dimer calculator using current mode as initial N
        calc_sp = mlmm(**self.calc_kwargs)

        n_atoms = len(self.geom.atomic_numbers)
        active_idx = _active_indices(n_atoms, self.freeze_atoms)
        rigid_masses = self.masses_au_t.detach().to(device="cpu", dtype=torch.float64)
        projection_mode = self.tr_projection
        rigid_basis, projection = full_cartesian_tr_basis(
            torch.as_tensor(self.geom.cart_coords.reshape(-1, 3), dtype=torch.float64),
            rigid_masses,
            active_idx,
            mode=projection_mode,
        )

        def rigid_basis_at(coords_flat: np.ndarray) -> np.ndarray:
            basis, _ = full_cartesian_tr_basis(
                torch.as_tensor(coords_flat.reshape(-1, 3), dtype=torch.float64),
                rigid_masses,
                active_idx,
                mode=projection_mode,
            )
            return basis.detach().cpu().numpy()

        self.rigid_projection_info.clear()
        self.rigid_projection_info.update(projection.as_dict())

        # Merge user dimer kwargs, while enforcing the geometry's exact Cartesian
        # constraints and its selected rigid-null treatment.
        dimer_kwargs = dict(self.dimer_kwargs)
        dimer_kwargs.update({
            "calculator": calc_sp,
            "N_raw": str(self.mode_path),
            "frozen_atoms": self.freeze_atoms,
            "rigid_basis": rigid_basis.detach().cpu().numpy(),
            "rigid_basis_getter": rigid_basis_at,
            "seed": 0,                    # runner override for determinism
            "mem": self.mem,              # accepted by Calculator base through **kwargs
            "out_dir": str(self.out_dir),
        })
        dimer = Dimer(**dimer_kwargs)

        self.geom.set_calculator(dimer)

        # LBFGS kwargs: enforce thresh/max_cycles/out_dir/dump; allow others
        lbfgs_kwargs = _force_ts_reject_uphill_off(self.lbfgs_kwargs)
        lbfgs_kwargs.update({
            "max_cycles": n_steps,
            "thresh": threshold,
            "out_dir": str(self.out_dir),
            "dump": self.dump,
        })
        opt = LBFGS(self.geom, **lbfgs_kwargs)
        opt.run()
        # pysisyphus uses 0-indexed cur_cycle; keep budget accounting strict by clamping
        # to the requested segment step count.
        steps = min(max(int(opt.cur_cycle) + 1, 1), int(n_steps))
        converged = opt.is_converged
        # Propagate an energy-plateau stall from the child LBFGS . A
        # stalled child is not converged; the caller stops all later segments.
        if getattr(opt, "is_stalled", False):
            self.is_stalled = True
            self.stop_reason = getattr(opt, "stop_reason", "") or self.stop_reason
        self.geom.set_calculator(None)

        # Free dimer/optimizer GPU resources before next Hessian computation
        del calc_sp, dimer, opt
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        # Append to concatenated trajectory if dump enabled
        if self.dump:
            _append_xyz_trajectory(self.optim_all_path, self.out_dir / "optimization_trj.xyz")
        return steps, converged

    # ----- Hessian caching for 0-step convergence -----
    def _cache_raw_hessian_cpu(
        self,
        H: torch.Tensor,
        calc_kwargs: Dict[str, Any],
    ) -> None:
        """Cache a raw Hessian with its geometry and evaluation identity."""
        from mlmm.io.hessian_cache import persistent_identity_from_context

        self._raw_hessian_cache_cpu = H.detach().cpu().clone()
        self._raw_hessian_coords_cpu = self.geom.cart_coords.copy()
        self._raw_hessian_identity = persistent_identity_from_context(
            self.geom,
            calc_kwargs,
            method="hessian_dimer_raw",
        )

    def _reuse_cached_hessian(
        self,
        calc_kwargs: Dict[str, Any],
    ) -> Optional[torch.Tensor]:
        """Return a cache hit only for identical coordinates and evaluator."""
        from mlmm.io.hessian_cache import persistent_identity_from_context

        if (
            self._raw_hessian_cache_cpu is None
            or self._raw_hessian_coords_cpu is None
            or self._raw_hessian_identity is None
        ):
            return None
        if not np.array_equal(self.geom.cart_coords, self._raw_hessian_coords_cpu):
            return None
        expected = persistent_identity_from_context(
            self.geom,
            calc_kwargs,
            method="hessian_dimer_raw",
        )
        if expected != self._raw_hessian_identity:
            return None
        H_dev = self._raw_hessian_cache_cpu.to(self.device)
        if self.device.type == "cpu":
            H_dev = H_dev.clone()
        return H_dev

    def _calc_full_hessian_cached(
        self, calc_kwargs: Dict[str, Any], allow_reuse: bool
    ) -> torch.Tensor:
        """Compute Hessian, caching on CPU. Reuse if allow_reuse and geometry unchanged."""
        if allow_reuse:
            cached = self._reuse_cached_hessian(calc_kwargs)
            if cached is not None:
                click.echo("[tsopt] Reusing cached raw Hessian (0-step convergence).")
                return cached
        H = _calc_full_hessian_torch(self.geom, calc_kwargs, self.device)
        H = self._compact_hessian_to_computed_coverage(H)
        self._cache_raw_hessian_cpu(H, calc_kwargs)
        return H

    def _compact_hessian_to_computed_coverage(
        self, hessian: torch.Tensor
    ) -> torch.Tensor:
        """Remove full-shape zero padding without inventing missing curvature."""
        coverage = _ordered_hessian_coverage_atoms(
            self.geom, len(self.geom.atomic_numbers)
        )
        if coverage is None:
            return hessian
        compact, _, _, _ = _reconcile_hessian_analysis_basis(
            hessian, self.geom, coverage
        )
        del hessian
        return compact

    def _resolve_hessian_active_subspace(self, H_t: torch.Tensor, N: int) -> Tuple[List[int], np.ndarray]:
        """
        Resolve active atoms/DOFs for a Hessian tensor.
        For partial Hessians, prefer geometry metadata populated by the calculator.
        """
        h_dim = int(H_t.size(0))
        full_dim = 3 * int(N)
        freeze = self.freeze_atoms if len(self.freeze_atoms) > 0 else []

        if h_dim == full_dim:
            active_idx = _active_indices(N, freeze)
            mask_dof = _active_mask_dof(N, freeze)
            self._last_active_idx = list(active_idx)
            self._last_active_mask_dof = mask_dof.copy()
            return active_idx, mask_dof

        def _norm_atoms(vals: Optional[Any]) -> np.ndarray:
            if vals is None:
                return np.zeros(0, dtype=int)
            arr = np.asarray(vals, dtype=int).reshape(-1)
            return arr[(arr >= 0) & (arr < N)]

        def _norm_dofs(vals: Optional[Any]) -> np.ndarray:
            if vals is None:
                return np.zeros(0, dtype=int)
            arr = np.asarray(vals, dtype=int).reshape(-1)
            return arr[(arr >= 0) & (arr < full_dim)]

        def _stable_unique(vals: np.ndarray) -> np.ndarray:
            seen = set()
            out: List[int] = []
            for v in vals.tolist():
                iv = int(v)
                if iv not in seen:
                    seen.add(iv)
                    out.append(iv)
            return np.asarray(out, dtype=int)

        candidates: List[Tuple[str, np.ndarray, np.ndarray]] = []
        try:
            candidates.append((
                "geom.hess_active_*",
                _norm_atoms(self.geom.hess_active_atom_indices),
                _norm_dofs(self.geom.hess_active_dof_indices),
            ))
        except Exception:
            logger.debug("Failed to read hess_active_* indices", exc_info=True)

        within = getattr(self.geom, "within_partial_hessian", None)
        if isinstance(within, dict):
            candidates.append((
                "geom.within_partial_hessian",
                _norm_atoms(within.get("active_atoms")),
                _norm_dofs(within.get("active_dofs")),
            ))

        candidates.append((
            "geom._hess_active_*_last",
            _norm_atoms(getattr(self.geom, "_hess_active_atoms_last", None)),
            _norm_dofs(getattr(self.geom, "_hess_active_dofs_last", None)),
        ))

        if self._last_active_idx is not None or self._last_active_mask_dof is not None:
            cached_atoms = _norm_atoms(self._last_active_idx)
            cached_dofs = np.flatnonzero(self._last_active_mask_dof).astype(int) \
                if self._last_active_mask_dof is not None else np.zeros(0, dtype=int)
            candidates.append(("cached_active_subspace", cached_atoms, _norm_dofs(cached_dofs)))

        fallback_atoms = _norm_atoms(_active_indices(N, freeze))
        fallback_dofs = np.flatnonzero(_active_mask_dof(N, freeze)).astype(int)
        candidates.append(("freeze_based", fallback_atoms, _norm_dofs(fallback_dofs)))

        for _, atoms_arr, dofs_arr in candidates:
            if dofs_arr.size > 0:
                mask_dof = np.zeros(full_dim, dtype=bool)
                mask_dof[dofs_arr] = True
            elif atoms_arr.size > 0:
                mask_dof = _mask_dof_from_active_idx(N, atoms_arr.tolist())
            else:
                continue

            if int(mask_dof.sum()) != h_dim:
                continue

            if dofs_arr.size > 0:
                atoms_arr = _stable_unique((dofs_arr // 3).astype(int))
            elif atoms_arr.size == 0:
                atoms_arr = _stable_unique((np.flatnonzero(mask_dof) // 3).astype(int))
            active_idx = [int(i) for i in atoms_arr.tolist()]
            self._last_active_idx = list(active_idx)
            self._last_active_mask_dof = mask_dof.copy()
            return active_idx, mask_dof

        raise RuntimeError(
            f"Failed to resolve active subspace for partial Hessian: "
            f"H_dim={h_dim}, full_dim={full_dim}, freeze_active_dof={int(fallback_dofs.size)}"
        )

    # ----- Loop dimer segments, updating mode from Hessian every interval -----
    def _dimer_loop(
        self,
        threshold: str,
        *,
        reserve_cycles: int = 0,
    ) -> Tuple[int, bool, bool]:
        """
        Run multiple LBFGS segments separated by periodic Hessian-based mode updates.
        Consumes from a *global* cycle budget self.max_total_cycles.

        Returns:
            (steps_in_this_call, zero_step_converged, loop_converged)
        where `zero_step_converged` is True iff the loop terminated by convergence
        without changing the geometry (i.e., 0-step convergence; a Hessian-reuse
        hint), and `loop_converged` is True iff the loop terminated because the
        underlying optimizer reported convergence rather than the global cycle
        budget being exhausted. `loop_converged` is the honest convergence signal
        threaded out to the result.json status.
        """
        steps_in_this_call = 0
        zero_step_converged = False
        loop_converged = False
        while True:
            remaining_global = (
                None
                if self.max_total_cycles is None
                else max(
                    0,
                    self.max_total_cycles
                    - self._cycles_spent
                    - int(reserve_cycles),
                )
            )
            if remaining_global == 0:
                break
            steps_this = min(self.update_interval_hessian, remaining_global)
            coords_before = self.geom.cart_coords.copy()
            steps, ok = self._dimer_segment(threshold, steps_this)
            self._cycles_spent += steps
            steps_in_this_call += steps
            # A stalled child stops all further segments in this loop.
            if self.is_stalled:
                break
            if ok:
                loop_converged = True
                # Check if geometry unchanged (0-step convergence)
                if np.array_equal(self.geom.cart_coords, coords_before):
                    zero_step_converged = True
                break
            # If budget exhausted after this segment, stop before doing a Hessian update
            if (
                self.max_total_cycles is not None
                and (self.max_total_cycles - self._cycles_spent) <= 0
            ):
                break
            # Update mode from Hessian (respect freeze atoms via PHVA)
            # Ensure VRAM is fully released after dimer segment before heavy Hessian computation
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            # Choose ML-only or full active-DOF Hessian for mode direction
            hess_kw = self.calc_kwargs_ml_only if self.ml_only_hessian_dimer else self.calc_kwargs_partial
            H_t = _calc_full_hessian_torch(self.geom, hess_kw, self.device)
            H_t = self._compact_hessian_to_computed_coverage(H_t)
            N = len(self.geom.atomic_numbers)
            coords_bohr_t = torch.as_tensor(self.geom.cart_coords.reshape(-1, 3),
                                            dtype=H_t.dtype, device=H_t.device)
            # full vs active-block Hessian
            if H_t.size(0) == 3 * N:
                mode_xyz = _mode_direction_by_root(
                    H_t, coords_bohr_t, self.masses_au_t,
                    root=self.root,
                    freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None,
                    tr_projection=self.tr_projection,
                    projection_info=self.rigid_projection_info,
                )
            else:
                # partial (active) Hessian returned by UMA
                active_idx, _ = self._resolve_hessian_active_subspace(H_t, N)
                mode_xyz = _mode_direction_by_root_from_Hact(
                    H_t, self.geom.cart_coords.reshape(-1, 3), self.geom.atomic_numbers,
                    self.masses_au_t, active_idx, self.device, root=self.root,
                    tr_projection=self.tr_projection,
                    projection_info=self.rigid_projection_info,
                )
            np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")
            del H_t, coords_bohr_t, mode_xyz
            _clear_cuda_cache()
        return steps_in_this_call, zero_step_converged, loop_converged

    def _flatten_once_with_modes(self, freqs_cm: np.ndarray, modes: torch.Tensor) -> bool:
        """
        Flatten using precomputed (approximate) modes (mass-weighted, embedded).

        Uses spatial separation (if flatten_sep_cutoff > 0) to select only modes
        whose representative atoms are well-separated from each other. This avoids
        applying conflicting displacements to nearby regions. Modes are applied
        sequentially, updating the reference position after each mode.
        """
        neg_idx_all = np.where(freqs_cm < -abs(self.neg_freq_thresh_cm))[0]
        if len(neg_idx_all) <= 1:
            return False

        # Use spatial separation if cutoff > 0, otherwise select all non-primary modes
        if self.flatten_sep_cutoff > 0:
            targets = _select_flatten_targets_for_geom(
                freqs_cm,
                modes,
                self.geom.cart_coords,
                self.neg_freq_thresh_cm,
                self.root,
                self.flatten_sep_cutoff,
                self.flatten_k,
            )
        else:
            # Legacy behavior: select all imaginary modes except primary
            order = np.argsort(freqs_cm[neg_idx_all])
            root_clamped = max(0, min(self.root, len(order) - 1))
            primary_idx = neg_idx_all[order[root_clamped]]
            targets = [i for i in neg_idx_all if i != primary_idx]

        if not targets:
            return False

        amp_bohr = self.flatten_amp_ang / BOHR2ANG

        # Get reference energy
        E_ref = _calc_energy(self.geom, self.calc_kwargs)

        # Apply modes sequentially
        for idx in targets:
            v_mw = modes[idx].detach().cpu().numpy().reshape(-1, 3)
            m3 = np.repeat(self.masses_amu, 3).reshape(-1, 3)
            v_cart = v_mw / np.sqrt(m3)
            v_cart /= np.linalg.norm(v_cart)

            disp = amp_bohr * v_cart
            ref = self.geom.cart_coords.reshape(-1, 3)

            plus = ref + disp
            minus = ref - disp

            _set_cartesian_flatten_coords(self.geom, plus)
            E_plus = _calc_energy(self.geom, self.calc_kwargs)

            _set_cartesian_flatten_coords(self.geom, minus)
            E_minus = _calc_energy(self.geom, self.calc_kwargs)

            # Keep lower-energy side and continue from there
            use_plus = E_plus <= E_minus
            _set_cartesian_flatten_coords(
                self.geom, plus if use_plus else minus
            )
            E_keep = E_plus if use_plus else E_minus
            delta_e = E_keep - E_ref
            click.echo(
                f"[Flatten] mode={idx} freq={freqs_cm[idx]:+.2f} cm^-1 "
                f"E_disp={E_keep:.8f} Ha ΔE={delta_e:+.8f} Ha"
            )

        _clear_cuda_cache()
        return True

    # ----- Run full procedure -----
    def run(self) -> None:
        if self.dump and self.optim_all_path.exists():
            self.optim_all_path.unlink()

        N = len(self.geom.atomic_numbers)
        H_final_reuse_cpu: Optional[torch.Tensor] = None
        H_final_reuse_coords: Optional[np.ndarray] = None

        # (1) Initial Hessian → pick direction by `root`
        hess_kw_init = self.calc_kwargs_ml_only if self.ml_only_hessian_dimer else self.calc_kwargs_partial
        if self.ml_only_hessian_dimer:
            click.echo("[tsopt] Using ML-only Hessian for dimer orientation.")
        H_t = _calc_full_hessian_torch(self.geom, hess_kw_init, self.device)
        H_t = self._compact_hessian_to_computed_coverage(H_t)
        coords_bohr_t = torch.as_tensor(self.geom.cart_coords.reshape(-1, 3),
                                        dtype=H_t.dtype, device=H_t.device)
        active_idx, mask_dof = self._resolve_hessian_active_subspace(H_t, N)
        if H_t.size(0) != 3 * N:
            click.echo(
                f"[tsopt] H_act={int(H_t.size(0))} active_atoms={len(active_idx)} "
                f"active_dofs={int(mask_dof.sum())} within={self.geom.within_partial_hessian is not None}"
            )

        if H_t.size(0) == 3 * N:
            # Skip heavy TR-projection residual check to conserve VRAM.
            click.echo("[tsopt] TR-projection residual check skipped to conserve VRAM.")
            mode_xyz = _mode_direction_by_root(
                H_t, coords_bohr_t, self.masses_au_t,
                root=self.root,
                freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None,
                tr_projection=self.tr_projection,
                projection_info=self.rigid_projection_info,
            )
        else:
            click.echo("[tsopt] Using active-block Hessian from UMA (partial Hessian). Skip full-space TR check.")
            mode_xyz = _mode_direction_by_root_from_Hact(
                H_t, self.geom.cart_coords.reshape(-1, 3), self.geom.atomic_numbers,
                self.masses_au_t, active_idx, self.device, root=self.root,
                tr_projection=self.tr_projection,
                projection_info=self.rigid_projection_info,
            )
        np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")
        del mode_xyz, coords_bohr_t, H_t
        _clear_cuda_cache()

        # (2) Loose loop
        if self.root!=0:
            click.echo("[tsopt] root != 0. Use this 'root' in first dimer loop", err=True)
            click.echo(f"[tsopt] Dimer Loop with initial direction from mode {self.root}...")
            self.root=0
            self.thresh_loose = self.thresh
        else:
            click.echo("[tsopt] Loose Dimer Loop...")

        thresholds_match = self.thresh_loose == self.thresh
        strict_reserve = 0 if thresholds_match else 1
        _, zero_step_loose, conv_loose = self._dimer_loop(
            self.thresh_loose,
            reserve_cycles=strict_reserve,
        )
        # A loose-threshold pass is phase progress, not terminal convergence.
        self.is_converged = bool(conv_loose and thresholds_match)

        zero_step_normal = False
        # A stalled loose loop stops all further optimization work :
        # skip the Hessian/mode update and the normal + flatten loops so a
        # stalled TS search is never retried.
        if self.is_stalled:
            click.echo("[tsopt] Optimization stalled (energy plateau); skipping the normal dimer loop.")
        elif thresholds_match and conv_loose:
            click.echo("[tsopt] Loose and final thresholds are identical; strict pass is complete.")
        elif (
            self.max_total_cycles is None
            or (self.max_total_cycles - self._cycles_spent) > 0
        ):
            # (3) Update mode & normal loop (reuse Hessian if 0-step converged)
            H_t = self._calc_full_hessian_cached(self.calc_kwargs_partial, allow_reuse=zero_step_loose)
            coords_bohr_t = torch.as_tensor(self.geom.cart_coords.reshape(-1, 3),
                                            dtype=H_t.dtype, device=H_t.device)
            if H_t.size(0) == 3 * N:
                click.echo("[tsopt] TR-projection residual check skipped to conserve VRAM.")
                mode_xyz = _mode_direction_by_root(
                    H_t, coords_bohr_t, self.masses_au_t,
                    root=self.root,
                    freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None,
                    tr_projection=self.tr_projection,
                    projection_info=self.rigid_projection_info,
                )
            else:
                click.echo("[tsopt] Using active-block Hessian from UMA (partial Hessian). Skip full-space TR check.")
                active_idx, mask_dof = self._resolve_hessian_active_subspace(H_t, N)
                mode_xyz = _mode_direction_by_root_from_Hact(
                    H_t, self.geom.cart_coords.reshape(-1, 3), self.geom.atomic_numbers,
                    self.masses_au_t, active_idx, self.device, root=self.root,
                    tr_projection=self.tr_projection,
                    projection_info=self.rigid_projection_info,
                )
            np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")
            del mode_xyz, coords_bohr_t, H_t
            _clear_cuda_cache()

            click.echo("[tsopt] Normal Dimer Loop...")
            _, zero_step_normal, conv_normal = self._dimer_loop(self.thresh)
            self.is_converged = conv_normal
        else:
            click.echo("[tsopt] Reached --max-cycles budget after loose loop; skipping normal dimer loop.")

        # A stalled optimization never enters the flatten/retry loop.
        if self.flatten_max_iter > 0 and self.is_stalled:
            self.flatten_skip_reason = "optimization stalled before flattening"
            click.echo("[tsopt] Optimization stalled (energy plateau); skipping the flatten loop.")
        elif self.flatten_max_iter > 0 and (
            self.max_total_cycles is None
            or (self.max_total_cycles - self._cycles_spent) > 0
        ):
            # (4) Flatten loop.  Bofill is an explicit approximation policy;
            # with the default off setting, refresh the exact Hessian after
            # each dimer segment.
            if self.flatten_loop_bofill:
                click.echo("[tsopt] Flatten loop with Bofill-updated active Hessian...")
            else:
                click.echo("[tsopt] Flatten loop with exact Hessian refreshes...")

            # (4.1) Evaluate one exact Hessian at the loop start and prepare the active block
            # (reuse Hessian if 0-step converged)
            H_any = self._calc_full_hessian_cached(self.calc_kwargs_full, allow_reuse=zero_step_normal)
            # Keep a CPU copy so we can skip the final Hessian recomputation
            # when the flatten loop leaves geometry unchanged.
            H_final_reuse_cpu = H_any.detach().cpu().clone()
            H_final_reuse_coords = self.geom.cart_coords.copy()
            if H_any.size(0) == 3 * N:
                # full → extract active
                H_act = _extract_active_block(H_any, mask_dof)  # torch (3N_act,3N_act)
            else:
                # UMA already returned active-block Hessian
                active_idx, mask_dof = self._resolve_hessian_active_subspace(H_any, N)
                H_act = H_any
            del H_any
            _clear_cuda_cache()

            # Gradient snapshots are only needed by the optional quasi-Newton
            # path; avoid both the calculation and retained arrays otherwise.
            g_prev = (
                _calc_gradient(self.geom, self.calc_kwargs).reshape(-1)
                if self.flatten_loop_bofill
                else None
            )

            # Flatten iterations with *approximate* Hessian updates
            for _it in range(self.flatten_max_iter):
                if (
                    self.max_total_cycles is not None
                    and (self.max_total_cycles - self._cycles_spent) <= 0
                ):
                    self.flatten_skip_reason = (
                        "max-cycles budget exhausted during flattening"
                    )
                    break

                # (a) Estimate current imaginary modes using the *active* Hessian
                freqs_est = _frequencies_from_Hact(H_act, self.geom.atomic_numbers,
                                                   self.geom.cart_coords.reshape(-1, 3), active_idx, self.device,
                                                   tr_projection=self.tr_projection,
                                                   projection_info=self.rigid_projection_info,
                                                   frequency_zero_cutoff_cm=self.neg_freq_thresh_cm)
                n_imag = int(np.sum(freqs_est < -abs(self.neg_freq_thresh_cm)))
                click.echo(f"[tsopt] n≈{n_imag}  (approx imag: {[float(x) for x in freqs_est if x < -abs(self.neg_freq_thresh_cm)]})")
                if n_imag <= 1:
                    break

                # (b) Get approximate modes for flattening (embedded, mass-weighted)
                freqs_cm_approx, modes_embedded = _modes_from_Hact_embedded(
                    H_act, self.geom.atomic_numbers, self.geom.cart_coords.reshape(-1, 3), active_idx, self.device,
                    tr_projection=self.tr_projection,
                    projection_info=self.rigid_projection_info,
                    frequency_zero_cutoff_cm=self.neg_freq_thresh_cm,
                )

                # (c) Do flatten step using the approximate modes
                if self.flatten_loop_bofill:
                    x_before_flat = self.geom.cart_coords.copy().reshape(-1)
                did_flatten = self._flatten_once_with_modes(freqs_cm_approx, modes_embedded)
                # Free GPU tensors from mode computation immediately after use
                del freqs_cm_approx, modes_embedded
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
                if not did_flatten:
                    self.flatten_skip_reason = "no eligible extra imaginary modes"
                    break

                # (d) Bofill update using UMA gradients across the flatten displacement
                if self.flatten_loop_bofill:
                    x_after_flat = self.geom.cart_coords.copy().reshape(-1)
                    g_after_flat = _calc_gradient(
                        self.geom, self.calc_kwargs
                    ).reshape(-1)
                    delta_flat_full = x_after_flat - x_before_flat
                    delta_flat_act = delta_flat_full[mask_dof]
                    g_old_act = g_prev[mask_dof]
                    g_new_act = g_after_flat[mask_dof]
                    H_act = _bofill_update_active(
                        H_act, delta_flat_act, g_new_act, g_old_act
                    )

                # (e) Refresh dimer direction from updated active Hessian
                mode_xyz = _mode_direction_by_root_from_Hact(
                    H_act, self.geom.cart_coords.reshape(-1, 3), self.geom.atomic_numbers,
                    self.masses_au_t, active_idx, self.device, root=self.root,
                    tr_projection=self.tr_projection,
                    projection_info=self.rigid_projection_info,
                )
                np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")
                del mode_xyz

                # (f) Re-optimize with Dimer (consumes global cycle budget)
                # Clear VRAM before dimer loop to ensure space for Hessian recomputation
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
                _, zero_step_flat, conv_flat = self._dimer_loop(self.thresh)
                self.is_converged = conv_flat

                # A stall inside the flatten loop stops the remaining iterations
                # : do not keep retrying a stalled optimization.
                if self.is_stalled:
                    self.flatten_skip_reason = (
                        "optimization stalled during flattening"
                    )
                    break

                if (
                    self.max_total_cycles is not None
                    and (self.max_total_cycles - self._cycles_spent) <= 0
                ):
                    self.flatten_skip_reason = (
                        "max-cycles budget exhausted during flattening"
                    )
                    break

                if self.flatten_loop_bofill:
                    # (g) Update across the optimization displacement.
                    x_after_opt = self.geom.cart_coords.copy().reshape(-1)
                    g_after_opt = _calc_gradient(
                        self.geom, self.calc_kwargs
                    ).reshape(-1)
                    delta_opt_full = x_after_opt - x_after_flat
                    delta_opt_act = delta_opt_full[mask_dof]
                    g_old_act2 = g_after_flat[mask_dof]
                    g_new_act2 = g_after_opt[mask_dof]
                    H_act = _bofill_update_active(
                        H_act, delta_opt_act, g_new_act2, g_old_act2
                    )
                    g_prev = g_after_opt
                else:
                    # The default path makes the next mode decision from the
                    # Hessian at the current coordinates, not from a stale
                    # pre-flatten approximation.
                    H_next = self._calc_full_hessian_cached(
                        self.calc_kwargs_full,
                        allow_reuse=zero_step_flat,
                    )
                    H_final_reuse_cpu = H_next.detach().cpu().clone()
                    H_final_reuse_coords = self.geom.cart_coords.copy()
                    if H_next.size(0) == 3 * N:
                        H_next_active = _extract_active_block(H_next, mask_dof)
                        del H_next
                    else:
                        active_idx, mask_dof = self._resolve_hessian_active_subspace(
                            H_next, N
                        )
                        H_next_active = H_next
                    old_H_act = H_act
                    H_act = H_next_active
                    del old_H_act
                    _clear_cuda_cache()
        elif self.flatten_max_iter > 0:
            self.flatten_skip_reason = (
                "max-cycles budget exhausted before flattening"
            )
            click.echo("[tsopt] Reached --max-cycles budget; skipping flatten loop.")

        # (5) Final outputs
        final_xyz = self.out_dir / "final_geometry.xyz"
        atoms_final = Atoms(self.geom.atoms, positions=(self.geom.coords3d * BOHR2ANG), pbc=False)
        write(final_xyz, atoms_final)

        if not self.is_converged:
            self.saddle_order_verified = False
            self.n_imaginary_modes = None
            self.imaginary_frequencies_cm = []
            self.hessian_status = "skipped"
            self.hessian_error = None
            return

        # Final Hessian → imaginary mode trajectory
        if self.skip_final_freq and not self.is_stalled:
            click.echo(
                "[tsopt] WARNING: TS saddle-point order is not verified "
                "(--skip-final-freq).",
                err=True,
            )
            self.saddle_order_verified = False
            self.n_imaginary_modes = None
            self.imaginary_frequencies_cm = []
            self.hessian_status = "skipped"
            self.hessian_error = None
            return

        try:
            reuse_final_hessian = (
                H_final_reuse_cpu is not None
                and H_final_reuse_coords is not None
                and np.array_equal(self.geom.cart_coords, H_final_reuse_coords)
            )
            if reuse_final_hessian:
                click.echo("[tsopt] Reusing flatten-start Hessian for final frequency analysis (geometry unchanged).")
                H_t = H_final_reuse_cpu.to(self.device)
            else:
                H_t = _calc_full_hessian_torch(self.geom, self.calc_kwargs_full, self.device)
            raw_hessian_shape = tuple(H_t.shape)
            H_analysis, active_idx_final, computed_atoms, storage = (
                _reconcile_hessian_analysis_basis(
                    H_t,
                    self.geom,
                    self.analysis_active_atoms,
                )
            )
            del H_t
            freqs_cm, modes = _modes_from_Hact_embedded(
                H_analysis,
                self.geom.atomic_numbers,
                self.geom.cart_coords.reshape(-1, 3),
                active_idx_final,
                self.device,
                tr_projection=self.tr_projection,
                projection_info=self.rigid_projection_info,
                frequency_zero_cutoff_cm=self.neg_freq_thresh_cm,
            )

            self.rigid_projection_info.update({
                "hessian_space": "full" if len(active_idx_final) == N else "active",
                "analysis_hessian_shape": list(H_analysis.shape),
                "raw_hessian_shape": list(raw_hessian_shape),
                "computed_atom_count": len(computed_atoms),
                "analysis_atom_count": len(active_idx_final),
                "storage": storage,
                "source": "tsopt_exact",
            })
            projection_block = pretty_block(
                "rigid_projection", self.rigid_projection_info
            )
            if projection_block:
                click.echo(projection_block)

            del H_analysis
            del H_final_reuse_cpu, H_final_reuse_coords
            _ref_pdb_light = (
                self.source_path
                if self.source_path is not None and self.source_path.suffix.lower() == ".pdb"
                else None
            )
            n_written = _write_all_imag_modes(
                self.geom,
                freqs_cm,
                modes,
                self.neg_freq_thresh_cm,
                self.vib_dir,
                ref_pdb=_ref_pdb_light,
            )
            _finalize_dimer_saddle_status(
                self, freqs_cm, self.neg_freq_thresh_cm
            )
            _n_imag = int(self.n_imaginary_modes or 0)
            mode_message, mode_message_is_diagnostic = _dimer_mode_export_message(
                n_written,
                _n_imag,
                self.neg_freq_thresh_cm,
                float(freqs_cm.min()),
            )
            click.echo(mode_message, err=mode_message_is_diagnostic)
            del modes, freqs_cm

            _clear_cuda_cache()
            self.hessian_status = "completed"
            self.hessian_error = None
        except Exception as exc:
            self.hessian_status = "failed"
            self.hessian_error = f"{type(exc).__name__}: {exc}"
            self.saddle_order_verified = False
            self.n_imaginary_modes = None
            self.imaginary_frequencies_cm = []
            click.echo("[tsopt] ERROR: Terminal PHVA failed.", err=True)
            emit(
                f"[tsopt] Terminal PHVA error: {self.hessian_error}",
                detail=True,
            )
            _clear_cuda_cache()
            click.echo(f"[tsopt] Saved final geometry → {final_xyz}")
            return
        click.echo(f"[tsopt] Saved final geometry → {final_xyz}")
        emit(f"[tsopt] Mode files → {self.vib_dir}", detail=True)




# Macro/micro alternation follows the Gaussian 16 ONIOM(QM:MM) algorithm
# described in the function docstring.
# CHEMISTRY-RULE:3 Macro/micro alternation (Gaussian 16 microiteration、ML+linkparent co-macro)。
def _run_microiter_tsopt(
    geometry,
    calc_cfg: Dict[str, Any],
    rsirfo_cfg: Dict[str, Any],
    lbfgs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    microiter_cfg: Dict[str, Any],
    out_dir_path: Path,
    *,
    dump: bool = False,
    thresh: Optional[str] = None,
    mode: str = "rsirfo",
    reference_mode: Optional[np.ndarray] = None,
    flatten_enabled: bool = False,
) -> Dict[str, Any]:
    """Run macro/micro alternating TS optimization (Gaussian 16-style microiteration).

    Macro step: 1 Hessian-based TS step (RS-I-RFO / RS-P-RFO / TRIM) moving ML atoms + link-atom MM parents (full ONIOM force).
    Micro step: LBFGS relaxing remaining MM atoms with MM-only forces until convergence.
    Link-atom MM parents are included in the macro step to maintain consistency
    of the link atom position across macro/micro boundaries.

    Convergence semantics (CHEMISTRY-RULE:3 macro/micro alternation):
      - Macro: full ONIOM gradient (= ML force + MM correction at link), converged
        when max ‖F_ML+linkparent‖ < `thresh` (default ``baker``).
      - Micro: MM-only on remaining movable_mm + hess_mm atoms, converged when
        max ‖F_MM‖ < ``micro_thresh`` (defaults to the macro preset).
      - Outer convergence: the first macro step that satisfies the
        macro thresholds terminates the loop (no consecutive-step
        debounce — the macro/micro cadence already filters spurious
        single-cycle convergence).

    Why ML+linkparent in macro (not ML-only):
      Link atom (capped H) position is a function of (ML host atom, MM parent atom).
      If MM parent moves in micro without ML host moving, link atom shifts incorrectly.
      Treating MM parent as "co-macro" (= moved together with ML in macro step)
      keeps link atom geometry consistent across macro/micro alternation cycles.
      Reference: Vreven & Morokuma 2003 (ONIOM microiteration scheme).

    GPU memory pattern (CHEMISTRY-RULE:1 + bofill_update CPU fallback):
      RS-I-RFO Bofill Hessian update is forced to CPU (= long-loop GPU peak
      avoidance); this function's macro step runs ML force on GPU + Hessian
      update on CPU + MM-only micro on CPU. No GPU residency growth across
      macro cycles.
    """
    # Resolve the immutable partition from a single accepted core
    # (loud on failure, never a swallowed empty set), consuming the exact
    # geometry freeze mask as the immutable original freeze. The user's freeze is
    # preserved in BOTH phase masks and restored exactly in ``finally``.
    entry_calculator = getattr(geometry, "calculator", None)
    temp_calc = mlmm(**dict(calc_cfg))
    calc_core = temp_calc.core if hasattr(temp_calc, "core") else temp_calc
    partition = resolve_partition_from_core(
        calc_core, len(geometry.atoms), geometry.freeze_atoms
    )
    del temp_calc, calc_core
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    if not partition.has_macro_active:
        # A TS microiteration without an ML macro region is uncomputable; fail
        # loudly (never return an unchanged geometry as an optimized result).
        raise PartitionError(
            "Microiteration requires at least one ML macro-active atom; the "
            "resolved partition has none."
        )

    ml_indices = list(partition.ml_atoms)
    movable_mm = list(partition.movable_mm_atoms)
    link_mm_parents = set(partition.link_parent_atoms)
    original_freeze = list(partition.original_freeze)
    frozen_mm = list(original_freeze)

    n_atoms = partition.n_atoms

    macro_freeze = list(partition.macro_freeze_atoms)
    micro_freeze = list(partition.micro_freeze_atoms)

    max_cycles = optional_positive_int(opt_cfg.get("max_cycles"), "opt.max_cycles")
    macro_thresh = thresh if thresh is not None else rsirfo_cfg.get("thresh", "baker")
    micro_thresh = microiter_cfg.get("micro_thresh") or macro_thresh
    micro_max_cycles = optional_positive_int(
        microiter_cfg.get("micro_max_cycles"), "microiter.micro_max_cycles"
    )

    click.echo(
        f"[microiter] ML atoms: {len(ml_indices)}, "
        f"Link MM parents: {len(link_mm_parents)}, "
        f"Movable MM atoms: {len(movable_mm)}, "
        f"Frozen MM atoms: {len(frozen_mm)}"
    )
    click.echo(f"[microiter] Macro thresh: {macro_thresh}, Micro thresh: {micro_thresh}")

    # Create ONIOM calculator (shared core for MM-only calc)
    macro_calc_cfg = dict(calc_cfg)
    macro_calc_cfg["freeze_atoms"] = macro_freeze
    macro_calc_cfg["hess_mm_atoms"] = sorted(link_mm_parents)  # ML + link MM parents in Hessian
    macro_calc = mlmm(**macro_calc_cfg)
    mm_calc = mlmm_mm_only(macro_calc.core, freeze_atoms=micro_freeze)
    # Full-constraint calculator restored on the normal exit path . Built
    # lazily just before returning so it does not add a second live core to the
    # macro/micro loop's VRAM footprint (mm_calc reuses macro_calc.core).
    base_calc = None

    # Ordered record of every micro (MM) relaxation outcome.
    micro_attempts: List[OptimizerOutcome] = []
    micro_cycles_total = 0

    def _relax_micro() -> Tuple[Any, int]:
        """Run one MM-only micro relaxation on a cart twin; copy coords back.

        Returns the LBFGS optimizer (for its explicit convergence bit) and the
        number of executed micro cycles.  A normal Python return is not evidence
        of convergence.
        """

        macro_coord_type = getattr(geometry, "coord_type", "cart")
        if macro_coord_type != "cart":
            from pysisyphus.Geometry import Geometry as _Geom
            micro_geom = _Geom(
                atoms=tuple(geometry.atoms),
                coords=geometry.coords3d.copy().flatten(),
                coord_type="cart",
                freeze_atoms=micro_freeze,
            )
        else:
            geometry.freeze_atoms = micro_freeze
            micro_geom = geometry
        micro_geom.set_calculator(mm_calc)

        micro_lbfgs_args = dict(lbfgs_cfg)
        micro_lbfgs_args["max_cycles"] = micro_max_cycles
        micro_lbfgs_args["thresh"] = micro_thresh
        micro_lbfgs_args["out_dir"] = str(out_dir_path)
        micro_lbfgs_args["dump"] = dump
        # The MM equilibration never uses the plateau stop. A flat energy while
        # its forces are still above threshold is a stalled optimizer, not MM
        # equilibrium, and stopping there ends the whole macro/micro alternation
        # with the environment unrelaxed. `micro_max_cycles` is the real bound.
        micro_lbfgs_args["energy_plateau"] = False

        _micro_opt = LBFGS(micro_geom, **micro_lbfgs_args)
        with contextlib.redirect_stdout(io.StringIO()):
            _micro_opt.run()
        _micro_steps = max(int(_micro_opt.cur_cycle) + 1, 1)
        if macro_coord_type != "cart":
            geometry.coords3d = micro_geom.coords3d.flatten()
            micro_geom.set_calculator(None)
            del micro_geom
        return _micro_opt, _micro_steps

    try:
        optim_all_path = out_dir_path / "optimization_all_trj.xyz"
        macro_trj_path = out_dir_path / "optimization_trj.xyz"
        total_macro_steps = 0
        run_macro = True
        latest_micro_stalled = False
        latest_micro_stop_reason = ""

        # Establish the MM equilibrium before resolving Hessian identity or
        # constructing the persistent macro optimizer.  Otherwise its first
        # RFO/TRIM step would use curvature from pre-relaxation coordinates.
        if partition.has_micro_active:
            _init_micro_opt, _init_micro_steps = _relax_micro()
            micro_cycles_total += _init_micro_steps
            _init_micro_out = OptimizerOutcome.from_optimizer(
                _init_micro_opt, max_cycles=micro_max_cycles
            )
            # keeps the macro step off unless the micro relaxation settled.
            # A plateau whose forces already meet the configured thresholds IS
            # MM equilibrium, though: refusing it costs the whole TS search
            # (zero macro steps) over step criteria the macro does not need.
            _init_micro_equilibrium = (
                _init_micro_out.converged is not True
                and _init_micro_out.stalled
                and micro_reached_force_equilibrium(_init_micro_opt)
            )
            if _init_micro_equilibrium:
                _init_micro_out = _init_micro_out.accept_force_equilibrium()
                emit(
                    "[microiter] Initial MM equilibration plateaued with its "
                    "force criteria met; accepting it as MM equilibrium.",
                    narrative=True,
                )
            micro_attempts.append(_init_micro_out)
            latest_micro_stalled = _init_micro_out.stalled
            latest_micro_stop_reason = _init_micro_out.stop_reason or ""
            if _init_micro_out.converged is not True and not _init_micro_equilibrium:
                run_macro = False
                if dump:
                    _append_xyz_trajectory(
                        optim_all_path,
                        out_dir_path / "optimization_trj.xyz",
                        reset=True,
                    )
                emit(
                    "[microiter] Initial MM equilibration did not converge "
                    f"(status={_init_micro_out.status}); no macro step is taken.",
                    narrative=True,
                )
            del _init_micro_opt
            _clear_cuda_cache()
            geometry.freeze_atoms = macro_freeze
            geometry.set_calculator(macro_calc)

        if not run_macro:
            macro_outcome = OptimizerOutcome.not_executed(
                max_cycles=max_cycles,
                reason=_init_micro_out.stop_reason or "initial_micro_not_converged",
            )
            aggregate = build_aggregate(
                macro_outcome, micro_attempts, max_cycles=max_cycles
            )
            micro_outcome = MicroiterationOutcome(
                aggregate=aggregate,
                macro=macro_outcome,
                micro_attempts=tuple(micro_attempts),
                macro_cycles=0,
                micro_cycles=micro_cycles_total,
                partition=partition,
            )
            return {
                "converged": False,
                "cycles": 0,
                "stop_requested": False,
                "stop_reason": aggregate.stop_reason,
                "is_stalled": bool(aggregate.stalled),
                "safeguards": {},
                "optimizer": None,
                "micro_cycles": micro_cycles_total,
                "outcome": micro_outcome,
            }

        # Seed initial Hessian for RS-I-RFO (with macro freeze)
        # Try TS Hessian cache first; fall back to full Hessian calculation.
        from mlmm.io.hessian_cache import (
            load_matching as _hess_load_matching,
            identity_from_context as _hess_identity,
            reconcile_active_hessian as _hess_reconcile_active,
        )
        hess_device = _torch_device(calc_cfg.get("ml_device", "auto"))

        # reuse a cached TS Hessian only on a full evaluation-identity match.
        cached_ts = _hess_load_matching(
            "ts",
            _hess_identity(geometry, calc_cfg, role="ts"),
            atol=1.1e-3,
        )
        macro_free_atoms = sorted(
            set(range(geometry.cart_coords.size // 3)) - set(macro_freeze)
        )
        macro_free_dofs = [
            3 * atom + axis
            for atom in macro_free_atoms
            for axis in range(3)
        ]
        _cache_used = False
        if cached_ts is not None:
            h_init = _hess_reconcile_active(
                cached_ts,
                macro_free_dofs,
                full_n_dof=geometry.cart_coords.size,
            )
            if h_init is not None:
                click.echo(
                    "[microiter] Reusing cached TS Hessian for the macro TS step."
                )
                geometry.freeze_atoms = macro_freeze
                geometry.set_calculator(macro_calc)
                geometry.within_partial_hessian = {
                    "active_n_dof": len(macro_free_dofs),
                    "full_n_dof": geometry.cart_coords.size,
                    "active_dofs": macro_free_dofs,
                    "active_atoms": macro_free_atoms,
                }
                geometry.cart_hessian = h_init
                click.echo(
                    "[microiter] Initial Hessian seeded from cache "
                    f"(shape={h_init.shape[0]}x{h_init.shape[1]})."
                )
                _cache_used = True
                del h_init
            else:
                click.echo(
                    "[microiter] Cached TS Hessian basis does not cover the "
                    "ordered macro DOFs. Falling back to a fresh Hessian."
                )
        if not _cache_used:
            click.echo("[microiter] Seeding initial Hessian for the macro TS step.")

            geometry.freeze_atoms = macro_freeze
            geometry.set_calculator(macro_calc)

            h_init = _calc_full_hessian_torch(geometry, macro_calc_cfg, hess_device)
            geometry.cart_hessian = h_init
            click.echo(f"[microiter] Initial Hessian seeded (shape={h_init.shape[0]}x{h_init.shape[1]}).")
            del h_init

        # Create the persistent macro TS optimizer once (LayerOpt pattern).
        # This preserves the BFGS Hessian update chain across macro iterations.
        # The macro optimizer is the resolved --opt-mode (RS-I-RFO / RS-P-RFO / TRIM);
        # all three are TSHessianOptimizer subclasses sharing the optimize/prepare_opt
        # + Bofill-update contract the macro loop drives.
        # NOTE: geometry already has macro_calc set (line above); do NOT call
        # set_calculator again as it clears the pre-computed cart_hessian.
        geometry.freeze_atoms = macro_freeze

        rsirfo_args = _build_rsirfo_kwargs(
            rsirfo_cfg,
            max_cycles=max_cycles,
            out_dir=out_dir_path,
            macro_thresh=macro_thresh,
            mode=mode,
            opt_cfg=opt_cfg,
            reference_mode=reference_mode,
            flatten_enabled=flatten_enabled,
        )

        macro_optimizer = TSOPT_CLASS_MAP[mode](geometry, **rsirfo_args)
        macro_optimizer.prepare_opt()  # initialize Hessian from geometry.cart_hessian

        # Microiteration progress table (pysisyphus-style with micro_steps column)
        micro_header = "cycle Δ(energy) max(|force|) rms(force) max(|step|) rms(step) micro_steps s/cycle".split()
        micro_col_fmts = "int float float float float float int float_short".split()
        micro_table = TablePrinter(micro_header, micro_col_fmts, width=12)
        micro_table.print_header()
        printed_macro_rows = 0

        macro_converged = False
        macro_iterable = (
            count()
            if run_macro and max_cycles is None
            else range(max_cycles if run_macro else 0)
        )
        for macro_iter in macro_iterable:
            # ---- Macro step: 1 RS-I-RFO step with ONIOM forces, MM frozen ----
            geometry.freeze_atoms = macro_freeze
            geometry.set_calculator(macro_calc)

            # Manually feed state to the persistent optimizer (cf. LayerOpt lines 358-364)
            macro_optimizer.coords.append(geometry.coords.copy())
            macro_optimizer.cart_coords.append(geometry.cart_coords.copy())
            macro_optimizer.cur_cycle = macro_iter

            t_start = time.time()
            step = macro_optimizer.optimize()  # housekeeping() triggers BFGS update
            macro_optimizer.steps.append(step)

            # Convergence check
            macro_converged, conv_info = macro_optimizer.check_convergence()
            total_macro_steps += 1

            if macro_optimizer.stop_requested:
                # A real macro energy-plateau stall  is a distinct outcome
                # from a generic stop; it is never convergence.
                if macro_optimizer.is_stalled:
                    click.echo(
                        "[microiter] Stalled (energy plateau; not converged): "
                        f"{macro_optimizer.stop_reason}",
                        err=True,
                    )
                else:
                    click.echo(
                        "[microiter] Stopped without convergence: "
                        f"{macro_optimizer.stop_reason}",
                        err=True,
                    )
                break

            if dump:
                with open(macro_trj_path, "a") as f:
                    f.write(geometry.as_xyz() + "\n")
                _append_xyz_trajectory(optim_all_path, macro_trj_path)

            if macro_converged:
                # Print final converged row (no micro steps)
                energy_diff = macro_optimizer.energies[-1] - macro_optimizer.energies[-2] if len(macro_optimizer.energies) >= 2 else float("nan")
                marks = [False, *conv_info.get_convergence()[:-1], False, False]
                if not np.all(np.isfinite(np.asarray(energy_diff))):
                    marks[1] = False
                cycle_time = time.time() - t_start
                if printed_macro_rows and printed_macro_rows % 10 == 0:
                    micro_table.print_sep()
                micro_table.print_row(
                    (macro_iter, energy_diff, macro_optimizer.max_forces[-1], macro_optimizer.rms_forces[-1],
                     macro_optimizer.max_steps[-1], macro_optimizer.rms_steps[-1], 0, cycle_time),
                    marks=marks,
                )
                printed_macro_rows += 1
                print()  # blank line closes the table (print() shares the table's stdout path)
                emit("[microiter] Converged!", detail=True)
                break

            # Apply step to geometry
            new_coords = geometry.coords.copy() + step
            geometry.coords = new_coords
            # Record actual step (may differ due to coordinate back-transformation)
            macro_optimizer.steps[-1] = geometry.coords - macro_optimizer.coords[-1]

            # ---- Micro step: LBFGS with MM-only forces, ML frozen ----
            # ``_relax_micro`` runs the MM relaxation on a cart-only twin and copies
            # the converged positions back via the coords3d setter (preserving the
            # macro RFO's persistent Bofill-updated Hessian); the macro chemistry
            # stays in DLC.
            if partition.has_micro_active:
                micro_opt, micro_steps = _relax_micro()
                micro_cycles_total += micro_steps
                _micro_out = OptimizerOutcome.from_optimizer(micro_opt, max_cycles=micro_max_cycles)
                # Evaluate before `del`: the predicate needs the live optimizer.
                _micro_equilibrium = (
                    _micro_out.converged is not True
                    and _micro_out.stalled
                    and micro_reached_force_equilibrium(micro_opt)
                )
                if _micro_equilibrium:
                    _micro_out = _micro_out.accept_force_equilibrium()
                    emit(
                        "[microiter] MM relaxation plateaued with its force "
                        "criteria met; accepting it as MM equilibrium.",
                        narrative=True,
                    )
                del micro_opt
                _clear_cuda_cache()
            else:
                # A validated partition with macro-active atoms but
                # ZERO micro-active MM atoms (e.g. the entire MM region is
                # user-frozen) has no movable MM coordinate to relax. Append the
                # zero-cycle vacuous micro success instead of building an LBFGS with
                # every atom frozen (mirrors the initial-equilibration guard above).
                _micro_out = OptimizerOutcome.vacuous_success()
                micro_steps = 0
                _micro_equilibrium = False
            micro_attempts.append(_micro_out)
            # remember whether THIS micro (MM) relaxation stalled on an
            # energy plateau. Coordinate copying alone is not evidence of
            # convergence, so a stalled/non-converged latest micro relaxation must
            # not later read as clean macro convergence.
            latest_micro_stalled = _micro_out.stalled
            latest_micro_stop_reason = _micro_out.stop_reason or ""

            if dump:
                _append_xyz_trajectory(optim_all_path, out_dir_path / "optimization_trj.xyz")

            # a required micro relaxation that did not explicitly converge stops
            # the macro/micro alternation; the aggregate cannot be converged and no
            # further macro step is taken (a normal LBFGS return on max-cycle
            # exhaustion is not convergence).
            if _micro_out.converged is not True and not _micro_equilibrium:
                emit(
                    "[microiter] Latest MM relaxation did not converge "
                    f"({describe_micro_stop(_micro_out, micro_opt)}); "
                    "stopping the macro/micro loop.",
                    narrative=True,
                )
                print()
                break

            # Print progress row with micro_steps
            cycle_time = time.time() - t_start
            energy_diff = macro_optimizer.energies[-1] - macro_optimizer.energies[-2] if len(macro_optimizer.energies) >= 2 else float("nan")
            marks = [False, *conv_info.get_convergence()[:-1], False, False]
            if not np.all(np.isfinite(np.asarray(energy_diff))):
                marks[1] = False
            macro_print_every = getattr(
                macro_optimizer, "print_every", opt_cfg.get("print_every", 1)
            )
            if macro_progress_due(macro_iter, macro_print_every):
                if printed_macro_rows and printed_macro_rows % 10 == 0:
                    micro_table.print_sep()
                micro_table.print_row(
                    (macro_iter, energy_diff, macro_optimizer.max_forces[-1], macro_optimizer.rms_forces[-1],
                     macro_optimizer.max_steps[-1], macro_optimizer.rms_steps[-1], micro_steps, cycle_time),
                    marks=marks,
                )
                printed_macro_rows += 1

        else:
            if run_macro:
                print()  # blank line closes the table (print() shares the table's stdout path)
                emit(f"[microiter] Reached max macro iterations ({max_cycles}).", detail=True)

        # a stalled latest micro (MM) relaxation must not masquerade as
        # clean macro convergence, and it must not be lost as a reasonless
        # not_converged when the macro merely ran out of cycles: surface it as a
        # stall (with its reason) in either case so the returned state is explicit.
        finalize_microiter_macro_convergence(
            macro_optimizer,
            macro_converged=macro_converged,
            latest_micro_stalled=latest_micro_stalled,
            latest_micro_stop_reason=latest_micro_stop_reason,
        )
        # Fold the macro state after stall demotion and the ordered micro
        # attempts into one fail-closed aggregate outcome.
        macro_outcome = OptimizerOutcome.from_optimizer(macro_optimizer, max_cycles=max_cycles)
        aggregate = build_aggregate(macro_outcome, micro_attempts, max_cycles=max_cycles)
        micro_outcome = MicroiterationOutcome(
            aggregate=aggregate,
            macro=macro_outcome,
            micro_attempts=tuple(micro_attempts),
            macro_cycles=total_macro_steps,
            micro_cycles=micro_cycles_total,
            partition=partition,
        )
        outcome = {
            "converged": bool(aggregate.converged is True),
            "cycles": total_macro_steps,
            "stop_requested": bool(macro_optimizer.stop_requested),
            "stop_reason": aggregate.stop_reason or (macro_optimizer.stop_reason or None),
            "is_stalled": bool(aggregate.stalled),
            "safeguards": _optimizer_safeguard_payload(macro_optimizer),
            "optimizer": macro_optimizer,
            "micro_cycles": micro_cycles_total,
            "outcome": micro_outcome,
        }
        _clear_cuda_cache()

        emit(f"[microiter] Total macro steps: {total_macro_steps}", detail=True)

        # Restore the full-constraint calculator for downstream use.
        base_calc = mlmm(**calc_cfg)
        return outcome
    finally:
        # restore the EXACT original freeze mask on every exit path
        # (success, micro non-convergence, macro stall, and raised exception);
        # restore the full-constraint calculator whenever one was built.
        geometry.freeze_atoms = list(original_freeze)
        if base_calc is not None:
            geometry.set_calculator(base_calc)
        else:
            geometry.set_calculator(entry_calculator)
        _clear_cuda_cache()



# Configuration defaults (imported from defaults.py)
GEOM_KW: Dict[str, Any] = deepcopy(GEOM_KW_DEFAULT)
CALC_KW: Dict[str, Any] = deepcopy(MLMM_CALC_KW)

# HessianDimer defaults - combine imported DIMER_KW and HESSIAN_DIMER_KW
hessian_dimer_KW = {
    **HESSIAN_DIMER_KW,
    "dimer": {**DIMER_KW},
    # The dimer rotates the physical gradient into an effective TS force, so
    # physical-energy uphill rejection is not meaningful for its inner LBFGS.
    "lbfgs": {
        **{k: v for k, v in LBFGS_KW.items() if k != "max_cycles"},
        "reject_uphill": False,
    },
}


def _load_reference_mode(path: Path, expected_size: int) -> np.ndarray:
    """Load the primary normalized mode through the shared cache reader.

    Direct callers therefore use the same validation as the multi-candidate
    TS workflow instead of maintaining a second file parser.
    """

    candidates, _labels, _metadata = read_reference_mode_candidates(
        Path(path), int(expected_size)
    )
    return candidates[0]

def _validate_reference_mode_optimizer(
    mode: str,
    reference_mode_path: Optional[Path],
) -> None:
    """Reject path-mode guidance for optimizers that do not consume it."""
    if reference_mode_path is not None and mode == "dimer":
        raise click.BadParameter(
            "--ref-mode requires a Hessian TS optimizer; use "
            "--opt-mode hess, rsirfo, rsprfo, or trim.",
            param_hint="--ref-mode",
        )


def _heavy_mode_label(mode: str) -> str:
    """Return the public label for one Hessian transition-state optimizer."""

    return {
        "rsirfo": "RS-I-RFO",
        "rsprfo": "RS-P-RFO",
        "trim": "TRIM",
    }[mode]


def _post_analysis_hessian_config(
    calc_cfg: Dict[str, Any],
    *,
    partial: bool,
) -> Dict[str, Any]:
    """Resolve the Hessian storage policy used by final analysis/flattening."""

    resolved = dict(calc_cfg)
    resolved["out_hess_torch"] = True
    resolved["return_partial_hessian"] = bool(partial)
    return resolved


class _TSOPTOutputCollisionError(click.UsageError):
    """A TSOPT output/input collision."""


def _tsopt_owned_output_paths(path: Path) -> List[Path]:
    """Return command-owned TS artifacts for one output directory."""

    resolved = Path(path).resolve()
    vib_dir = resolved / "vib"
    if vib_dir.is_symlink() or (vib_dir.exists() and not vib_dir.is_dir()):
        raise _TSOPTOutputCollisionError(
            f"TSOPT vib output must be a real directory: {vib_dir}."
        )
    return [
        *(
            resolved / name
            for name in (
                "final_geometry.xyz",
                "final_geometry.pdb",
                "final_geometry.gjf",
                "optimization_all_trj.xyz",
                "optimization_all.pdb",
                "optimization_trj.xyz",
                "optimization.pdb",
                "result.json",
                "summary.json",
            )
        ),
        *(
            candidate
            for pattern in ("imag_*.pdb", "imag_*_trj.xyz")
            for candidate in vib_dir.glob(pattern)
            if candidate.is_file()
        ),
    ]


def _reject_tsopt_output_collisions(
    path: Path,
    *,
    protected_inputs: Sequence[Optional[Path]] = (),
) -> None:
    """Reject inputs that occupy deterministic TSOPT output paths."""

    resolved = Path(path).resolve()
    reserved = {
        *(candidate.resolve() for candidate in _tsopt_owned_output_paths(resolved)),
        (resolved / "model_from_bfactor.pdb").resolve(),
    }
    for protected in protected_inputs:
        if protected is not None and Path(protected).resolve() in reserved:
            raise _TSOPTOutputCollisionError(
                f"Input {protected} collides with a reserved TSOPT output path "
                f"under {resolved}."
            )


def _prepare_tsopt_output_dir(
    path: Path,
    *,
    protected_inputs: Sequence[Optional[Path]] = (),
) -> Path:
    """Invalidate command-owned TS artifacts before a real generation."""

    resolved = Path(path).resolve()
    resolved.mkdir(parents=True, exist_ok=True)
    _reject_tsopt_output_collisions(
        resolved,
        protected_inputs=protected_inputs,
    )
    owned = _tsopt_owned_output_paths(resolved)
    for candidate in owned:
        candidate.unlink(missing_ok=True)
    return resolved


@click.command(
    help="TS optimization: grad (Dimer) or hess (RS-I-RFO) for the ML/MM calculator.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Starting geometry (PDB/mmCIF or XYZ). XYZ provides higher coordinate precision. "
         "If XYZ, use --ref-pdb to specify PDB topology for atom ordering and output conversion.",
)
@click.option(
    "--ref-mode",
    "reference_mode_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help=(
        "Advanced/internal Cartesian reference direction(s) for Hessian TS "
        "root selection and overlap tracking. Accepts .npz path-mode caches, "
        ".npy arrays, or whitespace text containing one 3N vector or a 2-D "
        "candidate table. This guides mode identity; it does not replace the "
        "Hessian and is not supported by Dimer. The all workflow supplies it "
        "from the MEP; standalone tsopt users normally leave it unset."
    ),
)
@click.option(
    "--ref-pdb",
    "ref_pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    show_default=False,
    help="Reference PDB topology when input is XYZ. XYZ coordinates are used (higher precision) "
         "while PDB provides atom ordering and residue information for output conversion.",
)
@click.option(
    "--parm",
    "real_parm7",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Amber parm7 topology for the whole enzyme (MM region).",
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
    help="Total charge of the ML region. Required unless --ligand-charge is provided.",
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
    help="Spin multiplicity (2S+1) for the ML region.",
)
@click.option(
    "--freeze-atoms",
    "freeze_atoms_text",
    type=str,
    default=None,
    show_default=False,
    help="Comma-separated 1-based indices to freeze (e.g., '1,3,5').",
)
@click.option(
    "--radius-hessian",
    "--hess-cutoff",
    "hess_cutoff",
    type=float,
    default=None,
    show_default="all movable MM atoms",
    help="Distance cutoff (Å) from ML region for MM atoms to include in Hessian calculation. "
         "Applied to movable MM atoms. Unset includes every required movable MM atom; "
         "0.0 requests an ML-only Hessian and should be paired with "
         "--active-dof-mode ml-only for final frequency validation.",
)
@click.option(
    "--movable-cutoff",
    "movable_cutoff",
    type=float,
    default=None,
    show_default="use freeze_atoms",
    help="Distance cutoff (Å) from ML region for movable MM atoms. "
         "MM atoms beyond this are frozen. "
         "Providing --movable-cutoff disables --detect-layer.",
)
@click.option(
    "--hessian-calc-mode",
    type=click.Choice(["Analytical", "FiniteDifference"], case_sensitive=False),
    default=None, show_default="FiniteDifference",
    help="How the ML backend builds the Hessian (Analytical or FiniteDifference); "
         "overrides calc.hessian_calc_mode from YAML. "
         "Runtime and memory depend on the backend and system; compare both "
         "modes on a representative pilot.",
)
@click.option(
    "--max-cycles",
    type=click.IntRange(min=1),
    default=None,
    show_default="100000",
    help="Maximum total optimization cycles.",
)
@click.option(
    "--dump/--no-dump",
    default=False,
    show_default=True,
    help="Write concatenated trajectory 'optimization_all_trj.xyz'.",
)
@click.option("-o", "--out-dir", type=str, default=OUT_DIR_TSOPT, show_default=True, help="Output directory.")
@click.option(
    "--thresh",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default=None,
    show_default="baker",
    help="Convergence preset.",
)
@click.option(
    "--opt-mode",
    type=click.Choice(
        ["grad", "hess", "light", "heavy", "dimer", "rsirfo", "trim", "rsprfo"],
        case_sensitive=False,
    ),
    default="hess",
    show_default=True,
    help=(
        "grad/dimer/light → Hessian Guided Dimer; "
        "hess/rsirfo/heavy → RS-I-RFO; trim → TRIM (Helgaker); "
        "rsprfo → RS-P-RFO (Banerjee). "
        "All three Hessian TS optimizers (rsirfo/rsprfo/trim) are microiter-capable."
    ),
)
@click.option(
    "--microiter/--no-microiter",
    "microiter",
    default=True,
    show_default=True,
    help="Enable microiteration: alternate a 1-step macro TS move (RS-I-RFO / RS-P-RFO / TRIM) and MM relaxation (L-BFGS with MM-only forces). "
         "Effective in any Hessian --opt-mode (hess/rsirfo/rsprfo/trim); ignored in grad/dimer mode.",
)
@click.option(
    "--partial-hessian-flatten/--full-hessian-flatten",
    "partial_hessian_flatten",
    default=True,
    show_default=True,
    help="Use partial (active-block) Hessian for imaginary mode detection in flatten loop.",
)
@click.option(
    "--flatten/--no-flatten",
    "flatten",
    default=None,
    show_default="no-flatten",
    help="Enable/disable extra imaginary-mode flattening loop. "
         "--flatten uses the default flatten_max_iter (50); --no-flatten forces it to 0. "
         "When not provided, the loop is disabled unless YAML/config enables it.",
)
@click.option(
    "--ml-only-hessian-dimer/--no-ml-only-hessian-dimer",
    "ml_only_hessian_dimer",
    default=False,
    show_default=True,
    help="Use ML-region-only Hessian (no MM Hessian contribution) for dimer orientation "
         "in grad mode. Faster but less accurate for mode direction.",
)
@click.option(
    "--active-dof-mode",
    type=click.Choice(["all", "ml-only", "partial", "unfrozen"], case_sensitive=False),
    default="partial",
    show_default=True,
    help="Active DOF selection for final frequency analysis: "
         "all (all atoms), ml-only (ML only), partial (ML + MovableMM, default), "
         "unfrozen (all except frozen layer).",
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
    help="Validate options and print the execution plan without running TS optimization.",
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
    help="ML backend for the ONIOM high-level region.",
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
    help="Link-atom position mode: scaled (g-factor) or fixed (legacy 1.09/1.01 Å).",
)
@click.option(
    "--mm-backend",
    "mm_backend",
    type=click.Choice(["hessian_ff", "openmm"], case_sensitive=False),
    default=None,
    show_default="hessian_ff",
    help="MM backend. MM Hessians use finite differences by default; set calc.mm_fd: false for the hessian_ff analytical path.",
)
@click.option(
    "--cmap/--no-cmap",
    "use_cmap",
    default=None,
    show_default="cmap",
    help="Preserve CMAP terms in both real and model MM layers when present in parm7.",
)
@click.option(
    "--skip-final-freq/--no-skip-final-freq",
    "skip_final_freq",
    default=False,
    show_default=True,
    help=(
        "Skip terminal PHVA/frequency analysis and imaginary-mode flattening. "
        "Standalone tsopt retains the final structure with unverified saddle "
        "order; mlmm all stops before IRC because no imaginary direction can "
        "be validated."
    ),
)
@click.option(
    "--out-json/--no-out-json",
    "out_json",
    default=False,
    show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
@add_ml_layer_detection_options()
@add_precision_option()
@add_workers_options()
@add_backend_model_option()
@add_calc_file_option()
@add_deterministic_option()
@add_coord_type_option()
@add_print_every_option()
@add_allow_charge_mult_mismatch_option()
@click.pass_context
@click.option(
    "--stop-plateau/--no-stop-plateau",
    "stop_plateau",
    default=False,
    show_default=True,
    help=(
        "Stop when the energy stops changing while the convergence criteria are "
        "still unmet, and report the run as stalled. It never signals "
        "convergence; --max-cycles remains the real bound. The MM micro "
        "iterations are never stopped this way."
    ),
)
@click.option(
    "--stop-plateau-thresh",
    "stop_plateau_thresh",
    type=float,
    default=None,
    show_default="1e-4",
    help="Energy range (hartree) below which --stop-plateau treats the window as flat.",
)
@click.option(
    "--stop-plateau-window",
    "stop_plateau_window",
    type=int,
    default=None,
    show_default="50",
    help="Number of consecutive cycles --stop-plateau inspects.",
)
def cli(
    ctx: click.Context,
    stop_plateau: bool,
    stop_plateau_thresh: Optional[float],
    stop_plateau_window: Optional[int],
    input_path: Path,
    reference_mode_path: Optional[Path],
    ref_pdb: Optional[Path],
    real_parm7: Path,
    model_pdb: Optional[Path],
    model_indices_str: Optional[str],
    model_indices_one_based: bool,
    detect_layer: bool,
    charge: Optional[int],
    ligand_charge: Optional[str],
    spin: Optional[int],
    freeze_atoms_text: Optional[str],
    hess_cutoff: Optional[float],
    movable_cutoff: Optional[float],
    hessian_calc_mode: Optional[str],
    max_cycles: Optional[int],
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    opt_mode: str,
    microiter: bool,
    partial_hessian_flatten: bool,
    flatten: Optional[bool],
    ml_only_hessian_dimer: bool,
    active_dof_mode: str,
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    convert_files: bool,
    backend: Optional[str],
    embedcharge: bool,
    embedcharge_cutoff: Optional[float],
    link_atom_method: Optional[str],
    mm_backend: Optional[str],
    use_cmap: Optional[bool],
    skip_final_freq: bool,
    out_json: bool,
    precision: Optional[str],
    workers: Optional[int],
    workers_per_node: Optional[int],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: Optional[str],
    cli_coord_type: Optional[str],
    print_every: Optional[int],
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

    # Handle PDB/mmCIF directly, or XYZ with --ref-pdb for topology.
    suffix = input_path.suffix.lower()
    if suffix in {".pdb", ".cif", ".mmcif"}:
        prepared_input = prepare_input_structure(input_path)
    elif suffix == ".xyz":
        # XYZ input: require --ref-pdb for topology
        if ref_pdb is None:
            click.echo("ERROR: XYZ/TRJ input requires --ref-pdb to specify PDB topology.", err=True)
            sys.exit(1)
        prepared_input = prepare_input_structure(input_path)
        apply_ref_pdb_override(prepared_input, ref_pdb)
        click.echo(f"[input] Using XYZ coordinates from {input_path.name}, PDB topology from {ref_pdb.name}")
    else:
        click.echo(f"ERROR: Unsupported input format: {suffix}. Use .pdb/.cif/.mmcif or .xyz (with --ref-pdb).", err=True)
        sys.exit(1)

    geom_input_path = prepared_input.geom_path
    source_path = prepared_input.source_path
    reference_modes: List[np.ndarray] = []
    reference_mode_labels: List[str] = []
    reference_mode_metadata: Dict[str, Any] = {}
    if reference_mode_path is not None:
        try:
            n_atoms_reference = len(ase_read(str(geom_input_path), index=0))
            (
                reference_modes,
                reference_mode_labels,
                reference_mode_metadata,
            ) = read_reference_mode_candidates(
                reference_mode_path, 3 * n_atoms_reference
            )
        except (OSError, RuntimeError, ValueError) as exc:
            prepared_input.cleanup()
            raise click.BadParameter(str(exc), param_hint="--ref-mode") from exc
    reference_mode = reference_modes[0] if reference_modes else None
    charge, spin = resolve_charge_spin_or_raise(
        prepared_input, charge, spin,
        ligand_charge=ligand_charge, prefix="[tsopt]",
        model_pdb=model_pdb,
        model_indices_spec=model_indices_str,
        detect_layer=detect_layer,
        yaml_cfg=merged_yaml_cfg,
    )

    try:
        freeze_atoms_cli = _parse_freeze_atoms_opt(freeze_atoms_text)
    except click.BadParameter as e:
        click.echo(f"ERROR: {e}", err=True)
        prepared_input.cleanup()
        sys.exit(1)

    model_indices: Optional[List[int]] = None
    if model_indices_str:
        try:
            model_indices = parse_indices_string(model_indices_str, one_based=model_indices_one_based)
        except click.BadParameter as e:
            click.echo(f"ERROR: {e}", err=True)
            prepared_input.cleanup()
            sys.exit(1)

    time_start = time.perf_counter()

    # Resolve optimizer mode (default is now hess/RS-I-RFO)
    mode_resolved = normalize_choice(
        opt_mode,
        param="--opt-mode",
        alias_groups=TSOPT_MODE_ALIASES,
        allowed_hint="grad|hess|dimer|rsirfo|trim|rsprfo",
    )
    _validate_reference_mode_optimizer(mode_resolved, reference_mode_path)
    # trim/rsprfo are Hessian-based TS opts like rsirfo (use the same non-dimer code path).
    # Flatten outcome, published in result.json: today the only signal that a
    # requested flatten never ran is one stderr line, which no downstream
    # consumer reads.
    _flatten_skip_reason = None
    use_heavy = (mode_resolved in ("rsirfo", "trim", "rsprfo"))

    config_layer_cfg = load_yaml_dict(config_yaml)
    override_layer_cfg = load_yaml_dict(override_yaml)
    geom_cfg: Dict[str, Any] = deepcopy(GEOM_KW)
    calc_cfg: Dict[str, Any] = deepcopy(CALC_KW)
    opt_cfg: Dict[str, Any] = dict(OPT_BASE_KW)
    opt_cfg["out_dir"] = OUT_DIR_TSOPT
    lbfgs_cfg: Dict[str, Any] = dict(LBFGS_KW)
    simple_cfg: Dict[str, Any] = deepcopy(hessian_dimer_KW)
    # Keep the flatten loop off unless enabled by YAML/config or explicit --flatten.
    simple_cfg["flatten_max_iter"] = 0
    rsirfo_cfg: Dict[str, Any] = dict(RSIRFO_KW)
    frequency_cfg = {"zero_cutoff_cm": FREQ_KW["zero_cutoff_cm"]}

    apply_yaml_overrides(
        config_layer_cfg,
        [
            (geom_cfg, (("geom",),)),
            (calc_cfg, (("calc",), ("mlmm",))),
            (opt_cfg, (("opt",),)),
            # The microiteration MM relaxation is an L-BFGS run, so it reads the
            # same `lbfgs` section `opt` does. Without this mapping its memory
            # and step controls were reachable from neither YAML nor CLI.
            (lbfgs_cfg, (("lbfgs",), ("opt", "lbfgs"))),
            (simple_cfg, (("hessian_dimer",),)),
            (rsirfo_cfg, (("rsirfo",),)),
            (frequency_cfg, (("freq",),)),
        ],
    )
    if _is_param_explicit("hessian_calc_mode") and hessian_calc_mode is not None:
        calc_cfg["hessian_calc_mode"] = str(hessian_calc_mode)
    if _is_param_explicit("print_every") and print_every is not None:
        opt_cfg["print_every"] = int(print_every)
    if _is_param_explicit("max_cycles") and max_cycles is not None:
        opt_cfg["max_cycles"] = int(max_cycles)
    if _is_param_explicit("dump"):
        opt_cfg["dump"] = bool(dump)
    if _is_param_explicit("out_dir"):
        opt_cfg["out_dir"] = out_dir
    if _is_param_explicit("thresh") and thresh is not None:
        opt_cfg["thresh"] = str(thresh)
    # --stop-plateau* rides the shared `opt` block, which the RS-I-RFO/RS-P-RFO
    # macro now reads. The dimer's inner L-BFGS takes `hessian_dimer.lbfgs`
    # alone and inherits nothing from `opt`, so it is written here too. The MM
    # micro relaxation is deliberately left out: it never takes the plateau stop.
    _cli_plateau: Dict[str, Any] = {}
    if _is_param_explicit("stop_plateau"):
        _cli_plateau["energy_plateau"] = bool(stop_plateau)
    if stop_plateau_thresh is not None:
        _cli_plateau["energy_plateau_thresh"] = float(stop_plateau_thresh)
    if stop_plateau_window is not None:
        _cli_plateau["energy_plateau_window"] = int(stop_plateau_window)
    if _cli_plateau:
        for _plateau_key, _plateau_val in _cli_plateau.items():
            opt_cfg[_plateau_key] = _plateau_val
    if _is_param_explicit("cli_coord_type") and cli_coord_type is not None:
        geom_cfg["coord_type"] = str(cli_coord_type).lower()
    # Handle --flatten/--no-flatten CLI toggle
    if flatten is not None:
        if flatten:
            # --flatten explicitly enables flattening even when defaults/config disable it.
            default_flatten_iter = int(HESSIAN_DIMER_KW.get("flatten_max_iter", 0))
            if int(simple_cfg.get("flatten_max_iter", 0)) <= 0 and default_flatten_iter > 0:
                simple_cfg["flatten_max_iter"] = default_flatten_iter
        else:
            simple_cfg["flatten_max_iter"] = 0
    if _is_param_explicit("detect_layer"):
        calc_cfg["use_bfactor_layers"] = bool(detect_layer)
    if _is_param_explicit("hess_cutoff") and hess_cutoff is not None:
        calc_cfg["hess_cutoff"] = float(hess_cutoff)
    if _is_param_explicit("movable_cutoff") and movable_cutoff is not None:
        calc_cfg["movable_cutoff"] = float(movable_cutoff)
        calc_cfg["use_bfactor_layers"] = False

    # CLI-resolved charge/spin (from -q / -l derivation, or -m / spin_default)
    # always wins over the CALC_KW default carried in calc_cfg.
    calc_cfg["model_charge"] = int(charge)
    calc_cfg["model_mult"] = int(spin)

    if model_pdb is not None:
        calc_cfg["model_pdb"] = str(model_pdb)
    calc_cfg["input_pdb"] = str(source_path)
    calc_cfg["real_parm7"] = str(real_parm7)

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
    if link_atom_method is not None:
        calc_cfg["link_atom_method"] = str(link_atom_method).lower()
    if mm_backend is not None:
        calc_cfg["mm_backend"] = str(mm_backend).lower()
    if use_cmap is not None:
        calc_cfg["use_cmap"] = use_cmap

    apply_yaml_overrides(
        override_layer_cfg,
        [
            (geom_cfg, (("geom",),)),
            (calc_cfg, (("calc",), ("mlmm",))),
            (opt_cfg, (("opt",),)),
            (lbfgs_cfg, (("lbfgs",), ("opt", "lbfgs"))),
            (simple_cfg, (("hessian_dimer",),)),
            (rsirfo_cfg, (("rsirfo",),)),
            (frequency_cfg, (("freq",),)),
        ],
    )
    def _yaml_has(paths: Tuple[Tuple[str, ...], ...], key: str) -> bool:
        return yaml_section_has_key(config_layer_cfg, paths, key) or yaml_section_has_key(
            override_layer_cfg, paths, key
        )

    cutoff = normalize_frequency_zero_cutoff_cm(
        frequency_cfg["zero_cutoff_cm"]
    )
    cutoff_source = "freq.zero_cutoff_cm"
    cutoff_explicit = _yaml_has((("freq",),), "zero_cutoff_cm")
    for paths, label, value in (
        ((("hessian_dimer",),), "hessian_dimer.neg_freq_thresh_cm", simple_cfg["neg_freq_thresh_cm"]),
        ((("rsirfo",),), "rsirfo.saddle_imaginary_threshold_cm", rsirfo_cfg["saddle_imaginary_threshold_cm"]),
    ):
        key = label.rsplit(".", 1)[1]
        if not _yaml_has(paths, key):
            continue
        alias_value = normalize_frequency_zero_cutoff_cm(value)
        if (cutoff_explicit or cutoff_source != "freq.zero_cutoff_cm") and not np.isclose(
            alias_value, cutoff, rtol=0.0, atol=1e-12
        ):
            raise click.BadParameter(
                f"freq.zero_cutoff_cm and {label} conflict; "
                "set only freq.zero_cutoff_cm."
            )
        cutoff = alias_value
        cutoff_source = label
    frequency_cfg["zero_cutoff_cm"] = cutoff
    simple_cfg["neg_freq_thresh_cm"] = cutoff
    rsirfo_cfg["saddle_imaginary_threshold_cm"] = cutoff

    if use_heavy:
        cli_shared = {
            "max_cycles": "max_cycles",
            "dump": "dump",
            "thresh": "thresh",
            "out_dir": "out_dir",
            "print_every": "print_every",
            "energy_plateau": "stop_plateau",
            "energy_plateau_thresh": "stop_plateau_thresh",
            "energy_plateau_window": "stop_plateau_window",
        }
        for key in sorted(OPT_BASE_KW.keys() & RSIRFO_KW.keys()):
            cli_name = cli_shared.get(key)
            _resolve_shared_optimizer_value(
                opt_cfg,
                rsirfo_cfg,
                key,
                opt_explicit=(
                    _yaml_has((("opt",),), key)
                    or (cli_name is not None and _is_param_explicit(cli_name))
                ),
                downstream_explicit=_yaml_has((("rsirfo",),), key),
                downstream_default=RSIRFO_KW[key],
                downstream_section="rsirfo",
            )
    else:
        _resolve_shared_optimizer_value(
            opt_cfg,
            simple_cfg,
            "thresh",
            opt_explicit=(
                _yaml_has((("opt",),), "thresh")
                or _is_param_explicit("thresh")
            ),
            downstream_explicit=_yaml_has((("hessian_dimer",),), "thresh"),
            downstream_default=HESSIAN_DIMER_KW["thresh"],
            downstream_section="hessian_dimer",
        )
        simple_lbfgs = dict(simple_cfg.get("lbfgs", {}))
        dimer_shared = {
            "print_every": "print_every",
            "energy_plateau": "stop_plateau",
            "energy_plateau_thresh": "stop_plateau_thresh",
            "energy_plateau_window": "stop_plateau_window",
        }
        for key, cli_name in dimer_shared.items():
            _resolve_shared_optimizer_value(
                opt_cfg,
                simple_lbfgs,
                key,
                opt_explicit=(
                    _yaml_has((("opt",),), key)
                    or _is_param_explicit(cli_name)
                ),
                downstream_explicit=_yaml_has(
                    (("hessian_dimer", "lbfgs"),), key
                ),
                downstream_default=hessian_dimer_KW["lbfgs"][key],
                downstream_section="hessian_dimer.lbfgs",
            )
        simple_cfg["lbfgs"] = simple_lbfgs

    if _yaml_has((("hessian_dimer", "lbfgs"),), "max_cycles"):
        raise click.BadParameter(
            "hessian_dimer.lbfgs.max_cycles is not configurable; "
            "use opt.max_cycles."
        )

    # A TS search follows a saddle-search direction, so physical energy is not
    # required to decrease. Keep this invariant after every YAML merge.
    opt_cfg["reject_uphill"] = False
    rsirfo_cfg["reject_uphill"] = False
    simple_cfg["lbfgs"] = _force_ts_reject_uphill_off(
        simple_cfg.get("lbfgs", {})
    )
    partial_hessian_flatten_effective = bool(
        partial_hessian_flatten
        if _is_param_explicit("partial_hessian_flatten")
        else simple_cfg.get("partial_hessian_flatten", True)
    )
    try:
        geom_cfg["tr_projection"] = normalize_tr_projection_mode(
            geom_cfg.get("tr_projection")
        )
    except ValueError as exc:
        prepared_input.cleanup()
        raise click.ClickException(str(exc)) from exc
    if not use_heavy and str(geom_cfg.get("coord_type", "cart")).lower() != "cart":
        click.echo(
            "[tsopt] Gradient/dimer mode uses Cartesian Hessian and mode kernels; "
            "using coord_type=cart."
        )
        geom_cfg["coord_type"] = "cart"
    calc_paths = (("calc",), ("mlmm",))
    partial_explicit = (
        yaml_section_has_key(config_layer_cfg, calc_paths, "return_partial_hessian")
        or yaml_section_has_key(override_layer_cfg, calc_paths, "return_partial_hessian")
    )
    if not partial_explicit:
        calc_cfg["return_partial_hessian"] = True

    # Resolve microiteration config from YAML
    microiter_cfg = dict(MICROITER_KW)
    apply_yaml_overrides(
        config_layer_cfg,
        [(microiter_cfg, (("microiter",),))],
    )
    apply_yaml_overrides(
        override_layer_cfg,
        [(microiter_cfg, (("microiter",),))],
    )

    # Microiteration drives one macro TS step per cycle and works with any of the
    # Hessian TS optimizers (RS-I-RFO / RS-P-RFO / TRIM), which share the
    # optimize/prepare_opt/Bofill-update contract. Only grad/dimer modes lack it.
    use_microiter = bool(microiter) and use_heavy
    if bool(microiter) and not use_heavy:
        click.echo("[microiter] --microiter needs a Hessian TS optimizer (hess/rsirfo/rsprfo/trim); ignoring for grad/dimer.")

    try:
        geom_freeze = _normalize_geom_freeze_opt(geom_cfg.get("freeze_atoms"))
    except click.BadParameter as e:
        click.echo(f"ERROR: {e}", err=True)
        prepared_input.cleanup()
        sys.exit(1)
    geom_cfg["freeze_atoms"] = geom_freeze
    _convert_yaml_layer_atoms_1to0(calc_cfg)
    if freeze_atoms_cli:
        merge_freeze_atom_indices(geom_cfg, freeze_atoms_cli)
    freeze_atoms_final = list(geom_cfg.get("freeze_atoms") or [])
    calc_cfg["freeze_atoms"] = freeze_atoms_final

    # Propagate opt.print_every only when it is explicitly different from the
    # base default. This avoids clobbering optimizer-specific YAML settings
    # (e.g. hessian_dimer.lbfgs.print_every / rsirfo.print_every) with the
    # inherited OPT_BASE default value.
    try:
        pe_opt = int(opt_cfg.get("print_every", OPT_BASE_KW.get("print_every", 100)))
        pe_base = int(OPT_BASE_KW.get("print_every", 100))
        if pe_opt >= 1 and pe_opt != pe_base:
            simple_cfg.setdefault("lbfgs", {})
            simple_cfg["lbfgs"]["print_every"] = pe_opt
            rsirfo_cfg["print_every"] = pe_opt
    except Exception:
        logger.debug("Failed to configure print_every", exc_info=True)

    out_dir_path = Path(opt_cfg["out_dir"]).resolve()

    # movable_cutoff implies full distance-based layer assignment.
    # hess_cutoff alone is allowed with detect-layer and is applied on movable MM atoms.
    detect_layer_enabled = bool(calc_cfg.get("use_bfactor_layers", True))
    model_pdb_cfg = calc_cfg.get("model_pdb")
    if calc_cfg.get("movable_cutoff") is not None:
        if detect_layer_enabled:
            click.echo("[layer] movable_cutoff is set; disabling --detect-layer.", err=True)
        detect_layer_enabled = False
        calc_cfg["use_bfactor_layers"] = False

    layer_source_pdb = source_path
    if detect_layer_enabled and layer_source_pdb.suffix.lower() != ".pdb":
        click.echo("ERROR: --detect-layer requires a PDB input (or --ref-pdb).", err=True)
        prepared_input.cleanup()
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

    if dry_run:
        if model_pdb_cfg is not None:
            model_region_source = "model_pdb"
        elif model_indices:
            model_region_source = "model_indices"
        elif detect_layer_enabled:
            model_region_source = "bfactor"
        else:
            click.echo("ERROR: Provide --model-pdb or --model-indices when B-factor layer detection is disabled in the configuration.", err=True)
            prepared_input.cleanup()
            sys.exit(1)
        if (
            not detect_layer_enabled
            and model_pdb_cfg is None
            and model_indices
            and layer_source_pdb.suffix.lower() != ".pdb"
        ):
            click.echo("ERROR: --model-indices requires a PDB input (or --ref-pdb).", err=True)
            prepared_input.cleanup()
            sys.exit(1)
        click.echo(
            pretty_block(
                "dry_run_plan",
                {
                    "input_geometry": str(geom_input_path),
                    "output_dir": str(out_dir_path),
                    "optimizer_mode": (
                        f"hess-{mode_resolved}" if use_heavy else "grad-dimer"
                    ),
                    "detect_layer": bool(detect_layer_enabled),
                    "model_region_source": model_region_source,
                    "model_indices_count": 0 if not model_indices else len(model_indices),
                    "hessian_calc_mode": calc_cfg.get("hessian_calc_mode"),
                    "partial_hessian_flatten": (
                        partial_hessian_flatten_effective
                    ),
                    "active_dof_mode": str(active_dof_mode),
                    "tr_projection": geom_cfg["tr_projection"],
                    "will_run_tsopt": True,
                    "will_write_summary": True,
                    "backend": calc_cfg.get("backend", "uma"),
                    "reference_mode": (
                        None if reference_mode_path is None else str(reference_mode_path)
                    ),
                    "embedcharge": bool(calc_cfg.get("embedcharge", False)),
                },
            )
        )
        click.echo("[dry-run] Validation complete. TS optimization execution was skipped.")
        prepared_input.cleanup()
        emit(
            format_elapsed("[time] Elapsed Time for TS Opt", time_start),
            narrative=True,
        )
        return

    tsopt_protected_inputs = (
        input_path,
        prepared_input.original_path,
        prepared_input.source_path,
        geom_input_path,
        ref_pdb,
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
        Path(model_pdb_cfg) if model_pdb_cfg is not None else None,
        config_yaml,
        override_yaml,
        (
            Path(calc_cfg["calc_file"])
            if calc_cfg.get("calc_file")
            else None
        ),
        reference_mode_path,
    )
    try:
        _reject_tsopt_output_collisions(
            out_dir_path,
            protected_inputs=tsopt_protected_inputs,
        )
    except _TSOPTOutputCollisionError:
        prepared_input.cleanup()
        raise

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
            protected_inputs=tsopt_protected_inputs,
            echo_fn=click.echo,
        )
    except click.ClickException as exc:
        click.echo(f"ERROR: {exc.message}", err=True)
        prepared_input.cleanup()
        sys.exit(1)
    freeze_atoms_final = apply_layer_freeze_constraints(
        geom_cfg,
        calc_cfg,
        layer_info,
        echo_fn=click.echo,
    )
    _align_three_layer_hessian_targets(calc_cfg, echo_fn=click.echo)

    for key in ("input_pdb", "real_parm7", "model_pdb", "mm_fd_dir"):
        val = calc_cfg.get(key)
        if val:
            calc_cfg[key] = str(Path(val).expanduser().resolve())

    # Pretty-print config summary (only non-default values for concise logging)
    heavy_mode_label = _heavy_mode_label(mode_resolved) if use_heavy else None
    mode_desc = (
        f"{heavy_mode_label} (hess)" if use_heavy else "Dimer (grad)"
    )
    if use_microiter:
        mode_desc += " + Microiteration"

    # Default-verbosity entry summary (skipped in child mode).
    from mlmm.core.utils import calculator_run_label, echo_run_summary
    echo_run_summary({
        "input": str(input_path),
        "backend": calculator_run_label(calc_cfg),
        "opt": mode_desc,
        "out": str(out_dir_path),
    })

    click.echo(f"\n[mode] TS Optimizer: {mode_desc}\n")
    click.echo(pretty_block("geom", format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")))
    echo_calc = format_freeze_atoms_for_echo(filter_calc_for_echo(calc_cfg), key="freeze_atoms")
    click.echo(pretty_block("calc", echo_calc))
    echo_opt = strip_inherited_keys({**opt_cfg, "out_dir": str(out_dir_path)}, OPT_BASE_KW, mode="same")
    click.echo(pretty_block("opt", echo_opt))
    # Show only optimizer-specific settings, not inherited from opt_cfg
    if use_heavy:
        echo_rsirfo = strip_inherited_keys(rsirfo_cfg, opt_cfg)
        click.echo(pretty_block("rsirfo", echo_rsirfo))
    else:
        sd_cfg_for_echo: Dict[str, Any] = {}
        sd_cfg_for_echo["dimer"] = dict(simple_cfg.get("dimer", {}))
        sd_cfg_for_echo["lbfgs"] = strip_inherited_keys(
            dict(simple_cfg.get("lbfgs", {})), opt_cfg
        )
        click.echo(pretty_block("hessian_dimer", sd_cfg_for_echo))

    geometry = None
    try:
        out_dir_path = _prepare_tsopt_output_dir(
            out_dir_path,
            protected_inputs=tsopt_protected_inputs,
        )
        if use_heavy:
            # Heavy mode: RS-I-RFO with full Hessian
            optim_all_path = out_dir_path / "optimization_all_trj.xyz"
            if bool(opt_cfg["dump"]) and optim_all_path.exists():
                optim_all_path.unlink()

            coord_type = geom_cfg.get("coord_type", "cart")
            coord_kwargs = dict(geom_cfg)
            coord_kwargs.pop("coord_type", None)
            geometry = geom_loader(
                geom_input_path,
                coord_type=coord_type,
                **coord_kwargs,
            )
            initial_ts_cart_coords = geometry.cart_coords.copy()

            echo_resolved_device()
            _heavy_optimizer_converged = False
            _heavy_safeguards: Dict[str, Any] = {}
            # additive microiteration serialization (populated on the
            # microiteration path only). Carries the macro/micro leaf outcomes and
            # the separate executed micro-cycle total.
            _heavy_microiteration_obj: Optional[Dict[str, Any]] = None
            _heavy_micro_cycles: Optional[int] = None
            user_max_cycles = optional_positive_int(
                opt_cfg.get("max_cycles"), "opt.max_cycles"
            )
            _heavy_cycle_ledger = _OptimizationCycleLedger(user_max_cycles)
            # Same construction as the microiteration macro step, so the two
            # paths cannot drift apart key by key.
            rsirfo_args = _build_rsirfo_kwargs(
                rsirfo_cfg,
                max_cycles=user_max_cycles,
                out_dir=out_dir_path,
                macro_thresh=thresh,
                mode=mode_resolved,
                opt_cfg=opt_cfg,
                dump=bool(opt_cfg["dump"]),
                reference_mode=reference_mode,
                flatten_enabled=bool(
                    int(simple_cfg.get("flatten_max_iter", 0)) > 0
                ),
            )

            if use_microiter:
                # --- Microiteration path ---
                microiter_outcome = _run_microiter_tsopt(
                    geometry,
                    calc_cfg,
                    rsirfo_cfg,
                    lbfgs_cfg,
                    opt_cfg,
                    microiter_cfg,
                    out_dir_path,
                    dump=bool(opt_cfg["dump"]),
                    thresh=thresh,
                    mode=mode_resolved,
                    reference_mode=reference_mode,
                    flatten_enabled=bool(
                        int(simple_cfg.get("flatten_max_iter", 0)) > 0
                    ),
                )
                _heavy_optimizer_converged = bool(microiter_outcome["converged"])
                _heavy_safeguards = dict(microiter_outcome["safeguards"])
                last_optimizer = microiter_outcome["optimizer"]
                _mi_outcome_obj = microiter_outcome.get("outcome")
                if _mi_outcome_obj is not None:
                    _heavy_microiteration_obj = _mi_outcome_obj.to_result_object()
                    _heavy_micro_cycles = int(microiter_outcome.get("micro_cycles", 0))

                # Write final geometry
                final_xyz = out_dir_path / "final_geometry.xyz"
                final_xyz.write_text(geometry.as_xyz(), encoding="utf-8")

                initial_run_cycles = int(microiter_outcome["cycles"])
                _heavy_cycle_ledger.debit(initial_run_cycles)
            else:
                # --- Standard RS-I-RFO path ---
                base_calc = mlmm(**calc_cfg)
                geometry.set_calculator(base_calc)

                click.echo("[tsopt] Seeding initial Hessian via shared freq backend.")
                hess_device = _torch_device(simple_cfg.get("device", calc_cfg.get("ml_device", "auto")))
                h_init = _calc_full_hessian_torch(geometry, calc_cfg, hess_device)
                geometry.cart_hessian = h_init
                click.echo(
                    f"[tsopt] Initial Hessian seeded (shape={h_init.shape[0]}x{h_init.shape[1]})."
                )
                del h_init

                calc_core = base_calc.core if hasattr(base_calc, "core") else base_calc
                hess_active_atoms = list(getattr(calc_core, "hess_active_atoms", []))

                optimizer = TSOPT_CLASS_MAP[mode_resolved](geometry, **rsirfo_args)
                optimizer.run()
                last_optimizer = optimizer
                _heavy_optimizer_converged = bool(optimizer.is_converged)
                _heavy_safeguards = _optimizer_safeguard_payload(optimizer)
                emit_optimizer_terminal_status(
                    "tsopt",
                    converged=getattr(optimizer, "is_converged", None),
                    cycles=optimizer_cycle_count(optimizer),
                    max_cycles=opt_cfg.get("max_cycles"),
                    stalled=getattr(optimizer, "is_stalled", False),
                    stop_reason=getattr(optimizer, "stop_reason", None) or None,
                    converged_message="Numerical optimization converged.",
                )
                if bool(opt_cfg["dump"]):
                    _append_xyz_trajectory(optim_all_path, out_dir_path / "optimization_trj.xyz")

                # --- Post-RSIRFO: count imaginary modes and optional flatten loop ---
                # Save cycle count before deleting optimizer for budget check.
                initial_run_cycles = int(optimizer_cycle_count(optimizer) or 0)
                _heavy_cycle_ledger.debit(initial_run_cycles)
                geometry.set_calculator(None)
                del calc_core
                del base_calc
                _clear_cuda_cache()
            hessian_postprocessing_ready = _hessian_postprocessing_is_ready(
                last_optimizer
            )
            _do_final_freq = hessian_postprocessing_ready and (
                not skip_final_freq or getattr(last_optimizer, "is_stalled", False)
            )
            if skip_final_freq and not getattr(last_optimizer, "is_stalled", False):
                click.echo(
                    "[tsopt] WARNING: TS saddle-point order is not verified "
                    "(--skip-final-freq).",
                    err=True,
                )
            mlmm_kwargs_for_heavy = _post_analysis_hessian_config(
                calc_cfg,
                partial=partial_hessian_flatten_effective,
            )
            device = _torch_device(simple_cfg.get("device", calc_cfg.get("ml_device", "auto")))

            # Determine active atoms for frequency analysis based on --active-dof-mode.
            active_atoms_freq = _resolve_validated_hessian_analysis_atoms(
                calc_cfg,
                len(geometry.atomic_numbers),
                active_dof_mode,
                freeze_atoms_final,
                validate_coverage=_do_final_freq,
            )
            n_atoms = len(geometry.atomic_numbers)
            requested_active_atoms = list(active_atoms_freq)

            rigid_projection_info: Dict[str, Any] = {}

            def _store_ts_hessian(
                H: Any,
                *,
                source: str,
                energy_ha: Optional[float] = None,
            ) -> None:
                """Publish the terminal exact Hessian for the following IRC.

                Optimizer-owned and workflow-owned PHVA data share this path,
                so reusing the terminal diagnostics cannot accidentally omit
                the raw cache artifact and force IRC to recompute it.
                """
                from mlmm.io.hessian_cache import (
                    store as _hess_store,
                    identity_from_context as _hess_identity,
                )

                computed_atoms = _ordered_hessian_coverage_atoms(
                    geometry, len(geometry.atomic_numbers)
                )
                active_dofs = None
                if computed_atoms is not None:
                    active_dofs = [
                        3 * atom + axis
                        for atom in computed_atoms
                        for axis in range(3)
                    ]
                meta = {
                    "cart_coords": geometry.cart_coords,
                    "source": source,
                }
                if energy_ha is not None:
                    meta["energy_ha"] = energy_ha
                _hess_store(
                    "ts",
                    H,
                    active_dofs=active_dofs,
                    meta=meta,
                    identity=_hess_identity(geometry, calc_cfg, role="ts"),
                )

            def _calc_freqs_and_modes() -> Tuple[np.ndarray, torch.Tensor]:
                H, _energy_ha = _freq_calc_full_hessian_torch(
                    geometry, mlmm_kwargs_for_heavy, device, refresh_geom_meta=True,
                )
                _store_ts_hessian(
                    H,
                    source="tsopt_exact",
                    energy_ha=_energy_ha,
                )
                _raw_shape = tuple(H.shape)
                H_analysis, _analysis_atoms, _coverage_atoms, _storage = (
                    _reconcile_hessian_analysis_basis(
                        H,
                        geometry,
                        requested_active_atoms,
                    )
                )
                del H
                freqs_local, modes_gpu = _modes_from_Hact_embedded(
                    H_analysis,
                    geometry.atomic_numbers,
                    geometry.cart_coords.reshape(-1, 3),
                    _analysis_atoms,
                    device,
                    tr_projection=geom_cfg["tr_projection"],
                    projection_info=rigid_projection_info,
                    frequency_zero_cutoff_cm=neg_freq_thresh_cm,
                )
                modes_local = modes_gpu.detach().cpu()
                del modes_gpu
                rigid_projection_info.update({
                    "hessian_space": (
                        "full" if len(_analysis_atoms) == n_atoms else "active"
                    ),
                    "analysis_hessian_shape": list(H_analysis.shape),
                    "raw_hessian_shape": list(_raw_shape),
                    "computed_atom_count": len(_coverage_atoms),
                    "analysis_atom_count": len(_analysis_atoms),
                    "storage": _storage,
                    "source": "tsopt_exact",
                })
                del H_analysis
                _clear_cuda_cache()
                return freqs_local, modes_local

            def _terminal_freqs_and_modes(
                current_optimizer: Any,
            ) -> Tuple[np.ndarray, torch.Tensor]:
                cached = _optimizer_exact_frequency_data(
                    current_optimizer, geometry
                )
                if cached is not None:
                    (
                        freqs_local,
                        modes_local,
                        projection_local,
                        exact_hessian,
                    ) = cached
                    if exact_hessian is not None:
                        _store_ts_hessian(
                            exact_hessian,
                            source="optimizer_terminal_exact_phva",
                        )
                    rigid_projection_info.clear()
                    rigid_projection_info.update(projection_local)
                    emit(
                        "[hessian] Reused terminal exact PHVA; no duplicate "
                        "Hessian calculation was performed.",
                        narrative=True,
                    )
                    return freqs_local, modes_local
                return _calc_freqs_and_modes()

            def _optimizer_safeguards(opt) -> Dict[str, Any]:
                return _optimizer_safeguard_payload(opt)

            def _run_path_restart(
                restart_reference: Optional[np.ndarray],
            ) -> Tuple[
                Any, bool, Dict[str, Any], int, Optional[Dict[str, Any]], Optional[int]
            ]:
                """Run the selected Hessian optimizer from current coordinates.

                Returns ``(optimizer, converged, safeguards, cycles,
                microiteration_obj, micro_cycles)``.  On the microiteration path the
                last two carry THIS restart's own additive ``microiteration`` result
                object and executed micro-cycle total so the serialized block
                describes the run that produced the final geometry (never a
                superseded initial run); on the ordinary path both are ``None``.
                """
                remaining_cycles = _heavy_cycle_ledger.remaining
                if remaining_cycles is not None and remaining_cycles <= 0:
                    raise OptimizationError(
                        "Command-level --max-cycles budget exhausted."
                    )
                if use_microiter:
                    restart_opt_cfg = dict(opt_cfg)
                    restart_opt_cfg["max_cycles"] = remaining_cycles
                    outcome = _run_microiter_tsopt(
                        geometry,
                        calc_cfg,
                        rsirfo_cfg,
                        lbfgs_cfg,
                        restart_opt_cfg,
                        microiter_cfg,
                        out_dir_path,
                        dump=bool(opt_cfg["dump"]),
                        thresh=thresh,
                        mode=mode_resolved,
                        reference_mode=restart_reference,
                        flatten_enabled=bool(
                            int(simple_cfg.get("flatten_max_iter", 0)) > 0
                        ),
                    )
                    restart_optimizer = outcome["optimizer"]
                    restart_micro_obj, restart_micro_cycles = (
                        _restart_microiteration_carry(outcome)
                    )
                    restart_cycles = int(outcome["cycles"])
                    _heavy_cycle_ledger.debit(restart_cycles)
                    return (
                        restart_optimizer,
                        bool(outcome["converged"]),
                        dict(outcome["safeguards"]),
                        restart_cycles,
                        restart_micro_obj,
                        restart_micro_cycles,
                    )

                restart_calc = mlmm(**calc_cfg)
                geometry.set_calculator(restart_calc)
                restart_hessian = _calc_full_hessian_torch(
                    geometry, calc_cfg, device
                )
                geometry.cart_hessian = restart_hessian
                del restart_hessian
                restart_args = dict(rsirfo_args)
                restart_args["max_cycles"] = remaining_cycles
                if restart_reference is not None:
                    restart_args["reference_mode"] = restart_reference
                restart_optimizer = TSOPT_CLASS_MAP[mode_resolved](
                    geometry, **restart_args
                )
                restart_optimizer.run()
                emit_optimizer_terminal_status(
                    "tsopt",
                    converged=getattr(restart_optimizer, "is_converged", None),
                    cycles=optimizer_cycle_count(restart_optimizer),
                    max_cycles=remaining_cycles,
                    stalled=getattr(restart_optimizer, "is_stalled", False),
                    stop_reason=getattr(restart_optimizer, "stop_reason", None) or None,
                    converged_message="Numerical optimization converged.",
                )
                cycles = int(optimizer_cycle_count(restart_optimizer) or 0)
                _heavy_cycle_ledger.debit(cycles)
                converged = bool(restart_optimizer.is_converged)
                safeguards = _optimizer_safeguards(restart_optimizer)
                geometry.set_calculator(None)
                del restart_calc
                _clear_cuda_cache()
                return restart_optimizer, converged, safeguards, cycles, None, None

            hessian_error: Optional[str] = getattr(
                last_optimizer, "_last_exact_failure_reason", None
            )
            if not _do_final_freq:
                freqs_cm, modes = None, None
            try:
                if _do_final_freq:
                    freqs_cm, modes = _terminal_freqs_and_modes(last_optimizer)
            except Exception as exc:
                hessian_error = f"{type(exc).__name__}: {exc}"
                click.echo("[tsopt] ERROR: Terminal PHVA failed.", err=True)
                emit(
                    f"[tsopt] Terminal PHVA error: {hessian_error}",
                    detail=True,
                )
                _clear_cuda_cache()
                freqs_cm, modes = None, None
            neg_freq_thresh_cm = float(simple_cfg.get("neg_freq_thresh_cm", 5.0))

            if freqs_cm is not None and modes is not None:
                neg_mask = freqs_cm < -abs(neg_freq_thresh_cm)
                n_imag = int(np.sum(neg_mask))
                ims = [float(x) for x in freqs_cm if x < -abs(neg_freq_thresh_cm)]
                emit(f"[Imaginary modes] n={n_imag} ({ims})", narrative=True)
                _warn_if_leading_imaginary_mode_is_soft(ims)

                saddle_multistart_attempts: List[Dict[str, Any]] = []
                target_mode_is_negative = getattr(
                    last_optimizer,
                    "_last_exact_target_mode_is_negative",
                    None,
                )
                if (
                    int(rsirfo_args.get("saddle_recovery_max_cycles", 0)) > 0
                    and not _heavy_optimizer_converged
                    and reference_mode is not None
                    and (n_imag <= 1 or target_mode_is_negative is False)
                ):
                    baseline = {
                        "coords": geometry.cart_coords.copy(),
                        "optimizer": last_optimizer,
                        "freqs": freqs_cm.copy(),
                        "modes": modes.detach().cpu().clone(),
                        "converged": _heavy_optimizer_converged,
                        "safeguards": dict(_heavy_safeguards),
                        "cycles": initial_run_cycles,
                        # The additive microiteration block that describes the
                        # initial run (the baseline geometry), carried so a
                        # baseline selection re-anchors it to the selected mode.
                        "microiteration_obj": _heavy_microiteration_obj,
                        "micro_cycles": _heavy_micro_cycles,
                    }
                    best_path_negative = None
                    multistart_success = False
                    multistart_budget_exhausted = False
                    for mode_source, restart_unit in _path_restart_mode_candidates(
                        last_optimizer,
                        geometry,
                        reference_modes,
                        reference_mode_labels,
                    ):
                        if mode_source == "initial-soft-root" and best_path_negative is not None:
                            break
                        for amplitude_ang in PATH_MODE_RESTART_AMPLITUDES_ANG:
                            if (
                                _heavy_cycle_ledger.remaining is not None
                                and _heavy_cycle_ledger.remaining <= 0
                            ):
                                multistart_budget_exhausted = True
                                click.echo(
                                    "[tsopt] Reached --max-cycles budget; "
                                    "stopping path-mode restarts."
                                )
                                break
                            geometry.cart_coords = (
                                initial_ts_cart_coords
                                + (amplitude_ang / BOHR2ANG) * restart_unit
                            )
                            geometry.set_calculator(None)
                            emit(
                                "[path-mode restart] "
                                f"source={mode_source}, displacement={amplitude_ang:+.2f} Å",
                                narrative=True,
                            )
                            (
                                restart_optimizer,
                                restart_converged,
                                restart_safeguards,
                                restart_cycles,
                                restart_micro_obj,
                                restart_micro_cycles,
                            ) = _run_path_restart(restart_unit)
                            last_optimizer = restart_optimizer
                            _heavy_optimizer_converged = restart_converged
                            _heavy_safeguards = restart_safeguards
                            # Re-anchor the additive microiteration block on the run
                            # that produced the current geometry (the last restart).
                            _heavy_microiteration_obj = restart_micro_obj
                            _heavy_micro_cycles = restart_micro_cycles
                            geometry.set_calculator(None)
                            restart_freqs, restart_modes = _terminal_freqs_and_modes(restart_optimizer)
                            restart_n_imag = int(
                                np.sum(
                                    restart_freqs
                                    < -abs(neg_freq_thresh_cm)
                                )
                            )
                            target_negative = getattr(
                                restart_optimizer,
                                "_last_exact_target_mode_is_negative",
                                None,
                            )
                            attempt = {
                                "mode_source": mode_source,
                                "displacement_ang": amplitude_ang,
                                "converged": restart_converged,
                                "n_imaginary": restart_n_imag,
                                "target_mode_negative": target_negative,
                            }
                            saddle_multistart_attempts.append(attempt)
                            emit(
                                "[path-mode restart] "
                                f"converged={restart_converged}, n_imag={restart_n_imag}",
                                narrative=True,
                            )
                            if restart_converged and restart_n_imag == 1:
                                freqs_cm, modes = restart_freqs, restart_modes
                                n_imag = restart_n_imag
                                multistart_success = True
                                break
                            if target_negative is True:
                                force = (
                                    restart_optimizer.forces[-1]
                                    if restart_optimizer.forces
                                    else geometry.cart_forces
                                )
                                if isinstance(force, torch.Tensor):
                                    force = force.detach().cpu().numpy()
                                score = (
                                    abs(restart_n_imag - 1),
                                    float(np.max(np.abs(np.asarray(force, dtype=float)))),
                                )
                                if best_path_negative is None or score < best_path_negative["score"]:
                                    best_path_negative = {
                                        "score": score,
                                        "coords": geometry.cart_coords.copy(),
                                        "optimizer": restart_optimizer,
                                        "freqs": restart_freqs.copy(),
                                        "modes": restart_modes.detach().cpu().clone(),
                                        "converged": restart_converged,
                                        "safeguards": dict(restart_safeguards),
                                        "cycles": restart_cycles,
                                        "microiteration_obj": restart_micro_obj,
                                        "micro_cycles": restart_micro_cycles,
                                    }
                        if multistart_success:
                            break
                        if multistart_budget_exhausted:
                            break

                    if not multistart_success:
                        selected = best_path_negative or baseline
                        geometry.cart_coords = selected["coords"]
                        last_optimizer = selected["optimizer"]
                        freqs_cm = selected["freqs"]
                        modes = selected["modes"]
                        _heavy_optimizer_converged = bool(selected["converged"])
                        _heavy_safeguards = dict(selected["safeguards"])
                        # The selected run (best-path-negative or baseline) owns the
                        # additive microiteration block for the final geometry.
                        _heavy_microiteration_obj = selected["microiteration_obj"]
                        _heavy_micro_cycles = selected["micro_cycles"]
                        n_imag = int(
                            np.sum(freqs_cm < -abs(neg_freq_thresh_cm))
                        )
                    (out_dir_path / "final_geometry.xyz").write_text(
                        geometry.as_xyz(), encoding="utf-8"
                    )
                if saddle_multistart_attempts:
                    _heavy_safeguards["path_mode_restarts"] = saddle_multistart_attempts

                flatten_max_iter = int(simple_cfg.get("flatten_max_iter", 0))
                target_mode_is_negative = getattr(
                    last_optimizer,
                    "_last_exact_target_mode_is_negative",
                    None,
                )
                flatten_max_iter, flatten_vetoed = _effective_flatten_iterations(
                    flatten_max_iter,
                    has_reference_mode=reference_mode is not None,
                    n_imag=n_imag,
                    target_mode_is_negative=target_mode_is_negative,
                )
                if flatten_vetoed:
                    _flatten_skip_reason = (
                        "target mode sign never determined"
                        if target_mode_is_negative is None
                        else "target mode is not negative"
                    )
                    click.echo(
                        "[flatten] Skipping extra-mode flattening: "
                        f"{_flatten_skip_reason}.",
                        err=True,
                    )
                budget_remaining = (
                    _heavy_cycle_ledger.remaining is None
                    or _heavy_cycle_ledger.remaining > 0
                )

                if flatten_max_iter > 0 and n_imag > 1 and not budget_remaining:
                    _flatten_skip_reason = (
                        "max-cycles budget exhausted before flattening"
                    )
                    click.echo("[tsopt] Reached --max-cycles budget; skipping flatten loop.")
                elif flatten_max_iter > 0 and n_imag > 1 and budget_remaining:
                    click.echo(
                        "[flatten] Extra imaginary modes detected; starting "
                        f"{heavy_mode_label} flatten loop."
                    )
                    masses_amu = np.array([atomic_masses[z] for z in geometry.atomic_numbers])
                    main_root = int(simple_cfg.get("root", 0))

                    def _run_flatten_branch(
                        start_coords: np.ndarray,
                        branch_reference: Optional[np.ndarray],
                        label: str,
                    ) -> Dict[str, Any]:
                        geometry.cart_coords = np.asarray(start_coords, dtype=float).copy()
                        (
                            branch_optimizer,
                            converged,
                            safeguards,
                            cycles,
                            branch_micro_obj,
                            branch_micro_cycles,
                        ) = _run_path_restart(branch_reference)
                        geometry.set_calculator(None)
                        branch_ready = _hessian_postprocessing_is_ready(
                            branch_optimizer
                        )
                        if branch_ready:
                            branch_freqs, branch_modes = _terminal_freqs_and_modes(
                                branch_optimizer
                            )
                        else:
                            branch_freqs = np.asarray([], dtype=float)
                            branch_modes = torch.empty(
                                (0, geometry.cart_coords.size), dtype=torch.float64
                            )
                        branch_n_imag = int(
                            np.sum(branch_freqs < -abs(neg_freq_thresh_cm))
                        )
                        emit(
                            f"[Imaginary modes:{label}] n={branch_n_imag}",
                            narrative=True,
                        )
                        return {
                            "label": label,
                            "optimizer": branch_optimizer,
                            "coords": geometry.cart_coords.copy(),
                            "freqs": branch_freqs,
                            "modes": branch_modes,
                            "n_imag": branch_n_imag,
                            "converged": converged,
                            "safeguards": safeguards,
                            "cycles": cycles,
                            "microiteration_obj": branch_micro_obj,
                            "micro_cycles": branch_micro_cycles,
                        }

                    def _flatten_branch_score(
                        result: Dict[str, Any],
                    ) -> Tuple[int, int, int, float, float]:
                        branch_optimizer = result["optimizer"]
                        target_negative = (
                            getattr(
                                branch_optimizer,
                                "_last_exact_target_mode_is_negative",
                                None,
                            )
                            is True
                        )
                        branch_freqs = np.asarray(result["freqs"], dtype=float)
                        negative_indices = np.flatnonzero(
                            branch_freqs < -abs(neg_freq_thresh_cm)
                        )
                        target_index = getattr(
                            branch_optimizer,
                            "_last_exact_target_mode_index",
                            None,
                        )
                        surplus_strength = sum(
                            abs(float(branch_freqs[int(mode_index)]))
                            for mode_index in negative_indices
                            if not (
                                target_negative
                                and target_index is not None
                                and int(mode_index) == int(target_index)
                            )
                        )
                        force = (
                            branch_optimizer.forces[-1]
                            if branch_optimizer.forces
                            else np.asarray([], dtype=float)
                        )
                        if isinstance(force, torch.Tensor):
                            force = force.detach().cpu().numpy()
                        force = np.asarray(force, dtype=float).reshape(-1)
                        return (
                            0 if result.get("converged", False) else 1,
                            0 if target_negative else 1,
                            abs(int(result["n_imag"]) - 1),
                            surplus_strength,
                            float(np.max(np.abs(force))) if force.size else float("inf"),
                        )

                    for it in range(flatten_max_iter):
                        if (
                            _heavy_cycle_ledger.remaining is not None
                            and _heavy_cycle_ledger.remaining <= 0
                        ):
                            _flatten_skip_reason = (
                                "max-cycles budget exhausted during flattening"
                            )
                            click.echo(
                                "[tsopt] Reached --max-cycles budget; "
                                "stopping flatten loop."
                            )
                            break
                        click.echo(
                            f"[flatten] {heavy_mode_label} iteration "
                            f"{it + 1}/{flatten_max_iter}"
                        )
                        flatten_reference_mode = _transported_path_mode_full(
                            last_optimizer, geometry, reference_mode
                        )
                        pre_flatten_coords = geometry.cart_coords.copy()
                        did_flatten = _flatten_once_with_modes_for_geom(
                            geometry,
                            masses_amu,
                            mlmm_kwargs_for_heavy,
                            freqs_cm,
                            modes,
                            neg_freq_thresh_cm,
                            float(simple_cfg.get("flatten_amp_ang", 0.10)),
                            float(simple_cfg.get("flatten_sep_cutoff", 0.0)),
                            int(simple_cfg.get("flatten_k", 10)),
                            main_root,
                            reference_mode=flatten_reference_mode,
                        )
                        if not did_flatten:
                            _flatten_skip_reason = (
                                "no eligible extra imaginary modes"
                            )
                            click.echo("[flatten] No eligible modes to flatten; stopping.")
                            break

                        try:
                            primary_start = geometry.cart_coords.copy()
                            primary_result = _run_flatten_branch(
                                primary_start,
                                flatten_reference_mode,
                                "primary",
                            )
                            selected_result = primary_result
                            if (
                                reference_mode is not None
                                and _flatten_branch_needs_alternate(primary_result)
                            ):
                                primary_score = _flatten_branch_score(primary_result)
                                if (
                                    _heavy_cycle_ledger.remaining is None
                                    or _heavy_cycle_ledger.remaining > 0
                                ):
                                    alternate_result = _run_flatten_branch(
                                        _mirrored_flatten_start(
                                            pre_flatten_coords, primary_start
                                        ),
                                        flatten_reference_mode,
                                        "alternate",
                                    )
                                    alternate_score = _flatten_branch_score(
                                        alternate_result
                                    )
                                    if alternate_score < primary_score:
                                        selected_result = alternate_result
                                    emit(
                                        "[flatten] Signed-branch probe selected "
                                        f"{selected_result['label']} "
                                        f"(primary={primary_score}, "
                                        f"alternate={alternate_score}).",
                                        narrative=True,
                                    )
                                else:
                                    click.echo(
                                        "[flatten] Skipping alternate signed branch: "
                                        "--max-cycles budget exhausted."
                                    )
                        except Exception as exc:
                            is_oom = isinstance(exc, torch.OutOfMemoryError) or ("cuda out of memory" in str(exc).lower())
                            if is_oom:
                                click.echo(
                                    "[tsopt] WARNING: CUDA OOM during final frequency analysis; "
                                    "stopping flatten loop.",
                                    err=True,
                                )
                                _clear_cuda_cache()
                                freqs_cm, modes = None, None
                                _flatten_skip_reason = (
                                    "GPU memory exhausted during flattening"
                                )
                                # Roll back to the last committed state before leaving the loop.
                                # `geometry` still holds this iteration's flatten displacement,
                                # which is an uncommitted, unvalidated branch candidate: branch
                                # selection only commits via `geometry.cart_coords =
                                # selected_result["coords"]` after the frequency check below.
                                # Without this, the reported energy and final_geometry.* would
                                # describe that uncommitted candidate.
                                geometry.cart_coords = pre_flatten_coords
                                break
                            raise

                        last_optimizer = selected_result["optimizer"]
                        hessian_postprocessing_ready = _hessian_postprocessing_is_ready(
                            last_optimizer
                        )
                        geometry.cart_coords = selected_result["coords"]
                        freqs_cm = selected_result["freqs"]
                        modes = selected_result["modes"]
                        n_imag = int(selected_result["n_imag"])
                        _heavy_optimizer_converged = bool(selected_result["converged"])
                        _heavy_safeguards = dict(selected_result["safeguards"])
                        # Re-anchor the additive microiteration block on the selected
                        # flatten branch (the run that produced the final geometry).
                        _heavy_microiteration_obj = selected_result["microiteration_obj"]
                        _heavy_micro_cycles = selected_result["micro_cycles"]
                        ims = [float(x) for x in freqs_cm if x < -abs(neg_freq_thresh_cm)]
                        emit(f"[Imaginary modes] n={n_imag} ({ims})", narrative=True)
                        _warn_if_leading_imaginary_mode_is_soft(ims)
                        (out_dir_path / "final_geometry.xyz").write_text(
                            geometry.as_xyz(), encoding="utf-8"
                        )
                        if not hessian_postprocessing_ready:
                            freqs_cm, modes = None, None
                            break
                        if (
                            reference_mode is not None
                            and getattr(
                                last_optimizer,
                                "_last_exact_target_mode_is_negative",
                                None,
                            )
                            is False
                        ):
                            click.echo(
                                "[flatten] Path-correlated mode was lost; stopping "
                                "without preserving an unrelated negative mode.",
                                err=True,
                            )
                            break
                        if n_imag <= 1:
                            break

            if freqs_cm is not None and modes is not None:
                # --- Write all final imaginary modes like light mode ---
                vib_dir = out_dir_path / "vib"
                vib_dir.mkdir(parents=True, exist_ok=True)
                _ref_pdb_for_modes = source_path if source_path.suffix.lower() == ".pdb" else None
                n_written = _write_all_imag_modes(
                    geometry,
                    freqs_cm,
                    modes,
                    neg_freq_thresh_cm,
                    vib_dir,
                    ref_pdb=_ref_pdb_for_modes,
                )
                _export_n_imag = _certified_saddle_order(
                    freqs_cm, neg_freq_thresh_cm
                )
                _export_message, _export_message_is_diagnostic = (
                    _dimer_mode_export_message(
                        n_written,
                        _export_n_imag,
                        neg_freq_thresh_cm,
                        float(np.min(freqs_cm)),
                    )
                )
                click.echo(_export_message, err=_export_message_is_diagnostic)
                if n_written:
                    emit(f"[DONE] Mode files → {vib_dir}", detail=True)
            elif _do_final_freq:
                emit(
                    "[tsopt] Skipped final imaginary-mode trajectory after "
                    "frequency-analysis fallback.",
                    detail=True,
                )

            # Capture freq/energy data for result.json BEFORE deleting
            _heavy_imag_freqs: Optional[list] = (
                None
                if not skip_final_freq and not hessian_postprocessing_ready
                else []
            )
            _heavy_n_imag: Optional[int] = None
            _heavy_energy = None
            if freqs_cm is not None:
                # Use the same configured magnitude gate for the console list,
                # exported modes, result metadata, and saddle-order verdict.
                _heavy_imag_freqs = _certified_negative_frequencies(
                    freqs_cm, neg_freq_thresh_cm
                )
                _heavy_n_imag = len(_heavy_imag_freqs)
                if _heavy_n_imag != 1:
                    click.echo(
                        _unexpected_saddle_order_message(_heavy_n_imag),
                        err=True,
                    )
            # a stall (energy-plateau outcome of the selected optimizer)
            # wins over every convergence/saddle-order verdict — it is never a
            # converged saddle. n_imag / stop_reason are recorded separately so
            # a stall does not hide saddle-order evidence.
            _heavy_stalled = bool(
                'last_optimizer' in dir()
                and getattr(last_optimizer, "is_stalled", False)
            )
            _heavy_status = _heavy_ts_terminal_status(
                optimizer_converged=_heavy_optimizer_converged,
                n_imag=_heavy_n_imag,
                stalled=_heavy_stalled,
            )
            _heavy_saddle_validation = _saddle_validation_from_count(
                _heavy_n_imag
            )
            _heavy_reaction_mode_index = None
            _heavy_reaction_mode_frequency = None
            _heavy_reaction_mode_overlap = None
            _heavy_candidate_mode_index = getattr(
                last_optimizer, "_last_exact_target_mode_index", None
            )
            if (
                _heavy_candidate_mode_index is not None
                and freqs_cm is not None
                and 0 <= int(_heavy_candidate_mode_index) < len(freqs_cm)
                and float(freqs_cm[int(_heavy_candidate_mode_index)]) < 0.0
            ):
                _heavy_reaction_mode_index = int(_heavy_candidate_mode_index)
                _heavy_reaction_mode_frequency = float(
                    freqs_cm[_heavy_reaction_mode_index]
                )
                _heavy_reaction_mode_overlap = getattr(
                    last_optimizer, "_last_exact_target_mode_overlap", None
                )
            elif freqs_cm is not None and len(freqs_cm):
                _negative_indices = np.flatnonzero(np.asarray(freqs_cm) < 0.0)
                if _negative_indices.size:
                    _heavy_reaction_mode_index = int(_negative_indices[0])
                    _heavy_reaction_mode_frequency = float(
                        freqs_cm[_heavy_reaction_mode_index]
                    )
            projection_block = pretty_block(
                "rigid_projection", rigid_projection_info
            )
            if projection_block:
                click.echo(projection_block)
            # Restart branches share the Hessian cache. Evaluate the selected
            # final coordinates directly so an unselected branch cannot leak
            # its cached energy into result.json.
            try:
                _heavy_energy = _calc_energy(geometry, calc_cfg)
            except (RuntimeError, ValueError, AttributeError, KeyError) as _ee:
                logger.warning(
                    "Heavy-mode final energy evaluation failed: %s. "
                    "Reporting status='energy_missing' with NaN energy.",
                    _ee,
                )
                _heavy_energy = float("nan")

            if modes is not None:
                del modes
            if freqs_cm is not None:
                del freqs_cm
            _clear_cuda_cache()

            # Ensure final_geometry.xyz exists (partial-micro path may not write it).
            final_xyz = out_dir_path / "final_geometry.xyz"
            if not final_xyz.exists():
                final_xyz.write_text(geometry.as_xyz(), encoding="utf-8")

        else:
            # Light mode: Partial Hessian guided Dimer
            light_n_atoms = len(ase_read(str(geom_input_path), index=0))
            light_active_atoms = _resolve_validated_hessian_analysis_atoms(
                calc_cfg,
                light_n_atoms,
                active_dof_mode,
                freeze_atoms_final,
                validate_coverage=not skip_final_freq,
            )
            runner = HessianDimer(
                fn=str(geom_input_path),
                out_dir=str(out_dir_path),
                thresh_loose=simple_cfg.get("thresh_loose", "gau_loose"),
                thresh=simple_cfg.get("thresh", "baker"),
                update_interval_hessian=int(simple_cfg.get("update_interval_hessian", 500)),
                neg_freq_thresh_cm=float(simple_cfg.get("neg_freq_thresh_cm", 5.0)),
                flatten_amp_ang=float(simple_cfg.get("flatten_amp_ang", 0.10)),
                flatten_max_iter=int(simple_cfg.get("flatten_max_iter", 50)),
                mem=int(simple_cfg.get("mem", 100000)),
                use_lobpcg=bool(simple_cfg.get("use_lobpcg", True)),
                calc_kwargs=dict(calc_cfg),
                device=str(simple_cfg.get("device", calc_cfg.get("ml_device", "auto"))),
                dump=bool(opt_cfg["dump"]),
                root=int(simple_cfg.get("root", 0)),
                dimer_kwargs=dict(simple_cfg.get("dimer", {})),
                lbfgs_kwargs=dict(simple_cfg.get("lbfgs", {})),
                max_total_cycles=opt_cfg.get("max_cycles"),
                geom_kwargs=dict(geom_cfg),
                # `simple_cfg` starts from HESSIAN_DIMER_KW, so these keys are ALWAYS present and
                # a `.get(key, <cli value>)` default can never be reached — the CLI flag would be
                # dead and only YAML would work. Apply the documented precedence explicitly.
                partial_hessian_flatten=partial_hessian_flatten_effective,
                flatten_sep_cutoff=float(simple_cfg.get("flatten_sep_cutoff", 0.0)),
                flatten_k=int(simple_cfg.get("flatten_k", 10)),
                flatten_loop_bofill=bool(simple_cfg.get("flatten_loop_bofill", False)),
                ml_only_hessian_dimer=bool(
                    ml_only_hessian_dimer
                    if _is_param_explicit("ml_only_hessian_dimer")
                    else simple_cfg.get("ml_only_hessian_dimer", False)
                ),
                analysis_active_atoms=light_active_atoms,
                source_path=source_path,
                skip_final_freq=skip_final_freq,
            )

            echo_resolved_device()

            runner.run()
            _flatten_skip_reason = runner.flatten_skip_reason
            emit_optimizer_terminal_status(
                "tsopt",
                converged=getattr(runner, "is_converged", None),
                cycles=getattr(runner, "_cycles_spent", None),
                max_cycles=opt_cfg.get("max_cycles"),
                stalled=getattr(runner, "is_stalled", False),
                stop_reason=getattr(runner, "stop_reason", None) or None,
                converged_message="Numerical optimization converged.",
            )

        if is_convert_file_enabled() and source_path.suffix.lower() == ".pdb":
            ref_pdb = source_path.resolve()
            final_xyz = out_dir_path / "final_geometry.xyz"
            final_pdb = out_dir_path / "final_geometry.pdb"

            # Get layer indices for B-factor annotation
            # For heavy mode, base_calc is available; for light mode, create temporary calc
            layer_indices = None
            if use_heavy and 'base_calc' in dir():
                calc_core = base_calc.core if hasattr(base_calc, 'core') else base_calc
                layer_indices = {
                    "ml": getattr(calc_core, 'ml_indices', None),
                    "hess_mm": getattr(calc_core, 'hess_mm_indices', None),
                    "movable_mm": getattr(calc_core, 'movable_mm_indices', None),
                    "frozen": getattr(calc_core, 'frozen_layer_indices', None),
                }
            else:
                # For light mode, create a temporary calculator to get layer indices
                try:
                    temp_calc = mlmm(**calc_cfg)
                    calc_core = temp_calc.core if hasattr(temp_calc, 'core') else temp_calc
                    layer_indices = {
                        "ml": getattr(calc_core, 'ml_indices', None),
                        "hess_mm": getattr(calc_core, 'hess_mm_indices', None),
                        "movable_mm": getattr(calc_core, 'movable_mm_indices', None),
                        "frozen": getattr(calc_core, 'frozen_layer_indices', None),
                    }
                    del temp_calc
                except Exception:
                    layer_indices = None

            try:
                convert_xyz_to_pdb(final_xyz, ref_pdb, final_pdb)
                click.echo(f"[convert] Wrote '{final_pdb}'.")

                # Annotate B-factors with layer-based encoding
                if layer_indices and layer_indices.get("ml") is not None:
                    update_pdb_bfactors_from_layers(
                        final_pdb,
                        ml_indices=layer_indices["ml"] or [],
                        hess_mm_indices=layer_indices.get("hess_mm"),
                        movable_mm_indices=layer_indices.get("movable_mm"),
                        frozen_indices=layer_indices.get("frozen"),
                    )
                    click.echo(
                        f"[annotate]   B-factors set in '{final_pdb}' "
                        f"(ML={BFACTOR_ML:.0f}, MovableMM={BFACTOR_MOVABLE_MM:.0f}, "
                        f"FrozenMM={BFACTOR_FROZEN:.0f})."
                    )
            except Exception as e:
                click.echo(f"[convert] WARNING: Failed to convert final geometry to PDB: {e}", err=True)

            all_trj = out_dir_path / "optimization_all_trj.xyz"
            if all_trj.exists():
                try:
                    opt_pdb = out_dir_path / "optimization_all.pdb"
                    convert_xyz_to_pdb(all_trj, ref_pdb, opt_pdb)
                    click.echo(f"[convert] Wrote '{opt_pdb}'.")

                    # Annotate B-factors with layer-based encoding
                    if layer_indices and layer_indices.get("ml") is not None:
                        update_pdb_bfactors_from_layers(
                            opt_pdb,
                            ml_indices=layer_indices["ml"] or [],
                            hess_mm_indices=layer_indices.get("hess_mm"),
                            movable_mm_indices=layer_indices.get("movable_mm"),
                            frozen_indices=layer_indices.get("frozen"),
                        )
                        click.echo(
                            f"[annotate]   B-factors set in '{opt_pdb}' "
                            f"(ML={BFACTOR_ML:.0f}, MovableMM={BFACTOR_MOVABLE_MM:.0f}, "
                            f"FrozenMM={BFACTOR_FROZEN:.0f})."
                        )
                except Exception as e:
                    click.echo(f"[convert] WARNING: Failed to convert optimization trajectory to PDB: {e}", err=True)
        else:
            final_xyz = out_dir_path / "final_geometry.xyz"

        if out_json:
            from mlmm.core.utils import calculator_provenance, write_result_json
            _tsopt_imag_freqs: Optional[list] = []
            _tsopt_n_imag: Optional[int] = None
            _tsopt_energy = None
            _tsopt_status = "unverified"
            _tsopt_saddle_validation = "unavailable"
            _tsopt_hessian_status = "unavailable"
            _tsopt_hessian_error = None
            _tsopt_reaction_mode_index = None
            _tsopt_reaction_mode_frequency = None
            _tsopt_reaction_mode_overlap = None

            if use_heavy:
                # Heavy mode: use captured data from before del
                _tsopt_status = _heavy_status
                _tsopt_imag_freqs = _heavy_imag_freqs
                _tsopt_n_imag = _heavy_n_imag
                _tsopt_energy = _heavy_energy
                _tsopt_saddle_validation = _heavy_saddle_validation
                _tsopt_hessian_status = (
                    "skipped"
                    if skip_final_freq and not _heavy_stalled
                    else "completed" if _heavy_n_imag is not None
                    else "failed" if hessian_error
                    else "unavailable"
                )
                _tsopt_hessian_error = hessian_error
                _tsopt_reaction_mode_index = _heavy_reaction_mode_index
                _tsopt_reaction_mode_frequency = _heavy_reaction_mode_frequency
                _tsopt_reaction_mode_overlap = _heavy_reaction_mode_overlap
                _tsopt_n_atoms = len(geometry.atomic_numbers) if 'geometry' in dir() and geometry is not None else None
                _tsopt_n_opt_cycles = (
                    _heavy_cycle_ledger.spent
                    if "_heavy_cycle_ledger" in dir()
                    else None
                )
            else:
                # Light mode: compute freq/energy from runner
                _light_optimizer_converged = bool(
                    'runner' in dir()
                    and hasattr(runner, 'is_converged')
                    and runner.is_converged
                )
                # a dimer runner stall (energy-plateau child) wins over
                # every convergence/saddle-order verdict — even under
                # --skip-final-freq — and is never a converged saddle.
                _light_stalled = bool(
                    'runner' in dir() and getattr(runner, "is_stalled", False)
                )
                _tsopt_n_atoms = len(runner.geom.atomic_numbers) if 'runner' in dir() and hasattr(runner, 'geom') else None
                _tsopt_n_opt_cycles = runner._cycles_spent if 'runner' in dir() and hasattr(runner, '_cycles_spent') else None
                _tsopt_status = (
                    "stalled"
                    if _light_stalled
                    else "converged" if _light_optimizer_converged
                    else "not_converged"
                )
                if not skip_final_freq and 'runner' in dir():
                    _tsopt_n_imag = getattr(runner, "n_imaginary_modes", None)
                    _tsopt_imag_freqs = list(
                        getattr(runner, "imaginary_frequencies_cm", [])
                    )
                _tsopt_saddle_validation = _saddle_validation_from_count(
                    _tsopt_n_imag
                )
                _tsopt_hessian_status = (
                    getattr(runner, "hessian_status", "unavailable")
                    if 'runner' in dir()
                    else "unavailable"
                )
                _tsopt_hessian_error = (
                    getattr(runner, "hessian_error", None)
                    if 'runner' in dir()
                    else None
                )
                if _tsopt_n_imag and _tsopt_n_imag > 0:
                    _tsopt_reaction_mode_index = 0
                    _tsopt_reaction_mode_frequency = float(
                        _tsopt_imag_freqs[0]
                    ) if _tsopt_imag_freqs else None
                if 'runner' in dir() and hasattr(runner, 'geom'):
                    try:
                        _tsopt_energy = _calc_energy(runner.geom, calc_cfg)
                    except (RuntimeError, ValueError, AttributeError, KeyError) as _ee2:
                        logger.warning(
                            "Light-mode energy fallback failed: %s. "
                            "Reporting status='energy_missing' with NaN energy.",
                            _ee2,
                        )
                        _tsopt_energy = float("nan")

            result_data = {
                "status": _tsopt_status,
                "optimization_status": _tsopt_status,
                "saddle_validation": _tsopt_saddle_validation,
                "saddle_order_verified": _tsopt_saddle_validation == "first_order",
                "hessian_status": _tsopt_hessian_status,
                "hessian_error": _tsopt_hessian_error,
                "reaction_mode_index": _tsopt_reaction_mode_index,
                "reaction_mode_frequency_cm": _tsopt_reaction_mode_frequency,
                "reaction_mode_overlap": _tsopt_reaction_mode_overlap,
                "reaction_mode_source": (
                    "mep-reference-overlap"
                    if _tsopt_reaction_mode_overlap is not None
                    else "lowest-imaginary" if _tsopt_reaction_mode_index is not None
                    else None
                ),
                "flatten_requested": bool(simple_cfg.get("flatten_max_iter", 0)),
                "flatten_enabled": bool(simple_cfg.get("flatten_max_iter", 0)),
                "flatten_skip_reason": _flatten_skip_reason,
                "energy_hartree": _tsopt_energy,
                "n_imaginary_modes": _tsopt_n_imag,
                "frequency_zero_cutoff_cm": float(
                    frequency_cfg["zero_cutoff_cm"]
                ),
                "imaginary_frequencies_cm": _tsopt_imag_freqs,
                "opt_mode": opt_mode,
                "n_atoms": _tsopt_n_atoms,
                "n_opt_cycles": _tsopt_n_opt_cycles,
                **calculator_provenance(calc_cfg),
                "charge": calc_cfg.get("model_charge"),
                "spin": calc_cfg.get("model_mult"),
                "n_freeze_atoms": len(geom_cfg.get("freeze_atoms", [])),
                "thresh": (
                    rsirfo_cfg.get("thresh", simple_cfg.get("thresh"))
                    if use_heavy
                    else simple_cfg.get("thresh")
                ),
                "max_cycles": opt_cfg.get("max_cycles"),
                "input_file": str(input_path),
                "reference_mode_file": (
                    None if reference_mode_path is None else str(reference_mode_path)
                ),
                "reference_mode_candidate_count": len(reference_modes),
                "reference_mode_candidate_labels": list(reference_mode_labels),
                "reference_mode_cache": dict(reference_mode_metadata),
                "files": {"final_geometry_xyz": "final_geometry.xyz"},
                "rigid_projection": dict(
                    rigid_projection_info
                    if use_heavy
                    else getattr(runner, "rigid_projection_info", {})
                ),
            }
            if use_heavy:
                result_data["safeguards"] = _heavy_safeguards
                # additive microiteration serialization (present only when the
                # microiteration path ran). Legacy keys are unchanged.
                if _heavy_microiteration_obj is not None:
                    if _heavy_micro_cycles is not None:
                        result_data["n_micro_cycles"] = int(_heavy_micro_cycles)
                    result_data["microiteration"] = _heavy_microiteration_obj
            # Additive stop_reason, present only for a non-converged stop
            # (stalled/stopped) so a converged TS run's JSON stays
            # byte-compatible.
            if use_heavy:
                _tsopt_stop_reason = (
                    getattr(last_optimizer, "stop_reason", "") or ""
                    if 'last_optimizer' in dir()
                    else ""
                )
            else:
                _tsopt_stop_reason = (
                    getattr(runner, "stop_reason", "") or ""
                    if 'runner' in dir()
                    else ""
                )
            if _tsopt_stop_reason:
                result_data["stop_reason"] = _tsopt_stop_reason
            for ext in (".pdb", ".gjf"):
                f = out_dir_path / f"final_geometry{ext}"
                if f.exists():
                    result_data["files"][f"final_geometry_{ext[1:]}"] = f.name
            # Add trajectory files if they exist
            for _trj_name in ("optimization_all_trj.xyz", "optimization_all.pdb", "optimization_trj.xyz", "optimization.pdb"):
                _tf = out_dir_path / _trj_name
                if _tf.exists():
                    _key = _trj_name.replace(".", "_").replace("-", "_")
                    result_data["files"][_key] = _trj_name
            # List imaginary mode vib files
            _vib_dir = out_dir_path / "vib"
            if _vib_dir.exists():
                result_data["files"]["imaginary_mode_files"] = sorted([
                    f"vib/{f.name}"
                    for pattern in ("imag_*.pdb", "imag_*_trj.xyz")
                    for f in _vib_dir.glob(pattern)
                    if f.is_file()
                ])
            write_result_json(
                out_dir_path, result_data,
                command="tsopt",
                elapsed_seconds=time.perf_counter() - time_start,
            )

        # summary.md and key_* outputs are disabled.
        emit(
            format_elapsed("[time] Elapsed Time for TS Opt", time_start),
            narrative=True,
        )

    except ZeroStepLength as e:
        _write_error_json(out_dir_path, "tsopt", e, "ZeroStepLength", time_start)
        click.echo("ERROR: Proposed step length dropped below the minimum allowed (ZeroStepLength).", err=True)
        sys.exit(2)
    except OptimizationError as e:
        _write_error_json(out_dir_path, "tsopt", e, "OptimizationError", time_start)
        click.echo(f"ERROR: Optimization failed — {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except _TSOPTOutputCollisionError:
        raise
    except Exception as e:
        render_cli_exception(e, label="TS optimization", out_dir=out_dir_path, command="tsopt", time_start=time_start)
    finally:
        prepared_input.cleanup()
        # Release GPU memory (model + Hessian) so subsequent stages don't OOM.
        # ``geometry`` is a closure cell for the post-optimization helpers, so
        # clear it without deleting the name. Deleting that cell makes the
        # helpers' otherwise valid references statically and dynamically unsafe.
        base_calc = geometry = optimizer = last_optimizer = None
        macro_calc = macro_optimizer = mm_calc = None
        del base_calc, optimizer, last_optimizer
        del macro_calc, macro_optimizer, mm_calc
        gc.collect()  # break cyclic refs inside torch.nn.Module
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


# Allow `python -m mlmm.tsopt` direct execution
if __name__ == "__main__":
    cli()
