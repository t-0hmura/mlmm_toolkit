"""ML/MM geometry optimization with L-BFGS or RFO."""

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import contextlib
import gc
import io
from itertools import count
import logging

import sys

logger = logging.getLogger(__name__)

import click
from mlmm.core.output import emit
import numpy as np
import torch
import time

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.intcoords.exceptions import RebuiltInternalsException
from pysisyphus.constants import ANG2BOHR, BOHR2ANG, AU2EV
from pysisyphus.tr_projection import normalize_tr_projection_mode
from mlmm.workflows.restraints import HarmonicBiasCalculator
from pysisyphus.TablePrinter import TablePrinter

from mlmm.backends.mlmm_calc import mlmm, mlmm_mm_only
from mlmm.core.defaults import (
    BIAS_KW,
    HESSIAN_DIMER_KW,
    FREQ_KW,
    OPT_BASE_KW,
    LBFGS_KW,
    RFO_KW,
    THRESH_CHOICES,
    OPT_MODE_ALIASES,
    MICROITER_KW,
    OUT_DIR_OPT,
    BFACTOR_ML,
    BFACTOR_MOVABLE_MM,
    BFACTOR_FROZEN,
)
from mlmm.workflows.charge_prep import resolve_charge_spin_or_raise
from mlmm.core.utils import (
    append_xyz_trajectory as _append_xyz_trajectory,
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
    is_scan_spec_file,
    parse_dist_freeze_list,
    parse_dist_freeze_spec,
    load_pdb_atom_metadata,
    echo_resolved_device,
    emit_optimizer_terminal_status,
    finalize_microiter_macro_convergence,
    optimizer_cycle_count,
    pdb_keys_from_line as _pdb_keys_from_line,
    collect_ml_atom_keys as _collect_ml_atom_keys,
    format_pdb_with_bfactor as _format_with_bfactor,
    unbiased_energy_hartree,
    optional_positive_int,
)
from pysisyphus.normal_modes import (
    normalize_frequency_zero_cutoff_cm,
    resolved_imaginary_mask,
)
from mlmm.cli.common_options import (
    add_ml_charge_spin_options,
    add_ml_layer_detection_options,
    add_coord_type_option,
    add_print_every_option,
    add_precision_option, add_backend_model_option, add_calc_file_option,
    add_deterministic_option, add_allow_charge_mult_mismatch_option,
    add_workers_options,
)
from mlmm.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, make_is_param_explicit, _write_error_json, render_cli_exception
from mlmm.workflows._microiteration import (
    MicroiterationOutcome,
    MicroiterationPartition,
    OptimizerOutcome,
    PartitionError,
    build_aggregate,
    resolve_partition_from_core,
    micro_reached_force_equilibrium,
    describe_micro_stop,
    macro_progress_due,
)

EV2AU = 1.0 / AU2EV                 # eV → Hartree


class _OptOutputCollisionError(click.UsageError):
    """An OPT-owned destination aliases a consumed input."""


def _paths_physically_alias(path_a: Path, path_b: Path) -> bool:
    """Return whether two paths identify the same filesystem object."""

    try:
        if path_a.exists() and path_b.exists() and path_a.samefile(path_b):
            return True
    except OSError:
        pass
    return path_a.expanduser().resolve() == path_b.expanduser().resolve()


def _reject_opt_output_collisions(
    out_dir: Path,
    protected_inputs: Sequence[Optional[Path]],
) -> None:
    """Reject fixed OPT outputs that would overwrite a consumed input."""

    destinations = (
        out_dir / "final_geometry.xyz",
        out_dir / "final_geometry.pdb",
    )
    for destination in destinations:
        for protected in protected_inputs:
            if protected is None:
                continue
            if _paths_physically_alias(destination, Path(protected)):
                raise _OptOutputCollisionError(
                    f"Output '{destination}' aliases consumed input '{protected}'."
                )


def _invalidate_opt_optional_outputs(out_dir: Path) -> None:
    """Remove optional artifacts owned by an earlier OPT generation."""

    for name in (
        "final_geometry.pdb",
        "final_geometry.gjf",
        "optimization_all_trj.xyz",
        "optimization_trj.xyz",
        "optimization_all.pdb",
        "optimization.pdb",
        "result.json",
        "summary.json",
    ):
        candidate = out_dir / name
        if candidate.is_file() or candidate.is_symlink():
            candidate.unlink()
H_EVAA_2_AU = EV2AU / (ANG2BOHR * ANG2BOHR)  # (eV/Å^2) → (Hartree/Bohr^2)

# Flatten-loop constants (sourced from defaults.py)
OPT_FLATTEN_AMP_ANG = HESSIAN_DIMER_KW["flatten_amp_ang"]
OPT_FLATTEN_MAX_ITER = HESSIAN_DIMER_KW["flatten_max_iter"]
# Guard: a structure near a stationary point has at most a handful of
# spurious imaginary modes. Many (e.g. 100+ on a debug-capped, unconverged
# input) means the flatten loop is being applied off its design regime and
# would iterate per-mode (displace + re-opt) effectively forever. Skip with
# a clear "pre-optimize first" warning instead.
OPT_FLATTEN_UNCONVERGED_GUARD = 25


# Default settings + shared layer helpers now live in the neutral
# ``_opt_freq_common`` module so ``freq`` does not import ``opt``.
# Re-exported here for opt's own use and for the many workflows that import
# GEOM_KW / CALC_KW / _normalize_geom_freeze / _convert_yaml_layer_atoms_1to0
# from ``mlmm.workflows.opt``.
from mlmm.workflows._opt_freq_common import (  # noqa: F401
    GEOM_KW,
    CALC_KW,
    _normalize_geom_freeze,
    _convert_yaml_layer_atoms_1to0,
)

# Note: OPT_BASE_KW, LBFGS_KW, RFO_KW are imported from defaults.py

# Canonical home moved to mlmm.core.utils so cross-subcommand callers (e.g. sp.py)
# can import the same parser without depending on workflows/opt.py.
from mlmm.core.utils import _parse_freeze_atoms  # re-export for backward compat


def _parse_dist_freeze_args(
    raw_args: Sequence[str],
    one_based: bool,
    atom_meta: Optional[Sequence[Dict[str, Any]]],
) -> List[Tuple[int, int, Optional[float]]]:
    """Parse all ``--dist-freeze`` arguments (inline literal or spec file).

    Accepts the same format as ``--scan-lists``: inline Python literal
    (e.g. ``'[(1,5,1.4)]'``) or a YAML/JSON spec file path.  String atom
    specs (e.g. ``'A:SER123:OG'``) are supported when *atom_meta* is
    available.  Target distance is optional — omit to freeze at the current
    distance.
    """
    all_pairs: List[Tuple[int, int, Optional[float]]] = []
    for raw in raw_args:
        if is_scan_spec_file(raw):
            all_pairs.extend(parse_dist_freeze_spec(
                Path(raw),
                one_based_default=one_based,
                atom_meta=atom_meta,
            ))
        else:
            all_pairs.extend(parse_dist_freeze_list(
                raw,
                one_based=one_based,
                atom_meta=atom_meta,
            ))
    return all_pairs


def _resolve_dist_freeze_targets(
    geometry,
    tuples: List[Tuple[int, int, Optional[float]]],
) -> List[Tuple[int, int, float]]:
    coords_bohr = np.array(geometry.coords3d, dtype=float).reshape(-1, 3)
    coords_ang = coords_bohr * BOHR2ANG
    n = coords_ang.shape[0]
    resolved: List[Tuple[int, int, float]] = []
    for (i, j, target) in tuples:
        if not (0 <= i < n and 0 <= j < n):
            raise click.BadParameter(
                f"--dist-freeze indices {(i, j)} are out of bounds for the loaded geometry (N={n})."
            )
        if target is None:
            vec = coords_ang[i] - coords_ang[j]
            dist = float(np.linalg.norm(vec))
        else:
            dist = float(target)
        resolved.append((i, j, dist))
    return resolved



def _annotate_b_factors_inplace(
    pdb_path: Path,
    model_pdb: Path,
    freeze_indices_0based: Sequence[int],
    beta_ml: float = 100.0,
    beta_frz: float = 50.0,
    beta_both: float = 150.0,
) -> None:
    """
    Overwrite B-factors in-place:
      - ML-region atoms: 100.00
      - frozen atoms: 50.00
      - ML ∩ frozen: 150.00
    Indexing for 'frozen' is 0-based and resets at each MODEL.
    """
    ml_full, ml_simple = _collect_ml_atom_keys(model_pdb)
    frozen_set = set(int(i) for i in (freeze_indices_0based or []))

    try:
        lines = pdb_path.read_text().splitlines(keepends=True)
    except Exception:
        logger.debug("Failed to read PDB file for B-factor annotation: %s", pdb_path, exc_info=True)
        return

    out_lines: List[str] = []
    atom_idx = 0  # resets per MODEL

    for line in lines:
        rec = line[:6]
        if rec.startswith("MODEL"):
            # reset atom counter for each model
            atom_idx = 0
            out_lines.append(line)
            continue
        if rec.startswith("ATOM  ") or rec.startswith("HETATM"):
            kf, ks = _pdb_keys_from_line(line)
            is_ml = (kf in ml_full) or (ks in ml_simple)
            is_frz = (atom_idx in frozen_set)
            if is_ml and is_frz:
                out_lines.append(_format_with_bfactor(line, beta_both))
            elif is_ml:
                out_lines.append(_format_with_bfactor(line, beta_ml))
            elif is_frz:
                out_lines.append(_format_with_bfactor(line, beta_frz))
            else:
                out_lines.append(line)
            atom_idx += 1
        else:
            out_lines.append(line)

    try:
        pdb_path.write_text("".join(out_lines))
    except Exception:
        logger.debug("Failed to write B-factor annotated PDB: %s", pdb_path, exc_info=True)


def _maybe_convert_outputs_to_pdb(
    input_path: Path,
    out_dir: Path,
    dump: bool,
    get_trj_fn,
    final_xyz_path: Path,
    model_pdb: Path,
    freeze_indices_0based: Sequence[int],
    ml_indices: Optional[List[int]] = None,
    hess_mm_indices: Optional[List[int]] = None,
    movable_mm_indices: Optional[List[int]] = None,
    frozen_layer_indices: Optional[List[int]] = None,
) -> None:
    """
    If the input is a PDB, convert outputs (final_geometry.xyz and, if dump, optimization_all_trj.xyz /
    optimization_trj.xyz) to PDB,
    and annotate B-factors for the 3-layer ML/MM system.

    B-factor encoding (3-layer system):
        ML atoms: 0.0
        Movable MM atoms: 10.0
        Frozen MM atoms: 20.0

    If layer indices are not provided, falls back to legacy encoding:
        ML atoms: 100.0
        Frozen atoms: 50.0
        ML ∩ frozen: 150.0
    """
    if not is_convert_file_enabled():
        return
    if input_path.suffix.lower() != ".pdb":
        return

    # Determine if we should use the layer-based B-factor encoding
    use_layer_bfactors = ml_indices is not None

    ref_pdb = input_path.resolve()
    # final_geometry.xyz → final_geometry.pdb
    final_pdb = out_dir / "final_geometry.pdb"
    try:
        convert_xyz_to_pdb(final_xyz_path, ref_pdb, final_pdb)
        click.echo(f"[convert] Wrote '{final_pdb}'.")

        if use_layer_bfactors:
            update_pdb_bfactors_from_layers(
                final_pdb,
                ml_indices=ml_indices or [],
                hess_mm_indices=hess_mm_indices,
                movable_mm_indices=movable_mm_indices,
                frozen_indices=frozen_layer_indices,
            )
            click.echo(
                f"[annotate]   B-factors set in '{final_pdb}' "
                f"(ML={BFACTOR_ML:.0f}, MovableMM={BFACTOR_MOVABLE_MM:.0f}, "
                f"FrozenMM={BFACTOR_FROZEN:.0f})."
            )
        else:
            # Fall back to legacy encoding
            _annotate_b_factors_inplace(
                final_pdb,
                model_pdb=model_pdb,
                freeze_indices_0based=freeze_indices_0based,
            )
            click.echo(f"[annotate]   B-factors set in '{final_pdb}' (ML=100, frozen=50, both=150).")
    except Exception as e:
        click.echo(f"[convert] WARNING: Failed to convert final geometry to PDB: {e}", err=True)

    # optimization_all_trj.xyz / optimization_trj.xyz → PDB (if dump)
    if dump:
        try:
            wrote_any = False
            all_trj_path = get_trj_fn("optimization_all_trj.xyz")
            if all_trj_path.exists():
                all_opt_pdb = out_dir / "optimization_all.pdb"
                convert_xyz_to_pdb(all_trj_path, ref_pdb, all_opt_pdb)
                click.echo(f"[convert] Wrote '{all_opt_pdb}'.")
                wrote_any = True

                if use_layer_bfactors:
                    update_pdb_bfactors_from_layers(
                        all_opt_pdb,
                        ml_indices=ml_indices or [],
                        hess_mm_indices=hess_mm_indices,
                        movable_mm_indices=movable_mm_indices,
                        frozen_indices=frozen_layer_indices,
                    )
                    click.echo(
                        f"[annotate]   B-factors set in '{all_opt_pdb}' "
                        f"(ML={BFACTOR_ML:.0f}, MovableMM={BFACTOR_MOVABLE_MM:.0f}, "
                        f"FrozenMM={BFACTOR_FROZEN:.0f})."
                    )
                else:
                    _annotate_b_factors_inplace(
                        all_opt_pdb,
                        model_pdb=model_pdb,
                        freeze_indices_0based=freeze_indices_0based,
                    )
                    click.echo(f"[annotate]   B-factors set in '{all_opt_pdb}' (ML=100, frozen=50, both=150).")

            trj_path = get_trj_fn("optimization_trj.xyz")
            if trj_path.exists():
                opt_pdb = out_dir / "optimization.pdb"
                convert_xyz_to_pdb(trj_path, ref_pdb, opt_pdb)
                click.echo(f"[convert] Wrote '{opt_pdb}'.")
                wrote_any = True

                if use_layer_bfactors:
                    update_pdb_bfactors_from_layers(
                        opt_pdb,
                        ml_indices=ml_indices or [],
                        hess_mm_indices=hess_mm_indices,
                        movable_mm_indices=movable_mm_indices,
                        frozen_indices=frozen_layer_indices,
                    )
                    click.echo(
                        f"[annotate]   B-factors set in '{opt_pdb}' "
                        f"(ML={BFACTOR_ML:.0f}, MovableMM={BFACTOR_MOVABLE_MM:.0f}, "
                        f"FrozenMM={BFACTOR_FROZEN:.0f})."
                    )
                else:
                    _annotate_b_factors_inplace(
                        opt_pdb,
                        model_pdb=model_pdb,
                        freeze_indices_0based=freeze_indices_0based,
                    )
                    click.echo(f"[annotate]   B-factors set in '{opt_pdb}' (ML=100, frozen=50, both=150).")

            if not wrote_any:
                click.echo(
                    "[convert] WARNING: neither 'optimization_all_trj.xyz' nor 'optimization_trj.xyz' was found; "
                    "skipping trajectory PDB conversion.",
                    err=True,
                )
        except Exception as e:
            click.echo(f"[convert] WARNING: Failed to convert optimization trajectory to PDB: {e}", err=True)




# Called in flatten-loop tight succession; per-call empty_cache
# prevents unbounded VRAM growth across iterations. The cache-clear is the
# load-bearing part of this helper. Lives in mlmm.core.calc_eval so the
# tsopt module can reuse the same implementation without re-duplicating it.
from mlmm.core.calc_eval import calc_energy as _calc_energy  # noqa: E402


def _set_cartesian_flatten_coords(geom, cart_coords: np.ndarray) -> None:
    """Install a Cartesian trial while accepting an internal-basis rebuild."""

    try:
        geom.cart_coords = np.asarray(cart_coords, dtype=float).reshape(-1)
    except RebuiltInternalsException:
        # Geometry has already installed the Cartesian coordinates and rebuilt
        # its primitive set before signalling this control-flow exception.
        # Flatten probes evaluate directly through the active calculator, so a
        # state clear is sufficient; no optimizer reset is involved here.
        geom.clear()


def _flatten_all_imag_modes_for_geom(
    geom,
    masses_amu: np.ndarray,
    calc_kwargs: dict,
    freqs_cm: np.ndarray,
    modes: torch.Tensor,
    neg_freq_thresh_cm: float,
    flatten_amp_ang: float,
    calculator=None,
) -> bool:
    """
    Flatten all imaginary modes for a geometry in a single pass.
    """
    neg_idx_all = np.where(freqs_cm < -abs(neg_freq_thresh_cm))[0]
    if len(neg_idx_all) == 0:
        return False

    if len(neg_idx_all) > OPT_FLATTEN_UNCONVERGED_GUARD:
        click.echo(
            f"[Flatten] WARNING: {len(neg_idx_all)} imaginary modes "
            f"(> {OPT_FLATTEN_UNCONVERGED_GUARD}; below "
            f"{-abs(neg_freq_thresh_cm):.1f} cm^-1) — the structure is far "
            f"from a stationary point, so the flatten loop is skipped (it "
            f"would displace + re-optimize per mode, effectively forever). "
            f"Pre-optimize first (e.g. `mlmm opt --thresh baker`) and rerun "
            f"with --flatten.",
            err=True,
        )
        return False

    order = np.argsort(freqs_cm[neg_idx_all])  # most negative first
    targets = [int(x) for x in neg_idx_all[order]]
    amp_bohr = float(flatten_amp_ang) / BOHR2ANG
    E_ref = _calc_energy(geom, calc_kwargs, calc=calculator)

    m3 = np.repeat(masses_amu, 3).reshape(-1, 3)
    for idx in targets:
        v_mw = modes[idx].detach().cpu().numpy().reshape(-1, 3)
        # A returned mode row is a mass-weighted eigenvector q of
        # M^(-1/2) H M^(-1/2).  The Cartesian normal-mode direction is
        # u = M^(-1/2) q / ||M^(-1/2) q||, computed once here (divide by
        # sqrt(m) then L2-normalize).  The flatten displacement is
        # amp_bohr * u so ||disp|| == amp_bohr.  A second per-atom mass
        # factor would rotate the direction toward M^(-1) q and change the
        # amplitude; there is no such factor.
        v_cart = v_mw / np.sqrt(m3)
        v_cart /= np.linalg.norm(v_cart)

        disp = amp_bohr * v_cart
        ref = geom.cart_coords.reshape(-1, 3)

        plus = ref + disp
        minus = ref - disp

        _set_cartesian_flatten_coords(geom, plus)
        E_plus = _calc_energy(geom, calc_kwargs, calc=calculator)

        _set_cartesian_flatten_coords(geom, minus)
        E_minus = _calc_energy(geom, calc_kwargs, calc=calculator)

        use_plus = E_plus <= E_minus
        _set_cartesian_flatten_coords(geom, plus if use_plus else minus)
        E_keep = E_plus if use_plus else E_minus
        delta_e = E_keep - E_ref
        click.echo(
            f"[Flatten] mode={idx} freq={freqs_cm[idx]:+.2f} cm^-1 "
            f"E_disp={E_keep:.8f} Ha \u0394E={delta_e:+.8f} Ha"
        )

    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return True


def _seed_rfo_initial_hessian(
    geometry,
    calc_cfg: Dict[str, Any],
    calculator,
    *,
    restraints_active: bool,
) -> str:
    """Seed RFO from the exact active PES, with safe cache reuse."""

    from mlmm.workflows.freq import (
        _calc_full_hessian_torch as _freq_calc_full_hessian_torch,
        _torch_device as _freq_torch_device,
    )

    hess_device = _freq_torch_device(calc_cfg.get("ml_device", "auto"))
    if restraints_active:
        click.echo(
            "[opt] Distance restraints are active; calculating "
            "the initial RFO Hessian on the restrained PES."
        )
        h_init, _ = _freq_calc_full_hessian_torch(
            geometry,
            calc_cfg,
            hess_device,
            refresh_geom_meta=True,
            calculator=calculator,
        )
        geometry.cart_hessian = h_init
        click.echo(
            f"[opt] Initial restrained Hessian seeded "
            f"(shape={h_init.shape[0]}x{h_init.shape[1]})."
        )
        return "restrained"

    from mlmm.io.hessian_cache import (
        load_matching as _hess_load_matching,
        identity_from_context as _hess_identity,
    )

    # reuse an IRC endpoint Hessian only on a full evaluation-identity
    # match (run/system/evaluator/active space/potential).
    cached = _hess_load_matching(
        "irc_endpoint",
        _hess_identity(geometry, calc_cfg, role="irc_endpoint"),
    )
    if cached is not None:
        click.echo("[opt] Reusing IRC endpoint Hessian for RFO seeding.")
        active_dofs = cached.get("active_dofs")
        h_raw = cached["hessian"]
        if isinstance(h_raw, torch.Tensor):
            h_init = h_raw.clone()
        else:
            h_init = torch.as_tensor(h_raw, dtype=torch.float64)
        if active_dofs is not None:
            geometry.within_partial_hessian = {
                "active_n_dof": len(active_dofs),
                "full_n_dof": geometry.cart_coords.size,
                "active_dofs": active_dofs,
                "active_atoms": sorted(set(d // 3 for d in active_dofs)),
            }
        geometry.cart_hessian = h_init
        click.echo(
            f"[opt] Initial Hessian seeded "
            f"(shape={h_init.shape[0]}x{h_init.shape[1]})."
        )
        return "irc_cache"

    click.echo("[opt] Seeding initial Hessian via shared freq backend.")
    h_init, _ = _freq_calc_full_hessian_torch(
        geometry,
        calc_cfg,
        hess_device,
        refresh_geom_meta=True,
        calculator=calculator,
    )
    geometry.cart_hessian = h_init
    click.echo(
        f"[opt] Initial Hessian seeded "
        f"(shape={h_init.shape[0]}x{h_init.shape[1]})."
    )
    return "fresh"


def _opt_terminal_converged(
    use_microiter: bool,
    microiter_result: Optional[Dict[str, Any]],
    optimizer: Any,
) -> Optional[bool]:
    """Source the opt terminal convergence flag from whichever runner produced it.

    On the QM/MM microiteration path there is no standalone ``optimizer`` in
    scope, so the convergence truth lives in ``microiter_result['converged']``
    (mirroring how ``is_stalled`` / ``stop_reason`` are already sourced).  A
    genuinely converged microiter run must therefore report ``converged`` in
    result.json, matching the ``[opt] Converged!`` console emit -- reading a
    missing ``optimizer`` here was the bug that mislabeled it ``not_converged``.
    Returns ``None`` (unknown) when neither runner is available.
    """
    if use_microiter:
        if microiter_result is not None:
            return bool(microiter_result.get("converged"))
        return None
    if optimizer is not None and hasattr(optimizer, "is_converged"):
        return bool(optimizer.is_converged)
    return None


def _run_microiter_opt(
    geometry,
    base_calc,
    calc_cfg: Dict[str, Any],
    rfo_cfg: Dict[str, Any],
    lbfgs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    microiter_cfg: Dict[str, Any],
    out_dir_path: Path,
    *,
    partition: MicroiterationPartition,
    dump: bool = False,
) -> Dict[str, Any]:
    """Run macro/micro alternating optimization (Gaussian 16-style microiteration).

    Macro step: 1 RFO step moving ML atoms + link-atom MM parents (full ONIOM force).
    Micro step: LBFGS relaxing remaining MM atoms with MM-only forces until convergence.
    Link-atom MM parents are included in the macro step to maintain consistency
    of the link atom position across macro/micro boundaries.

    The macro/micro phase masks and the user's original freeze mask come from the
    single immutable :class:`MicroiterationPartition` resolved by the caller from
    the accepted calculator core: the original freeze is preserved in
    both phases and restored exactly on every exit path.  A micro relaxation that
    does not explicitly converge fails closed and never reads as macro
    convergence.
    """
    # a valid empty-ML partition is a documented caller-side fallback to a
    # real standard optimization; reaching this driver with no macro-active atoms
    # is a contract violation, raised loudly (never an unchanged-geometry return).
    if not partition.has_macro_active:
        raise PartitionError(
            "microiteration reached the macro/micro driver with no ML macro-active "
            "atoms; the caller must dispatch a standard optimization instead."
        )

    # consume the single immutable partition. The user's original freeze
    # mask is preserved in BOTH phase masks and restored exactly in ``finally``.
    ml_indices = list(partition.ml_atoms)
    movable_mm = list(partition.movable_mm_atoms)
    link_mm_parents = set(partition.link_parent_atoms)
    original_freeze = list(partition.original_freeze)
    frozen_mm = list(original_freeze)

    n_atoms = partition.n_atoms

    macro_freeze = list(partition.macro_freeze_atoms)
    micro_freeze = list(partition.micro_freeze_atoms)

    max_cycles = optional_positive_int(opt_cfg.get("max_cycles"), "opt.max_cycles")
    thresh = opt_cfg.get("thresh", "gau")
    micro_thresh = microiter_cfg.get("micro_thresh") or thresh
    micro_max_cycles = optional_positive_int(
        microiter_cfg.get("micro_max_cycles"), "microiter.micro_max_cycles"
    )

    click.echo(
        f"[microiter] ML atoms: {len(ml_indices)}, "
        f"Link MM parents: {len(link_mm_parents)}, "
        f"Movable MM atoms: {len(movable_mm)}, "
        f"Frozen MM atoms: {len(frozen_mm)}"
    )
    click.echo(f"[microiter] Macro thresh: {thresh}, Micro thresh: {micro_thresh}")

    # Reuse the caller's accepted ONIOM core for MM-only work.
    mm_calc = mlmm_mm_only(base_calc.core, freeze_atoms=micro_freeze)

    # Ordered record of every micro (MM) relaxation outcome.
    micro_attempts: List[OptimizerOutcome] = []
    micro_cycles_total = 0

    def _relax_micro() -> Tuple[Any, int]:
        """Run one MM-only micro relaxation on a cart twin; copy coords back.

        Returns the LBFGS optimizer (for its explicit convergence bit) and the
        number of executed micro cycles.  A normal Python return is not, by
        itself, evidence of convergence: the caller reads
        ``micro_opt.is_converged`` via :class:`OptimizerOutcome`.
        """

        macro_coord_type = getattr(geometry, "coord_type", "cart")
        if macro_coord_type != "cart":
            from pysisyphus.Geometry import Geometry as _Geometry
            micro_geom = _Geometry(
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
        # The MM equilibration never uses the plateau stop; see the twin in
        # `tsopt.py`. `micro_max_cycles` is the real bound.
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

        # Resolve the initial MM equilibrium before matching or computing the
        # macro Hessian.  The persistent RFO model must start from the same
        # coordinates and gradients that its first macro step sees.
        if partition.has_micro_active:
            _init_micro_opt, _init_micro_steps = _relax_micro()
            micro_cycles_total += _init_micro_steps
            _init_micro_out = OptimizerOutcome.from_optimizer(
                _init_micro_opt, max_cycles=micro_max_cycles
            )
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
                emit(
                    "[microiter] Initial MM equilibration did not converge "
                    f"(status={_init_micro_out.status}); no macro step is taken.",
                    narrative=True,
                )
            del _init_micro_opt
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

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
                "optimizer": None,
                "micro_cycles": micro_cycles_total,
                "outcome": micro_outcome,
            }

        # Seed initial Hessian for RFO (with macro freeze)
        # Try IRC endpoint cache first; fall back to full Hessian calculation.
        from mlmm.io.hessian_cache import (
            load_matching as _hess_load_matching,
            identity_from_context as _hess_identity,
            reconcile_active_hessian as _hess_reconcile_active,
        )
        from mlmm.workflows.freq import (
            _calc_full_hessian_torch as _freq_calc_full_hessian_torch,
            _torch_device as _freq_torch_device,
        )
        hess_device = _freq_torch_device(calc_cfg.get("ml_device", "auto"))

        # Always create macro calculator (needed for optimization loop below)
        macro_calc_cfg = dict(calc_cfg)
        macro_calc_cfg["freeze_atoms"] = macro_freeze
        macro_calc_cfg["hess_mm_atoms"] = sorted(link_mm_parents)  # ML + link MM parents in Hessian
        macro_calc = mlmm(
            **macro_calc_cfg,
            _high_level_backend=base_calc.core._ml_backend,
        )

        # reuse an IRC endpoint Hessian only on a full evaluation-identity
        # match (run/system/evaluator/active space/potential).
        cached = _hess_load_matching(
            "irc_endpoint",
            _hess_identity(geometry, calc_cfg, role="irc_endpoint"),
        )
        _cache_used = False
        macro_free_atoms = sorted(
            set(range(geometry.cart_coords.size // 3)) - set(macro_freeze)
        )
        macro_free_dofs = [
            3 * atom + axis
            for atom in macro_free_atoms
            for axis in range(3)
        ]
        if cached is not None:
            h_init = _hess_reconcile_active(
                cached,
                macro_free_dofs,
                full_n_dof=geometry.cart_coords.size,
            )
            if h_init is not None:
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
                    "[microiter] Reusing IRC endpoint Hessian for RFO macro "
                    f"step (shape={h_init.shape[0]}x{h_init.shape[1]})."
                )
                _cache_used = True
            else:
                click.echo(
                    "[microiter] IRC endpoint Hessian basis does not cover the "
                    "ordered macro DOFs. Falling back to a fresh Hessian."
                )
            if h_init is not None:
                del h_init
        if not _cache_used:
            click.echo("[microiter] Seeding initial Hessian for RFO macro step.")
            geometry.freeze_atoms = macro_freeze
            geometry.set_calculator(macro_calc)

            h_init, _ = _freq_calc_full_hessian_torch(
                geometry,
                macro_calc_cfg,
                hess_device,
                refresh_geom_meta=True,
                calculator=macro_calc,
            )
            geometry.cart_hessian = h_init
            click.echo(f"[microiter] Initial Hessian seeded (shape={h_init.shape[0]}x{h_init.shape[1]}).")
            del h_init

        # Create persistent RFOptimizer once (LayerOpt pattern).
        # This preserves the BFGS Hessian update chain across macro iterations.
        # NOTE: geometry already has macro_calc set (line above); do NOT call
        # set_calculator again as it clears the pre-computed cart_hessian.
        geometry.freeze_atoms = macro_freeze

        # The macro step IS this run's optimizer, so it honours the shared `opt`
        # block exactly as an ordinary (non-microiter) RFO run does -- same merge
        # rule, so the two paths cannot drift apart key by key. Only values the
        # user actually changed are passed, which keeps optimizer-specific
        # `rfo.*` settings authoritative for untouched defaults. The micro step
        # below deliberately takes none of this: it has its own `microiter.*`
        # threshold and cycle bound.
        rfo_args = {
            **rfo_cfg,
            **strip_inherited_keys(dict(opt_cfg), OPT_BASE_KW, mode="same"),
        }
        rfo_args["max_cycles"] = max_cycles
        rfo_args["out_dir"] = str(out_dir_path)
        rfo_args["dump"] = False  # trajectory dumping handled externally
        rfo_args["thresh"] = thresh

        macro_optimizer = RFOptimizer(geometry, **rfo_args)
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
            # ---- Macro step: 1 RFO step with ONIOM forces, MM frozen ----
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

            # A real macro energy-plateau stall  stops the loop BEFORE the
            # step is applied or a micro relaxation is launched; it is never
            # convergence.
            if macro_optimizer.stop_requested:
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
                print()  # blank line closes the table
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
            rebuilt_internals = False
            try:
                geometry.coords = new_coords
                # Record actual step (may differ due to coordinate back-transformation)
                macro_optimizer.steps[-1] = (
                    geometry.coords - macro_optimizer.coords[-1]
                )
            except RebuiltInternalsException:
                click.echo(
                    "[microiter] Internal coordinates were rebuilt; resetting "
                    "the macro optimizer after MM relaxation.",
                    err=True,
                )
                geometry.clear()
                rebuilt_internals = True

            # ---- Micro step: MM relaxation on a cart-only twin geometry ----
            # ``_relax_micro`` runs the MM relaxation on a cart twin and copies the
            # converged positions back via the coords3d setter; the macro chemistry
            # stays in DLC (mirrors tsopt's _run_microiter_tsopt).
            if partition.has_micro_active:
                micro_opt, micro_steps = _relax_micro()
                micro_cycles_total += micro_steps
                _micro_out = OptimizerOutcome.from_optimizer(micro_opt, max_cycles=micro_max_cycles)
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
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
            else:
                # A validated partition with macro-active atoms but
                # ZERO micro-active MM atoms (e.g. the entire MM region is
                # user-frozen) has no movable MM coordinate to relax. Append the
                # zero-cycle vacuous micro success instead of building an LBFGS with
                # every atom frozen (mirrors the initial-equilibration guard above).
                _micro_out = OptimizerOutcome.vacuous_success()
                micro_steps = 0
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
            if _micro_out.converged is not True:
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
            if rebuilt_internals:
                macro_optimizer = RFOptimizer(geometry, **rfo_args)
                macro_optimizer.prepare_opt()

        else:
            if run_macro:
                print()  # blank line closes the table (print() shares the table's stdout path)
                emit(f"[microiter] Reached max macro iterations ({max_cycles}).", detail=True)

        # terminal outcome. A stalled latest micro (MM) relaxation must not
        # masquerade as clean macro convergence, and it must not be lost as a
        # reasonless not_converged when the macro merely ran out of cycles: surface
        # it as a stall (with its reason) in either case.
        finalize_microiter_macro_convergence(
            macro_optimizer,
            macro_converged=macro_converged,
            latest_micro_stalled=latest_micro_stalled,
            latest_micro_stop_reason=latest_micro_stop_reason,
        )
        # Fold the macro state after stall demotion and the ordered micro
        # attempts into one fail-closed outcome.
        # into a single aggregate: converged only when the macro is converged AND the
        # latest required micro relaxation is converged.
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
            "optimizer": macro_optimizer,
            "micro_cycles": micro_cycles_total,
            "outcome": micro_outcome,
        }

        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        emit(f"[microiter] Total macro steps: {total_macro_steps}", detail=True)

        return outcome
    finally:
        # restore the EXACT original freeze mask + base calculator on
        # every exit path (success, micro non-convergence, macro stall, and
        # raised exception), so no user freeze constraint leaks past the
        # microiteration driver.
        geometry.freeze_atoms = list(original_freeze)
        geometry.set_calculator(base_calc)
        if torch.cuda.is_available():
            torch.cuda.empty_cache()



@click.command(
    help="ML/MM geometry optimization with L-BFGS or RFO.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (PDB/mmCIF, or XYZ). XYZ provides higher coordinate precision. "
         "If XYZ, use --ref-pdb to specify PDB topology for atom ordering and output conversion.",
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
    help="Amber parm7 topology covering the whole enzyme complex.",
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
    show_default="all movable MM atoms",
    help="Distance cutoff (Å) from ML region for MM atoms to include in Hessian calculation. "
         "Applied to movable MM atoms and can be combined with --detect-layer.",
)
@click.option(
    "--movable-cutoff",
    "movable_cutoff",
    type=float,
    default=None,
    show_default="use freeze_atoms",
    help="Distance cutoff (Å) from ML region for movable MM atoms. "
     "MM atoms beyond this are frozen. "
         "Providing --movable-cutoff disables --detect-layer and uses distance-based layer assignment.",
)
@click.option(
    "--dist-freeze",
    "dist_freeze_raw",
    type=str,
    multiple=True,
    default=(),
    show_default=False,
    help="Distance restraints: inline Python literal (e.g. '[(1,5,1.4)]') or a YAML/JSON spec file path. "
         "Format: (i,j,target_Å) triples. "
         "Target may be omitted to freeze at the current distance: (i,j).",
)
@click.option(
    "--one-based/--zero-based",
    "one_based",
    default=True,
    show_default=True,
    help="Interpret --dist-freeze indices as 1-based or 0-based.",
)
@click.option(
    "--bias-k",
    type=float,
    default=None,
    show_default="300.0",
    help=(
        "Harmonic restraint strength k [eV/Å^2] for --dist-freeze. "
        "YAML bias.k applies when this option is omitted; explicit CLI wins."
    ),
)
@click.option("--max-cycles", type=click.IntRange(min=1), default=None, show_default="100000", help="Maximum number of optimization cycles.")
@click.option(
    "--dump/--no-dump",
    default=False,
    show_default=True,
    help="Write optimization trajectories ('optimization_trj.xyz' and 'optimization_all_trj.xyz').",
)
@click.option("-o", "--out-dir", type=str, default=OUT_DIR_OPT, show_default=True, help="Output directory.")
@click.option(
    "--thresh",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default=None,
    show_default="gau",
    help="Convergence preset.",
)
@click.option(
    "--opt-mode",
    type=click.Choice(["grad", "hess", "lbfgs", "rfo"], case_sensitive=False),
    default="grad",
    show_default=True,
    help="Optimization mode: grad/lbfgs or hess/rfo.",
)
@click.option(
    "--microiter/--no-microiter",
    "microiter",
    default=True,
    show_default=True,
    help="Enable microiteration: alternate ML 1-step (RFO) and MM relaxation (L-BFGS with MM-only forces). "
         "Only effective in --opt-mode hess (RFO). Ignored in grad mode.",
)
@click.option(
    "--flatten/--no-flatten",
    "flatten",
    default=False,
    show_default=True,
    help="Enable/disable imaginary-mode flatten loop after optimization.",
)
@click.option(
    "--reject-uphill/--no-reject-uphill",
    "reject_uphill",
    default=False,
    show_default=True,
    help=(
        "Opt in to rejecting uphill RFO trials in hess mode (tolerance: "
        "1e-4 Hartree) and final-check the retained geometry at the emergency "
        "floor. Ignored in grad/lbfgs mode."
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
    help="Validate options and print the execution plan without running optimization.",
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
    "--mm-only/--no-mm-only",
    "mm_only",
    default=False,
    show_default=True,
    help="Skip the MLIP component entirely and minimize using only the MM "
         "force field on the full system. Layers (movable/frozen) are still "
         "honored via B-factor encoding or --movable-cutoff. Only "
         "--opt-mode grad (L-BFGS) is supported in this mode; microiteration "
         "is automatically disabled.",
)
@click.option(
    "--cmap/--no-cmap",
    "use_cmap",
    default=None,
    show_default="cmap",
    help="Preserve CMAP terms in both real and model MM layers when present in parm7.",
)
@click.option(
    "--out-json/--no-out-json",
    "out_json",
    default=False,
    show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
@add_ml_layer_detection_options()
@add_ml_charge_spin_options()
@add_coord_type_option()
@add_print_every_option()
@add_precision_option()
@add_workers_options()
@add_backend_model_option()
@add_calc_file_option()
@add_deterministic_option()
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
        "convergence; an explicit --max-cycles remains the hard bound. The MM micro "
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
    dist_freeze_raw: Sequence[str],
    one_based: bool,
    bias_k: Optional[float],
    max_cycles: int,
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    opt_mode: str,
    microiter: bool,
    flatten: bool,
    reject_uphill: bool,
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    convert_files: bool,
    backend: Optional[str],
    embedcharge: bool,
    embedcharge_cutoff: Optional[float],
    link_atom_method: Optional[str],
    mm_backend: Optional[str],
    mm_only: bool,
    use_cmap: Optional[bool],
    out_json: bool,
    cli_coord_type: Optional[str],
    print_every: Optional[int],
    precision: Optional[str],
    workers: Optional[int],
    workers_per_node: Optional[int],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: Optional[str],
) -> None:
    set_convert_file_enabled(convert_files)
    time_start = time.perf_counter()
    error_out_dir = Path(out_dir).resolve()
    prepared_input = None

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
            raise click.UsageError(
                "XYZ/TRJ input requires --ref-pdb to specify PDB topology."
            )
        prepared_input = prepare_input_structure(input_path)
        apply_ref_pdb_override(prepared_input, ref_pdb)
        click.echo(f"[input] Using XYZ coordinates from {input_path.name}, PDB topology from {ref_pdb.name}")
    else:
        click.echo(f"ERROR: Unsupported input format: {suffix}. Use .pdb/.cif/.mmcif or .xyz (with --ref-pdb).", err=True)
        sys.exit(1)

    geom_input_path = prepared_input.geom_path
    charge, spin = resolve_charge_spin_or_raise(
        prepared_input, charge, spin,
        ligand_charge=ligand_charge, prefix="[opt]",
        model_pdb=model_pdb,
        model_indices_spec=model_indices_str,
        detect_layer=detect_layer,
        yaml_cfg=merged_yaml_cfg,
    )

    try:
        freeze_atoms_cli = _parse_freeze_atoms(freeze_atoms_text)
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

    pdb_atom_meta: List[Dict[str, Any]] = []
    if prepared_input.source_path.suffix.lower() == ".pdb":
        pdb_atom_meta = load_pdb_atom_metadata(prepared_input.source_path)

    try:
        dist_freeze = _parse_dist_freeze_args(
            dist_freeze_raw, one_based=bool(one_based), atom_meta=pdb_atom_meta,
        )
    except click.BadParameter as e:
        click.echo(f"ERROR: {e}", err=True)
        prepared_input.cleanup()
        sys.exit(1)

    # Resolve optimizer mode
    mode_resolved = normalize_choice(
        opt_mode,
        param="--opt-mode",
        alias_groups=OPT_MODE_ALIASES,
        allowed_hint="grad|hess|lbfgs|rfo",
    )
    use_rfo = (mode_resolved == "rfo")

    try:
        config_layer_cfg = load_yaml_dict(config_yaml)
        override_layer_cfg = load_yaml_dict(override_yaml)
        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        opt_cfg = dict(OPT_BASE_KW)
        lbfgs_cfg = dict(LBFGS_KW)
        rfo_cfg = dict(RFO_KW)
        frequency_cfg = {"zero_cutoff_cm": FREQ_KW["zero_cutoff_cm"]}

        apply_yaml_overrides(
            config_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",), ("mlmm",))),
                (opt_cfg, (("opt",),)),
                (lbfgs_cfg, (("lbfgs",), ("opt", "lbfgs"))),
                (rfo_cfg, (("rfo",), ("opt", "rfo"))),
                (frequency_cfg, (("freq",),)),
            ],
        )

        if _is_param_explicit("max_cycles"):
            opt_cfg["max_cycles"] = int(max_cycles)
        if _is_param_explicit("dump"):
            opt_cfg["dump"] = bool(dump)
        if _is_param_explicit("out_dir"):
            opt_cfg["out_dir"] = out_dir
        if _is_param_explicit("thresh") and thresh is not None:
            opt_cfg["thresh"] = str(thresh)
        # --stop-plateau* rides the shared `opt` block, which the macro LBFGS/RFO
        # inherit. The MM micro relaxation never takes it (see _run_macro_micro).
        if _is_param_explicit("stop_plateau"):
            opt_cfg["energy_plateau"] = bool(stop_plateau)
        if stop_plateau_thresh is not None:
            opt_cfg["energy_plateau_thresh"] = float(stop_plateau_thresh)
        if stop_plateau_window is not None:
            opt_cfg["energy_plateau_window"] = int(stop_plateau_window)
        if _is_param_explicit("print_every") and print_every is not None:
            opt_cfg["print_every"] = int(print_every)
        if _is_param_explicit("cli_coord_type") and cli_coord_type is not None:
            geom_cfg["coord_type"] = str(cli_coord_type).lower()
        if _is_param_explicit("reject_uphill"):
            rfo_cfg["reject_uphill"] = bool(reject_uphill)

        if _is_param_explicit("detect_layer"):
            calc_cfg["use_bfactor_layers"] = bool(detect_layer)
        if _is_param_explicit("hess_cutoff") and hess_cutoff is not None:
            calc_cfg["hess_cutoff"] = float(hess_cutoff)
        if _is_param_explicit("movable_cutoff") and movable_cutoff is not None:
            calc_cfg["movable_cutoff"] = float(movable_cutoff)
            calc_cfg["use_bfactor_layers"] = False

        # CLI-resolved charge/spin (from -q / -l derivation in resolve_charge_spin_or_raise,
        # or -m / spin_default) always wins over the CALC_KW default carried in calc_cfg.
        # YAML calc.model_charge would have been merged earlier; the explicit user CLI
        # intent (or -l-derived total) supersedes it.
        calc_cfg["model_charge"] = int(charge)
        calc_cfg["model_mult"] = int(spin)
        if model_pdb is not None:
            calc_cfg["model_pdb"] = str(model_pdb)
        calc_cfg["input_pdb"] = str(prepared_input.source_path)
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
                (rfo_cfg, (("rfo",), ("opt", "rfo"))),
                (frequency_cfg, (("freq",),)),
            ],
        )
        frequency_zero_cutoff_cm = normalize_frequency_zero_cutoff_cm(
            frequency_cfg["zero_cutoff_cm"]
        )
        model_pdb_cfg = calc_cfg.get("model_pdb")
        # Revalidate after the highest-precedence YAML layer so dry-run and
        # real execution share the constructor's strict method vocabulary.
        apply_workers_to_calc_cfg(calc_cfg, None, None)
        try:
            geom_cfg["tr_projection"] = normalize_tr_projection_mode(
                geom_cfg.get("tr_projection")
            )
        except ValueError as exc:
            prepared_input.cleanup()
            raise click.ClickException(str(exc)) from exc

        # DLC is only meaningful with Hessian-based microiteration (ML region in
        # internal coordinates, MM as a Cartesian twin). Under plain L-BFGS
        # (--opt-mode grad) it would build delocalized internals over the whole
        # ML/MM system, which is needlessly slow with no benefit; fall back to cart.
        if not use_rfo and str(geom_cfg.get("coord_type", "cart")).lower() == "dlc":
            click.echo(
                "[opt] --coord-type dlc needs Hessian-based optimization "
                "(--opt-mode hess); L-BFGS runs in Cartesian — falling back to cart."
            )
            geom_cfg["coord_type"] = "cart"

        calc_paths = (("calc",), ("mlmm",))
        partial_explicit = (
            yaml_section_has_key(config_layer_cfg, calc_paths, "return_partial_hessian")
            or yaml_section_has_key(override_layer_cfg, calc_paths, "return_partial_hessian")
        )
        if not partial_explicit:
            calc_cfg["return_partial_hessian"] = True

        try:
            geom_freeze = _normalize_geom_freeze(geom_cfg.get("freeze_atoms"))
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

        out_dir_path = Path(opt_cfg["out_dir"]).resolve()
        error_out_dir = out_dir_path
        _reject_opt_output_collisions(
            out_dir_path,
            (
                input_path,
                prepared_input.original_path,
                prepared_input.source_path,
                geom_input_path,
                ref_pdb,
                real_parm7,
                Path(model_pdb_cfg) if model_pdb_cfg else None,
                config_yaml,
                override_yaml,
                Path(calc_cfg["calc_file"]) if calc_cfg.get("calc_file") else None,
            ),
        )

        # movable_cutoff implies full distance-based layer assignment.
        # hess_cutoff alone can be combined with --detect-layer.
        detect_layer_enabled = bool(calc_cfg.get("use_bfactor_layers", True))
        if movable_cutoff is not None:
            if detect_layer_enabled:
                click.echo("[layer] --movable-cutoff provided; disabling --detect-layer.", err=True)
            detect_layer_enabled = False
            calc_cfg["use_bfactor_layers"] = False

        layer_source_pdb = prepared_input.source_path
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

        mode_str = "RFO (hess)" if use_rfo else "LBFGS (grad)"

        if dry_run:
            from mlmm.core.utils import calculator_run_label, echo_run_summary
            echo_run_summary({
                "input": str(input_path),
                "backend": calculator_run_label(calc_cfg),
                "opt": f"{mode_str}, max_cycles={opt_cfg.get('max_cycles', '?')}",
                "out": str(out_dir_path),
            })
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
                        "optimizer_mode": "rfo" if use_rfo else "lbfgs",
                        "detect_layer": bool(detect_layer_enabled),
                        "model_region_source": model_region_source,
                        "model_indices_count": 0 if not model_indices else len(model_indices),
                        "tr_projection": geom_cfg["tr_projection"],
                        "will_run_optimization": True,
                        "will_convert_outputs": True,
                        "backend": calc_cfg.get("backend", "uma"),
                        "embedcharge": bool(calc_cfg.get("embedcharge", False)),
                    },
                )
            )
            click.echo("[dry-run] Validation complete. Optimization execution was skipped.")
            emit(
                format_elapsed("[time] Elapsed Time for Opt", time_start),
                narrative=True,
            )
            return

        _invalidate_opt_optional_outputs(out_dir_path)

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
                protected_inputs=(
                    input_path,
                    prepared_input.original_path,
                    prepared_input.source_path,
                    geom_input_path,
                    ref_pdb,
                    real_parm7,
                    Path(model_pdb_cfg) if model_pdb_cfg else None,
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
            prepared_input.cleanup()
            sys.exit(1)

        # When layer detection is enabled, also freeze frozen-layer atoms at the
        # optimizer geometry level (not only inside the calculator).
        # Otherwise LBFGS may still move those coordinates through coupled
        # inverse-Hessian updates, even if raw forces are zeroed there.
        if layer_info is not None:
            frozen_from_layer = [int(i) for i in layer_info.get("frozen_indices", [])]
            if frozen_from_layer:
                before = set(freeze_atoms_final)
                merged = sorted(before | set(frozen_from_layer))
                added = len(set(merged) - before)
                freeze_atoms_final = merged
                geom_cfg["freeze_atoms"] = freeze_atoms_final
                calc_cfg["freeze_atoms"] = freeze_atoms_final
                click.echo(
                    f"[layer] Applied optimizer freeze constraints: "
                    f"total={len(freeze_atoms_final)} (added_from_layer={added})"
                )

        # Distance-based overrides for Hessian-target and movable MM selection.
        hess_cutoff_final = calc_cfg.get("hess_cutoff")
        movable_cutoff_final = calc_cfg.get("movable_cutoff")
        if hess_cutoff_final is not None or movable_cutoff_final is not None:
            click.echo(
                f"[layer] Applied distance cutoffs: "
                f"hess={hess_cutoff_final} Å, freeze={movable_cutoff_final} Å"
            )
        from mlmm.workflows.freq import _align_three_layer_hessian_targets as _freq_align_three_layer_hessian_targets
        _freq_align_three_layer_hessian_targets(calc_cfg, echo_fn=click.echo)

        for key in ("input_pdb", "real_parm7", "model_pdb", "mm_fd_dir"):
            val = calc_cfg.get(key)
            if val:
                calc_cfg[key] = str(Path(val).expanduser().resolve())

        # Default-verbosity entry summary (skipped in child mode).
        from mlmm.core.utils import calculator_run_label, echo_run_summary
        echo_run_summary({
            "input": str(input_path),
            "backend": calculator_run_label(calc_cfg),
            "opt": f"{mode_str}, max_cycles={opt_cfg.get('max_cycles', '?')}",
            "out": str(out_dir_path),
        })

        click.echo(f"\n[mode] Optimizer: {mode_str}\n")
        click.echo(pretty_block("geom", format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")))
        echo_calc = format_freeze_atoms_for_echo(filter_calc_for_echo(calc_cfg), key="freeze_atoms")
        click.echo(pretty_block("calc", echo_calc))
        # Show only non-default opt settings
        echo_opt = strip_inherited_keys({**opt_cfg, "out_dir": str(out_dir_path)}, OPT_BASE_KW, mode="same")
        click.echo(pretty_block("opt", echo_opt))
        # Show only optimizer-specific settings, not inherited from opt_cfg
        if use_rfo:
            echo_rfo = strip_inherited_keys(rfo_cfg, opt_cfg)
            click.echo(pretty_block("rfo", echo_rfo))
        else:
            echo_lbfgs = strip_inherited_keys(lbfgs_cfg, opt_cfg)
            click.echo(pretty_block("lbfgs", echo_lbfgs))
        # Resolve effective bias_k: CLI value wins, else BIAS_KW['k'] default.
        # CLI default flipped to None so that downstream consumers (--dist-freeze
        # harmonic restraint) fall back to the single source of truth in
        # defaults.py rather than a hardcoded 300 in the @click.option.
        bias_k_eff = float(bias_k) if bias_k is not None else float(BIAS_KW["k"])

        if dist_freeze:
            display_pairs = []
            for (i, j, target) in dist_freeze:
                label = (f"{target:.4f}" if target is not None else "<current>")
                display_pairs.append((int(i) + 1, int(j) + 1, label))
            click.echo(
                pretty_block(
                    "dist_freeze (input)",
                    {
                        "k (eV/Å^2)": bias_k_eff,
                        "pairs_1based": display_pairs,
                    },
                )
            )

        out_dir_path.mkdir(parents=True, exist_ok=True)
        coord_type = geom_cfg.get("coord_type", "cart")
        coord_kwargs = dict(geom_cfg)
        coord_kwargs.pop("coord_type", None)
        geometry = geom_loader(
            geom_input_path,
            coord_type=coord_type,
            **coord_kwargs,
        )

        if mm_only:
            if use_rfo:
                click.echo(
                    "ERROR: --mm-only is incompatible with --opt-mode hess (RFO needs a Hessian, "
                    "but the MM-only calculator does not provide one). Use --opt-mode grad.",
                    err=True,
                )
                sys.exit(1)
            if microiter:
                click.echo("[opt] --mm-only: microiteration disabled (no ML component to alternate with).")
                microiter = False
            mm_core_calc = mlmm(
                **calc_cfg,
                _skip_high_level_backend=True,
            )
            base_calc = mlmm_mm_only(
                mm_core_calc.core,
                freeze_atoms=freeze_atoms_final,
            )
            click.echo("[opt] --mm-only: MLIP component skipped; minimizing on MM force field only.")
        else:
            base_calc = mlmm(**calc_cfg)
        geometry.set_calculator(base_calc)

        echo_resolved_device()

        resolved_dist_freeze: List[Tuple[int, int, float]] = []
        active_calc = base_calc
        if dist_freeze:
            try:
                resolved_dist_freeze = _resolve_dist_freeze_targets(geometry, dist_freeze)
            except click.BadParameter as e:
                click.echo(f"ERROR: {e}", err=True)
                sys.exit(1)
            click.echo(
                pretty_block(
                    "dist_freeze (active)",
                    {
                        "k (eV/Å^2)": bias_k_eff,
                        "pairs_1based": [
                            (int(i) + 1, int(j) + 1, float(f"{t:.4f}"))
                            for (i, j, t) in resolved_dist_freeze
                        ],
                    },
                )
            )
            bias_calc = HarmonicBiasCalculator(base_calc, k=bias_k_eff)
            bias_calc.set_pairs(resolved_dist_freeze)
            active_calc = bias_calc
            geometry.set_calculator(active_calc)

        # Pass only opt-level values that differ from OPT_BASE defaults, so
        # optimizer-specific YAML (e.g. rfo.print_every / lbfgs.print_every)
        # is not overwritten by inherited defaults such as opt.print_every=100.
        common_kwargs = strip_inherited_keys(dict(opt_cfg), OPT_BASE_KW, mode="same")
        common_kwargs["out_dir"] = str(out_dir_path)

        def _build_optimizer(run_kind: str):
            if run_kind == "lbfgs":
                lbfgs_args = {**lbfgs_cfg, **common_kwargs}
                return LBFGS(geometry, **lbfgs_args)
            if run_kind == "rfo":
                rfo_args = {**rfo_cfg, **common_kwargs}
                return RFOptimizer(geometry, **rfo_args)
            raise click.BadParameter(f"Unknown optimizer kind '{run_kind}'.")

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

        # Serialize state through one set of names (no dir-based discovery).
        optimizer = None
        microiter_result = None
        microiter_partition = None
        microiter_fallback_reason = None

        use_microiter = bool(microiter) and use_rfo and not dist_freeze
        if bool(microiter) and not use_rfo:
            microiter_fallback_reason = "requires_hessian_mode"
            click.echo("[microiter] --microiter is only effective with --opt-mode hess (RFO). Ignoring.")
        if bool(microiter) and use_rfo and dist_freeze:
            microiter_fallback_reason = "distance_restraints"
            click.echo("[microiter] --microiter is not compatible with --dist-freeze. Falling back to standard RFO.")

        if use_microiter:
            # resolve the ONE immutable partition strictly from the accepted
            # calculator core BEFORE dispatch. A construction/partition failure
            # raises PartitionError (loud, ordinary error envelope); a VALID
            # empty-ML partition is a documented fallback to a REAL standard RFO
            # run, never an unchanged geometry returned as a completed result.
            _mi_core = base_calc.core if hasattr(base_calc, "core") else base_calc
            microiter_partition = resolve_partition_from_core(
                _mi_core, len(geometry.atoms), geometry.freeze_atoms
            )
            if not microiter_partition.has_macro_active:
                microiter_fallback_reason = "no_ml_atoms"
                click.echo(
                    "[microiter] No ML macro-active atoms in the resolved "
                    "partition; running a standard RFO optimization "
                    "(fallback_reason=no_ml_atoms)."
                )
                use_microiter = False

        if use_microiter:
            emit("\n====== Optimization (RFO + Microiteration) ======\n", narrative=True)
            microiter_result = _run_microiter_opt(
                geometry,
                base_calc,
                calc_cfg,
                rfo_cfg,
                lbfgs_cfg,
                opt_cfg,
                microiter_cfg,
                out_dir_path,
                partition=microiter_partition,
                dump=bool(opt_cfg["dump"]),
            )
            if microiter_result is not None:
                emit_optimizer_terminal_status(
                    "opt",
                    converged=microiter_result.get("converged"),
                    cycles=microiter_result.get("cycles"),
                    max_cycles=opt_cfg.get("max_cycles"),
                    stalled=bool(microiter_result.get("is_stalled")),
                    stop_reason=microiter_result.get("stop_reason") or None,
                )

            # Write final geometry
            from ase import Atoms as _Atoms
            from ase.io import write as _write
            final_xyz_path = out_dir_path / "final_geometry.xyz"
            final_coords_ang = geometry.coords3d * BOHR2ANG
            atoms_final = _Atoms(geometry.atoms, positions=final_coords_ang, pbc=False)
            _write(final_xyz_path, atoms_final)

        else:
            main_kind = "rfo" if use_rfo else "lbfgs"
            if use_rfo:
                _seed_rfo_initial_hessian(
                    geometry,
                    calc_cfg,
                    active_calc,
                    restraints_active=bool(resolved_dist_freeze),
                )

            main_label = "RFO" if use_rfo else "LBFGS"
            optimizer = _build_optimizer(main_kind)
            emit(f"\n====== Optimization ({main_label}) ======\n", narrative=True)
            optimizer.run()
            emit_optimizer_terminal_status(
                "opt",
                converged=getattr(optimizer, "is_converged", None),
                cycles=optimizer_cycle_count(optimizer),
                max_cycles=opt_cfg.get("max_cycles"),
                stalled=getattr(optimizer, "is_stalled", False),
                stop_reason=getattr(optimizer, "stop_reason", None) or None,
            )

            # Get final geometry path
            final_xyz_path = optimizer.final_fn if isinstance(optimizer.final_fn, Path) else Path(optimizer.final_fn)

            if bool(opt_cfg["dump"]):
                optim_all_path = out_dir_path / "optimization_all_trj.xyz"
                if not optim_all_path.exists():
                    trj_path = optimizer.get_path_for_fn("optimization_trj.xyz")
                    _append_xyz_trajectory(optim_all_path, trj_path, reset=True)

        # The terminal owner is the run that produced the coordinates ultimately
        # written to final_geometry. A flatten retry supersedes this initial owner.
        terminal_optimizer = optimizer
        terminal_microiter_result = microiter_result
        terminal_use_microiter = use_microiter

        rigid_projection_info: Dict[str, Any] = {}

        # track a real energy-plateau stall from whichever optimizer
        # ran (standard, microiteration macro/latest-micro, or a flatten retry
        # below). A stall stops further flatten/retry work and is reported as a
        # distinct, non-converged outcome (never converged).
        _opt_stalled = False
        _opt_stop_reason = ""
        if use_microiter:
            if microiter_result is not None:
                _opt_stalled = bool(microiter_result.get("is_stalled"))
                _opt_stop_reason = microiter_result.get("stop_reason") or ""
        elif optimizer is not None:
            _opt_stalled = bool(getattr(optimizer, "is_stalled", False))
            _opt_stop_reason = getattr(optimizer, "stop_reason", "") or ""

        # Flatten loop (all imaginary modes).  A stalled optimization is
        # precisely when this is wanted: it rebuilds the Hessian and displaces
        # along the remaining imaginary modes to leave the plateau.  's
        # no-retry rule belongs inside the loop (a flatten *retry* that stalls
        # again stops there and sets ``_opt_stalled``), not in front of it.
        if flatten:
            from mlmm.workflows.freq import (
                _torch_device,
                _calc_full_hessian_torch,
                _frequencies_cm_and_modes,
                _safe_masses_amu,
                _active_atoms_from_partial_hessian_metadata,
            )

            emit("\n====== Optimization (Flatten loop) ======\n", narrative=True)

            geometry.set_calculator(None)
            calc_kwargs_for_flatten = dict(calc_cfg)
            calc_kwargs_for_flatten["out_hess_torch"] = True
            device = _torch_device(calc_cfg.get("ml_device", "auto"))
            freeze_idx = list(geom_cfg.get("freeze_atoms", [])) if len(geom_cfg.get("freeze_atoms", [])) > 0 else None
            masses_amu = _safe_masses_amu(geometry.atomic_numbers)

            def _attach_opt_calc() -> None:
                geometry.set_calculator(active_calc)

            def _calc_freqs_and_modes() -> Tuple[np.ndarray, torch.Tensor]:
                # Refresh the active-DOF metadata used by PHVA routing.
                H, _e = _calc_full_hessian_torch(
                    geometry,
                    calc_kwargs_for_flatten,
                    device,
                    refresh_geom_meta=True,
                    calculator=active_calc,
                )
                effective_freeze_idx = freeze_idx
                if H.shape[0] != 3 * len(geometry.atomic_numbers):
                    active_atoms = _active_atoms_from_partial_hessian_metadata(
                        geometry, int(H.shape[0])
                    )
                    if active_atoms is None:
                        raise RuntimeError(
                            "Partial Hessian metadata does not identify its active atoms."
                        )
                    active_set = set(active_atoms)
                    effective_freeze_idx = [
                        i for i in range(len(geometry.atomic_numbers))
                        if i not in active_set
                    ]
                freqs_local, modes_local = _frequencies_cm_and_modes(
                    H,
                    geometry.atomic_numbers,
                    geometry.cart_coords.reshape(-1, 3),
                    device,
                    freeze_idx=effective_freeze_idx,
                    tr_projection=geom_cfg["tr_projection"],
                    projection_info=rigid_projection_info,
                    frequency_zero_cutoff_cm=frequency_zero_cutoff_cm,
                )
                rigid_projection_info.update({
                    "hessian_space": (
                        "full" if H.shape[0] == 3 * len(geometry.atomic_numbers)
                        else "active"
                    ),
                    "raw_hessian_shape": list(H.shape),
                    "source": "opt_flatten",
                })
                del H
                return freqs_local, modes_local

            freqs_cm, modes = _calc_freqs_and_modes()
            neg_mask = resolved_imaginary_mask(
                freqs_cm, frequency_zero_cutoff_cm
            )
            n_imag = int(np.sum(neg_mask))
            ims = [float(x) for x in freqs_cm[neg_mask]]
            emit(f"[Imaginary modes] n={n_imag} ({ims})", narrative=True)

            flatten_kind = mode_resolved  # reuse same optimizer type
            for it in range(OPT_FLATTEN_MAX_ITER):
                if n_imag == 0:
                    break
                click.echo(f"[flatten] iteration {it + 1}/{OPT_FLATTEN_MAX_ITER}")
                did_flatten = _flatten_all_imag_modes_for_geom(
                    geometry,
                    masses_amu,
                    calc_kwargs_for_flatten,
                    freqs_cm,
                    modes,
                    frequency_zero_cutoff_cm,
                    OPT_FLATTEN_AMP_ANG,
                    calculator=active_calc,
                )
                if not did_flatten:
                    click.echo("[flatten] No eligible imaginary modes to flatten; stopping.")
                    break

                _attach_opt_calc()
                opt_restart = _build_optimizer(flatten_kind)
                restart_label = "LBFGS" if flatten_kind == "lbfgs" else "RFO"
                emit(f"\n====== Optimization ({restart_label}, flatten retry) ======\n", narrative=True)
                opt_restart.run()
                terminal_optimizer = opt_restart
                terminal_microiter_result = None
                terminal_use_microiter = False
                emit_optimizer_terminal_status(
                    "opt",
                    converged=getattr(opt_restart, "is_converged", None),
                    cycles=optimizer_cycle_count(opt_restart),
                    max_cycles=opt_cfg.get("max_cycles"),
                    stalled=getattr(opt_restart, "is_stalled", False),
                    stop_reason=getattr(opt_restart, "stop_reason", None) or None,
                )

                # Stop retrying a stalled optimization : a flatten
                # retry that stalled is not making progress, so re-running it
                # would only repeat the stall.
                _opt_stalled = bool(getattr(opt_restart, "is_stalled", False))
                _opt_stop_reason = (
                    getattr(opt_restart, "stop_reason", "") if _opt_stalled else ""
                )
                if _opt_stalled:
                    click.echo(
                        "[flatten] Optimization stalled (energy plateau); "
                        "stopping the flatten loop."
                    )
                    break

                geometry.set_calculator(None)
                freqs_cm, modes = _calc_freqs_and_modes()
                neg_mask = resolved_imaginary_mask(
                    freqs_cm, frequency_zero_cutoff_cm
                )
                n_imag = int(np.sum(neg_mask))
                ims = [float(x) for x in freqs_cm[neg_mask]]
                emit(f"[Imaginary modes] n={n_imag} ({ims})", narrative=True)

            if n_imag > 0:
                click.echo(
                    f"[flatten] WARNING: Remaining imaginary modes after {OPT_FLATTEN_MAX_ITER} iterations: {n_imag}",
                    err=True,
                )
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            projection_block = pretty_block(
                "rigid_projection", rigid_projection_info
            )
            if projection_block:
                click.echo(projection_block)

            # Update final geometry after flatten
            final_xyz_path = out_dir_path / "final_geometry.xyz"
            from ase import Atoms as _Atoms
            from ase.io import write as _write
            final_coords_ang = geometry.coords3d * BOHR2ANG
            atoms_final = _Atoms(geometry.atoms, positions=final_coords_ang, pbc=False)
            _write(final_xyz_path, atoms_final)

        # Extract layer indices from calculator for layer-based B-factor encoding
        calc_core = base_calc.core if hasattr(base_calc, 'core') else base_calc
        ml_indices = getattr(calc_core, 'ml_indices', None)
        hess_mm_indices = getattr(calc_core, 'hess_mm_indices', None)
        movable_mm_indices = getattr(calc_core, 'movable_mm_indices', None)
        frozen_layer_indices = getattr(calc_core, 'frozen_layer_indices', None)

        _maybe_convert_outputs_to_pdb(
            input_path=prepared_input.source_path,  # Use PDB topology for conversion
            out_dir=out_dir_path,
            dump=bool(opt_cfg["dump"]),
            get_trj_fn=(
                (lambda fn: out_dir_path / fn)
                if terminal_use_microiter
                else terminal_optimizer.get_path_for_fn
            ),
            final_xyz_path=final_xyz_path,
            model_pdb=Path(calc_cfg["model_pdb"]),
            freeze_indices_0based=freeze_atoms_final,
            ml_indices=ml_indices,
            hess_mm_indices=hess_mm_indices,
            movable_mm_indices=movable_mm_indices,
            frozen_layer_indices=frozen_layer_indices,
        )

        if out_json:
            from mlmm.core.utils import calculator_provenance, write_result_json
            _opt_converged = _opt_terminal_converged(
                terminal_use_microiter,
                terminal_microiter_result,
                terminal_optimizer,
            )
            # n_opt_cycles is the EXECUTED macro cycle count, never the
            # configured budget. On the microiteration path there is no standalone
            # optimizer in scope, so the executed macro cycles come from the
            # driver's outcome. Preserve the actual microiteration cycle count.
            # The ordinary path reports ``cur_cycle + 1`` (executed cycles). A
            # one-cycle converged optimization
            # therefore reported n_opt_cycles=0 in JSON but "Total cycles: 1" in the
            # log. Use ``optimizer_cycle_count`` so JSON == log == tsopt.
            if terminal_use_microiter and terminal_microiter_result is not None:
                _opt_cycles = int(terminal_microiter_result.get("cycles", 0))
            elif terminal_optimizer is not None and hasattr(terminal_optimizer, "cur_cycle"):
                _opt_cycles = optimizer_cycle_count(terminal_optimizer)
            else:
                _opt_cycles = None
            final_energy_hartree = unbiased_energy_hartree(geometry, base_calc)
            # an energy-plateau stall is a distinct, additive outcome
            # that is never reported as converged.  ``converged`` / ``not_converged``
            # remain byte-compatible; only ``stalled`` is new.
            provenance = calculator_provenance(calc_cfg)
            if mm_only:
                provenance.update(
                    {
                        "mlip_backend": None,
                        "mlip_model": None,
                        "mlip_precision": None,
                    }
                )
            result_data = {
                "status": "stalled" if _opt_stalled else ("converged" if _opt_converged else "not_converged"),
                "energy_hartree": final_energy_hartree,
                "n_opt_cycles": _opt_cycles,
                "opt_mode": opt_cfg.get("opt_mode", opt_mode),
                **provenance,
                "charge": calc_cfg.get("model_charge"),
                "spin": calc_cfg.get("model_mult"),
                "n_atoms": len(geometry.atoms),
                "n_freeze_atoms": len(geom_cfg.get("freeze_atoms", [])),
                "thresh": opt_cfg.get("thresh", "gau"),
                "max_cycles": opt_cfg.get("max_cycles"),
                "input_file": str(prepared_input.source_path),
                "files": {
                    "final_geometry_xyz": str(final_xyz_path.name),
                },
            }
            # Additive stop_reason, present only for a non-converged stop
            # (stalled/stopped) so a genuinely converged run's JSON stays
            # byte-compatible.
            if _opt_stop_reason:
                result_data["stop_reason"] = _opt_stop_reason
            # additive microiteration serialization. Executed macro cycles
            # already populate n_opt_cycles; n_micro_cycles and the microiteration
            # object carry the separate micro totals + macro/micro leaf outcomes.
            # These are additive: a converged run's legacy keys are unchanged.
            if terminal_use_microiter and terminal_microiter_result is not None:
                _mi_outcome = terminal_microiter_result.get("outcome")
                if _mi_outcome is not None:
                    result_data["n_micro_cycles"] = int(
                        terminal_microiter_result.get("micro_cycles", 0)
                    )
                    result_data["microiteration"] = _mi_outcome.to_result_object()
            elif bool(microiter) and microiter_fallback_reason:
                result_data["microiteration"] = {
                    "requested": True,
                    "used": False,
                    "fallback_reason": microiter_fallback_reason,
                }
            if rigid_projection_info:
                result_data["rigid_projection"] = dict(rigid_projection_info)
            # Final force convergence values
            if (
                terminal_optimizer is not None
                and hasattr(terminal_optimizer, "max_forces")
                and terminal_optimizer.max_forces
            ):
                result_data["final_max_force"] = float(terminal_optimizer.max_forces[-1])
                result_data["final_rms_force"] = float(terminal_optimizer.rms_forces[-1])
            # Convergence thresholds (numeric values for the named preset)
            if (
                terminal_optimizer is not None
                and hasattr(terminal_optimizer, "convergence")
                and terminal_optimizer.convergence
            ):
                result_data["convergence_thresholds"] = {
                    k: float(v)
                    for k, v in terminal_optimizer.convergence.items()
                }
            # Final step convergence values
            if (
                terminal_optimizer is not None
                and hasattr(terminal_optimizer, "max_steps")
                and terminal_optimizer.max_steps
            ):
                result_data["final_max_step"] = float(terminal_optimizer.max_steps[-1])
                result_data["final_rms_step"] = float(terminal_optimizer.rms_steps[-1])
            # Add PDB/GJF if generated
            for ext in (".pdb", ".gjf"):
                f = out_dir_path / f"final_geometry{ext}"
                if f.exists():
                    result_data["files"][f"final_geometry_{ext[1:]}"] = f.name
            # Add trajectory files if they exist
            for name in ("optimization_trj.xyz", "optimization.pdb"):
                _tf = out_dir_path / name
                if _tf.exists():
                    key = name.replace(".", "_").replace("-", "_")
                    result_data["files"][key] = name
            write_result_json(
                out_dir_path, result_data,
                command="opt",
                elapsed_seconds=time.perf_counter() - time_start,
            )

        emit(
            format_elapsed("[time] Elapsed Time for Opt", time_start),
            narrative=True,
        )

    except ZeroStepLength as e:
        _write_error_json(Path(out_dir).resolve(), "opt", e, "ZeroStepLength", time_start)
        click.echo("ERROR: Step length fell below the minimum allowed (ZeroStepLength).", err=True)
        sys.exit(2)
    except OptimizationError as e:
        _write_error_json(Path(out_dir).resolve(), "opt", e, "OptimizationError", time_start)
        click.echo(f"ERROR: Optimization failed - {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except _OptOutputCollisionError:
        raise
    except Exception as e:
        render_cli_exception(e, label="optimization", out_dir=error_out_dir, command="opt", time_start=time_start)
    finally:
        if prepared_input is not None:
            prepared_input.cleanup()
        # Release GPU memory so subsequent pipeline stages don't OOM.
        # `= None` decref's the heavy refs; `del` then removes names from
        # the local frame so torch.nn.Module hooks / closures cannot retain.
        base_calc = bias_calc = active_calc = geometry = optimizer = mm_calc = macro_calc = macro_optimizer = None
        del base_calc, bias_calc, active_calc, geometry, optimizer, mm_calc, macro_calc, macro_optimizer
        gc.collect()  # break cyclic refs inside torch.nn.Module
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


# Allow `python -m mlmm.opt` direct execution
if __name__ == "__main__":
    cli()
