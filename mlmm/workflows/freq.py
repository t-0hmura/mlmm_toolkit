"""ML/MM vibrational analysis with PHVA and thermochemistry."""
# DOMAIN_PURE

from __future__ import annotations

import gc
import logging
import sys
import textwrap
import time

logger = logging.getLogger(__name__)
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import click
from mlmm.core.output import emit
import numpy as np
import torch
import ase.units as units
import yaml
from ase import Atoms
from ase.data import atomic_masses
from ase.io import write

from pysisyphus.constants import AMU2AU, ANG2BOHR, AU2EV, BOHR2ANG
from pysisyphus._array import active_square
from pysisyphus.helpers import geom_loader
from pysisyphus.tr_projection import active_tr_basis, project_hessian_inplace

# Compatibility re-exports: the pure normal-mode kernel lives in the
# lower bundled-engine module ``pysisyphus.normal_modes`` (a sibling of
# ``pysisyphus.tr_projection``), which does not import ``mlmm``. It is
# re-exported here so existing callers of
# ``mlmm.workflows.freq`` (opt/tsopt/tests) keep working unchanged and resolve to
# the same function objects as the lower implementation.
#
# CHEMISTRY-RULE:6 PHVA + MLIP active-block: mass-weighted Hessian only;
# TR projection is applied separately downstream. The kernel lives in
# ``pysisyphus.normal_modes``, but the rule stays this workflow's
# responsibility, so the marker is placed at the consuming site.
from pysisyphus.normal_modes import (  # noqa: F401
    _safe_masses_amu,
    _mw_projected_hessian,
    _mass_weighted_hessian,
    _frequencies_cm_and_modes,
    _mw_mode_to_cart,
)

from mlmm.backends.mlmm_calc import mlmm
from mlmm.core.defaults import FREQ_KW, THERMO_KW
# Import shared layer/keyword helpers from the neutral ``_opt_freq_common``
# module and ``_parse_freeze_atoms`` from its canonical home in core.utils.
from mlmm.workflows._opt_freq_common import (
    CALC_KW as OPT_CALC_KW,
    GEOM_KW as OPT_GEOM_KW,
    _normalize_geom_freeze as _normalize_geom_freeze_opt,
    _convert_yaml_layer_atoms_1to0,
)
from mlmm.core.utils import _parse_freeze_atoms as _parse_freeze_atoms_opt
from mlmm.workflows.charge_prep import resolve_charge_spin_or_raise
from mlmm.core.utils import (
    apply_ref_pdb_override,
    apply_layer_freeze_constraints,
    apply_yaml_overrides,
    convert_xyz_to_pdb,
    set_convert_file_enabled,
    is_convert_file_enabled,
    filter_calc_for_echo,
    format_elapsed,
    format_freeze_atoms_for_echo,
    load_yaml_dict,
    merge_freeze_atom_indices,
    normalize_choice,
    prepare_input_structure,
    pretty_block,
    parse_indices_string,
    resolve_ml_layer_assignment,
    strip_inherited_keys,
    yaml_section_has_key,
    echo_resolved_device,
)
from mlmm.cli.common_options import add_ml_charge_spin_options, add_ml_layer_detection_options, add_precision_option, add_workers_options, add_backend_model_option, add_calc_file_option, add_deterministic_option, add_allow_charge_mult_mismatch_option
from mlmm.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, make_is_param_explicit, render_cli_exception


def _torch_device(auto: str = "auto") -> torch.device:
    if auto == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(auto)


def _calc_full_hessian_torch(
    geom,
    calc_kwargs: Dict[str, Any],
    device: torch.device,
    *,
    refresh_geom_meta: bool = False,
    calculator=None,
) -> Tuple[torch.Tensor, float]:
    """Return (Hessian torch tensor, energy Hartree) for the active PES.

    ``calculator`` may supply an already-resolved evaluator, including a
    restraint wrapper. When omitted, a temporary ML/MM calculator is built
    from ``calc_kwargs``.
    """

    kw = dict(calc_kwargs or {})
    kw["out_hess_torch"] = True
    owns_calc = calculator is None
    calc = mlmm(**kw) if owns_calc else calculator
    # Pass Cartesian coords explicitly (geom.coords returns the active coord
    # space which is INTERNAL for coord_type ∈ {redund, dlc, tric}, breaking
    # mlmm_calc._run_core's `np.asarray(coords).reshape(-1, 3)` assumption).
    # cart_coords is always 3N Cartesian regardless of geom.coord_type.
    result = calc.get_hessian(geom.atoms, geom.cart_coords)

    if refresh_geom_meta:
        within = result.get("within_partial_hessian")
        if within is None and kw.get("return_partial_hessian"):
            try:
                core = getattr(calc, "core", None)
                if core is not None and hasattr(core, "_build_within_partial_hessian"):
                    within = core._build_within_partial_hessian()
            except Exception:
                within = None
        if within is not None:
            geom.within_partial_hessian = within
        elif "hessian" in result:
            geom.within_partial_hessian = None

        try:
            core = getattr(calc, "core", None)
            if core is not None and hasattr(core, "hess_active_atoms"):
                active_atoms = np.asarray(core.hess_active_atoms, dtype=int)
                geom._hess_active_atoms_last = active_atoms
                if active_atoms.size:
                    active_dofs = np.empty(active_atoms.size * 3, dtype=int)
                    for i, a in enumerate(active_atoms):
                        base = 3 * int(a)
                        active_dofs[3 * i:3 * i + 3] = (base, base + 1, base + 2)
                else:
                    active_dofs = np.zeros(0, dtype=int)
                geom._hess_active_dofs_last = active_dofs
        except Exception:
            logger.debug("Failed to extract active DOF info from calculator", exc_info=True)

    H = result["hessian"]
    if not isinstance(H, torch.Tensor):
        H = torch.as_tensor(H)
    H = H.to(device=device)
    if "energy" not in result:
        raise KeyError("Hessian result is missing 'energy'.")
    energy = float(result["energy"])
    if not np.isfinite(energy):
        raise ValueError("Hessian energy must be finite for thermochemistry.")

    del result
    if owns_calc:
        del calc
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    return H, energy


def _ordered_hessian_coverage_atoms(
    geom: Any,
    n_atoms: int,
) -> Optional[List[int]]:
    """Return the ordered atoms whose curvature was actually evaluated.

    Coverage metadata is independent of storage shape: ML/MM may scatter a
    partial Hessian into a zero-padded ``3N x 3N`` tensor.  Treating that tensor
    as fully evaluated solely because of its shape would certify artificial
    zero-curvature rows.
    """

    def _atoms_from_dofs(raw: Any) -> Optional[List[int]]:
        if raw is None:
            return None
        dofs = [int(d) for d in np.asarray(raw, dtype=int).reshape(-1).tolist()]
        if not dofs or len(dofs) % 3:
            return None
        atoms: List[int] = []
        for pos in range(0, len(dofs), 3):
            triplet = dofs[pos:pos + 3]
            atom = triplet[0] // 3
            if atom < 0 or atom >= int(n_atoms):
                return None
            if triplet != [3 * atom, 3 * atom + 1, 3 * atom + 2]:
                return None
            atoms.append(atom)
        if len(set(atoms)) != len(atoms):
            return None
        return atoms

    def _atoms_from_atoms(raw: Any) -> Optional[List[int]]:
        if raw is None:
            return None
        atoms = [int(a) for a in np.asarray(raw, dtype=int).reshape(-1).tolist()]
        if not atoms or len(set(atoms)) != len(atoms):
            return None
        if any(a < 0 or a >= int(n_atoms) for a in atoms):
            return None
        return atoms

    within = getattr(geom, "within_partial_hessian", None)
    if isinstance(within, dict):
        for candidate in (
            _atoms_from_dofs(within.get("active_dofs")),
            _atoms_from_atoms(within.get("active_atoms")),
        ):
            if candidate is not None:
                return candidate

    for candidate in (
        _atoms_from_dofs(getattr(geom, "_hess_active_dofs_last", None)),
        _atoms_from_atoms(getattr(geom, "_hess_active_atoms_last", None)),
    ):
        if candidate is not None:
            return candidate
    return None


def _active_atoms_from_partial_hessian_metadata(
    geom: Any,
    hessian_dim: int,
) -> Optional[List[int]]:
    """Return atoms spanning a compact Hessian, retaining the legacy helper API."""
    if hessian_dim <= 0 or hessian_dim % 3:
        return None
    n_atoms = len(getattr(geom, "atomic_numbers", []))
    if n_atoms == 0:
        within = getattr(geom, "within_partial_hessian", None)
        if isinstance(within, dict):
            raw_atoms = np.asarray(within.get("active_atoms", []), dtype=int).reshape(-1)
            raw_dofs = np.asarray(within.get("active_dofs", []), dtype=int).reshape(-1)
            maxima = []
            if raw_atoms.size:
                maxima.append(int(raw_atoms.max()) + 1)
            if raw_dofs.size:
                maxima.append(int(raw_dofs.max()) // 3 + 1)
            if maxima:
                n_atoms = max(maxima)
    coverage = _ordered_hessian_coverage_atoms(geom, n_atoms)
    if coverage is None or 3 * len(coverage) != int(hessian_dim):
        return None
    return coverage


def _reconcile_hessian_analysis_basis(
    hessian: torch.Tensor,
    geom: Any,
    requested_atoms: List[int],
) -> Tuple[torch.Tensor, List[int], List[int], str]:
    """Slice an evaluated Hessian onto the exact requested atom basis.

    A computed superset is sliced in its recorded order.  Missing requested
    atoms are a hard error: silently shrinking a PHVA basis changes the
    frequencies and can change transition-state order.
    """
    n_atoms = len(getattr(geom, "atomic_numbers", []))
    full_dim = 3 * n_atoms
    if hessian.ndim != 2 or hessian.shape[0] != hessian.shape[1]:
        raise click.ClickException(
            f"Hessian must be square; got shape={tuple(hessian.shape)}."
        )

    requested: List[int] = []
    seen = set()
    for raw in requested_atoms:
        atom = int(raw)
        if atom < 0 or atom >= n_atoms:
            raise click.ClickException(
                f"Requested Hessian atom index {atom} is outside 0..{n_atoms - 1}."
            )
        if atom not in seen:
            seen.add(atom)
            requested.append(atom)
    requested.sort()
    if not requested:
        raise click.ClickException("Frequency analysis requires at least one active atom.")

    coverage = _ordered_hessian_coverage_atoms(geom, n_atoms)
    h_dim = int(hessian.shape[0])
    if coverage is None:
        if h_dim != full_dim:
            raise click.ClickException(
                "A compact Hessian was returned without ordered coverage metadata."
            )
        coverage = list(range(n_atoms))

    missing = sorted(set(requested) - set(coverage))
    if missing:
        preview = ", ".join(str(i + 1) for i in missing[:12])
        suffix = " …" if len(missing) > 12 else ""
        raise click.ClickException(
            "The requested frequency-analysis basis is wider than the evaluated "
            f"Hessian coverage; missing 1-based atom indices: {preview}{suffix}. "
            "Increase the Hessian region or choose a narrower --active-dof-mode."
        )

    within = getattr(geom, "within_partial_hessian", None)
    declared_compact = (
        isinstance(within, dict)
        and h_dim == 3 * len(coverage)
        and (
            within.get("active_dofs") is not None
            or within.get("active_atoms") is not None
        )
    )
    if declared_compact:
        local = {atom: pos for pos, atom in enumerate(coverage)}
        source_dofs = [
            3 * local[atom] + axis for atom in requested for axis in range(3)
        ]
        storage = "compact"
    elif h_dim == full_dim:
        source_dofs = [3 * atom + axis for atom in requested for axis in range(3)]
        storage = "full"
    elif h_dim == 3 * len(coverage):
        local = {atom: pos for pos, atom in enumerate(coverage)}
        source_dofs = [
            3 * local[atom] + axis for atom in requested for axis in range(3)
        ]
        storage = "compact"
    else:
        raise click.ClickException(
            "Hessian shape is inconsistent with its coverage metadata: "
            f"shape={tuple(hessian.shape)}, coverage_atoms={len(coverage)}, "
            f"full_atoms={n_atoms}."
        )

    if source_dofs == list(range(h_dim)):
        analysis_hessian = hessian
    else:
        index = torch.as_tensor(
            source_dofs, dtype=torch.long, device=hessian.device
        )
        analysis_hessian = active_square(hessian, index)
        del index
    return analysis_hessian, requested, coverage, storage


def _record_hessian_result_path(files: Dict[str, str], path: Path) -> Dict[str, str]:
    """Record the exact Hessian artifact path under its legacy JSON key."""

    files["hessian_npz"] = str(Path(path))
    return files


def _collect_layer_atom_sets(calc_cfg: Dict[str, Any]) -> Dict[str, set[int]]:
    """Collect ML/MM layer index sets from a temporary calculator instance."""
    empty = {"ml": set(), "hess_mm": set(), "movable_mm": set(), "frozen_mm": set()}
    try:
        temp_calc = mlmm(**dict(calc_cfg))
        calc_core = temp_calc.core if hasattr(temp_calc, "core") else temp_calc
        layer_sets = {
            "ml": set(getattr(calc_core, "ml_indices", []) or []),
            "hess_mm": set(getattr(calc_core, "hess_mm_indices", []) or []),
            "movable_mm": set(getattr(calc_core, "movable_mm_indices", []) or []),
            "frozen_mm": set(getattr(calc_core, "frozen_layer_indices", []) or []),
        }
        del temp_calc
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        return layer_sets
    except Exception as exc:
        logger.debug(
            "_collect_layer_atom_sets: failed to construct temp mlmm(); "
            "returning empty layer sets",
            exc_info=True,
        )
        # Empty layer sets make every --active-dof-mode collapse to ALL atoms.
        # A debug-level line is invisible on a normal run, so the user would
        # read a full-system Hessian as the ml-only one they asked for.
        click.echo(
            f"[active-dof] WARNING: could not resolve ML/MM layer sets ({exc}); "
            "the frequency analysis falls back to ALL atoms regardless of "
            "--active-dof-mode.",
            err=True,
        )
        return empty


def _align_three_layer_hessian_targets(
    calc_cfg: Dict[str, Any],
    *,
    echo_fn=None,
) -> bool:
    """
    In 3-layer detect-layer mode, align Hessian targets to MovableMM by default.

    Returns True when a default policy was applied.
    """
    if calc_cfg.get("hess_cutoff") is not None:
        return False

    explicit_layer_lists = any(
        calc_cfg.get(key) is not None
        for key in ("hess_mm_atoms", "movable_mm_atoms", "frozen_mm_atoms")
    )
    if explicit_layer_lists:
        return False

    movable_cutoff = calc_cfg.get("movable_cutoff")
    if movable_cutoff is not None:
        calc_cfg["hess_cutoff"] = float(movable_cutoff)
        if echo_fn is not None:
            echo_fn(
                "[layer] Using all atoms within movable_cutoff as Hessian targets."
            )
        return True

    if not bool(calc_cfg.get("use_bfactor_layers", True)):
        return False

    calc_cfg["hess_cutoff"] = float("inf")
    if echo_fn is not None:
        echo_fn("[layer] 3-layer mode: using MovableMM atoms as Hessian targets.")
    return True


def _resolve_active_atom_indices(
    calc_cfg: Dict[str, Any],
    n_atoms: int,
    active_dof_mode: str,
) -> Tuple[Optional[set[int]], Dict[str, set[int]]]:
    """Resolve active atom indices for active_dof_mode from calculator layer sets."""
    layer_sets = _collect_layer_atom_sets(calc_cfg)
    mode = str(active_dof_mode).lower()
    if mode == "all":
        return None, layer_sets

    ml_indices = layer_sets["ml"]
    hess_mm_indices = layer_sets["hess_mm"]
    movable_mm_indices = layer_sets["movable_mm"]
    frozen_mm_indices = layer_sets["frozen_mm"]
    # For freq, hess_cutoff=None is normalized to hess_cutoff=inf so every
    # MovableMM atom becomes a Hessian target. Keep both sets active: hess_mm
    # is not frozen; it is the movable subset currently covered by the Hessian.
    partial_mm_indices = set(hess_mm_indices) | set(movable_mm_indices)

    if mode == "ml-only":
        active_indices = set(ml_indices)
    elif mode == "partial":
        active_indices = set(ml_indices) | set(partial_mm_indices)
    elif mode == "unfrozen":
        if ml_indices or hess_mm_indices or movable_mm_indices:
            active_indices = set(ml_indices) | set(hess_mm_indices) | set(movable_mm_indices)
        elif frozen_mm_indices:
            active_indices = set(range(int(n_atoms))) - set(frozen_mm_indices)
        else:
            return None, layer_sets
    else:
        active_indices = set(ml_indices) | set(partial_mm_indices)

    if not active_indices:
        return None, layer_sets
    return active_indices, layer_sets


def _validate_hessian_basis_coverage(
    active_indices: set[int],
    layer_sets: Dict[str, set[int]],
    frozen_indices: set[int],
) -> None:
    """Reject an analysis basis wider than the configured Hessian target."""

    requested = set(active_indices) - set(frozen_indices)
    coverage = set(layer_sets["ml"]) | set(layer_sets["hess_mm"])
    missing = sorted(requested - coverage)
    if missing:
        preview = ", ".join(str(index + 1) for index in missing[:12])
        suffix = " …" if len(missing) > 12 else ""
        raise click.ClickException(
            "The requested frequency-analysis basis is wider than the configured "
            "Hessian target before evaluation; missing 1-based atom indices: "
            f"{preview}{suffix}. Increase hess_cutoff or choose a narrower "
            "--active-dof-mode."
        )


def _write_mode_trj_and_pdb(geom,
                            mode_vec_3N: np.ndarray,
                            out_trj: Path,
                            out_pdb: Path,
                            amplitude_ang: float = 0.8,
                            n_frames: int = 20,
                            comment: str = "mode",
                            ref_pdb: Optional[Path] = None) -> None:
    """Write a single mode animation as _trj.xyz (XYZ-like) and .pdb.

    If `ref_pdb` is provided and is a .pdb file, the .pdb is generated by
    converting the _trj.xyz using the input PDB as the template (same as path_opt).
    """
    ref_ang = geom.coords3d * BOHR2ANG
    mode = mode_vec_3N.reshape(-1, 3).copy()
    mode /= np.linalg.norm(mode)

    # _trj.xyz (concatenated XYZ-like trajectory)
    if ref_pdb is not None and ref_pdb.suffix.lower() == ".pdb":
        # Emit a simple XYZ-like trajectory in Å for the converter
        with out_trj.open("w", encoding="utf-8") as f:
            for i in range(n_frames):
                phase = np.sin(2.0 * np.pi * i / n_frames)
                coords = ref_ang + phase * amplitude_ang * mode  # Å
                f.write(f"{len(geom.atoms)}\n{comment} frame={i+1}/{n_frames}\n")
                for sym, (x, y, z) in zip(geom.atoms, coords):
                    f.write(f"{sym:2s} {x: .8f} {y: .8f} {z: .8f}\n")
        # Generate PDB using the input PDB as template (respects convert-files toggle)
        if is_convert_file_enabled():
            try:
                convert_xyz_to_pdb(out_trj, ref_pdb, out_pdb)
            except Exception as exc:
                # Fallback: generate MODEL/ENDMDL using ASE. Say so — the file
                # still appears, but without the reference topology its atom
                # names and residues are not the input's.
                click.echo(
                    "[convert] WARNING: mode PDB fell back to plain ASE output "
                    f"without the reference topology: {exc}",
                    err=True,
                )
                atoms0 = Atoms(geom.atoms, positions=ref_ang, pbc=False)
                for i in range(n_frames):
                    phase = np.sin(2.0 * np.pi * i / n_frames)
                    ai = atoms0.copy()
                    ai.set_positions(ref_ang + phase * amplitude_ang * mode)
                    write(out_pdb, ai, append=(i != 0))
        return

    # If no ref_pdb is given, use the legacy behavior (use pysisyphus.make_trj_str if available)
    try:
        from pysisyphus.xyzloader import make_trj_str  # type: ignore
        amp_ang = amplitude_ang
        steps = np.sin(2.0 * np.pi * np.arange(n_frames) / n_frames)[:, None, None] * (amp_ang * mode[None, :, :])
        traj_ang = ref_ang[None, :, :] + steps  # (T,N,3) in Å
        traj_bohr = traj_ang.reshape(n_frames, -1, 3) * ANG2BOHR
        comments = [f"{comment}  frame={i+1}/{n_frames}" for i in range(n_frames)]
        trj_str = make_trj_str(geom.atoms, traj_bohr, comments=comments)
        out_trj.write_text(trj_str, encoding="utf-8")
    except Exception:
        with out_trj.open("w", encoding="utf-8") as f:
            for i in range(n_frames):
                phase = np.sin(2.0 * np.pi * i / n_frames)
                coords = ref_ang + phase * amplitude_ang * mode
                f.write(f"{len(geom.atoms)}\n{comment} frame={i+1}/{n_frames}\n")
                for sym, (x, y, z) in zip(geom.atoms, coords):
                    f.write(f"{sym:2s} {x: .8f} {y: .8f} {z: .8f}\n")

    # .pdb (MODEL/ENDMDL via ASE)
    atoms0 = Atoms(geom.atoms, positions=ref_ang, pbc=False)
    for i in range(n_frames):
        phase = np.sin(2.0 * np.pi * i / n_frames)
        ai = atoms0.copy()
        ai.set_positions(ref_ang + phase * amplitude_ang * mode)
        write(out_pdb, ai, append=(i != 0))



# Geometry defaults — shared with opt.py
GEOM_KW: Dict[str, Any] = deepcopy(OPT_GEOM_KW)

# ML/MM calculator defaults — shared with opt.py
CALC_KW: Dict[str, Any] = deepcopy(OPT_CALC_KW)

# FREQ_KW and THERMO_KW are imported from .defaults


def _validated_symmetry_number(value: object) -> Optional[int]:
    """Return an external rotational symmetry number accepted by thermoanalysis."""
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise click.UsageError(
            "thermo.symmetry_number must be an integer greater than or equal to 1."
        )
    return int(value)


def _symmetry_number_source(
    *, config_has_value: bool, override_has_value: bool, resolved_value: object
) -> str:
    """Describe which configuration layer supplied the resolved symmetry number."""
    if resolved_value is None:
        return "auto"
    if override_has_value:
        return "override"
    if config_has_value:
        return "config"
    return "auto"


def _validated_thermo_condition(value: object, *, name: str) -> float:
    """Return a finite, strictly positive thermochemistry state variable."""
    if isinstance(value, bool):
        raise click.UsageError(
            f"thermo.{name} must be a finite number greater than zero."
        )
    try:
        resolved = float(value)
    except (TypeError, ValueError) as exc:
        raise click.UsageError(
            f"thermo.{name} must be a finite number greater than zero."
        ) from exc
    if not np.isfinite(resolved) or resolved <= 0.0:
        raise click.UsageError(
            f"thermo.{name} must be a finite number greater than zero."
        )
    return resolved


class _FrequencyOutputCollisionError(click.UsageError):
    """A frequency output/input collision."""


def _prepare_thermo_output_paths(
    out_dir: Path,
    *,
    protected_inputs: Tuple[Optional[Path], ...] = (),
) -> Tuple[Path, Path]:
    """Invalidate a prior thermochemistry generation before real work starts."""
    thermo_yaml = Path(out_dir) / "thermoanalysis.yaml"
    thermo_yaml_tmp = Path(out_dir) / "thermoanalysis.yaml.tmp"
    reserved = {thermo_yaml.resolve(), thermo_yaml_tmp.resolve()}
    for protected in protected_inputs:
        if protected is not None and Path(protected).resolve() in reserved:
            raise _FrequencyOutputCollisionError(
                f"Configuration input {protected} collides with a reserved "
                f"frequency output path under {out_dir}."
            )
    thermo_yaml.unlink(missing_ok=True)
    thermo_yaml_tmp.unlink(missing_ok=True)
    return thermo_yaml, thermo_yaml_tmp


def _prepare_frequency_output_paths(
    out_dir: Path,
    *,
    protected_inputs: Tuple[Optional[Path], ...] = (),
) -> Tuple[Path, Path]:
    """Invalidate every public artifact owned by one real frequency run."""
    out_dir = Path(out_dir)
    owned = [
        out_dir / "frequencies_cm-1.txt",
        out_dir / "result.json",
        out_dir / "summary.json",
        *out_dir.glob("mode_*cm-1_trj.xyz"),
        *out_dir.glob("mode_*cm-1.pdb"),
    ]
    reserved = {path.resolve() for path in owned}
    for protected in protected_inputs:
        if protected is not None and Path(protected).resolve() in reserved:
            raise _FrequencyOutputCollisionError(
                f"Input {protected} collides with a reserved frequency output "
                f"path under {out_dir}."
            )
    thermo_paths = _prepare_thermo_output_paths(
        out_dir,
        protected_inputs=protected_inputs,
    )
    for path in owned:
        path.unlink(missing_ok=True)
    return thermo_paths



@click.command(
    help="ML/MM vibrational frequency analysis (PHVA-compatible).",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Enzyme complex PDB used by both geom_loader and the ML/MM calculator.",
)
@click.option(
    "--parm",
    "real_parm7",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Amber parm7 topology for the full enzyme complex.",
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
    "--hessian-calc-mode",
    type=click.Choice(["Analytical", "FiniteDifference"], case_sensitive=False),
    default=None,
    help="How the ML backend builds the Hessian (Analytical or FiniteDifference); "
         "overrides calc.hessian_calc_mode from YAML. Default: 'FiniteDifference'. "
         "Runtime and memory depend on the backend and system; compare both "
         "modes on a representative pilot.",
)
@click.option("--max-write", type=int, default=FREQ_KW["max_write"], show_default=True,
              help="Maximum number of modes to export.")
@click.option("--amplitude-ang", type=float, default=FREQ_KW["amplitude_ang"], show_default=True,
              help="Mode animation amplitude (Å).")
@click.option("--n-frames", type=int, default=FREQ_KW["n_frames"], show_default=True,
              help="Frames per vibrational mode animation.")
@click.option(
    "--sort",
    type=click.Choice(["value", "abs"]),
    default=FREQ_KW["sort"],
    show_default=True,
    help="Sort modes by signed value or absolute value.",
)
@click.option("--temperature", type=float, default=THERMO_KW["temperature"], show_default=True,
              help="Temperature (K) for thermochemistry summary.")
@click.option("--pressure", "pressure_atm",
              type=float, default=THERMO_KW["pressure_atm"], show_default=True,
              help="Pressure (atm) for thermochemistry summary.")
@click.option(
    "--dump/--no-dump",
    default=THERMO_KW["dump"],
    show_default=True,
    help="Write 'thermoanalysis.yaml' alongside the console summary.",
)
@click.option("-o", "--out-dir", type=str, default=FREQ_KW["out_dir"], show_default=True, help="Output directory.")
@click.option(
    "--active-dof-mode",
    type=click.Choice(["all", "ml-only", "partial", "unfrozen"], case_sensitive=False),
    default="partial",
    show_default=True,
    help="Active DOF selection for frequency analysis: "
         "all (all atoms), ml-only (ML only), partial (ML + MovableMM, default), "
         "unfrozen (all non-frozen atoms).",
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
    help="Validate options and print the execution plan without running frequency analysis.",
)
@click.option(
    "--ref-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Reference PDB topology to use when --input is XYZ (keeps XYZ coordinates).",
)
@click.option(
    "--convert-files/--no-convert-files",
    "convert_files",
    default=True,
    show_default=True,
    help="Convert XYZ/TRJ outputs into PDB companions based on the input format.",
)
@click.option(
    "--hess-device",
    "hess_device",
    type=click.Choice(["auto", "cuda", "cpu"], case_sensitive=False),
    default="auto",
    show_default=True,
    help="Device for Hessian assembly and diagonalization (auto/cuda/cpu). "
         "Use 'cpu' to avoid VRAM issues with large unfrozen systems. "
         "ML model inference always uses ml_device (typically GPU).",
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
@click.option(
    "--dump-hess",
    "dump_hess",
    type=click.Path(dir_okay=False),
    default=None,
    show_default=False,
    help="Save the computed Hessian and geometry/active-basis identity to a "
         "compressed .npz file for a matching 'mlmm irc --read-hess' run. "
         "The file also identifies model charge and multiplicity.",
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
@add_precision_option()
@add_workers_options()
@add_backend_model_option()
@add_calc_file_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path,
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
    max_write: int,
    amplitude_ang: float,
    n_frames: int,
    sort: str,
    temperature: float,
    pressure_atm: float,
    dump: bool,
    out_dir: str,
    active_dof_mode: str,
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    ref_pdb: Optional[Path],
    convert_files: bool,
    hess_device: str,
    backend: Optional[str],
    embedcharge: bool,
    embedcharge_cutoff: Optional[float],
    link_atom_method: Optional[str],
    mm_backend: Optional[str],
    use_cmap: Optional[bool],
    dump_hess: Optional[str],
    out_json: bool,
    precision: Optional[str],
    workers: Optional[int],
    workers_per_node: Optional[int],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: Optional[str],
) -> None:
    set_convert_file_enabled(convert_files)
    time_start = time.perf_counter()
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

    # Validate input format: PDB/mmCIF directly, or XYZ with --ref-pdb.
    suffix = input_path.suffix.lower()
    if suffix not in (".pdb", ".cif", ".mmcif", ".xyz"):
        click.echo("ERROR: --input must be a PDB, mmCIF, or XYZ file.", err=True)
        sys.exit(1)
    if suffix == ".xyz" and ref_pdb is None:
        click.echo("ERROR: --ref-pdb is required when --input is an XYZ file.", err=True)
        sys.exit(1)

    prepared_input = prepare_input_structure(input_path)
    try:
        apply_ref_pdb_override(prepared_input, ref_pdb)
    except click.BadParameter as e:
        click.echo(f"ERROR: {e}", err=True)
        prepared_input.cleanup()
        sys.exit(1)

    geom_input_path = prepared_input.geom_path
    source_path = prepared_input.source_path  # PDB topology for output conversion
    charge, spin = resolve_charge_spin_or_raise(
        prepared_input, charge, spin,
        ligand_charge=ligand_charge, prefix="[freq]",
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

    try:
        config_layer_cfg = load_yaml_dict(config_yaml)
        override_layer_cfg = load_yaml_dict(override_yaml)
    except ValueError as e:
        click.echo(f"ERROR: {e}", err=True)
        prepared_input.cleanup()
        sys.exit(1)

    geom_cfg = deepcopy(GEOM_KW)
    calc_cfg = deepcopy(CALC_KW)
    freq_cfg = dict(FREQ_KW)
    thermo_cfg = dict(THERMO_KW)
    # Keep the command-level detect-layer default unless YAML or an explicit
    # CLI option overrides it.
    calc_cfg["use_bfactor_layers"] = bool(detect_layer)

    apply_yaml_overrides(
        config_layer_cfg,
        [
            (geom_cfg, (("geom",),)),
            (calc_cfg, (("calc",), ("mlmm",))),
            (freq_cfg, (("freq",),)),
            (thermo_cfg, (("thermo",), ("freq", "thermo"))),
        ],
    )
    thermo_paths = (("thermo",), ("freq", "thermo"))
    _config_has_symmetry_number = (
        yaml_section_has_key(config_layer_cfg, thermo_paths, "symmetry_number")
        and thermo_cfg.get("symmetry_number") is not None
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
    if link_atom_method is not None:
        calc_cfg["link_atom_method"] = str(link_atom_method).lower()
    if mm_backend is not None:
        calc_cfg["mm_backend"] = str(mm_backend).lower()
    if use_cmap is not None:
        calc_cfg["use_cmap"] = use_cmap

    if _is_param_explicit("hessian_calc_mode") and hessian_calc_mode is not None:
        calc_cfg["hessian_calc_mode"] = str(hessian_calc_mode)
    # The first routing pass precedes this explicit CLI override. Re-run the
    # compatibility guard so ``--workers >1 --hessian-calc-mode Analytical``
    # cannot evade validation.
    apply_workers_to_calc_cfg(calc_cfg, None, None)

    if _is_param_explicit("max_write"):
        freq_cfg["max_write"] = int(max_write)
    if _is_param_explicit("amplitude_ang"):
        freq_cfg["amplitude_ang"] = float(amplitude_ang)
    if _is_param_explicit("n_frames"):
        freq_cfg["n_frames"] = int(n_frames)
    if _is_param_explicit("sort"):
        freq_cfg["sort"] = str(sort)
    if _is_param_explicit("out_dir"):
        freq_cfg["out_dir"] = out_dir
    if _is_param_explicit("active_dof_mode"):
        freq_cfg["active_dof_mode"] = str(active_dof_mode)
    if _is_param_explicit("temperature"):
        thermo_cfg["temperature"] = float(temperature)
    if _is_param_explicit("pressure_atm"):
        thermo_cfg["pressure_atm"] = float(pressure_atm)
    if _is_param_explicit("dump"):
        thermo_cfg["dump"] = bool(dump)

    if _is_param_explicit("hess_cutoff") and hess_cutoff is not None:
        calc_cfg["hess_cutoff"] = float(hess_cutoff)
    if _is_param_explicit("movable_cutoff") and movable_cutoff is not None:
        calc_cfg["movable_cutoff"] = float(movable_cutoff)
    if _is_param_explicit("detect_layer"):
        calc_cfg["use_bfactor_layers"] = bool(detect_layer)

    # CLI-resolved charge/spin (from -q / -l derivation in resolve_charge_spin_or_raise,
    # or -m / spin_default) always wins over the CALC_KW default carried in calc_cfg.
    # Same pattern as the opt command: -l 'IN2:-1' must propagate through
    # to model_charge instead of being silently overridden by MLMM_CALC_KW default 0.
    calc_cfg["model_charge"] = int(charge)
    calc_cfg["model_mult"] = int(spin)

    calc_cfg["input_pdb"] = str(source_path)
    calc_cfg["real_parm7"] = str(real_parm7)
    if model_pdb is not None:
        calc_cfg["model_pdb"] = str(model_pdb)

    apply_yaml_overrides(
        override_layer_cfg,
        [
            (geom_cfg, (("geom",),)),
            (calc_cfg, (("calc",), ("mlmm",))),
            (freq_cfg, (("freq",),)),
            (thermo_cfg, (("thermo",), ("freq", "thermo"))),
        ],
    )
    _override_has_symmetry_number = (
        yaml_section_has_key(
            override_layer_cfg, thermo_paths, "symmetry_number"
        )
        and thermo_cfg.get("symmetry_number") is not None
    )
    thermo_cfg["symmetry_number"] = _validated_symmetry_number(
        thermo_cfg.get("symmetry_number")
    )
    symmetry_number_source = _symmetry_number_source(
        config_has_value=_config_has_symmetry_number,
        override_has_value=_override_has_symmetry_number,
        resolved_value=thermo_cfg["symmetry_number"],
    )
    thermo_cfg["temperature"] = _validated_thermo_condition(
        thermo_cfg.get("temperature"), name="temperature"
    )
    thermo_cfg["pressure_atm"] = _validated_thermo_condition(
        thermo_cfg.get("pressure_atm"), name="pressure_atm"
    )
    active_dof_mode_value = normalize_choice(
        str(freq_cfg.get("active_dof_mode", active_dof_mode)),
        param="--active-dof-mode",
        alias_groups=(
            (("all",), "all"),
            (("ml-only",), "ml-only"),
            (("partial",), "partial"),
            (("unfrozen",), "unfrozen"),
        ),
        allowed_hint="all|ml-only|partial|unfrozen",
    )
    freq_cfg["active_dof_mode"] = active_dof_mode_value
    for cutoff_key in ("hess_cutoff", "movable_cutoff"):
        cutoff_value = calc_cfg.get(cutoff_key)
        if cutoff_value is None:
            continue
        try:
            cutoff_float = float(cutoff_value)
        except (TypeError, ValueError) as exc:
            raise click.BadParameter(
                f"{cutoff_key} must be a finite non-negative distance.",
                param_hint=f"--{cutoff_key.replace('_', '-')}",
            ) from exc
        if not np.isfinite(cutoff_float) or cutoff_float < 0.0:
            raise click.BadParameter(
                f"{cutoff_key} must be a finite non-negative distance.",
                param_hint=f"--{cutoff_key.replace('_', '-')}",
            )
        calc_cfg[cutoff_key] = cutoff_float
    from pysisyphus.tr_projection import normalize_tr_projection_mode
    geom_cfg["tr_projection"] = normalize_tr_projection_mode(
        geom_cfg.get("tr_projection")
    )
    calc_paths = (("calc",), ("mlmm",))
    partial_explicit = (
        yaml_section_has_key(config_layer_cfg, calc_paths, "return_partial_hessian")
        or yaml_section_has_key(override_layer_cfg, calc_paths, "return_partial_hessian")
    )
    if not partial_explicit:
        calc_cfg["return_partial_hessian"] = True

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

    out_dir_path = Path(freq_cfg.get("out_dir", FREQ_KW["out_dir"])).resolve()
    layer_source_pdb = source_path
    detect_layer_enabled = bool(calc_cfg.get("use_bfactor_layers", True))
    model_pdb_cfg = calc_cfg.get("model_pdb")
    movable_cutoff_value = calc_cfg.get("movable_cutoff")
    if movable_cutoff_value is not None:
        if detect_layer_enabled:
            click.echo("[layer] movable_cutoff is set; disabling detect-layer mode.", err=True)
        detect_layer_enabled = False
        calc_cfg["use_bfactor_layers"] = False

    from mlmm.core.embedcharge_policy import reject_retired_embedcharge_cli

    reject_retired_embedcharge_cli(
        calc_cfg,
        cutoff_requested=_is_param_explicit("embedcharge_cutoff"),
    )

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
        if detect_layer_enabled and layer_source_pdb.suffix.lower() != ".pdb":
            click.echo("ERROR: --detect-layer requires a PDB input (or --ref-pdb).", err=True)
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
                    "detect_layer": bool(detect_layer_enabled),
                    "model_region_source": model_region_source,
                    "model_indices_count": 0 if not model_indices else len(model_indices),
                    "active_dof_mode": str(freq_cfg.get("active_dof_mode", active_dof_mode)),
                    "tr_projection": geom_cfg["tr_projection"],
                    "will_run_frequency_analysis": True,
                    "will_write_modes": True,
                    "will_dump_thermo_yaml": bool(thermo_cfg.get("dump", False)),
                    "backend": calc_cfg.get("backend", "uma"),
                    "embedcharge": bool(calc_cfg.get("embedcharge", False)),
                },
            )
        )
        click.echo("[dry-run] Validation complete. Frequency execution was skipped.")
        emit(
            format_elapsed("[time] Elapsed Time for Freq", time_start),
            narrative=True,
        )
        return

    out_dir_path.mkdir(parents=True, exist_ok=True)
    # A real invocation owns this exact output generation.  Invalidate both
    # published and staged thermochemistry before layer preparation, calculator,
    # or Hessian work can fail, including under --no-dump.
    frequency_protected_inputs = (
        input_path,
        prepared_input.original_path,
        source_path,
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
    )
    try:
        _thermo_yaml, _thermo_yaml_tmp = _prepare_frequency_output_paths(
            out_dir_path,
            protected_inputs=frequency_protected_inputs,
        )
    except _FrequencyOutputCollisionError:
        prepared_input.cleanup()
        raise

    if detect_layer_enabled and layer_source_pdb.suffix.lower() != ".pdb":
        click.echo("ERROR: --detect-layer requires a PDB input (or --ref-pdb).", err=True)
        prepared_input.cleanup()
        sys.exit(1)

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
            protected_inputs=frequency_protected_inputs,
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

    # Default-verbosity entry summary (skipped in child mode).
    from mlmm.core.utils import calculator_run_label, echo_run_summary
    echo_run_summary({
        "input": str(input_path),
        "backend": calculator_run_label(calc_cfg),
        "out": str(out_dir_path),
    })

    click.echo(pretty_block("geom", format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")))
    echo_calc = format_freeze_atoms_for_echo(filter_calc_for_echo(calc_cfg), key="freeze_atoms")
    click.echo(pretty_block("calc", echo_calc))
    echo_freq = strip_inherited_keys({**freq_cfg, "out_dir": str(out_dir_path)}, FREQ_KW, mode="same")
    click.echo(pretty_block("freq", echo_freq))
    echo_thermo = strip_inherited_keys(thermo_cfg, THERMO_KW, mode="same")
    click.echo(pretty_block("thermo", echo_thermo))

    # freq has no microiteration and computes a cartesian Hessian; internal
    # coordinates (dlc/redund/tric) give no benefit and the result is identical.
    # Fixed to cartesian (also ignores any coord_type injected by `mlmm all`).
    coord_type = "cart"
    coord_kwargs = dict(geom_cfg)
    coord_kwargs.pop("coord_type", None)
    geometry = geom_loader(geom_input_path, coord_type=coord_type, **coord_kwargs)

    masses_amu = np.array([atomic_masses[z] for z in geometry.atomic_numbers])
    # Resolve Hessian assembly/diagonalization device separately from ML inference device.
    # --hess-device=cpu allows large Hessians to be assembled on CPU while ML model stays on GPU.
    if hess_device.lower() == "auto":
        device = _torch_device(calc_cfg.get("ml_device", "auto"))
    else:
        device = _torch_device(hess_device.lower())
    if device.type == "cpu":
        click.echo("[device] Hessian assembly and diagonalization will run on CPU.")
    masses_au_t = torch.as_tensor(masses_amu * AMU2AU, dtype=torch.float32, device=device)

    n_atoms = len(geometry.atoms)
    all_indices = set(range(n_atoms))

    # Determine active atoms based on mode
    active_dof_mode_lower = str(freq_cfg.get("active_dof_mode", active_dof_mode)).lower()
    active_indices, layer_sets = _resolve_active_atom_indices(calc_cfg, n_atoms, active_dof_mode_lower)
    ml_indices = layer_sets["ml"]
    movable_mm_indices = layer_sets["movable_mm"]
    partial_mm_indices = movable_mm_indices  # hess_cutoff is microiter-only; freq/tsopt active DOF stays full-movable

    if active_dof_mode_lower == "all" or active_indices is None:
        active_indices = all_indices
        emit("[active-dof] Using all atoms for frequency analysis.", detail=True)
    elif active_dof_mode_lower == "ml-only":
        emit(f"[active-dof] Using ML atoms only for frequency analysis (n={len(active_indices)}).", detail=True)
    elif active_dof_mode_lower == "partial":
        emit(f"[active-dof] Using ML + MovableMM atoms for frequency analysis (n={len(active_indices)}).", detail=True)
    elif active_dof_mode_lower == "unfrozen":
        emit(f"[active-dof] Using all non-frozen atoms for frequency analysis (n={len(active_indices)}).", detail=True)
    else:
        active_indices = set(ml_indices) | set(partial_mm_indices)
        emit(f"[active-dof] Defaulting to partial mode (n={len(active_indices)}).", detail=True)

    # Atoms not in active_indices become frozen for frequency analysis
    freeze_for_freq = sorted(all_indices - active_indices)
    # Also include any explicitly frozen atoms from config
    explicit_freeze = set(calc_cfg.get("freeze_atoms") or [])
    freeze_list = sorted(set(freeze_for_freq) | explicit_freeze)
    _validate_hessian_basis_coverage(
        set(active_indices),
        layer_sets,
        set(freeze_list),
    )

    try:
        from mlmm.io.hessian_cache import (
            load_matching as _hess_load_matching,
            identity_from_context as _hess_identity,
        )
        # reuse a cached TS Hessian only on a full evaluation-identity
        # match (run/system/evaluator/active space/potential).  The all
        # workflow may round-trip the TS through a three-decimal PDB, so the
        # coordinate field keeps the wider bohr tolerance.
        _cached_ts = _hess_load_matching(
            "ts",
            _hess_identity(geometry, calc_cfg, role="ts"),
            atol=1.1e-3,
        )
        if _cached_ts is not None:
            emit("[freq] Reusing cached TS Hessian.", narrative=True)
            H_t = _cached_ts["hessian"]
            if isinstance(H_t, torch.Tensor):
                H_t = H_t.to(device=device)
            else:
                H_t = torch.as_tensor(H_t, device=device)
            energy_ha = _cached_ts.get("meta", {}).get("energy_ha")
            if energy_ha is None:
                energy_ha = float(geometry.energy)
            # Restore active-DOF metadata so cached partial Hessians follow the
            # same PHVA route as freshly evaluated Hessians.
            _cached_active_dofs = _cached_ts.get("active_dofs")
            if _cached_active_dofs is not None and len(_cached_active_dofs) > 0:
                _active_atoms_from_cache = sorted({d // 3 for d in _cached_active_dofs})
                geometry._hess_active_atoms_last = np.asarray(
                    _active_atoms_from_cache, dtype=int
                )
                geometry._hess_active_dofs_last = np.asarray(
                    _cached_active_dofs, dtype=int
                )
                if getattr(geometry, "within_partial_hessian", None) is None:
                    geometry.within_partial_hessian = {
                        "active_n_dof": len(_cached_active_dofs),
                        "full_n_dof": int(geometry.cart_coords.size),
                        "active_dofs": list(_cached_active_dofs),
                        "active_atoms": _active_atoms_from_cache,
                    }
        else:
            # Populate active-DOF metadata for downstream PHVA routing.
            H_t, energy_ha = _calc_full_hessian_torch(
                geometry, calc_cfg, device, refresh_geom_meta=True
            )

        echo_resolved_device()

        _raw_hessian_shape = tuple(H_t.shape)
        _requested_atoms = [
            i for i in range(len(geometry.atomic_numbers))
            if i not in set(int(j) for j in freeze_list)
        ]
        H_analysis, _analysis_atoms, _computed_atoms, _hessian_storage = (
            _reconcile_hessian_analysis_basis(
                H_t,
                geometry,
                _requested_atoms,
            )
        )
        del H_t

        # --dump-hess: save Hessian to compressed .npz
        if dump_hess:
            _h_np = H_analysis.detach().cpu().numpy()
            # The opt-in dump needs a host copy; H_analysis remains live for
            # the frequency calculation immediately below.
            _dump_path = Path(dump_hess)
            # Persist partial-Hessian metadata so `mlmm irc --read-hess` can
            # restore geometry.within_partial_hessian. Without it the loaded
            # partial Hessian (active_n_dof != 3N) tripped a Geometry shape
            # assertion (the npz was not a usable round-trip).
            from mlmm.io.hessian_file import save_hessian_file
            from mlmm.io.hessian_cache import persistent_identity_from_context

            _analysis_dofs = [
                3 * atom + axis for atom in _analysis_atoms for axis in range(3)
            ]
            _analysis_metadata = None
            if len(_analysis_atoms) != len(geometry.atomic_numbers):
                _analysis_metadata = {
                    "active_n_dof": len(_analysis_dofs),
                    "full_n_dof": int(geometry.cart_coords.size),
                    "active_dofs": _analysis_dofs,
                    "active_atoms": list(_analysis_atoms),
                }
            _dump_path = save_hessian_file(
                _dump_path,
                hessian=_h_np,
                energy_ha=float(energy_ha) if energy_ha is not None else 0.0,
                cart_coords_bohr=geometry.cart_coords,
                atomic_numbers=geometry.atomic_numbers,
                model_charge=int(calc_cfg["model_charge"]),
                model_mult=int(calc_cfg["model_mult"]),
                potential_identity=persistent_identity_from_context(
                    geometry,
                    calc_cfg,
                ),
                partial_metadata=_analysis_metadata,
            )
            emit(f"[freq] Hessian saved → {_dump_path} (shape={_h_np.shape})", narrative=True)
            del _h_np

        coords_bohr = geometry.coords3d
        _effective_freeze = sorted(
            i for i in range(len(geometry.atomic_numbers))
            if i not in set(_analysis_atoms)
        )
        _n_atoms = len(geometry.atomic_numbers)
        _n_frozen = len(_effective_freeze)
        _n_active = len(_analysis_atoms)
        emit(
            f"[freq] Hessian ready: raw_shape={_raw_hessian_shape}, "
            f"analysis_shape={tuple(H_analysis.shape)}, "
            f"active_atoms={_n_active}/{_n_atoms}, frozen_atoms={_n_frozen}, "
            f"computed_atoms={len(_computed_atoms)}, active_dof={3 * _n_active}",
            detail=True,
        )
        _rigid_projection = {}
        freqs_cm, modes_mw = _frequencies_cm_and_modes(
            H_analysis,
            geometry.atomic_numbers,
            coords_bohr,
            device,
            freeze_idx=_effective_freeze if _effective_freeze else None,
            tr_projection=geom_cfg["tr_projection"],
            projection_info=_rigid_projection,
        )
        _rigid_projection.update(
            {
                "hessian_space": (
                    "full" if not _effective_freeze else "active"
                ),
                "hessian_shape": list(H_analysis.shape),
                "raw_hessian_shape": list(_raw_hessian_shape),
                "computed_atom_count": len(_computed_atoms),
                "analysis_atom_count": len(_analysis_atoms),
                "storage": _hessian_storage,
                "hessian_source": "cache" if _cached_ts is not None else "fresh",
                "hessian_representation": "cartesian-unweighted-unprojected",
            }
        )
        click.echo(
            "[freq] Rigid projection: "
            f"treatment={_rigid_projection['treatment']}, "
            f"rank={_rigid_projection['effective_rank']}, "
            f"full_rigid_rank={_rigid_projection['full_rigid_rank']}."
        )

        del H_analysis
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        order = (
            np.argsort(np.abs(freqs_cm))
            if freq_cfg["sort"] == "abs"
            else np.argsort(freqs_cm)
        )
        n_write = int(min(freq_cfg["max_write"], len(order)))
        _imag = [f for f in freqs_cm if f < 0]
        _preview_n = min(20, len(order))
        _freq_preview = ", ".join(f"{float(freqs_cm[j]):+.1f}" for j in order[:_preview_n])
        _suffix = ", ..." if len(order) > _preview_n else ""
        emit(
            f"[freq] {len(freqs_cm)} modes ({len(_imag)} imaginary); "
            f"first {_preview_n} by {freq_cfg['sort']}: [{_freq_preview}{_suffix}] cm⁻¹; "
            f"full list: {out_dir_path / 'frequencies_cm-1.txt'}",
            narrative=True,
        )
        if len(order) > _preview_n:
            _freq_str = ", ".join(f"{float(freqs_cm[j]):+.1f}" for j in order)
            click.echo(f"[freq:all] [{_freq_str}] cm⁻¹")
        emit(f"[INFO] Writing {n_write} mode(s) ({freq_cfg['sort']} ordering).", detail=True)

        ref_pdb_for_modes = source_path if source_path.suffix.lower() == ".pdb" else None
        _mode_output_files: list[str] = []
        for k, idx in enumerate(order[:n_write], start=1):
            freq_val = float(freqs_cm[idx])
            mode_cart_3N = _mw_mode_to_cart(modes_mw[idx], masses_au_t)
            out_trj = out_dir_path / f"mode_{k:04d}_{freq_val:+.2f}cm-1_trj.xyz"
            out_pdb = out_dir_path / f"mode_{k:04d}_{freq_val:+.2f}cm-1.pdb"
            _write_mode_trj_and_pdb(
                geometry,
                mode_cart_3N,
                out_trj,
                out_pdb,
                amplitude_ang=freq_cfg["amplitude_ang"],
                n_frames=freq_cfg["n_frames"],
                comment=f"mode {k}  {freq_val:+.2f} cm-1",
                ref_pdb=ref_pdb_for_modes,
            )
            _mode_output_files.append(out_trj.name)
            if out_pdb.is_file():
                _mode_output_files.append(out_pdb.name)

        (out_dir_path / "frequencies_cm-1.txt").write_text(
            "\n".join(f"{i+1:4d}  {float(freqs_cm[j]):+12.4f}" for i, j in enumerate(order)),
            encoding="utf-8",
        )

        del modes_mw
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        _thermo_data = None
        try:
            from thermoanalysis.QCData import (
                QCData,
                detect_point_group_and_symmetry_number,
            )
            from thermoanalysis.constants import J2AU, J2CAL, NA
            from thermoanalysis.thermo import thermochemistry
            from thermoanalysis.config import WORKFLOW_THERMO_POLICY

            qc_data = {
                "coords3d": geometry.coords3d * BOHR2ANG,
                "wavenumbers": freqs_cm,
                "scf_energy": float(energy_ha),
                "masses": masses_amu,
                "mult": int(calc_cfg["model_mult"]),
            }
            point_group, detected_symmetry_number, point_group_source = (
                detect_point_group_and_symmetry_number(
                    geometry.atomic_numbers,
                    qc_data["coords3d"],
                )
            )
            configured_symmetry_number = thermo_cfg["symmetry_number"]
            symmetry_number = (
                detected_symmetry_number
                if configured_symmetry_number is None
                else int(configured_symmetry_number)
            )
            if configured_symmetry_number is None:
                symmetry_number_source = point_group_source
            qc = QCData(
                qc_data,
                point_group=point_group,
                mult=int(calc_cfg["model_mult"]),
            )
            qc.symmetry_number = symmetry_number

            T = float(thermo_cfg["temperature"])
            p_atm = float(thermo_cfg["pressure_atm"])
            p_pa = p_atm * 101325.0  # Pa

            # The standalone-freq policy is library-default QRRHO with no
            # imaginary inversion and NO positive-frequency floor. Pass every value
            # explicitly and serialize the effective policy below; this reproduces
            # the historical bare thermochemistry(qc, T, pressure=p) numbers.
            _thermo_policy = WORKFLOW_THERMO_POLICY
            tr = thermochemistry(
                qc, T, pressure=p_pa, **_thermo_policy.thermochemistry_kwargs()
            )  # default: QRRHO

            # Converters
            au2CalMol = (1.0 / J2AU) * NA * J2CAL
            to_cal_per_mol = lambda x: float(x) * au2CalMol
            J_per_Kmol_to_cal_per_Kmol = lambda j: float(j) * J2CAL

            # Counts
            n_imag = int(np.sum(freqs_cm < 0.0))

            # Compose summary
            EE = float(tr.U_el)
            ZPE = float(tr.ZPE)
            dE_therm = float(tr.U_therm)               # Thermal correction to Energy (includes ZPE)
            dH_therm = float(tr.H - tr.U_el)           # Thermal correction to Enthalpy (= U_therm + kBT)
            dG_therm = float(tr.dG)                    # Thermal correction to Free Energy (= G - EE)

            sum_EE_ZPE = EE + ZPE
            sum_EE_thermal_E = float(tr.U_tot)         # = EE + U_therm
            sum_EE_thermal_H = float(tr.H)             # = H
            sum_EE_thermal_G = float(tr.G)             # = G

            E_thermal_cal = to_cal_per_mol(tr.U_therm)               # cal/mol
            Cv_cal_per_Kmol = J_per_Kmol_to_cal_per_Kmol(tr.c_tot)   # cal/(mol*K)
            S_cal_per_Kmol  = to_cal_per_mol(tr.S_tot)               # cal/(mol*K)

            # Echo summary (Gaussian-like)
            click.echo("\nThermochemistry Summary")
            click.echo("------------------------")
            click.echo(f"Structure               = {input_path}")
            click.echo(f"Temperature (K)         = {T:.2f}")
            click.echo(f"Pressure    (atm)       = {p_atm:.4f}")
            click.echo(
                f"Molecular point group   = {point_group} "
                f"({point_group_source})"
            )
            click.echo(
                f"Rotational symmetry no. = {symmetry_number:d} "
                f"({symmetry_number_source})"
            )
            if freeze_list:
                emit("[NOTE] Thermochemistry uses active DOF (PHVA) due to frozen atoms.", narrative=True)
            click.echo(f"Number of Imaginary Freq = {n_imag:d}\n")

            def _ha(x): return f"{float(x): .6f} Ha"
            def _cal(x): return f"{float(x): .2f} cal/mol"
            def _calK(x): return f"{float(x): .2f} cal/(mol*K)"

            click.echo(f"Electronic Energy (E)                  = {_ha(EE)}")
            click.echo(f"Zero-point Energy Correction           = {_ha(ZPE)}")
            click.echo(f"Thermal Correction to Energy           = {_ha(dE_therm)}")
            click.echo(f"Thermal Correction to Enthalpy         = {_ha(dH_therm)}")
            click.echo(f"Gibbs Free Energy Correction (G_corr)  = {_ha(dG_therm)}")
            click.echo(f"EE + Zero-point Energy                 = {_ha(sum_EE_ZPE)}")
            click.echo(f"EE + Thermal Energy Correction         = {_ha(sum_EE_thermal_E)}")
            click.echo(f"EE + Thermal Enthalpy Correction       = {_ha(sum_EE_thermal_H)}")
            click.echo(f"Gibbs Free Energy (G = E + G_corr)      = {_ha(sum_EE_thermal_G)}")
            click.echo("")
            click.echo(f"E (Thermal)                            = {_cal(E_thermal_cal)}")
            click.echo(f"Heat Capacity (Cv)                     = {_calK(Cv_cal_per_Kmol)}")
            click.echo(f"Entropy (S)                            = {_calK(S_cal_per_Kmol)}")
            click.echo("")

            # Dump YAML when requested
            if bool(thermo_cfg["dump"]):
                payload = {
                    "structure": str(input_path),
                    "temperature_K": T,
                    "pressure_atm": p_atm,
                    "point_group": point_group,
                    "point_group_source": point_group_source,
                    "symmetry_number": symmetry_number,
                    "symmetry_number_source": symmetry_number_source,
                    "num_imag_freq": n_imag,
                    "n_freeze_atoms": int(_n_frozen),
                    "thermo_policy": _thermo_policy.as_dict(),
                    "rigid_projection": _rigid_projection,
                    "electronic_energy_ha": EE,
                    "zpe_correction_ha": ZPE,
                    "thermal_correction_energy_ha": dE_therm,
                    "thermal_correction_enthalpy_ha": dH_therm,
                    "thermal_correction_free_energy_ha": dG_therm,
                    "sum_EE_and_ZPE_ha": sum_EE_ZPE,
                    "sum_EE_and_thermal_energy_ha": sum_EE_thermal_E,
                    "sum_EE_and_thermal_enthalpy_ha": sum_EE_thermal_H,
                    "sum_EE_and_thermal_free_energy_ha": sum_EE_thermal_G,
                    "E_thermal_cal_per_mol": E_thermal_cal,
                    "Cv_cal_per_mol_K": Cv_cal_per_Kmol,
                    "S_cal_per_mol_K": S_cal_per_Kmol,
                }
                try:
                    with _thermo_yaml_tmp.open("w", encoding="utf-8") as f:
                        yaml.safe_dump(
                            payload, f, sort_keys=False, allow_unicode=True
                        )
                    _thermo_yaml_tmp.replace(_thermo_yaml)
                finally:
                    _thermo_yaml_tmp.unlink(missing_ok=True)
                emit(
                    f"[dump] Wrote thermoanalysis summary → {_thermo_yaml}",
                    detail=True,
                )

            _thermo_data = {
                "thermo_policy": _thermo_policy.as_dict(),
                "temperature_K": T,
                "pressure_atm": p_atm,
                "point_group": point_group,
                "point_group_source": point_group_source,
                "symmetry_number": symmetry_number,
                "symmetry_number_source": symmetry_number_source,
                # The E of the reported "E + G_corr = G" identity, under the same key name
                # thermoanalysis.yaml uses, so a consumer can check the identity from
                # result.json alone.
                "electronic_energy_ha": EE,
                "zpe_ha": ZPE,
                "thermal_correction_energy_ha": dE_therm,
                "thermal_correction_enthalpy_ha": dH_therm,
                "thermal_correction_free_energy_ha": dG_therm,
                "sum_EE_and_ZPE_ha": sum_EE_ZPE,
                "sum_EE_and_thermal_energy_ha": sum_EE_thermal_E,
                "sum_EE_and_thermal_enthalpy_ha": sum_EE_thermal_H,
                "sum_EE_and_thermal_free_energy_ha": sum_EE_thermal_G,
                "E_thermal_cal_per_mol": E_thermal_cal,
                "Cv_cal_per_mol_K": Cv_cal_per_Kmol,
                "S_cal_per_mol_K": S_cal_per_Kmol,
            }

        except ImportError as e:
            raise click.ClickException(
                "Thermochemistry failed because 'thermoanalysis' is unavailable."
            ) from e
        except Exception as e:
            raise click.ClickException(f"Thermochemistry failed: {e}") from e

        # summary.md and key_* outputs are disabled.
        emit(f"[DONE] Wrote modes and list → {out_dir_path}", detail=True)

        if out_json:
            from mlmm.core.utils import calculator_provenance, write_result_json
            _all_freqs = [float(f) for f in freqs_cm]
            _imag_freqs = [f for f in _all_freqs if f < 0.0]
            result_data = {
                "status": "completed",
                "n_modes": len(_all_freqs),
                "n_imaginary": len(_imag_freqs),
                "frequencies_cm": _all_freqs,
                "imaginary_frequencies_cm": _imag_freqs,
                "thermochemistry": _thermo_data,
                "rigid_projection": _rigid_projection,
                **calculator_provenance(calc_cfg),
                "charge": calc_cfg.get("model_charge"),
                "spin": calc_cfg.get("model_mult"),
                "n_atoms": len(geometry.atomic_numbers),
                "n_freeze_atoms": int(_n_frozen),
                "input_file": str(input_path),
                "files": {
                    "frequencies_txt": "frequencies_cm-1.txt",
                    "mode_files": _mode_output_files,
                },
            }
            if _thermo_data is not None and bool(thermo_cfg.get("dump", False)):
                if _thermo_yaml.exists():
                    result_data["files"]["thermoanalysis_yaml"] = "thermoanalysis.yaml"
            if dump_hess:
                _record_hessian_result_path(result_data["files"], _dump_path)
            write_result_json(
                out_dir_path, result_data,
                command="freq",
                elapsed_seconds=time.perf_counter() - time_start,
            )

        emit(
            format_elapsed("[time] Elapsed Time for Freq", time_start),
            narrative=True,
        )

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        render_cli_exception(e, label="frequency analysis", out_dir=out_dir, command="freq", time_start=time_start)
    finally:
        prepared_input.cleanup()
        # Release GPU memory so subsequent pipeline stages don't OOM.
        # `= None` decref's the heavy refs; `del` then removes names from
        # the local frame so torch.nn.Module hooks / closures cannot retain.
        geometry = H_t = H_analysis = modes = modes_mw = None
        del geometry, H_t, H_analysis, modes, modes_mw
        gc.collect()  # break cyclic refs inside torch.nn.Module
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


# Allow `python -m mlmm.freq` direct execution
if __name__ == "__main__":
    cli()
