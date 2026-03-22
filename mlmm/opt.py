# mlmm/opt.py

"""
ML/MM geometry optimization (LBFGS or RFO) with UMA + hessian_ff calculator.

Example:
    mlmm opt -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0

For detailed documentation, see: docs/opt.md
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Set

import ast
import contextlib
import gc
import io
import logging

import sys
import textwrap
import traceback

logger = logging.getLogger(__name__)

import click
import numpy as np
import torch
import time

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import ANG2BOHR, BOHR2ANG, AU2EV
from pysisyphus.TablePrinter import TablePrinter

from .mlmm_calc import mlmm, mlmm_mm_only
from .defaults import (
    GEOM_KW_DEFAULT,
    HESSIAN_DIMER_KW,
    MLMM_CALC_KW,
    OPT_BASE_KW,
    LBFGS_KW,
    RFO_KW,
    OPT_MODE_ALIASES,
    MICROITER_KW,
    BFACTOR_ML,
    BFACTOR_MOVABLE_MM,
    BFACTOR_FROZEN,
)
from .utils import (
    append_xyz_trajectory as _append_xyz_trajectory,
    convert_xyz_to_pdb,
    set_convert_file_enabled,
    is_convert_file_enabled,
    convert_xyz_like_outputs,
    deep_update,
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
    resolve_charge_spin_or_raise,
    parse_indices_string,
    build_model_pdb_from_bfactors,
    build_model_pdb_from_indices,
    update_pdb_bfactors_from_layers,
    normalize_choice,
    yaml_section_has_key,
    is_scan_spec_file,
    parse_dist_freeze_list,
    parse_dist_freeze_spec,
    load_pdb_atom_metadata,
)
from .cli_utils import resolve_yaml_sources, load_merged_yaml_cfg, make_is_param_explicit

EV2AU = 1.0 / AU2EV                 # eV → Hartree
H_EVAA_2_AU = EV2AU / (ANG2BOHR * ANG2BOHR)  # (eV/Å^2) → (Hartree/Bohr^2)

# Flatten-loop constants (sourced from defaults.py)
OPT_FLATTEN_NEG_FREQ_THRESH_CM = HESSIAN_DIMER_KW["neg_freq_thresh_cm"]
OPT_FLATTEN_AMP_ANG = HESSIAN_DIMER_KW["flatten_amp_ang"]
OPT_FLATTEN_MAX_ITER = HESSIAN_DIMER_KW["flatten_max_iter"]


# -----------------------------------------------
# Default settings (imported from defaults.py, aliased for compatibility)
# -----------------------------------------------

GEOM_KW: Dict[str, Any] = dict(GEOM_KW_DEFAULT)
CALC_KW: Dict[str, Any] = dict(MLMM_CALC_KW)

# Note: OPT_BASE_KW, LBFGS_KW, RFO_KW are imported from defaults.py


class HarmonicBiasCalculator:
    """Wrap a base UMA calculator with harmonic distance restraints."""

    def __init__(self, base_calc, k: float = 10.0, pairs: Optional[List[Tuple[int, int, float]]] = None):
        self.base = base_calc
        self.k_evAA = float(k)
        self.k_au_bohr2 = self.k_evAA * H_EVAA_2_AU
        self._pairs: List[Tuple[int, int, float]] = list(pairs or [])

    def set_pairs(self, pairs: List[Tuple[int, int, float]]) -> None:
        self._pairs = [(int(i), int(j), float(t)) for (i, j, t) in pairs]

    def _bias_energy_forces_bohr(self, coords_bohr: np.ndarray) -> Tuple[float, np.ndarray]:
        coords = np.array(coords_bohr, dtype=float).reshape(-1, 3)
        n = coords.shape[0]
        E_bias = 0.0
        F_bias = np.zeros((n, 3), dtype=float)
        k = self.k_au_bohr2
        for (i, j, target_ang) in self._pairs:
            if not (0 <= i < n and 0 <= j < n):
                continue
            rij_vec = coords[i] - coords[j]
            rij = float(np.linalg.norm(rij_vec))
            if rij < 1e-14:
                continue
            target_bohr = float(target_ang) * ANG2BOHR
            diff_bohr = rij - target_bohr
            E_bias += 0.5 * k * diff_bohr * diff_bohr
            u = rij_vec / max(rij, 1e-14)
            Fi = -k * diff_bohr * u
            F_bias[i] += Fi
            F_bias[j] -= Fi
        return E_bias, F_bias.reshape(-1)

    def get_forces(self, elem, coords):
        coords_bohr = np.asarray(coords, dtype=float).reshape(-1, 3)
        base = self.base.get_forces(elem, coords_bohr)
        E0 = float(base["energy"])
        F0 = np.asarray(base["forces"], dtype=float).reshape(-1)
        Ebias, Fbias = self._bias_energy_forces_bohr(coords_bohr)
        return {"energy": E0 + Ebias, "forces": F0 + Fbias}

    def get_energy(self, elem, coords):
        coords_bohr = np.asarray(coords, dtype=float).reshape(-1, 3)
        E0 = float(self.base.get_energy(elem, coords_bohr)["energy"])
        Ebias, _ = self._bias_energy_forces_bohr(coords_bohr)
        return {"energy": E0 + Ebias}

    def get_energy_and_forces(self, elem, coords):
        res = self.get_forces(elem, coords)
        return res["energy"], res["forces"]

    def get_energy_and_gradient(self, elem, coords):
        res = self.get_forces(elem, coords)
        return res["energy"], -np.asarray(res["forces"], dtype=float).reshape(-1)

    def __getattr__(self, name: str):
        return getattr(self.base, name)


def _parse_freeze_atoms(arg: Optional[str]) -> List[int]:
    """Parse comma-separated 1-based indices (e.g., "1,3,5") into 0-based ints."""
    if arg is None:
        return []

    items = [chunk.strip() for chunk in str(arg).split(",")]
    indices: List[int] = []
    for idx, chunk in enumerate(items, start=1):
        if not chunk:
            continue
        try:
            value = int(chunk)
        except ValueError as exc:
            raise click.BadParameter(
                f"Invalid integer in --freeze-atoms entry #{idx}: '{chunk}'"
            ) from exc
        if value <= 0:
            raise click.BadParameter(
                f"--freeze-atoms expects 1-based positive indices; got {value}"
            )
        indices.append(value - 1)
    return sorted(set(indices))


def _normalize_geom_freeze(value: Any) -> List[int]:
    """Normalize YAML-provided geom.freeze_atoms to a sorted 0-based list."""
    if value is None:
        return []
    if isinstance(value, str):
        tokens = [tok.strip() for tok in value.split(",") if tok.strip()]
        try:
            return sorted({int(tok) for tok in tokens})
        except ValueError as exc:
            raise click.BadParameter(
                "geom.freeze_atoms must contain integers (string form)."
            ) from exc
    try:
        return sorted({int(idx) for idx in value})
    except TypeError as exc:
        raise click.BadParameter("geom.freeze_atoms must be iterable of integers.") from exc


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


# -------------------------------
# PDB helpers for B-factor patch
# -------------------------------

def _pdb_keys_from_line(line: str) -> Tuple[Tuple, Tuple]:
    """
    Extract robust keys from a PDB ATOM/HETATM record.

    Returns:
        key_full: (chain, resseq, icode, resname, atomname, altloc)
        key_simple: (chain, resseq, icode, atomname)
    """
    atom_name = line[12:16].strip()
    altloc = line[16:17].strip()
    resname = line[17:20].strip()
    chain = line[21:22].strip()
    resseq_str = line[22:26].strip()
    try:
        resseq = int(resseq_str)
    except ValueError:
        resseq = -10**9  # unlikely sentinel when missing
    icode = line[26:27].strip()
    key_full = (chain, resseq, icode, resname, atom_name, altloc)
    key_simple = (chain, resseq, icode, atom_name)
    return key_full, key_simple


def _collect_ml_atom_keys(model_pdb: Path) -> Tuple[Set[Tuple], Set[Tuple]]:
    """Collect ML-region atom keys from model_pdb."""
    keys_full: Set[Tuple] = set()
    keys_simple: Set[Tuple] = set()
    try:
        with model_pdb.open("r") as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    kf, ks = _pdb_keys_from_line(line)
                    keys_full.add(kf)
                    keys_simple.add(ks)
    except Exception:
        logger.debug("Failed to collect ML atom keys from %s", model_pdb, exc_info=True)
    return keys_full, keys_simple


def _format_with_bfactor(line: str, b: float) -> str:
    """Return PDB line with B-factor field (cols 61-66) set to b (6.2f)."""
    if len(line) < 66:
        line = line.rstrip("\n")
        line = line + " " * max(0, 66 - len(line))
        line = line + "\n"
    bf_str = f"{b:6.2f}"
    # Preserve occupancy (cols 55-60), overwrite tempFactor (61-66).
    new_line = line[:60] + bf_str + line[66:]
    return new_line


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
                f"[annot]   B-factors set in '{final_pdb}' "
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
            click.echo(f"[annot]   B-factors set in '{final_pdb}' (ML=100, frozen=50, both=150).")
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
                        f"[annot]   B-factors set in '{all_opt_pdb}' "
                        f"(ML={BFACTOR_ML:.0f}, MovableMM={BFACTOR_MOVABLE_MM:.0f}, "
                        f"FrozenMM={BFACTOR_FROZEN:.0f})."
                    )
                else:
                    _annotate_b_factors_inplace(
                        all_opt_pdb,
                        model_pdb=model_pdb,
                        freeze_indices_0based=freeze_indices_0based,
                    )
                    click.echo(f"[annot]   B-factors set in '{all_opt_pdb}' (ML=100, frozen=50, both=150).")

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
                        f"[annot]   B-factors set in '{opt_pdb}' "
                        f"(ML={BFACTOR_ML:.0f}, MovableMM={BFACTOR_MOVABLE_MM:.0f}, "
                        f"FrozenMM={BFACTOR_FROZEN:.0f})."
                    )
                else:
                    _annotate_b_factors_inplace(
                        opt_pdb,
                        model_pdb=model_pdb,
                        freeze_indices_0based=freeze_indices_0based,
                    )
                    click.echo(f"[annot]   B-factors set in '{opt_pdb}' (ML=100, frozen=50, both=150).")

            if not wrote_any:
                click.echo(
                    "[convert] WARNING: neither 'optimization_all_trj.xyz' nor 'optimization_trj.xyz' was found; "
                    "skipping trajectory PDB conversion.",
                    err=True,
                )
        except Exception as e:
            click.echo(f"[convert] WARNING: Failed to convert optimization trajectory to PDB: {e}", err=True)


# -----------------------------------------------
# Flatten helpers
# -----------------------------------------------


def _calc_energy(geom, calc_kwargs: dict, calc=None) -> float:
    """Compute energy (Hartree) from ML/MM calculator."""
    owns_calc = calc is None
    if owns_calc:
        kw = dict(calc_kwargs or {})
        kw["out_hess_torch"] = False
        calc = mlmm(**kw)
    result = calc.get_energy(geom.atoms, geom.coords)
    energy = float(result.get("energy", 0.0))
    del result
    if owns_calc:
        del calc
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return energy


def _flatten_all_imag_modes_for_geom(
    geom,
    masses_amu: np.ndarray,
    calc_kwargs: dict,
    freqs_cm: np.ndarray,
    modes: torch.Tensor,
    neg_freq_thresh_cm: float,
    flatten_amp_ang: float,
) -> bool:
    """
    Flatten all imaginary modes for a geometry in a single pass.
    """
    neg_idx_all = np.where(freqs_cm < -abs(neg_freq_thresh_cm))[0]
    if len(neg_idx_all) == 0:
        return False

    order = np.argsort(freqs_cm[neg_idx_all])  # most negative first
    targets = [int(x) for x in neg_idx_all[order]]
    mass_scale = np.sqrt(12.011 / masses_amu)[:, None]
    amp_bohr = float(flatten_amp_ang) / BOHR2ANG
    E_ref = _calc_energy(geom, calc_kwargs)

    m3 = np.repeat(masses_amu, 3).reshape(-1, 3)
    for idx in targets:
        v_mw = modes[idx].detach().cpu().numpy().reshape(-1, 3)
        v_cart = v_mw / np.sqrt(m3)
        v_cart /= np.linalg.norm(v_cart)

        disp = amp_bohr * mass_scale * v_cart
        ref = geom.cart_coords.reshape(-1, 3)

        plus = ref + disp
        minus = ref - disp

        geom.coords = plus.reshape(-1)
        E_plus = _calc_energy(geom, calc_kwargs)

        geom.coords = minus.reshape(-1)
        E_minus = _calc_energy(geom, calc_kwargs)

        use_plus = E_plus <= E_minus
        geom.coords = (plus if use_plus else minus).reshape(-1)
        E_keep = E_plus if use_plus else E_minus
        delta_e = E_keep - E_ref
        click.echo(
            f"[Flatten] mode={idx} freq={freqs_cm[idx]:+.2f} cm^-1 "
            f"E_disp={E_keep:.8f} Ha \u0394E={delta_e:+.8f} Ha"
        )

    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return True


# -----------------------------------------------
# Microiteration optimizer
# -----------------------------------------------


def _run_microiter_opt(
    geometry,
    calc_cfg: Dict[str, Any],
    rfo_cfg: Dict[str, Any],
    lbfgs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    microiter_cfg: Dict[str, Any],
    out_dir_path: Path,
    *,
    dump: bool = False,
) -> None:
    """Run macro/micro alternating optimization (Gaussian 16-style microiteration).

    Macro step: 1 RFO step moving ML atoms + link-atom MM parents (full ONIOM force).
    Micro step: LBFGS relaxing remaining MM atoms with MM-only forces until convergence.
    Link-atom MM parents are included in the macro step to maintain consistency
    of the link atom position across macro/micro boundaries.
    """
    from .freq import _collect_layer_atom_sets

    # Resolve layer atom sets
    layer_sets = _collect_layer_atom_sets(calc_cfg)
    ml_indices = sorted(layer_sets["ml"])
    movable_mm = sorted(layer_sets["movable_mm"] | layer_sets["hess_mm"])
    frozen_mm = sorted(layer_sets["frozen_mm"])

    if not ml_indices:
        click.echo("[microiter] WARNING: No ML atoms found. Falling back to standard optimization.")
        return None

    # Identify link-atom MM parent atoms (boundary atoms that should move
    # with ML atoms in the macro step to keep link atom positions consistent).
    link_mm_parents: set[int] = set()
    temp_calc = mlmm(**dict(calc_cfg))
    calc_core = temp_calc.core if hasattr(temp_calc, "core") else temp_calc
    for _ml_1, mm_1 in getattr(calc_core, "mlmm_links", []):
        link_mm_parents.add(mm_1 - 1)  # mlmm_links uses 1-based indices
    del temp_calc
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    n_atoms = len(geometry.atoms)
    all_indices = list(range(n_atoms))
    mm_indices = sorted(set(all_indices) - set(ml_indices))

    # Macro step: optimize ML atoms + link MM parents; freeze remaining MM.
    # Micro step: optimize movable MM except link MM parents; freeze ML + link MM parents.
    macro_active = set(ml_indices) | link_mm_parents
    macro_freeze = sorted(set(all_indices) - macro_active)
    micro_freeze = sorted(set(ml_indices) | link_mm_parents | set(frozen_mm))

    max_cycles = int(opt_cfg.get("max_cycles", 10000))
    thresh = opt_cfg.get("thresh", "gau")
    micro_thresh = microiter_cfg.get("micro_thresh") or thresh
    micro_max_cycles = int(microiter_cfg.get("micro_max_cycles", 10000))

    click.echo(
        f"[microiter] ML atoms: {len(ml_indices)}, "
        f"Link MM parents: {len(link_mm_parents)}, "
        f"Movable MM atoms: {len(movable_mm)}, "
        f"Frozen MM atoms: {len(frozen_mm)}"
    )
    click.echo(f"[microiter] Macro thresh: {thresh}, Micro thresh: {micro_thresh}")

    # Create ONIOM calculator (shared core for MM-only calc)
    base_calc = mlmm(**calc_cfg)
    mm_calc = mlmm_mm_only(base_calc.core, freeze_atoms=micro_freeze)

    # Seed initial Hessian for RFO (with macro freeze)
    # Try IRC endpoint cache first; fall back to full Hessian calculation.
    from .hessian_cache import load as _hess_load
    from .freq import (
        _calc_full_hessian_torch as _freq_calc_full_hessian_torch,
        _torch_device as _freq_torch_device,
    )
    hess_device = _freq_torch_device(calc_cfg.get("ml_device", "auto"))

    # Always create macro calculator (needed for optimization loop below)
    macro_calc_cfg = dict(calc_cfg)
    macro_calc_cfg["freeze_atoms"] = macro_freeze
    macro_calc_cfg["hess_mm_atoms"] = sorted(link_mm_parents)  # ML + link MM parents in Hessian
    macro_calc = mlmm(**macro_calc_cfg)

    cached = _hess_load("irc_endpoint")
    _cache_used = False
    if cached is not None:
        active_dofs = cached.get("active_dofs")
        h_raw = cached["hessian"]
        if isinstance(h_raw, torch.Tensor):
            h_init = h_raw.clone()
        else:
            h_init = torch.as_tensor(h_raw, dtype=torch.float64)

        # Macro step freezes MM atoms → only ML DOFs are free.
        # The cached IRC Hessian covers ML+MovableMM DOFs and is generally
        # larger.  Extract the ML-only sub-block when active_dofs are known;
        # otherwise fall back to a fresh Hessian calculation.
        n_free = geometry.cart_coords.size - 3 * len(macro_freeze)
        if h_init.shape[0] == n_free:
            # Size already matches (e.g. non-microiter or same freeze set)
            geometry.set_calculator(macro_calc)
            if active_dofs is not None:
                geometry.within_partial_hessian = {
                    "active_n_dof": len(active_dofs),
                    "full_n_dof": geometry.cart_coords.size,
                    "active_dofs": active_dofs,
                    "active_atoms": sorted(set(d // 3 for d in active_dofs)),
                }
            geometry.cart_hessian = h_init
            click.echo(f"[microiter] Reusing IRC endpoint Hessian for RFO macro step (shape={h_init.shape[0]}x{h_init.shape[1]}).")
            _cache_used = True
        elif active_dofs is not None:
            # Extract ML-only sub-block from the larger cached Hessian.
            macro_free_atoms = sorted(set(range(geometry.cart_coords.size // 3)) - set(macro_freeze))
            macro_free_dofs = []
            for a in macro_free_atoms:
                macro_free_dofs.extend([3 * a, 3 * a + 1, 3 * a + 2])
            # Map macro free DOFs to indices within the cached active_dofs
            cached_dof_set = set(active_dofs)
            sub_indices = []
            for d in macro_free_dofs:
                if d in cached_dof_set:
                    sub_indices.append(active_dofs.index(d))
            if len(sub_indices) == n_free:
                idx = torch.tensor(sub_indices, dtype=torch.long)
                h_sub = h_init[idx][:, idx]
                macro_active_dofs = macro_free_dofs
                geometry.set_calculator(macro_calc)
                geometry.within_partial_hessian = {
                    "active_n_dof": len(macro_active_dofs),
                    "full_n_dof": geometry.cart_coords.size,
                    "active_dofs": macro_active_dofs,
                    "active_atoms": macro_free_atoms,
                }
                geometry.cart_hessian = h_sub
                click.echo(
                    f"[microiter] Reusing IRC endpoint Hessian (sub-block) for RFO macro step "
                    f"(cached {h_init.shape[0]}x{h_init.shape[1]} → extracted {h_sub.shape[0]}x{h_sub.shape[1]})."
                )
                _cache_used = True
                del h_sub
            else:
                click.echo(
                    f"[microiter] IRC endpoint Hessian sub-block extraction failed "
                    f"(expected {n_free}, got {len(sub_indices)}). Falling back to fresh Hessian."
                )
        else:
            click.echo(
                f"[microiter] IRC endpoint Hessian size mismatch "
                f"(cached={h_init.shape[0]}, needed={n_free}). Falling back to fresh Hessian."
            )
        del h_init
    if not _cache_used:
        click.echo("[microiter] Seeding initial Hessian for RFO macro step.")
        geometry.set_calculator(macro_calc)

        h_init, _ = _freq_calc_full_hessian_torch(
            geometry, macro_calc_cfg, hess_device, refresh_geom_meta=True,
        )
        geometry.cart_hessian = h_init
        click.echo(f"[microiter] Initial Hessian seeded (shape={h_init.shape[0]}x{h_init.shape[1]}).")
        del h_init

    optim_all_path = out_dir_path / "optimization_all_trj.xyz"
    macro_trj_path = out_dir_path / "optimization_trj.xyz"
    total_macro_steps = 0

    # Create persistent RFOptimizer once (LayerOpt pattern).
    # This preserves the BFGS Hessian update chain across macro iterations.
    # NOTE: geometry already has macro_calc set (line above); do NOT call
    # set_calculator() again as it clears the pre-computed cart_hessian.
    geometry.freeze_atoms = macro_freeze

    rfo_args = dict(rfo_cfg)
    rfo_args["max_cycles"] = max_cycles
    rfo_args["out_dir"] = str(out_dir_path)
    rfo_args["dump"] = False  # trajectory dumping handled externally
    rfo_args["thresh"] = thresh

    macro_optimizer = RFOptimizer(geometry, **rfo_args)
    macro_optimizer.prepare_opt()  # initialise Hessian from geometry.cart_hessian

    # Microiteration progress table (pysisyphus-style with micro_steps column)
    micro_header = "cycle Δ(energy) max(|force|) rms(force) max(|step|) rms(step) micro_steps s/cycle".split()
    micro_col_fmts = "int float float float float float int float_short".split()
    micro_table = TablePrinter(micro_header, micro_col_fmts, width=12)
    micro_table.print_header()

    for macro_iter in range(max_cycles):
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

        if dump:
            with open(macro_trj_path, "a") as f:
                f.write(geometry.as_xyz() + "\n")
            _append_xyz_trajectory(optim_all_path, macro_trj_path)

        if macro_converged:
            # Print final converged row (no micro steps)
            energy_diff = macro_optimizer.energies[-1] - macro_optimizer.energies[-2] if len(macro_optimizer.energies) >= 2 else float("nan")
            marks = [False, *conv_info.get_convergence()[:-1], False, False]
            cycle_time = time.time() - t_start
            micro_table.print_row(
                (macro_iter, energy_diff, macro_optimizer.max_forces[-1], macro_optimizer.rms_forces[-1],
                 macro_optimizer.max_steps[-1], macro_optimizer.rms_steps[-1], 0, cycle_time),
                marks=marks,
            )
            click.echo(f"[microiter] Macro convergence reached at iteration {macro_iter + 1}.")
            break

        # Apply step to geometry
        new_coords = geometry.coords.copy() + step
        geometry.coords = new_coords
        # Record actual step (may differ due to coordinate back-transformation)
        macro_optimizer.steps[-1] = geometry.coords - macro_optimizer.coords[-1]

        # ---- Micro step: LBFGS with MM-only forces, ML frozen ----
        geometry.freeze_atoms = micro_freeze
        geometry.set_calculator(mm_calc)

        micro_lbfgs_args = dict(lbfgs_cfg)
        micro_lbfgs_args["max_cycles"] = micro_max_cycles
        micro_lbfgs_args["thresh"] = micro_thresh
        micro_lbfgs_args["out_dir"] = str(out_dir_path)
        micro_lbfgs_args["dump"] = dump

        micro_opt = LBFGS(geometry, **micro_lbfgs_args)
        with contextlib.redirect_stdout(io.StringIO()):
            micro_opt.run()
        micro_steps = max(int(micro_opt.cur_cycle) + 1, 1)

        if dump:
            _append_xyz_trajectory(optim_all_path, out_dir_path / "optimization_trj.xyz")

        del micro_opt
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        # Print progress row with micro_steps
        cycle_time = time.time() - t_start
        energy_diff = macro_optimizer.energies[-1] - macro_optimizer.energies[-2] if len(macro_optimizer.energies) >= 2 else float("nan")
        marks = [False, *conv_info.get_convergence()[:-1], False, False]
        if (macro_iter > 1) and (macro_iter % 10 == 0):
            micro_table.print_sep()
        micro_table.print_row(
            (macro_iter, energy_diff, macro_optimizer.max_forces[-1], macro_optimizer.rms_forces[-1],
             macro_optimizer.max_steps[-1], macro_optimizer.rms_steps[-1], micro_steps, cycle_time),
            marks=marks,
        )

    else:
        click.echo(f"[microiter] Reached max macro iterations ({max_cycles}).")

    del macro_optimizer
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    click.echo(f"[microiter] Total macro steps: {total_macro_steps}")
    # Restore full calculator
    geometry.freeze_atoms = list(set(frozen_mm))
    geometry.set_calculator(base_calc)

    return geometry


# -----------------------------------------------
# CLI
# -----------------------------------------------

@click.command(
    help="ML/MM geometry optimization with LBFGS (light) or RFO (heavy).",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (PDB, XYZ). XYZ provides higher coordinate precision. "
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
    help="PDB defining atoms that belong to the ML (high-level) region. "
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
    help="Detect ML/MM layers from input PDB B-factors (ML=0, MovableMM=10, FrozenMM=20). "
         "If disabled, you must provide --model-pdb or --model-indices.",
)
@click.option("-q", "--charge", type=int, required=False,
              help="ML region charge. Required unless --ligand-charge is provided.")
@click.option("-l", "--ligand-charge", type=str, default=None, show_default=False,
              help="Total charge or per-resname mapping (e.g., GPP:-3,SAM:1) used to derive "
                   "charge when -q is omitted (requires PDB input or --ref-pdb).")
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=None,
    show_default=False,
    help="Spin multiplicity (2S+1) for the ML region. Defaults to 1 when omitted.",
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
    "--radius-partial-hessian",
    "--hess-cutoff",
    "radius_partial_hessian",
    type=float,
    default=None,
    show_default=False,
    help="Distance cutoff (Å) from ML region for MM atoms to include in Hessian calculation. "
         "Applied to movable MM atoms and can be combined with --detect-layer. "
         "`--hess-cutoff` is a compatibility alias.",
)
@click.option(
    "--radius-freeze",
    "--movable-cutoff",
    "radius_freeze",
    type=float,
    default=None,
    show_default=False,
    help="Distance cutoff (Å) from ML region for movable MM atoms. "
     "MM atoms beyond this are frozen. "
         "Providing --radius-freeze disables --detect-layer and uses distance-based layer assignment. "
         "`--movable-cutoff` is a compatibility alias.",
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
    help="Interpret --dist-freeze indices as 1-based (default) or 0-based.",
)
@click.option(
    "--bias-k",
    type=float,
    default=300.0,
    show_default=True,
    help="Harmonic restraint strength k [eV/Å^2] for --dist-freeze.",
)
@click.option("--max-cycles", type=int, default=10000, show_default=True, help="Maximum number of optimization cycles.")
@click.option(
    "--dump/--no-dump",
    default=False,
    show_default=True,
    help="Write optimization trajectories ('optimization_trj.xyz' and 'optimization_all_trj.xyz').",
)
@click.option("-o", "--out-dir", type=str, default="./result_opt/", show_default=True, help="Output directory.")
@click.option(
    "--thresh",
    type=click.Choice(["gau_loose", "gau", "gau_tight", "gau_vtight", "baker", "never"], case_sensitive=False),
    default=None,
    help="Convergence preset.",
)
@click.option(
    "--opt-mode",
    type=click.Choice(["grad", "hess", "light", "heavy", "lbfgs", "rfo"], case_sensitive=False),
    default="grad",
    show_default=True,
    help="Optimization mode: grad (lbfgs) or hess (rfo). Aliases light/heavy and lbfgs/rfo are accepted.",
)
@click.option(
    "--microiter/--no-microiter",
    "microiter",
    default=True,
    show_default=True,
    help="Enable microiteration: alternate ML 1-step (RFO) and MM relaxation (LBFGS with MM-only forces). "
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
    show_default=False,
    help="ML backend for the ONIOM high-level region (default: uma).",
)
@click.option(
    "--embedcharge/--no-embedcharge",
    "embedcharge",
    default=False,
    show_default=True,
    help="Enable xTB point-charge embedding correction for MM→ML environmental effects.",
)
@click.option(
    "--embedcharge-cutoff",
    "embedcharge_cutoff",
    type=float,
    default=None,
    show_default=False,
    help="Distance cutoff (Å) from ML region for MM point charges in xTB embedding. "
         "Default: 12.0 Å. Only used when --embedcharge is enabled.",
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
    help="MM backend: hessian_ff (analytical Hessian, default) or openmm (finite-difference Hessian, slower).",
)
@click.option(
    "--cmap/--no-cmap",
    "use_cmap",
    default=None,
    show_default=False,
    help="Enable CMAP (backbone cross-map) terms in model parm7. Default: disabled (Gaussian ONIOM-compatible).",
)
@click.pass_context
def cli(
    ctx: click.Context,
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
    radius_partial_hessian: Optional[float],
    radius_freeze: Optional[float],
    dist_freeze_raw: Sequence[str],
    one_based: bool,
    bias_k: float,
    max_cycles: int,
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    opt_mode: str,
    microiter: bool,
    flatten: bool,
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
) -> None:
    set_convert_file_enabled(convert_files)
    time_start = time.perf_counter()
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

    # Handle input: PDB directly, or XYZ with --ref-pdb for topology
    suffix = input_path.suffix.lower()
    if suffix == ".pdb":
        # PDB input: use directly
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
        click.echo(f"ERROR: Unsupported input format: {suffix}. Use .pdb or .xyz (with --ref-pdb).", err=True)
        sys.exit(1)

    geom_input_path = prepared_input.geom_path
    charge, spin = resolve_charge_spin_or_raise(
        prepared_input, charge, spin,
        ligand_charge=ligand_charge, prefix="[opt]",
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

        apply_yaml_overrides(
            config_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",), ("mlmm",))),
                (opt_cfg, (("opt",),)),
                (lbfgs_cfg, (("lbfgs",), ("opt", "lbfgs"))),
                (rfo_cfg, (("rfo",), ("opt", "rfo"))),
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

        if _is_param_explicit("detect_layer"):
            calc_cfg["use_bfactor_layers"] = bool(detect_layer)
        if _is_param_explicit("radius_partial_hessian") and radius_partial_hessian is not None:
            calc_cfg["hess_cutoff"] = float(radius_partial_hessian)
        if _is_param_explicit("radius_freeze") and radius_freeze is not None:
            calc_cfg["movable_cutoff"] = float(radius_freeze)
            calc_cfg["use_bfactor_layers"] = False

        model_charge_value = calc_cfg.get("model_charge", charge)
        if model_charge_value is None:
            model_charge_value = charge
        calc_cfg["model_charge"] = int(model_charge_value)
        if _is_param_explicit("charge"):
            calc_cfg["model_charge"] = int(charge)
        model_mult_value = calc_cfg.get("model_mult", spin)
        if model_mult_value is None:
            model_mult_value = spin
        calc_cfg["model_mult"] = int(model_mult_value)
        if _is_param_explicit("spin"):
            calc_cfg["model_mult"] = int(spin)
        if model_pdb is not None:
            calc_cfg["model_pdb"] = str(model_pdb)
        calc_cfg["input_pdb"] = str(prepared_input.source_path)
        calc_cfg["real_parm7"] = str(real_parm7)
        if backend is not None:
            calc_cfg["backend"] = str(backend).lower()
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
            ],
        )
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
        if freeze_atoms_cli:
            merge_freeze_atom_indices(geom_cfg, freeze_atoms_cli)
        freeze_atoms_final = list(geom_cfg.get("freeze_atoms") or [])
        calc_cfg["freeze_atoms"] = freeze_atoms_final

        out_dir_path = Path(opt_cfg["out_dir"]).resolve()

        # radius_freeze implies full distance-based layer assignment.
        # radius_partial_hessian alone can be combined with --detect-layer.
        detect_layer_enabled = bool(calc_cfg.get("use_bfactor_layers", True))
        model_pdb_cfg = calc_cfg.get("model_pdb")
        if radius_freeze is not None:
            if detect_layer_enabled:
                click.echo("[layer] --radius-freeze provided; disabling --detect-layer.", err=True)
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
                )
            )

        if dry_run:
            model_region_source = "bfactor"
            if not detect_layer_enabled:
                if model_pdb_cfg is not None:
                    model_region_source = "model_pdb"
                elif model_indices:
                    model_region_source = "model_indices"
                else:
                    click.echo("ERROR: Provide --model-pdb or --model-indices when --no-detect-layer.", err=True)
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
                        "will_run_optimization": True,
                        "will_convert_outputs": True,
                        "backend": calc_cfg.get("backend", "uma"),
                        "embedcharge": bool(calc_cfg.get("embedcharge", False)),
                    },
                )
            )
            click.echo("[dry-run] Validation complete. Optimization execution was skipped.")
            return

        model_pdb_path: Optional[Path] = None
        layer_info: Optional[Dict[str, List[int]]] = None

        if detect_layer_enabled:
            try:
                model_pdb_path, layer_info = build_model_pdb_from_bfactors(layer_source_pdb, out_dir_path)
                calc_cfg["use_bfactor_layers"] = True
                click.echo(
                    f"[layer] Detected B-factor layers: ML={len(layer_info.get('ml_indices', []))}, "
                    f"MovableMM={len(layer_info.get('movable_mm_indices', []))}, "
                    f"FrozenMM={len(layer_info.get('frozen_indices', []))}"
                )
            except Exception as e:
                if model_pdb_cfg is None and not model_indices:
                    click.echo(f"ERROR: {e}", err=True)
                    prepared_input.cleanup()
                    sys.exit(1)
                click.echo(f"[layer] WARNING: {e} Falling back to explicit ML region.", err=True)
                detect_layer_enabled = False

        if not detect_layer_enabled:
            if model_pdb_cfg is None and not model_indices:
                click.echo("ERROR: Provide --model-pdb or --model-indices when --no-detect-layer.", err=True)
                prepared_input.cleanup()
                sys.exit(1)
            if model_pdb_cfg is not None:
                model_pdb_path = Path(model_pdb_cfg)
            else:
                if layer_source_pdb.suffix.lower() != ".pdb":
                    click.echo("ERROR: --model-indices requires a PDB input (or --ref-pdb).", err=True)
                    prepared_input.cleanup()
                    sys.exit(1)
                try:
                    model_pdb_path = build_model_pdb_from_indices(layer_source_pdb, out_dir_path, model_indices or [])
                except Exception as e:
                    click.echo(f"ERROR: {e}", err=True)
                    prepared_input.cleanup()
                    sys.exit(1)
            calc_cfg["use_bfactor_layers"] = False

        if model_pdb_path is None:
            click.echo("ERROR: Failed to resolve model PDB for the ML region.", err=True)
            prepared_input.cleanup()
            sys.exit(1)

        calc_cfg["model_pdb"] = str(model_pdb_path)

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
        from .freq import _align_three_layer_hessian_targets as _freq_align_three_layer_hessian_targets
        _freq_align_three_layer_hessian_targets(calc_cfg, echo_fn=click.echo)

        for key in ("input_pdb", "real_parm7", "model_pdb", "mm_fd_dir"):
            val = calc_cfg.get(key)
            if val:
                calc_cfg[key] = str(Path(val).expanduser().resolve())

        mode_str = "RFO (hess)" if use_rfo else "LBFGS (grad)"
        click.echo(f"\n[mode] Optimizer: {mode_str}\n")
        click.echo(pretty_block("geom", format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")))
        echo_calc = format_freeze_atoms_for_echo(filter_calc_for_echo(calc_cfg), key="freeze_atoms")
        click.echo(pretty_block("calc", echo_calc, defaults=MLMM_CALC_KW))
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
        if dist_freeze:
            display_pairs = []
            for (i, j, target) in dist_freeze:
                label = (f"{target:.4f}" if target is not None else "<current>")
                display_pairs.append((int(i) + 1, int(j) + 1, label))
            click.echo(
                pretty_block(
                    "dist_freeze (input)",
                    {
                        "k (eV/Å^2)": float(bias_k),
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

        base_calc = mlmm(**calc_cfg)
        geometry.set_calculator(base_calc)

        resolved_dist_freeze: List[Tuple[int, int, float]] = []
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
                        "k (eV/Å^2)": float(bias_k),
                        "pairs_1based": [
                            (int(i) + 1, int(j) + 1, float(f"{t:.4f}"))
                            for (i, j, t) in resolved_dist_freeze
                        ],
                    },
                )
            )
            bias_calc = HarmonicBiasCalculator(base_calc, k=float(bias_k))
            bias_calc.set_pairs(resolved_dist_freeze)
            geometry.set_calculator(bias_calc)

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

        def _seed_rfo_hessian():
            """Seed initial Hessian via shared freq backend for RFO."""
            from .hessian_cache import load as _hess_load
            cached = _hess_load("irc_endpoint")
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
                click.echo(f"[opt] Initial Hessian seeded (shape={h_init.shape[0]}x{h_init.shape[1]}).")
                del h_init
                return
            click.echo("[opt] Seeding initial Hessian via shared freq backend.")
            from .freq import (
                _calc_full_hessian_torch as _freq_calc_full_hessian_torch,
                _torch_device as _freq_torch_device,
            )
            hess_device = _freq_torch_device(calc_cfg.get("ml_device", "auto"))
            h_init, _ = _freq_calc_full_hessian_torch(
                geometry,
                calc_cfg,
                hess_device,
                refresh_geom_meta=True,
            )
            geometry.cart_hessian = h_init
            click.echo(f"[opt] Initial Hessian seeded (shape={h_init.shape[0]}x{h_init.shape[1]}).")
            del h_init

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

        use_microiter = bool(microiter) and use_rfo and not dist_freeze
        if bool(microiter) and not use_rfo:
            click.echo("[microiter] --microiter is only effective with --opt-mode hess (RFO). Ignoring.")
        if bool(microiter) and use_rfo and dist_freeze:
            click.echo("[microiter] --microiter is not compatible with --dist-freeze. Falling back to standard RFO.")

        if use_microiter:
            click.echo("\n====== Optimization (RFO + Microiteration) started ======\n")
            _run_microiter_opt(
                geometry,
                calc_cfg,
                rfo_cfg,
                lbfgs_cfg,
                opt_cfg,
                microiter_cfg,
                out_dir_path,
                dump=bool(opt_cfg["dump"]),
            )
            click.echo("\n====== Optimization (RFO + Microiteration) finished ======\n")

            # Write final geometry
            from ase import Atoms as _Atoms
            from ase.io import write as _write
            final_xyz_path = out_dir_path / "final_geometry.xyz"
            final_coords_ang = geometry.coords.reshape(-1, 3) * BOHR2ANG
            atoms_final = _Atoms(geometry.atoms, positions=final_coords_ang, pbc=False)
            _write(final_xyz_path, atoms_final)

        else:
            main_kind = "rfo" if use_rfo else "lbfgs"
            if use_rfo:
                _seed_rfo_hessian()

            main_label = "RFO" if use_rfo else "LBFGS"
            optimizer = _build_optimizer(main_kind)
            click.echo(f"\n====== Optimization ({main_label}) started ======\n")
            optimizer.run()
            click.echo(f"\n====== Optimization ({main_label}) finished ======\n")

            # Get final geometry path
            final_xyz_path = optimizer.final_fn if isinstance(optimizer.final_fn, Path) else Path(optimizer.final_fn)

            if bool(opt_cfg["dump"]):
                optim_all_path = out_dir_path / "optimization_all_trj.xyz"
                if not optim_all_path.exists():
                    trj_path = optimizer.get_path_for_fn("optimization_trj.xyz")
                    _append_xyz_trajectory(optim_all_path, trj_path, reset=True)

        # --------------------------
        # Flatten loop (all imaginary modes)
        # --------------------------
        if flatten:
            from .freq import (
                _torch_device,
                _calc_full_hessian_torch,
                _frequencies_cm_and_modes,
                _safe_masses_amu,
            )

            click.echo("\n====== Optimization (Flatten loop) started ======\n")

            geometry.set_calculator(None)
            uma_kwargs_for_flatten = dict(calc_cfg)
            uma_kwargs_for_flatten["out_hess_torch"] = True
            device = _torch_device(calc_cfg.get("ml_device", "auto"))
            freeze_idx = list(geom_cfg.get("freeze_atoms", [])) if len(geom_cfg.get("freeze_atoms", [])) > 0 else None
            masses_amu = _safe_masses_amu(geometry.atomic_numbers)

            def _attach_opt_calc() -> None:
                geometry.set_calculator(
                    bias_calc if resolved_dist_freeze else base_calc
                )

            def _calc_freqs_and_modes() -> Tuple[np.ndarray, torch.Tensor]:
                H, _e = _calc_full_hessian_torch(geometry, uma_kwargs_for_flatten, device)
                freqs_local, modes_local = _frequencies_cm_and_modes(
                    H,
                    geometry.atomic_numbers,
                    geometry.cart_coords.reshape(-1, 3),
                    device,
                    freeze_idx=freeze_idx,
                )
                del H
                return freqs_local, modes_local

            freqs_cm, modes = _calc_freqs_and_modes()
            neg_mask = freqs_cm < -abs(OPT_FLATTEN_NEG_FREQ_THRESH_CM)
            n_imag = int(np.sum(neg_mask))
            ims = [float(x) for x in freqs_cm if x < -abs(OPT_FLATTEN_NEG_FREQ_THRESH_CM)]
            click.echo(f"[Imaginary modes] n={n_imag}  ({ims})")

            flatten_kind = mode_resolved  # reuse same optimizer type
            for it in range(OPT_FLATTEN_MAX_ITER):
                if n_imag == 0:
                    break
                click.echo(f"[flatten] iteration {it + 1}/{OPT_FLATTEN_MAX_ITER}")
                did_flatten = _flatten_all_imag_modes_for_geom(
                    geometry,
                    masses_amu,
                    uma_kwargs_for_flatten,
                    freqs_cm,
                    modes,
                    OPT_FLATTEN_NEG_FREQ_THRESH_CM,
                    OPT_FLATTEN_AMP_ANG,
                )
                if not did_flatten:
                    click.echo("[flatten] No eligible imaginary modes to flatten; stopping.")
                    break

                _attach_opt_calc()
                opt_restart = _build_optimizer(flatten_kind)
                restart_label = "LBFGS" if flatten_kind == "lbfgs" else "RFO"
                click.echo(f"\n====== Optimization ({restart_label}) restarted ======\n")
                opt_restart.run()
                click.echo(f"\n====== Optimization ({restart_label}) finished ======\n")

                geometry.set_calculator(None)
                freqs_cm, modes = _calc_freqs_and_modes()
                neg_mask = freqs_cm < -abs(OPT_FLATTEN_NEG_FREQ_THRESH_CM)
                n_imag = int(np.sum(neg_mask))
                ims = [float(x) for x in freqs_cm if x < -abs(OPT_FLATTEN_NEG_FREQ_THRESH_CM)]
                click.echo(f"[Imaginary modes] n={n_imag}  ({ims})")

            if n_imag > 0:
                click.echo(
                    f"[flatten] WARNING: Remaining imaginary modes after {OPT_FLATTEN_MAX_ITER} iterations: {n_imag}",
                    err=True,
                )
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            click.echo("\n====== Optimization (Flatten loop) finished ======\n")

            # Update final geometry after flatten
            final_xyz_path = out_dir_path / "final_geometry.xyz"
            from ase import Atoms as _Atoms
            from ase.io import write as _write
            final_coords_ang = geometry.coords.reshape(-1, 3) * BOHR2ANG
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
            get_trj_fn=(lambda fn: out_dir_path / fn) if use_microiter else optimizer.get_path_for_fn,
            final_xyz_path=final_xyz_path,
            model_pdb=Path(calc_cfg["model_pdb"]),
            freeze_indices_0based=freeze_atoms_final,
            ml_indices=ml_indices,
            hess_mm_indices=hess_mm_indices,
            movable_mm_indices=movable_mm_indices,
            frozen_layer_indices=frozen_layer_indices,
        )

        # summary.md and key_* outputs are disabled.
        click.echo(format_elapsed("[time] Elapsed Time for Opt", time_start))

    except ZeroStepLength:
        click.echo("ERROR: Step length fell below the minimum allowed (ZeroStepLength).", err=True)
        sys.exit(2)
    except OptimizationError as e:
        click.echo(f"ERROR: Optimization failed - {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled exception during optimization:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)
    finally:
        if prepared_input is not None:
            prepared_input.cleanup()
        # Release GPU memory so subsequent pipeline stages don't OOM
        base_calc = bias_calc = geometry = optimizer = mm_calc = macro_calc = macro_optimizer = None
        gc.collect()  # break cyclic refs inside torch.nn.Module
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


# Allow `python -m mlmm.opt` direct execution
if __name__ == "__main__":
    cli()
