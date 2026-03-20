"""
ML/MM vibrational frequency analysis with PHVA support and thermochemistry.

Example:
    mlmm freq -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0

For detailed documentation, see: docs/freq.md
"""

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
import numpy as np
import torch
import ase.units as units
import yaml
from ase import Atoms
from ase.data import atomic_masses
from ase.io import write

from pysisyphus.constants import AMU2AU, ANG2BOHR, AU2EV, BOHR2ANG
from pysisyphus.helpers import geom_loader

from .mlmm_calc import mlmm
from .defaults import FREQ_KW, THERMO_KW
from .opt import (
    CALC_KW as OPT_CALC_KW,
    GEOM_KW as OPT_GEOM_KW,
    _normalize_geom_freeze as _normalize_geom_freeze_opt,
    _parse_freeze_atoms as _parse_freeze_atoms_opt,
)
from .utils import (
    apply_ref_pdb_override,
    apply_layer_freeze_constraints,
    apply_yaml_overrides,
    convert_xyz_to_pdb,
    set_convert_file_enabled,
    is_convert_file_enabled,
    convert_xyz_like_outputs,
    deep_update,
    filter_calc_for_echo,
    format_elapsed,
    format_freeze_atoms_for_echo,
    load_yaml_dict,
    merge_freeze_atom_indices,
    prepare_input_structure,
    pretty_block,
    resolve_charge_spin_or_raise,
    parse_indices_string,
    build_model_pdb_from_bfactors,
    build_model_pdb_from_indices,
    strip_inherited_keys,
    yaml_section_has_key,
)
from .cli_utils import resolve_yaml_sources, load_merged_yaml_cfg, make_is_param_explicit


def _safe_masses_amu(atomic_numbers) -> np.ndarray:
    """Look up atomic masses with a clear error for unknown atomic numbers."""
    max_z = len(atomic_masses) - 1
    bad = [z for z in atomic_numbers if z < 0 or z > max_z or atomic_masses[z] == 0.0]
    if bad:
        raise ValueError(
            f"Unknown or unsupported atomic number(s): {sorted(set(bad))}. "
            "Check that all elements in the input structure are valid."
        )
    return np.array([atomic_masses[z] for z in atomic_numbers])


def _torch_device(auto: str = "auto") -> torch.device:
    if auto == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(auto)


# ===================================================================
#          Mass-weighted TR projection & vibrational analysis
# ===================================================================

def _build_tr_basis(coords_bohr_t: torch.Tensor,
                    masses_au_t: torch.Tensor) -> torch.Tensor:
    """
    Mass-weighted translation/rotation basis (Tx, Ty, Tz, Rx, Ry, Rz), shape (3N, r<=6).
    """
    device, dtype = coords_bohr_t.device, coords_bohr_t.dtype
    N = coords_bohr_t.shape[0]
    m_au = masses_au_t.to(dtype=dtype, device=device)
    m_sqrt = torch.sqrt(m_au).reshape(-1, 1)

    com = (m_au.reshape(-1, 1) * coords_bohr_t).sum(0) / m_au.sum()
    x = coords_bohr_t - com

    eye3 = torch.eye(3, dtype=dtype, device=device)
    cols = []
    for i in range(3):
        cols.append((eye3[i].repeat(N, 1) * m_sqrt).reshape(-1, 1))
    for i in range(3):
        rot = torch.cross(x, eye3[i].expand_as(x), dim=1) * m_sqrt
        cols.append(rot.reshape(-1, 1))
    return torch.cat(cols, dim=1)


def _tr_orthonormal_basis(coords_bohr_t: torch.Tensor,
                          masses_au_t: torch.Tensor,
                          rtol: float = 1e-12) -> Tuple[torch.Tensor, int]:
    """
    Orthonormalize TR basis in mass-weighted space by SVD. Returns (Q, rank).
    """
    B = _build_tr_basis(coords_bohr_t, masses_au_t)
    U, S, Vh = torch.linalg.svd(B, full_matrices=False)
    r = int((S > rtol * S.max()).sum().item())
    Q = U[:, :r]
    del B, S, Vh, U
    return Q, r


def _mw_projected_hessian(H_t: torch.Tensor,
                          coords_bohr_t: torch.Tensor,
                          masses_au_t: torch.Tensor) -> torch.Tensor:
    """
    Project out translations/rotations in mass-weighted space:
    Hmw = M^{-1/2} H M^{-1/2};  P = I - QQ^T;  Hmw_proj = P Hmw P

    To save memory, update **H_t in-place** (no clone) and return it.
    The output is explicitly symmetrized after TR projection.
    """
    if H_t.dtype != torch.float64:
        H_t = H_t.to(dtype=torch.float64)
    dtype, device = H_t.dtype, H_t.device
    with torch.no_grad():
        masses_amu_t = (masses_au_t / AMU2AU).to(dtype=dtype, device=device)
        m3 = torch.repeat_interleave(masses_amu_t, 3)
        # Use a single base vector for inverse sqrt mass and create views (no extra large allocations)
        inv_sqrt_m = torch.sqrt(1.0 / m3)
        inv_sqrt_m_col = inv_sqrt_m.view(1, -1)
        inv_sqrt_m_row = inv_sqrt_m.view(-1, 1)

        # In-place mass-weighting on input Hessian
        H_t.mul_(inv_sqrt_m_row)
        H_t.mul_(inv_sqrt_m_col)

        Q, _ = _tr_orthonormal_basis(coords_bohr_t, masses_au_t)  # (3N, r)
        Q = Q.to(dtype=dtype, device=device)
        Qt = Q.T

        QtH = Qt @ H_t                   # (r,3N)
        H_t.addmm_(Q, QtH, beta=1.0, alpha=-1.0)

        HQ = QtH.T                       # (3N,r)
        H_t.addmm_(HQ, Qt, beta=1.0, alpha=-1.0)

        QtHQ = QtH @ Q                   # (r,r)
        tmp = Q @ QtHQ                   # (3N,r)
        H_t.addmm_(tmp, Qt, beta=1.0, alpha=1.0)

        # Explicit symmetrization: H = (H + H^T) / 2
        H_sym = H_t.T.clone()
        H_t.add_(H_sym).mul_(0.5)
        del H_sym

        del masses_amu_t, m3, inv_sqrt_m, inv_sqrt_m_col, inv_sqrt_m_row
        del Q, Qt, QtH, HQ, QtHQ, tmp

        if torch.cuda.is_available() and device.type == "cuda":
            torch.cuda.empty_cache()
        return H_t


# ---- PHVA helper: mass-weighted Hessian without TR projection (for active subspace) ----
def _mass_weighted_hessian(H_t: torch.Tensor,
                           masses_au_t: torch.Tensor) -> torch.Tensor:
    """
    Return Hmw = M^{-1/2} H M^{-1/2} (no symmetrization/TR projection; in-place).
    """
    dtype, device = H_t.dtype, H_t.device
    with torch.no_grad():
        masses_amu_t = (masses_au_t / AMU2AU).to(dtype=dtype, device=device)
        m3 = torch.repeat_interleave(masses_amu_t, 3)
        inv_sqrt_m = torch.sqrt(1.0 / m3)
        inv_sqrt_m_col = inv_sqrt_m.view(1, -1)
        inv_sqrt_m_row = inv_sqrt_m.view(-1, 1)
        # In-place mass-weighting on input Hessian
        H_t.mul_(inv_sqrt_m_row)
        H_t.mul_(inv_sqrt_m_col)
        del masses_amu_t, m3, inv_sqrt_m, inv_sqrt_m_col, inv_sqrt_m_row
        return H_t


def _frequencies_cm_and_modes(H_t: torch.Tensor,
                              atomic_numbers: List[int],
                              coords_bohr: np.ndarray,
                              device: torch.device,
                              tol: float = 1e-6,
                              freeze_idx: Optional[List[int]] = None) -> Tuple[np.ndarray, torch.Tensor]:
    """
    Diagonalize a (possibly PHVA/active-subspace) TR-projected mass-weighted Hessian
    to obtain frequencies (cm^-1) and mass-weighted eigenvectors (modes).

    If `freeze_idx` is provided (list of 0-based atom indices), perform
    Partial Hessian Vibrational Analysis (PHVA). Supports two cases:

      A) Full Hessian given (3N×3N):
         1) build Hmw = M^{-1/2} H M^{-1/2}
         2) take the active subspace by removing DOF of frozen atoms
         3) perform TR projection **only in the active subspace** (always applied)
         4) diagonalize and embed eigenvectors back to 3N by zero-filling frozen DOF

      B) Already-reduced (active-block) Hessian given (3N_act×3N_act), e.g.
         when UMA is called with return_partial_hessian=True:
         1) mass-weight with **active** masses only
         2) TR projection in the active space
         3) diagonalize and embed back to 3N by zero-filling frozen DOF

    Returns:
      freqs_cm : (nmode,) numpy, negatives are imaginary
      modes    : (nmode, 3N) torch (mass-weighted eigenvectors)
    """
    with torch.no_grad():
        if H_t.dtype != torch.float64:
            H_t = H_t.to(dtype=torch.float64)
        Z = np.array(atomic_numbers, dtype=int)
        N = int(len(Z))
        masses_amu = np.array([atomic_masses[z] for z in Z])  # amu
        masses_au_t = torch.as_tensor(masses_amu * AMU2AU, dtype=H_t.dtype, device=device)
        coords_bohr_t = torch.as_tensor(coords_bohr.reshape(-1, 3), dtype=H_t.dtype, device=device)

        # --------------------------------------------
        # PHVA path (active DOF subspace with TR-proj)
        # --------------------------------------------
        if freeze_idx is not None and len(freeze_idx) > 0:
            # Active atom indices
            frozen_set = set(int(i) for i in freeze_idx if 0 <= int(i) < N)
            active_idx = [i for i in range(N) if i not in frozen_set]
            n_active = len(active_idx)
            if n_active == 0:
                # All atoms are frozen → no modes
                freqs_cm = np.zeros((0,), dtype=float)
                modes = torch.zeros((0, 3 * N), dtype=H_t.dtype, device=H_t.device)
                return freqs_cm, modes

            # Determine whether the provided Hessian is already the active block (3N_act×3N_act).
            expected_act_dim = 3 * n_active
            is_partial = (H_t.shape[0] == expected_act_dim and H_t.shape[1] == expected_act_dim)

            if is_partial:
                # --- Case B: Active-subspace Hessian supplied ---
                # Mass-weight using only active atoms → project TR modes in the active space
                # → diagonalise → embed back into the full space.
                masses_act = masses_au_t[active_idx]
                coords_act = coords_bohr_t[active_idx, :]

                # in-place mass-weight (active masses)
                Hmw_act = _mass_weighted_hessian(H_t, masses_act)

                # TR basis and projection in the active space
                Q, _ = _tr_orthonormal_basis(coords_act, masses_act)  # (3N_act, r)
                Qt = Q.T
                QtH = Qt @ Hmw_act
                Hmw_act.addmm_(Q, QtH, beta=1.0, alpha=-1.0)
                Hmw_act.addmm_(QtH.T, Qt, beta=1.0, alpha=-1.0)
                QtHQ = QtH @ Q
                Hmw_act.addmm_(Q @ QtHQ, Qt, beta=1.0, alpha=1.0)

                # Explicit symmetrization before eigendecomposition
                _t = Hmw_act.T.clone()
                Hmw_act.add_(_t).mul_(0.5)
                del _t
                omega2, Vsub = torch.linalg.eigh(Hmw_act, UPLO="U")

                # Free the (only) Hessian ASAP
                del Hmw_act
                del H_t
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

                sel = torch.abs(omega2) > tol
                omega2 = omega2[sel]
                Vsub = Vsub[:, sel]  # (3N_act, nsel)

                # Embed to full 3N (mass-weighted eigenvectors)
                modes = torch.zeros((Vsub.shape[1], 3 * N), dtype=Vsub.dtype, device=Vsub.device)
                mask_dof = torch.ones(3 * N, dtype=torch.bool, device=Vsub.device)
                for i in frozen_set:
                    mask_dof[3 * i:3 * i + 3] = False
                modes[:, mask_dof] = Vsub.T
                del Q, Qt, QtH, QtHQ, mask_dof

            else:
                # --- Case A: Full Hessian (3N×3N) supplied ---
                # Apply full mass-weighting → extract the active block → project TR modes in the active space.
                H_t = _mass_weighted_hessian(H_t, masses_au_t)

                # Build active mask (boolean) and immediately carve out the active block
                mask_dof = torch.ones(3 * N, dtype=torch.bool, device=H_t.device)
                for i in frozen_set:
                    mask_dof[3 * i:3 * i + 3] = False

                # Create the reduced Hessian; free the full one immediately to keep only one in VRAM
                H_act = H_t[mask_dof][:, mask_dof]
                del H_t
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
                H_t = H_act
                del H_act

                coords_act = coords_bohr_t[active_idx, :]
                masses_act = masses_au_t[active_idx]
                Q, _ = _tr_orthonormal_basis(coords_act, masses_act)  # (3N_act, r)
                Qt = Q.T

                QtH = Qt @ H_t
                H_t.addmm_(Q, QtH, beta=1.0, alpha=-1.0)

                H_t.addmm_(QtH.T, Qt, beta=1.0, alpha=-1.0)

                QtH = QtH @ Q
                H_t.addmm_(Q @ QtH, Qt, beta=1.0, alpha=1.0)

                # Explicit symmetrization before eigendecomposition
                _t = H_t.T.clone()
                H_t.add_(_t).mul_(0.5)
                del _t
                omega2, Vsub = torch.linalg.eigh(H_t, UPLO="U")

                # Free the (only) Hessian ASAP
                del H_t
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

                sel = torch.abs(omega2) > tol
                omega2 = omega2[sel]
                Vsub = Vsub[:, sel]  # (3N_act, nsel)

                modes = torch.zeros((Vsub.shape[1], 3 * N), dtype=Vsub.dtype, device=Vsub.device)
                modes[:, mask_dof] = Vsub.T  # (nsel, 3N_act) → place into active DOF
                del Vsub, mask_dof, Q, Qt, QtH

        else:
            # Legacy behavior: TR-projection in full DOF → diagonalization (both in-place)
            H_t = _mw_projected_hessian(H_t, coords_bohr_t, masses_au_t)
            # Explicit symmetrization before eigendecomposition
            _t = H_t.T.clone()
            H_t.add_(_t).mul_(0.5)
            del _t
            omega2, V = torch.linalg.eigh(H_t, UPLO="U")

            # Free the (only) Hessian ASAP
            del H_t
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

            sel = torch.abs(omega2) > tol
            omega2 = omega2[sel]
            modes = V[:, sel].T
            del V

        # Convert to frequencies (cm^-1)
        s_new = (units._hbar * 1e10 / np.sqrt(units._e * units._amu) * np.sqrt(AU2EV) / BOHR2ANG)
        hnu = s_new * torch.sqrt(torch.abs(omega2))
        hnu = torch.where(omega2 < 0, -hnu, hnu)
        freqs_cm = (hnu / units.invcm).detach().cpu().numpy()

        del omega2, hnu, sel
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        return freqs_cm, modes


def _mw_mode_to_cart(mode_mw_3N_t: torch.Tensor,
                     masses_au_t: torch.Tensor) -> np.ndarray:
    """
    Convert one mass-weighted eigenvector (3N,) to Cartesian (3N,) and L2-normalize.
    """
    with torch.no_grad():
        masses_amu_t = (masses_au_t / AMU2AU).to(dtype=mode_mw_3N_t.dtype, device=mode_mw_3N_t.device)
        m3 = torch.repeat_interleave(masses_amu_t, 3)
        v_cart = torch.sqrt(1.0 / m3) * mode_mw_3N_t
        v_cart.div_(torch.linalg.norm(v_cart))
        arr = v_cart.detach().cpu().numpy()
        del masses_amu_t, m3, v_cart
        return arr


def _calc_full_hessian_torch(
    geom,
    calc_kwargs: Dict[str, Any],
    device: torch.device,
    *,
    refresh_geom_meta: bool = False,
) -> Tuple[torch.Tensor, float]:
    """Return (Hessian torch tensor, energy Hartree) from the ML/MM calculator."""

    kw = dict(calc_kwargs or {})
    kw["out_hess_torch"] = True
    calc = mlmm(**kw)
    result = calc.get_hessian(geom.atoms, geom.coords)

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
    energy = float(result.get("energy", 0.0))

    del calc, result
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    return H, energy


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
    except Exception:
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
    if not bool(calc_cfg.get("use_bfactor_layers", True)):
        return False
    if calc_cfg.get("movable_cutoff") is not None:
        return False
    if calc_cfg.get("hess_cutoff") is not None:
        return False

    explicit_layer_lists = any(
        calc_cfg.get(key) is not None
        for key in ("hess_mm_atoms", "movable_mm_atoms", "frozen_mm_atoms")
    )
    if explicit_layer_lists:
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
    partial_mm_indices = hess_mm_indices if hess_mm_indices else movable_mm_indices

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
    ref_ang = geom.coords.reshape(-1, 3) * BOHR2ANG
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
            except Exception:
                # Fallback: generate MODEL/ENDMDL using ASE
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


# ===================================================================
#                         Defaults for CLI
# ===================================================================

# Geometry defaults — shared with opt.py
GEOM_KW: Dict[str, Any] = deepcopy(OPT_GEOM_KW)

# ML/MM calculator defaults — shared with opt.py
CALC_KW: Dict[str, Any] = deepcopy(OPT_CALC_KW)

# FREQ_KW and THERMO_KW are imported from .defaults


# ===================================================================
#                            CLI
# ===================================================================

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
    help="PDB defining atoms belonging to the ML region. Optional when --detect-layer is enabled.",
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
         "overrides calc.hessian_calc_mode from YAML. "
         "Default: 'FiniteDifference'. Use 'Analytical' when VRAM is sufficient.",
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

    # Validate input format: PDB directly, or XYZ with --ref-pdb
    suffix = input_path.suffix.lower()
    if suffix not in (".pdb", ".xyz"):
        click.echo("ERROR: --input must be a PDB or XYZ file.", err=True)
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
    # Keep command-level default for detect-layer unless YAML/explicit CLI overrides it.
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

    # CLI explicit overrides (after config YAML, before override YAML)
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

    if _is_param_explicit("hessian_calc_mode") and hessian_calc_mode is not None:
        calc_cfg["hessian_calc_mode"] = str(hessian_calc_mode)

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
                    "will_run_frequency_analysis": True,
                    "will_write_modes": True,
                    "will_dump_thermo_yaml": bool(thermo_cfg.get("dump", False)),
                    "backend": calc_cfg.get("backend", "uma"),
                    "embedcharge": bool(calc_cfg.get("embedcharge", False)),
                },
            )
        )
        click.echo("[dry-run] Validation complete. Frequency execution was skipped.")
        return

    out_dir_path.mkdir(parents=True, exist_ok=True)

    if detect_layer_enabled and layer_source_pdb.suffix.lower() != ".pdb":
        click.echo("ERROR: --detect-layer requires a PDB input (or --ref-pdb).", err=True)
        prepared_input.cleanup()
        sys.exit(1)

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

    click.echo(pretty_block("geom", format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")))
    echo_calc = format_freeze_atoms_for_echo(filter_calc_for_echo(calc_cfg), key="freeze_atoms")
    click.echo(pretty_block("calc", echo_calc))
    echo_freq = strip_inherited_keys({**freq_cfg, "out_dir": str(out_dir_path)}, FREQ_KW, mode="same")
    click.echo(pretty_block("freq", echo_freq))
    echo_thermo = strip_inherited_keys(thermo_cfg, THERMO_KW, mode="same")
    click.echo(pretty_block("thermo", echo_thermo))

    coord_type = geom_cfg.get("coord_type", "cart")
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
    three_layer_policy_applied = _align_three_layer_hessian_targets(calc_cfg, echo_fn=click.echo)

    # Determine active atoms based on mode
    active_dof_mode_lower = str(freq_cfg.get("active_dof_mode", active_dof_mode)).lower()
    active_indices, layer_sets = _resolve_active_atom_indices(calc_cfg, n_atoms, active_dof_mode_lower)
    ml_indices = layer_sets["ml"]
    hess_mm_indices = layer_sets["hess_mm"]
    movable_mm_indices = layer_sets["movable_mm"]
    partial_mm_indices = hess_mm_indices if hess_mm_indices else movable_mm_indices

    if active_dof_mode_lower == "all" or active_indices is None:
        active_indices = all_indices
        click.echo("[active-dof] Using all atoms for frequency analysis.")
    elif active_dof_mode_lower == "ml-only":
        click.echo(f"[active-dof] Using ML atoms only for frequency analysis (n={len(active_indices)}).")
    elif active_dof_mode_lower == "partial":
        if three_layer_policy_applied:
            click.echo(f"[active-dof] Using ML + MovableMM atoms for frequency analysis (n={len(active_indices)}).")
        elif hess_mm_indices:
            click.echo(f"[active-dof] Using ML + Hessian-target MM atoms for frequency analysis (n={len(active_indices)}).")
        else:
            click.echo(f"[active-dof] Using ML + MovableMM atoms for frequency analysis (n={len(active_indices)}).")
    elif active_dof_mode_lower == "unfrozen":
        click.echo(f"[active-dof] Using all non-frozen atoms for frequency analysis (n={len(active_indices)}).")
    else:
        active_indices = set(ml_indices) | set(partial_mm_indices)
        click.echo(f"[active-dof] Defaulting to partial mode (n={len(active_indices)}).")

    # Atoms not in active_indices become frozen for frequency analysis
    freeze_for_freq = sorted(all_indices - active_indices)
    # Also include any explicitly frozen atoms from config
    explicit_freeze = set(calc_cfg.get("freeze_atoms") or [])
    freeze_list = sorted(set(freeze_for_freq) | explicit_freeze)

    try:
        from .hessian_cache import load as _hess_load
        _cached_ts = _hess_load("ts")
        if _cached_ts is not None:
            click.echo("[freq] Reusing cached TS Hessian.")
            H_t = _cached_ts["hessian"]
            if isinstance(H_t, torch.Tensor):
                H_t = H_t.to(device=device)
            else:
                H_t = torch.as_tensor(H_t, device=device)
            energy_ha = _cached_ts.get("meta", {}).get("energy_ha")
            if energy_ha is None:
                energy_ha = float(geometry.energy)
        else:
            H_t, energy_ha = _calc_full_hessian_torch(geometry, calc_cfg, device)
        coords_bohr = geometry.coords.reshape(-1, 3)
        freqs_cm, modes_mw = _frequencies_cm_and_modes(
            H_t,
            geometry.atomic_numbers,
            coords_bohr,
            device,
            freeze_idx=freeze_list if freeze_list else None,
        )

        del H_t
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        order = (
            np.argsort(np.abs(freqs_cm))
            if freq_cfg["sort"] == "abs"
            else np.argsort(freqs_cm)
        )
        n_write = int(min(freq_cfg["max_write"], len(order)))
        click.echo(
            f"[INFO] Total modes: {len(freqs_cm)} → writing {n_write} mode(s) ({freq_cfg['sort']} ordering)."
        )

        ref_pdb_for_modes = source_path if source_path.suffix.lower() == ".pdb" else None
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

        (out_dir_path / "frequencies_cm-1.txt").write_text(
            "\n".join(f"{i+1:4d}  {float(freqs_cm[j]):+12.4f}" for i, j in enumerate(order)),
            encoding="utf-8",
        )

        del modes_mw
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        try:
            from thermoanalysis.QCData import QCData
            from thermoanalysis.constants import J2AU, J2CAL, NA
            from thermoanalysis.thermo import thermochemistry

            qc_data = {
                "coords3d": geometry.coords.reshape(-1, 3) * BOHR2ANG,
                "wavenumbers": freqs_cm,
                "scf_energy": float(energy_ha),
                "masses": masses_amu,
                "mult": int(calc_cfg["model_mult"]),
            }
            qc = QCData(qc_data, point_group="c1", mult=int(calc_cfg["model_mult"]))

            T = float(thermo_cfg["temperature"])
            p_atm = float(thermo_cfg["pressure_atm"])
            p_pa = p_atm * 101325.0  # Pa

            tr = thermochemistry(qc, T, pressure=p_pa)  # default: QRRHO

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
            click.echo(f"Temperature (K)         = {T:.2f}")
            click.echo(f"Pressure    (atm)       = {p_atm:.4f}")
            if freeze_list:
                click.echo("[NOTE] Thermochemistry uses active DOF (PHVA) due to frozen atoms.")
            click.echo(f"Number of Imaginary Freq = {n_imag:d}\n")

            def _ha(x): return f"{float(x): .6f} Ha"
            def _cal(x): return f"{float(x): .2f} cal/mol"
            def _calK(x): return f"{float(x): .2f} cal/(mol*K)"

            click.echo(f"Electronic Energy (EE)                 = {_ha(EE)}")
            click.echo(f"Zero-point Energy Correction           = {_ha(ZPE)}")
            click.echo(f"Thermal Correction to Energy           = {_ha(dE_therm)}")
            click.echo(f"Thermal Correction to Enthalpy         = {_ha(dH_therm)}")
            click.echo(f"Thermal Correction to Free Energy      = {_ha(dG_therm)}")
            click.echo(f"EE + Zero-point Energy                 = {_ha(sum_EE_ZPE)}")
            click.echo(f"EE + Thermal Energy Correction         = {_ha(sum_EE_thermal_E)}")
            click.echo(f"EE + Thermal Enthalpy Correction       = {_ha(sum_EE_thermal_H)}")
            click.echo(f"EE + Thermal Free Energy Correction    = {_ha(sum_EE_thermal_G)}")
            click.echo("")
            click.echo(f"E (Thermal)                            = {_cal(E_thermal_cal)}")
            click.echo(f"Heat Capacity (Cv)                     = {_calK(Cv_cal_per_Kmol)}")
            click.echo(f"Entropy (S)                            = {_calK(S_cal_per_Kmol)}")
            click.echo("")

            # Dump YAML when requested
            if bool(thermo_cfg["dump"]):
                out_yaml = out_dir_path / "thermoanalysis.yaml"
                payload = {
                    "temperature_K": T,
                    "pressure_atm": p_atm,
                    "num_imag_freq": n_imag,
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
                with out_yaml.open("w", encoding="utf-8") as f:
                    yaml.safe_dump(payload, f, sort_keys=False, allow_unicode=True)
                click.echo(f"[dump] Wrote thermoanalysis summary → {out_yaml}")

        except ImportError:
            click.echo("[thermo] WARNING: 'thermoanalysis' package not found; skipped thermochemistry summary.", err=True)
        except Exception as e:
            import traceback
            tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
            click.echo("Unhandled error during thermochemistry summary:\n" + textwrap.indent(tb, "  "), err=True)

        # summary.md and key_* outputs are disabled.
        click.echo(f"[DONE] Wrote modes and list → {out_dir_path}")

        click.echo(format_elapsed("[time] Elapsed Time for Freq", time_start))

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        import traceback
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled error during frequency analysis:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)
    finally:
        prepared_input.cleanup()
        # Release GPU memory so subsequent pipeline stages don't OOM
        geometry = H_t = modes = None
        gc.collect()  # break cyclic refs inside torch.nn.Module
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


# Allow `python -m mlmm.freq` direct execution
if __name__ == "__main__":
    cli()
