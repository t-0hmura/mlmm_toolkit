# mlmm_toolkit/tsopt.py

"""Partial Hessian guided Dimer / RS-I-RFO transition-state search with ML/MM.

Example:
    mlmm tsopt -i ts_guess.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
        -q 0 -m 1 --max-cycles 8000

For detailed documentation, see: docs/tsopt.md
"""

from __future__ import annotations

import sys
import textwrap
import os
import shutil
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Any, Optional, Tuple, List, Sequence

import click
import numpy as np
import torch
from click.core import ParameterSource
from ase import Atoms
from ase.io import write
from ase.data import atomic_masses
import ase.units as units
import time

# ---------------- pysisyphus / mlmm imports ----------------
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AMU2AU, AU2EV
from pysisyphus.calculators.Dimer import Dimer  # Dimer calculator (orientation-projected forces)

# RS-I-RFO optimizer for heavy mode
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer

# local helpers from mlmm
from .mlmm_calc import mlmm
from .defaults import OUT_DIR_TSOPT
from .defaults import (
    GEOM_KW_DEFAULT,
    MLMM_CALC_KW,
    OPT_BASE_KW,
    LBFGS_KW,
    DIMER_KW,
    HESSIAN_DIMER_KW,
    RSIRFO_KW,
    LAYEROPT_KW,
    TSOPT_MODE_ALIASES,
    BFACTOR_ML,
    BFACTOR_MOVABLE_MM,
    BFACTOR_FROZEN,
)
from .opt import (
    _parse_freeze_atoms as _parse_freeze_atoms_opt,
    _normalize_geom_freeze as _normalize_geom_freeze_opt,
    PartialHessianMicroIterationOptimizer,
)
from .utils import (
    append_xyz_trajectory as _append_xyz_trajectory,
    apply_layer_freeze_constraints,
    convert_xyz_to_pdb,
    deep_update,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    strip_inherited_keys,
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
)
from .freq import (
    _calc_full_hessian_torch as _freq_calc_full_hessian_torch,
    _torch_device,
    _build_tr_basis,
    _tr_orthonormal_basis,
    _mass_weighted_hessian,
)


# ===================================================================
#               Mass-weighted projection & vib analysis
# ===================================================================


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


def _calc_energy(geom, calc_kwargs: Dict[str, Any]) -> float:
    kw = dict(calc_kwargs or {})
    kw["out_hess_torch"] = False
    calc = mlmm(**kw)
    result = calc.get_energy(geom.atoms, geom.coords)
    energy = float(result.get("energy", 0.0))
    del result, calc
    _clear_cuda_cache()
    return energy


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


def _first_existing_artifact(out_dir: Path, patterns: Sequence[str]) -> Optional[Path]:
    """Resolve the first existing artifact for a list of relative patterns."""
    for pattern in patterns:
        if any(ch in pattern for ch in "*?[]"):
            for candidate in sorted(out_dir.glob(pattern)):
                if candidate.is_file():
                    return candidate.resolve()
            continue
        candidate = out_dir / pattern
        if candidate.is_file():
            return candidate.resolve()
    return None


def _link_or_copy_file(src: Path, dst: Path) -> bool:
    """Create a symlink when possible; fall back to copy."""
    try:
        if dst.exists() or dst.is_symlink():
            if dst.is_dir():
                return False
            dst.unlink()
        rel = os.path.relpath(src, start=dst.parent)
        dst.symlink_to(rel)
        return True
    except Exception:
        try:
            shutil.copy2(src, dst)
            return True
        except Exception:
            return False


def _write_output_summary_md(out_dir: Path) -> None:
    """summary.md and key_* outputs are disabled."""
    return None

def _resolve_yaml_sources(
    config_yaml: Optional[Path],
    override_yaml: Optional[Path],
    args_yaml_legacy: Optional[Path],
) -> Tuple[Optional[Path], Optional[Path], bool]:
    if override_yaml is not None and args_yaml_legacy is not None:
        raise click.BadParameter(
            "Use a single YAML source option."
        )
    if args_yaml_legacy is not None:
        return config_yaml, args_yaml_legacy, True
    return config_yaml, override_yaml, False


def _load_merged_yaml_cfg(
    config_yaml: Optional[Path],
    override_yaml: Optional[Path],
) -> Dict[str, Any]:
    merged: Dict[str, Any] = {}
    deep_update(merged, load_yaml_dict(config_yaml))
    deep_update(merged, load_yaml_dict(override_yaml))
    return merged




def _mw_projected_hessian_inplace(H_t: torch.Tensor,
                                  coords_bohr_t: torch.Tensor,
                                  masses_au_t: torch.Tensor,
                                  freeze_idx: Optional[List[int]] = None) -> torch.Tensor:
    """
    Mass-weight H in-place, optionally restrict to active DOF subspace (PHVA) and
    project out TR motions (in that subspace), also in-place. No explicit symmetrization.
    Returns the (possibly reduced) Hessian to be diagonalized.
    """
    dtype, device = H_t.dtype, H_t.device
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
            H_t = H_t[mask_dof][:, mask_dof]
            # TR basis and projection in active subspace (in-place)
            coords_act = coords_bohr_t[active_idx, :]
            masses_act = masses_au_t[active_idx]
            Q, _ = _tr_orthonormal_basis(coords_act, masses_act)  # (3N_act, r)
            Qt = Q.T
            QtH = Qt @ H_t
            H_t.addmm_(Q, QtH, beta=1.0, alpha=-1.0)
            H_t.addmm_((QtH.T), Qt, beta=1.0, alpha=-1.0)
            QtHQ = QtH @ Q
            H_t.addmm_(Q @ QtHQ, Qt, beta=1.0, alpha=1.0)
            del Q, Qt, QtH, QtHQ, mask_dof, coords_act, masses_act, active_idx, frozen
        else:
            # Full DOF: mass-weight + TR projection in-place
            H_t = _mass_weighted_hessian(H_t, masses_au_t)
            Q, _ = _tr_orthonormal_basis(coords_bohr_t, masses_au_t)  # (3N, r)
            Qt = Q.T
            QtH = Qt @ H_t
            H_t.addmm_(Q, QtH, beta=1.0, alpha=-1.0)
            H_t.addmm_(QtH.T, Qt, beta=1.0, alpha=-1.0)
            QtHQ = QtH @ Q
            H_t.addmm_(Q @ QtHQ, Qt, beta=1.0, alpha=1.0)
            del Q, Qt, QtH, QtHQ
        _clear_cuda_cache()
        return H_t


def _self_check_tr_projection(H_t: torch.Tensor,
                              coords_bohr_t: torch.Tensor,
                              masses_au_t: torch.Tensor,
                              freeze_idx: Optional[List[int]] = None) -> Tuple[float, float, int]:
    """
    (Low-VRAM) Skipped heavy clone-based check. Keep signature for compatibility.
    """
    # To conserve VRAM, do not allocate a full clone for checks.
    # Return zeros and rank estimate based on current TR basis only.
    with torch.no_grad():
        N = coords_bohr_t.shape[0]
        if freeze_idx:
            frozen = set(int(i) for i in freeze_idx if 0 <= int(i) < N)
            active_idx = [i for i in range(N) if i not in frozen]
            if len(active_idx) == 0:
                return 0.0, 0.0, 0
            coords_act = coords_bohr_t[active_idx, :]
            masses_act = masses_au_t[active_idx]
            _, r = _tr_orthonormal_basis(coords_act, masses_act)
            return 0.0, 0.0, r
        else:
            _, r = _tr_orthonormal_basis(coords_bohr_t, masses_au_t)
            return 0.0, 0.0, r


def _mode_direction_by_root(H_t: torch.Tensor,
                            coords_bohr_t: torch.Tensor,
                            masses_au_t: torch.Tensor,
                            root: int = 0,
                            freeze_idx: Optional[List[int]] = None) -> np.ndarray:
    """
    Get the eigenvector (Cartesian space) corresponding to the `root`-th most negative
    eigenvalue (root=0: most negative) of the mass-weighted, TR-projected Hessian.
    PHVA (active-subspace) is applied if freeze_idx is provided: frozen DOFs are zero.
    root==0 prefers torch.lobpcg; fallback to eigh (UPLO='U').
    """
    with torch.no_grad():
        # In-place: mass weight + (active-subspace) TR projection
        Hmw_proj = _mw_projected_hessian_inplace(H_t, coords_bohr_t, masses_au_t, freeze_idx=freeze_idx)

        # Solve eigenproblem in the (possibly reduced) space
        if int(root) == 0:
            try:
                w, v_mw_sub = torch.lobpcg(Hmw_proj, k=1, largest=False)
                u_mw_sub = v_mw_sub[:, 0]
            except Exception:
                evals_f, evecs_f = torch.linalg.eigh(Hmw_proj, UPLO="U")
                u_mw_sub = evecs_f[:, torch.argmin(evals_f)]
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
            u_mw_sub = evecs_mw[:, pick]
            del evals, evecs_mw

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
        m3 = torch.repeat_interleave(masses_amu_t, 3)
        inv_sqrt_m = torch.sqrt(1.0 / m3)
        v = inv_sqrt_m * u_mw
        v = v / torch.linalg.norm(v)
        mode = v.reshape(-1, 3).detach().cpu().numpy()

        del masses_amu_t, m3, inv_sqrt_m, v, u_mw, u_mw_sub
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
    g = np.array(geom.gradient, dtype=float).reshape(-1)
    geom.set_calculator(None)
    del calc
    _clear_cuda_cache()
    return g


def _frequencies_cm_and_modes(H_t: torch.Tensor,
                              atomic_numbers: List[int],
                              coords_bohr: np.ndarray,
                              device: torch.device,
                              tol: float = 1e-6,
                              freeze_idx: Optional[List[int]] = None) -> Tuple[np.ndarray, torch.Tensor]:
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

        # in-place mass-weight + (active-subspace) TR projection
        Hmw = _mw_projected_hessian_inplace(H_t, coords_bohr_t, masses_au_t, freeze_idx=freeze_idx)

        # eigensolve (upper triangle)
        omega2, Vsub = torch.linalg.eigh(Hmw, UPLO="U")

        sel = torch.abs(omega2) > tol
        omega2 = omega2[sel]
        Vsub = Vsub[:, sel]  # (3N_act or 3N, nsel)

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

        del omega2, Vsub, sel, masses_amu, masses_au_t, coords_bohr_t, Hmw
        _clear_cuda_cache(H_t)
        return freqs_cm, modes


def _frequencies_cm_only(H_t: torch.Tensor,
                         atomic_numbers: List[int],
                         coords_bohr: np.ndarray,
                         device: torch.device,
                         tol: float = 1e-6,
                         freeze_idx: Optional[List[int]] = None) -> np.ndarray:
    """
    Frequencies only (PHVA/TR in-place; no eigenvectors) for quick checks.
    """
    with torch.no_grad():
        Z = np.array(atomic_numbers, dtype=int)
        masses_amu = np.array([atomic_masses[z] for z in Z])  # amu
        masses_au_t = torch.as_tensor(masses_amu * AMU2AU, dtype=H_t.dtype, device=device)
        coords_bohr_t = torch.as_tensor(coords_bohr.reshape(-1, 3), dtype=H_t.dtype, device=device)

        Hmw = _mw_projected_hessian_inplace(H_t, coords_bohr_t, masses_au_t, freeze_idx=freeze_idx)
        omega2 = torch.linalg.eigvalsh(Hmw, UPLO="U")

        sel = torch.abs(omega2) > tol
        omega2 = omega2[sel]

        freqs_cm = _omega2_to_freqs_cm(omega2)

        del omega2, sel
        _clear_cuda_cache(H_t)
        return freqs_cm


def _write_mode_trj_and_pdb(geom,
                            mode_vec_3N: np.ndarray,
                            out_trj: Path,
                            out_pdb: Path,
                            amplitude_ang: float = 0.25,
                            n_frames: int = 20,
                            comment: str = "imag mode") -> None:
    """
    Write a single imaginary mode animation both as _trj.xyz (XYZ-like) and .pdb.
    """
    ref_ang = geom.coords.reshape(-1, 3) * BOHR2ANG
    mode = mode_vec_3N.reshape(-1, 3).copy()
    mode /= np.linalg.norm(mode)

    # _trj.xyz (XYZ-like concatenation)
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

    # .pdb (MODEL/ENDMDL)
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
    filename_prefix: str = "final_imag_mode",
    amplitude_ang: float = 0.8,
    n_frames: int = 20,
) -> int:
    """
    Write all imaginary modes (freq < -|threshold|) to vib_dir.

    Returns:
        Number of mode trajectories written.
    """
    neg_idx = np.where(freqs_cm < -abs(neg_freq_thresh_cm))[0]
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

        stem = f"{filename_prefix}_{rank:02d}_mode{mode_idx:04d}_{freq:+.2f}cm-1"
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
        )
        del mode_mw, v_cart
        written += 1

    del masses_amu, sqrt_m3, order, neg_idx
    _clear_cuda_cache()
    return written


# ===================================================================
#            Active-subspace helpers & Bofill update
# ===================================================================

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
    device = H_full.device
    m = torch.as_tensor(mask_dof, device=device, dtype=torch.bool)
    return H_full[m][:, m].clone()


def _embed_active_vector(vec_act: torch.Tensor,
                         mask_dof: np.ndarray,
                         total_3N: int) -> torch.Tensor:
    """
    Embed a (3N_act,) vector back to full (3N,) with zeros on frozen DOFs.
    """
    device = vec_act.device
    dtype = vec_act.dtype
    full = torch.zeros(total_3N, device=device, dtype=dtype)
    m = torch.as_tensor(mask_dof, device=device, dtype=torch.bool)
    full[m] = vec_act
    return full


def _mw_tr_project_active_inplace(H_act: torch.Tensor,
                                  coords_act_t: torch.Tensor,
                                  masses_act_au_t: torch.Tensor) -> torch.Tensor:
    """
    Mass-weight & project TR in the *active* subspace (in-place; no explicit symmetrization).
    """
    with torch.no_grad():
        # mass-weight
        masses_amu_t = (masses_act_au_t / AMU2AU).to(dtype=H_act.dtype, device=H_act.device)
        m3 = torch.repeat_interleave(masses_amu_t, 3)
        inv_sqrt_m_col = torch.sqrt(1.0 / m3).view(1, -1)
        inv_sqrt_m_row = inv_sqrt_m_col.view(-1, 1)
        H_act.mul_(inv_sqrt_m_row)
        H_act.mul_(inv_sqrt_m_col)
        # TR basis & projection
        Q, _ = _tr_orthonormal_basis(coords_act_t, masses_act_au_t)  # (3N_act, r)
        Qt = Q.T
        QtH = Qt @ H_act
        H_act.addmm_(Q, QtH, beta=1.0, alpha=-1.0)
        H_act.addmm_(QtH.T, Qt, beta=1.0, alpha=-1.0)
        QtHQ = QtH @ Q
        H_act.addmm_(Q @ QtHQ, Qt, beta=1.0, alpha=1.0)
        del masses_amu_t, m3, inv_sqrt_m_col, inv_sqrt_m_row, Q, Qt, QtH, QtHQ
        return H_act


def _frequencies_from_Hact(H_act: torch.Tensor,
                           atomic_numbers: List[int],
                           coords_bohr: np.ndarray,
                           active_idx: List[int],
                           device: torch.device,
                           tol: float = 1e-6) -> np.ndarray:
    """
    Frequencies (cm^-1) computed from active-block Hessian with active-space TR projection.
    """
    with torch.no_grad():
        coords_act = torch.as_tensor(coords_bohr.reshape(-1, 3)[active_idx, :], dtype=H_act.dtype, device=device)
        masses_act_au = torch.as_tensor([atomic_masses[int(z)] * AMU2AU
                                         for z in np.array(atomic_numbers, int)[active_idx]],
                                        dtype=H_act.dtype, device=device)
        Hmw = H_act.clone()
        _mw_tr_project_active_inplace(Hmw, coords_act, masses_act_au)
        omega2 = torch.linalg.eigvalsh(Hmw, UPLO="U")
        sel = torch.abs(omega2) > tol
        omega2 = omega2[sel]
        freqs_cm = _omega2_to_freqs_cm(omega2)
        del coords_act, masses_act_au, Hmw, omega2, sel
        _clear_cuda_cache(H_act)
        return freqs_cm


def _modes_from_Hact_embedded(H_act: torch.Tensor,
                              atomic_numbers: List[int],
                              coords_bohr: np.ndarray,
                              active_idx: List[int],
                              device: torch.device,
                              tol: float = 1e-6) -> Tuple[np.ndarray, torch.Tensor]:
    """
    Diagonalize active-block Hessian with mass-weight/TR in active space and return:
      freqs_cm : (nmode,)
      modes    : (nmode, 3N) mass-weighted eigenvectors embedded to full 3N (torch)
    """
    with torch.no_grad():
        N = len(atomic_numbers)
        coords_act = torch.as_tensor(coords_bohr.reshape(-1, 3)[active_idx, :], dtype=H_act.dtype, device=device)
        masses_act_au = torch.as_tensor([atomic_masses[int(z)] * AMU2AU
                                         for z in np.array(atomic_numbers, int)[active_idx]],
                                        dtype=H_act.dtype, device=device)
        Hmw = H_act.clone()
        _mw_tr_project_active_inplace(Hmw, coords_act, masses_act_au)
        omega2, Vsub = torch.linalg.eigh(Hmw, UPLO="U")
        sel = torch.abs(omega2) > tol
        omega2 = omega2[sel]
        Vsub = Vsub[:, sel]  # (3N_act, nsel)

        # Embed to full 3N (mass-weighted eigenvectors)
        modes_full = torch.zeros((Vsub.shape[1], 3 * N), dtype=Hmw.dtype, device=device)
        mask_dof = _active_mask_dof(N, list(set(range(N)) - set(active_idx)))  # give frozen list
        mask_t = torch.as_tensor(mask_dof, dtype=torch.bool, device=device)
        modes_full[:, mask_t] = Vsub.T
        # frequencies
        freqs_cm = _omega2_to_freqs_cm(omega2)

        del coords_act, masses_act_au, Hmw, omega2, Vsub, mask_t
        _clear_cuda_cache(H_act)
        return freqs_cm, modes_full


def _mode_direction_by_root_from_Hact(H_act: torch.Tensor,
                                      coords_bohr: np.ndarray,
                                      atomic_numbers: List[int],
                                      masses_au_t: torch.Tensor,
                                      active_idx: List[int],
                                      device: torch.device,
                                      root: int = 0) -> np.ndarray:
    """
    TS direction from the *active* Hessian block. Mass-weighting/TR are done in the
    active space. Result is embedded back to full 3N in Cartesian space.
    """
    with torch.no_grad():
        N = len(atomic_numbers)
        coords_act = torch.as_tensor(coords_bohr.reshape(-1, 3)[active_idx, :], dtype=H_act.dtype, device=device)
        masses_act_au = masses_au_t[active_idx].to(device=device, dtype=H_act.dtype)
        # mass-weight + TR in active space
        Hmw = H_act.clone()
        _mw_tr_project_active_inplace(Hmw, coords_act, masses_act_au)

        # eigenvector for requested root
        if int(root) == 0:
            try:
                w, V = torch.lobpcg(Hmw, k=1, largest=False)
                u_mw = V[:, 0]
            except Exception:
                vals, vecs = torch.linalg.eigh(Hmw, UPLO="U")
                u_mw = vecs[:, torch.argmin(vals)]
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
            u_mw = vecs[:, pick]
            del vals, vecs

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

        del coords_act, masses_act_au, masses_act_amu, m3, v_cart_act, full, mask_t, Hmw, u_mw
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
) -> List[int]:
    """
    Select a subset of imaginary modes to flatten for a geometry.
    """
    neg_idx_all = np.where(freqs_cm < -abs(neg_freq_thresh_cm))[0]
    if len(neg_idx_all) <= 1:
        return []

    order = np.argsort(freqs_cm[neg_idx_all])
    sorted_neg = neg_idx_all[order]
    root_clamped = max(0, min(int(root), len(order) - 1))
    primary_idx = sorted_neg[root_clamped]
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
) -> bool:
    """
    Flatten extra imaginary modes for a geometry (single pass).
    """
    neg_idx_all = np.where(freqs_cm < -abs(neg_freq_thresh_cm))[0]
    if len(neg_idx_all) <= 1:
        return False

    targets = _select_flatten_targets_for_geom(
        freqs_cm,
        modes,
        geom.cart_coords,
        neg_freq_thresh_cm,
        root,
        flatten_sep_cutoff,
        flatten_k,
    )
    if not targets:
        return False

    mass_scale = np.sqrt(12.011 / masses_amu)[:, None]
    amp_bohr = float(flatten_amp_ang) / BOHR2ANG

    for idx in targets:
        v_mw = modes[idx].detach().cpu().numpy().reshape(-1, 3)
        m3 = np.repeat(masses_amu, 3).reshape(-1, 3)
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

        # Move towards lower energy
        if E_plus <= E_minus:
            geom.coords = plus.reshape(-1)
        else:
            geom.coords = minus.reshape(-1)

    return True


def _get_active_dof_indices(
    calc_cfg: Dict[str, Any],
    n_atoms: int,
    active_dof_mode: str,
    freeze_atoms_final: List[int],
) -> Optional[List[int]]:
    """
    Determine active DOF indices based on active_dof_mode.

    Returns:
        None if all atoms should be active, otherwise list of atom indices (0-based).
    """
    if active_dof_mode.lower() == "all":
        return None  # All atoms active

    # Get layer indices from calculator
    try:
        temp_calc = mlmm(**calc_cfg)
        calc_core = temp_calc.core if hasattr(temp_calc, 'core') else temp_calc
        ml_indices = set(getattr(calc_core, 'ml_indices', []))
        hess_mm_indices = set(getattr(calc_core, 'hess_mm_indices', []))
        movable_mm_indices = set(getattr(calc_core, 'movable_mm_indices', []))
        del temp_calc
    except Exception:
        ml_indices = set()
        hess_mm_indices = set()
        movable_mm_indices = set()

    mode_lower = active_dof_mode.lower()
    if mode_lower == "ml-only":
        # ML only
        active_indices = sorted(ml_indices)
    elif mode_lower == "partial":
        # ML + Hessian-target MM (default)
        active_indices = sorted(ml_indices | hess_mm_indices)
    elif mode_lower == "unfrozen":
        # All non-frozen atoms
        active_indices = sorted(ml_indices | hess_mm_indices | movable_mm_indices)
    else:
        return None

    if not active_indices:
        return None

    return active_indices


def _bofill_update_active(H_act: torch.Tensor,
                          delta_act: np.ndarray,
                          g_new_act: np.ndarray,
                          g_old_act: np.ndarray,
                          eps: float = 1e-12) -> torch.Tensor:
    """
    Memory-efficient Bofill update on the *active* Cartesian Hessian block.
    Apply symmetric rank-1/2 updates directly **in place** using only the **upper triangle**
    index set (and mirror to the lower) to avoid allocating large NxN temporaries.
    No explicit (H+H^T)/2 symmetrization step is performed.
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
    denom_ms = d_dot_xi if torch.abs(d_dot_xi) > eps else torch.sign(d_dot_xi + 0.0) * eps
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
        H_act[idx, idx].add_(alpha * xi[idx] * xi[idx]
                             + beta * d[idx] * d[idx]
                             + 2.0 * gamma * d[idx] * xi[idx])

    # Off-diagonal (i < j): symmetric update
    if off.any():
        i = iu0[off]; j = iu1[off]
        inc = (alpha * xi[i] * xi[j]
               + beta * d[i] * d[j]
               + gamma * (d[i] * xi[j] + xi[i] * d[j]))
        H_act[i, j].add_(inc)
        H_act[j, i].add_(inc)

    return H_act


# ===================================================================
#                        HessianDimer (extended)
# ===================================================================

class HessianDimer:
    """
    Dimer-based TS search with periodic Hessian updates.

    Extensions in this implementation:
      - `root` parameter: choose which imaginary mode to follow (0 = most negative).
      - Pass-through kwargs: `dimer_kwargs` and `lbfgs_kwargs` to tune internals.
      - Hard cap on total LBFGS steps across segments: `max_total_cycles`.
      - PHVA (active DOF subspace) + TR projection for mode picking,
        respecting ``freeze_atoms``, with in-place operations. When ``root == 0`` the
        implementation prefers LOBPCG.
      - The flatten loop uses a *Bofill*-updated active Hessian block, so the
        expensive exact Hessian is evaluated only once before the flatten loop and
        once at the end for the final frequency analysis.
      - UMA calculator kwargs accept ``freeze_atoms`` and ``hessian_calc_mode`` and
        default to ``return_partial_hessian=True`` (active-block Hessian when frozen).
    """

    def __init__(self,
                 fn: str,
                 out_dir: str = "./result_dimer",
                 thresh_loose: str = "gau_loose",
                 thresh: str = "baker",
                 update_interval_hessian: int = 50,
                 neg_freq_thresh_cm: float = 5.0,
                 flatten_amp_ang: float = 0.20,
                 flatten_max_iter: int = 20,
                 mem: int = 100000,
                 use_lobpcg: bool = True,  # kept for backward compat (not used when root!=0)
                 calc_kwargs: Optional[dict] = None,
                 device: str = "auto",
                 dump: bool = False,
                 #
                 # New:
                 root: int = 0,
                 dimer_kwargs: Optional[Dict[str, Any]] = None,
                 lbfgs_kwargs: Optional[Dict[str, Any]] = None,
                 max_total_cycles: int = 10000,
                 #
                # Pass geom kwargs so freeze-atoms and YAML geometry overrides apply on the light path (fix #1)
                 geom_kwargs: Optional[Dict[str, Any]] = None,
                 # New: Use partial Hessian for imaginary mode detection in flatten loop
                 partial_hessian_flatten: bool = True,
                 # Spatial separation for flatten mode selection (from pdb2reaction)
                 flatten_sep_cutoff: float = 0.0,
                 flatten_k: int = 10,
                 flatten_loop_bofill: bool = False,
                 ) -> None:

        self.fn = fn
        self.out_dir = Path(out_dir); self.out_dir.mkdir(parents=True, exist_ok=True)
        self.vib_dir = self.out_dir / "vib"; self.vib_dir.mkdir(parents=True, exist_ok=True)

        self.thresh_loose = thresh_loose
        self.thresh = thresh
        self.update_interval_hessian = int(update_interval_hessian)
        self.neg_freq_thresh_cm = float(neg_freq_thresh_cm)
        self.flatten_amp_ang = float(flatten_amp_ang)
        self.flatten_max_iter = int(flatten_max_iter)
        self.mem = int(mem)
        self.use_lobpcg = bool(use_lobpcg)  # used only when root==0 shortcut
        self.root = int(root)
        self.dimer_kwargs = dict(dimer_kwargs or {})
        self.lbfgs_kwargs = dict(lbfgs_kwargs or {})
        self.max_total_cycles = int(max_total_cycles)
        self.partial_hessian_flatten = bool(partial_hessian_flatten)
        # Spatial separation for flatten mode selection
        self.flatten_sep_cutoff = float(flatten_sep_cutoff)
        self.flatten_k = int(flatten_k)
        self.flatten_loop_bofill = bool(flatten_loop_bofill)

        # Track total cycles globally across ALL loops/segments (fix #2)
        self._cycles_spent = 0

        # Hessian caching for 0-step convergence (avoid redundant recalculation)
        self._raw_hessian_cache_cpu: Optional[torch.Tensor] = None
        self._raw_hessian_coords_cpu: Optional[np.ndarray] = None
        self._last_active_idx: Optional[List[int]] = None
        self._last_active_mask_dof: Optional[np.ndarray] = None

        # ML/MM calculator settings
        self.calc_kwargs = dict(calc_kwargs or {})
        self.calc_kwargs.setdefault("out_hess_torch", False)

        # Geometry & masses (use provided geom kwargs so freeze_atoms etc. apply)
        gkw = dict(geom_kwargs or {})
        coord_type = gkw.pop("coord_type", "cart")
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
        self.calc_kwargs_partial["mm_fd"] = False
        self.calc_kwargs_partial["return_partial_hessian"] = False
        self.calc_kwargs_partial["out_hess_torch"] = True
        self.calc_kwargs_full = dict(self.calc_kwargs)
        self.calc_kwargs_full.setdefault("mm_fd", True)
        self.calc_kwargs_full["return_partial_hessian"] = False
        self.calc_kwargs_full["out_hess_torch"] = True
        self.geom = geom_loader(fn, coord_type=coord_type, **gkw)
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
        self.freeze_atoms: List[int] = list(gkw.get("freeze_atoms", [])) if "freeze_atoms" in gkw else []

        # Device
        self.device = _torch_device(device)
        self.masses_au_t = self.masses_au_t.to(self.device)

        # temp file for Dimer orientation (N_raw)
        self.mode_path = self.out_dir / ".dimer_mode.dat"

        self.dump = bool(dump)
        self.optim_all_path = self.out_dir / "optimization_all_trj.xyz"

    # ----- One dimer segment for up to n_steps; returns (steps_done, converged) -----
    def _dimer_segment(self, threshold: str, n_steps: int) -> Tuple[int, bool]:
        # Dimer calculator using current mode as initial N
        calc_sp = mlmm(**self.calc_kwargs)

        # Merge user dimer kwargs (but enforce N_raw & write_orientations)
        dimer_kwargs = dict(self.dimer_kwargs)
        dimer_kwargs.update({
            "calculator": calc_sp,
            "N_raw": str(self.mode_path),
            "write_orientations": False,  # runner override to reduce IO
            "seed": 0,                    # runner override for determinism
            "mem": self.mem,              # accepted by Calculator base through **kwargs
        })
        dimer = Dimer(**dimer_kwargs)

        self.geom.set_calculator(dimer)

        # LBFGS kwargs: enforce thresh/max_cycles/out_dir/dump; allow others
        lbfgs_kwargs = dict(self.lbfgs_kwargs)
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
    def _cache_raw_hessian_cpu(self, H: torch.Tensor) -> None:
        """Cache the raw Hessian on CPU for the current geometry."""
        self._raw_hessian_cache_cpu = H.detach().cpu().clone()
        self._raw_hessian_coords_cpu = self.geom.cart_coords.copy()

    def _reuse_cached_hessian(self) -> Optional[torch.Tensor]:
        """If the cached geometry matches current, return cached Hessian on device."""
        if self._raw_hessian_cache_cpu is None or self._raw_hessian_coords_cpu is None:
            return None
        if not np.array_equal(self.geom.cart_coords, self._raw_hessian_coords_cpu):
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
            cached = self._reuse_cached_hessian()
            if cached is not None:
                click.echo("[tsopt] Reusing cached raw Hessian (0-step convergence).")
                return cached
        H = _calc_full_hessian_torch(self.geom, calc_kwargs, self.device)
        self._cache_raw_hessian_cpu(H)
        return H

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
            pass

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
    def _dimer_loop(self, threshold: str) -> Tuple[int, bool]:
        """
        Run multiple LBFGS segments separated by periodic Hessian-based mode updates.
        Consumes from a *global* cycle budget self.max_total_cycles.

        Returns:
            (steps_in_this_call, zero_step_converged)
        where `zero_step_converged` is True iff the loop terminated by convergence
        without changing the geometry (i.e., 0-step convergence).
        """
        steps_in_this_call = 0
        zero_step_converged = False
        while True:
            remaining_global = max(0, self.max_total_cycles - self._cycles_spent)
            if remaining_global == 0:
                break
            steps_this = min(self.update_interval_hessian, remaining_global)
            coords_before = self.geom.cart_coords.copy()
            steps, ok = self._dimer_segment(threshold, steps_this)
            self._cycles_spent += steps
            steps_in_this_call += steps
            if ok:
                # Check if geometry unchanged (0-step convergence)
                if np.array_equal(self.geom.cart_coords, coords_before):
                    zero_step_converged = True
                break
            # If budget exhausted after this segment, stop before doing a Hessian update
            if (self.max_total_cycles - self._cycles_spent) <= 0:
                break
            # Update mode from Hessian (respect freeze atoms via PHVA)
            # Ensure VRAM is fully released after dimer segment before heavy Hessian computation
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            H_t = _calc_full_hessian_torch(self.geom, self.calc_kwargs_partial, self.device)
            N = len(self.geom.atomic_numbers)
            coords_bohr_t = torch.as_tensor(self.geom.coords.reshape(-1, 3),
                                            dtype=H_t.dtype, device=H_t.device)
            # full vs active-block Hessian
            if H_t.size(0) == 3 * N:
                mode_xyz = _mode_direction_by_root(
                    H_t, coords_bohr_t, self.masses_au_t,
                    root=self.root, freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None
                )
            else:
                # partial (active) Hessian returned by UMA
                active_idx, _ = self._resolve_hessian_active_subspace(H_t, N)
                mode_xyz = _mode_direction_by_root_from_Hact(
                    H_t, self.geom.coords.reshape(-1, 3), self.geom.atomic_numbers,
                    self.masses_au_t, active_idx, self.device, root=self.root
                )
            np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")
            del H_t, coords_bohr_t, mode_xyz
            _clear_cuda_cache()
        return steps_in_this_call, zero_step_converged

    # ----- Flatten (resolve multiple imaginary modes) -----
    def _flatten_once(self) -> bool:
        """
        Legacy: exact-Hessian-based flattening (kept for reference / fallback).
        """
        H_t = _calc_full_hessian_torch(self.geom, self.calc_kwargs_full, self.device)
        freqs_cm, modes = _frequencies_cm_and_modes(
            H_t, self.geom.atomic_numbers, self.geom.coords.reshape(-1, 3), self.device,
            freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None
        )
        del H_t
        neg_idx_all = np.where(freqs_cm < -abs(self.neg_freq_thresh_cm))[0]
        if len(neg_idx_all) <= 1:
            del modes
            return False

        # Identify the "primary" imaginary by root among negative modes
        order = np.argsort(freqs_cm[neg_idx_all])  # ascending (more negative first)
        root_clamped = max(0, min(self.root, len(order) - 1))
        primary_idx = neg_idx_all[order[root_clamped]]

        targets = [i for i in neg_idx_all if i != primary_idx]
        if not targets:
            del modes
            return False

        # Reference structure and energy
        ref = self.geom.coords.reshape(-1, 3).copy()
        _ = _calc_energy(self.geom, self.calc_kwargs)  # E_ref (unused, but keeps semantics)

        # mass scaling so that carbon ~ amplitude
        mass_scale = np.sqrt(12.011 / self.masses_amu)[:, None]
        amp_bohr = self.flatten_amp_ang / BOHR2ANG

        disp_total = np.zeros_like(ref)
        for idx in targets:
            v_mw = modes[idx].detach().cpu().numpy().reshape(-1, 3)  # mass-weighted eigenvector embedded to 3N
            # Convert to Cartesian step direction already done downstream in writer,
            # but for flattening we only need a normalized direction in Cartesian:
            # use masses to unweight:
            m3 = np.repeat(self.masses_amu, 3).reshape(-1, 3)
            v_cart = v_mw / np.sqrt(m3)
            v_cart /= np.linalg.norm(v_cart)
            disp0 = amp_bohr * mass_scale * v_cart

            self.geom.coords = (ref +  disp0).reshape(-1)
            E_plus  = _calc_energy(self.geom, self.calc_kwargs)
            self.geom.coords = (ref -  disp0).reshape(-1)
            E_minus = _calc_energy(self.geom, self.calc_kwargs)
            self.geom.coords = ref.reshape(-1)

            disp_total += (disp0 if E_plus <= E_minus else -disp0)

        del modes
        _clear_cuda_cache()

        self.geom.coords = (ref + disp_total).reshape(-1)
        return True

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

        # Mass scaling (carbon moves exactly flatten_amp_ang Å)
        mass_scale = np.sqrt(12.011 / self.masses_amu)[:, None]
        amp_bohr = self.flatten_amp_ang / BOHR2ANG

        # Get reference energy
        E_ref = _calc_energy(self.geom, self.calc_kwargs)

        # Apply modes sequentially (like pdb2reaction)
        for idx in targets:
            v_mw = modes[idx].detach().cpu().numpy().reshape(-1, 3)
            m3 = np.repeat(self.masses_amu, 3).reshape(-1, 3)
            v_cart = v_mw / np.sqrt(m3)
            v_cart /= np.linalg.norm(v_cart)

            disp = amp_bohr * mass_scale * v_cart
            ref = self.geom.coords.reshape(-1, 3)

            plus = ref + disp
            minus = ref - disp

            self.geom.coords = plus.reshape(-1)
            E_plus = _calc_energy(self.geom, self.calc_kwargs)

            self.geom.coords = minus.reshape(-1)
            E_minus = _calc_energy(self.geom, self.calc_kwargs)

            # Keep lower-energy side and continue from there
            use_plus = E_plus <= E_minus
            self.geom.coords = (plus if use_plus else minus).reshape(-1)
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
        H_t = _calc_full_hessian_torch(self.geom, self.calc_kwargs_partial, self.device)
        coords_bohr_t = torch.as_tensor(self.geom.coords.reshape(-1, 3),
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
                root=self.root, freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None
            )
        else:
            click.echo("[tsopt] Using active-block Hessian from UMA (partial Hessian). Skip full-space TR check.")
            mode_xyz = _mode_direction_by_root_from_Hact(
                H_t, self.geom.coords.reshape(-1, 3), self.geom.atomic_numbers,
                self.masses_au_t, active_idx, self.device, root=self.root
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

        _, zero_step_loose = self._dimer_loop(self.thresh_loose)

        zero_step_normal = False
        if (self.max_total_cycles - self._cycles_spent) > 0:
            # (3) Update mode & normal loop (reuse Hessian if 0-step converged)
            H_t = self._calc_full_hessian_cached(self.calc_kwargs_partial, allow_reuse=zero_step_loose)
            coords_bohr_t = torch.as_tensor(self.geom.coords.reshape(-1, 3),
                                            dtype=H_t.dtype, device=H_t.device)
            if H_t.size(0) == 3 * N:
                click.echo("[tsopt] TR-projection residual check skipped to conserve VRAM.")
                mode_xyz = _mode_direction_by_root(
                    H_t, coords_bohr_t, self.masses_au_t,
                    root=self.root, freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None
                )
            else:
                click.echo("[tsopt] Using active-block Hessian from UMA (partial Hessian). Skip full-space TR check.")
                active_idx, mask_dof = self._resolve_hessian_active_subspace(H_t, N)
                mode_xyz = _mode_direction_by_root_from_Hact(
                    H_t, self.geom.coords.reshape(-1, 3), self.geom.atomic_numbers,
                    self.masses_au_t, active_idx, self.device, root=self.root
                )
            np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")
            del mode_xyz, coords_bohr_t, H_t
            _clear_cuda_cache()

            click.echo("[tsopt] Normal Dimer Loop...")
            _, zero_step_normal = self._dimer_loop(self.thresh)
        else:
            click.echo("[tsopt] Reached --max-cycles budget after loose loop; skipping normal dimer loop.")

        if self.flatten_max_iter > 0 and (self.max_total_cycles - self._cycles_spent) > 0:
            # (4) Flatten Loop — *reduced* exact Hessian calls via Bofill updates (active DOF only)
            click.echo("[tsopt] Flatten Loop with Bofill-updated active Hessian...")

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

            # Gradient & coordinates snapshot for quasi-Newton updates
            x_prev = self.geom.coords.copy().reshape(-1)             # (3N,)
            g_prev = _calc_gradient(self.geom, self.calc_kwargs).reshape(-1)  # (3N,)

            # Flatten iterations with *approximate* Hessian updates
            for it in range(self.flatten_max_iter):
                if (self.max_total_cycles - self._cycles_spent) <= 0:
                    break

                # (a) Estimate current imaginary modes using the *active* Hessian
                freqs_est = _frequencies_from_Hact(H_act, self.geom.atomic_numbers,
                                                   self.geom.coords.reshape(-1, 3), active_idx, self.device)
                n_imag = int(np.sum(freqs_est < -abs(self.neg_freq_thresh_cm)))
                click.echo(f"[tsopt] n≈{n_imag}  (approx imag: {[float(x) for x in freqs_est if x < -abs(self.neg_freq_thresh_cm)]})")
                if n_imag <= 1:
                    break

                # (b) Get approximate modes for flattening (embedded, mass-weighted)
                freqs_cm_approx, modes_embedded = _modes_from_Hact_embedded(
                    H_act, self.geom.atomic_numbers, self.geom.coords.reshape(-1, 3), active_idx, self.device
                )

                # (c) Do flatten step using the approximate modes
                x_before_flat = self.geom.coords.copy().reshape(-1)
                did_flatten = self._flatten_once_with_modes(freqs_cm_approx, modes_embedded)
                # Free GPU tensors from mode computation immediately after use
                del freqs_cm_approx, modes_embedded
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
                if not did_flatten:
                    break
                x_after_flat = self.geom.coords.copy().reshape(-1)

                # (d) Bofill update using UMA gradients across the flatten displacement
                g_after_flat = _calc_gradient(self.geom, self.calc_kwargs).reshape(-1)
                delta_flat_full = x_after_flat - x_before_flat
                delta_flat_act = delta_flat_full[mask_dof]
                g_old_act = g_prev[mask_dof]
                g_new_act = g_after_flat[mask_dof]
                H_act = _bofill_update_active(H_act, delta_flat_act, g_new_act, g_old_act)

                # (e) Refresh dimer direction from updated active Hessian
                mode_xyz = _mode_direction_by_root_from_Hact(
                    H_act, self.geom.coords.reshape(-1, 3), self.geom.atomic_numbers,
                    self.masses_au_t, active_idx, self.device, root=self.root
                )
                np.savetxt(self.mode_path, mode_xyz, fmt="%.12f")
                del mode_xyz

                # (f) Re-optimize with Dimer (consumes global cycle budget)
                # Clear VRAM before dimer loop to ensure space for Hessian recomputation
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
                _, zero_step_flat = self._dimer_loop(self.thresh)

                # (g) Bofill update again across the optimization displacement
                x_after_opt = self.geom.coords.copy().reshape(-1)
                g_after_opt = _calc_gradient(self.geom, self.calc_kwargs).reshape(-1)
                delta_opt_full = x_after_opt - x_after_flat
                delta_opt_act = delta_opt_full[mask_dof]
                g_old_act2 = g_after_flat[mask_dof]
                g_new_act2 = g_after_opt[mask_dof]
                H_act = _bofill_update_active(H_act, delta_opt_act, g_new_act2, g_old_act2)

                # (h) Prepare for next iteration
                x_prev = x_after_opt
                g_prev = g_after_opt
        elif self.flatten_max_iter > 0:
            click.echo("[tsopt] Reached --max-cycles budget; skipping flatten loop.")

        # (5) Final outputs
        final_xyz = self.out_dir / "final_geometry.xyz"
        atoms_final = Atoms(self.geom.atoms, positions=(self.geom.coords.reshape(-1, 3) * BOHR2ANG), pbc=False)
        write(final_xyz, atoms_final)

        # Final Hessian → imaginary mode animation
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
        if H_t.size(0) == 3 * N:
            freqs_cm, modes = _frequencies_cm_and_modes(
                H_t, self.geom.atomic_numbers, self.geom.coords.reshape(-1, 3), self.device,
                freeze_idx=self.freeze_atoms if len(self.freeze_atoms) > 0 else None
            )
        else:
            active_idx_final, _ = self._resolve_hessian_active_subspace(H_t, N)
            freqs_cm, modes = _modes_from_Hact_embedded(
                H_t, self.geom.atomic_numbers, self.geom.coords.reshape(-1, 3),
                active_idx_final, self.device
            )

        del H_t
        del H_final_reuse_cpu, H_final_reuse_coords
        n_written = _write_all_imag_modes(
            self.geom,
            freqs_cm,
            modes,
            self.neg_freq_thresh_cm,
            self.vib_dir,
        )
        if n_written == 0:
            click.echo(
                "[tsopt] No imaginary mode found at the end (nu_min = %.2f cm^-1)." % (float(freqs_cm.min()),),
                err=True,
            )
        else:
            click.echo(f"[tsopt] Wrote {n_written} final imaginary mode(s).")
        del modes, freqs_cm

        _clear_cuda_cache()
        click.echo(f"[tsopt] Saved final geometry → {final_xyz}")
        click.echo(f"[tsopt] Mode files → {self.vib_dir}")


# ===================================================================
#                         Defaults for CLI
# ===================================================================

# Configuration defaults (imported from defaults.py)
GEOM_KW: Dict[str, Any] = deepcopy(GEOM_KW_DEFAULT)
CALC_KW: Dict[str, Any] = deepcopy(MLMM_CALC_KW)

# HessianDimer defaults - combine imported DIMER_KW and HESSIAN_DIMER_KW
hessian_dimer_KW = {
    **HESSIAN_DIMER_KW,
    "dimer": {**DIMER_KW},
    "lbfgs": {**LBFGS_KW},
}

# ===================================================================
#                            CLI
# ===================================================================

@click.command(
    help="TS optimization: light (Dimer) or heavy (RS-I-RFO) for the ML/MM calculator.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Starting geometry (PDB or XYZ). XYZ provides higher coordinate precision. "
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
    "--real-parm7",
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
    help="PDB containing the ML-region atoms. Optional when --detect-layer is enabled.",
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
    help="Total charge of the ML region.",
)
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=None,
    show_default=False,
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
    default=0.0,
    show_default=True,
    help="Distance cutoff (Å) from ML region for MM atoms to include in Hessian calculation. "
         "Applied to movable MM atoms. Default 0.0 means ML-only partial Hessian.",
)
@click.option(
    "--movable-cutoff",
    "movable_cutoff",
    type=float,
    default=None,
    show_default=False,
    help="Distance cutoff (Å) from ML region for movable MM atoms. "
         "MM atoms beyond this are frozen. "
         "Providing --movable-cutoff disables --detect-layer.",
)
@click.option(
    "--hessian-calc-mode",
    type=click.Choice(["Analytical", "FiniteDifference"], case_sensitive=False),
    default=None,
    help="How UMA builds the ML Hessian (Analytical or FiniteDifference); "
         "overrides calc.hessian_calc_mode from YAML.",
)
@click.option("--max-cycles", type=int, default=10000, show_default=True, help="Maximum total LBFGS cycles.")
@click.option(
    "--dump/--no-dump",
    default=False,
    show_default=True,
    help="Write concatenated trajectory 'optimization_all_trj.xyz'.",
)
@click.option("--out-dir", type=str, default=OUT_DIR_TSOPT, show_default=True, help="Output directory.")
@click.option(
    "--thresh",
    type=str,
    default=None,
    help="Convergence preset (gau_loose|gau|gau_tight|gau_vtight|baker|never).",
)
@click.option(
    "--opt-mode",
    type=click.Choice(["light", "heavy"], case_sensitive=False),
    default="heavy",
    show_default=True,
    help="TS optimizer mode: light (Dimer) or heavy (RS-I-RFO with full Hessian).",
)
@click.option(
    "--partial-hessian-flatten/--full-hessian-flatten",
    "partial_hessian_flatten",
    default=True,
    show_default=True,
    help="Use partial Hessian (ML region only) for imaginary mode detection in flatten loop.",
)
@click.option(
    "--active-dof-mode",
    type=click.Choice(["all", "ml-only", "partial", "unfrozen"], case_sensitive=False),
    default="partial",
    show_default=True,
    help="Active DOF selection for final frequency analysis: "
         "all (all atoms), ml-only (ML only), partial (ML + Hessian-target MM, default), "
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
    spin: Optional[int],
    freeze_atoms_text: Optional[str],
    hess_cutoff: Optional[float],
    movable_cutoff: Optional[float],
    hessian_calc_mode: Optional[str],
    max_cycles: int,
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    opt_mode: str,
    partial_hessian_flatten: bool,
    active_dof_mode: str,
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
) -> None:
    def _is_param_explicit(name: str) -> bool:
        try:
            source = ctx.get_parameter_source(name)
            return source not in (None, ParameterSource.DEFAULT)
        except Exception:
            return False

    config_yaml, override_yaml, used_legacy_yaml = _resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )
    merged_yaml_cfg = _load_merged_yaml_cfg(
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
    source_path = prepared_input.source_path
    charge, spin = resolve_charge_spin_or_raise(prepared_input, charge, spin)

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

    # Resolve optimizer mode (default is now heavy/RS-I-RFO)
    mode_resolved = normalize_choice(
        opt_mode,
        param="--opt-mode",
        alias_groups=TSOPT_MODE_ALIASES,
        allowed_hint="light, heavy",
    )
    use_heavy = (mode_resolved == "heavy")

    config_layer_cfg = load_yaml_dict(config_yaml)
    override_layer_cfg = load_yaml_dict(override_yaml)
    geom_cfg: Dict[str, Any] = deepcopy(GEOM_KW)
    calc_cfg: Dict[str, Any] = deepcopy(CALC_KW)
    opt_cfg: Dict[str, Any] = dict(OPT_BASE_KW)
    lbfgs_cfg: Dict[str, Any] = dict(LBFGS_KW)
    simple_cfg: Dict[str, Any] = dict(hessian_dimer_KW)
    rsirfo_cfg: Dict[str, Any] = dict(RSIRFO_KW)

    apply_yaml_overrides(
        config_layer_cfg,
        [
            (geom_cfg, (("geom",),)),
            (calc_cfg, (("calc",), ("mlmm",))),
            (opt_cfg, (("opt",),)),
            (simple_cfg, (("hessian_dimer",),)),
            (rsirfo_cfg, (("rsirfo",),)),
        ],
    )
    if _is_param_explicit("hessian_calc_mode") and hessian_calc_mode is not None:
        calc_cfg["hessian_calc_mode"] = str(hessian_calc_mode)
    if _is_param_explicit("max_cycles"):
        opt_cfg["max_cycles"] = int(max_cycles)
    if _is_param_explicit("dump"):
        opt_cfg["dump"] = bool(dump)
    if _is_param_explicit("out_dir"):
        opt_cfg["out_dir"] = out_dir
    if _is_param_explicit("thresh") and thresh is not None:
        opt_cfg["thresh"] = str(thresh)
        simple_cfg["thresh"] = str(thresh)
        rsirfo_cfg["thresh"] = str(thresh)
    if _is_param_explicit("detect_layer"):
        calc_cfg["use_bfactor_layers"] = bool(detect_layer)
    if _is_param_explicit("hess_cutoff") and hess_cutoff is not None:
        calc_cfg["hess_cutoff"] = float(hess_cutoff)
    if _is_param_explicit("movable_cutoff") and movable_cutoff is not None:
        calc_cfg["movable_cutoff"] = float(movable_cutoff)
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
    calc_cfg["input_pdb"] = str(source_path)
    calc_cfg["real_parm7"] = str(real_parm7)

    apply_yaml_overrides(
        override_layer_cfg,
        [
            (geom_cfg, (("geom",),)),
            (calc_cfg, (("calc",), ("mlmm",))),
            (opt_cfg, (("opt",),)),
            (simple_cfg, (("hessian_dimer",),)),
            (rsirfo_cfg, (("rsirfo",),)),
        ],
    )

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

    # Propagate opt.print_every to TS optimizers (LBFGS in light mode, RS-I-RFO in heavy mode).
    if "print_every" in opt_cfg:
        try:
            pe = int(opt_cfg["print_every"])
            if pe >= 1:
                simple_cfg.setdefault("lbfgs", {})
                simple_cfg["lbfgs"]["print_every"] = pe
                rsirfo_cfg["print_every"] = pe
        except Exception:
            pass

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
                    "optimizer_mode": "heavy-rsirfo" if use_heavy else "light-dimer",
                    "detect_layer": bool(detect_layer_enabled),
                    "model_region_source": model_region_source,
                    "model_indices_count": 0 if not model_indices else len(model_indices),
                    "hessian_calc_mode": calc_cfg.get("hessian_calc_mode"),
                    "partial_hessian_flatten": bool(partial_hessian_flatten),
                    "active_dof_mode": str(active_dof_mode),
                    "will_run_tsopt": True,
                    "will_write_summary": True,
                },
            )
        )
        click.echo("[dry-run] Validation complete. TS optimization execution was skipped.")
        prepared_input.cleanup()
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

    # Pretty-print config summary (only non-default values for concise logging)
    click.echo(f"\n[mode] TS Optimizer: {'RS-I-RFO (heavy)' if use_heavy else 'Dimer (light)'}\n")
    click.echo(pretty_block("geom", format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")))
    echo_calc = strip_inherited_keys(calc_cfg, CALC_KW, mode="same")
    echo_calc = format_freeze_atoms_for_echo(echo_calc, key="freeze_atoms")
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

    # --------------------------
    # 2) Prepare geometry dir
    # --------------------------
    out_dir_path.mkdir(parents=True, exist_ok=True)

    # --------------------------
    # 3) Run
    # --------------------------
    try:
        if use_heavy:
            # Heavy mode: RS-I-RFO with full Hessian
            click.echo("\n=== TS optimization (RS-I-RFO heavy mode) started ===\n")
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

            base_calc = mlmm(**calc_cfg)
            geometry.set_calculator(base_calc)

            rsirfo_args = {**rsirfo_cfg}
            rsirfo_args["out_dir"] = str(out_dir_path)
            rsirfo_args["max_cycles"] = int(opt_cfg["max_cycles"])
            rsirfo_args["dump"] = bool(opt_cfg["dump"])
            if thresh is not None:
                rsirfo_args["thresh"] = str(thresh)

            # If non-Hessian movable atoms exist, use partial-Hessian microiterations.
            calc_core = base_calc.core if hasattr(base_calc, "core") else base_calc
            movable_mm_indices = list(getattr(calc_core, "movable_mm_indices", []))
            hess_active_atoms = list(getattr(calc_core, "hess_active_atoms", []))
            use_partial_micro = bool(movable_mm_indices) and bool(hess_active_atoms)

            if use_partial_micro:
                layeropt_cfg = dict(LAYEROPT_KW)
                outer_lbfgs_kwargs = {**lbfgs_cfg, **opt_cfg}
                outer_lbfgs_kwargs["out_dir"] = str(out_dir_path)
                outer_lbfgs_kwargs["thresh"] = layeropt_cfg.get("outer_thresh", "gau_loose")
                outer_max_cycles = int(layeropt_cfg.get("outer_max_cycles", 1500))
                outer_max_cycles = min(outer_max_cycles, int(opt_cfg.get("max_cycles", outer_max_cycles)))

                optimizer = PartialHessianMicroIterationOptimizer(
                    geometry=geometry,
                    hess_indices=hess_active_atoms,
                    outer_indices=movable_mm_indices,
                    outer_lbfgs_kwargs=outer_lbfgs_kwargs,
                    inner_rfo_kwargs=rsirfo_args,
                    max_macro_cycles=int(opt_cfg.get("max_cycles", 10000)),
                    outer_max_cycles=outer_max_cycles,
                    outer_thresh=layeropt_cfg.get("outer_thresh", "gau_loose"),
                    echo_fn=click.echo,
                    inner_opt_cls=RSIRFOptimizer,
                    dump=bool(opt_cfg["dump"]),
                )

                # Keep Hessian dimensionality consistent with inner RS-I-RFO active DOF.
                if hasattr(calc_core, "return_partial_hessian"):
                    calc_core.return_partial_hessian = True

                click.echo("\n=== Partial-Hessian MicroIteration (RS-I-RFO) started ===\n")
                optimizer.run()
                click.echo("\n=== Partial-Hessian MicroIteration (RS-I-RFO) finished ===\n")
            else:
                optimizer = RSIRFOptimizer(geometry, **rsirfo_args)
                optimizer.run()
                if bool(opt_cfg["dump"]):
                    _append_xyz_trajectory(optim_all_path, out_dir_path / "optimization_trj.xyz")

            click.echo("\n=== TS optimization (RS-I-RFO heavy mode) finished ===\n")

            # --- Post-RSIRFO: count imaginary modes and optional flatten loop ---
            # Save cycle count before deleting optimizer for budget check.
            _rsirfo_cycles_spent = getattr(optimizer, "cur_cycle", 0) + 1
            geometry.set_calculator(None)
            del optimizer
            del calc_core
            del base_calc
            _clear_cuda_cache()
            mlmm_kwargs_for_heavy = dict(calc_cfg)
            mlmm_kwargs_for_heavy["out_hess_torch"] = True
            if use_partial_micro:
                mlmm_kwargs_for_heavy["return_partial_hessian"] = True
            device = _torch_device(simple_cfg.get("device", calc_cfg.get("ml_device", "auto")))

            # Determine active atoms for frequency analysis based on --active-dof-mode.
            active_atoms_freq = _get_active_dof_indices(
                calc_cfg, len(geometry.atomic_numbers), active_dof_mode, freeze_atoms_final
            )
            n_atoms = len(geometry.atomic_numbers)
            if active_atoms_freq is not None:
                active_set = set(active_atoms_freq)
                freeze_atoms_freq = [i for i in range(n_atoms) if i not in active_set]
            else:
                freeze_atoms_freq = freeze_atoms_final if freeze_atoms_final else None

            def _calc_freqs_and_modes() -> Tuple[np.ndarray, torch.Tensor]:
                H = _calc_full_hessian_torch(geometry, mlmm_kwargs_for_heavy, device)
                n_full = len(geometry.atomic_numbers)
                if H.shape[0] != 3 * n_full:
                    # Partial Hessian: use Hessian-target atoms only and embed modes back.
                    active_atoms = None
                    if getattr(geometry, "within_partial_hessian", None) is not None:
                        active_atoms = geometry.within_partial_hessian.get("active_atoms")
                    if active_atoms is None:
                        active_atoms = hess_active_atoms
                    if active_atoms is None:
                        active_atoms = active_atoms_freq
                    if active_atoms is None:
                        active_atoms = []
                    else:
                        active_atoms = [int(i) for i in np.asarray(active_atoms, dtype=int).reshape(-1).tolist()]
                    if not active_atoms:
                        raise RuntimeError(
                            "No active atoms available for partial Hessian frequency analysis."
                        )
                    coords_act = geometry.cart_coords.reshape(-1, 3)[active_atoms]
                    nums_act = np.asarray(geometry.atomic_numbers)[active_atoms]
                    freqs_local, modes_act = _frequencies_cm_and_modes(
                        H,
                        nums_act,
                        coords_act,
                        device,
                        freeze_idx=None,
                    )
                    # Embed modes into full 3N on CPU to reduce VRAM peak.
                    modes_local = torch.zeros((modes_act.shape[0], 3 * n_full), dtype=modes_act.dtype, device="cpu")
                    mask_dof = torch.as_tensor(_mask_dof_from_active_idx(n_full, active_atoms), dtype=torch.bool)
                    modes_local[:, mask_dof] = modes_act.detach().cpu()
                    del coords_act, nums_act, modes_act, mask_dof
                else:
                    freqs_local, modes_gpu = _frequencies_cm_and_modes(
                        H,
                        geometry.atomic_numbers,
                        geometry.cart_coords.reshape(-1, 3),
                        device,
                        freeze_idx=freeze_atoms_freq,
                    )
                    modes_local = modes_gpu.detach().cpu()
                    del modes_gpu
                del H
                _clear_cuda_cache()
                return freqs_local, modes_local

            try:
                freqs_cm, modes = _calc_freqs_and_modes()
            except Exception as exc:
                is_oom = isinstance(exc, torch.OutOfMemoryError) or ("cuda out of memory" in str(exc).lower())
                if is_oom:
                    click.echo(
                        "[tsopt] WARNING: CUDA OOM during final frequency analysis; "
                        "skipping imaginary-mode analysis/flatten loop.",
                        err=True,
                    )
                    _clear_cuda_cache()
                    freqs_cm, modes = None, None
                else:
                    raise
            neg_freq_thresh_cm = float(simple_cfg.get("neg_freq_thresh_cm", 5.0))

            if freqs_cm is not None and modes is not None:
                neg_mask = freqs_cm < -abs(neg_freq_thresh_cm)
                n_imag = int(np.sum(neg_mask))
                ims = [float(x) for x in freqs_cm if x < -abs(neg_freq_thresh_cm)]
                click.echo(f"[Imaginary modes] n={n_imag}  ({ims})")

                flatten_max_iter = int(simple_cfg.get("flatten_max_iter", 0))
                user_max_cycles = int(opt_cfg.get("max_cycles", 10000))
                budget_remaining = user_max_cycles - _rsirfo_cycles_spent > 0

                if flatten_max_iter > 0 and n_imag > 1 and not budget_remaining:
                    click.echo("[tsopt] Reached --max-cycles budget; skipping flatten loop.")
                elif flatten_max_iter > 0 and n_imag > 1 and budget_remaining:
                    click.echo("[flatten] Extra imaginary modes detected; starting RS-I-RFO flatten loop.")
                    masses_amu = np.array([atomic_masses[z] for z in geometry.atomic_numbers])
                    main_root = int(simple_cfg.get("root", 0))

                    for it in range(flatten_max_iter):
                        click.echo(f"[flatten] RS-I-RFO iteration {it + 1}/{flatten_max_iter}")
                        did_flatten = _flatten_once_with_modes_for_geom(
                            geometry,
                            masses_amu,
                            mlmm_kwargs_for_heavy,
                            freqs_cm,
                            modes,
                            neg_freq_thresh_cm,
                            float(simple_cfg.get("flatten_amp_ang", 0.10)),
                            float(simple_cfg.get("flatten_sep_cutoff", 2.0)),
                            int(simple_cfg.get("flatten_k", 10)),
                            main_root,
                        )
                        if not did_flatten:
                            click.echo("[flatten] No eligible modes to flatten; stopping.")
                            break

                        del freqs_cm, modes
                        _clear_cuda_cache()
                        base_calc = mlmm(**calc_cfg)
                        geometry.set_calculator(base_calc)
                        optimizer = RSIRFOptimizer(geometry, **rsirfo_args)
                        click.echo("\n=== TS optimization (RS-I-RFO) restarted ===\n")
                        optimizer.run()
                        click.echo("\n=== TS optimization (RS-I-RFO) finished ===\n")
                        if bool(opt_cfg["dump"]):
                            _append_xyz_trajectory(optim_all_path, out_dir_path / "optimization_trj.xyz")
                        geometry.set_calculator(None)
                        del optimizer, base_calc
                        _clear_cuda_cache()

                        try:
                            freqs_cm, modes = _calc_freqs_and_modes()
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
                                break
                            raise

                        neg_mask = freqs_cm < -abs(neg_freq_thresh_cm)
                        n_imag = int(np.sum(neg_mask))
                        ims = [float(x) for x in freqs_cm if x < -abs(neg_freq_thresh_cm)]
                        click.echo(f"[Imaginary modes] n={n_imag}  ({ims})")
                        if n_imag <= 1:
                            break

            if freqs_cm is not None and modes is not None:
                # --- Write all final imaginary modes like light mode ---
                vib_dir = out_dir_path / "vib"
                vib_dir.mkdir(parents=True, exist_ok=True)
                n_written = _write_all_imag_modes(
                    geometry,
                    freqs_cm,
                    modes,
                    neg_freq_thresh_cm,
                    vib_dir,
                )
                if n_written == 0:
                    click.echo("[INFO] No imaginary mode found at the end for RS-I-RFO.")
                else:
                    click.echo(f"[DONE] Wrote {n_written} final imaginary mode(s).")
                    click.echo(f"[DONE] Mode files → {vib_dir}")
            else:
                click.echo("[INFO] Skipped final imaginary-mode export due to frequency-analysis fallback.")

            if modes is not None:
                del modes
            if freqs_cm is not None:
                del freqs_cm
            _clear_cuda_cache()

            # Ensure final_geometry.xyz exists (partial-micro path may not write it).
            final_xyz = out_dir_path / "final_geometry.xyz"
            if use_partial_micro or not final_xyz.exists():
                final_xyz.write_text(geometry.as_xyz(), encoding="utf-8")

        else:
            # Light mode: Partial Hessian guided Dimer
            runner = HessianDimer(
                fn=str(geom_input_path),
                out_dir=str(out_dir_path),
                thresh_loose=simple_cfg.get("thresh_loose", "gau_loose"),
                thresh=simple_cfg.get("thresh", "gau"),
                update_interval_hessian=int(simple_cfg.get("update_interval_hessian", 200)),
                neg_freq_thresh_cm=float(simple_cfg.get("neg_freq_thresh_cm", 5.0)),
                flatten_amp_ang=float(simple_cfg.get("flatten_amp_ang", 0.10)),
                flatten_max_iter=int(simple_cfg.get("flatten_max_iter", 20)),
                mem=int(simple_cfg.get("mem", 100000)),
                use_lobpcg=bool(simple_cfg.get("use_lobpcg", True)),
                calc_kwargs=dict(calc_cfg),
                device=str(simple_cfg.get("device", calc_cfg.get("ml_device", "auto"))),
                dump=bool(opt_cfg["dump"]),
                root=int(simple_cfg.get("root", 0)),
                dimer_kwargs=dict(simple_cfg.get("dimer", {})),
                lbfgs_kwargs=dict(simple_cfg.get("lbfgs", {})),
                max_total_cycles=int(opt_cfg["max_cycles"]),
                geom_kwargs=dict(geom_cfg),
                partial_hessian_flatten=partial_hessian_flatten,
                flatten_sep_cutoff=float(simple_cfg.get("flatten_sep_cutoff", 0.0)),
                flatten_k=int(simple_cfg.get("flatten_k", 10)),
                flatten_loop_bofill=bool(simple_cfg.get("flatten_loop_bofill", False)),
            )

            click.echo("\n=== TS optimization (Partial Hessian Dimer) started ===\n")
            runner.run()
            click.echo("\n=== TS optimization (Partial Hessian Dimer) finished ===\n")

        if source_path.suffix.lower() == ".pdb":
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
                        f"[annot]   B-factors set in '{final_pdb}' "
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
                            f"[annot]   B-factors set in '{opt_pdb}' "
                            f"(ML={BFACTOR_ML:.0f}, MovableMM={BFACTOR_MOVABLE_MM:.0f}, "
                            f"FrozenMM={BFACTOR_FROZEN:.0f})."
                        )
                except Exception as e:
                    click.echo(f"[convert] WARNING: Failed to convert optimization trajectory to PDB: {e}", err=True)
        else:
            final_xyz = out_dir_path / "final_geometry.xyz"

        # summary.md and key_* outputs are disabled.
        click.echo(format_elapsed("[time] Elapsed Time for TS Opt", time_start))

    except ZeroStepLength:
        click.echo("ERROR: Proposed step length dropped below the minimum allowed (ZeroStepLength).", err=True)
        sys.exit(2)
    except OptimizationError as e:
        click.echo(f"ERROR: Optimization failed — {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        import traceback
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled error during optimization:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)
    finally:
        prepared_input.cleanup()


# Allow `python -m mlmm.tsopt` direct execution
if __name__ == "__main__":
    cli()
