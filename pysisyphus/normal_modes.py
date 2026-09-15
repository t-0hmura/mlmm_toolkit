# pysisyphus/normal_modes.py

"""Pure normal-mode kernel: mass weighting, rigid projection, diagonalization.

This is a lower bundled-engine module, alongside
``pysisyphus.tr_projection``. It does not import the product layer, so
bundled-engine code (e.g.
``pysisyphus.tsoptimizers.TSHessianOptimizer``) can consume the mass/mode kernel
without an upward product dependency.

The bounded-peak Hessian symmetrizer lives here as well so the module is
self-contained; ``mlmm.core.utils`` re-exports it for backward compatibility.
"""

from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np
import torch
import ase.units as units
from ase.data import atomic_masses

from pysisyphus.constants import BOHR2ANG, AMU2AU, AU2EV
from pysisyphus._array import active_square
from pysisyphus.tr_projection import (
    active_tr_basis,
    compact_project_hessian,
    project_hessian_inplace,
)


DEFAULT_FREQUENCY_ZERO_CUTOFF_CM = 5.00


def normalize_frequency_zero_cutoff_cm(value) -> float:
    """Validate a configurable non-negative zero-mode cutoff."""
    cutoff = float(value)
    if not np.isfinite(cutoff) or cutoff < 0.0:
        raise ValueError("frequency zero cutoff must be finite and non-negative")
    return cutoff


def resolved_frequency_mask(
    freqs_cm, cutoff_cm=DEFAULT_FREQUENCY_ZERO_CUTOFF_CM
) -> np.ndarray:
    """Select modes outside the configured symmetric zero window."""
    cutoff = normalize_frequency_zero_cutoff_cm(cutoff_cm)
    return np.abs(np.asarray(freqs_cm, dtype=float)) > cutoff


def resolved_imaginary_mask(
    freqs_cm, cutoff_cm=DEFAULT_FREQUENCY_ZERO_CUTOFF_CM
) -> np.ndarray:
    """Select resolved imaginary modes using the configured zero window."""
    cutoff = normalize_frequency_zero_cutoff_cm(cutoff_cm)
    return np.asarray(freqs_cm, dtype=float) < -cutoff


def _strict_negative_count(freqs_cm, projection_info) -> Optional[int]:
    """Count all finite negative modes in an explicitly complete PHVA partition.

    Display/eligibility filtering is unchanged. Missing near-zero metadata or
    incomplete/nonfinite partitions provide no strict curvature certificate.
    """
    if not isinstance(projection_info, dict):
        return None
    try:
        frequencies = np.asarray(freqs_cm, dtype=float)
        near = np.asarray(projection_info["near_zero_frequencies_cm"], dtype=float)
        raw = projection_info["raw_mode_count"]
    except (KeyError, TypeError, ValueError):
        return None
    if (
        frequencies.ndim != 1 or near.ndim != 1
        or not np.isfinite(frequencies).all() or not np.isfinite(near).all()
        or isinstance(raw, (bool, np.bool_))
        or not isinstance(raw, (int, np.integer)) or raw < 0
        or frequencies.size + near.size != raw
        or projection_info.get("resolved_mode_count", frequencies.size) != frequencies.size
        or projection_info.get("near_zero_mode_count", near.size) != near.size
    ):
        return None
    return int(np.count_nonzero(frequencies < 0.0) + np.count_nonzero(near < 0.0))


def filter_resolved_modes(
    freqs_cm,
    modes,
    cutoff_cm=DEFAULT_FREQUENCY_ZERO_CUTOFF_CM,
    *,
    filter_info=None,
):
    """Remove |ν| <= cutoff while preserving frequency/mode alignment."""
    frequencies = np.asarray(freqs_cm, dtype=float)
    cutoff = normalize_frequency_zero_cutoff_cm(cutoff_cm)
    keep = resolved_frequency_mask(frequencies, cutoff)
    if filter_info is not None:
        near_zero = frequencies[~keep]
        filter_info.clear()
        filter_info.update(
            {
                "frequency_zero_cutoff_cm": cutoff,
                "raw_mode_count": int(frequencies.size),
                "resolved_mode_count": int(np.count_nonzero(keep)),
                "near_zero_mode_count": int(near_zero.size),
                "near_zero_frequencies_cm": [float(value) for value in near_zero],
            }
        )
    if isinstance(modes, torch.Tensor):
        mode_keep = torch.as_tensor(keep, dtype=torch.bool, device=modes.device)
        filtered_modes = modes[mode_keep]
    else:
        filtered_modes = np.asarray(modes)[keep]
    return frequencies[keep], filtered_modes


def symmetrize_inplace(H, chunk: int = 512):
    """Symmetrize a square Hessian-like tensor in place with bounded peak VRAM.

    Replaces the 2x-peak idiom ``_t = H.T.clone(); H.add_(_t).mul_(0.5); del _t``
    with a chunked average that writes BOTH triangles symmetrically (no
    upper-triangle-only tricks). Peak extra allocation is bounded by
    ``chunk * chunk`` elements (vs ``N * N`` for the naive form).

    Parameters
    ----------
    H : torch.Tensor
        2-D square tensor of any floating dtype (fp32 / fp64) on any device.
        Modified in place; partial Hessian semantics preserved (no shape change).
    chunk : int, optional
        Block edge length for the off-diagonal averaging loop. Peak extra VRAM
        is bounded by ``chunk * chunk * dtype.itemsize`` bytes. Default 512.

    Returns
    -------
    torch.Tensor
        The same ``H`` object, modified in place, for chainability.
    """

    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise ValueError(
            f"symmetrize_inplace expects a square 2-D tensor, got shape {tuple(H.shape)}"
        )
    N = H.shape[0]
    if N == 0:
        return H

    # Degenerate fast path: whole matrix fits inside one chunk — naive in-place
    # average on a single chunk-sized temp (still bounded by chunk*chunk).
    if N <= chunk:
        tmp = H.T.contiguous()
        H.add_(tmp).mul_(0.5)
        del tmp
        return H

    # Chunked loop: bound peak extra alloc to chunk*chunk; writes BOTH triangles.
    for i in range(0, N, chunk):
        ie = min(i + chunk, N)
        diag = H[i:ie, i:ie]
        diag_tmp = diag.T.contiguous()
        # Out-of-place add then assign: an in-place `diag.add_(diag_tmp)` on a
        # strided self-overlapping view raises RuntimeError ("some elements of
        # the input tensor and the written-to tensor refer to a single memory
        # location") for some block sizes. Mirror the off-diagonal assign below.
        H[i:ie, i:ie] = diag.add(diag_tmp).mul(0.5)
        del diag_tmp
        for j in range(ie, N, chunk):
            je = min(j + chunk, N)
            upper = H[i:ie, j:je]
            lower_T = H[j:je, i:ie].T
            avg = upper.add(lower_T).mul_(0.5)
            upper.copy_(avg)
            H[j:je, i:ie].copy_(avg.T)
            del avg
    return H


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


def _mw_projected_hessian(H_t: torch.Tensor,
                          coords_bohr_t: torch.Tensor,
                          masses_au_t: torch.Tensor,
                          tr_projection: str = "constrained",
                          projection_info: Optional[dict] = None,
                          compact: bool = False):
    """
    Project out translations/rotations in mass-weighted space:
    Hmw = M^{-1/2} H M^{-1/2};  P = I - QQ^T;  Hmw_proj = P Hmw P

    To save memory, update **H_t in-place** (no clone) and return it.
    The output is explicitly symmetrized after TR projection.

    With ``compact=True``, return ``(H_compact, lift)`` for the true orthogonal
    complement of the rigid basis instead. The reduced matrix carries exactly
    one root per remaining degree of freedom, so no root has to be discarded by
    magnitude afterwards, and ``lift`` maps reduced eigenvectors back to the
    mass-weighted Cartesian space.
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

        Q, info = active_tr_basis(
            coords_bohr_t,
            masses_au_t,
            list(range(int(coords_bohr_t.shape[0]))),
            mode=tr_projection,
        )
        if compact:
            H_t, lift = compact_project_hessian(H_t, Q)
        else:
            project_hessian_inplace(H_t, Q)
            lift = None
        if projection_info is not None:
            projection_info.clear()
            projection_info.update(info.as_dict())

        # Bounded-peak symmetrization (writes BOTH triangles; peak temp <= chunk^2
        # instead of full N×N clone). Do not rely on upper-triangle-only tricks
        # like `eigh(UPLO='U')` to skip symmetrization.
        symmetrize_inplace(H_t)

        del masses_amu_t, m3, inv_sqrt_m, inv_sqrt_m_col, inv_sqrt_m_row
        del Q

        if torch.cuda.is_available() and device.type == "cuda":
            torch.cuda.empty_cache()
        return (H_t, lift) if compact else H_t


# CHEMISTRY-RULE:6 PHVA + MLIP active-block: mass-weighted Hessian only;
# TR projection is applied separately downstream.
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


# Keep the PHVA/MLIP active-block derivation in the function docstring beside
# the numerical implementation.
def _frequencies_cm_and_modes(H_t: torch.Tensor,
                              atomic_numbers: List[int],
                              coords_bohr: np.ndarray,
                              device: torch.device,
                              freeze_idx: Optional[List[int]] = None,
                              tr_projection: str = "constrained",
                              projection_info: Optional[dict] = None,
                              frequency_zero_cutoff_cm: float = DEFAULT_FREQUENCY_ZERO_CUTOFF_CM) -> Tuple[np.ndarray, torch.Tensor]:
    """
    Diagonalize a (possibly PHVA/active-subspace) TR-projected mass-weighted Hessian
    to obtain frequencies (cm^-1) and mass-weighted eigenvectors (modes).

    If `freeze_idx` is provided (list of 0-based atom indices), perform
    Partial Hessian Vibrational Analysis (PHVA). Supports two cases:

      A) Full Hessian given (3N×3N):
         1) build Hmw = M^{-1/2} H M^{-1/2}
         2) take the active subspace by removing DOF of frozen atoms
         3) remove only constrained-system rigid null modes represented in the active subspace
         4) diagonalize and embed eigenvectors back to 3N by zero-filling frozen DOF

      B) Already-reduced (active-block) Hessian given (3N_act×3N_act), e.g.
         when UMA is called with return_partial_hessian=True:
         1) mass-weight with **active** masses only
         2) apply the same constrained-system rigid-null treatment in active space
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
        masses_amu = _safe_masses_amu(Z)
        masses_au_t = torch.as_tensor(masses_amu * AMU2AU, dtype=H_t.dtype, device=device)
        coords_bohr_t = torch.as_tensor(coords_bohr.reshape(-1, 3), dtype=H_t.dtype, device=device)

        # PHVA path (active DOF subspace with TR-proj)
        if freeze_idx is not None and len(freeze_idx) > 0:
            # Active atom indices
            frozen_set = set(int(i) for i in freeze_idx if 0 <= int(i) < N)
            active_idx = [i for i in range(N) if i not in frozen_set]
            n_active = len(active_idx)
            if n_active == 0:
                raise ValueError("PHVA requires at least one active atom")

            # Determine whether the provided Hessian is already the active block (3N_act×3N_act).
            expected_act_dim = 3 * n_active
            is_partial = (H_t.shape[0] == expected_act_dim and H_t.shape[1] == expected_act_dim)
            is_full = (H_t.shape[0] == 3 * N and H_t.shape[1] == 3 * N)
            if not (is_partial or is_full):
                raise ValueError(
                    "Hessian shape is inconsistent with the full and active PHVA spaces: "
                    f"got {tuple(H_t.shape)}, expected {(3 * N, 3 * N)} or "
                    f"{(expected_act_dim, expected_act_dim)}"
                )
            Q, info = active_tr_basis(
                coords_bohr_t,
                masses_au_t,
                active_idx,
                mode=tr_projection,
            )
            if projection_info is not None:
                projection_info.clear()
                projection_info.update(info.as_dict())

            if is_partial:
                # --- Case B: Active-subspace Hessian supplied ---
                # Mass-weight using only active atoms → project TR modes in the active space
                # → diagonalize → embed back into the full space.
                masses_act = masses_au_t[active_idx]
                # in-place mass-weight (active masses)
                Hmw_act = _mass_weighted_hessian(H_t, masses_act)
                # Diagonalize the orthogonal complement of the rigid basis, so
                # exactly the constrained rigid rank is removed and every genuine
                # low-magnitude root survives for the imaginary count and for
                # thermochemistry.
                Hmw_act, lift = compact_project_hessian(Hmw_act, Q)

                # Bounded-peak symmetrization (helper writes both triangles).
                symmetrize_inplace(Hmw_act)
                omega2, Vsub = torch.linalg.eigh(Hmw_act)

                # Free the (only) Hessian ASAP
                del Hmw_act
                del H_t
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

                if lift is not None:
                    # Lift the reduced eigenvectors back to the active DOF space.
                    Vsub = lift.T @ Vsub  # (3N_act, nmode)
                    del lift

                # Embed to full 3N (mass-weighted eigenvectors)
                modes = torch.zeros((Vsub.shape[1], 3 * N), dtype=Vsub.dtype, device=Vsub.device)
                mask_dof = torch.ones(3 * N, dtype=torch.bool, device=Vsub.device)
                for i in frozen_set:
                    mask_dof[3 * i:3 * i + 3] = False
                modes[:, mask_dof] = Vsub.T
                del Q, mask_dof

            else:
                # --- Case A: Full Hessian (3N×3N) supplied ---
                # Apply full mass-weighting → extract the active block → project TR modes in the active space.
                H_t = _mass_weighted_hessian(H_t, masses_au_t)

                # Build active mask (boolean) and immediately carve out the active block
                mask_dof = torch.ones(3 * N, dtype=torch.bool, device=H_t.device)
                for i in frozen_set:
                    mask_dof[3 * i:3 * i + 3] = False

                # Create the reduced Hessian; free the full one immediately to keep only one in VRAM
                active_dof = torch.nonzero(mask_dof, as_tuple=False).flatten()
                H_act = active_square(H_t, active_dof)
                del active_dof
                del H_t
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
                H_t = H_act
                del H_act

                # Same compact complement treatment as in the active-block case.
                H_t, lift = compact_project_hessian(H_t, Q)

                # Bounded-peak symmetrization (helper writes both triangles).
                symmetrize_inplace(H_t)
                omega2, Vsub = torch.linalg.eigh(H_t)

                # Free the (only) Hessian ASAP
                del H_t
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

                if lift is not None:
                    Vsub = lift.T @ Vsub  # (3N_act, nmode)
                    del lift

                modes = torch.zeros((Vsub.shape[1], 3 * N), dtype=Vsub.dtype, device=Vsub.device)
                modes[:, mask_dof] = Vsub.T  # (nmode, 3N_act) → place into active DOF
                del Vsub, mask_dof, Q

        else:
            H_t, lift = _mw_projected_hessian(
                H_t,
                coords_bohr_t,
                masses_au_t,
                tr_projection=tr_projection,
                projection_info=projection_info,
                compact=True,
            )
            omega2, V = torch.linalg.eigh(H_t)

            # Free the (only) Hessian ASAP
            del H_t
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

            if lift is not None:
                V = lift.T @ V
                del lift
            modes = V.T
            del V

        # Convert to frequencies (cm^-1)
        s_new = (units._hbar * 1e10 / np.sqrt(units._e * units._amu) * np.sqrt(AU2EV) / BOHR2ANG)
        hnu = s_new * torch.sqrt(torch.abs(omega2))
        hnu = torch.where(omega2 < 0, -hnu, hnu)
        freqs_cm = (hnu / units.invcm).detach().cpu().numpy()
        zero_filter_info = {}
        freqs_cm, modes = filter_resolved_modes(
            freqs_cm,
            modes,
            frequency_zero_cutoff_cm,
            filter_info=zero_filter_info,
        )
        if projection_info is not None:
            projection_info.update(zero_filter_info)

        del omega2, hnu
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
