"""Rigid-body projection for unconstrained and Cartesian-constrained systems."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import torch


TR_PROJECTION_MODES = ("constrained", "legacy-active")

# Single source of truth for the default rigid-mode treatment.
#
# ``constrained`` projects only an actual null space of the constrained problem:
# it builds the full-system rigid basis, keeps the part that leaves the frozen
# atoms fixed, and skips the projection entirely when that rank is 0 -- which it
# is for every real frozen active-site cluster.
#
# ``legacy-active`` is the superseded treatment: it projects the ACTIVE
# FRAGMENT's own translations/rotations out of the Hessian.  With frozen atoms
# those are not null modes (K_BB*t_B = -K_BC*t_C != 0), so it is a Rayleigh-Ritz
# compression onto a different, more constrained model -- structurally biased
# toward a SMALLER n_imag, i.e. it can hide a real imaginary mode but never
# create one.  It is retained ONLY as a named opt-in, for reproducing output
# produced before the fix.  Do not make it the default again.
DEFAULT_TR_PROJECTION = "constrained"


def allows_saddle_certification(
    mode: str | None,
    frozen_atoms: Sequence[int],
) -> bool:
    """Return whether this projection can certify a frozen-system saddle."""
    value = DEFAULT_TR_PROJECTION if mode is None else str(mode).strip().lower()
    return value != "legacy-active" or len(tuple(frozen_atoms)) == 0


@dataclass(frozen=True)
class TRProjectionInfo:
    """Provenance for a rigid-body projection in the active Cartesian space."""

    treatment: str
    effective_rank: int
    full_rigid_rank: int
    frozen_constraint_rank: int
    svd_rtol: float
    active_atoms: tuple[int, ...]
    frozen_atoms: tuple[int, ...]

    def as_dict(self) -> dict[str, object]:
        algorithm = (
            "constrained-rigid-null-v1"
            if self.treatment == "constrained"
            else "legacy-active-fragment-v1"
        )
        return {
            "treatment": self.treatment,
            "algorithm": algorithm,
            "effective_rank": self.effective_rank,
            "full_rigid_rank": self.full_rigid_rank,
            "frozen_constraint_rank": self.frozen_constraint_rank,
            "svd_rtol": self.svd_rtol,
            "active_atom_count": len(self.active_atoms),
            "frozen_atom_count": len(self.frozen_atoms),
            "active_atoms": list(self.active_atoms),
            "frozen_atoms": list(self.frozen_atoms),
        }


def normalize_tr_projection_mode(mode: str | None) -> str:
    """Return a validated projection mode name."""

    value = DEFAULT_TR_PROJECTION if mode is None else str(mode).strip().lower()
    if value not in TR_PROJECTION_MODES:
        choices = ", ".join(TR_PROJECTION_MODES)
        raise ValueError(f"Unknown TR projection mode {mode!r}; choose one of: {choices}.")
    if value == "legacy-active":
        import warnings

        warnings.warn(
            "tr_projection='legacy-active' is deprecated and scheduled for removal: it "
            "projects the active fragment's own rigid modes, which is not a null space of "
            "the frozen system and can hide a real imaginary mode on frozen-boundary "
            "systems. Use the default 'constrained' treatment; do not use 'legacy-active' "
            "for pass/HOSP transition-state certification. To reproduce pre-fix results "
            "bitwise, install the pinned pre-fix release instead.",
            DeprecationWarning,
            stacklevel=2,
        )
    return value


def _rigid_basis(coords: torch.Tensor, masses: torch.Tensor) -> torch.Tensor:
    """Return the six mass-weighted rigid-motion columns for ``coords``."""

    coords = coords.reshape(-1, 3)
    masses = masses.reshape(-1).to(dtype=coords.dtype, device=coords.device)
    if coords.shape[0] != masses.numel():
        raise ValueError("coords and masses must describe the same number of atoms")
    if coords.shape[0] == 0:
        return torch.zeros((0, 0), dtype=coords.dtype, device=coords.device)
    if not bool(torch.all(masses > 0)):
        raise ValueError("all atomic masses must be positive")

    sqrt_m = torch.sqrt(masses).reshape(-1, 1)
    center = (masses[:, None] * coords).sum(dim=0) / masses.sum()
    centered = coords - center
    eye = torch.eye(3, dtype=coords.dtype, device=coords.device)
    columns = [(eye[i].expand_as(centered) * sqrt_m).reshape(-1, 1) for i in range(3)]
    columns.extend(
        (torch.cross(centered, eye[i].expand_as(centered), dim=1) * sqrt_m).reshape(-1, 1)
        for i in range(3)
    )
    return torch.cat(columns, dim=1)


def _orthonormal_columns(
    matrix: torch.Tensor,
    *,
    rtol: float,
) -> tuple[torch.Tensor, int]:
    """Return an orthonormal basis for the numerical column space."""

    if matrix.numel() == 0 or matrix.shape[1] == 0:
        return matrix.new_zeros((matrix.shape[0], 0)), 0
    u, singular, _ = torch.linalg.svd(matrix, full_matrices=False)
    scale = float(singular[0].item()) if singular.numel() else 0.0
    if scale == 0.0:
        return matrix.new_zeros((matrix.shape[0], 0)), 0
    rank = int((singular > float(rtol) * scale).sum().item())
    return u[:, :rank], rank


def active_tr_basis(
    coords: torch.Tensor,
    masses: torch.Tensor,
    active_atoms: Sequence[int] | torch.Tensor,
    *,
    mode: str = "constrained",
    rtol: float = 1.0e-10,
) -> tuple[torch.Tensor, TRProjectionInfo]:
    """Build rigid null modes represented on the active Cartesian DOFs.

    ``constrained`` retains only full-system rigid motions that leave every
    inactive atom fixed.  With zero, one, two, or at least three non-collinear
    frozen anchors this has the generic ranks 6, 3, 1, and 0.  ``legacy-active``
    is an isolated-active comparison treatment: it removes the active
    fragment's own translations and rotations using this common kernel.
    """

    mode = normalize_tr_projection_mode(mode)
    if not np.isfinite(rtol) or float(rtol) <= 0.0:
        raise ValueError("rtol must be a finite positive number")
    coords = coords.reshape(-1, 3)
    masses = masses.reshape(-1).to(device=coords.device)
    # The rigid problem has at most six columns, so carrying it out in float64
    # is effectively free and prevents the allowed two-anchor rotation from
    # disappearing under float32 roundoff.
    coords = coords.to(dtype=torch.float64)
    masses = masses.to(dtype=torch.float64)
    n_atoms = int(coords.shape[0])
    if masses.numel() != n_atoms:
        raise ValueError("coords and masses must describe the same number of atoms")

    active = torch.as_tensor(active_atoms, dtype=torch.long, device=coords.device).reshape(-1)
    if active.numel():
        if int(active.min().item()) < 0 or int(active.max().item()) >= n_atoms:
            raise ValueError("active atom index is outside the coordinate array")
        if torch.unique(active).numel() != active.numel():
            raise ValueError("active atom indices must be unique")
    active_tuple = tuple(int(i) for i in active.detach().cpu().tolist())
    active_set = set(active_tuple)
    frozen_tuple = tuple(i for i in range(n_atoms) if i not in active_set)
    n_active = len(active_tuple)

    full_q, full_rank = _orthonormal_columns(_rigid_basis(coords, masses), rtol=rtol)
    if n_active == 0:
        info = TRProjectionInfo(
            mode, 0, full_rank, full_rank, float(rtol), active_tuple, frozen_tuple
        )
        return coords.new_zeros((0, 0)), info

    if mode == "legacy-active":
        q_active, rank = _orthonormal_columns(
            _rigid_basis(coords.index_select(0, active), masses.index_select(0, active)),
            rtol=rtol,
        )
        info = TRProjectionInfo(
            mode, rank, full_rank, 0, float(rtol), active_tuple, frozen_tuple
        )
        return q_active, info

    atom_mask = torch.zeros(n_atoms, dtype=torch.bool, device=coords.device)
    atom_mask[active] = True
    frozen_dofs = (~atom_mask)[:, None].expand(-1, 3).reshape(-1)
    active_dof_indices = torch.as_tensor(
        [3 * atom + axis for atom in active_tuple for axis in range(3)],
        dtype=torch.long,
        device=coords.device,
    )
    q_active_rows = full_q.index_select(0, active_dof_indices)

    if not bool(frozen_dofs.any()):
        frozen_rank = 0
        q_active, rank = _orthonormal_columns(q_active_rows, rtol=rtol)
    else:
        q_frozen_rows = full_q[frozen_dofs]
        if full_rank == 0:
            coefficient_null = full_q.new_zeros((0, 0))
        else:
            _, singular, vh = torch.linalg.svd(q_frozen_rows, full_matrices=True)
            scale = float(singular[0].item()) if singular.numel() else 0.0
            frozen_rank = (
                int((singular > float(rtol) * scale).sum().item())
                if scale > 0.0
                else 0
            )
            coefficient_null = vh[frozen_rank:, :].T
        candidate = q_active_rows @ coefficient_null
        q_active, rank = _orthonormal_columns(candidate, rtol=rtol)

    info = TRProjectionInfo(
        mode,
        rank,
        full_rank,
        frozen_rank,
        float(rtol),
        active_tuple,
        frozen_tuple,
    )
    return q_active, info


def full_cartesian_tr_basis(
    coords: torch.Tensor,
    masses: torch.Tensor,
    active_atoms: Sequence[int] | torch.Tensor,
    *,
    mode: str = "constrained",
    rtol: float = 1.0e-10,
) -> tuple[torch.Tensor, TRProjectionInfo]:
    """Return feasible rigid displacements in full Cartesian coordinates.

    :func:`active_tr_basis` represents rigid motions in the active,
    mass-weighted space used for Hessian diagonalization.  Dimer orientations
    are ordinary Cartesian displacements, so they require the corresponding
    mass-unweighted basis embedded in the full ``3N`` vector.  Frozen rows are
    exactly zero and the returned columns are Euclidean-orthonormal.
    """

    coords = torch.as_tensor(coords).reshape(-1, 3)
    masses = torch.as_tensor(masses, device=coords.device).reshape(-1)
    basis_mw, info = active_tr_basis(
        coords,
        masses,
        active_atoms,
        mode=mode,
        rtol=rtol,
    )
    n_atoms = int(coords.shape[0])
    full = coords.new_zeros((3 * n_atoms, info.effective_rank), dtype=torch.float64)
    if info.effective_rank == 0:
        return full, info

    active = torch.as_tensor(
        info.active_atoms, dtype=torch.long, device=coords.device
    )
    masses64 = masses.to(dtype=torch.float64, device=coords.device)
    active_masses = masses64.index_select(0, active)
    scale = torch.sqrt(torch.repeat_interleave(active_masses, 3)).reshape(-1, 1)
    basis_cart = basis_mw.to(dtype=torch.float64, device=coords.device) / scale
    basis_cart, cart_rank = _orthonormal_columns(basis_cart, rtol=rtol)
    if cart_rank != info.effective_rank:
        raise RuntimeError(
            "Mass-unweighting changed the numerical rank of the rigid basis "
            f"({info.effective_rank} -> {cart_rank})."
        )

    active_dofs = torch.as_tensor(
        [3 * atom + axis for atom in info.active_atoms for axis in range(3)],
        dtype=torch.long,
        device=coords.device,
    )
    full.index_copy_(0, active_dofs, basis_cart)
    return full, info


def project_hessian_inplace(hessian: torch.Tensor, basis: torch.Tensor) -> torch.Tensor:
    """Apply ``(I - QQ.T) H (I - QQ.T)`` without forming the projector."""

    if basis.shape[1] == 0:
        return hessian
    basis = basis.to(dtype=hessian.dtype, device=hessian.device)
    qt = basis.T
    # Compute both sides from the unmodified input.  H@Q is not generally
    # equal to (Q.T@H).T for a finite-difference Hessian with a small
    # antisymmetric residual.
    qt_h = qt @ hessian
    h_q = hessian @ basis
    qt_h_q = qt_h @ basis
    hessian.addmm_(basis, qt_h, beta=1.0, alpha=-1.0)
    hessian.addmm_(h_q, qt, beta=1.0, alpha=-1.0)
    hessian.addmm_(basis @ qt_h_q, qt, beta=1.0, alpha=1.0)
    return hessian


def compact_project_hessian(hessian, basis: torch.Tensor):
    """Project into the orthogonal complement and return ``(H, lift)``.

    ``lift`` contains orthonormal rows that map reduced eigenvectors back to
    the active mass-weighted Cartesian space.  A zero-rank projection returns
    the original Hessian and ``None`` without allocating an identity matrix.
    """

    rank = int(basis.shape[1])
    if rank == 0:
        return hessian, None
    if isinstance(hessian, torch.Tensor):
        basis = basis.to(dtype=hessian.dtype, device=hessian.device)
        complete, _ = torch.linalg.qr(basis, mode="complete")
        lift = complete[:, rank:].T
        projected = lift @ hessian @ lift.T
        projected = 0.5 * (projected + projected.T)
        return projected, lift

    basis_np = basis.detach().cpu().numpy().astype(np.asarray(hessian).dtype, copy=False)
    complete, _ = np.linalg.qr(basis_np, mode="complete")
    lift = complete[:, rank:].T
    projected = lift.dot(np.asarray(hessian)).dot(lift.T)
    projected = 0.5 * (projected + projected.T)
    return projected, lift
