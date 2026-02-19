"""Tests: Hessian correctness of hessian_ff (analytical vs autograd vs FD vs OpenMM FD)."""

from __future__ import annotations

import numpy as np
import pytest
import torch

from hessian_ff.analytical_hessian import build_analytical_hessian

FD_DELTA = 1.0e-4  # Angstrom
H_RELFRO_TOL = 1.0e-3  # relative Frobenius norm tolerance


def _symmetrize(h: torch.Tensor) -> torch.Tensor:
    h = h.detach().cpu().to(torch.float64)
    return 0.5 * (h + h.T)


def _relative_frobenius(a: torch.Tensor, b: torch.Tensor) -> float:
    d = b - a
    return float(torch.linalg.norm(d) / (torch.linalg.norm(a) + 1e-20))


def _analytical_hessian(system, coords):
    """Compute full analytical Hessian via build_analytical_hessian."""
    all_atoms = list(range(system.natom))
    h, _ = build_analytical_hessian(
        system=system, coords=coords, active_atoms=all_atoms,
        mpi_rank=0, mpi_size=1,
    )
    return _symmetrize(h)


def _autograd_hessian(ff, coords):
    """Compute full Hessian via torch.autograd.functional.hessian."""
    x = coords.clone().detach().requires_grad_(True)
    h4 = torch.autograd.functional.hessian(
        lambda x_: ff(x_)["E_total"], x, vectorize=False,
    )
    return _symmetrize(h4.reshape(x.numel(), x.numel()))


def _fd_hessian(ff, coords, delta=FD_DELTA, force_calc_mode="Analytical"):
    """Compute full Hessian via finite-difference on forces."""
    ndof = 3 * coords.shape[0]
    h = torch.empty((ndof, ndof), dtype=torch.float64, device="cpu")
    for col in range(ndof):
        atom, comp = col // 3, col % 3
        xp, xm = coords.clone(), coords.clone()
        xp[atom, comp] += delta
        xm[atom, comp] -= delta
        _, fp = ff.energy_force(xp, force_calc_mode=force_calc_mode)
        _, fm = ff.energy_force(xm, force_calc_mode=force_calc_mode)
        fp_flat = fp.detach().cpu().to(torch.float64).reshape(-1)
        fm_flat = fm.detach().cpu().to(torch.float64).reshape(-1)
        h[:, col] = -(fp_flat - fm_flat) / (2.0 * delta)
    return _symmetrize(h)


def _openmm_fd_hessian(context, base_pos, delta=FD_DELTA):
    """Compute full Hessian via finite-difference on OpenMM forces."""
    from openmm import unit

    natom = base_pos.shape[0]
    ndof = 3 * natom
    h = np.empty((ndof, ndof), dtype=np.float64)
    force_unit = unit.kilocalories_per_mole / unit.angstrom

    for col in range(ndof):
        atom, comp = col // 3, col % 3
        xp, xm = base_pos.copy(), base_pos.copy()
        xp[atom, comp] += delta
        xm[atom, comp] -= delta

        context.setPositions(xp * unit.angstrom)
        fp = np.asarray(
            context.getState(getForces=True).getForces(asNumpy=True)
            .value_in_unit(force_unit), dtype=np.float64,
        ).reshape(-1)

        context.setPositions(xm * unit.angstrom)
        fm = np.asarray(
            context.getState(getForces=True).getForces(asNumpy=True)
            .value_in_unit(force_unit), dtype=np.float64,
        ).reshape(-1)

        h[:, col] = -(fp - fm) / (2.0 * delta)

    # Restore original positions
    context.setPositions(base_pos * unit.angstrom)
    return _symmetrize(torch.from_numpy(h))


class TestHessianConsistency:
    """Verify Hessian methods produce consistent results on the small system."""

    def test_analytical_vs_autograd(
        self, small_system, small_ff_autograd, small_coords,
    ):
        h_ana = _analytical_hessian(small_system, small_coords)
        h_auto = _autograd_hessian(small_ff_autograd, small_coords)
        rel = _relative_frobenius(h_ana, h_auto)
        assert rel <= H_RELFRO_TOL, (
            f"analytical vs autograd: relative Frobenius = {rel:.6e} > {H_RELFRO_TOL}"
        )

    def test_analytical_vs_fd(
        self, small_system, small_ff, small_coords,
    ):
        h_ana = _analytical_hessian(small_system, small_coords)
        h_fd = _fd_hessian(small_ff, small_coords)
        rel = _relative_frobenius(h_ana, h_fd)
        assert rel <= H_RELFRO_TOL, (
            f"analytical vs FD: relative Frobenius = {rel:.6e} > {H_RELFRO_TOL}"
        )

    def test_analytical_vs_openmm_fd(
        self, small_system, small_coords, openmm_context_and_pos,
    ):
        context, base_pos = openmm_context_and_pos
        h_ana = _analytical_hessian(small_system, small_coords)
        h_omm = _openmm_fd_hessian(context, base_pos)
        rel = _relative_frobenius(h_ana, h_omm)
        assert rel <= H_RELFRO_TOL, (
            f"analytical vs OpenMM FD: relative Frobenius = {rel:.6e} > {H_RELFRO_TOL}"
        )

    def test_analytical_hessian_shape(self, small_system, small_coords):
        ndof = 3 * small_system.natom
        h_ana = _analytical_hessian(small_system, small_coords)
        assert h_ana.shape == (ndof, ndof)

    def test_analytical_hessian_symmetric(self, small_system, small_coords):
        h_ana = _analytical_hessian(small_system, small_coords)
        diff = float(torch.max(torch.abs(h_ana - h_ana.T)))
        assert diff < 1e-10, f"Hessian not symmetric: max asymmetry = {diff:.6e}"
