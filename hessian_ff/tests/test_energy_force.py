"""Tests: energy and force correctness of hessian_ff vs OpenMM."""

from __future__ import annotations

import numpy as np
import torch

from hessian_ff.terms.angle import AngleTerm
from hessian_ff.terms.dihedral import DihedralTerm
from hessian_ff.terms.nonbonded import NonbondedTerm

E_ABS_TOL = 1.0e-3  # kcal/mol
F_RMS_TOL = 1.0e-3  # kcal/mol/A


def _openmm_energy_force(context):
    """Extract energy (kcal/mol) and forces (kcal/mol/A) from OpenMM context."""
    from openmm import unit

    state = context.getState(getEnergy=True, getForces=True)
    e = float(state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole))
    f = state.getForces(asNumpy=True).value_in_unit(
        unit.kilocalories_per_mole / unit.angstrom,
    )
    return e, torch.as_tensor(np.asarray(f), dtype=torch.float64)


class TestEnergyForceVsOpenMM:
    """Verify hessian_ff energy/force match OpenMM on the small system."""

    def test_energy_matches_openmm(
        self, small_ff, small_coords, openmm_context_and_pos,
    ):
        context, _ = openmm_context_and_pos
        e_hff = float(small_ff(small_coords)["E_total"].detach().cpu())
        e_omm, _ = _openmm_energy_force(context)
        assert abs(e_hff - e_omm) <= E_ABS_TOL, (
            f"Energy mismatch: hff={e_hff:.6f}, omm={e_omm:.6f}, "
            f"diff={abs(e_hff - e_omm):.6e}"
        )

    def test_force_matches_openmm(
        self, small_ff, small_coords, openmm_context_and_pos,
    ):
        context, _ = openmm_context_and_pos
        _, f_hff = small_ff.energy_force(small_coords, force_calc_mode="Analytical")
        f_hff = f_hff.detach().cpu().to(torch.float64)
        _, f_omm = _openmm_energy_force(context)
        diff = (f_hff - f_omm).numpy()
        rms = float(np.sqrt(np.mean(diff ** 2)))
        assert rms <= F_RMS_TOL, (
            f"Force RMS mismatch: {rms:.6e} > tol={F_RMS_TOL}"
        )

    def test_python_bonded_force_path_matches_openmm(
        self, small_ff, small_coords, openmm_context_and_pos,
    ):
        context, _ = openmm_context_and_pos
        coords = small_coords.clone().detach().requires_grad_(True)
        _, f_hff = small_ff.energy_force(coords, force_calc_mode="Analytical")
        _, f_omm = _openmm_energy_force(context)
        diff = f_hff.detach().cpu().to(torch.float64) - f_omm
        rms = float(torch.sqrt(torch.mean(diff * diff)))
        assert rms <= F_RMS_TOL


def test_exact_linear_angle_energy_force_is_finite_and_consistent() -> None:
    term = AngleTerm(
        i=torch.tensor([0]),
        j=torch.tensor([1]),
        k=torch.tensor([2]),
        k_theta=torch.tensor([50.0], dtype=torch.float64),
        theta0=torch.tensor([2.0], dtype=torch.float64),
    )
    coords = torch.tensor(
        [[-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        dtype=torch.float64,
        requires_grad=True,
    )

    energy, force = term.energy_force(coords)
    autograd_gradient, = torch.autograd.grad(term(coords), coords)

    assert torch.isfinite(energy)
    assert torch.isfinite(force).all()
    assert torch.isfinite(autograd_gradient).all()
    torch.testing.assert_close(energy, term(coords))
    hessian = torch.autograd.functional.hessian(term, coords)
    assert torch.isfinite(hessian).all()


def test_planar_dihedral_autograd_hessian_is_finite() -> None:
    term = DihedralTerm(
        i=torch.tensor([0]),
        j=torch.tensor([1]),
        k=torch.tensor([2]),
        l=torch.tensor([3]),
        force=torch.tensor([1.0], dtype=torch.float64),
        period=torch.tensor([2.0], dtype=torch.float64),
        phase=torch.tensor([0.0], dtype=torch.float64),
    )
    coords = torch.tensor(
        [
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
        ],
        dtype=torch.float64,
    )

    hessian = torch.autograd.functional.hessian(term, coords)

    assert torch.isfinite(hessian).all()


def test_nonbonded_autograd_preserves_force_and_hessian_graph() -> None:
    empty_i = torch.empty(0, dtype=torch.int64)
    empty_f = torch.empty(0, dtype=torch.float64)
    term = NonbondedTerm(
        natom=2,
        charge=torch.tensor([1.0, -1.0], dtype=torch.float64),
        atom_type=torch.tensor([0, 0], dtype=torch.int64),
        lj_acoef=torch.tensor([1.0], dtype=torch.float64),
        lj_bcoef=torch.tensor([0.5], dtype=torch.float64),
        hb_acoef=empty_f,
        hb_bcoef=empty_f,
        nb_index=torch.tensor([[1]], dtype=torch.int64),
        pair_i=torch.tensor([0], dtype=torch.int64),
        pair_j=torch.tensor([1], dtype=torch.int64),
        pair14_i=empty_i,
        pair14_j=empty_i,
        pair14_inv_scee=empty_f,
        pair14_inv_scnb=empty_f,
    )
    coords = torch.tensor(
        [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
        dtype=torch.float64,
        requires_grad=True,
    )

    energy = term(coords).coulomb + term(coords).lj
    gradient, = torch.autograd.grad(energy, coords, create_graph=True)
    hessian = torch.autograd.functional.hessian(
        lambda value: term(value).coulomb + term(value).lj,
        coords,
    )

    assert energy.grad_fn is not None
    assert torch.isfinite(gradient).all()
    assert torch.isfinite(hessian).all()
