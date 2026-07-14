"""Regression tests for RFO trial acceptance and minimization BFGS safety."""

from __future__ import annotations

import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


class _QuadraticCalculator(Calculator):
    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    @staticmethod
    def _results(coords):
        coords = np.asarray(coords, dtype=float)
        return {
            "energy": float(coords @ coords),
            "forces": -2.0 * coords,
        }

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        results = self._results(coords)
        results["hessian"] = 2.0 * np.eye(len(coords))
        return results


def _optimizer(tmp_path):
    geom = Geometry(["H"], np.array([1.0, 0.0, 0.0]), coord_type="cart")
    geom.set_calculator(_QuadraticCalculator(tmp_path))
    opt = RFOptimizer(
        geom,
        hessian_init="unit",
        trust_radius=0.1,
        trust_min=1e-4,
        trust_max=0.1,
        line_search=False,
        gdiis=False,
        out_dir=tmp_path,
    )
    opt.H = np.eye(3)
    return geom, opt


@pytest.mark.parametrize("predicted_change", [-0.5, 0.5])
@pytest.mark.parametrize(
    ("initial_radius", "expected_radius"),
    [(0.1, 0.025), (1e-4, 2.5e-5)],
)
def test_rfo_rejects_uphill_trial_and_rolls_back_all_histories(
    tmp_path,
    predicted_change,
    initial_radius,
    expected_radius,
) -> None:
    geom, opt = _optimizer(tmp_path)
    accepted = geom.coords.copy()
    accepted_cart = geom.cart_coords.copy()
    accepted_energy = geom.energy
    accepted_force = geom.forces.copy()

    # Reproduce Optimizer.run's state at the start of the next cycle: a trial
    # coordinate has already been installed and its step/model prediction are
    # the final entries in their respective histories.
    trial_step = np.array([1.0, 0.0, 0.0])
    geom.coords = accepted + trial_step
    opt.coords = [accepted.copy(), geom.coords.copy()]
    opt.cart_coords = [accepted_cart.copy(), geom.cart_coords.copy()]
    opt.energies = [accepted_energy]
    opt.forces = [accepted_force]
    opt.steps = [trial_step.copy()]
    opt.predicted_energy_changes = [predicted_change]
    opt.image_inds = [[0], [0]]
    opt.image_nums = [1, 1]
    opt.cur_cycle = 1
    opt.trust_radius = initial_radius

    energy, gradient, hessian, _, _, resetted = opt.housekeeping()

    np.testing.assert_allclose(geom.coords, accepted)
    np.testing.assert_allclose(geom.cart_coords, accepted_cart)
    assert energy == accepted_energy
    np.testing.assert_allclose(gradient, -accepted_force)
    np.testing.assert_allclose(hessian, np.eye(3))
    assert resetted is True
    assert opt.trust_radius == expected_radius
    assert opt.rejected_uphill_steps == 1
    assert len(opt.coords) == len(opt.cart_coords) == 1
    assert len(opt.energies) == len(opt.forces) == 1
    assert opt.steps == []
    assert opt.predicted_energy_changes == []
    assert len(opt.image_inds) == len(opt.image_nums) == 1


def test_rfo_skips_bfgs_update_when_curvature_condition_fails(tmp_path) -> None:
    _, opt = _optimizer(tmp_path)
    original_hessian = opt.H.copy()
    opt.coords = [np.zeros(3), np.array([1.0, 0.0, 0.0])]
    opt.steps = [np.array([1.0, 0.0, 0.0])]
    # dg = -(f_new - f_old) = (-1, 0, 0), hence s·y = -1.
    opt.forces = [np.zeros(3), np.array([1.0, 0.0, 0.0])]

    opt.update_hessian()

    np.testing.assert_array_equal(opt.H, original_hessian)
    assert opt.skipped_bfgs_updates == 1


def test_rfo_retains_lower_energy_point_after_rejections_at_floor(tmp_path) -> None:
    geom, opt = _optimizer(tmp_path)
    accepted = geom.coords.copy()
    accepted_cart = geom.cart_coords.copy()
    opt.coords = [accepted.copy()]
    opt.cart_coords = [accepted_cart.copy()]
    opt.energies = [geom.energy]
    opt.forces = [geom.forces.copy()]
    opt.image_inds = [[0]]
    opt.image_nums = [1]
    opt.trust_radius = opt.rejection_trust_floor

    trial_step = np.array([1.0, 0.0, 0.0])
    for cycle in range(1, opt.max_rejections_at_floor + 1):
        geom.coords = accepted + trial_step
        opt.coords.append(geom.coords.copy())
        opt.cart_coords.append(geom.cart_coords.copy())
        opt.steps.append(trial_step.copy())
        opt.predicted_energy_changes.append(-0.5)
        opt.image_inds.append([0])
        opt.image_nums.append(1)
        opt.cur_cycle = cycle
        opt.housekeeping()

    assert opt.uphill_rejection_stalled is True
    assert opt.rejections_at_floor == opt.max_rejections_at_floor
    assert opt.stop_requested is True
    assert "emergency trust floor" in opt.stop_reason
    np.testing.assert_allclose(geom.coords, accepted)
    # The ordinary force/step thresholds fail at this test point.  Retaining
    # the last accepted geometry is safe, but it must not be called converged.
    converged, _ = opt.check_convergence(np.full(3, 0.01))
    assert converged is False
