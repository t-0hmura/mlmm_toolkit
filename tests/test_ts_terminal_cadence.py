"""Terminal Hessian cadence on a counted constrained quadratic."""

import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer


class CountedQuadratic(Calculator):
    """One movable H relative to three frozen noncollinear H anchors."""

    def __init__(self, out_dir, diagonal):
        super().__init__(out_dir=out_dir, check_mem=False)
        self.diagonal = np.asarray(diagonal, dtype=float)
        self.hessian_calls = []
        self.optimizer = None

    def get_forces(self, atoms, coords, **kwargs):
        displacement = np.asarray(coords)[-3:] - [0., 0., 1.]
        gradient = np.zeros(12)
        gradient[-3:] = self.diagonal * displacement
        return {"energy": float(.5 * np.dot(gradient[-3:], displacement)),
                "forces": -gradient}  # no hidden Hessian in E/g results

    get_energy = get_forces

    def get_hessian(self, atoms, coords, **kwargs):
        self.hessian_calls.append({
            "cycle": int(self.optimizer.cur_cycle),
            "coords": np.asarray(coords).copy(),
        })
        hessian = np.zeros((12, 12))
        hessian[-3:, -3:] = np.diag(self.diagonal)
        return {**self.get_forces(atoms, coords), "hessian": hessian}


def run_quadratic(
    tmp_path, cadence, *, first_order=False, flatten=False, zero_proposal=False,
):
    diagonal = [-.04, .02 if first_order else -.02, .10]
    coords = np.array([0., 0., 0., 2., 0., 0., 0., 2., 0., 0., 0., 1.])
    geom = Geometry(["H"] * 4, coords, coord_type="cart", freeze_atoms=[0, 1, 2])
    calculator = CountedQuadratic(tmp_path, diagonal)
    geom.set_calculator(calculator)
    optimizer = RSPRFOptimizer(
        geom, hessian_init="calc", hessian_update="bofill", hessian_recalc=cadence,
        hessian_recalc_adapt=None, hessian_xtb=False, thresh="baker", max_cycles=6,
        trust_radius=1e-4, trust_min=1e-4, trust_max=1e-4, trust_update=False,
        min_line_search=False, max_line_search=False, verify_saddle=True,
        saddle_imaginary_threshold_cm=5., saddle_recovery_max_cycles=0,
        flatten_enabled=flatten, energy_plateau=False,
        reference_mode=np.r_[np.zeros(9), 1., 0., 0.], out_dir=tmp_path, dump=False,
    )
    calculator.optimizer = optimizer
    images, terminals = {}, []
    original_image = optimizer._image_trust_step
    original_terminal = optimizer._validate_terminal_exact_saddle

    def record_image():
        step, gradient = original_image()
        # The last image call in a cycle owns the ultimately applied proposal.
        images[optimizer.cur_cycle] = (step.copy(), gradient.copy(), optimizer.cur_H.copy())
        return step, gradient

    def record_terminal():
        before = len(calculator.hessian_calls)
        result = original_terminal()
        terminals.append({"cycle": optimizer.cur_cycle,
                          "new_calculator_H": len(calculator.hessian_calls) - before,
                          "n_imaginary": optimizer._last_exact_n_imaginary})
        return result

    optimizer._image_trust_step = record_image
    optimizer._validate_terminal_exact_saddle = record_terminal
    if zero_proposal:
        # Force a zero ordinary proposal after the first exact HOSP. The
        # deferred gate must still select the physical image step.
        original_line_search = optimizer.step_and_grad_from_line_search

        def zero_after_first(energy, gradient, *args):
            if optimizer.cur_cycle > 0:
                return np.zeros_like(gradient), np.zeros_like(gradient)
            return original_line_search(energy, gradient, *args)

        optimizer.step_and_grad_from_line_search = zero_after_first
    optimizer.run()
    np.testing.assert_array_equal(geom.cart_coords[:9], coords[:9])
    for recorded in calculator.hessian_calls:
        np.testing.assert_array_equal(recorded["coords"][:9], coords[:9])
    return geom, optimizer, calculator, images, terminals


def assert_persistent_hosp_dynamics(geom, optimizer, calculator, images):
    assert len(optimizer.steps) == len(optimizer.predicted_energy_changes) == 6
    assert not optimizer.is_converged and not optimizer.stopped
    assert optimizer.stop_reason == ""
    assert optimizer._last_exact_n_imaginary == 2
    assert not optimizer._last_exact_saddle_verified
    assert not optimizer._exact_terminal_candidate_matches_current_geometry()
    assert not optimizer._exact_phva_matches_current_geometry()  # final step moved
    assert optimizer._last_rigid_projection_info["effective_rank"] == 0
    assert sorted(images) == list(range(6)), "cadence screening lost image continuation"
    for cycle, full_step in enumerate(optimizer.steps):
        step = np.asarray(full_step)[-3:]
        assert np.isfinite(full_step).all()
        assert 0. < np.linalg.norm(full_step) <= 1e-4 * (1 + 1e-10)
        np.testing.assert_array_equal(np.asarray(full_step)[:9], np.zeros(9))
        image_step, gradient, physical_h = images[cycle]
        np.testing.assert_allclose(step, image_step, rtol=0., atol=1e-14)
        # Independent E/g/H oracle, not an RFO/PHVA helper used as its own answer.
        displacement = np.asarray(optimizer.coords[cycle])[-3:] - [0., 0., 1.]
        true_gradient = calculator.diagonal * displacement
        np.testing.assert_allclose(gradient, true_gradient, rtol=1e-10, atol=1e-14)
        np.testing.assert_allclose(physical_h, np.diag(calculator.diagonal), atol=1e-12)
        expected = true_gradient @ step + .5 * np.dot(calculator.diagonal * step, step)
        # At tiny trust RFO differs by ~1e-8 relative; tolerance still detects it.
        assert optimizer.predicted_energy_changes[cycle] == pytest.approx(
            expected, rel=1e-10, abs=1e-22,
        )


@pytest.mark.parametrize("cadence, expected_cycles, countdown", [
    (500, [0], 495), (3, [0, 3], 1),
])
@pytest.mark.parametrize("zero_proposal", [False, True])
def test_hosp_model_preserves_terminal_hessian_cadence(
    tmp_path, cadence, expected_cycles, countdown, zero_proposal,
):
    geom, optimizer, calculator, images, terminals = run_quadratic(
        tmp_path, cadence, zero_proposal=zero_proposal,
    )
    assert_persistent_hosp_dynamics(geom, optimizer, calculator, images)
    assert terminals and terminals[0]["cycle"] == 0
    assert all(item["n_imaginary"] == 2 for item in terminals)
    actual_cycles = [item["cycle"] for item in calculator.hessian_calls]
    # Deliberately last: old-source failure retains all preceding safety checks.
    assert (actual_cycles, optimizer.hessian_recalc_in) == (expected_cycles, countdown), (
        "Repeated terminal acquisition/reset after measured HOSP", actual_cycles,
        optimizer.hessian_recalc_in, terminals,
    )


@pytest.mark.parametrize("cadence", [None, np.inf], ids=["none", "infinite"])
def test_no_finite_physical_cadence_retains_exact_rechecks(tmp_path, cadence):
    geom, optimizer, calculator, images, _ = run_quadratic(tmp_path, cadence)
    assert_persistent_hosp_dynamics(geom, optimizer, calculator, images)
    assert [item["cycle"] for item in calculator.hessian_calls] == list(range(6))


@pytest.mark.parametrize("first_order, flatten, count", [
    (True, False, 1), (False, True, 2),
], ids=["exact-fosp", "explicit-flatten-hosp"])
def test_legitimate_exact_terminal_candidate_is_not_deferred(
    tmp_path, first_order, flatten, count,
):
    geom, optimizer, calculator, _, terminals = run_quadratic(
        tmp_path, 500, first_order=first_order, flatten=flatten,
    )
    assert optimizer.is_converged and len(optimizer.steps) == 1
    assert optimizer._last_exact_n_imaginary == count
    assert optimizer._exact_phva_matches_current_geometry()
    assert optimizer._exact_terminal_candidate_matches_current_geometry()
    assert optimizer._exact_saddle_matches_current_geometry() is first_order
    assert [item["cycle"] for item in calculator.hessian_calls] == [0]
    assert len(terminals) == 1
    before = len(calculator.hessian_calls)
    optimizer.validate_terminal_saddle_for_step(np.zeros(3))
    assert len(calculator.hessian_calls) == before  # unchanged exact coordinates
