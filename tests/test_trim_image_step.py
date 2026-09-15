"""Actual one-cycle TRIM steps, without terminal physical certification.

Run with the selected source tree and its tests directory on PYTHONPATH.
The constrained analytic calculator is reused from test_ts_terminal_cadence.
These tests examine the optional TRIM path, not default RS-P-RFO behavior.
"""

import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.tsoptimizers.TRIM import TRIM
from test_ts_terminal_cadence import CountedQuadratic


def make_trim(tmp_path, diagonal, displacement):
    diagonal = np.asarray(diagonal, dtype=float)
    displacement = np.asarray(displacement, dtype=float)
    coords = np.array([
        0., 0., 0.,
        2., 0., 0.,
        0., 2., 0.,
        0., 0., 1.,
    ])
    coords[-3:] += displacement

    geom = Geometry(
        ["H"] * 4,
        coords.copy(),
        coord_type="cart",
        freeze_atoms=[0, 1, 2],
    )
    calculator = CountedQuadratic(tmp_path, diagonal)
    geom.set_calculator(calculator)
    optimizer = TRIM(
        geom,
        hessian_init="calc",
        hessian_update="bofill",
        hessian_recalc=500,
        hessian_recalc_adapt=None,
        hessian_xtb=False,
        roots=[0],
        track_mode_by_overlap=False,
        reject_mode_loss=False,
        trust_radius=.1,
        trust_min=.1,
        trust_max=.1,
        trust_update=False,
        min_line_search=False,
        max_line_search=False,
        verify_saddle=False,
        saddle_recovery_max_cycles=0,
        flatten_enabled=False,
        thresh="never",
        assert_min_step=False,
        energy_plateau=False,
        max_cycles=1,
        out_dir=tmp_path,
        dump=False,
    )
    calculator.optimizer = optimizer
    return geom, optimizer, calculator, coords


def assert_physical_prediction(geom, optimizer, diagonal, initial_coords):
    diagonal = np.asarray(diagonal, dtype=float)
    assert len(optimizer.steps) == len(optimizer.predicted_energy_changes) == 1
    assert not optimizer.is_converged
    assert optimizer.trust_radius == pytest.approx(.1)

    full_step = np.asarray(optimizer.steps[0])
    assert full_step.shape == initial_coords.shape
    assert np.isfinite(full_step).all()
    step = full_step[-3:]
    displacement = initial_coords[-3:] - [0., 0., 1.]
    gradient = diagonal * displacement

    np.testing.assert_array_equal(
        geom.cart_coords[:9], initial_coords[:9]
    )
    np.testing.assert_array_equal(full_step[:9], np.zeros(9))
    np.testing.assert_allclose(
        optimizer.cur_H, np.diag(diagonal), atol=1e-12
    )
    np.testing.assert_allclose(
        geom.cart_coords[-3:], initial_coords[-3:] + step,
        rtol=0., atol=1e-12,
    )

    # Independent physical quadratic oracle, not the optimizer's predictor.
    expected = (
        gradient @ step + .5 * np.dot(diagonal * step, step)
    )
    assert float(optimizer.predicted_energy_changes[0]) == pytest.approx(
        expected, rel=1e-11, abs=1e-13
    )
    return step, gradient


@pytest.mark.parametrize("diagonal, displacement, pd_image", [
    # Root 0 is physical y; complementary negative curvature is physical z.
    ([3., -2., -1.], [0., 0., .01], False),
    ([3., -2., -1.], [0., 0., 0.], False),
    # Nonidentity eigenvector order detects a double basis transformation.
    ([1., -2., 3.], [.01, .02, -.01], True),
], ids=["extra-negative", "hard-case", "pd-image-interior"])
def test_trim_actual_image_subproblem(
    tmp_path, monkeypatch, diagonal, displacement, pd_image
):
    monkeypatch.chdir(tmp_path)
    geom, optimizer, calculator, initial = make_trim(
        tmp_path, diagonal, displacement
    )
    optimizer.run()

    step, gradient = assert_physical_prediction(
        geom, optimizer, diagonal, initial
    )
    assert len(calculator.hessian_calls) == 1
    assert not optimizer._saddle_recovery_active
    assert optimizer.saddle_recovery_steps == 0
    radius = optimizer.trust_radius
    norm = np.linalg.norm(step)
    assert norm <= radius * (1 + 1e-11)

    # Independent image construction in Cartesian active-coordinate order.
    image_diagonal = np.asarray(diagonal, dtype=float).copy()
    image_gradient = gradient.copy()
    target_axis = int(np.argmin(image_diagonal))
    image_diagonal[target_axis] *= -1
    image_gradient[target_axis] *= -1

    if norm < radius * (1 - 1e-9):
        shift = 0.
    else:
        shift = -np.dot(
            step, image_diagonal * step + image_gradient
        ) / np.dot(step, step)

    # Radius alone is insufficient: enforce the independent quadratic KKT
    # conditions, including positive-semidefinite shifted image curvature.
    assert shift >= -1e-10
    assert np.min(image_diagonal + shift) >= -1e-10
    np.testing.assert_allclose(
        (image_diagonal + shift) * step + image_gradient,
        0., atol=1e-10,
    )
    image_change = (
        image_gradient @ step
        + .5 * np.dot(image_diagonal * step, step)
    )
    assert image_change <= 1e-12

    if pd_image:
        np.testing.assert_allclose(
            step, -np.asarray(displacement), rtol=0., atol=1e-12
        )
    else:
        assert norm == pytest.approx(radius, rel=1e-10)
        assert image_change < 0.


def test_trim_recovery_order_and_physical_prediction(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    diagonal = [1., -2., 3.]
    geom, optimizer, calculator, initial = make_trim(
        tmp_path, diagonal, [.01, .02, -.01]
    )
    events = []
    original_terminal = optimizer.validate_terminal_saddle_for_step
    original_recovery = optimizer.apply_saddle_recovery_step
    original_prediction = optimizer.quadratic_model

    def terminal_then_arm(step):
        events.append("terminal")
        original_terminal(step)
        # Inject only the terminal-to-recovery interface state. This does not
        # test physical-Hessian activation of recovery. The real recovery
        # implementation below is exercised without alteration.
        optimizer._saddle_recovery_active = True
        optimizer._saddle_recovery_mode = np.array([0., 1., 0.])
        optimizer._saddle_recovery_sign = 1.
        optimizer.saddle_recovery_step = .03

    def record_recovery(step):
        events.append("recovery")
        return original_recovery(step)

    def record_prediction(gradient, hessian, step):
        events.append("prediction")
        return original_prediction(gradient, hessian, step)

    monkeypatch.setattr(
        optimizer, "validate_terminal_saddle_for_step", terminal_then_arm
    )
    monkeypatch.setattr(
        optimizer, "apply_saddle_recovery_step", record_recovery
    )
    monkeypatch.setattr(optimizer, "quadratic_model", record_prediction)
    optimizer.run()

    step, _ = assert_physical_prediction(
        geom, optimizer, diagonal, initial
    )
    assert len(calculator.hessian_calls) == 1
    assert events == ["terminal", "recovery", "prediction"]
    assert optimizer.saddle_recovery_steps == 1
    np.testing.assert_allclose(
        step, [-.01, .03, .01], rtol=0., atol=1e-12
    )
    assert np.linalg.norm(step) <= optimizer.trust_radius
