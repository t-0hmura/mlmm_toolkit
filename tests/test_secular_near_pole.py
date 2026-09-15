"""Secular RFO roots must give accurate steps or request dense solving."""
from decimal import Decimal, localcontext

import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


@pytest.fixture
def optimizer(tmp_path):
    geometry = Geometry(["H"], np.array([1.0, 0.0, 0.0]), coord_type="cart")
    return RFOptimizer(
        geometry, hessian_init="unit", max_cycles=1, dump=False,
        line_search=False, gdiis=False, out_dir=tmp_path,
    )  # Real constructor initializes all fields; no prepare_opt/run call.


def analytic_1d(lam, gradient, alpha, kind):
    """90-digit reference for [[lam/alpha,g/sqrt(alpha)],[g/sqrt(alpha),0]]."""
    with localcontext() as context:
        context.prec = 90
        h, g, a = (Decimal.from_float(float(x)) for x in (lam, gradient, alpha))
        diagonal = h/a
        discriminant = (diagonal*diagonal+4*g*g/a).sqrt()
        root = (diagonal+(discriminant if kind == "max" else -discriminant))/2
        # The second original eigen-equation gives g*s=mu. No near-pole
        # denominator subtraction is needed to obtain the reference step.
        return float(root/g), float(root)


def assert_accepted_1d(result, expected_step, expected_root, lam, gradient, alpha):
    step, root, nu, vector = result
    step, root = float(np.asarray(step)[0]), float(root)
    vector = np.asarray(vector)
    assert np.isfinite([step, root, float(nu)]).all() and np.isfinite(vector).all()
    assert float(nu) != 0.0 and vector.shape == (2,) and vector[-1] != 0.0
    np.testing.assert_allclose(step, expected_step, rtol=1e-8, atol=0.0)
    np.testing.assert_allclose(root, expected_root, rtol=1e-8, atol=0.0)
    np.testing.assert_allclose(step, vector[0]/vector[1], rtol=1e-8, atol=0.0)
    np.testing.assert_allclose(np.linalg.norm(vector), 1.0, rtol=1e-12, atol=0.0)
    # Scale EACH original eigen-equation separately. An absolute 1e-12
    # residual would miss a large relative step error when the root is near 1e-12.
    left, right = lam*step+gradient, alpha*root*step
    scale = abs(lam*step)+abs(gradient)+abs(right)
    assert abs(left-right) <= 1e-8*scale, (left, right, scale)
    left = gradient*step
    scale = abs(left)+abs(root)
    assert abs(left-root) <= 1e-8*scale, (left, root, scale)


@pytest.mark.parametrize("kind", ["max", "min"])
@pytest.mark.parametrize(
    "lam,gradient,alpha",
    [
        pytest.param(.057354192967919285, 5.255878300403e-9, 29425517334.491104,
                     id="near-pole-fallback-allowed"),
        pytest.param(.057354192967919285, 5.255878300403e-9, 50927225059.16852,
                     id="near-pole-wrong-small-step-forbidden"),
        pytest.param(.2, .12, 1.0, id="benign-coupled"),
        pytest.param(.2, .12, 1000.0, id="benign-moderate-alpha"),
        pytest.param(0.0, .02, 4.0, id="coupled-zero-curvature"),
    ],
)
def test_secular_result_matches_stable_reference_or_declines(optimizer, kind, lam, gradient, alpha):
    lam = lam if kind == "max" else -lam  # Mirror max pole into min pole.
    values, forces = np.array([lam]), np.array([gradient])
    expected_step, expected_root = analytic_1d(lam, gradient, alpha, kind)
    matrix = optimizer.get_augmented_hessian(values, forces, alpha)
    before = matrix.copy()
    dense = optimizer.solve_rfo(matrix, kind, alpha=alpha)
    assert_accepted_1d(dense, expected_step, expected_root, lam, gradient, alpha)
    np.testing.assert_array_equal(matrix, before)
    secular = optimizer.solve_rfo_secular(values, forces, alpha, kind=kind)
    if secular is not None:
        assert_accepted_1d(secular, expected_step, expected_root, lam, gradient, alpha)
        np.testing.assert_allclose(secular[0], dense[0], rtol=1e-8, atol=0.0)
    np.testing.assert_array_equal(values, [lam])
    np.testing.assert_array_equal(forces, [gradient])


@pytest.mark.parametrize("kind,lam", [("min", 1.0), ("max", -1.0), ("min", 0.0), ("max", 0.0)])
def test_decoupled_finite_zero_step_remains_valid(optimizer, kind, lam):
    values, gradient, alpha = np.array([lam]), np.array([0.0]), 4.0
    dense = optimizer.solve_rfo(optimizer.get_augmented_hessian(values, gradient, alpha), kind, alpha=alpha)
    assert_accepted_1d(dense, 0.0, 0.0, lam, 0.0, alpha)
    secular = optimizer.solve_rfo_secular(values, gradient, alpha, kind=kind)
    if secular is not None:
        assert_accepted_1d(secular, 0.0, 0.0, lam, 0.0, alpha)


@pytest.mark.parametrize("kind,sign", [("max", 1.0), ("min", -1.0)])
@pytest.mark.parametrize("partial_coupling", [False, True])
def test_uncoupled_extremal_root_does_not_fabricate_finite_step(optimizer, kind, sign, partial_coupling):
    values = sign*np.array([1.0, .2] if partial_coupling else [1.0])
    gradient = np.array([0.0, .1] if partial_coupling else [0.0])
    # The true extremal root belongs to the uncoupled first coordinate;
    # its augmented component is zero, even if another coordinate is coupled.
    assert optimizer.solve_rfo_secular(values, gradient, 4.0, kind=kind) is None
    with pytest.raises(ZeroDivisionError, match="zero augmented component"):
        optimizer.solve_rfo(optimizer.get_augmented_hessian(values, gradient, 4.0), kind, alpha=4.0)
