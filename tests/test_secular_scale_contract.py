"""Homogeneous RFO scaling and NumPy/Torch return-value contracts."""
import numpy as np
import pytest
import torch

from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

BACKENDS = ["numpy", "cpu"] + (["cuda"] if torch.cuda.is_available() else [])
BASE_LAMBDA = np.array([-.4, .2, .7], dtype=np.float64)
BASE_GRADIENT = np.array([.11, -.08, .05], dtype=np.float64)


@pytest.fixture
def optimizer(tmp_path):
    geometry = Geometry(["H"], np.array([1., 0., 0.]), coord_type="cart")
    return RFOptimizer(geometry, hessian_init="unit", max_cycles=1, dump=False,
                       line_search=False, gdiis=False, out_dir=tmp_path)


def host(value):
    return value.detach().cpu().numpy() if isinstance(value, torch.Tensor) else np.asarray(value)


def backend_array(values, backend):
    if backend == "numpy":
        return np.asarray(values, dtype=np.float64).copy()
    return torch.tensor(values, dtype=torch.float64, device=backend)


def assert_eigen_equations(result, values, gradient, alpha):
    step, root, nu, vector = result
    step, root, vector = host(step), float(root), host(vector)
    assert step.shape == (3,) and vector.shape == (4,)
    assert np.isfinite(step).all() and np.isfinite(vector).all() and np.isfinite(root)
    assert np.isfinite(float(nu)) and float(nu) != 0. and vector[-1] != 0.
    np.testing.assert_allclose(np.linalg.norm(vector), 1., rtol=1e-12, atol=0.)
    np.testing.assert_allclose(step, vector[:-1]/vector[-1], rtol=1e-8, atol=0.)
    # Componentwise original equations, not a residual divided by the largest
    # matrix entry: every small row and the g.s=mu equation must be accurate.
    diagonal_term, rhs = values*step, alpha*root*step
    residual = diagonal_term+gradient-rhs
    scales = np.abs(diagonal_term)+np.abs(gradient)+np.abs(rhs)
    assert np.all(np.abs(residual) <= 1e-8*scales), (residual, scales)
    residual = float(gradient@step)-root
    scale = float(np.sum(np.abs(gradient*step)))+abs(root)
    assert abs(residual) <= 1e-8*scale, (residual, scale)


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize(
    "scale,alpha,kind,require_secular",
    [
        pytest.param(1e-12, 1., "min", False, id="tiny-scale-normal-alpha"),
        pytest.param(1., 1., "max", True, id="ordinary-benign-must-use-secular"),
        pytest.param(1e12, 1., "min", False, id="large-scale-normal-alpha"),
        pytest.param(1e-12, 1e10, "max", False, id="tiny-scale-large-alpha"),
        pytest.param(1., 1e10, "min", False, id="unit-scale-large-alpha"),
        pytest.param(1e12, 1e10, "max", False, id="large-scale-large-alpha"),
    ],
)
def test_multidimensional_secular_scale_contract(optimizer, backend, scale, alpha, kind, require_secular):
    values, gradient = BASE_LAMBDA*scale, BASE_GRADIENT*scale
    lam, g = backend_array(values, backend), backend_array(gradient, backend)
    # These all-coupled examples have a well-separated extremal augmented
    # root. The unscaled dense reference avoids judging tiny roots against an
    # absolute error floor and explicitly checks homogeneous scaling.
    reference = optimizer.solve_rfo(
        optimizer.get_augmented_hessian(BASE_LAMBDA, BASE_GRADIENT, alpha), kind, alpha=alpha)
    assert_eigen_equations(reference, BASE_LAMBDA, BASE_GRADIENT, alpha)
    dense = optimizer.solve_rfo(optimizer.get_augmented_hessian(lam, g, alpha), kind, alpha=alpha)
    assert_eigen_equations(dense, values, gradient, alpha)
    np.testing.assert_allclose(host(dense[0]), host(reference[0]), rtol=1e-8, atol=0.)
    np.testing.assert_allclose(float(dense[1])/scale, float(reference[1]), rtol=1e-8, atol=0.)
    result = optimizer.solve_rfo_secular(lam, g, alpha, kind=kind)
    if require_secular:
        assert result is not None, "Ordinary coupled case must not unconditionally fall back"
    if result is not None:
        assert_eigen_equations(result, values, gradient, alpha)
        np.testing.assert_allclose(host(result[0]), host(reference[0]), rtol=1e-8, atol=0.)
        np.testing.assert_allclose(float(result[1])/scale, float(reference[1]), rtol=1e-8, atol=0.)
        if backend != "numpy":
            assert result[0].dtype == lam.dtype and result[0].device == lam.device
            assert result[3].dtype == lam.dtype and result[3].device == lam.device
    np.testing.assert_array_equal(host(lam), values)
    np.testing.assert_array_equal(host(g), gradient)
