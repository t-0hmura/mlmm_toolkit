"""Raw RS-RFO feasibility before any caller clamp; fixed diagonal witnesses only."""
import numpy as np
import pytest
import torch
from scipy.optimize import brentq

from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer


LAM = np.array([-3e-5, .01, 1.])
GRAD = np.array([6e-8, 6e-6, 1e-7])


@pytest.fixture(params=['numpy', 'torch'])
def backend(request):
    return request.param


def optimizer(tmp_path, radius, *, micro=50, cls=RFOptimizer):
    geom = Geometry(['H'], np.array([1., 0., 0.]), coord_type='cart')
    options = dict(gdiis=False) if cls is RFOptimizer else dict(roots=[0], verify_saddle=False)
    return cls(geom, hessian_init='unit', dump=False, max_cycles=1,
               trust_radius=radius, trust_min=1e-8, trust_max=1.,
               alpha0=1., max_micro_cycles=micro, rfo_overlaps=False,
               line_search=False, out_dir=tmp_path, **options)


def arguments(lam, gradient, backend):
    values = [np.asarray(lam, dtype=float), np.eye(3), np.asarray(gradient, dtype=float)]
    return [torch.tensor(v, dtype=torch.float64) for v in values] if backend == 'torch' else values


def rfo_value(lam, gradient, step):
    return (gradient @ step + .5 * step @ (lam * step)) / (1. + step @ step)


def boundary_reference(lam, gradient, radius):
    # Independent scalar boundary solve. This reference is certified for these
    # fixtures, not asserted to describe every bounded RFO problem.
    lower = np.nextafter(max(0., -float(np.min(lam))), np.inf)
    shift = brentq(lambda t: np.linalg.norm(gradient / (lam + t)) - radius,
                   lower, lower + 2. * np.linalg.norm(gradient) / radius,
                   xtol=np.finfo(float).smallest_subnormal, rtol=4. * np.finfo(float).eps)
    step = -gradient / (lam + shift)
    value = rfo_value(lam, gradient, step)
    # L = numerator - value*(1+s.s) + mu/2*(s.s-radius**2).
    # PSD Hessian and mu>=0 certify the constrained rational objective.
    mu = shift + 2. * value
    assert mu >= 0. and np.min(lam + shift) >= 0.
    np.testing.assert_allclose((lam + shift) * step + gradient, 0., atol=1e-18)
    assert np.linalg.norm(step) == pytest.approx(radius, rel=1e-12, abs=0.)
    return step, value


def assert_small_boundary(opt, backend):
    reference, value = boundary_reference(LAM, GRAD, opt.trust_radius)
    step = opt.get_rs_step(*arguments(LAM, GRAD, backend))
    assert np.isfinite(step).all()
    assert np.linalg.norm(step) <= opt.trust_radius * (1. + 1e-12)
    np.testing.assert_allclose(step, reference, rtol=1e-7, atol=1e-12)
    assert rfo_value(LAM, GRAD, step) == pytest.approx(value, rel=1e-7, abs=1e-17)


@pytest.mark.parametrize('radius', [1e-4, 2e-4])
def test_small_radius_raw_step_and_certified_objective(tmp_path, backend, radius):
    assert_small_boundary(optimizer(tmp_path, radius), backend)


def test_legitimate_interior_step_is_not_expanded(tmp_path, backend):
    lam, g = np.array([1., 2., 3.]), np.array([.01, 0., 0.])
    step = optimizer(tmp_path, .1).get_rs_step(*arguments(lam, g, backend))
    expected = np.array([-2. * g[0] / (1. + np.sqrt(1. + 4. * g[0]**2)), 0., 0.])
    np.testing.assert_allclose(step, expected, rtol=1e-10, atol=1e-13)
    assert np.linalg.norm(step) < .1


def test_ordinary_point_one_boundary(tmp_path, backend):
    lam, g = np.array([1., 2., 3.]), np.array([1., 2., 3.])
    step = optimizer(tmp_path, .1).get_rs_step(*arguments(lam, g, backend))
    assert np.isfinite(step).all() and np.linalg.norm(step) <= .1 * (1. + 1e-12)
    assert np.linalg.norm(step) == pytest.approx(.1, rel=1e-7, abs=0.)
    assert rfo_value(lam, g, step) < 0.


def test_micro_cap_preserves_native_fallback(tmp_path, backend):
    assert_small_boundary(optimizer(tmp_path, 1e-4, micro=1), backend)


def test_zero_gradient_negative_hessian_keeps_hard_case(tmp_path, backend):
    lam, g, radius = np.array([-.01, .02, 1.]), np.zeros(3), 1e-4
    step = optimizer(tmp_path, radius).get_rs_step(*arguments(lam, g, backend))
    assert np.isfinite(step).all()
    assert np.linalg.norm(step) <= radius * (1. + 1e-12)
    assert abs(step[0]) == pytest.approx(radius, rel=1e-12, abs=0.)
    np.testing.assert_allclose(step[1:], 0., atol=1e-15)
    assert rfo_value(lam, g, step) == pytest.approx(-.005 * radius**2 / (1. + radius**2), rel=1e-12)


def test_rs_irfo_actual_image_path_is_raw_feasible(tmp_path, monkeypatch, backend):
    # Reflect the actual negative root0; the complementary negative direction
    # makes the image problem the small-radius witness above, up to permutation.
    opt = optimizer(tmp_path, 1e-4, cls=RSIRFOptimizer)
    physical_lam, physical_g = np.array([-.01, -3e-5, 1.]), np.array([-6e-6, 6e-8, 1e-7])
    lam, vectors, g = arguments(physical_lam, physical_g, backend)
    H = torch.diag(lam) if backend == 'torch' else np.diag(lam)
    opt.H = opt.cur_H = H
    monkeypatch.setattr(opt, 'housekeeping', lambda: (0., g, H, lam, vectors, False))
    step = opt.optimize()  # Real root selection/projection/eigh/RS step/prediction.
    assert np.isfinite(step).all() and np.linalg.norm(step) <= opt.trust_radius * (1. + 1e-12)
    assert len(opt.predicted_energy_changes) == 1
    assert float(opt.predicted_energy_changes[0]) == pytest.approx(rfo_value(physical_lam, physical_g, step))
    assert not opt._saddle_recovery_active and opt.geometry.calculator is None
