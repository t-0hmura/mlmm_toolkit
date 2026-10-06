"""RS-P-RFO restricted steps at weakly coupled negative roots."""
from types import SimpleNamespace

import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from test_rfo_symmetric import rs_optimizer


def model(values, gradient, cycles=25):
    values, gradient = np.asarray(values, dtype=np.float32), np.asarray(gradient, dtype=np.float32)
    opt = rs_optimizer(cycles)
    opt.geometry = Geometry(['H'], np.zeros(3), coord_type='cart')
    opt.H = opt.cur_H = np.diag(values)
    opt.forces = [-gradient.copy()]
    opt.small_eigval_thresh = 1e-8
    opt.table = SimpleNamespace(print=lambda *_: None)
    opt.housekeeping = lambda: (0., gradient, opt.cur_H, values, np.eye(3, dtype=np.float32), False)
    return opt


@pytest.mark.parametrize('values,gradient', [
    ([-.1, -.00384848, .4], [.001, 1e-7, .01]),
    ([-.1, -.00384848, .4], [.001, 1e-10, .01]),
])
def test_near_pole_rsprfo_returns_a_finite_bounded_rfo_step(values, gradient):
    opt = model(values, gradient)
    def forbidden():
        raise AssertionError('A coupled finite RFO proposal used image recovery.')
    opt._image_trust_step = forbidden
    step = opt.optimize()
    assert np.isfinite(step).all()
    assert np.linalg.norm(step) <= opt.trust_radius * (1 + 1e-12)
    physical = -opt.forces[-1] @ step + .5 * step @ opt.cur_H @ step
    expected = physical / (1 + step @ step)
    assert len(opt.predicted_energy_changes) == 1
    assert opt.predicted_energy_changes[0] == pytest.approx(expected)


def test_successful_rsprfo_proposal_retains_the_rfo_prediction():
    opt = model([-.2, .3, .7], [.1, .12, .07])
    step = opt.optimize()
    assert np.linalg.norm(step) <= opt.trust_radius * (1 + 1e-12)
    physical = -opt.forces[-1] @ step + .5 * step @ opt.cur_H @ step
    assert len(opt.predicted_energy_changes) == 1
    assert opt.predicted_energy_changes[0] == pytest.approx(physical / (1 + step @ step))


def test_uncoupled_adverse_root_retains_existing_bounded_recovery():
    opt = model([-.1, -.00384848, .4], [.001, 0., .01])
    step = opt.optimize()
    assert np.isfinite(step).all()
    assert np.linalg.norm(step) <= opt.trust_radius * (1 + 1e-12)
    assert len(opt.predicted_energy_changes) == 1
    assert opt.predicted_energy_changes[0] == pytest.approx(
        -opt.forces[-1] @ step + .5 * step @ opt.cur_H @ step)


# Minimizing partition of a recorded TS cycle (6-atom test system, alpha = 1):
# a second negative root (-0.00385) couples to the gradient only at 7e-10.
WEAK_ROOT_EIGVALS = np.array([
    -0.0038501204962899313, -1.5722653248683023e-07, 2.6950288638850595e-08,
    2.2529665962161209e-07, 0.0018986704102999048, 0.0019009172892921133,
    0.010790207661807936, 0.010793101472913368, 0.04246346240604662,
    0.04247310251041955, 0.05783863083073537, 0.08365661003569581,
    0.10371178260570468, 0.10372192237815704, 0.4337802182945495,
    1.0590315417260516, 1.059173735684847,
])
WEAK_ROOT_GRADIENT = np.array([
    6.789405691452958e-10, 7.145981972253986e-11, -1.7696133942333863e-07,
    -5.213124320546961e-08, -4.0594148675673313e-07, -1.9699333957658966e-07,
    -2.6632272423195018e-06, 4.154277920053116e-06, 1.0402383425909612e-05,
    -1.2402168526087384e-06, -0.007544491999937432, -0.0015886418180242107,
    -4.355036006592267e-06, 1.9875919170432025e-06, 5.385294792141848e-05,
    -4.792773977281868e-07, -2.468935840178442e-07,
])
# Eq. 18 at the lowest root of the secular equation, evaluated with 60 digits.
WEAK_ROOT_DSTEP2_DALPHA = -4.84201848253e13


def test_eq18_derivative_is_accurate_at_a_weakly_coupled_root():
    opt = rs_optimizer()
    augmented = opt.get_augmented_hessian(WEAK_ROOT_EIGVALS, WEAK_ROOT_GRADIENT, 1.0)
    step, eigval, _, _ = opt.solve_rfo(augmented, "min", alpha=1.0)
    derivative = RSPRFOptimizer._partition_dstep2_dalpha(
        1.0, eigval, step, WEAK_ROOT_EIGVALS, WEAK_ROOT_GRADIENT
    )
    assert derivative == pytest.approx(WEAK_ROOT_DSTEP2_DALPHA, rel=1e-6)


def test_eq18_ignores_roundoff_in_uncoupled_components():
    eigvals = np.array([-.1, .3, .5])
    gradient = np.array([.01, .02, 0.])
    step = np.array([.05, -.04, 1e-20])
    coupled = RSPRFOptimizer._partition_dstep2_dalpha(
        2., -.01, step[:2], eigvals[:2], gradient[:2]
    )
    derivative = RSPRFOptimizer._partition_dstep2_dalpha(2., -.01, step, eigvals, gradient)
    assert np.isfinite(derivative)
    assert derivative == pytest.approx(coupled)


def test_positive_alpha_update_below_one_ulp_still_advances():
    opt = rs_optimizer()
    alpha0 = opt.alpha0
    def noisy_boundary(values, gradient, alpha, kind, **kwargs):
        # The step norm is outside the trust radius only at the first alpha.
        excess = 1e-9 if alpha == alpha0 else -1e-9
        length = opt.trust_radius * (1 + excess) if kind == "max" else 0.0
        return np.array([length]), -.01, 1.0, np.array([0.0, 1.0])
    opt.solve_rfo_secular = noisy_boundary
    # A steep derivative makes the Newton update smaller than half an ulp.
    opt._partition_dstep2_dalpha = lambda *_: -5e5
    step = opt.optimize()
    assert np.linalg.norm(step) <= opt.trust_radius
    assert len(opt.predicted_energy_changes) == 1
