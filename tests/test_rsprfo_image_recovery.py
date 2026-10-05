"""Bounded physical-model recovery for failed restricted RFO proposals."""
from types import SimpleNamespace

import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
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
def test_near_pole_rsprfo_returns_a_finite_bounded_image_step(values, gradient):
    opt = model(values, gradient)
    step = opt.optimize()
    assert np.isfinite(step).all()
    assert np.linalg.norm(step) <= opt.trust_radius * (1 + 1e-12)
    assert len(opt.predicted_energy_changes) == 1
    assert opt.predicted_energy_changes[0] == pytest.approx(
        opt.quadratic_model(-opt.forces[-1], opt.cur_H, step))


def test_successful_rsprfo_proposal_retains_the_rfo_prediction():
    opt = model([-.2, .3, .7], [.1, .12, .07])
    step = opt.optimize()
    assert np.linalg.norm(step) <= opt.trust_radius * (1 + 1e-12)
    assert opt.predicted_energy_changes == [opt.rfo_model(-opt.forces[-1], opt.cur_H, step)]


def test_exhausted_rsprfo_recovers_using_its_available_physical_model():
    opt = model([-.2, .3, .7], [.1, .12, .07], cycles=2)
    step = opt.optimize()
    assert np.linalg.norm(step) <= opt.trust_radius * (1 + 1e-12)
    assert opt.predicted_energy_changes == [opt.quadratic_model(-opt.forces[-1], opt.cur_H, step)]
