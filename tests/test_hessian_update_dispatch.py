"""Honor the configured single-pair Hessian update without an implicit method switch."""
import importlib
from types import SimpleNamespace

import numpy as np
import pytest
import torch

from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers import hessian_updates

owner = importlib.import_module("pysisyphus.optimizers.HessianOptimizer")
BACKENDS = ["numpy", "torch_cpu", pytest.param("torch_cuda", marks=pytest.mark.skipif(
    not torch.cuda.is_available(), reason="CUDA required for Hessian routing parity",
))]


def _numpy(value):
    return value.detach().cpu().numpy() if isinstance(value, torch.Tensor) else np.asarray(value)


def _forbidden(*_args, **_kwargs):
    pytest.fail("This branch must not call the new routing decision or a PES")


@pytest.fixture(params=BACKENDS)
def make_case(request, tmp_path):
    def make(diagonal=(-1., 2.), *, y=(1., 2.), family="bfgs", window=1):
        def backend(value):
            value = np.asarray(value, dtype=float)
            if request.param == "numpy":
                return value.copy()
            device = "cuda" if request.param == "torch_cuda" else "cpu"
            return torch.as_tensor(value.copy(), dtype=torch.float64, device=device)

        # Exercise the native compact 2x2 update after selecting full-vector DOFs.
        active = np.array([0, 2])
        geom = Geometry(["H"], np.zeros(3), coord_type="cart")
        geom.set_calculator(SimpleNamespace(
            get_energy=_forbidden, get_forces=_forbidden, get_hessian=_forbidden,
        ))
        geom.within_partial_hessian = dict(
            active_dofs=active, active_atoms=np.array([0]), active_n_dof=2, full_n_dof=3,
        )
        opt = RFOptimizer(
            geom, hessian_init="unit", hessian_update=family, hessian_update_window=window,
            hessian_recalc=500, trust_radius=.1, trust_min=1e-4, trust_max=.1,
            thresh="baker", dump=False, out_dir=tmp_path,
        )
        opt._set_active_dofs(True)
        H = backend(np.diag(diagonal)); opt.H = H
        step = np.zeros(3); step[active] = 1.
        force = np.zeros(3); force[active] = -np.asarray(y)
        opt.steps = [step.copy()]; opt.forces = [np.zeros(3), force.copy()]
        opt.cur_cycle = 1; opt.hessian_recalc_in = 500
        logs = []; opt.log = logs.append
        return SimpleNamespace(opt=opt, H=H, H_before=_numpy(H).copy(), s=np.ones(2),
                               y=np.asarray(y), step=step, force=force, logs=logs, backend=backend)
    return make


def _unchanged_inputs(case):
    np.testing.assert_array_equal(_numpy(case.H), case.H_before)
    np.testing.assert_array_equal(case.opt.steps[-1], case.step)
    np.testing.assert_array_equal(case.opt.forces[-1], case.force)
    assert case.opt.hessian_recalc_in == 499
    assert case.opt.H.dtype == case.H.dtype
    if isinstance(case.H, torch.Tensor):
        assert case.opt.H.device == case.H.device


@pytest.mark.parametrize("diagonal", [(-1., 1.001), (-1e-12, 2.)])
def test_indefinite_model_retains_selected_bfgs(make_case, monkeypatch, diagonal):
    case = make_case(diagonal)
    expected_delta, _ = hessian_updates.bfgs_update(case.H, case.s, case.y)
    alternative_delta, _ = hessian_updates.ts_bfgs_update(case.H, case.s, case.y)
    assert not np.allclose(_numpy(expected_delta), _numpy(alternative_delta), rtol=1e-8, atol=1e-14)
    monkeypatch.setattr(owner, "get_xp", _forbidden)
    monkeypatch.setattr(owner, "ts_bfgs_update", _forbidden)

    case.opt.update_hessian()

    np.testing.assert_allclose(_numpy(case.opt.H), _numpy(case.H + expected_delta), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(_numpy(case.opt.H) @ case.s, case.y, rtol=1e-12, atol=1e-12)
    assert "Did BFGS Hessian update." in case.logs
    assert case.opt.skipped_bfgs_updates == 0
    _unchanged_inputs(case)


@pytest.mark.parametrize("diagonal", [(1., 2.), (0., 2.)])
def test_nonnegative_model_retains_exact_bfgs(make_case, monkeypatch, diagonal):
    case = make_case(diagonal)
    expected_delta, _ = hessian_updates.bfgs_update(case.H, case.s, case.y)
    expected = _numpy(case.H + expected_delta).copy()
    monkeypatch.setattr(owner, "ts_bfgs_update", _forbidden)

    case.opt.update_hessian()

    np.testing.assert_array_equal(_numpy(case.opt.H), expected)
    np.testing.assert_allclose(_numpy(case.opt.H) @ case.s, case.y, atol=1e-12)
    assert "Did BFGS Hessian update." in case.logs
    _unchanged_inputs(case)


@pytest.mark.parametrize("y", [(-1., 0.), (-1., 1.)])
def test_nonpositive_sTy_still_skips_before_model_routing(make_case, monkeypatch, y):
    case = make_case(y=y)
    monkeypatch.setattr(owner, "get_xp", _forbidden)
    monkeypatch.setattr(owner, "ts_bfgs_update", _forbidden)

    case.opt.update_hessian()

    assert case.opt.H is case.H
    assert case.opt.skipped_bfgs_updates == 1
    assert "Skipped unsafe BFGS Hessian update." in case.logs
    assert case.opt._sy_buffer_S == case.opt._sy_buffer_Y == []
    _unchanged_inputs(case)


@pytest.mark.parametrize("family", [
    "none", "damped_bfgs", "flowchart", "bofill", "ts_bfgs", "ts_bfgs_org", "ts_bfgs_rev",
])
def test_other_families_do_not_enter_new_routing(make_case, monkeypatch, family):
    case = make_case(family=family)
    expected_delta, _ = case.opt.hessian_update_func(case.H, case.s, case.y)
    expected = _numpy(case.H + expected_delta).copy()
    monkeypatch.setattr(owner, "get_xp", _forbidden)

    case.opt.update_hessian()

    np.testing.assert_allclose(_numpy(case.opt.H), expected, rtol=1e-12, atol=1e-12)
    assert case.opt.skipped_bfgs_updates == 0
    _unchanged_inputs(case)


def test_window_two_keeps_existing_multistep_owner(make_case, monkeypatch):
    case = make_case(window=2)
    old_s = np.array([1., 0.]); old_y = np.array([2., 0.])
    case.opt._sy_buffer_S = [old_s.copy()]; case.opt._sy_buffer_Y = [old_y.copy()]
    S = np.column_stack((old_s, case.s)); Y = np.column_stack((old_y, case.y))
    expected_delta, _ = hessian_updates.multistep_ts_bfgs_update(case.H_before, S, Y)
    expected = _numpy(case.H + case.backend(expected_delta)).copy()
    monkeypatch.setattr(owner, "get_xp", _forbidden)
    monkeypatch.setattr(owner, "ts_bfgs_update", _forbidden)

    case.opt.update_hessian()

    np.testing.assert_allclose(_numpy(case.opt.H), expected, rtol=1e-12, atol=1e-12)
    np.testing.assert_array_equal(np.column_stack(case.opt._sy_buffer_S), S)
    np.testing.assert_array_equal(np.column_stack(case.opt._sy_buffer_Y), Y)
    assert "Did MS-TS-BFGS Hessian update (window=2)." in case.logs
    _unchanged_inputs(case)
