"""Current-model safeguard through native RFO.optimize, without a PES.

Only housekeeping and the interpolation result are controlled. Native step
selection, composition, prediction, active-space expansion and the new guard
run unchanged. No physical Hessian is requested by this model safeguard.
"""
from types import SimpleNamespace
import importlib

import numpy as np
import pytest
import torch

from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

rf_module = importlib.import_module("pysisyphus.optimizers.RFOptimizer")
BACKENDS = [
    "numpy", "torch_cpu",
    pytest.param("torch_cuda", marks=pytest.mark.skipif(
        not torch.cuda.is_available(), reason="CUDA required for RFO model-guard parity"
    )),
]


def _array(value):
    if isinstance(value, torch.Tensor):
        return value.detach().cpu().numpy()
    return np.asarray(value)


@pytest.fixture(params=BACKENDS)
def make_case(request, tmp_path, monkeypatch):
    def make(offset, *, branch="gdiis", hard_case=False, active_order=(0, 1, 2)):
        def backend(value):
            value = np.asarray(value, dtype=float)
            if request.param == "numpy":
                return value
            device = "cuda" if request.param == "torch_cuda" else "cpu"
            return torch.as_tensor(value, dtype=torch.float64, device=device)

        active = np.asarray(active_order)
        geom = Geometry(["H", "H"], np.zeros(6), coord_type="cart", freeze_atoms=[1])

        def no_pes(*_args, **_kwargs):
            pytest.fail("Model selection must not request a PES evaluation")

        geom.set_calculator(SimpleNamespace(
            get_energy=no_pes, get_forces=no_pes, get_hessian=no_pes,
        ))
        geom.within_partial_hessian = dict(
            active_dofs=active, active_atoms=np.array([0]),
            active_n_dof=3, full_n_dof=6,
        )
        opt = RFOptimizer(
            geom, hessian_init="unit", trust_radius=1e-3, trust_min=1e-4, trust_max=1e-3,
            thresh="baker", dump=False, out_dir=tmp_path, adapt_step_func=True,
            gdiis=True, gediis=False, line_search=True,
        )
        opt._set_active_dofs(True)
        values = np.array([-1., 1., 2.]) if hard_case else np.full(3, .01)
        gradient_np = np.zeros(3) if hard_case else np.array([5e-5, 0., 0.])
        gradient = backend(gradient_np); H = backend(np.diag(values))
        opt.H = opt.cur_H = H
        force = np.zeros(6); force[active] = -gradient_np
        previous_step = np.zeros(6); previous_step[active[0]] = -1e-3
        opt.coords = [-previous_step, np.zeros(6)]
        opt.cart_coords = [v.copy() for v in opt.coords]
        opt.forces = [force.copy(), force.copy()]
        opt.steps = [previous_step]
        opt.energies = [1e-5, 0.]
        opt.cur_cycle = 1
        monkeypatch.setattr(opt, "housekeeping", lambda: (
            0., gradient, H, backend(values), backend(np.eye(3)), False,
        ))
        seen = {}
        native_accept = opt._accept_accelerated_step

        def accept(step, ip_step, reference):
            result = native_accept(step, ip_step, reference)
            seen["reference"] = _array(reference).copy()
            seen["composed"] = _array(result).copy()
            return result

        monkeypatch.setattr(opt, "_accept_accelerated_step", accept)

        def diis(_errors, coords, _forces, _reference, **_kwargs):
            if branch != "gdiis":
                return None
            return SimpleNamespace(coords=coords[-1]+np.asarray(offset),
                                   forces=backend(np.zeros(3)))

        def line_search(*_args, **kwargs):
            assert kwargs["cubic_max_x"] == -1 and kwargs["quartic_max_x"] == 2
            if branch == "poly":
                return -1e-5, backend(np.zeros(3)), backend(offset)
            assert branch == "no_fit"
            return None, None, None

        monkeypatch.setattr(rf_module, "gdiis", diis)
        monkeypatch.setattr(rf_module, "poly_line_search", line_search)
        return SimpleNamespace(opt=opt, active=active, gradient=gradient, H=H,
                               seen=seen, backend=backend)
    return make


@pytest.mark.parametrize("offset", [[4e-4, 0., 0.], [0., 0., 0.]])
@pytest.mark.parametrize("branch", ["gdiis", "poly"])
def test_model_uphill_or_zero_acceleration_keeps_descending_reference(make_case, offset, branch):
    case = make_case(offset, branch=branch)
    step = case.opt.optimize()
    np.testing.assert_allclose(step[case.active], case.seen["reference"], rtol=1e-12, atol=1e-14)
    np.testing.assert_array_equal(step[3:], 0.)
    assert case.opt.quadratic_model(case.gradient, case.H, case.seen["composed"]) >= 0
    assert case.opt.predicted_energy_changes[-1] < 0


def test_descending_acceleration_is_retained_even_when_reference_predicts_more_descent(make_case):
    case = make_case([-4e-4, 0., 0.])
    step = case.opt.optimize()
    np.testing.assert_array_equal(step[case.active], case.seen["composed"])
    qref = case.opt.quadratic_model(case.gradient, case.H, case.seen["reference"])
    assert qref < case.opt.predicted_energy_changes[-1] < 0


def test_non_descending_reference_does_not_add_a_new_rejection_rule(make_case, monkeypatch):
    case = make_case([4e-4, 0., 0.])
    calls = []

    def controlled_reference(_values, _vectors, _gradient):
        calls.append(None)
        return np.array([1e-3, 0., 0.]) if len(calls) == 1 else np.zeros(3)

    monkeypatch.setattr(case.opt, "get_step_func", lambda *_: (
        controlled_reference, case.opt.quadratic_model,
    ))
    step = case.opt.optimize()
    np.testing.assert_array_equal(step[case.active], case.seen["composed"])
    assert case.opt.quadratic_model(case.gradient, case.H, case.seen["reference"]) > 0
    assert case.opt.predicted_energy_changes[-1] > 0


def test_no_fit_keeps_the_native_reference_and_prediction(make_case):
    case = make_case([0., 0., 0.], branch="no_fit")
    step = case.opt.optimize()
    assert "composed" not in case.seen
    np.testing.assert_allclose(step[case.active], [-1e-3, 0., 0.], atol=1e-14)
    assert case.opt.predicted_energy_changes[-1] < 0


def test_zero_gradient_negative_curvature_keeps_descending_acceleration(make_case):
    case = make_case([0., 2e-4, 0.], hard_case=True)
    step = case.opt.optimize()
    np.testing.assert_array_equal(step[case.active], case.seen["composed"])
    assert np.linalg.norm(step) == pytest.approx(1e-3)
    assert case.opt.predicted_energy_changes[-1] < 0
    assert np.dot(_array(case.gradient), step[case.active]) == 0


def test_guard_uses_compact_model_order_before_expanding_frozen_dofs(make_case):
    case = make_case([4e-4, 0., 0.], active_order=(2, 0, 1))
    step = case.opt.optimize()
    np.testing.assert_allclose(step, [0., 0., -1e-3, 0., 0., 0.], atol=1e-14)


@pytest.mark.parametrize(("target", "value"), [
    (target, value) for target in ("reference", "candidate")
    for value in (np.nan, np.inf, -np.inf)
] + [("reference", 0.)])
def test_guard_requires_finite_predictions_and_strictly_descending_reference(
    make_case, monkeypatch, target, value,
):
    case = make_case([4e-4, 0., 0.])
    native_model = case.opt.quadratic_model

    def model(gradient, H, step):
        is_reference = _array(step)[0] < 0
        if is_reference == (target == "reference"):
            return value
        return native_model(gradient, H, step)

    monkeypatch.setattr(case.opt, "quadratic_model", model)
    step = case.opt.optimize()
    np.testing.assert_array_equal(step[case.active], case.seen["composed"])
