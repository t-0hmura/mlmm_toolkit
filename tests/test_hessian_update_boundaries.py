"""Finite secant boundaries and coordinate-basis history regressions."""

from __future__ import annotations

import importlib
from types import SimpleNamespace

import numpy as np
import pytest
import torch

from pysisyphus.optimizers import hessian_updates
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


TS_UPDATES = (
    "ts_bfgs_update",
    "ts_bfgs_update_org",
    "ts_bfgs_update_revised",
)


@pytest.fixture(params=[
    ("numpy", np.float32), ("numpy", np.float64),
    ("cpu", np.float32), ("cpu", np.float64),
    ("cuda", np.float32), ("cuda", np.float64),
])
def array_backend(request):
    device, dtype = request.param
    if device == "cuda" and not torch.cuda.is_available():
        pytest.skip("CUDA is unavailable")

    def convert(value):
        value = np.asarray(value, dtype=dtype)
        if device == "numpy":
            return value
        return torch.as_tensor(value, device=device)

    return convert, 5e-5 if dtype == np.float32 else 1e-11


def _numpy(value):
    if isinstance(value, torch.Tensor):
        return value.detach().cpu().numpy()
    return np.asarray(value)


def _checked_update(name, hessian, step, gradient_delta):
    before = _numpy(hessian).copy()
    update, label = getattr(hessian_updates, name)(hessian, step, gradient_delta)
    assert update.dtype == hessian.dtype
    if isinstance(hessian, torch.Tensor):
        assert update.device == hessian.device
    np.testing.assert_array_equal(_numpy(hessian), before)
    assert np.isfinite(_numpy(update)).all()
    return _numpy(update), label


@pytest.mark.parametrize("name", TS_UPDATES)
@pytest.mark.parametrize("case", ["exact", "zero_step", "exact_null"])
def test_ts_update_without_new_secant_information(array_backend, name, case):
    convert, _ = array_backend
    hessian = convert(np.diag([0.0 if case == "exact_null" else -1.0, 2.0]))
    step = convert([0.25, 0.0 if case == "exact_null" else 0.5])
    gradient_delta = hessian @ step
    if case == "zero_step":
        step = convert([0.0, 0.0])
        gradient_delta = convert([1.0, 2.0])

    update, _ = _checked_update(name, hessian, step, gradient_delta)

    np.testing.assert_array_equal(update, np.zeros((2, 2)))


@pytest.mark.parametrize("name", TS_UPDATES)
def test_nonzero_null_direction_is_not_silently_skipped(array_backend, name):
    convert, _ = array_backend
    hessian = convert(np.diag([0.0, 2.0]))
    with pytest.raises(ValueError, match="singular.*nonzero secant"):
        getattr(hessian_updates, name)(
            hessian, convert([1.0, 0.0]), convert([0.0, 1.0])
        )
    np.testing.assert_array_equal(_numpy(hessian), np.diag([0.0, 2.0]))


@pytest.mark.parametrize("name", TS_UPDATES)
def test_singular_hessian_with_informative_direction_is_updated(array_backend, name):
    convert, tol = array_backend
    hessian = convert(np.diag([0.0, 2.0]))
    step = convert([0.0, 0.5])
    gradient_delta = convert([0.25, 1.5])

    update, _ = _checked_update(name, hessian, step, gradient_delta)

    np.testing.assert_allclose(
        (_numpy(hessian) + update) @ _numpy(step),
        _numpy(gradient_delta), rtol=tol, atol=tol,
    )
    assert np.linalg.norm(update) > 0.1


def _dense_ts_reference(name, hessian, step, gradient_delta):
    """Evaluate the published weight-matrix form in float64."""
    s = np.asarray(step, dtype=np.float64)
    y = np.asarray(gradient_delta, dtype=np.float64)
    hessian = np.asarray(hessian, dtype=np.float64)
    residual = y - hessian @ s
    eigenvalues, vectors = np.linalg.eigh(hessian)
    absolute_hessian = (vectors * np.abs(eigenvalues)) @ vectors.T
    if name == "ts_bfgs_update":
        absolute_product = absolute_hessian @ s
        weight = np.outer(y, y) + np.outer(absolute_product, absolute_product)
    else:
        phi = (residual @ s)**2 / ((s @ s) * (residual @ residual))
        if name == "ts_bfgs_update_org":
            coefficient = 1.0
        else:
            coefficient = abs(
                (y @ y - y @ hessian @ s) / ((y @ y) * (y @ s))
            )
        weight = (1.0 - phi) * absolute_hessian + coefficient * phi * np.outer(y, y)
    u = weight @ s / (s @ weight @ s)
    return (
        np.outer(residual, u) + np.outer(u, residual)
        - (residual @ s) * np.outer(u, u)
    )


@pytest.mark.parametrize("name", TS_UPDATES)
@pytest.mark.parametrize("scale", [1.0, 1e-9, 1e-12])
@pytest.mark.parametrize("case", ["mixed", "negative_curvature", "parallel"])
def test_ts_update_keeps_formula_and_tiny_secants(array_backend, name, scale, case):
    convert, tol = array_backend
    hessian = convert(np.diag([-1.0, 2.0]))
    if case == "mixed":
        step = convert(scale * np.array([0.25, 0.5]))
        gradient_delta = hessian @ step + convert(scale * np.array([0.125, 0.125]))
    elif case == "negative_curvature":
        step = convert(scale * np.array([0.5, 0.125]))
        gradient_delta = convert(scale * np.array([-0.75, 0.5]))
    else:
        step = convert(scale * np.array([0.5, 0.0]))
        gradient_delta = convert(scale * np.array([0.25, 0.0]))
    expected = _dense_ts_reference(
        name, _numpy(hessian), _numpy(step), _numpy(gradient_delta)
    )

    update, _ = _checked_update(name, hessian, step, gradient_delta)

    np.testing.assert_allclose(update, expected, rtol=tol, atol=tol)
    np.testing.assert_allclose(update, update.T, rtol=tol, atol=tol)
    # Scale the secant residual back to order one: a loose absolute tolerance
    # on the tiny original vectors would allow a zero correction to pass.
    np.testing.assert_allclose(
        (_numpy(hessian) + update) @ (_numpy(step) / scale),
        _numpy(gradient_delta) / scale, rtol=tol, atol=tol,
    )


@pytest.mark.parametrize("name", ["ts_bfgs_update", "ts_bfgs_update_org"])
def test_zero_curvature_ts_secant_is_valid(array_backend, name):
    convert, tol = array_backend
    hessian = convert(np.eye(2))
    step = convert([1.0, 0.0])
    gradient_delta = convert([0.0, 1.0])
    update, _ = _checked_update(name, hessian, step, gradient_delta)
    np.testing.assert_allclose(
        (_numpy(hessian) + update) @ _numpy(step),
        _numpy(gradient_delta), rtol=tol, atol=tol,
    )


def test_revised_zero_curvature_weight_is_explicit(array_backend):
    convert, _ = array_backend
    with pytest.raises(ValueError, match="singular.*nonzero secant"):
        hessian_updates.ts_bfgs_update_revised(
            convert(np.eye(2)), convert([1.0, 0.0]), convert([0.0, 1.0])
        )


@pytest.mark.parametrize("gradient_delta", [[2.0, -2.0], [0.0, 0.0]])
def test_revised_zero_mixing_weight_has_removable_singularity(
    array_backend, gradient_delta
):
    convert, tol = array_backend
    hessian = convert(np.diag([-1.0, 1.0]))
    step = convert([1.0, 1.0])
    gradient_delta = convert(gradient_delta)
    update, _ = _checked_update(
        "ts_bfgs_update_revised", hessian, step, gradient_delta
    )
    np.testing.assert_allclose(
        (_numpy(hessian) + update) @ _numpy(step),
        _numpy(gradient_delta), rtol=tol, atol=tol,
    )


@pytest.mark.parametrize("name", ("bfgs_update",) + TS_UPDATES)
@pytest.mark.parametrize("field", ["hessian", "step", "gradient_delta"])
def test_nonfinite_update_input_is_rejected(array_backend, name, field):
    convert, _ = array_backend
    values = {
        "hessian": np.diag([-1.0, 2.0]),
        "step": np.array([0.25, 0.5]),
        "gradient_delta": np.array([0.125, 1.125]),
    }
    values[field].flat[0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        getattr(hessian_updates, name)(
            *(convert(values[key]) for key in ("hessian", "step", "gradient_delta"))
        )


@pytest.mark.parametrize("case", ["indefinite_model", "zero_step", "zero_curvature"])
def test_bfgs_skips_zero_denominators(array_backend, case):
    convert, _ = array_backend
    step, gradient_delta = [1.0, 1.0], [1.0, 2.0]
    if case == "zero_step":
        step = [0.0, 0.0]
    elif case == "zero_curvature":
        step, gradient_delta = [1.0, 0.0], [0.0, 1.0]
    update, label = _checked_update(
        "bfgs_update", convert(np.diag([-1.0, 1.0])),
        convert(step), convert(gradient_delta),
    )
    np.testing.assert_array_equal(update, np.zeros((2, 2)))
    assert "skipped" in label


@pytest.mark.parametrize("diagonal", [[-1.0, 1.001], [1e-6, 2e-6]])
@pytest.mark.parametrize("scale", [1.0, 1e-12])
def test_bfgs_retains_nonzero_small_denominator(array_backend, diagonal, scale):
    convert, tol = array_backend
    hessian = convert(np.diag(diagonal))
    step = convert(scale * np.array([1.0, 1.0]))
    gradient_delta = convert(scale * np.array([1.0, 2.0]))
    s = _numpy(step).astype(np.float64)
    y = _numpy(gradient_delta).astype(np.float64)
    hs = _numpy(hessian).astype(np.float64) @ s
    expected = np.outer(y, y) / (s @ y) - np.outer(hs, hs) / (s @ hs)

    update, label = _checked_update("bfgs_update", hessian, step, gradient_delta)

    assert "skipped" not in label
    # The indefinite near-cancellation case is ill-conditioned in fp32.
    np.testing.assert_allclose(update, expected, rtol=tol * 5, atol=tol)


@pytest.mark.parametrize("old_size", [2, 3])
def test_basis_reset_discards_history_before_next_housekeeping(monkeypatch, old_size):
    module = importlib.import_module("pysisyphus.optimizers.HessianOptimizer")
    exact_hessian = np.diag([2.0, 3.0, 4.0])
    monkeypatch.setattr(
        module, "get_guess_hessian",
        lambda _geometry, _init: (exact_hessian.copy(), "exact"),
    )
    opt = object.__new__(RFOptimizer)
    opt.geometry = SimpleNamespace(
        coord_type="cart", is_analytical_2d=False, internal=None,
        cart_coords=np.zeros(3), gradient=np.zeros(3), energy=0.0,
        freeze_atoms=[], within_partial_hessian=None,
    )
    opt.hessian_init = "calc"
    opt.hessian_recalc_reset = True
    opt.hessian_recalc_adapt = None
    opt.hessian_recalc = None
    opt.hessian_recalc_in = None
    opt.adapt_norm = None
    opt.hessian_update = "bofill"
    opt.hessian_update_window = 2
    opt._using_active_dofs = False
    opt.trust_update = False
    opt.reject_uphill = False
    opt.small_eigval_thresh = 1e-8
    opt.cur_cycle = 1
    opt.log = lambda *_args: None
    old_step = np.eye(old_size)[0]
    opt._sy_buffer_S = [old_step.copy()]
    opt._sy_buffer_Y = [-old_step.copy()]

    opt.reset()

    assert opt._sy_buffer_S == []
    assert opt._sy_buffer_Y == []
    # First evaluation in the rebuilt basis cannot update yet.
    opt.coords = [np.zeros(3)]
    opt.energies = []
    opt.forces = []
    opt.steps = []
    opt.housekeeping()
    step = np.array([0.0, 0.25, 0.0])
    opt.coords.append(step.copy())
    opt.steps.append(step.copy())
    opt.geometry.cart_coords = step.copy()
    opt.geometry.gradient = exact_hessian @ step
    opt.geometry.energy = 0.5 * step @ opt.geometry.gradient

    opt.housekeeping()

    assert len(opt._sy_buffer_S) == len(opt._sy_buffer_Y) == 1
    np.testing.assert_array_equal(opt._sy_buffer_S[0], step)
    np.testing.assert_allclose(opt.H, exact_hessian, atol=1e-12)
