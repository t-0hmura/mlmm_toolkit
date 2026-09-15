"""Structured RFO eigenproblem and restricted-step failure contracts."""
import numpy as np
import pytest
import torch

from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer

BACKENDS = ["numpy", "cpu"] + (["cuda"] if torch.cuda.is_available() else [])


def optimizer():
    opt = RFOptimizer.__new__(RFOptimizer)
    opt.log = lambda *_: None
    return opt


def array(values, backend, dtype=np.float64):
    if backend == "numpy":
        return np.asarray(values, dtype=dtype)
    dt = torch.float64 if dtype == np.float64 else torch.float32
    return torch.tensor(values, dtype=dt, device=backend)


def numpy(value):
    return value.detach().cpu().numpy() if isinstance(value, torch.Tensor) else value


def original_matrix(lam, gradient, alpha):
    result = np.zeros((len(lam) + 1, len(lam) + 1))
    result[:-1, :-1] = np.diag(lam / alpha)
    result[:-1, -1] = gradient / alpha
    result[-1, :-1] = gradient
    return result


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("dtype", [np.float64, np.float32])
@pytest.mark.parametrize("alpha", [0.2, 1.0, 4.0, 1000.0])
@pytest.mark.parametrize("kind", ["min", "max"])
def test_symmetric_rfo_matches_original_problem(backend, dtype, alpha, kind):
    opt = optimizer()
    lam = array([-0.2, 0.3, 0.7], backend, dtype)
    gradient = array([0.12, -0.03, 0.07], backend, dtype)
    matrix = opt.get_augmented_hessian(lam, gradient, alpha)
    before = numpy(matrix).copy()
    step, root, nu, vector = opt.solve_rfo(matrix, kind, alpha=alpha)
    original = original_matrix(numpy(lam), numpy(gradient), alpha)
    roots, vectors = np.linalg.eig(original)
    index = np.argmin(roots) if kind == "min" else np.argmax(roots)
    expected = vectors[:-1, index] / vectors[-1, index]
    tol = 3e-5 if dtype == np.float32 else 1e-11
    if dtype == np.float32:
        # Eigensolvers control normwise error, not relative error in each
        # near-zero component of an eigenvector divided by its last entry.
        assert np.linalg.norm(numpy(step) - expected) / np.linalg.norm(expected) <= 32 * np.finfo(dtype).eps
    else:
        np.testing.assert_allclose(numpy(step), expected, rtol=tol, atol=tol * 0.01)
    np.testing.assert_allclose(float(root), roots[index], rtol=tol, atol=tol * 0.01)
    residual = original @ numpy(vector) - float(root) * numpy(vector)
    assert np.linalg.norm(residual) / np.linalg.norm(original) <= 8 * np.finfo(dtype).eps
    np.testing.assert_allclose(np.linalg.norm(numpy(vector)), 1.0, rtol=tol)
    np.testing.assert_array_equal(numpy(matrix), before)
    assert float(nu) != 0.0
    if backend != "numpy":
        assert step.dtype == lam.dtype and step.device == lam.device


@pytest.mark.parametrize("backend", BACKENDS)
def test_small_coupling_uses_exact_symmetric_problem(backend):
    opt = optimizer()
    lam = array([0.1, 0.2, 0.5], backend)
    gradient = array([1e-9, -2e-9, 1e-9], backend)
    alpha = 4.0
    matrix = opt.get_augmented_hessian(lam, gradient, alpha)
    step = opt.solve_rfo(matrix, alpha=alpha)[0]
    # The correction to the Newton limit is O(g**3).
    np.testing.assert_allclose(numpy(step), -numpy(gradient) / numpy(lam),
                               rtol=1e-6, atol=1e-15)


@pytest.mark.parametrize("backend", BACKENDS)
def test_overlap_uses_backmapped_normalized_vectors(backend):
    opt = optimizer()
    lam = array([-0.2, 0.3, 0.7], backend)
    gradient = array([0.12, -0.03, 0.07], backend)
    alpha = 100.0
    original = original_matrix(numpy(lam), numpy(gradient), alpha)
    roots, vectors = np.linalg.eig(original)
    # Arbitrary previous directions test actual overlap selection, not only
    # recovering an eigenvector that was selected in advance.
    for previous in np.random.default_rng(8).normal(size=(12, 4)):
        previous /= np.linalg.norm(previous)
        selected = np.argmax(np.abs(previous @ vectors))
        result = opt.solve_rfo(
            opt.get_augmented_hessian(lam, gradient, alpha), "min",
            prev_eigvec=array(previous, backend), alpha=alpha,
        )
        assert float(result[1]) == pytest.approx(roots[selected], rel=1e-10, abs=1e-12)
        assert opt.solve_rfo_secular(lam, gradient, alpha,
                                    prev_eigvec=array(previous, backend)) is None


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("kind,lam", [("min", [-1.0]), ("max", [1.0])])
def test_uncoupled_adverse_root_is_not_false_convergence(backend, kind, lam):
    opt = optimizer()
    values, gradient = array(lam, backend), array([0.0], backend)
    assert opt.solve_rfo_secular(values, gradient, 4.0, kind=kind) is None
    with pytest.raises(ZeroDivisionError, match="zero augmented component"):
        opt.solve_rfo(opt.get_augmented_hessian(values, gradient, 4.0),
                      kind, alpha=4.0)


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("kind,lam", [("min", [1.0]), ("max", [-1.0]),
                                     ("min", [0.0]), ("max", [0.0])])
def test_zero_step_and_degenerate_finite_root_are_valid(backend, kind, lam):
    opt = optimizer()
    values, gradient = array(lam, backend), array([0.0], backend)
    secular = opt.solve_rfo_secular(values, gradient, 4.0, kind=kind)
    dense = opt.solve_rfo(opt.get_augmented_hessian(values, gradient, 4.0),
                         kind, alpha=4.0)
    assert secular is not None
    np.testing.assert_array_equal(numpy(dense[0]), [0.0])


@pytest.mark.parametrize("alpha", [0.0, -1.0, np.inf, np.nan])
def test_invalid_alpha_is_an_explicit_error(alpha):
    opt = optimizer()
    with pytest.raises(ValueError, match="alpha"):
        opt.get_augmented_hessian(np.array([1.0]), np.array([0.1]), alpha)
    with pytest.raises(ValueError, match="alpha"):
        opt.solve_rfo(np.eye(2), alpha=alpha)
    with pytest.raises(ValueError, match="alpha"):
        opt.solve_rfo_secular(np.array([1.0]), np.array([0.1]), alpha)


def test_nonfinite_and_nonsymmetric_inputs_are_not_sanitized():
    opt = optimizer()
    with pytest.raises(ValueError, match="NaN/inf"):
        opt.solve_rfo(np.array([[1.0, np.nan], [np.nan, 0.0]]))
    with pytest.raises(ValueError, match="symmetric"):
        opt.solve_rfo(np.array([[1.0, 1e-10], [2e-10, 0.0]]))


def test_dense_solver_does_not_use_general_eig(monkeypatch):
    opt = optimizer()
    def forbidden(*_args, **_kwargs):
        raise AssertionError("general eig was called")
    monkeypatch.setattr(np.linalg, "eig", forbidden)
    opt.solve_rfo(opt.get_augmented_hessian(np.array([1.0]), np.array([0.1]), 4.0),
                  alpha=4.0)
    def failed(*_args, **_kwargs):
        raise np.linalg.LinAlgError("eigh failed")
    monkeypatch.setattr(np.linalg, "eigh", failed)
    with pytest.raises(np.linalg.LinAlgError, match="eigh failed"):
        opt.solve_rfo(np.eye(2))


def test_near_pole_secular_failure_has_dense_solution():
    opt = optimizer()
    lam, gradient = np.array([-1.0]), np.array([1e-5])
    assert opt.solve_rfo_secular(lam, gradient, 4.0) is None
    step, root, _, vector = opt.solve_rfo(
        opt.get_augmented_hessian(lam, gradient, 4.0), alpha=4.0)
    assert np.isfinite(step).all()
    np.testing.assert_allclose(original_matrix(lam, gradient, 4.0) @ vector,
                               root * vector, atol=1e-15)


@pytest.mark.parametrize("cycles", [0, 25])
def test_minimum_rfo_preserves_existing_negative_curvature_step(cycles):
    opt = optimizer()
    opt.alpha0 = 1.0
    opt.max_micro_cycles = cycles
    opt.trust_radius = 0.1
    opt._prev_eigvec_min = None
    step = opt.get_rs_step(np.array([-1.0, 2.0]), np.eye(2), np.zeros(2))
    assert np.isfinite(step).all()
    assert np.linalg.norm(step) == pytest.approx(opt.trust_radius)
    assert abs(step[0]) == pytest.approx(opt.trust_radius)


def rs_optimizer(cycles=25):
    opt = RSPRFOptimizer.__new__(RSPRFOptimizer)
    opt.log = lambda *_: None
    opt.alpha0 = 1.0
    opt.max_micro_cycles = cycles
    opt.trust_radius = 0.1
    opt.rfo_overlaps = False
    opt._prev_eigvec_min = opt._prev_eigvec_max = None
    opt._physical_ts_mode = None
    opt.stop_requested = False
    opt.flatten_enabled = False
    opt.verify_saddle = False  # This fixture isolates the restricted step solver.
    opt._last_exact_cart_coords = None
    opt._last_exact_n_imaginary = None
    opt._last_exact_n_negative = None
    opt._last_exact_frequencies_cm = None
    opt.roots = [0]
    opt.cur_H = np.diag([-0.2, 0.3])
    opt.housekeeping = lambda: (0.0, np.array([0.1, 0.12]), opt.cur_H,
                               np.array([-0.2, 0.3]), np.eye(2), False)
    opt.update_ts_mode = lambda *_: None
    opt.step_and_grad_from_line_search = lambda e, g, *args: (np.zeros_like(g), g)
    opt.validate_terminal_saddle_for_step = lambda *_: None
    opt.apply_saddle_recovery_step = lambda step: step
    opt.full_from_active = lambda step: step
    opt.predicted_energy_changes = []
    return opt


def test_restricted_step_converges_inside_trust():
    opt = rs_optimizer()
    step = opt.optimize()
    assert np.isfinite(step).all()
    assert np.linalg.norm(step) <= opt.trust_radius * (1 + 1e-12)
    assert len(opt.predicted_energy_changes) == 1


def test_restricted_step_accepts_roundoff_at_trust_boundary():
    opt = rs_optimizer()
    def at_boundary(values, gradient, alpha, kind, **kwargs):
        length = np.nextafter(opt.trust_radius, np.inf) if kind == "max" else 0.0
        return np.array([length]), 0.0, 1.0, np.array([0.0, 1.0])
    opt.solve_rfo_secular = at_boundary
    step = opt.optimize()
    assert np.linalg.norm(step) == np.nextafter(opt.trust_radius, np.inf)


@pytest.mark.parametrize("values,gradient,trust", [
    ([-4.548177703484925, 1.7539536727302055e-7],
     [-2.4395511966045625e-12, 1.2141663190912835e-5], 0.007296478720777277),
    ([-4.018913834147417, 0.039037139803876374],
     [0.7598286922625115, 0.002919111026902365], 0.00010459214082846993),
    ([-3.6298388411081137, 5.988358245562869e-7, 2.2724789392745173],
     [-5.1581915837608926e-8, -6.598007880889851e-6, -0.004636282165313381],
     0.000800031945589288),
])
def test_restricted_step_preserves_near_boundary_cases(values, gradient, trust):
    opt = rs_optimizer(cycles=50)
    values, gradient = np.asarray(values), np.asarray(gradient)
    opt.cur_H = np.diag(values)
    opt.housekeeping = lambda: (0.0, gradient, opt.cur_H, values,
                               np.eye(values.size), False)
    opt.trust_radius = trust
    step = opt.optimize()
    assert np.isfinite(step).all()
    assert np.linalg.norm(step) <= trust * (1 + 1e-12)


def test_one_cycle_prfo_keeps_its_explicit_scaling_mode():
    opt = rs_optimizer(cycles=1)
    step = opt.optimize()
    assert np.linalg.norm(step) == pytest.approx(opt.trust_radius)


def test_zero_cycle_rsprfo_is_rejected():
    with pytest.raises(ValueError, match="at least one"):
        rs_optimizer(cycles=0).optimize()


def test_exhausted_rsprfo_does_not_accept_outside_step():
    opt = rs_optimizer(cycles=2)
    with pytest.raises(ValueError, match="exhausted"):
        opt.optimize()
    assert opt.predicted_energy_changes == []


@pytest.mark.parametrize("derivative", [0.0, np.nan, np.inf, 0.01])
def test_invalid_rsprfo_update_does_not_advance_coordinates(derivative):
    opt = rs_optimizer()
    opt._partition_dstep2_dalpha = lambda *_: derivative
    with pytest.raises(ValueError, match="derivative|alpha update"):
        opt.optimize()
    assert opt.predicted_energy_changes == []
