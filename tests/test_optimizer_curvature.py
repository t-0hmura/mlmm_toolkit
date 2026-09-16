"""Ordinary curvature steps and separately evaluated saddle diagnostics."""
import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer


class DoubleWell(Calculator):
    """Two double wells on one atom relative to three frozen anchors."""

    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    def results_at(self, coords):
        x, y, z = np.asarray(coords)[-3:] - [0.0, 0.0, 1.0]
        gradient = np.zeros(12)
        gradient[-3:] = [4*x*(x*x-1), 4*y*(y*y-1), 10*z]
        hessian = np.zeros((12, 12))
        hessian[-3:, -3:] = np.diag([12*x*x-4, 12*y*y-4, 10])
        return {
            "energy": float((x*x-1)**2 + (y*y-1)**2 + 5*z*z),
            "forces": -gradient,
            "hessian": hessian,
        }

    def get_energy(self, atoms, coords, **kwargs):
        return self.results_at(coords)

    def get_forces(self, atoms, coords, **kwargs):
        return self.results_at(coords)

    def get_hessian(self, atoms, coords, **kwargs):
        return self.results_at(coords)


def make_optimizer(tmp_path, kind, position, **kwargs):
    coords = np.r_[0., 0., 0., 2., 0., 0., 0., 2., 0., position]
    geom = Geometry(["H"]*4, coords, coord_type="cart", freeze_atoms=[0, 1, 2])
    geom.set_calculator(DoubleWell(tmp_path))
    settings = dict(hessian_init="calc", hessian_recalc=1, trust_radius=.1,
                    trust_max=.1, thresh="baker", max_cycles=100,
                    out_dir=tmp_path, dump=False)
    settings.update(kwargs)
    if kind is RSPRFOptimizer:
        settings["reference_mode"] = np.r_[np.zeros(9), 1., 0., 0.]
    return geom, kind(geom, **settings)


def test_minimum_ordinary_step_uses_available_negative_curvature(tmp_path):
    geom, opt = make_optimizer(
        tmp_path, RFOptimizer, [0., 0., 1.], hessian_recalc=500,
        max_cycles=1, trust_update=False, gdiis=False, line_search=False,
    )
    before = geom.cart_coords.copy()
    energy = geom.calculator.results_at(before)["energy"]
    opt.run()
    assert not opt.is_converged
    np.testing.assert_array_equal(geom.cart_coords[:9], before[:9])
    assert 0. < np.linalg.norm(opt.steps[0]) <= opt.trust_radius * (1 + 1e-12)
    assert geom.calculator.results_at(geom.cart_coords)["energy"] < energy


def test_full_hessian_reference_root_uses_active_eigenspace(tmp_path):
    geom, opt = make_optimizer(tmp_path, RSPRFOptimizer, [1., 1., 1.])
    opt.prepare_opt()
    assert opt.using_active_dofs
    assert opt.cur_H.shape == (3, 3)
    assert 0 <= opt.root < 3
    assert opt.ts_modes.shape == (1, 3)


def test_old_saddle_mode_does_not_reflect_new_physical_negative_direction(tmp_path):
    geom, opt = make_optimizer(
        tmp_path, RSPRFOptimizer, [0., 1., 1.], trust_update=False,
        min_line_search=False, max_line_search=False,
    )
    opt.prepare_opt()
    opt.coords = [geom.coords.copy()]
    opt.cart_coords = [geom.cart_coords.copy()]
    opt.housekeeping()
    opt._validate_terminal_exact_saddle()
    assert opt._last_exact_n_imaginary == 1
    assert opt._physical_ts_mode is not None
    certificate = opt._last_exact_cart_coords.copy()
    exact_checks = opt.exact_saddle_checks

    # Move the fixture to a different point, then take the ordinary cadence
    # refresh. This is a lifecycle boundary test, not a saved optimizer restart.
    current = geom.coords.copy()
    current[-3:] = [.005, 0., 1.]
    geom.coords = current
    opt.steps.append(geom.coords - opt.coords[-1])
    opt.coords.append(geom.coords.copy())
    opt.cart_coords.append(geom.cart_coords.copy())
    opt.cur_cycle += 1
    step = opt.optimize()

    np.testing.assert_array_equal(opt._last_exact_cart_coords, certificate)
    assert not opt._exact_phva_matches_current_geometry()
    assert opt.exact_saddle_checks == exact_checks
    assert np.count_nonzero(np.linalg.eigvalsh(opt.cur_H) < 0.) == 2
    # Three noncollinear frozen anchors leave no feasible rigid displacement.
    # The second negative direction is a true double-well maximum, with zero
    # gradient: it must retain the existing uncoupled image-step escape.
    assert abs(step[-2]) > .01
    assert np.linalg.norm(step) <= opt.trust_radius * (1. + 1e-12)
    np.testing.assert_array_equal(step[:9], np.zeros(9))


def test_saved_model_does_not_replace_an_independent_physical_hessian(tmp_path):
    from pysisyphus.optimizers.guess_hessians import get_guess_hessian

    seed = tmp_path / "positive_model.dat"
    np.savetxt(seed, np.eye(12))
    geom, _ = make_optimizer(tmp_path, RFOptimizer, [0., 0., 1.])
    model, _ = get_guess_hessian(geom, str(seed))
    np.testing.assert_array_equal(model, np.eye(12))
    assert geom._hessian is None
    physical = geom.hessian
    np.testing.assert_array_equal(physical[-3:, -3:], np.diag([-4., -4., 10.]))
    np.testing.assert_array_equal(model, np.eye(12))


def test_explicit_ts_flatten_receives_hosp_candidate(tmp_path):
    _, opt = make_optimizer(
        tmp_path, RSPRFOptimizer, [0., 0., 1.],
        flatten_enabled=True, max_cycles=1,
    )
    opt.run()
    assert opt.is_converged
    assert opt._last_exact_n_imaginary == 2
    assert not opt._exact_saddle_matches_current_geometry()


def test_explicit_minimum_flatten_owns_curvature_check(tmp_path):
    _, opt = make_optimizer(
        tmp_path, RFOptimizer, [0., 0., 1.], hessian_init="unit",
        hessian_recalc=500, flatten_enabled=True, max_cycles=1,
    )
    opt.run()
    assert opt.is_converged
    np.testing.assert_array_equal(opt.steps[-1], np.zeros(12))


@pytest.mark.parametrize("position", [[.005, 0., 1.], [1., .005, 1.]])
@pytest.mark.parametrize("line_offset", [0., .012])
def test_uncoupled_partition_uses_physical_image_step(tmp_path, position, line_offset):
    _, opt = make_optimizer(
        tmp_path, RSPRFOptimizer, position, max_cycles=1,
    )
    opt.step_and_grad_from_line_search = (
        lambda energy, gradient, *args:
        (np.full_like(gradient, line_offset), gradient)
    )
    captured = []
    original = opt._image_trust_step

    def image_step():
        step, gradient = original()
        captured.append((step.copy(), gradient.copy(), opt.cur_H.copy()))
        return step, gradient

    opt._image_trust_step = image_step
    opt.run()
    assert not opt.is_converged
    assert len(opt.steps) == len(opt.predicted_energy_changes) == len(captured) == 1
    step, gradient, hessian = captured[0]
    np.testing.assert_allclose(opt.steps[0][-3:], step, atol=1e-14)
    assert 0 < np.linalg.norm(step) <= .1 * (1 + 1e-12)
    expected = gradient @ step + .5 * step @ hessian @ step
    assert opt.predicted_energy_changes[0] == pytest.approx(expected)


def _strict_phva_packet_optimizer(kind, values, metadata=True):
    """A complete filtered PHVA packet; no calculator or trajectory is fabricated."""
    from types import SimpleNamespace
    from pysisyphus.normal_modes import filter_resolved_modes

    info = {}
    frequencies, modes = filter_resolved_modes(
        np.array(values), np.eye(len(values)), 5., filter_info=info,
    )
    opt = kind.__new__(kind)
    opt._mw_frequencies_and_modes = lambda: (frequencies, modes)
    opt._last_rigid_projection_info = info if metadata else {}
    opt._recovery_mode_from_mw = lambda _modes, index: _modes[index].copy()
    opt.geometry = SimpleNamespace(cart_coords=np.zeros(len(values)))
    opt.reference_mode = None
    opt.roots = [0]
    opt.saddle_imaginary_threshold_cm = 5.
    opt.cur_cycle = 7
    opt.table = SimpleNamespace(print=lambda *_a, **_kw: None)
    opt.request_stop = lambda *_a: None
    opt._record_exact_saddle_candidate = lambda: None
    opt._last_exact_saddle_verified = True  # failed packets must clear a past proof
    return opt


@pytest.mark.parametrize("values,resolved,strict,accepted", [
    ([-100., -3.2, 12.], 1, 2, True),
    ([-100., 3.2, 12.], 1, 1, True),
    ([-100., -0., 12.], 1, 1, True),
    ([100., -3.2, 12.], 0, 1, False),
])
def test_exact_phva_reports_raw_count_and_uses_selected_eligibility(values, resolved, strict, accepted):
    opt = _strict_phva_packet_optimizer(RSPRFOptimizer, values)
    has_modes, _, performed = opt._verify_exact_vibrational_structure(None, None)
    assert opt._last_exact_saddle_verified is accepted
    assert opt._last_exact_n_imaginary == resolved
    assert opt._last_exact_n_negative == strict
    assert has_modes is (resolved >= 1)
    assert performed is True
    assert opt._last_exact_validation == ("first_order" if accepted else "higher_order" if resolved > 1 else "no_imaginary")


def test_invalid_partition_clears_exact_saddle_proof():
    opt = _strict_phva_packet_optimizer(RSPRFOptimizer, [-100., 3.2, 12.], metadata=False)
    opt._verify_exact_vibrational_structure(None, None)
    assert opt._last_exact_saddle_verified is False
    assert opt._last_exact_n_negative is None
    assert opt._last_exact_n_imaginary is None
    assert opt._last_exact_validation == "unavailable"
