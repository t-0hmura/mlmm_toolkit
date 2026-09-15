"""Curvature-certified minimum and saddle searches on an exact constrained PES."""
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


@pytest.mark.parametrize("initial", ["calc", "unit"])
def test_minimum_resolves_real_curvature_with_exact_or_stale_model(tmp_path, initial):
    geom, opt = make_optimizer(
        tmp_path, RFOptimizer, [0., 0., 1.],
        hessian_init=initial, hessian_recalc=500,
    )
    frozen = geom.cart_coords[:9].copy()
    opt.run()
    assert opt.is_converged
    assert opt._minimum_curvature_valid
    np.testing.assert_array_equal(geom.cart_coords[:9], frozen)
    hessian = geom.calculator.results_at(geom.cart_coords)["hessian"][-3:, -3:]
    assert np.linalg.eigvalsh(hessian).min() > 0.
    np.testing.assert_allclose(np.abs(geom.cart_coords[-3:-1]), 1., atol=1e-4)


@pytest.mark.parametrize("y", [0., -.005, .005, 1.])
def test_rsprfo_reaches_first_order_without_flatten(tmp_path, y):
    geom, opt = make_optimizer(tmp_path, RSPRFOptimizer, [0., y, 1.])
    frozen = geom.cart_coords[:9].copy()
    opt.run()
    assert opt.is_converged
    assert opt._last_exact_n_imaginary == 1
    assert opt._exact_saddle_matches_current_geometry()
    np.testing.assert_array_equal(geom.cart_coords[:9], frozen)
    hessian = geom.calculator.results_at(geom.cart_coords)["hessian"][-3:, -3:]
    assert np.count_nonzero(np.linalg.eigvalsh(hessian) < 0.) == 1
    assert abs(geom.cart_coords[-3]) < 1e-4
    assert abs(abs(geom.cart_coords[-2])-1.) < 1e-4


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


def test_saved_positive_model_cannot_certify_a_physical_hosp(tmp_path):
    class ForceOnlyDoubleWell(DoubleWell):
        hessian_calls = 0

        def get_forces(self, atoms, coords, **kwargs):
            result = self.results_at(coords)
            result.pop("hessian")
            return result

        get_energy = get_forces

        def get_hessian(self, atoms, coords, **kwargs):
            self.hessian_calls += 1
            return self.results_at(coords)

    seed = tmp_path / "positive_model.dat"
    np.savetxt(seed, np.eye(12))
    geom, opt = make_optimizer(
        tmp_path, RFOptimizer, [0., 0., 1.],
        hessian_init=str(seed), hessian_recalc=500,
    )
    calculator = ForceOnlyDoubleWell(tmp_path)
    geom.set_calculator(calculator)
    opt.run()
    assert calculator.hessian_calls > 0
    assert opt.is_converged
    hessian = calculator.results_at(geom.cart_coords)["hessian"][-3:, -3:]
    assert np.linalg.eigvalsh(hessian).min() > 0.


def test_refreshed_adaptive_step_keeps_its_predictor(tmp_path):
    _, opt = make_optimizer(
        tmp_path, RFOptimizer, [0., 0., 1.], hessian_init="unit",
        hessian_recalc=500, adapt_step_func=True, max_cycles=1,
    )
    opt.run()
    step = opt.steps[0][-3:]
    expected = opt.rfo_model(np.zeros(3), np.diag([-4., -4., 10.]), step)
    assert opt.predicted_energy_changes[0] == pytest.approx(expected)


@pytest.mark.parametrize(
    "cadence, expected_calls, countdown",
    [(500, 2, 496), (3, 3, 2), (None, 6, None), (np.inf, 6, np.inf)],
)
@pytest.mark.parametrize("zero_acceleration", [False, True])
def test_negative_model_preserves_refresh_schedule(
    tmp_path, cadence, expected_calls, countdown, zero_acceleration
):
    class SoftDoubleWell(DoubleWell):
        hessian_calls = 0

        def results_at(self, coords):
            return {k: v * .01 for k, v in super().results_at(coords).items()}

        def get_forces(self, atoms, coords, **kwargs):
            result = self.results_at(coords)
            result.pop("hessian")
            return result

        get_energy = get_forces

        def get_hessian(self, atoms, coords, **kwargs):
            self.hessian_calls += 1
            return self.results_at(coords)

    geom, opt = make_optimizer(
        tmp_path, RFOptimizer, [0., 0., 1.], hessian_recalc=cadence,
        trust_radius=1e-4, trust_min=1e-4, trust_max=1e-4,
        trust_update=False, max_cycles=6, gdiis=False, line_search=False,
    )
    calculator = SoftDoubleWell(tmp_path)
    geom.set_calculator(calculator)
    if zero_acceleration:
        original = opt._verify_minimum_step

        def verify(step, *args):
            if opt.cur_cycle >= 2:
                step = np.zeros_like(step)
            return original(step, *args)

        opt._verify_minimum_step = verify
    opt.run()
    assert not opt.is_converged
    assert calculator.hessian_calls == expected_calls
    assert opt.hessian_recalc_in == countdown
    assert all(np.linalg.norm(step) > 0 for step in opt.steps)
    assert not opt._minimum_matches_current_geometry()


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
    assert opt._minimum_exact_coords is None


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


@pytest.mark.parametrize("soft,expected", [(-3.2, 1), (3.2, 0), (-0.0, 0), (0.0, 0)])
def test_opt_strict_count_rejects_soft_negative(soft, expected):
    opt = _strict_phva_packet_optimizer(RFOptimizer, [10., soft, 20.])
    assert opt._minimum_negative_count(np.ones(3)) == expected


@pytest.mark.parametrize("values,resolved,strict,accepted", [
    ([-100., -3.2, 12.], 1, 2, False),
    ([-100., 3.2, 12.], 1, 1, True),
    ([-100., -0., 12.], 1, 1, True),
    ([100., -3.2, 12.], 0, 1, False),
])
def test_exact_phva_strict_order_keeps_resolved_eligibility(values, resolved, strict, accepted):
    opt = _strict_phva_packet_optimizer(RSPRFOptimizer, values)
    has_modes, _, performed = opt._verify_exact_vibrational_structure(None, None)
    assert opt._last_exact_saddle_verified is accepted
    assert opt._last_exact_n_imaginary == resolved
    assert opt._last_exact_n_negative == strict
    assert has_modes is (resolved >= 1)
    assert performed is True
    assert opt._last_exact_validation == ("first_order" if accepted else "higher_order" if strict > 1 else "no_imaginary")


def test_invalid_partition_clears_exact_proof_and_cannot_certify_minimum():
    opt = _strict_phva_packet_optimizer(RSPRFOptimizer, [-100., 3.2, 12.], metadata=False)
    opt._verify_exact_vibrational_structure(None, None)
    assert opt._last_exact_saddle_verified is False
    assert opt._last_exact_n_negative is None
    assert opt._last_exact_n_imaginary is None
    assert opt._last_exact_validation == "unavailable"
    minimum = _strict_phva_packet_optimizer(RFOptimizer, [10., 3.2, 12.], metadata=False)
    with pytest.raises(ValueError, match="partition"):
        minimum._minimum_negative_count(np.ones(3))
