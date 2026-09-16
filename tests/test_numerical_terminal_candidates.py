"""Numerical termination never implicitly launches curvature recovery.

These small force-only oracles distinguish an ordinary numerical candidate
from its physical Hessian classification. Explicit recovery is a separate,
caller-selected feature and is tested as such.
"""
import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from pysisyphus.tsoptimizers.TRIM import TRIM


CENTER = np.array([0., 0., 0., 2., 0., 0., 0., 2., 0., 0., 0., 1.])


class Quadratic(Calculator):
    def __init__(self, diagonal, out_dir, *, forbid_hessian=False):
        super().__init__(out_dir=out_dir, check_mem=False)
        self.diagonal = np.asarray(diagonal, dtype=float)
        self.forbid_hessian = forbid_hessian
        self.hessian_points = []
        self.force_points = []

    def get_forces(self, atoms, coords, **kwargs):
        self.force_points.append(np.asarray(coords).copy())
        delta = np.asarray(coords)[-3:] - CENTER[-3:]
        gradient = np.zeros(12)
        gradient[-3:] = self.diagonal * delta
        return dict(energy=float(.5 * delta @ gradient[-3:]), forces=-gradient)

    get_energy = get_forces

    def get_hessian(self, atoms, coords, **kwargs):
        self.hessian_points.append(np.asarray(coords).copy())
        if self.forbid_hessian:
            raise AssertionError("force-only RFO must not request a Hessian")
        hessian = np.zeros((12, 12))
        hessian[-3:, -3:] = np.diag(self.diagonal)
        return dict(self.get_forces(atoms, coords), hessian=hessian)


def make_geometry(calculator, displacement=0.):
    coords = CENTER.copy()
    coords[-3] += displacement
    geom = Geometry(["H"] * 4, coords, coord_type="cart", freeze_atoms=[0, 1, 2])
    geom.set_calculator(calculator)
    return geom


@pytest.mark.parametrize("displacement", [0., .02])
@pytest.mark.parametrize("flatten", [False, True])
def test_force_only_rfo_reaches_numerical_candidate_without_exact_hessian(
    tmp_path, displacement, flatten
):
    calc = Quadratic([1., 1., 1.], tmp_path, forbid_hessian=True)
    geom = make_geometry(calc, displacement)
    opt = RFOptimizer(
        geom, hessian_init="unit", hessian_recalc=None, max_cycles=20,
        trust_radius=.1, trust_max=.1, thresh="baker", gdiis=False,
        line_search=False, flatten_enabled=flatten, dump=False, out_dir=tmp_path,
    )
    opt.run()
    assert opt.is_converged
    assert calc.hessian_points == []
    np.testing.assert_array_equal(geom.cart_coords[:9], CENTER[:9])
    assert abs(geom.cart_coords[-3] - CENTER[-3]) < 1e-4


def test_stale_rfo_model_returns_candidate_without_implicit_curvature_escape(tmp_path):
    calc = Quadratic([-4., -2., 10.], tmp_path, forbid_hessian=True)
    geom = make_geometry(calc)
    opt = RFOptimizer(
        geom, hessian_init="unit", hessian_recalc=500, max_cycles=5,
        thresh="baker", dump=False, out_dir=tmp_path,
    )
    opt.run()
    assert opt.is_converged  # Numerical candidate, not a certified minimum.
    assert calc.hessian_points == []
    np.testing.assert_array_equal(geom.cart_coords, CENTER)
    assert len(opt.steps) == 1
    np.testing.assert_array_equal(opt.steps[0], np.zeros(12))


def make_ts(tmp_path, cls, diagonal, *, flatten=False, recovery=0):
    calc = Quadratic(diagonal, tmp_path)
    geom = make_geometry(calc)
    model = np.eye(12)
    model[-3, -3] = -1.
    path = tmp_path / "initial_model.dat"
    np.savetxt(path, model)
    opt = cls(
        geom, hessian_init=str(path), hessian_recalc=500, max_cycles=5,
        reference_mode=np.r_[np.zeros(9), 1., 0., 0.],
        trust_radius=.1, trust_max=.1, thresh="baker",
        flatten_enabled=flatten, saddle_recovery_max_cycles=recovery,
        dump=False, out_dir=tmp_path,
    )
    return geom, opt, calc


@pytest.mark.parametrize("cls", [RSPRFOptimizer, RSIRFOptimizer, TRIM])
@pytest.mark.parametrize("diagonal", [(-4., -2., 10.), (-4., -1e-7, 10.), (-4., 1e-7, 10.), (-4., 2., 10.)])
@pytest.mark.parametrize("flatten", [False, True])
def test_terminal_ts_phva_does_not_displace_numerical_candidate(
    tmp_path, cls, diagonal, flatten
):
    geom, opt, calc = make_ts(tmp_path, cls, diagonal, flatten=flatten)
    opt.run()
    assert opt.is_converged
    assert opt.saddle_recovery_max_cycles == 0
    assert opt.saddle_recovery_steps == 0
    assert not opt._saddle_recovery_active
    assert len(calc.hessian_points) == 1
    np.testing.assert_array_equal(calc.hessian_points[0], CENTER)
    np.testing.assert_array_equal(geom.cart_coords, CENTER)
    for point in calc.force_points:
        np.testing.assert_array_equal(point, CENTER)
    assert len(opt.steps) == 1
    np.testing.assert_array_equal(opt.steps[0], np.zeros(12))
    if diagonal[1] == -2.:
        assert opt._last_exact_n_imaginary == 2
        assert not opt._exact_saddle_matches_current_geometry()
    elif diagonal[1] > 0.:
        assert opt._exact_saddle_matches_current_geometry()
    # Retain the signed spectrum independently of the criterion used to
    # classify the soft root; neither criterion may initiate implicit motion.
    assert np.count_nonzero(opt._last_exact_frequencies_cm < 0.) == (2 if diagonal[1] < 0. else 1)


@pytest.mark.parametrize("cls", [RSPRFOptimizer, RSIRFOptimizer, TRIM])
def test_default_ts_stops_without_motion_when_exact_hessian_has_no_negative_mode(tmp_path, cls):
    geom, opt, calc = make_ts(tmp_path, cls, [1., 2., 3.])
    opt.run()
    assert not opt.is_converged
    assert opt.stopped
    assert opt.saddle_recovery_max_cycles == 0
    assert opt.saddle_recovery_steps == 0
    assert len(calc.hessian_points) == 1
    np.testing.assert_array_equal(geom.cart_coords, CENTER)


def test_explicit_legacy_recovery_budget_remains_an_opt_in(tmp_path):
    geom, opt, calc = make_ts(tmp_path, RSPRFOptimizer, [1., 2., 3.], recovery=2)
    opt.max_cycles = 1
    opt.run()
    assert not opt.is_converged
    assert opt.saddle_recovery_max_cycles == 2
    assert opt.saddle_recovery_steps == 1
    assert opt._saddle_recovery_active
    assert len(calc.hessian_points) == 1
    assert np.linalg.norm(geom.cart_coords - CENTER) > 0.
