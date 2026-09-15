"""Near-zero physical curvature must be resolved by motion, not only rejected.

The force-only oracle has two independent wells on one movable H and three
noncollinear frozen H anchors. A supported, explicitly stale initial-Hessian
file makes the numerical terminal gate acquire the true soft-negative H.
No exact certificate, optimizer history, reference swap or extra-H hook is used.
These trajectories discriminate the demonstrated resolved1/strict2 acceptance
gap; they do not establish convergence on every molecular HOSP.
"""

import numpy as np
import pytest
import torch

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.normal_modes import _frequencies_cm_and_modes
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer


CENTER = np.array([0., 0., 0., 2., 0., 0., 0., 2., 0., 0., 0., 1.])
FROZEN = [0, 1, 2]
MAX_CYCLES = 200
CRITERIA = ("energy_converged", "max_force_converged", "rms_force_converged",
            "max_step_converged", "rms_step_converged")


def true_phva(hessian, coords):
    """Use the native mass/conversion kernel, independently of acceptance state."""
    info = {}
    resolved, modes = _frequencies_cm_and_modes(
        torch.tensor(hessian, dtype=torch.float64), [1] * 4,
        np.asarray(coords).reshape(-1, 3), torch.device("cpu"),
        freeze_idx=FROZEN, tr_projection="constrained", projection_info=info,
    )
    assert info["frequency_zero_cutoff_cm"] == 5.
    assert info["effective_rank"] == 0 and info["raw_mode_count"] == 3
    assert modes.shape == (len(resolved), 12)
    assert torch.count_nonzero(modes[:, :9]).item() == 0
    near = np.asarray(info["near_zero_frequencies_cm"])
    complete = np.sort(np.r_[resolved, near])
    assert complete.shape == (3,) and np.isfinite(complete).all()
    return complete, resolved, near


@pytest.fixture
def soft_curvature():
    # Calibrate from the actual H mass and native units, not a guessed Hessian
    # threshold. Frequency scales with sqrt(curvature) for this diagonal mode.
    unit_h = np.zeros((12, 12))
    unit_h[-3:, -3:] = np.diag([-1., 1., 2.])
    frequencies, _, _ = true_phva(unit_h, CENTER)
    curvature = (3.2 / abs(frequencies[0])) ** 2
    assert curvature > 1e-8  # retained by the native raw-H small-eigenvalue filter
    return curvature


class SoftDoubleWell(Calculator):
    """E=.01(x²-1)² + k/4(y²-1)² + 5z²; positive control replaces y well."""

    def __init__(self, out_dir, curvature, positive_near=False):
        super().__init__(out_dir=out_dir, check_mem=False)
        self.curvature = curvature
        self.positive_near = positive_near
        self.hessian_calls = []
        self.optimizer = None

    def get_forces(self, atoms, coords, **kwargs):
        x, y, z = np.asarray(coords)[-3:] - CENTER[-3:]
        if self.positive_near:
            soft_energy, soft_gradient = .5 * self.curvature * y*y, self.curvature * y
        else:
            soft_energy = .25 * self.curvature * (y*y - 1.)**2
            soft_gradient = self.curvature * y * (y*y - 1.)
        gradient = np.zeros(12)
        gradient[-3:] = [.04 * x * (x*x - 1.), soft_gradient, 10. * z]
        # Deliberately no hidden H in either ordinary E/F entry point.
        return {"energy": float(.01 * (x*x - 1.)**2 + soft_energy + 5. * z*z),
                "forces": -gradient}

    get_energy = get_forces

    def hessian_at(self, coords):
        x, y, _ = np.asarray(coords)[-3:] - CENTER[-3:]
        soft = self.curvature * (1. if self.positive_near else 3.*y*y - 1.)
        hessian = np.zeros((12, 12))
        hessian[-3:, -3:] = np.diag([.04 * (3.*x*x - 1.), soft, 10.])
        return hessian

    def get_hessian(self, atoms, coords, **kwargs):
        self.hessian_calls.append((int(self.optimizer.cur_cycle), np.asarray(coords).copy()))
        return {**self.get_forces(atoms, coords),
                "hessian": self.hessian_at(coords)}


def run_case(tmp_path, curvature, *, minimum=False, y=0., initial="stale",
             positive_near=False):
    coords = CENTER.copy()
    coords[-3:] += [1. if minimum else 0., y, 0.]
    geom = Geometry(["H"] * 4, coords, coord_type="cart", freeze_atoms=FROZEN)
    calc = SoftDoubleWell(tmp_path, curvature, positive_near)
    geom.set_calculator(calc)
    for method in (calc.get_forces, calc.get_energy):
        assert set(method(geom.atoms, coords)) == {"energy", "forces"}

    initial_h = calc.hessian_at(coords)
    frequencies, resolved, near = true_phva(initial_h, coords)
    assert len(near) == 1 and 0. < abs(near[0]) < 5.
    assert np.count_nonzero(frequencies < 0.) == (int(not minimum) + int(not positive_near))
    assert np.count_nonzero(resolved < 0.) == int(not minimum)
    if y == 0.:
        assert abs(near[0]) == pytest.approx(3.2, rel=1e-12)

    hessian_init = "calc"
    if initial == "stale":
        # Keep the correct target curvature, but hide the soft instability in
        # the initial MODEL only. The analytic PES/terminal H are unchanged.
        stale = initial_h.copy()
        stale[-2, -2] = .04
        seed = tmp_path / "stale_initial_B.txt"
        np.savetxt(seed, stale)
        hessian_init = str(seed)
    else:
        assert initial == "calc"

    settings = dict(
        hessian_init=hessian_init, hessian_recalc=500,
        hessian_recalc_adapt=None, hessian_xtb=False, hessian_update_window=1,
        trust_radius=.1, trust_min=1e-4, trust_max=.1, trust_update=True,
        thresh="baker", max_cycles=MAX_CYCLES, flatten_enabled=False,
        out_dir=tmp_path, dump=False,
    )
    # Start at .1, not a tiny radius chosen to manufacture convergence.
    # Stale-B cases reach the terminal gate through their small model step.
    if minimum:
        opt = RFOptimizer(geom, **settings)
    else:
        opt = RSPRFOptimizer(
            geom, hessian_update="bofill", saddle_recovery_max_cycles=0,
            reference_mode=np.r_[np.zeros(9), 1., 0., 0.], **settings,
        )
        assert opt.saddle_imaginary_threshold_cm == 5.
    calc.optimizer = opt
    assert curvature > opt.small_eigval_thresh
    assert opt.hessian_recalc == 500 and opt.hessian_update_window == 1

    convergence = []
    native_check = opt.check_convergence

    def observe_check(*args, **kwargs):
        result = native_check(*args, **kwargs)
        convergence.append(result)
        return result

    opt.check_convergence = observe_check  # observation only, no altered result
    opt.run()
    # Separate analytic endpoint oracle; no calculator call/cadence mutation.
    endpoint_h = calc.hessian_at(geom.cart_coords)
    complete, resolved, _ = true_phva(endpoint_h, geom.cart_coords)
    print("soft-curvature endpoint", dict(
        minimum=minimum, initial=initial, initial_y=y, positive_near=positive_near,
        cycle=opt.cur_cycle, converged=opt.is_converged, stop=opt.stop_reason,
        coordinates=geom.cart_coords[-3:].tolist(), frequencies_cm=complete.tolist(),
        hessian_cycles=[cycle for cycle, _ in calc.hessian_calls],
    ))
    assert opt.is_converged, (opt.cur_cycle, opt.stop_reason, geom.cart_coords[-3:])
    assert opt.cur_cycle < MAX_CYCLES
    assert convergence[-1][0]
    assert all(getattr(convergence[-1][1], key) for key in CRITERIA)
    np.testing.assert_array_equal(geom.cart_coords[:9], coords[:9])
    for point in opt.coords:
        np.testing.assert_array_equal(np.asarray(point)[:9], coords[:9])
    for _, point in calc.hessian_calls:
        np.testing.assert_array_equal(point[:9], coords[:9])
    assert calc.hessian_calls
    assert opt.hessian_recalc == 500  # not a recalc=1 trajectory in disguise
    assert not opt.flatten_enabled

    # Neither opt.H nor an acceptance helper supplies the expected index.
    expected = int(not minimum)
    assert np.count_nonzero(np.diag(endpoint_h)[-3:] < 0.) == expected
    assert np.count_nonzero(complete < 0.) == expected
    assert np.count_nonzero(resolved < 0.) == expected
    fresh_force = calc.get_forces(geom.atoms, geom.cart_coords)["forces"][-3:]
    assert np.max(abs(fresh_force)) <= opt.max_force_thresh
    assert np.sqrt(np.mean(fresh_force**2)) <= opt.rms_force_thresh
    if not positive_near:
        # Rejection alone is not a pass: the surplus well must actually become
        # positively curved, with H requests absent from at least some cycles.
        assert endpoint_h[-2, -2] > 0.
        assert abs(geom.cart_coords[-2] - coords[-2]) > .5
        # A stale model can also reach the well before its first terminal H;
        # do not require a redundant acquisition in that legitimate case.
        assert 1 <= len(calc.hessian_calls) < len(opt.coords)
    else:
        np.testing.assert_array_equal(geom.cart_coords, coords)
    if minimum:
        assert opt._minimum_curvature_valid
    else:
        assert opt._exact_saddle_matches_current_geometry()
        mode = opt.ts_modes[0]
        if isinstance(mode, torch.Tensor):
            mode = mode.detach().cpu().numpy()
        mode = np.asarray(mode).reshape(-1)[-3:]
        assert abs(mode[0]) == pytest.approx(1., abs=1e-12)
        assert abs(geom.cart_coords[-3]) < 1e-8


@pytest.mark.parametrize("initial,y", [("stale", 0.), ("stale", .005), ("calc", 0.)],
                         ids=["stale-zero-gradient", "stale-nonzero-gradient", "exact-H-control"])
def test_soft_hosp_reaches_strict_first_order(tmp_path, soft_curvature, initial, y):
    run_case(tmp_path, soft_curvature, initial=initial, y=y)


def test_positive_near_mode_keeps_legitimate_ts(tmp_path, soft_curvature):
    run_case(tmp_path, soft_curvature, initial="calc", positive_near=True)


def test_soft_negative_minimum_search_reaches_strict_zero(tmp_path, soft_curvature):
    run_case(tmp_path, soft_curvature, minimum=True)
