"""Regression tests for shared rollback and first-order-saddle safeguards."""

from __future__ import annotations

import numpy as np
import pytest
import torch

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.intcoords.exceptions import RebuiltInternalsException
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from pysisyphus.tsoptimizers.TRIM import TRIM
from mlmm.workflows.tsopt import (
    FLATTEN_RETRY_HIGHER_ORDER_CHECKS,
    HessianDimer,
    PATH_MODE_RESTART_AMPLITUDES_ANG,
    _effective_flatten_iterations,
    _flatten_branch_needs_alternate,
    _flatten_once_with_modes_for_geom,
    _mirrored_flatten_start,
    _path_restart_mode_candidates,
    _select_flatten_targets_for_geom,
    _set_cartesian_flatten_coords,
    _transported_path_mode_full,
)


@pytest.mark.parametrize(
    ("configured", "has_reference", "n_imag", "target_negative", "expected"),
    [
        (0, False, 2, None, (0, False)),
        (0, True, 2, True, (0, False)),
        (50, False, 2, None, (50, False)),
        (50, True, 1, False, (50, False)),
        (50, True, 2, True, (50, False)),
        (50, True, 2, False, (0, True)),
        (50, True, 2, None, (0, True)),
    ],
)
def test_flatten_budget_is_explicit_and_path_safe(
    configured: int,
    has_reference: bool,
    n_imag: int,
    target_negative: bool | None,
    expected: tuple[int, bool],
) -> None:
    assert _effective_flatten_iterations(
        configured,
        has_reference_mode=has_reference,
        n_imag=n_imag,
        target_mode_is_negative=target_negative,
    ) == expected


def test_path_mode_restarts_use_two_bounded_displacement_shells() -> None:
    assert PATH_MODE_RESTART_AMPLITUDES_ANG == (-0.10, 0.10, -0.20, 0.20)


def test_flatten_retry_keeps_full_higher_order_persistence_window() -> None:
    assert FLATTEN_RETRY_HIGHER_ORDER_CHECKS == 3


def test_flatten_alternate_start_mirrors_primary_about_saddle() -> None:
    saddle = np.array([1.0, -2.0, 3.0])
    primary = np.array([1.2, -2.3, 3.4])

    alternate = _mirrored_flatten_start(saddle, primary)

    np.testing.assert_allclose(alternate, [0.8, -1.7, 2.6])
    np.testing.assert_allclose(
        0.5 * (primary + alternate),
        saddle,
    )


@pytest.mark.parametrize(
    ("n_imag", "target_negative", "expected"),
    [
        (0, False, True),
        (1, False, True),
        (2, True, True),
        (1, True, False),
    ],
)
def test_flatten_skips_alternate_only_for_requested_first_order_branch(
    n_imag: int,
    target_negative: bool,
    expected: bool,
) -> None:
    class _Optimizer:
        _last_exact_target_mode_is_negative = target_negative

    assert (
        _flatten_branch_needs_alternate(
            {"optimizer": _Optimizer(), "n_imag": n_imag}
        )
        is expected
    )


def test_flatten_restart_carries_full_transported_path_mode() -> None:
    class _Geometry:
        coord_type = "cart"
        cart_coords = np.zeros(6)

    class _Optimizer:
        _last_exact_physical_mode = np.array([3.0, 4.0, 0.0])
        _physical_ts_mode = np.array([1.0, 0.0, 0.0])
        _saddle_recovery_mode = None
        ts_modes = np.empty((0, 3))

        @staticmethod
        def full_from_active(mode):
            full = np.zeros(6)
            full[3:] = mode
            return full

    mode = _transported_path_mode_full(
        _Optimizer(),
        _Geometry(),
        fallback=np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    )

    np.testing.assert_allclose(mode, [0.0, 0.0, 0.0, 0.6, 0.8, 0.0])


def test_path_restart_adds_distinct_initial_soft_root_shell() -> None:
    class _Geometry:
        coord_type = "cart"
        cart_coords = np.zeros(6)

    class _Optimizer:
        _initial_reference_root_mode = np.array([0.0, 1.0, 0.0])

        @staticmethod
        def full_from_active(mode):
            full = np.zeros(6)
            full[3:] = mode
            return full

    modes = _path_restart_mode_candidates(
        _Optimizer(),
        _Geometry(),
        np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    )

    assert [source for source, _ in modes] == [
        "mep-tangent",
        "initial-soft-root",
    ]
    np.testing.assert_allclose(modes[1][1], [0.0, 0.0, 0.0, 0.0, 1.0, 0.0])


class _QuadraticCalculator(Calculator):
    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    @staticmethod
    def _results(coords):
        coords = np.asarray(coords, dtype=float)
        return {
            "energy": float(coords @ coords),
            "forces": -2.0 * coords,
            "hessian": 2.0 * np.eye(len(coords)),
        }

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)


@pytest.mark.parametrize("optimizer_cls", [LBFGS, RFOptimizer])
@pytest.mark.parametrize("tolerance", [-1.0, np.nan, np.inf])
def test_uphill_tolerance_must_be_finite_and_nonnegative(
    tmp_path, optimizer_cls, tolerance
) -> None:
    geom = Geometry(["H"], np.zeros(3), coord_type="cart")
    geom.set_calculator(_QuadraticCalculator(tmp_path))
    kwargs = {"out_dir": tmp_path, "uphill_tolerance": tolerance}
    if optimizer_cls is RFOptimizer:
        kwargs.update(hessian_init="calc", hessian_recalc=1)

    with pytest.raises(
        ValueError,
        match="must be finite and non-negative",
    ):
        optimizer_cls(geom, **kwargs)


class _DoubleWellCalculator(Calculator):
    """One saddle coordinate plus two stiff minimizing coordinates."""

    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    @staticmethod
    def _results(coords):
        x, y, z = np.asarray(coords, dtype=float)
        return {
            "energy": float((x * x - 1.0) ** 2 + 5.0 * y * y + 5.0 * z * z),
            "forces": np.array([-4.0 * x * (x * x - 1.0), -10.0 * y, -10.0 * z]),
            "hessian": np.diag([12.0 * x * x - 4.0, 10.0, 10.0]),
        }

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)


class _TranslationalArtifactCalculator(Calculator):
    """Stationary diatomic with a spurious negative translation only."""

    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    @staticmethod
    def _results(coords):
        coords = np.asarray(coords, dtype=float)
        translation_x = np.array([1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
        translation_x /= np.linalg.norm(translation_x)
        hessian = 2.0 * np.eye(6) - 3.0 * np.outer(
            translation_x, translation_x
        )
        return {
            "energy": 0.0,
            "forces": np.zeros_like(coords),
            "hessian": hessian,
        }

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)


class _SpuriousTranslationDoubleWellCalculator(Calculator):
    """Diatomic double well whose most-negative raw root is always TR noise."""

    well_position = 0.1

    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    @classmethod
    def _results(cls, coords):
        coords = np.asarray(coords, dtype=float)
        stretch = np.array([-1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
        stretch /= np.linalg.norm(stretch)
        translation = np.array([1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
        translation /= np.linalg.norm(translation)
        q = float(stretch @ coords)
        a = cls.well_position
        d_energy_dq = 4.0 * q * (q * q - a * a)
        d2_energy_dq2 = 12.0 * q * q - 4.0 * a * a

        gradient = d_energy_dq * stretch
        # Keep the four transverse directions stiff. The x translation is a
        # deliberately inconsistent negative numerical artifact.
        transverse = np.diag([0.0, 2.0, 2.0, 0.0, 2.0, 2.0])
        hessian = (
            d2_energy_dq2 * np.outer(stretch, stretch)
            - np.outer(translation, translation)
            + transverse
        )
        return {
            "energy": float((q * q - a * a) ** 2),
            "forces": -gradient,
            "hessian": hessian,
        }

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)


class _TorchSpuriousTranslationDoubleWellCalculator(
    _SpuriousTranslationDoubleWellCalculator
):
    def get_hessian(self, atoms, coords, **prepare_kwargs):
        results = self._results(coords)
        results["hessian"] = torch.as_tensor(
            results["hessian"], dtype=torch.float64
        )
        return results


def test_lbfgs_rejects_uphill_trial_through_shared_rollback(tmp_path) -> None:
    geom = Geometry(["H"], np.array([1.0, 0.0, 0.0]), coord_type="cart")
    geom.set_calculator(_QuadraticCalculator(tmp_path))
    opt = LBFGS(
        geom,
        max_step=0.1,
        line_search=False,
        reject_uphill=True,
        out_dir=tmp_path,
    )
    assert opt.reject_uphill is True

    accepted = geom.coords.copy()
    accepted_cart = geom.cart_coords.copy()
    accepted_energy = geom.energy
    accepted_force = geom.forces.copy()
    rejected_step = np.array([1.0, 0.0, 0.0])
    geom.coords = accepted + rejected_step
    opt.coords = [accepted.copy(), geom.coords.copy()]
    opt.cart_coords = [accepted_cart.copy(), geom.cart_coords.copy()]
    opt.energies = [accepted_energy]
    opt.forces = [accepted_force]
    opt.steps = [rejected_step.copy()]
    opt.image_inds = [[0], [0]]
    opt.image_nums = [1, 1]
    opt.cur_cycle = 1

    retry_step = opt.optimize()

    np.testing.assert_allclose(geom.coords, accepted)
    np.testing.assert_allclose(geom.cart_coords, accepted_cart)
    assert len(opt.coords) == len(opt.cart_coords) == 1
    assert len(opt.energies) == len(opt.forces) == 1
    assert opt.steps == []
    assert opt.rejected_uphill_steps == 1
    assert opt._trial_max_step == 0.05
    assert np.abs(retry_step).max() <= 0.05


@pytest.mark.parametrize(
    ("energy_change", "max_step", "expected"),
    [
        (0.0, 3.56e-4, True),
        (2.0e-6, 2.00e-4, True),
        (2.0e-6, 3.00e-4, True),
        (2.0e-6, 3.56e-4, False),
    ],
)
def test_baker_uses_energy_or_max_step_after_max_force(
    tmp_path, energy_change, max_step, expected
) -> None:
    geom = Geometry(["H"], np.zeros(3), coord_type="cart")
    geom.set_calculator(_QuadraticCalculator(tmp_path))
    opt = LBFGS(
        geom,
        thresh="baker",
        overachieve_factor=0.0,
        out_dir=tmp_path,
    )
    opt.cur_cycle = opt.last_cycle
    opt.forces = [np.array([1.0e-5, 0.0, 0.0])]
    opt.energies = [0.0, energy_change]
    opt.steps = [np.array([max_step, 0.0, 0.0])]

    converged, conv_info = opt.check_convergence()

    assert converged is expected
    assert bool(conv_info.max_force_converged)
    assert bool(conv_info.energy_converged) is (energy_change < 1.0e-6)


def test_emergency_stop_accepts_converged_retained_rfo_geometry(
    tmp_path, capsys
) -> None:
    class _EmergencyStopRF(RFOptimizer):
        def optimize(self):
            step = super().optimize()
            if self.cur_cycle >= 0:
                self.uphill_rejection_stalled = True
                self.request_stop(
                    "repeated uphill RFO trials at the emergency trust floor"
                )
            return step

    geom = Geometry(["H"], np.array([1.0e-4, 0.0, 0.0]), coord_type="cart")
    geom.set_calculator(_QuadraticCalculator(tmp_path))
    opt = _EmergencyStopRF(
        geom,
        hessian_init="calc",
        hessian_recalc=1,
        thresh="baker",
        max_cycles=3,
        out_dir=tmp_path,
    )

    opt.run()

    stdout = capsys.readouterr().out
    assert opt.is_converged, stdout
    assert not opt.stopped
    assert opt.termination_status == "converged"
    assert "Converged!" in stdout
    assert "retained lower-energy geometry satisfies" in stdout
    assert "without claiming convergence" not in stdout


def _ts_optimizer(tmp_path, x, **kwargs):
    geom = Geometry(["H"], np.array([x, 0.0, 0.0]), coord_type="cart")
    geom.set_calculator(_DoubleWellCalculator(tmp_path))
    opt_kwargs = dict(
        hessian_init="calc",
        hessian_recalc=1,
        trust_radius=0.1,
        trust_min=1e-5,
        trust_max=0.1,
        thresh="baker",
        check_eigval_structure=True,
        out_dir=tmp_path,
    )
    opt_kwargs.update(kwargs)
    opt = RSPRFOptimizer(geom, **opt_kwargs)
    return geom, opt


def test_exact_saddle_output_follows_the_cycle_row(
    tmp_path, monkeypatch, capsys
) -> None:
    _, opt = _ts_optimizer(tmp_path, 0.0, energy_plateau=False)

    def _refresh(_gradient):
        print("[hessian] Completed FiniteDifference Hessian: 1.00 s")
        return (
            np.zeros(3),
            np.eye(3),
            np.array([-1.0, 1.0, 2.0]),
            np.eye(3),
        )

    def _verify(_eigvals, _eigvecs):
        print("Exact PHVA saddle validation: n_imag=1, lowest=-500.00 cm^-1.")
        return True, np.array([1.0, 0.0, 0.0]), True

    monkeypatch.setattr(opt, "_refresh_exact_saddle_model", _refresh)
    monkeypatch.setattr(opt, "_verify_exact_vibrational_structure", _verify)

    opt._refresh_and_verify_exact_saddle_model(np.zeros(3))
    opt._print_or_defer_cycle_message(
        "Exact saddle validation found no physical imaginary mode."
    )
    assert capsys.readouterr().out == ""
    assert opt.has_deferred_cycle_output()

    opt.cur_cycle = 7
    opt.energies = [0.0, 0.0]
    opt.max_forces = [2.0e-4]
    opt.rms_forces = [1.0e-4]
    opt.max_steps = [3.0e-4]
    opt.rms_steps = [2.0e-4]
    opt.cycle_times = [1.25]

    class _ConvInfo:
        @staticmethod
        def get_convergence():
            return True, True, True, True, True, True

    opt.print_opt_progress(_ConvInfo())

    lines = capsys.readouterr().out.splitlines()
    assert lines[0].split()[0] == "7"
    assert lines[1].startswith("[hessian] Completed")
    assert lines[2].lstrip().startswith("Exact PHVA saddle validation:")
    assert lines[3].lstrip().startswith("Exact saddle validation found")
    assert not opt.has_deferred_cycle_output()


def test_undefined_initial_energy_change_has_no_convergence_mark(
    tmp_path, capsys
) -> None:
    _, opt = _ts_optimizer(tmp_path, 0.0, energy_plateau=False)
    opt.cur_cycle = 0
    opt.energies = [0.0]
    opt.max_forces = [2.0e-4]
    opt.rms_forces = [1.0e-4]
    opt.max_steps = [3.0e-4]
    opt.rms_steps = [2.0e-4]
    opt.cycle_times = [0.0]

    class _ConvInfo:
        @staticmethod
        def get_convergence():
            return True, True, True, True, True, True

    opt.print_opt_progress(_ConvInfo())

    row = capsys.readouterr().out.splitlines()[0]
    assert "nan" in row.lower()
    assert "nan*" not in row.lower()


def test_ts_mode_loss_rejects_trial_and_restores_hessian(tmp_path) -> None:
    geom, opt = _ts_optimizer(tmp_path, 0.4, reject_mode_loss=True)
    opt.prepare_opt()
    accepted = geom.coords.copy()
    accepted_cart = geom.cart_coords.copy()
    accepted_energy = geom.energy
    accepted_force = geom.forces.copy()
    accepted_hessian = opt.H.copy()

    # x=0.4 has negative curvature; x=0.8 is already inside the convex basin.
    trial_step = np.array([0.4, 0.0, 0.0])
    geom.coords = accepted + trial_step
    opt.coords = [accepted.copy(), geom.coords.copy()]
    opt.cart_coords = [accepted_cart.copy(), geom.cart_coords.copy()]
    opt.energies = [accepted_energy]
    opt.forces = [accepted_force]
    opt.steps = [trial_step.copy()]
    opt.predicted_energy_changes = [-0.1]
    opt.image_inds = [[0], [0]]
    opt.image_nums = [1, 1]
    opt.cur_cycle = 1
    opt.hessian_recalc_in = 0

    _, _, _, eigvals, _, resetted = opt.housekeeping()

    np.testing.assert_allclose(geom.coords, accepted)
    np.testing.assert_allclose(opt.H, accepted_hessian)
    assert eigvals[0] < 0.0
    assert resetted is True
    assert opt.rejected_mode_loss_steps == 1
    assert opt.trust_radius == 0.025
    assert opt.steps == []
    assert opt.predicted_energy_changes == []


def test_energy_plateau_cannot_bypass_ts_curvature_requirement(tmp_path) -> None:
    geom, opt = _ts_optimizer(tmp_path, 1.0, energy_plateau_window=5)
    opt.cur_cycle = 0
    opt.last_cycle = 0
    opt.forces = [np.zeros(3)]
    opt.energies = [0.0] * 5
    opt.ts_mode_eigvals = np.array([1.0])

    converged, conv_info = opt.check_convergence(np.zeros(3))

    assert converged is False
    assert conv_info.desired_eigval_structure is False


def test_verified_phva_mode_rechecks_when_force_and_energy_are_terminal(tmp_path) -> None:
    _, opt = _ts_optimizer(tmp_path, 0.0, energy_plateau=False)
    opt._physical_ts_mode = np.array([1.0, 0.0, 0.0])
    opt._last_exact_n_imaginary = 1
    opt.forces = [np.zeros(3)]
    opt.energies = [0.0, 0.0]
    opt.steps = [np.array([0.01, 0.0, 0.0])]

    # The outgoing step used by check_convergence() does not exist yet in
    # housekeeping().  A nonterminal incoming step must not suppress the exact
    # check once the current force and energy are terminal.
    assert opt._near_terminal_without_eigval_check() is True


def test_stale_exact_saddle_cannot_authorize_rolled_back_minimum(tmp_path) -> None:
    geom, opt = _ts_optimizer(tmp_path, 1.0, energy_plateau=False)
    opt.cur_cycle = 3
    opt.last_cycle = 0
    opt.forces = [np.zeros(3), np.zeros(3)]
    opt.energies = [0.0, 0.0]
    opt.steps = [np.zeros(3)]
    opt.ts_mode_eigvals = np.array([-1.0])
    opt._last_exact_n_imaginary = 1
    opt._last_exact_saddle_verified = True
    # Simulate exact verification of a saddle trial followed by rollback to
    # the current x=1 local minimum.
    opt._last_exact_cart_coords = np.array([0.0, 0.0, 0.0])

    converged, _ = opt.check_convergence()

    assert converged is False
    assert not opt._exact_saddle_matches_current_geometry()


def test_current_coordinate_exact_saddle_can_authorize_convergence(tmp_path) -> None:
    geom, opt = _ts_optimizer(tmp_path, 0.0, energy_plateau=False)
    opt.cur_cycle = 3
    opt.last_cycle = 0
    opt.forces = [np.zeros(3), np.zeros(3)]
    opt.energies = [0.0, 0.0]
    opt.steps = [np.zeros(3)]
    opt.ts_mode_eigvals = np.array([-4.0])
    opt._last_exact_n_imaginary = 1
    opt._last_exact_saddle_verified = True
    opt._last_exact_cart_coords = geom.cart_coords.copy()

    assert opt._near_terminal_without_eigval_check() is True
    converged, conv_info = opt.check_convergence()

    assert converged is True
    assert conv_info.energy_converged is True


def test_exact_non_baker_saddle_ignores_outgoing_step(tmp_path) -> None:
    geom, opt = _ts_optimizer(
        tmp_path, 0.0, thresh="gau_loose", energy_plateau=False
    )
    opt.cur_cycle = 3
    opt.last_cycle = 0
    opt.forces = [np.full(3, 9.0e-4)]
    opt.energies = [0.0, 1.0e-3]
    opt.steps = [np.full(3, 4.0e-2)]
    opt.ts_mode_eigvals = np.array([-4.0])
    opt._last_exact_n_imaginary = 1
    opt._last_exact_saddle_verified = True
    opt._last_exact_cart_coords = geom.cart_coords.copy()

    converged, conv_info = opt.check_convergence()

    assert bool(conv_info.max_force_converged)
    assert bool(conv_info.rms_force_converged)
    assert not bool(conv_info.max_step_converged)
    assert not bool(conv_info.rms_step_converged)
    assert converged is True


def test_exact_saddle_does_not_override_never_threshold(tmp_path) -> None:
    geom, opt = _ts_optimizer(
        tmp_path, 0.0, thresh="never", energy_plateau=False
    )
    opt.cur_cycle = 3
    opt.last_cycle = 0
    opt.forces = [np.zeros(3)]
    opt.energies = [0.0, 0.0]
    opt.steps = [np.zeros(3)]
    opt.ts_mode_eigvals = np.array([-4.0])
    opt._last_exact_n_imaginary = 1
    opt._last_exact_saddle_verified = True
    opt._last_exact_cart_coords = geom.cart_coords.copy()

    converged, _ = opt.check_convergence()

    assert converged is False


def test_energy_plateau_with_exact_saddle_but_failed_force_stalls(
    tmp_path,
) -> None:
    # Required curvature is present (a verified first-order saddle at these
    # exact coordinates) and the energy is flat over the window, but the
    # configured current force/step criteria fail: this is a stall,
    # never a converged saddle.
    geom, opt = _ts_optimizer(
        tmp_path, 0.0, energy_plateau=True, energy_plateau_window=3
    )
    opt.cur_cycle = 3
    opt.last_cycle = 0
    opt.forces = [np.full(3, 1.0)]      # large current force -> not converged
    opt.energies = [0.0, 0.0, 0.0]      # flat window -> energy plateau
    opt.steps = [np.full(3, 1.0)]
    opt.ts_mode_eigvals = np.array([-4.0])
    opt._last_exact_n_imaginary = 1
    opt._last_exact_saddle_verified = True
    opt._last_exact_cart_coords = geom.cart_coords.copy()

    converged, conv_info = opt.check_convergence()

    assert converged is False
    assert opt.is_stalled is True
    assert opt.termination_status == "stalled"
    assert "energy plateau" in opt.stop_reason
    assert not bool(conv_info.max_force_converged)


def test_energy_plateau_provisional_check_does_not_stall_ts_optimizer(
    tmp_path,
) -> None:
    # A provisional probe (allow_stall=False) must never arm the stall, even
    # with an exact saddle and a full flat window.
    geom, opt = _ts_optimizer(
        tmp_path, 0.0, energy_plateau=True, energy_plateau_window=3
    )
    opt.cur_cycle = 3
    opt.last_cycle = 0
    opt.forces = [np.full(3, 1.0)]
    opt.energies = [0.0, 0.0, 0.0]
    opt.steps = [np.full(3, 1.0)]
    opt.ts_mode_eigvals = np.array([-4.0])
    opt._last_exact_n_imaginary = 1
    opt._last_exact_saddle_verified = True
    opt._last_exact_cart_coords = geom.cart_coords.copy()

    converged, _ = opt.check_convergence(allow_stall=False)

    assert converged is False
    assert opt.is_stalled is False


def test_exact_current_saddle_still_requires_baker_energy_or_step(
    tmp_path
) -> None:
    geom, opt = _ts_optimizer(tmp_path, 0.0, energy_plateau=False)
    opt.cur_cycle = 3
    opt.last_cycle = 0
    opt.forces = [np.zeros(3), np.zeros(3)]
    opt.energies = [0.0, 1.0e-3]
    opt.steps = [np.full(3, 0.1)]
    opt.ts_mode_eigvals = np.array([1.0])
    opt._last_exact_n_imaginary = 1
    opt._last_exact_saddle_verified = True
    opt._last_exact_cart_coords = geom.cart_coords.copy()

    converged, conv_info = opt.check_convergence()

    assert not bool(conv_info.max_step_converged)
    assert conv_info.desired_eigval_structure is False
    assert converged is False


def test_ts_baker_keeps_rms_force_diagnostic_only(tmp_path) -> None:
    geom, opt = _ts_optimizer(tmp_path, 0.0, energy_plateau=False)
    opt.cur_cycle = 3
    opt.last_cycle = 0
    opt.forces = [np.full(3, 2.5e-4), np.full(3, 2.5e-4)]
    opt.energies = [0.0, 0.0]
    opt.steps = [np.zeros(3)]
    opt._last_exact_n_imaginary = 1
    opt._last_exact_saddle_verified = True
    opt._last_exact_cart_coords = geom.cart_coords.copy()

    assert opt._near_terminal_without_eigval_check() is True
    converged, conv_info = opt.check_convergence()

    assert bool(conv_info.max_force_converged)
    assert not bool(conv_info.rms_force_converged)
    assert converged is True


def test_exact_higher_order_saddle_cannot_authorize_first_order_convergence(
    tmp_path
) -> None:
    geom, opt = _ts_optimizer(tmp_path, 0.0, energy_plateau=False)
    opt.cur_cycle = 3
    opt.last_cycle = 0
    opt.forces = [np.zeros(3), np.zeros(3)]
    opt.energies = [0.0, 0.0]
    opt.steps = [np.zeros(3)]
    opt.ts_mode_eigvals = np.array([-4.0])
    opt._last_exact_n_imaginary = 2
    opt._last_exact_saddle_verified = False
    opt._last_exact_cart_coords = geom.cart_coords.copy()

    converged, _ = opt.check_convergence()

    assert converged is False
    assert not opt._exact_saddle_matches_current_geometry()


def test_exact_verifier_retains_curvature_but_rejects_higher_order_status(
    tmp_path, monkeypatch
) -> None:
    _, opt = _ts_optimizer(tmp_path, 0.0)
    opt.forces = [np.zeros(3)]
    opt.cur_H = np.diag([-2.0, -1.0, 3.0])
    monkeypatch.setattr(
        opt,
        "_mw_frequencies_and_modes",
        lambda: (np.array([-100.0, -20.0, 30.0]), torch.eye(3)),
    )
    monkeypatch.setattr(
        opt,
        "_recovery_mode_from_mw",
        lambda modes, mode_index: np.array([1.0, 0.0, 0.0]),
    )

    has_negative, _, _ = opt._verify_exact_vibrational_structure(
        np.array([-2.0, -1.0, 3.0]), np.eye(3)
    )

    assert has_negative is True
    assert opt._last_exact_n_imaginary == 2
    assert opt._last_exact_saddle_verified is False
    assert opt._best_exact_saddle is None


def test_exact_noncart_verifier_uses_first_order_saddle_not_selected_root(
    tmp_path, monkeypatch
) -> None:
    _, opt = _ts_optimizer(tmp_path, 0.0)
    opt.roots = [1]
    opt.forces = [np.zeros(3)]
    opt.cur_H = np.diag([-2.0, 1.0, 3.0])
    monkeypatch.setattr(opt, "_mw_frequencies_and_modes", lambda: None)

    has_negative, _, _ = opt._verify_exact_vibrational_structure(
        np.array([-2.0, 1.0, 3.0]), np.eye(3)
    )

    assert has_negative is True
    assert opt._last_exact_n_imaginary == 1
    assert opt._last_exact_saddle_verified is True


@pytest.mark.parametrize(
    "device",
    ["cpu"] + (["cuda"] if torch.cuda.is_available() else []),
)
def test_hessian_reference_accepts_torch_eigenvectors(
    tmp_path, monkeypatch, device
) -> None:
    _, opt = _ts_optimizer(tmp_path, 0.0)
    opt.hessian_ref = np.diag([-3.0, 1.0, 2.0])

    def seed_torch_hessian(self, *_args, **_kwargs):
        self.H = torch.diag(
            torch.tensor([-2.0, 1.0, 3.0], device=device)
        )

    monkeypatch.setattr(HessianOptimizer, "prepare_opt", seed_torch_hessian)

    opt.prepare_opt()

    assert opt.roots.tolist() == [0]
    assert isinstance(opt.ts_modes, torch.Tensor)


def test_exact_higher_order_saddle_keeps_path_correlated_negative_mode(
    tmp_path, monkeypatch
) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        0.0,
        reference_mode=np.array([0.0, 1.0, 0.0]),
    )
    # Continuous single-root transport has drifted to the other member of the
    # two-dimensional negative subspace.  The immutable MEP tangent must decide
    # which mode survives a flatten restart.
    opt.ts_modes = np.array([[1.0, 0.0, 0.0]])
    opt.forces = [np.zeros(3)]
    opt.cur_H = np.diag([-2.0, -1.0, 3.0])
    monkeypatch.setattr(
        opt,
        "_mw_frequencies_and_modes",
        lambda: (np.array([-100.0, -20.0, 30.0]), torch.eye(3)),
    )
    exact_modes = {
        0: np.array([1.0, 0.0, 0.0]),
        1: np.array([0.0, 1.0, 0.0]),
    }
    monkeypatch.setattr(
        opt,
        "_recovery_mode_from_mw",
        lambda modes, mode_index: exact_modes.get(mode_index),
    )

    has_negative, physical_mode, _ = opt._verify_exact_vibrational_structure(
        np.array([-2.0, -1.0, 3.0]), np.eye(3)
    )

    assert has_negative is True
    np.testing.assert_allclose(physical_mode, [0.0, 1.0, 0.0])
    assert opt._last_exact_target_mode_index == 1
    assert opt._last_exact_target_mode_overlap == pytest.approx(1.0)
    assert opt._last_exact_target_mode_reanchored is True
    np.testing.assert_allclose(
        opt._last_exact_physical_mode, [0.0, 1.0, 0.0]
    )


def test_reference_mismatch_is_diagnostic_for_exact_first_order_saddle(
    tmp_path, monkeypatch
) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        0.0,
        reference_mode=np.array([0.0, 1.0, 0.0]),
    )
    opt.forces = [np.zeros(3)]
    opt.cur_H = np.diag([-2.0, 1.0, 3.0])
    monkeypatch.setattr(
        opt,
        "_mw_frequencies_and_modes",
        lambda: (np.array([-100.0, 20.0, 30.0]), torch.eye(3)),
    )
    monkeypatch.setattr(
        opt,
        "_recovery_mode_from_mw",
        lambda modes, mode_index: np.eye(3)[mode_index],
    )

    has_saddle, recovery_mode, _ = opt._verify_exact_vibrational_structure(
        np.array([-2.0, 1.0, 3.0]), np.eye(3)
    )

    assert has_saddle is True
    assert opt._last_exact_n_imaginary == 1
    assert opt._last_exact_saddle_verified is True
    assert opt._last_exact_target_mode_index == 1
    assert opt._last_exact_target_mode_overlap == pytest.approx(1.0)
    assert opt._last_exact_target_mode_is_negative is False
    np.testing.assert_allclose(recovery_mode, [0.0, 1.0, 0.0])


def test_first_path_recovery_keeps_complete_multimode_tangent(
    tmp_path, monkeypatch
) -> None:
    reference = np.array([1.0, 1.0, 0.0])
    reference /= np.linalg.norm(reference)
    _, opt = _ts_optimizer(
        tmp_path,
        0.0,
        reference_mode=reference,
    )
    opt.ts_modes = np.array([[1.0, 0.0, 0.0]])
    opt.negative_mode_seen = False
    opt.forces = [np.zeros(3)]
    opt.cur_H = np.diag([1.0, 2.0, 3.0])
    monkeypatch.setattr(
        opt,
        "_mw_frequencies_and_modes",
        lambda: (np.array([10.0, 20.0, 30.0]), torch.eye(3)),
    )
    monkeypatch.setattr(
        opt,
        "_recovery_mode_from_mw",
        lambda modes, mode_index: np.eye(3)[mode_index],
    )

    has_saddle, recovery_mode, _ = opt._verify_exact_vibrational_structure(
        np.array([1.0, 2.0, 3.0]), np.eye(3)
    )

    assert has_saddle is False
    np.testing.assert_allclose(recovery_mode, reference)


def test_recovery_uses_transported_mode_after_target_was_negative(
    tmp_path, monkeypatch
) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        0.0,
        reference_mode=np.array([0.0, 1.0, 0.0]),
    )
    opt.ts_modes = np.array([[1.0, 0.0, 0.0]])
    opt.negative_mode_seen = True
    opt.forces = [np.zeros(3)]
    opt.cur_H = np.diag([1.0, 2.0, 3.0])
    monkeypatch.setattr(
        opt,
        "_mw_frequencies_and_modes",
        lambda: (np.array([10.0, 20.0, 30.0]), torch.eye(3)),
    )
    monkeypatch.setattr(
        opt,
        "_recovery_mode_from_mw",
        lambda modes, mode_index: np.eye(3)[mode_index],
    )

    has_saddle, recovery_mode, _ = opt._verify_exact_vibrational_structure(
        np.array([1.0, 2.0, 3.0]), np.eye(3)
    )

    assert has_saddle is False
    np.testing.assert_allclose(recovery_mode, [1.0, 0.0, 0.0])


def test_higher_order_saddle_is_rejected_independent_of_reference_mode(
    tmp_path, monkeypatch
) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        0.0,
        reference_mode=np.array([0.0, 0.0, 1.0]),
        max_higher_order_checks=1,
    )
    opt.forces = [np.zeros(3)]
    opt.cur_H = np.diag([-2.0, -1.0, 3.0])
    monkeypatch.setattr(
        opt,
        "_mw_frequencies_and_modes",
        lambda: (np.array([-100.0, -20.0, 30.0]), torch.eye(3)),
    )
    monkeypatch.setattr(
        opt,
        "_recovery_mode_from_mw",
        lambda modes, mode_index: np.eye(3)[mode_index],
    )

    has_saddle, recovery_mode, _ = opt._verify_exact_vibrational_structure(
        np.array([-2.0, -1.0, 3.0]), np.eye(3)
    )

    assert has_saddle is True
    assert opt._last_exact_n_imaginary == 2
    assert opt._last_exact_saddle_verified is False
    assert opt._last_exact_target_mode_is_negative is False
    assert opt.higher_order_saddle_checks == 1
    assert opt.stop_requested is True
    np.testing.assert_allclose(recovery_mode, [0.0, 0.0, 1.0])


def test_single_imaginary_path_mode_is_exact_first_order_saddle(
    tmp_path, monkeypatch
) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        0.0,
        reference_mode=np.array([1.0, 0.0, 0.0]),
    )
    opt.forces = [np.zeros(3)]
    opt.cur_H = np.diag([-2.0, 1.0, 3.0])
    monkeypatch.setattr(
        opt,
        "_mw_frequencies_and_modes",
        lambda: (np.array([-100.0, 20.0, 30.0]), torch.eye(3)),
    )
    monkeypatch.setattr(
        opt,
        "_recovery_mode_from_mw",
        lambda modes, mode_index: np.eye(3)[mode_index],
    )

    has_saddle, physical_mode, _ = opt._verify_exact_vibrational_structure(
        np.array([-2.0, 1.0, 3.0]), np.eye(3)
    )

    assert has_saddle is True
    assert opt._last_exact_saddle_verified is True
    assert opt._last_exact_target_mode_index == 0
    assert opt._last_exact_target_mode_is_negative is True
    np.testing.assert_allclose(physical_mode, [1.0, 0.0, 0.0])
    np.testing.assert_allclose(
        opt._last_exact_physical_mode, [1.0, 0.0, 0.0]
    )


def test_exact_identity_uses_overlap_transported_mode_after_rotation(
    tmp_path, monkeypatch
) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        0.0,
        reference_mode=np.array([0.0, 1.0, 0.0]),
    )
    # The initial path tangent was y, but continuous overlap tracking rotated
    # the same reaction mode to x after exact PHVA had first verified it.
    opt.ts_modes = np.array([[1.0, 0.0, 0.0]])
    opt.negative_mode_seen = True
    opt.forces = [np.zeros(3)]
    opt.cur_H = np.diag([-2.0, 1.0, 3.0])
    monkeypatch.setattr(
        opt,
        "_mw_frequencies_and_modes",
        lambda: (np.array([-100.0, 20.0, 30.0]), torch.eye(3)),
    )
    monkeypatch.setattr(
        opt,
        "_recovery_mode_from_mw",
        lambda modes, mode_index: np.eye(3)[mode_index],
    )

    has_saddle, physical_mode, _ = opt._verify_exact_vibrational_structure(
        np.array([-2.0, 1.0, 3.0]), np.eye(3)
    )

    assert has_saddle is True
    assert opt._last_exact_target_mode_index == 0
    assert opt._last_exact_target_mode_is_negative is True
    np.testing.assert_allclose(physical_mode, [1.0, 0.0, 0.0])


def test_exact_identity_keeps_full_path_until_first_physical_crossing(
    tmp_path, monkeypatch
) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        0.0,
        reference_mode=np.array([0.0, 1.0, 0.0]),
    )
    # The numerical uphill root drifted to x while the requested path mode has
    # never been physically negative. Exact validation must still identify y.
    opt.ts_modes = np.array([[1.0, 0.0, 0.0]])
    opt.negative_mode_seen = False
    opt.forces = [np.zeros(3)]
    opt.cur_H = np.diag([-2.0, 1.0, 3.0])
    monkeypatch.setattr(
        opt,
        "_mw_frequencies_and_modes",
        lambda: (np.array([-100.0, 20.0, 30.0]), torch.eye(3)),
    )
    monkeypatch.setattr(
        opt,
        "_recovery_mode_from_mw",
        lambda modes, mode_index: np.eye(3)[mode_index],
    )

    has_saddle, recovery_mode, _ = opt._verify_exact_vibrational_structure(
        np.array([-2.0, 1.0, 3.0]), np.eye(3)
    )

    assert has_saddle is True
    assert opt._last_exact_saddle_verified is True
    assert opt._last_exact_target_mode_index == 1
    assert opt._last_exact_target_mode_is_negative is False
    np.testing.assert_allclose(recovery_mode, [0.0, 1.0, 0.0])


def test_path_mode_overlap_tracking_can_follow_mode_into_positive_spectrum(
    tmp_path
) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        0.0,
        reference_mode=np.array([0.0, 1.0, 0.0]),
    )
    opt.ts_modes = np.array([[0.0, 1.0, 0.0]])
    opt.ts_mode_eigvals = np.array([-1.0])

    opt.update_ts_mode(
        np.array([-2.0, 1.0, 3.0]),
        np.eye(3),
    )

    assert int(opt.roots[0]) == 1
    assert float(opt.ts_mode_eigvals[0]) == pytest.approx(1.0)
    np.testing.assert_allclose(opt.ts_modes[0], [0.0, 1.0, 0.0])


def test_raw_root_tracking_preserves_separate_exact_phva_snapshot(tmp_path) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        0.0,
        reference_mode=np.array([0.0, 1.0, 0.0]),
    )
    physical_mode = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    opt._physical_ts_mode = physical_mode.copy()
    opt._last_exact_physical_mode = physical_mode.copy()

    opt.update_ts_mode(
        np.array([-2.0, 1.0, 3.0]),
        np.eye(3),
    )

    np.testing.assert_allclose(opt._physical_ts_mode, [1.0, 0.0, 0.0])
    np.testing.assert_allclose(opt._last_exact_physical_mode, physical_mode)


def test_initial_reference_root_prefers_softer_near_tied_mode(tmp_path) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        0.0,
        reference_mode=np.array([1.0, 1.0, 0.0]),
    )

    selected = opt._select_initial_reference_root(
        np.array([0.01, 0.50, 1.00]),
        np.array([0.69, 0.72, 0.10]),
    )

    assert selected == 0


def test_initial_reference_root_keeps_clearly_dominant_overlap(tmp_path) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        0.0,
        reference_mode=np.array([1.0, 1.0, 0.0]),
    )

    selected = opt._select_initial_reference_root(
        np.array([0.01, 0.50, 1.00]),
        np.array([0.30, 0.80, 0.10]),
    )

    assert selected == 1


def test_unrelated_initial_negative_root_does_not_mark_path_mode_seen(
    tmp_path,
) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        0.0,
        reference_mode=np.array([0.0, 1.0, 0.0]),
    )

    opt.prepare_opt()

    assert int(opt.roots[0]) == 1
    assert opt._initial_reference_root_index == 1
    assert opt._initial_reference_root_overlap == pytest.approx(1.0)
    assert np.isfinite(opt._initial_reference_root_eigenvalue)
    assert opt.ts_mode_eigvals[0] > 0.0
    assert opt.negative_mode_seen is False
    np.testing.assert_allclose(
        opt._initial_reference_root_mode, [0.0, 1.0, 0.0]
    )


def test_path_initial_raw_negative_waits_for_exact_phva(tmp_path) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        0.0,
        reference_mode=np.array([1.0, 0.0, 0.0]),
    )

    opt.prepare_opt()

    assert int(opt.roots[0]) == 0
    assert opt.ts_mode_eigvals[0] < 0.0
    assert opt.negative_mode_seen is False


def test_path_quasi_negative_does_not_arm_mode_loss_before_exact_check(
    tmp_path,
) -> None:
    geom, opt = _ts_optimizer(
        tmp_path,
        1.0,
        reference_mode=np.array([1.0, 0.0, 0.0]),
        hessian_recalc=None,
        energy_plateau=False,
    )
    opt.prepare_opt()
    assert opt.negative_mode_seen is False

    # Simulate a transient negative quasi-Newton root away from terminal force
    # criteria. It may guide the current step but must not latch rollback for
    # all subsequent trials until exact PHVA or recovery verifies the target.
    opt.H = np.diag([-1.0, 10.0, 10.0])
    opt.coords = [geom.coords.copy()]
    opt.cart_coords = [geom.cart_coords.copy()]
    opt.energies = [geom.energy]
    opt.forces = [np.array([1.0, 0.0, 0.0])]
    opt.image_inds = [[0]]
    opt.image_nums = [1]
    opt.cur_cycle = 1

    opt.housekeeping()

    assert opt.negative_mode_seen is False


def test_cartesian_reference_is_transformed_to_internal_hessian_space(
    tmp_path
) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        0.0,
        reference_mode=np.array([1.0, 2.0, 3.0]),
    )

    class _Internal:
        B = np.array(
            [
                [1.0, 0.0, 1.0],
                [0.0, 2.0, 0.0],
            ]
        )

    opt.geometry.internal = _Internal()
    mapped = opt._reference_mode_for_eigenspace(np.eye(2))

    expected = np.array([4.0, 4.0])
    expected /= np.linalg.norm(expected)
    np.testing.assert_allclose(mapped, expected)


def test_internal_exact_check_treats_reference_mismatch_as_diagnostic(
    tmp_path, monkeypatch
) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        0.0,
        reference_mode=np.array([0.0, 1.0, 0.0]),
    )
    opt.forces = [np.zeros(3)]
    monkeypatch.setattr(opt, "_mw_frequencies_and_modes", lambda: None)

    has_saddle, recovery_mode, physical_checked = (
        opt._verify_exact_vibrational_structure(
            np.array([-2.0, 1.0, 3.0]), np.eye(3)
        )
    )

    assert has_saddle is True
    assert physical_checked is True
    assert opt._last_exact_n_imaginary == 1
    assert opt._last_exact_saddle_verified is True
    assert opt._last_exact_target_mode_index == 1
    np.testing.assert_allclose(recovery_mode, [0.0, 1.0, 0.0])


def test_persistent_higher_order_saddle_requests_explicit_flatten_restart(
    tmp_path, monkeypatch
) -> None:
    _, opt = _ts_optimizer(tmp_path, 0.0, max_higher_order_checks=3)
    opt.forces = [np.zeros(3)]
    opt.cur_H = np.diag([-2.0, -1.0, 3.0])
    monkeypatch.setattr(
        opt,
        "_mw_frequencies_and_modes",
        lambda: (np.array([-100.0, -20.0, 30.0]), torch.eye(3)),
    )
    monkeypatch.setattr(
        opt,
        "_recovery_mode_from_mw",
        lambda modes, mode_index: np.eye(3)[mode_index],
    )

    for _ in range(3):
        opt._verify_exact_vibrational_structure(
            np.array([-2.0, -1.0, 3.0]), np.eye(3)
        )

    assert opt.stop_requested is True
    assert "higher-order saddle" in opt.stop_reason


def test_flatten_target_can_preserve_path_correlated_nonlowest_mode() -> None:
    freqs = np.array([-100.0, -20.0, 30.0, 40.0, 50.0, 60.0])
    modes = torch.eye(6, dtype=torch.float64)

    targets = _select_flatten_targets_for_geom(
        freqs,
        modes,
        np.zeros(6),
        neg_freq_thresh_cm=5.0,
        root=0,
        flatten_sep_cutoff=0.0,
        flatten_k=1,
        primary_idx=1,
    )

    assert targets == [0]


def test_flatten_displacement_preserves_unweighted_normal_mode(
    monkeypatch,
) -> None:
    class _Geometry:
        def __init__(self) -> None:
            self.cart_coords = np.zeros(6, dtype=float)

        @property
        def coords(self):
            return self.cart_coords

        @coords.setter
        def coords(self, value) -> None:
            self.cart_coords = np.asarray(value, dtype=float).copy()

    geom = _Geometry()
    masses = np.array([12.011, 1.008])
    modes = torch.zeros((2, 6), dtype=torch.float64)
    # Keep mode 0 as the path-correlated saddle direction. Mode 1 mixes C and
    # H so an erroneous second mass scaling would visibly rotate it.
    modes[0, 0] = 1.0
    modes[1, 1] = 1.0
    modes[1, 3] = 1.0
    reference = np.zeros(6)
    reference[0] = 1.0

    energies = iter((0.0, -1.0, -0.5))
    monkeypatch.setattr(
        "mlmm.workflows.tsopt._calc_energy",
        lambda geometry, kwargs: next(energies),
    )

    flattened = _flatten_once_with_modes_for_geom(
        geom,
        masses,
        {},
        np.array([-100.0, -40.0]),
        modes,
        neg_freq_thresh_cm=5.0,
        flatten_amp_ang=0.10,
        flatten_sep_cutoff=0.0,
        flatten_k=1,
        root=0,
        reference_mode=reference,
    )

    expected = modes[1].numpy() / np.sqrt(np.repeat(masses, 3))
    expected /= np.linalg.norm(expected)
    actual = geom.cart_coords / np.linalg.norm(geom.cart_coords)
    assert flattened is True
    np.testing.assert_allclose(actual, expected, atol=1.0e-12)


@pytest.mark.parametrize("implementation", ["module", "class"])
def test_ts_flatten_uses_cartesian_boundary_for_redundant_geometry(
    implementation: str,
    monkeypatch,
) -> None:
    atoms = ["C", "C", "O", "C", "N", "H"]
    cart_coords = np.array(
        [
            0.0, 0.0, 0.0,
            1.4, 0.0, 0.0,
            2.1, 1.2, 0.0,
            3.5, 1.2, 0.1,
            4.2, 0.1, 0.2,
            5.6, 0.1, 0.0,
        ],
        dtype=float,
    )
    geom = Geometry(atoms, cart_coords, coord_type="redund")
    assert geom.coords.size != geom.cart_coords.size

    masses = np.array([12.011, 12.011, 15.999, 12.011, 14.007, 1.008])
    modes = torch.zeros((2, cart_coords.size), dtype=torch.float64)
    modes[0, 0] = 1.0
    modes[1, 4] = 1.0
    probes = []
    energies = iter((0.0, -1.0, -2.0))

    class _CartesianProbeCalculator:
        def get_energy(self, elements, coords):
            probe = np.asarray(coords, dtype=float).reshape(-1).copy()
            assert probe.size == cart_coords.size
            probes.append(probe)
            return {"energy": next(energies)}

    calculator = _CartesianProbeCalculator()
    monkeypatch.setattr(
        "mlmm.workflows.tsopt._calc_energy",
        lambda geometry, kwargs: calculator.get_energy(
            geometry.atoms, geometry.cart_coords
        )["energy"],
    )

    if implementation == "module":
        flattened = _flatten_once_with_modes_for_geom(
            geom,
            masses,
            {},
            np.array([-100.0, -40.0]),
            modes,
            neg_freq_thresh_cm=5.0,
            flatten_amp_ang=0.01,
            flatten_sep_cutoff=0.0,
            flatten_k=1,
            root=0,
            reference_mode=modes[0].numpy(),
        )
    else:
        runner = object.__new__(HessianDimer)
        runner.geom = geom
        runner.masses_amu = masses
        runner.neg_freq_thresh_cm = 5.0
        runner.flatten_amp_ang = 0.01
        runner.flatten_sep_cutoff = 0.0
        runner.root = 0
        runner.calc_kwargs = {}
        flattened = HessianDimer._flatten_once_with_modes(
            runner, np.array([-100.0, -40.0]), modes
        )

    assert flattened is True
    assert len(probes) == 3
    assert all(probe.shape == (cart_coords.size,) for probe in probes)
    np.testing.assert_allclose(geom.cart_coords, probes[-1], atol=1.0e-12)
    np.testing.assert_allclose(
        geom.internal.coords3d.reshape(-1), geom.cart_coords, atol=1.0e-12
    )
    np.testing.assert_allclose(geom.coords, geom.internal.coords, atol=1.0e-12)


def test_cartesian_flatten_setter_accepts_only_completed_rebuild_signal() -> None:
    class _RebuiltGeometry:
        def __init__(self) -> None:
            self.installed = None
            self.clear_count = 0

        @property
        def cart_coords(self):
            return self.installed

        @cart_coords.setter
        def cart_coords(self, value) -> None:
            self.installed = np.asarray(value, dtype=float).copy()
            raise RebuiltInternalsException()

        def clear(self) -> None:
            self.clear_count += 1

    rebuilt = _RebuiltGeometry()
    _set_cartesian_flatten_coords(rebuilt, np.arange(6.0).reshape(2, 3))
    np.testing.assert_allclose(rebuilt.installed, np.arange(6.0))
    assert rebuilt.clear_count == 1

    class _InvalidGeometry(_RebuiltGeometry):
        @property
        def cart_coords(self):
            return self.installed

        @cart_coords.setter
        def cart_coords(self, value) -> None:
            raise ValueError("invalid Cartesian trial")

    with pytest.raises(ValueError, match="invalid Cartesian trial"):
        _set_cartesian_flatten_coords(_InvalidGeometry(), np.zeros(6))


def test_path_reference_mode_overrides_unrelated_lowest_mode_for_recovery(
    tmp_path, monkeypatch
) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        1.0,
        reference_mode=np.array([0.0, 1.0, 0.0]),
    )
    opt.forces = [np.zeros(3)]
    opt.cur_H = np.diag([1.0, 2.0, 3.0])
    monkeypatch.setattr(
        opt,
        "_mw_frequencies_and_modes",
        lambda: (np.array([10.0, 20.0, 30.0]), torch.eye(3)),
    )
    monkeypatch.setattr(
        opt,
        "_recovery_mode_from_mw",
        lambda modes, mode_index: np.array([1.0, 0.0, 0.0]),
    )

    has_saddle, recovery_mode, _ = opt._verify_exact_vibrational_structure(
        np.array([1.0, 2.0, 3.0]), np.eye(3)
    )

    assert has_saddle is False
    np.testing.assert_allclose(recovery_mode, [0.0, 1.0, 0.0])
    assert opt.track_mode_by_overlap is True


def test_guarded_failure_restores_best_exact_saddle_candidate(tmp_path) -> None:
    geom, opt = _ts_optimizer(tmp_path, 0.0)
    opt.forces = [np.array([1.0e-4, 0.0, 0.0])]
    opt.cur_cycle = 7
    saddle_coords = geom.cart_coords.copy()
    opt._record_exact_saddle_candidate()

    geom.cart_coords = np.array([1.0, 0.0, 0.0])
    assert opt._restore_best_exact_saddle() is True

    np.testing.assert_allclose(geom.cart_coords, saddle_coords)


def test_exact_hessian_recovery_displaces_back_along_incoming_mode(tmp_path) -> None:
    geom, opt = _ts_optimizer(
        tmp_path,
        1.0,
        hessian_recalc=None,
        energy_plateau_window=5,
        saddle_recovery_max_cycles=200,
    )
    # Simulate a stale quasi-Newton model that still claims negative curvature
    # at the true local minimum.  The exact Hessian must override it.
    opt.H = np.diag([-1.0, 10.0, 10.0])
    opt.negative_mode_seen = True
    opt.coords = [geom.coords.copy()]
    opt.cart_coords = [geom.cart_coords.copy()]
    opt.energies = [0.0] * 4
    opt.forces = [geom.forces.copy()] * 4
    opt.steps = [np.array([0.1, 0.0, 0.0])]
    opt.image_inds = [[0]]
    opt.image_nums = [1]
    opt.cur_cycle = 4

    step = opt.optimize()

    assert opt.exact_saddle_checks == 1
    assert opt._saddle_recovery_active is True
    assert opt.saddle_recovery_steps == 1
    assert opt.ts_mode_eigvals[0] > 0.0
    # Reverse the incoming motion toward x=+1, back toward the nearby saddle x=0.
    assert step[0] == -opt.saddle_recovery_step
    np.testing.assert_allclose(step[1:], 0.0)


def test_initial_near_stationary_minimum_gets_finite_recovery_step(tmp_path) -> None:
    geom, opt = _ts_optimizer(
        tmp_path, 1.0, hessian_recalc=None, saddle_recovery_max_cycles=200
    )
    opt.prepare_opt()
    opt.coords = [geom.coords.copy()]
    opt.cart_coords = [geom.cart_coords.copy()]
    opt.image_inds = [[0]]
    opt.image_nums = [1]
    opt.cur_cycle = opt.last_cycle

    step = opt.optimize()

    assert opt.exact_saddle_checks == 1
    assert opt._saddle_recovery_active is True
    assert np.linalg.norm(step) == opt.saddle_recovery_step

    # A minimizing proposal on the next cycle would point back toward the
    # minimum. Recovery must retain its chosen uphill direction instead of
    # oscillating across the minimum.
    next_step = opt.apply_saddle_recovery_step(-step)
    np.testing.assert_allclose(next_step, step)


def test_recovery_limit_rechecks_exact_curvature_before_stopping(
    tmp_path, monkeypatch
) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        1.0,
        hessian_recalc=None,
        saddle_recovery_max_cycles=2,
    )
    opt._saddle_recovery_active = True
    opt._saddle_recovery_mode = np.array([1.0, 0.0, 0.0])
    opt.saddle_recovery_cycles = 1
    opt.forces = [np.zeros(3)]

    model_hessian = np.diag([1.0, 2.0, 3.0])
    exact_hessian = np.diag([-1.0, 2.0, 3.0])
    monkeypatch.setattr(
        HessianOptimizer,
        "housekeeping",
        lambda self: (
            0.0,
            np.zeros(3),
            model_hessian,
            np.diag(model_hessian),
            np.eye(3),
            False,
        ),
    )
    monkeypatch.setattr(
        opt,
        "_refresh_exact_saddle_model",
        lambda gradient: (
            np.zeros(3),
            exact_hessian,
            np.diag(exact_hessian),
            np.eye(3),
        ),
    )
    monkeypatch.setattr(
        opt,
        "_verify_exact_vibrational_structure",
        lambda eigvals, eigvecs: (
            True,
            np.array([1.0, 0.0, 0.0]),
            True,
        ),
    )

    _, _, returned_hessian, returned_eigvals, _, resetted = opt.housekeeping()

    np.testing.assert_allclose(returned_hessian, exact_hessian)
    np.testing.assert_allclose(returned_eigvals, [-1.0, 2.0, 3.0])
    assert resetted is True
    assert opt.stop_requested is False
    assert opt.negative_mode_seen is True
    assert opt._saddle_recovery_active is False
    np.testing.assert_allclose(opt._physical_ts_mode, [1.0, 0.0, 0.0])


def test_recovery_checkpoint_continues_bounded_search_without_crossing(
    tmp_path, monkeypatch
) -> None:
    _, opt = _ts_optimizer(
        tmp_path,
        1.0,
        hessian_recalc=None,
        saddle_recovery_check_interval=2,
        saddle_recovery_max_cycles=4,
    )
    opt._saddle_recovery_active = True
    opt._saddle_recovery_mode = np.array([1.0, 0.0, 0.0])
    opt.saddle_recovery_cycles = 1
    opt.forces = [np.zeros(3)]

    positive_hessian = np.diag([1.0, 2.0, 3.0])
    monkeypatch.setattr(
        HessianOptimizer,
        "housekeeping",
        lambda self: (
            0.0,
            np.zeros(3),
            positive_hessian,
            np.diag(positive_hessian),
            np.eye(3),
            False,
        ),
    )
    monkeypatch.setattr(
        opt,
        "_refresh_exact_saddle_model",
        lambda gradient: (
            np.zeros(3),
            positive_hessian,
            np.diag(positive_hessian),
            np.eye(3),
        ),
    )
    monkeypatch.setattr(
        opt,
        "_verify_exact_vibrational_structure",
        lambda eigvals, eigvecs: (
            False,
            np.array([-1.0, 0.0, 0.0]),
            True,
        ),
    )

    opt.housekeeping()

    assert opt.stop_requested is False
    assert opt._saddle_recovery_active is True
    assert opt.saddle_recovery_cycles == 2
    # Eigenvector signs are arbitrary; checkpoint transport aligns the new
    # vector to the previous recovery direction instead of reversing motion.
    np.testing.assert_allclose(opt._saddle_recovery_mode, [1.0, 0.0, 0.0])


def test_exact_phva_ignores_negative_translation_and_uses_vibration(tmp_path) -> None:
    coords = np.array([-0.7, 0.0, 0.0, 0.7, 0.0, 0.0])
    geom = Geometry(["H", "H"], coords, coord_type="cart")
    geom.set_calculator(_TranslationalArtifactCalculator(tmp_path))
    opt = RSPRFOptimizer(
        geom,
        hessian_init="calc",
        hessian_recalc=None,
        trust_radius=0.1,
        trust_max=0.1,
        thresh="baker",
        check_eigval_structure=True,
        out_dir=tmp_path,
    )
    opt.prepare_opt()
    opt.coords = [geom.coords.copy()]
    opt.cart_coords = [geom.cart_coords.copy()]
    opt.image_inds = [[0]]
    opt.image_nums = [1]
    opt.cur_cycle = opt.last_cycle

    step = opt.optimize()

    # The unprojected Hessian still contains the artificial negative root.
    assert opt.ts_mode_eigvals[0] < 0.0
    # PHVA removes that translation and correctly identifies n_imag=0.
    assert opt._last_exact_n_imaginary == 0
    assert opt._saddle_recovery_active is False
    assert opt.stop_requested is True
    assert opt.stop_reason == "exact PHVA found no physical imaginary mode"


@pytest.mark.parametrize(
    ("optimizer_cls", "calculator_cls"),
    [
        (RSPRFOptimizer, _SpuriousTranslationDoubleWellCalculator),
        (RSIRFOptimizer, _TorchSpuriousTranslationDoubleWellCalculator),
        (TRIM, _TorchSpuriousTranslationDoubleWellCalculator),
    ],
)
def test_recovery_tracks_physical_mode_past_more_negative_tr_artifact(
    tmp_path, optimizer_cls, calculator_cls
) -> None:
    a = _SpuriousTranslationDoubleWellCalculator.well_position
    coords = np.array([-a / np.sqrt(2.0), 0.0, 0.0, a / np.sqrt(2.0), 0.0, 0.0])
    geom = Geometry(["H", "H"], coords, coord_type="cart")
    geom.set_calculator(calculator_cls(tmp_path))
    opt = optimizer_cls(
        geom,
        hessian_init="calc",
        hessian_recalc=1,
        trust_radius=0.05,
        trust_update=False,
        trust_min=1e-5,
        trust_max=0.05,
        max_cycles=40,
        saddle_recovery_max_cycles=200,
        thresh="baker",
        check_eigval_structure=True,
        out_dir=tmp_path,
    )

    opt.run()

    stretch = np.array([-1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
    stretch /= np.linalg.norm(stretch)
    final_q = float(stretch @ geom.coords)
    assert opt.is_converged is True
    assert abs(final_q) < 1e-3
    assert opt.saddle_recovery_steps > 0
    assert opt._last_exact_n_imaginary == 1
    # Root 0 remains the artificial translation (-1); the physical reaction
    # vibration is root 1 and must be followed after recovery.
    assert int(opt.roots[0]) == 1
