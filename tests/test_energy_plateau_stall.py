"""An energy-only plateau is an additive 'stalled' outcome,
never convergence, and it stops further retry work.

These bind to the production optimizer/renderer code paths (no reimplemented
gate logic): the bundled ``pysisyphus`` Optimizer/RFOptimizer/TSHessianOptimizer
state machine and the product-local status helpers in ``mlmm.core.utils`` /
``mlmm.workflows.tsopt``.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.Optimizer import CONV_THRESHS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

import mlmm.workflows.tsopt as tsopt_mod
from mlmm.core.utils import (
    optimizer_terminal_status,
    emit_optimizer_terminal_status,
)
from mlmm.workflows.tsopt import HessianDimer, _tsopt_terminal_status


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


class _ConstantPlateauCalculator(Calculator):
    """Constant energy (an energy plateau) with a constant large force so the
    configured force/step criteria can never be met."""

    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    @staticmethod
    def _results(coords):
        coords = np.asarray(coords, dtype=float)
        return {
            "energy": 1.0,
            "forces": np.array([0.5, 0.0, 0.0]),
            "hessian": np.eye(len(coords)),
        }

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)


def _lbfgs(tmp_path, **kwargs):
    geom = Geometry(["H"], np.zeros(3), coord_type="cart")
    geom.set_calculator(_QuadraticCalculator(tmp_path))
    opt_kwargs = dict(
        thresh="gau",
        energy_plateau=True,
        energy_plateau_window=3,
        out_dir=tmp_path,
    )
    opt_kwargs.update(kwargs)
    return geom, LBFGS(geom, **opt_kwargs)


def _seed(opt, *, force, step, energies):
    opt.cur_cycle = 5
    opt.last_cycle = 0
    opt.forces = [np.asarray(force, dtype=float)]
    opt.steps = [np.asarray(step, dtype=float)]
    opt.energies = list(energies)


_HIGH = np.full(3, 1.0)
_LOW = np.full(3, 1.0e-8)
_BETWEEN_GAU_AND_OLD = np.full(3, 7.0e-4)
_FLAT = [1.0, 1.0, 1.0]  # range 0 over window 3 -> a plateau


# ---- Base state machine: force/step truth table under a flat energy window ----

@pytest.mark.parametrize(
    "force, step, expect_converged",
    [
        (_HIGH, _HIGH, False),   # high force, high step   -> stalled
        (_HIGH, _LOW, False),    # high force, low step    -> stalled
        (_LOW, _HIGH, False),    # low force, high step    -> stalled
        (_BETWEEN_GAU_AND_OLD, _LOW, False),
        (_LOW, _LOW, True),      # all below thresh        -> converged (wins)
    ],
)
def test_energy_plateau_truth_table(tmp_path, force, step, expect_converged):
    _, opt = _lbfgs(tmp_path)
    _seed(opt, force=force, step=step, energies=_FLAT)

    converged, conv_info = opt.check_convergence(step)

    assert bool(converged) is expect_converged
    if expect_converged:
        # A real convergence wins over the plateau; no stall. (check_convergence
        # returns the bool; the run loop is what assigns self.is_converged.)
        assert opt.is_stalled is False
    else:
        # A plateau without real convergence stalls, never converges.  The
        # stall side effect (request_stall) sets the committed terminal state.
        assert opt.is_stalled is True
        assert opt.is_converged is False
        assert opt.stop_requested is True
        assert opt.termination_status == "stalled"
        assert "energy plateau" in opt.stop_reason
        # A plateau does not alter the ConvInfo force/step fields.
        max_force_thresh = CONV_THRESHS["gau"][0]
        assert bool(conv_info.max_force_converged) == bool(
            np.max(np.abs(force)) <= max_force_thresh
        )


def test_energy_range_just_above_threshold_neither_converges_nor_stalls(tmp_path):
    _, opt = _lbfgs(tmp_path)
    # range 2e-5 over the window > energy_plateau_thresh (1e-5): not a plateau.
    _seed(opt, force=_HIGH, step=_HIGH, energies=[1.0, 1.0, 1.0 + 2.0e-5])

    converged, _ = opt.check_convergence(_HIGH)

    assert converged is False
    assert opt.is_stalled is False
    assert opt.termination_status == "not_converged"


def test_energy_plateau_disabled_does_not_stall(tmp_path):
    _, opt = _lbfgs(tmp_path, energy_plateau=False)
    _seed(opt, force=_HIGH, step=_HIGH, energies=_FLAT)

    converged, _ = opt.check_convergence(_HIGH)

    assert converged is False
    assert opt.is_stalled is False


def test_too_few_energy_samples_does_not_stall(tmp_path):
    _, opt = _lbfgs(tmp_path, energy_plateau_window=5)
    _seed(opt, force=_HIGH, step=_HIGH, energies=[1.0, 1.0])  # < window

    converged, _ = opt.check_convergence(_HIGH)

    assert converged is False
    assert opt.is_stalled is False


def test_thresh_never_neither_converges_nor_stalls(tmp_path):
    _, opt = _lbfgs(tmp_path, thresh="never")
    _seed(opt, force=_LOW, step=_LOW, energies=_FLAT)

    converged, _ = opt.check_convergence(_LOW)

    assert converged is False
    assert opt.is_stalled is False


# ---- RFOptimizer provisional probe must not stall; the final check must ------

def test_rfo_provisional_probe_suppresses_stall_but_final_check_stalls(tmp_path):
    geom = Geometry(["H"], np.zeros(3), coord_type="cart")
    geom.set_calculator(_QuadraticCalculator(tmp_path))
    opt = RFOptimizer(
        geom,
        thresh="gau",
        energy_plateau=True,
        energy_plateau_window=3,
        out_dir=tmp_path,
    )
    _seed(opt, force=_HIGH, step=_HIGH, energies=_FLAT)

    # Provisional probe (RFOptimizer's ref_step probe uses this) must NOT stall.
    provisional, _ = opt.check_convergence(_HIGH, allow_stall=False)
    assert provisional is False
    assert opt.is_stalled is False

    # The subsequent run-loop final check (allow_stall default True) stalls.
    final, _ = opt.check_convergence(_HIGH)
    assert final is False
    assert opt.is_stalled is True
    assert opt.termination_status == "stalled"


# ---- Run-loop integration: stall, no further step, terminal says stalled -----

def test_run_loop_stalls_on_energy_plateau(tmp_path, capsys):
    geom = Geometry(["H"], np.zeros(3), coord_type="cart")
    geom.set_calculator(_ConstantPlateauCalculator(tmp_path))
    opt = LBFGS(
        geom,
        thresh="gau",
        energy_plateau=True,
        energy_plateau_window=2,
        max_cycles=20,
        out_dir=tmp_path,
    )
    opt.run()

    assert opt.is_converged is False
    assert opt.is_stalled is True
    assert opt.stop_requested is True
    assert opt.stopped is True
    assert opt.termination_status == "stalled"
    # Exited well before max_cycles (stalled after the window filled).
    assert opt.cur_cycle < 19
    out = capsys.readouterr().out
    assert "Stalled" in out
    assert "Converged!" not in out


# ---- Product-local status renderers ------------------------------------------

class _FakeOpt:
    def __init__(self, *, is_converged=False, is_stalled=False, stop_reason=""):
        self.is_converged = is_converged
        self.is_stalled = is_stalled
        self.stop_reason = stop_reason

    @property
    def termination_status(self):
        if self.is_stalled:
            return "stalled"
        if self.is_converged:
            return "converged"
        return "not_converged"


def test_optimizer_terminal_status_maps_stalled_and_converged():
    assert optimizer_terminal_status(_FakeOpt(is_stalled=True)) == "stalled"
    assert optimizer_terminal_status(_FakeOpt(is_converged=True)) == "converged"
    assert optimizer_terminal_status(_FakeOpt()) == "not_converged"


def test_stalled_optimizer_never_reports_converged():
    stalled = _FakeOpt(is_stalled=True, is_converged=False, stop_reason="energy plateau: ...")
    assert optimizer_terminal_status(stalled) == "stalled"
    assert optimizer_terminal_status(stalled) != "converged"
    # A stalled TS optimizer is never a converged saddle, even at n_imag==1.
    assert _tsopt_terminal_status(stalled, saddle_verified=True) == "stalled"


def test_tsopt_terminal_status_composition():
    converged = _FakeOpt(is_converged=True)
    assert _tsopt_terminal_status(converged, saddle_verified=True) == "converged"
    # Numerical convergence is independent of saddle order.
    assert _tsopt_terminal_status(converged, saddle_verified=False) == "converged"
    assert _tsopt_terminal_status(_FakeOpt(), saddle_verified=True) == "not_converged"


def test_emit_terminal_status_stalled_and_converged_are_distinct(capsys):
    emit_optimizer_terminal_status(
        "opt", converged=False, cycles=3, max_cycles=20,
        stalled=True, stop_reason="energy plateau: range=0.00e+00",
    )
    stalled_out = capsys.readouterr().out
    assert "Stalled" in stalled_out
    assert "Converged!" not in stalled_out

    emit_optimizer_terminal_status(
        "opt", converged=True, cycles=7, max_cycles=20,
    )
    conv_out = capsys.readouterr().out
    assert "Converged!" in conv_out
    assert "Stalled" not in conv_out


def test_tsopt_terminal_status_labels_numerical_convergence(capsys):
    emit_optimizer_terminal_status(
        "tsopt",
        converged=True,
        cycles=7,
        max_cycles=20,
        converged_message="Numerical optimization converged.",
    )
    out = capsys.readouterr().out
    assert "[tsopt] Numerical optimization converged." in out
    assert "[tsopt] Converged!" not in out


# ---- HessianDimer wrapper: a stalled child stops all further work -------------

@pytest.mark.parametrize("cadence", [0, -1])
def test_hessian_dimer_rejects_nonpositive_cadence_before_output_mutation(
    tmp_path, cadence
):
    out_dir = tmp_path / "not-created"

    with pytest.raises(ValueError, match="must be at least 1"):
        HessianDimer(
            fn=tmp_path / "missing.pdb",
            out_dir=out_dir,
            update_interval_hessian=cadence,
        )

    assert not out_dir.exists()


def test_hessian_dimer_stops_after_child_stall(tmp_path, monkeypatch):
    """A stalled child LBFGS makes the runner stalled and stops the segment
    loop before any further segment or Hessian update."""
    runner = HessianDimer.__new__(HessianDimer)
    runner.max_total_cycles = 100
    runner._cycles_spent = 0
    runner.update_interval_hessian = 5
    runner.is_stalled = False
    runner.is_converged = False
    runner.stop_reason = ""

    class _Geom:
        cart_coords = np.zeros(6)

    runner.geom = _Geom()

    calls = {"segments": 0}

    def _fake_segment(threshold, n_steps):
        calls["segments"] += 1
        # The child LBFGS stalled on an energy plateau.
        runner.is_stalled = True
        runner.stop_reason = "energy plateau: range=0.00e+00 au over 2 steps"
        return 3, False  # (steps, converged=False)

    def _raise_hessian(*args, **kwargs):
        raise AssertionError("no Hessian update should run after a stall")

    runner._dimer_segment = _fake_segment
    # The mlmm dimer loop updates the mode via the module-level Hessian helper;
    # it must never be reached after a child stall.
    monkeypatch.setattr(tsopt_mod, "_calc_full_hessian_torch", _raise_hessian)

    steps, zero_step_converged, loop_converged = runner._dimer_loop("gau")

    assert runner.is_stalled is True
    assert runner.termination_status == "stalled"
    assert loop_converged is False
    assert calls["segments"] == 1          # only the stalled segment ran
    assert steps == 3
    # And the public status mapper reports stalled, never converged.
    assert _tsopt_terminal_status(runner, saddle_verified=True) == "stalled"


def test_terminal_saddle_certification_uses_magnitude_threshold():
    from mlmm.workflows.tsopt import (
        _certified_negative_frequencies,
        _certified_saddle_order,
        _finalize_dimer_saddle_status,
    )

    freqs_cm = np.array([-450.0, -3.2, 12.0, 640.0])
    assert _certified_saddle_order(freqs_cm, 5.0) == 1
    assert _certified_negative_frequencies(freqs_cm, 5.0) == [-450.0]

    runner = _FakeOpt(is_converged=True)
    export_idx = _finalize_dimer_saddle_status(runner, freqs_cm, 5.0)
    assert runner.n_imaginary_modes == 1
    assert runner.imaginary_frequencies_cm == [-450.0]
    assert runner.saddle_order_verified is True
    assert runner.is_converged is True
    assert export_idx.tolist() == [0]

    soft = _FakeOpt(is_converged=True)
    soft_export = _finalize_dimer_saddle_status(
        soft, np.array([-3.2, 12.0, 640.0]), 5.0
    )
    assert soft.n_imaginary_modes == 0
    assert soft.imaginary_frequencies_cm == []
    assert soft.saddle_order_verified is False
    assert soft_export.tolist() == []


def test_exact_phva_validation_ignores_soft_negative_roots():
    from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer

    modes = np.eye(3)
    printed: list[str] = []
    optimizer = RSIRFOptimizer.__new__(RSIRFOptimizer)
    optimizer.saddle_imaginary_threshold_cm = 5.0
    optimizer.small_eigval_thresh = 1e-8
    optimizer.roots = np.array([0])
    optimizer.reference_mode = None
    optimizer.cur_cycle = 7
    optimizer.higher_order_saddle_checks = 0
    optimizer.max_higher_order_checks = 99
    optimizer.forces = []
    optimizer.geometry = SimpleNamespace(cart_coords=np.zeros(3))
    optimizer.table = SimpleNamespace(print=printed.append)
    optimizer._mw_frequencies_and_modes = lambda: (
        np.array([-450.0, -3.2, 12.0]),
        modes,
    )
    optimizer._recovery_mode_from_mw = (
        lambda _modes, index: modes[:, int(index)]
    )
    optimizer._record_exact_saddle_candidate = lambda: None
    optimizer.request_stop = lambda *_a, **_k: None

    has_saddle_modes, _mode, _has_mode = (
        optimizer._verify_exact_vibrational_structure(
            np.array([-0.1, -1.0e-7, 0.2]), np.eye(3)
        )
    )
    assert optimizer._last_exact_n_imaginary == 1
    assert optimizer._last_exact_saddle_verified is True
    assert has_saddle_modes is True
    assert any("n_imag=1" in message for message in printed)


def test_dimer_final_message_separates_no_mode_from_write_failure():
    from mlmm.workflows.tsopt import (
        _dimer_mode_export_message,
        _unexpected_saddle_order_message,
    )

    positive, positive_is_diagnostic = _dimer_mode_export_message(0, 0, 5.0, 12.0)
    assert positive == "[tsopt] No imaginary mode detected. Try all --refine-path."
    assert positive_is_diagnostic is True

    failed, failed_is_diagnostic = _dimer_mode_export_message(0, 1, 5.0, -100.0)
    assert failed == "[tsopt] ERROR: Failed to write imaginary mode trajectory."
    assert failed_is_diagnostic is True

    assert _unexpected_saddle_order_message(2) == (
        "[tsopt] WARNING: Higher-order stationary point (n_imag=2). "
        "Try --flatten or all --refine-path."
    )


# ---- Microiteration terminal outcome (opt + tsopt share the helper) ----------

def _real_terminal_optimizer(tmp_path, **state):
    """A real pysisyphus optimizer forced into a chosen terminal state.

    Binds the microiteration finalizer to the production ``request_stall`` state
    machine (no reimplemented gate logic).
    """
    geom = Geometry(["H"], np.zeros(3), coord_type="cart")
    geom.set_calculator(_QuadraticCalculator(tmp_path))
    opt = LBFGS(geom, thresh="gau", max_cycles=5, out_dir=tmp_path)
    for key, value in state.items():
        setattr(opt, key, value)
    return opt


def test_microiter_surfaces_latest_micro_stall_on_max_cycles(tmp_path):
    """A macro that merely ran out of cycles while the latest micro
    (MM) relaxation stalled is surfaced as ``stalled`` with the micro reason,
    not a reasonless ``not_converged``."""
    from mlmm.core.utils import finalize_microiter_macro_convergence

    opt = _real_terminal_optimizer(tmp_path, is_converged=False)
    macro_conv = finalize_microiter_macro_convergence(
        opt,
        macro_converged=False,           # macro reached max_cycles
        latest_micro_stalled=True,       # final MM relaxation stalled
        latest_micro_stop_reason="energy plateau: range=0.00e+00 au over 2 steps",
    )
    assert macro_conv is False
    assert opt.is_stalled is True
    assert opt.stop_requested is True
    assert opt.stop_reason == "energy plateau: range=0.00e+00 au over 2 steps"
    assert opt.termination_status == "stalled"


def test_microiter_demotes_converged_macro_on_latest_micro_stall(tmp_path):
    """A would-be macro convergence with
    a stalled latest micro relaxation is demoted to ``stalled``."""
    from mlmm.core.utils import finalize_microiter_macro_convergence

    opt = _real_terminal_optimizer(tmp_path, is_converged=True)
    macro_conv = finalize_microiter_macro_convergence(
        opt,
        macro_converged=True,
        latest_micro_stalled=True,
        latest_micro_stop_reason="mm plateau",
    )
    assert macro_conv is False
    assert opt.is_stalled is True
    assert opt.stop_reason == "mm plateau"


def test_microiter_converges_when_latest_micro_did_not_stall(tmp_path):
    """A converged macro with a non-stalled latest micro stays ``converged``."""
    from mlmm.core.utils import finalize_microiter_macro_convergence

    opt = _real_terminal_optimizer(tmp_path, is_converged=True)
    macro_conv = finalize_microiter_macro_convergence(
        opt,
        macro_converged=True,
        latest_micro_stalled=False,
    )
    assert macro_conv is True
    assert opt.is_stalled is False


def test_microiter_keeps_macro_own_stop_reason_over_micro(tmp_path):
    """A macro that made its own clean stop keeps its more specific reason and is
    not overwritten by the (also-stalled) latest micro relaxation."""
    from mlmm.core.utils import finalize_microiter_macro_convergence

    opt = _real_terminal_optimizer(
        tmp_path,
        is_converged=False,
        stop_requested=True,
        stop_reason="repeated uphill RFO trials at the emergency trust floor",
    )
    macro_conv = finalize_microiter_macro_convergence(
        opt,
        macro_converged=False,
        latest_micro_stalled=True,
        latest_micro_stop_reason="mm plateau",
    )
    assert macro_conv is False
    assert opt.is_stalled is False
    assert opt.stop_reason == "repeated uphill RFO trials at the emergency trust floor"


def test_opt_result_converged_sourced_from_microiter_result():
    """On the microiteration path the opt result.json convergence flag
    is sourced from ``microiter_result['converged']`` (there is no standalone
    optimizer in scope), so a converged microiter run is no longer mislabeled
    ``not_converged``. Feeding this non-stalled ``True`` into run_opt's inline
    status ternary yields ``status='converged'`` (console + JSON agree); a
    stalled/​not-converged microiter result yields ``stalled`` / ``not_converged``.
    """
    from mlmm.workflows.opt import _opt_terminal_converged

    # Microiteration path: converged truth lives in microiter_result.
    assert _opt_terminal_converged(True, {"converged": True}, None) is True
    assert _opt_terminal_converged(True, {"converged": False, "is_stalled": True}, None) is False
    assert _opt_terminal_converged(True, {"converged": False}, None) is False
    # Missing runner -> unknown (None), never a spurious converged.
    assert _opt_terminal_converged(True, None, None) is None
    # Standard (non-microiter) path still reads the optimizer object.
    assert _opt_terminal_converged(False, None, _FakeOpt(is_converged=True)) is True
    assert _opt_terminal_converged(False, None, _FakeOpt(is_converged=False)) is False
