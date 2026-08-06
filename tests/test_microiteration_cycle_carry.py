"""Microiteration cycle-count and vacuous-partition contracts.

The selected restart owns the serialized microiteration outcome, ordinary
optimization reports executed rather than zero-based cycle counts, and a
partition with no micro-active MM atoms records a zero-cycle micro success
without constructing an all-frozen optimizer.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11), reason="mlmm requires Python >= 3.11"
)

from mlmm.workflows._microiteration import (
    MicroiterationOutcome,
    OptimizerOutcome,
    build_aggregate,
    build_partition,
)


# ---------------------------------------------------------------------------
# Restart selection re-anchors the microiteration block
# ---------------------------------------------------------------------------


def _converged_micro_outcome(macro_cycles: int, micro_cycles: int) -> MicroiterationOutcome:
    macro = OptimizerOutcome(
        status="converged", executed=True, converged=True, cycles=macro_cycles
    )
    micro = OptimizerOutcome(
        status="converged", executed=True, converged=True, cycles=micro_cycles
    )
    agg = build_aggregate(macro, [micro], max_cycles=100)
    return MicroiterationOutcome(
        aggregate=agg,
        macro=macro,
        micro_attempts=(micro,),
        macro_cycles=macro_cycles,
        micro_cycles=micro_cycles,
    )


def _notconverged_micro_outcome() -> MicroiterationOutcome:
    macro = OptimizerOutcome(
        status="not_converged", executed=True, converged=False, cycles=10
    )
    agg = build_aggregate(macro, [], max_cycles=10)
    return MicroiterationOutcome(
        aggregate=agg, macro=macro, micro_attempts=(), macro_cycles=10, micro_cycles=0
    )


def test_restart_carry_serializes_selected_run_outcome():
    from mlmm.workflows.tsopt import _restart_microiteration_carry

    restart = _converged_micro_outcome(macro_cycles=4, micro_cycles=3)
    outcome = {"outcome": restart, "micro_cycles": 3}
    obj, micro_cycles = _restart_microiteration_carry(outcome)
    assert obj["aggregate"]["status"] == "converged"
    assert obj["micro_cycles"] == 3
    assert micro_cycles == 3


def test_restart_carry_ordinary_path_carries_nothing():
    from mlmm.workflows.tsopt import _restart_microiteration_carry

    # An ordinary (non-microiter) restart's outcome has no MicroiterationOutcome.
    assert _restart_microiteration_carry({"outcome": None}) == (None, None)
    assert _restart_microiteration_carry({}) == (None, None)


def test_selected_converged_restart_overrides_initial_notconverged_leaf():
    """An initial microiteration that did not converge,
    then a converging restart, must serialize the SELECTED (converged) run's
    aggregate status -- never the superseded initial not_converged leaf."""
    from mlmm.workflows.tsopt import _restart_microiteration_carry

    initial_obj = _notconverged_micro_outcome().to_result_object()
    assert initial_obj["aggregate"]["status"] == "not_converged"

    restart_outcome = {"outcome": _converged_micro_outcome(4, 1), "micro_cycles": 1}
    selected_obj, selected_micro = _restart_microiteration_carry(restart_outcome)

    # The carried block describes the converged restart, not the initial leaf.
    assert selected_obj["aggregate"]["status"] == "converged"
    assert selected_obj["aggregate"]["status"] != initial_obj["aggregate"]["status"]
    assert selected_micro == 1


# ---------------------------------------------------------------------------
# Ordinary-opt n_opt_cycles is the executed count (cur_cycle + 1)
# ---------------------------------------------------------------------------


def test_ordinary_opt_cycle_count_is_executed_not_zero_based():
    """A 1-cycle converged opt executed ONE cycle. The ordinary-opt JSON
    ``n_opt_cycles`` now uses ``optimizer_cycle_count`` (== console log ==
    tsopt), never the raw zero-based ``cur_cycle`` (which reported 0)."""
    from mlmm.core.utils import optimizer_cycle_count

    class _Opt:
        cur_cycle = 0

    class _Opt5:
        cur_cycle = 5

    assert optimizer_cycle_count(_Opt()) == 1
    assert optimizer_cycle_count(_Opt5()) == 6


# ---------------------------------------------------------------------------
# A zero-micro-active partition uses the vacuous micro, not an
# all-frozen LBFGS (bound end-to-end in BOTH drivers)
# ---------------------------------------------------------------------------


class _FakeCore:
    _ml_backend = object()


class _FakeCalc:
    def __init__(self, *args, **kwargs):
        self.core = _FakeCore()
        self.freeze_atoms = []


class _FakeConvInfo:
    def get_convergence(self):
        # >= number of table columns so the marks list is long enough.
        return [False] * 8


class _FakeMacroOptimizer:
    """A macro optimizer that never converges (so the in-loop micro step is
    reached), with just enough surface for the microiteration driver loop."""

    def __init__(self, geom, **kwargs):
        self.geom = geom
        self.coords = []
        self.cart_coords = []
        self.steps = []
        self.cur_cycle = 0
        self.stop_requested = False
        self.is_stalled = False
        self.is_converged = False
        self.stop_reason = ""
        self.energies = [-1.0]
        self.max_forces = [0.1]
        self.rms_forces = [0.05]
        self.max_steps = [0.01]
        self.rms_steps = [0.005]

    def prepare_opt(self):
        return None

    def optimize(self):
        return np.zeros_like(self.geom.coords)

    def check_convergence(self):
        return False, _FakeConvInfo()


class _FakeGeom:
    def __init__(self, n_atoms=2):
        self.atoms = tuple("H" for _ in range(n_atoms))
        self.atomic_numbers = [1] * n_atoms
        self._coords = np.zeros(3 * n_atoms, dtype=float)
        self.freeze_atoms = []
        self.cart_hessian = None
        self.within_partial_hessian = None
        self.energy = -1.0
        self._calc = None

    @property
    def coords(self):
        return self._coords

    @coords.setter
    def coords(self, value):
        self._coords = np.asarray(value, dtype=float)

    @property
    def cart_coords(self):
        return self._coords

    @cart_coords.setter
    def cart_coords(self, value):
        self._coords = np.asarray(value, dtype=float)

    def set_calculator(self, calc):
        self._calc = calc
        if hasattr(calc, "freeze_atoms"):
            calc.freeze_atoms = list(self.freeze_atoms)

    def as_xyz(self):
        return ""


def _no_micro_active_partition():
    """macro-active {0}, micro-active {} (atom 1 user-frozen)."""
    part = build_partition(
        2, ml=[0], link_parents=[], hess_mm=[], movable_mm=[], frozen_mm=[],
        original_freeze=[1],
    )
    assert part.has_macro_active is True
    assert part.has_micro_active is False
    return part


def _micro_active_partition():
    """macro-active {0}, micro-active {1}."""
    part = build_partition(
        2,
        ml=[0],
        link_parents=[],
        hess_mm=[1],
        movable_mm=[1],
        frozen_mm=[],
        original_freeze=[],
    )
    assert part.has_macro_active is True
    assert part.has_micro_active is True
    return part


def _install_common_fakes(monkeypatch, driver_mod, tripwire):
    import torch

    monkeypatch.setattr(driver_mod, "mlmm", lambda *a, **k: _FakeCalc())
    monkeypatch.setattr(driver_mod, "mlmm_mm_only", lambda *a, **k: _FakeCalc())
    monkeypatch.setattr(driver_mod, "LBFGS", tripwire)

    import mlmm.io.hessian_cache as _hc

    monkeypatch.setattr(_hc, "load_matching", lambda *a, **k: None)
    monkeypatch.setattr(_hc, "identity_from_context", lambda *a, **k: {})

    import mlmm.workflows.freq as _freq

    monkeypatch.setattr(
        _freq, "_calc_full_hessian_torch",
        lambda *a, **k: (torch.zeros((3, 3), dtype=torch.float64), None),
    )


def test_opt_zero_micro_active_uses_vacuous_micro_not_all_frozen_lbfgs(tmp_path, monkeypatch):
    import mlmm.workflows.opt as opt_mod

    def _lbfgs_tripwire(*args, **kwargs):
        raise AssertionError(
            "no LBFGS may be built when the partition has no micro-active DOF"
        )

    _install_common_fakes(monkeypatch, opt_mod, _lbfgs_tripwire)
    monkeypatch.setattr(opt_mod, "RFOptimizer", _FakeMacroOptimizer)

    geom = _FakeGeom(n_atoms=2)
    part = _no_micro_active_partition()

    outcome = opt_mod._run_microiter_opt(
        geom,
        _FakeCalc(),
        calc_cfg={},
        rfo_cfg={},
        lbfgs_cfg={},
        opt_cfg={"max_cycles": 1},
        microiter_cfg={},
        out_dir_path=tmp_path,
        partition=part,
        dump=False,
    )

    assert outcome["micro_cycles"] == 0
    mi = outcome["outcome"]
    assert len(mi.micro_attempts) >= 1
    for micro in mi.micro_attempts:
        assert micro.status == "converged"
        assert micro.cycles == 0
        assert micro.stalled is False
        assert micro.stop_reason is None
    assert mi.to_result_object()["micro_cycles"] == 0


def test_opt_initial_hessian_is_resolved_after_initial_mm_relaxation(
    tmp_path, monkeypatch
) -> None:
    import torch
    import mlmm.workflows.opt as opt_mod
    import mlmm.workflows.freq as freq_mod

    class _MovingMicro:
        cur_cycle = 0
        is_converged = True
        is_stalled = False
        stop_reason = None
        calls = 0

        def __init__(self, geom, **_kwargs):
            self.geom = geom

        def run(self):
            type(self).calls += 1
            self.geom.coords = np.full_like(
                self.geom.coords, 7.0 if type(self).calls == 1 else 9.0
            )

    _install_common_fakes(monkeypatch, opt_mod, _MovingMicro)
    monkeypatch.setattr(opt_mod, "RFOptimizer", _FakeMacroOptimizer)
    hessian_coords = []
    hessian_freeze_masks = []

    def _capture_hessian(geom, *_args, **kwargs):
        hessian_coords.append(np.asarray(geom.coords, dtype=float).copy())
        hessian_freeze_masks.append(
            (
                tuple(geom.freeze_atoms),
                tuple(kwargs["calculator"].freeze_atoms),
            )
        )
        return torch.zeros((3, 3), dtype=torch.float64), None

    monkeypatch.setattr(freq_mod, "_calc_full_hessian_torch", _capture_hessian)

    partition = _micro_active_partition()
    opt_mod._run_microiter_opt(
        _FakeGeom(n_atoms=2),
        _FakeCalc(),
        calc_cfg={},
        rfo_cfg={},
        lbfgs_cfg={},
        opt_cfg={"max_cycles": 1},
        microiter_cfg={"micro_max_cycles": 1},
        out_dir_path=tmp_path,
        partition=partition,
        dump=False,
    )

    assert len(hessian_coords) == 1
    np.testing.assert_allclose(hessian_coords[0], 7.0)
    assert hessian_freeze_masks == [
        (partition.macro_freeze_atoms, partition.macro_freeze_atoms)
    ]


def test_opt_initial_micro_failure_skips_hessian_and_macro(
    tmp_path, monkeypatch
) -> None:
    import mlmm.workflows.opt as opt_mod
    import mlmm.workflows.freq as freq_mod

    class _NonconvergedMicro:
        cur_cycle = 0
        is_converged = False
        is_stalled = False
        stop_reason = "maximum cycles reached"

        def __init__(self, _geom, **_kwargs):
            pass

        def run(self):
            return None

    _install_common_fakes(monkeypatch, opt_mod, _NonconvergedMicro)
    monkeypatch.setattr(
        freq_mod,
        "_calc_full_hessian_torch",
        lambda *a, **k: (_ for _ in ()).throw(
            AssertionError("Hessian must not run after initial micro failure")
        ),
    )
    monkeypatch.setattr(
        opt_mod,
        "RFOptimizer",
        lambda *a, **k: (_ for _ in ()).throw(
            AssertionError("macro optimizer must not be constructed")
        ),
    )

    outcome = opt_mod._run_microiter_opt(
        _FakeGeom(n_atoms=2),
        _FakeCalc(),
        calc_cfg={},
        rfo_cfg={},
        lbfgs_cfg={},
        opt_cfg={"max_cycles": 1},
        microiter_cfg={"micro_max_cycles": 1},
        out_dir_path=tmp_path,
        partition=_micro_active_partition(),
        dump=False,
    )

    assert outcome["cycles"] == 0
    assert outcome["optimizer"] is None
    assert outcome["outcome"].macro.executed is False


def test_tsopt_zero_micro_active_uses_vacuous_micro_not_all_frozen_lbfgs(tmp_path, monkeypatch):
    import mlmm.workflows.tsopt as tsopt_mod
    from mlmm.core.defaults import RSIRFO_KW

    def _lbfgs_tripwire(*args, **kwargs):
        raise AssertionError(
            "no LBFGS may be built when the partition has no micro-active DOF"
        )

    _install_common_fakes(monkeypatch, tsopt_mod, _lbfgs_tripwire)
    monkeypatch.setattr(
        tsopt_mod, "_calc_full_hessian_torch",
        lambda *a, **k: __import__("torch").zeros((3, 3), dtype=__import__("torch").float64),
    )
    part = _no_micro_active_partition()
    monkeypatch.setattr(tsopt_mod, "resolve_partition_from_core", lambda *a, **k: part)
    monkeypatch.setitem(tsopt_mod.TSOPT_CLASS_MAP, "rsirfo", _FakeMacroOptimizer)

    geom = _FakeGeom(n_atoms=2)

    outcome = tsopt_mod._run_microiter_tsopt(
        geom,
        {},
        dict(RSIRFO_KW),
        {},
        {"max_cycles": 1, "dump": False},
        {},
        tmp_path,
        dump=False,
        thresh=None,
        mode="rsirfo",
        reference_mode=None,
    )

    assert outcome["micro_cycles"] == 0
    mi = outcome["outcome"]
    assert len(mi.micro_attempts) >= 1
    for micro in mi.micro_attempts:
        assert micro.status == "converged"
        assert micro.cycles == 0
        assert micro.stalled is False
        assert micro.stop_reason is None
    assert mi.to_result_object()["micro_cycles"] == 0


def test_tsopt_dump_keeps_initial_micro_trajectory_when_macro_never_runs(
    tmp_path, monkeypatch
):
    import mlmm.workflows.tsopt as tsopt_mod
    from mlmm.core.defaults import RSIRFO_KW

    class _NonconvergedMicro:
        cur_cycle = 0
        is_converged = False
        is_stalled = False
        stop_reason = "maximum cycles reached"

        def __init__(self, _geom, **kwargs):
            self.out_dir = Path(kwargs["out_dir"])

        def run(self):
            (self.out_dir / "optimization_trj.xyz").write_text(
                "2\ninitial micro\nH 0 0 0\nH 0 0 1\n",
                encoding="utf-8",
            )

    _install_common_fakes(monkeypatch, tsopt_mod, _NonconvergedMicro)
    monkeypatch.setattr(
        tsopt_mod,
        "_calc_full_hessian_torch",
        lambda *a, **k: (_ for _ in ()).throw(
            AssertionError("Hessian must not run after initial micro failure")
        ),
    )
    monkeypatch.setattr(
        tsopt_mod,
        "resolve_partition_from_core",
        lambda *a, **k: _micro_active_partition(),
    )
    monkeypatch.setitem(
        tsopt_mod.TSOPT_CLASS_MAP,
        "rsirfo",
        lambda *a, **k: (_ for _ in ()).throw(
            AssertionError("macro optimizer must not be constructed")
        ),
    )

    outcome = tsopt_mod._run_microiter_tsopt(
        _FakeGeom(n_atoms=2),
        {},
        dict(RSIRFO_KW),
        {},
        {"max_cycles": 1, "dump": True},
        {"micro_max_cycles": 1},
        tmp_path,
        dump=True,
        thresh=None,
        mode="rsirfo",
        reference_mode=None,
    )

    assert outcome["cycles"] == 0
    assert outcome["optimizer"] is None
    assert outcome["outcome"].macro.executed is False
    assert (tmp_path / "optimization_all_trj.xyz").read_text(
        encoding="utf-8"
    ) == (tmp_path / "optimization_trj.xyz").read_text(encoding="utf-8")


def test_tsopt_dump_does_not_duplicate_converged_initial_micro_trajectory(
    tmp_path, monkeypatch
):
    import mlmm.workflows.tsopt as tsopt_mod
    from mlmm.core.defaults import RSIRFO_KW

    class _ConvergedMicro:
        cur_cycle = 0
        is_converged = True
        is_stalled = False
        stop_reason = None
        calls = 0

        def __init__(self, _geom, **kwargs):
            self.out_dir = Path(kwargs["out_dir"])

        def run(self):
            type(self).calls += 1
            label = f"micro {type(self).calls}"
            (self.out_dir / "optimization_trj.xyz").write_text(
                f"2\n{label}\nH 0 0 0\nH 0 0 1\n",
                encoding="utf-8",
            )

    _install_common_fakes(monkeypatch, tsopt_mod, _ConvergedMicro)
    monkeypatch.setattr(
        tsopt_mod,
        "_calc_full_hessian_torch",
        lambda *a, **k: __import__("torch").zeros(
            (3, 3), dtype=__import__("torch").float64
        ),
    )
    monkeypatch.setattr(
        tsopt_mod,
        "resolve_partition_from_core",
        lambda *a, **k: _micro_active_partition(),
    )
    monkeypatch.setitem(tsopt_mod.TSOPT_CLASS_MAP, "rsirfo", _FakeMacroOptimizer)

    tsopt_mod._run_microiter_tsopt(
        _FakeGeom(n_atoms=2),
        {},
        dict(RSIRFO_KW),
        {},
        {"max_cycles": 1, "dump": True},
        {"micro_max_cycles": 1},
        tmp_path,
        dump=True,
        thresh=None,
        mode="rsirfo",
        reference_mode=None,
    )

    combined = (tmp_path / "optimization_all_trj.xyz").read_text(encoding="utf-8")
    assert combined.count("micro 1") == 1
    assert combined.count("micro 2") == 1


def test_tsopt_initial_hessian_is_resolved_after_initial_mm_relaxation(
    tmp_path, monkeypatch
) -> None:
    import torch
    import mlmm.workflows.tsopt as tsopt_mod
    from mlmm.core.defaults import RSIRFO_KW

    class _MovingMicro:
        cur_cycle = 0
        is_converged = True
        is_stalled = False
        stop_reason = None
        calls = 0

        def __init__(self, geom, **_kwargs):
            self.geom = geom

        def run(self):
            type(self).calls += 1
            self.geom.coords = np.full_like(
                self.geom.coords, 7.0 if type(self).calls == 1 else 9.0
            )

    _install_common_fakes(monkeypatch, tsopt_mod, _MovingMicro)
    monkeypatch.setattr(
        tsopt_mod,
        "resolve_partition_from_core",
        lambda *a, **k: _micro_active_partition(),
    )
    monkeypatch.setitem(tsopt_mod.TSOPT_CLASS_MAP, "rsirfo", _FakeMacroOptimizer)
    hessian_coords = []

    def _capture_hessian(geom, *_args, **_kwargs):
        hessian_coords.append(np.asarray(geom.coords, dtype=float).copy())
        return torch.zeros((3, 3), dtype=torch.float64)

    monkeypatch.setattr(tsopt_mod, "_calc_full_hessian_torch", _capture_hessian)

    tsopt_mod._run_microiter_tsopt(
        _FakeGeom(n_atoms=2),
        {},
        dict(RSIRFO_KW),
        {},
        {"max_cycles": 1, "dump": False},
        {"micro_max_cycles": 1},
        tmp_path,
        dump=False,
        mode="rsirfo",
    )

    assert len(hessian_coords) == 1
    np.testing.assert_allclose(hessian_coords[0], 7.0)


def test_tsopt_exception_restores_entry_calculator_and_freeze_mask(
    tmp_path,
    monkeypatch,
) -> None:
    import mlmm.workflows.tsopt as tsopt_mod
    from mlmm.core.defaults import RSIRFO_KW

    class _ConvergedMicro:
        cur_cycle = 0
        is_converged = True
        is_stalled = False
        stop_reason = None

        def __init__(self, geom, **_kwargs):
            self.geom = geom

        def run(self):
            return None

    _install_common_fakes(monkeypatch, tsopt_mod, _ConvergedMicro)
    partition = _micro_active_partition()
    monkeypatch.setattr(
        tsopt_mod,
        "resolve_partition_from_core",
        lambda *a, **k: partition,
    )

    def fail_hessian(*args, **kwargs):
        raise RuntimeError("hessian sentinel")

    monkeypatch.setattr(tsopt_mod, "_calc_full_hessian_torch", fail_hessian)
    geometry = _FakeGeom(n_atoms=2)
    entry_calculator = object()
    geometry.calculator = entry_calculator
    geometry._calc = entry_calculator

    with pytest.raises(RuntimeError, match="hessian sentinel"):
        tsopt_mod._run_microiter_tsopt(
            geometry,
            {},
            dict(RSIRFO_KW),
            {},
            {"max_cycles": 1, "dump": False},
            {"micro_max_cycles": 1},
            tmp_path,
            dump=False,
            mode="rsirfo",
        )

    assert geometry.freeze_atoms == list(partition.original_freeze)
    assert geometry._calc is entry_calculator


# ---------------------------------------------------------------------------
# The MM micro relaxation never inherits the energy-plateau stop
# ---------------------------------------------------------------------------


class _CapturingMicro:
    """Records the kwargs the driver hands to the micro LBFGS."""

    captured: List[Dict[str, Any]] = []

    cur_cycle = 0
    is_converged = True
    is_stalled = False
    stop_reason = None

    def __init__(self, _geom, **kwargs):
        type(self).captured.append(dict(kwargs))

    def run(self):
        return None


_PLATEAU_ON = {
    "energy_plateau": True,
    "energy_plateau_thresh": 1.0e-4,
    "energy_plateau_window": 50,
}


def test_opt_micro_never_inherits_the_energy_plateau_stop(tmp_path, monkeypatch):
    """A flat MM energy above the force threshold is a stalled micro, not MM
    equilibrium: stopping there ends the macro/micro alternation with the
    environment unrelaxed. `micro_max_cycles` is the only micro bound, so the
    driver must strip the plateau stop even when every incoming config sets it.
    """
    import mlmm.workflows.opt as opt_mod

    _CapturingMicro.captured = []
    _install_common_fakes(monkeypatch, opt_mod, _CapturingMicro)
    monkeypatch.setattr(opt_mod, "RFOptimizer", _FakeMacroOptimizer)

    opt_mod._run_microiter_opt(
        _FakeGeom(n_atoms=2),
        _FakeCalc(),
        calc_cfg={},
        rfo_cfg={},
        lbfgs_cfg=dict(_PLATEAU_ON),
        opt_cfg={"max_cycles": 1, **_PLATEAU_ON},
        microiter_cfg={"micro_max_cycles": 1},
        out_dir_path=tmp_path,
        partition=_micro_active_partition(),
        dump=False,
    )

    assert _CapturingMicro.captured, "the micro LBFGS was never constructed"
    for kwargs in _CapturingMicro.captured:
        assert kwargs["energy_plateau"] is False


def test_tsopt_micro_never_inherits_the_energy_plateau_stop(tmp_path, monkeypatch):
    """Same contract on the TS driver: a plateau stop in the micro ended the
    whole macro search early and left extra imaginary modes behind.
    """
    import torch
    import mlmm.workflows.tsopt as tsopt_mod
    from mlmm.core.defaults import RSIRFO_KW

    _CapturingMicro.captured = []
    _install_common_fakes(monkeypatch, tsopt_mod, _CapturingMicro)
    monkeypatch.setattr(
        tsopt_mod,
        "_calc_full_hessian_torch",
        lambda *a, **k: torch.zeros((3, 3), dtype=torch.float64),
    )
    monkeypatch.setattr(
        tsopt_mod,
        "resolve_partition_from_core",
        lambda *a, **k: _micro_active_partition(),
    )
    monkeypatch.setitem(tsopt_mod.TSOPT_CLASS_MAP, "rsirfo", _FakeMacroOptimizer)

    tsopt_mod._run_microiter_tsopt(
        _FakeGeom(n_atoms=2),
        {},
        {**RSIRFO_KW, **_PLATEAU_ON},
        dict(_PLATEAU_ON),
        {"max_cycles": 1, "dump": False, **_PLATEAU_ON},
        {"micro_max_cycles": 1},
        tmp_path,
        dump=False,
        thresh=None,
        mode="rsirfo",
        reference_mode=None,
    )

    assert _CapturingMicro.captured, "the micro LBFGS was never constructed"
    for kwargs in _CapturingMicro.captured:
        assert kwargs["energy_plateau"] is False


def test_opt_micro_macro_takes_the_shared_plateau_setting(tmp_path, monkeypatch):
    """The macro step IS the run's optimizer, so an opted-in plateau stop must
    reach it even though the micro step is exempt."""
    import mlmm.workflows.opt as opt_mod

    captured: List[Dict[str, Any]] = []

    class _CapturingMacro(_FakeMacroOptimizer):
        def __init__(self, geom, **kwargs):
            captured.append(dict(kwargs))
            super().__init__(geom, **kwargs)

    _CapturingMicro.captured = []
    _install_common_fakes(monkeypatch, opt_mod, _CapturingMicro)
    monkeypatch.setattr(opt_mod, "RFOptimizer", _CapturingMacro)

    opt_mod._run_microiter_opt(
        _FakeGeom(n_atoms=2),
        _FakeCalc(),
        calc_cfg={},
        rfo_cfg={},
        lbfgs_cfg={},
        opt_cfg={"max_cycles": 1, **_PLATEAU_ON},
        microiter_cfg={"micro_max_cycles": 1},
        out_dir_path=tmp_path,
        partition=_micro_active_partition(),
        dump=False,
    )

    assert captured, "the macro RFO was never constructed"
    assert captured[0]["energy_plateau"] is True
    assert captured[0]["energy_plateau_window"] == 50
