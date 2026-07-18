"""C9 production-bound falsifiers for the shared microiteration contract.

These exercise the real ``mlmm.workflows._microiteration`` helpers that opt and
tsopt both consume:

* M44/M45 — one immutable partition that preserves the user freeze mask in BOTH
  phases and fails LOUDLY (never a swallowed empty set);
* M46 — a fail-closed macro/micro aggregate where a micro plateau/stall/max-cycle
  (or a missing convergence signal) never reads as macro convergence;
* M47 — one field-isomorphic optimizer outcome across ordinary/micro/restart/
  multistart/flatten, with executed (never configured) cycle counts.
"""

from __future__ import annotations

import sys

import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11), reason="mlmm requires Python >= 3.11"
)

from mlmm.workflows._microiteration import (
    MicroiterationOutcome,
    OptimizerOutcome,
    PartitionError,
    build_aggregate,
    build_partition,
    resolve_partition_from_core,
)


class _FakeOptimizer:
    """Minimal pysisyphus-optimizer stand-in for the outcome factory."""

    def __init__(self, *, is_converged, cur_cycle, is_stalled=False, stop_reason=""):
        self.is_converged = is_converged
        self.cur_cycle = cur_cycle
        self.is_stalled = is_stalled
        self.stop_reason = stop_reason


class _FakeCore:
    def __init__(self, ml, hess_mm, movable_mm, frozen_mm, links):
        self.ml_indices = list(ml)
        self.hess_mm_indices = list(hess_mm)
        self.movable_mm_indices = list(movable_mm)
        self.frozen_layer_indices = list(frozen_mm)
        self.mlmm_links = list(links)


# --------------------------------------------------------------------------
# M45 — immutable partition preserves the user freeze mask in BOTH phases
# --------------------------------------------------------------------------


def test_partition_preserves_user_freeze_in_both_phases():
    # atoms: 0,1 = ML; 2 = link parent; 3 = HessMM; 4 = MovableMM; 5 = frozen layer
    # user freezes one ML atom (1), the link parent (2), the HessMM atom (3) and
    # the MovableMM atom (4) -- one of each functional class.
    part = build_partition(
        6,
        ml=[0, 1],
        link_parents=[2],
        hess_mm=[3],
        movable_mm=[4],
        frozen_mm=[5],
        original_freeze=[1, 2, 3, 4],
    )
    for a in (1, 2, 3, 4):
        assert a in part.macro_freeze_atoms, a
        assert a in part.micro_freeze_atoms, a
        assert a not in part.macro_active_atoms, a
        assert a not in part.micro_active_atoms, a
    # The one unfrozen ML atom is macro-active; nothing MM is macro-active.
    assert part.macro_active_atoms == (0,)
    # original_freeze is the exact accepted mask (ascending, de-duplicated).
    assert part.original_freeze == (1, 2, 3, 4)


def test_partition_keeps_movable_mm_in_micro_when_not_frozen():
    # No user freeze: the MovableMM + HessMM atoms are the micro-active set and a
    # link parent is excluded from micro (it co-moves with ML in macro).
    part = build_partition(
        6,
        ml=[0],
        link_parents=[1],
        hess_mm=[2],
        movable_mm=[3, 4],
        frozen_mm=[5],
        original_freeze=[],
    )
    assert part.macro_active_atoms == (0, 1)
    assert part.micro_active_atoms == (2, 3, 4)
    assert 1 not in part.micro_active_atoms  # link parent never micro-active
    assert 5 in part.macro_freeze_atoms and 5 in part.micro_freeze_atoms


def test_partition_phase_masks_are_disjoint_and_cover():
    part = build_partition(
        7, ml=[0, 1], link_parents=[2], hess_mm=[3], movable_mm=[4, 5],
        frozen_mm=[6], original_freeze=[],
    )
    for active, freeze in (
        (part.macro_active_atoms, part.macro_freeze_atoms),
        (part.micro_active_atoms, part.micro_freeze_atoms),
    ):
        assert set(active).isdisjoint(freeze)
        assert sorted([*active, *freeze]) == list(range(7))
    # ordered tuples derived from range(n_atoms) (ascending)
    assert list(part.macro_active_atoms) == sorted(part.macro_active_atoms)
    assert list(part.micro_active_atoms) == sorted(part.micro_active_atoms)


def test_partition_flags_reflect_active_sets():
    empty_ml = build_partition(
        3, ml=[], link_parents=[], hess_mm=[1], movable_mm=[2], frozen_mm=[0],
        original_freeze=[],
    )
    assert empty_ml.has_macro_active is False
    assert empty_ml.has_micro_active is True

    no_micro = build_partition(
        2, ml=[0], link_parents=[], hess_mm=[], movable_mm=[], frozen_mm=[],
        original_freeze=[1],
    )
    assert no_micro.has_macro_active is True
    assert no_micro.has_micro_active is False


# --------------------------------------------------------------------------
# M44 — loud partition failure (never a swallowed empty set)
# --------------------------------------------------------------------------


def test_partition_out_of_range_index_raises_loudly():
    with pytest.raises(PartitionError):
        build_partition(
            3, ml=[0, 9], link_parents=[], hess_mm=[], movable_mm=[],
            frozen_mm=[], original_freeze=[],
        )


def test_resolve_from_none_core_raises_loudly():
    with pytest.raises(PartitionError):
        resolve_partition_from_core(None, 3, [])


def test_resolve_from_core_reads_layers_and_link_parents():
    core = _FakeCore(
        ml=[0, 1], hess_mm=[3], movable_mm=[4], frozen_mm=[5],
        links=[(1, 3)],  # 1-based MM parent atom 3 -> 0-based atom 2
    )
    part = resolve_partition_from_core(core, 6, original_freeze=[])
    assert part.link_parent_atoms == (2,)
    assert part.macro_active_atoms == (0, 1, 2)
    assert part.has_macro_active is True


def test_resolve_from_core_empty_ml_is_valid_fallback_not_error():
    # A genuine empty-ML region is a VALID partition (documented fallback),
    # distinguishable from a swallowed construction error which would raise.
    core = _FakeCore(ml=[], hess_mm=[1], movable_mm=[2], frozen_mm=[0], links=[])
    part = resolve_partition_from_core(core, 3, original_freeze=[])
    assert part.has_macro_active is False


# --------------------------------------------------------------------------
# M46 — fail-closed macro/micro aggregate
# --------------------------------------------------------------------------


def _macro(converged, *, cycles=3, stalled=False, stop_reason=None):
    status = "stalled" if stalled else ("converged" if converged else "not_converged")
    return OptimizerOutcome(
        status=status, executed=True, converged=(False if stalled else converged),
        cycles=cycles, max_cycles=100, stalled=stalled, stop_reason=stop_reason,
    )


def _micro(converged, *, cycles=2, stalled=False, stop_reason=None):
    status = "stalled" if stalled else ("converged" if converged is True else "not_converged")
    return OptimizerOutcome(
        status=status, executed=True, converged=(False if stalled else converged),
        cycles=cycles, max_cycles=50, stalled=stalled, stop_reason=stop_reason,
    )


def test_aggregate_converges_only_when_macro_and_latest_micro_converge():
    agg = build_aggregate(_macro(True), [_micro(True), _micro(True)], max_cycles=100)
    assert agg.converged is True and agg.status == "converged"


def test_aggregate_macro_converged_but_micro_maxcycle_is_not_converged():
    agg = build_aggregate(_macro(True), [_micro(True), _micro(False)], max_cycles=100)
    assert agg.converged is False
    assert agg.status == "not_converged"
    assert agg.stop_reason == "micro_max_cycles"


def test_aggregate_missing_micro_signal_fails_closed():
    agg = build_aggregate(_macro(True), [_micro(None)], max_cycles=100)
    assert agg.converged is False
    assert agg.stop_reason == "micro_convergence_unknown"


def test_aggregate_latest_micro_stall_is_stalled():
    agg = build_aggregate(
        _macro(True), [_micro(True), _micro(False, stalled=True, stop_reason="plateau")],
        max_cycles=100,
    )
    assert agg.stalled is True and agg.status == "stalled"
    assert agg.converged is False


def test_aggregate_macro_stall_is_stalled():
    agg = build_aggregate(_macro(False, stalled=True, stop_reason="macro plateau"), [_micro(True)])
    assert agg.status == "stalled" and agg.converged is False


def test_aggregate_macro_not_converged_reports_macro_reason():
    agg = build_aggregate(_macro(False, stop_reason="max_macro"), [_micro(True)])
    assert agg.converged is False
    assert agg.stop_reason == "max_macro"


def test_aggregate_vacuous_no_micro_follows_macro():
    # No micro-active DOFs -> no micro attempts -> aggregate tracks the macro.
    assert build_aggregate(_macro(True), []).converged is True
    assert build_aggregate(_macro(False), []).converged is False


# --------------------------------------------------------------------------
# M47 — one field-isomorphic outcome, executed cycle counts
# --------------------------------------------------------------------------


def test_from_optimizer_counts_executed_cycles_not_zero_based():
    # A fake optimizer that converges at cur_cycle == 0 reports ONE executed cycle.
    out = OptimizerOutcome.from_optimizer(
        _FakeOptimizer(is_converged=True, cur_cycle=0), max_cycles=10
    )
    assert out.cycles == 1
    assert out.status == "converged" and out.converged is True


def test_from_optimizer_stall_is_not_converged():
    out = OptimizerOutcome.from_optimizer(
        _FakeOptimizer(is_converged=False, cur_cycle=4, is_stalled=True, stop_reason="plateau"),
        max_cycles=10,
    )
    assert out.stalled is True and out.status == "stalled"
    assert out.converged is False


def test_outcome_shape_is_identical_across_paths():
    ordinary = OptimizerOutcome.from_optimizer(_FakeOptimizer(is_converged=True, cur_cycle=2))
    restart = OptimizerOutcome.from_optimizer(_FakeOptimizer(is_converged=False, cur_cycle=7))
    micro = OptimizerOutcome.not_executed()
    vacuous = OptimizerOutcome.vacuous_success()
    keys = set(ordinary.to_dict())
    for other in (restart, micro, vacuous):
        assert set(other.to_dict()) == keys
    assert keys == {
        "status", "executed", "converged", "cycles", "max_cycles", "stalled", "stop_reason",
    }


def test_microiteration_outcome_serializes_separate_macro_and_micro_cycles():
    macro = _macro(True, cycles=2)
    micro_attempts = [_micro(True, cycles=3), _micro(True, cycles=4), _micro(True, cycles=0)]
    agg = build_aggregate(macro, micro_attempts, max_cycles=100)
    part = build_partition(3, ml=[0], link_parents=[], hess_mm=[1], movable_mm=[2],
                           frozen_mm=[], original_freeze=[])
    outcome = MicroiterationOutcome(
        aggregate=agg, macro=macro, micro_attempts=tuple(micro_attempts),
        macro_cycles=2, micro_cycles=7, partition=part,
    )
    obj = outcome.to_result_object()
    assert obj["macro_cycles"] == 2
    assert obj["micro_cycles"] == 7
    assert len(obj["micro_attempts"]) == 3
    # The three ordered child outcomes are retained (not flattened to a boolean).
    assert [m["cycles"] for m in obj["micro_attempts"]] == [3, 4, 0]
    assert obj["aggregate"]["status"] == "converged"
    assert "partition" in obj


def test_microiteration_outcome_retains_stalled_child_not_bare_false():
    # A selected micro-stalled branch keeps its stalled status/reason.
    macro = _macro(False)
    stalled_micro = _micro(False, stalled=True, stop_reason="micro plateau")
    agg = build_aggregate(macro, [stalled_micro], max_cycles=100)
    outcome = MicroiterationOutcome(
        aggregate=agg, macro=macro, micro_attempts=(stalled_micro,),
        macro_cycles=1, micro_cycles=5,
    )
    obj = outcome.to_result_object()
    assert obj["aggregate"]["status"] == "stalled"
    assert obj["micro_attempts"][0]["stalled"] is True
    assert obj["micro_attempts"][0]["stop_reason"] == "micro plateau"
