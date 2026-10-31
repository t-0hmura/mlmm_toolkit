"""A plateaued MM relaxation whose forces are converged IS equilibrium.

The micro stage exists to settle the MM subsystem before the macro step reads
curvature. It inherits the TS threshold, so under the shipped default `baker`
it can plateau with its forces already under 3e-4/2e-4 while the step criteria
are still unmet. A force-converged micro stage is accepted as equilibrium;
every other stall remains fatal.
"""

from mlmm.workflows._microiteration import (
    OptimizerOutcome,
    micro_reached_force_equilibrium,
)


class _Opt:
    def __init__(self, max_force, rms_force, max_thresh=3.0e-4, rms_thresh=2.0e-4):
        self.max_forces = [max_force]
        self.rms_forces = [rms_force]
        self.convergence = {}
        if max_thresh is not None:
            self.convergence["max_force_thresh"] = max_thresh
        if rms_thresh is not None:
            self.convergence["rms_force_thresh"] = rms_thresh


def test_forces_under_baker_thresholds_are_equilibrium() -> None:
    # The observed run: stop_reason named only (max_step, rms_step), i.e. the
    # force criteria were already met.
    assert micro_reached_force_equilibrium(_Opt(2.9e-4, 1.9e-4)) is True


def test_boundary_is_inclusive() -> None:
    assert micro_reached_force_equilibrium(_Opt(3.0e-4, 2.0e-4)) is True


def test_max_force_above_threshold_is_not_equilibrium() -> None:
    assert micro_reached_force_equilibrium(_Opt(3.1e-4, 1.9e-4)) is False


def test_rms_force_above_threshold_is_not_equilibrium() -> None:
    assert micro_reached_force_equilibrium(_Opt(2.9e-4, 2.1e-4)) is False


def test_missing_threshold_is_not_equilibrium() -> None:
    # Fail closed: no configured tolerance means we cannot claim equilibrium.
    assert micro_reached_force_equilibrium(_Opt(1e-9, 1e-9, max_thresh=None)) is False


def test_no_force_history_is_not_equilibrium() -> None:
    opt = _Opt(1e-9, 1e-9)
    opt.max_forces = []
    assert micro_reached_force_equilibrium(opt) is False


def test_rms_absent_falls_back_to_max_force_only() -> None:
    opt = _Opt(2.9e-4, 1.9e-4, rms_thresh=None)
    assert micro_reached_force_equilibrium(opt) is True


def test_verdict_follows_the_optimizers_own_thresholds() -> None:
    """The predicate has no tolerance of its own: the same forces flip with the thresholds."""
    forces = (2.9e-4, 1.9e-4)
    assert micro_reached_force_equilibrium(_Opt(*forces)) is True
    assert micro_reached_force_equilibrium(_Opt(*forces, max_thresh=1.0e-4, rms_thresh=1.0e-4)) is False
    assert micro_reached_force_equilibrium(_Opt(5.0e-4, 4.0e-4, max_thresh=1.0e-3, rms_thresh=1.0e-3)) is True


def test_accepted_force_equilibrium_is_serialized_as_converged() -> None:
    stalled = OptimizerOutcome(
        status="stalled",
        executed=True,
        converged=False,
        cycles=4,
        max_cycles=20,
        stalled=True,
        stop_reason="energy plateau",
    )

    accepted = stalled.accept_force_equilibrium()

    assert accepted.to_dict() == {
        "status": "converged",
        "executed": True,
        "converged": True,
        "cycles": 4,
        "max_cycles": 20,
        "stalled": False,
        "stop_reason": None,
    }
