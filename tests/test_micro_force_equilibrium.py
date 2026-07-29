"""A plateaued MM relaxation whose forces are converged IS equilibrium.

The micro stage exists to settle the MM subsystem before the macro step reads
curvature. It inherits the TS threshold, so under the shipped default `baker`
it can plateau with its forces already under 3e-4/2e-4 while the step criteria
are still unmet. A force-converged micro stage is accepted as equilibrium;
every other stall remains fatal.
"""

import pytest

from mlmm.workflows._microiteration import micro_reached_force_equilibrium


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


def test_introduces_no_new_tolerance() -> None:
    """The predicate must read the optimizer's own configured thresholds."""
    import inspect

    src = inspect.getsource(micro_reached_force_equilibrium)
    assert "max_force_thresh" in src and "rms_force_thresh" in src
    # No numeric literal tolerance of its own.
    assert "e-" not in src.replace("1e-", "")


def test_flatten_outcome_is_published_not_only_on_stderr() -> None:
    """A requested flatten that never ran must be visible in result.json.

    Before this, the only signal was one stderr line; the run's own
    result.json, summary.json and summary.log said nothing, so a consumer
    could not tell a flattened result from a vetoed one. That is how I
    initially misread a correct veto as a silently ignored flag.
    """
    from pathlib import Path

    import mlmm.workflows.tsopt as tsopt_mod

    src = Path(tsopt_mod.__file__).read_text(encoding="utf-8")
    assert '"flatten_requested": bool(simple_cfg.get("flatten_max_iter", 0)),' in src
    assert '"flatten_skip_reason": _flatten_skip_reason,' in src
    # Initialised before any branch, so the key is always present.
    assert src.index("_flatten_skip_reason = None") < src.index('"flatten_skip_reason"')


def test_veto_message_distinguishes_undetermined_from_positive() -> None:
    """`None` (never determined) and `False` (determined positive) differ.

    The observed run took zero macro steps, so the sign was never determined,
    yet the message claimed the mode "is not negative" -- a stronger statement
    than the code had evidence for.
    """
    from pathlib import Path

    import mlmm.workflows.tsopt as tsopt_mod

    src = Path(tsopt_mod.__file__).read_text(encoding="utf-8")
    assert "target mode sign never determined" in src
    assert "target mode is not negative" in src
    assert "if target_mode_is_negative is None" in src
