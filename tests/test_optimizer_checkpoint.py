"""safe-primitive, atomic, explicitly-bounded optimizer checkpoints.

Supported classes (LBFGS, non-TS RF/Hessian optimizers) round-trip through a
YAML-safe envelope and resume the same trajectory; unsupported classes
(TS mode-following) fail loud without writing a partial checkpoint; a failed
write leaves any existing checkpoint byte-identical; and a corrupted checkpoint
is rejected before any optimizer/geometry state is mutated.
"""

from __future__ import annotations

import copy

import numpy as np
import pytest
import yaml

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from pysisyphus.optimizers import checkpoint


class _QuadraticCalculator(Calculator):
    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    @staticmethod
    def _results(coords):
        coords = np.asarray(coords, dtype=float)
        return {"energy": float(coords @ coords), "forces": -2.0 * coords}

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        results = self._results(coords)
        results["hessian"] = 2.0 * np.eye(len(coords))
        return results


def _rfo(out_dir, start, **kwargs):
    out_dir.mkdir(parents=True, exist_ok=True)
    geom = Geometry(["H"], np.array(start, dtype=float), coord_type="cart")
    geom.set_calculator(_QuadraticCalculator(out_dir))
    opt = RFOptimizer(
        geom,
        hessian_init="unit",
        trust_radius=0.1,
        trust_min=1e-4,
        trust_max=0.1,
        line_search=False,
        gdiis=False,
        out_dir=out_dir,
        dump=False,
        thresh="gau_loose",
        **kwargs,
    )
    return geom, opt


def _rfo_trust(out_dir, start, **kwargs):
    """RFO on the exact-Hessian quadratic with head-room for the trust radius
    to grow (``trust_radius`` 0.05 < ``trust_max`` 1.0).

    Every far-from-minimum step is trust-radius-limited, so a resumed run only
    reproduces the uninterrupted trajectory when the *adapted* trust radius is
    restored; a fresh optimizer would fall back to ``min(0.05, 1.0)``.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    geom = Geometry(["H"], np.array(start, dtype=float), coord_type="cart")
    geom.set_calculator(_QuadraticCalculator(out_dir))
    opt = RFOptimizer(
        geom,
        hessian_init="calc",
        trust_radius=0.05,
        trust_min=1e-4,
        trust_max=1.0,
        line_search=False,
        gdiis=False,
        out_dir=out_dir,
        dump=False,
        thresh="gau_loose",
        **kwargs,
    )
    return geom, opt


def _lbfgs(out_dir, start, **kwargs):
    out_dir.mkdir(parents=True, exist_ok=True)
    geom = Geometry(["H"], np.array(start, dtype=float), coord_type="cart")
    geom.set_calculator(_QuadraticCalculator(out_dir))
    opt = LBFGS(geom, out_dir=out_dir, dump=False, thresh="gau_loose", **kwargs)
    return geom, opt


def _tsopt(out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    geom = Geometry(["H"], np.array([0.5, 0.1, 0.0], dtype=float), coord_type="cart")
    geom.set_calculator(_QuadraticCalculator(out_dir))
    return geom, RSPRFOptimizer(geom, out_dir=out_dir, dump=False)


def _only_primitives(obj) -> bool:
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return True
    if isinstance(obj, dict):
        return all(isinstance(k, str) and _only_primitives(v) for k, v in obj.items())
    if isinstance(obj, list):
        return all(_only_primitives(v) for v in obj)
    return False


def test_supported_optimizer_dump_is_yaml_safe_and_atomic(tmp_path) -> None:
    _, opt = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt.run()

    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt, ck)

    loaded = yaml.safe_load(ck.read_text())
    assert _only_primitives(loaded)
    assert loaded["schema"] == checkpoint.CHECKPOINT_SCHEMA
    assert loaded["phase"] == checkpoint.RESUMABLE_PHASE
    # No temporary sibling remains after a successful atomic replace.
    assert not list(tmp_path.glob("restart.yaml.*.tmp"))


def test_supported_rfo_round_trips_and_resumes_uninterrupted_trace(tmp_path) -> None:
    geom_full, opt_full = _rfo(tmp_path / "full", [0.5, 0.0, 0.0], max_cycles=200)
    opt_full.run()

    geom_a, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)

    geom_b, opt_b = _rfo(tmp_path / "b", [9.9, 9.9, 9.9], max_cycles=200)
    checkpoint.load_and_apply(opt_b, ck)

    # Restored accepted-state histories match the interrupted optimizer exactly.
    np.testing.assert_allclose(opt_b.coords[-1], opt_a.coords[-1])
    np.testing.assert_allclose(opt_b.H, opt_a.H)
    np.testing.assert_allclose(opt_b.energies, opt_a.energies)

    opt_b.run()
    # The resumed optimizer converges to the same minimum as the uninterrupted run.
    np.testing.assert_allclose(geom_b.cart_coords, geom_full.cart_coords, atol=1e-7)


def test_supported_lbfgs_round_trips(tmp_path) -> None:
    geom_a, opt_a = _lbfgs(tmp_path / "la", [0.5, 0.0, 0.0], max_cycles=3)
    opt_a.run()
    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)

    assert _only_primitives(yaml.safe_load(ck.read_text()))

    geom_b, opt_b = _lbfgs(tmp_path / "lb", [9.9, 9.9, 9.9], max_cycles=200)
    checkpoint.load_and_apply(opt_b, ck)
    np.testing.assert_allclose(opt_b.coords[-1], opt_a.coords[-1])
    assert len(opt_b.coord_diffs) == len(opt_a.coord_diffs)
    for restored, original in zip(opt_b.coord_diffs, opt_a.coord_diffs):
        np.testing.assert_allclose(restored, original)


def test_unsupported_optimizer_fails_loud_without_writing(tmp_path) -> None:
    _, opt = _tsopt(tmp_path / "ts")
    ck = tmp_path / "restart.yaml"
    with pytest.raises(checkpoint.CheckpointUnsupportedError):
        checkpoint.save_checkpoint(opt, ck)
    assert not ck.exists()


def test_serialization_rejects_unsafe_state() -> None:
    with pytest.raises(checkpoint.CheckpointSerializationError):
        checkpoint.to_safe_primitives({"factory": lambda: 1})
    with pytest.raises(checkpoint.CheckpointSerializationError):
        checkpoint.to_safe_primitives(np.array([object()], dtype=object))


def test_write_failure_leaves_existing_checkpoint_byte_identical(tmp_path, monkeypatch) -> None:
    _, opt = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt.run()
    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt, ck)
    original = ck.read_bytes()

    def _boom(src, dst):
        raise OSError("simulated replace failure")

    monkeypatch.setattr(checkpoint.os, "replace", _boom)
    with pytest.raises(OSError):
        checkpoint.save_checkpoint(opt, ck)

    assert ck.read_bytes() == original
    assert not list(tmp_path.glob("restart.yaml.*.tmp"))


def test_corrupt_checkpoint_is_rejected_before_any_mutation(tmp_path) -> None:
    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)
    payload = checkpoint.load_payload(ck)

    _, opt_b = _rfo(tmp_path / "b", [9.9, 9.9, 9.9], max_cycles=200)

    # A clean payload validates against the sibling geometry (identity fields
    # match; only coordinates differ, which are not part of the identity).
    assert checkpoint.validate_payload(payload, opt_b)

    corruptions = []
    bad = dict(payload); bad["schema"] = "not-a-checkpoint"; corruptions.append(bad)
    bad = dict(payload); bad["version"] = 999; corruptions.append(bad)
    bad = dict(payload); bad["phase"] = "terminal"; corruptions.append(bad)
    bad = dict(payload); bad["optimizer_id"] = "some.other.Optimizer"; corruptions.append(bad)
    bad = copy.deepcopy(payload); bad["geometry_identity"]["atoms"] = ["He"]; corruptions.append(bad)
    bad = copy.deepcopy(payload)
    bad["restart_info"]["forces"] = bad["restart_info"]["forces"][:-1]
    corruptions.append(bad)
    bad = copy.deepcopy(payload)
    del bad["restart_info"]["coords"]
    corruptions.append(bad)

    for corrupt in corruptions:
        with pytest.raises(checkpoint.CheckpointValidationError):
            checkpoint.validate_payload(corrupt, opt_b)

    # A legacy unversioned/unsafe payload is rejected, never partially applied.
    with pytest.raises(checkpoint.CheckpointValidationError):
        checkpoint.validate_payload({"coords": [[0.0, 0.0, 0.0]]}, opt_b)


def test_unsupported_target_rejected_on_load(tmp_path) -> None:
    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)
    payload = checkpoint.load_payload(ck)

    _, ts_opt = _tsopt(tmp_path / "ts")
    with pytest.raises(checkpoint.CheckpointUnsupportedError):
        checkpoint.validate_payload(payload, ts_opt)


def test_supported_rfo_round_trips_adapted_trust_radius(tmp_path) -> None:
    """An RF/Hessian optimizer whose per-cycle ``trust_radius``
    was adapted away from its default round-trips that exact radius, and the
    resumed run reproduces the uninterrupted trajectory.
    """
    geom_full, opt_full = _rfo_trust(tmp_path / "full", [9.9, 9.9, 9.9], max_cycles=200)
    opt_full.run()
    assert opt_full.is_converged

    # Interrupt once the trust radius has grown past its 0.05 initial value.
    geom_a, opt_a = _rfo_trust(tmp_path / "a", [9.9, 9.9, 9.9], max_cycles=4)
    opt_a.run()
    assert opt_a.trust_radius != pytest.approx(0.05)  # genuinely adapted

    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)

    geom_b, opt_b = _rfo_trust(tmp_path / "b", [0.0, 0.0, 0.0], max_cycles=200)
    # A fresh optimizer initialises trust_radius to min(0.05, 1.0) = 0.05.
    assert opt_b.trust_radius == pytest.approx(0.05)
    checkpoint.load_and_apply(opt_b, ck)
    # The exact adapted trust radius round-trips through the checkpoint.
    assert opt_b.trust_radius == pytest.approx(opt_a.trust_radius)

    opt_b.run()
    # Resume-equivalence: same cycle count and cycle-by-cycle trajectory as the
    # uninterrupted reference; only the restored trust radius makes this hold.
    assert opt_b.cur_cycle == opt_full.cur_cycle
    assert len(opt_b.coords) == len(opt_full.coords)
    for resumed, reference in zip(opt_b.coords, opt_full.coords):
        np.testing.assert_allclose(resumed, reference, atol=1e-9)
    np.testing.assert_allclose(geom_b.cart_coords, geom_full.cart_coords, atol=1e-7)


def test_supported_lbfgs_round_trips_uphill_rejection_state(tmp_path) -> None:
    """An LBFGS interrupted at an uphill-rejection boundary
    preserves the adaptive state (``_trial_max_step`` and the rejection
    counters), so the resumed run continues the uninterrupted trajectory.

    ``beta > 1`` overshoots the quadratic minimum, firing the uphill-rejection
    safeguard.
    """
    geom_full, opt_full = _lbfgs(
        tmp_path / "lf", [0.3, 0.0, 0.0],
        beta=2.0, max_step=2.0, double_damp=False, max_cycles=60,
        reject_uphill=True,
    )
    opt_full.run()
    assert opt_full.is_converged

    geom_a, opt_a = _lbfgs(
        tmp_path / "la", [0.3, 0.0, 0.0],
        beta=2.0, max_step=2.0, double_damp=False, max_cycles=3,
        reject_uphill=True,
    )
    opt_a.run()
    # The scenario actually triggers a rejection and drives the adaptive state
    # away from the __init__ defaults (_trial_max_step=2.0, counters=0).
    assert opt_a.rejected_uphill_steps >= 1
    assert opt_a._trial_max_step != pytest.approx(2.0)

    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)

    geom_b, opt_b = _lbfgs(
        tmp_path / "lb", [0.0, 0.0, 0.0],
        beta=2.0, max_step=2.0, double_damp=False, max_cycles=60,
        reject_uphill=True,
    )
    # Fresh optimizer carries only the __init__ defaults.
    assert opt_b._trial_max_step == pytest.approx(2.0)
    assert opt_b.rejected_uphill_steps == 0
    assert opt_b.rejections_at_floor == 0

    checkpoint.load_and_apply(opt_b, ck)
    # The uphill-rejection adaptive state round-trips exactly.
    assert opt_b._trial_max_step == pytest.approx(opt_a._trial_max_step)
    assert opt_b.rejected_uphill_steps == opt_a.rejected_uphill_steps
    assert opt_b.rejections_at_floor == opt_a.rejections_at_floor

    opt_b.run()
    # The resumed run still converges to the same minimum as the reference.
    assert opt_b.is_converged
    np.testing.assert_allclose(geom_b.cart_coords, geom_full.cart_coords, atol=1e-7)


def test_missing_subclass_key_rejected_before_any_base_mutation(tmp_path) -> None:
    """A checkpoint missing a subclass-required restart key (H for Hessian/RF,
    coord_diffs for LBFGS) is rejected with the TYPED validation error and
    leaves the optimizer's base history untouched.

    Validation must complete before ``set_restart_info`` can overwrite
    coordinates, energies, forces, or steps.
    """
    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)
    payload = checkpoint.load_payload(ck)

    bad = copy.deepcopy(payload)
    del bad["restart_info"]["H"]

    _, opt_b = _rfo(tmp_path / "b", [9.9, 9.9, 9.9], max_cycles=200)
    # Typed rejection, not a bare KeyError.
    with pytest.raises(checkpoint.CheckpointValidationError):
        checkpoint.validate_payload(bad, opt_b)

    # A full load of the same corrupt checkpoint leaves the base history empty
    # (nothing was partially applied before the error).
    assert opt_b.coords == [] and opt_b.energies == []
    assert opt_b.forces == [] and opt_b.steps == []
    bad_ck = tmp_path / "bad.yaml"
    bad_ck.write_text(yaml.safe_dump(bad, default_flow_style=False, sort_keys=True))
    with pytest.raises(checkpoint.CheckpointValidationError):
        checkpoint.load_and_apply(opt_b, bad_ck)
    assert opt_b.coords == [] and opt_b.energies == []
    assert opt_b.forces == [] and opt_b.steps == []

    # The same guarantee holds for an LBFGS-specific missing key.
    _, opt_la = _lbfgs(tmp_path / "la", [0.5, 0.0, 0.0], max_cycles=3)
    opt_la.run()
    lck = tmp_path / "lbfgs.yaml"
    checkpoint.save_checkpoint(opt_la, lck)
    lbad = copy.deepcopy(checkpoint.load_payload(lck))
    del lbad["restart_info"]["coord_diffs"]

    _, opt_lb = _lbfgs(tmp_path / "lb", [9.9, 9.9, 9.9], max_cycles=200)
    with pytest.raises(checkpoint.CheckpointValidationError):
        checkpoint.validate_payload(lbad, opt_lb)


def test_backward_tolerant_load_of_pre_adaptive_checkpoint(tmp_path) -> None:
    """A checkpoint written before the adaptive keys existed (no trust_radius /
    no _trial_max_step / no rejection counters) still loads: the newly added
    keys are restored tolerantly, keeping the optimizer at its __init__ state.
    """
    # RFO: drop the newly added trust_radius key from an otherwise valid payload.
    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)
    payload = checkpoint.load_payload(ck)
    legacy = copy.deepcopy(payload)
    legacy["restart_info"].pop("trust_radius", None)
    legacy_ck = tmp_path / "legacy.yaml"
    legacy_ck.write_text(yaml.safe_dump(legacy, default_flow_style=False, sort_keys=True))

    _, opt_b = _rfo(tmp_path / "b", [9.9, 9.9, 9.9], max_cycles=200)
    init_trust = opt_b.trust_radius
    checkpoint.load_and_apply(opt_b, legacy_ck)  # must not raise
    # Absent key -> the __init__ trust radius is retained.
    assert opt_b.trust_radius == pytest.approx(init_trust)

    # LBFGS: drop the adaptive-state keys.
    _, opt_la = _lbfgs(tmp_path / "la", [0.5, 0.0, 0.0], max_cycles=3)
    opt_la.run()
    lck = tmp_path / "lbfgs.yaml"
    checkpoint.save_checkpoint(opt_la, lck)
    llegacy = copy.deepcopy(checkpoint.load_payload(lck))
    for key in ("_trial_max_step", "rejected_uphill_steps", "rejections_at_floor"):
        llegacy["restart_info"].pop(key, None)
    llegacy_ck = tmp_path / "lbfgs_legacy.yaml"
    llegacy_ck.write_text(yaml.safe_dump(llegacy, default_flow_style=False, sort_keys=True))

    _, opt_lb = _lbfgs(tmp_path / "lb", [9.9, 9.9, 9.9], max_cycles=200)
    checkpoint.load_and_apply(opt_lb, llegacy_ck)  # must not raise
    assert opt_lb._trial_max_step == pytest.approx(float(opt_lb.max_step))
    assert opt_lb.rejected_uphill_steps == 0
    assert opt_lb.rejections_at_floor == 0


# ---------------------------------------------------------------------------
# Complete HessianOptimizer/RFO adaptive state and base cart_coords.
#
# Most tests below verify resume equivalence against uninterrupted execution.
# engine (the attribute is not serialized, so a resumed optimizer silently
# reverts to its __init__ default) and passes once the attribute round-trips.
# RFOptimizer defines no restart override, so exercising it here also covers
# HessianOptimizer._get/_set_opt_restart_info (RFO inherits them unchanged).
# ---------------------------------------------------------------------------


class _AnisotropicCalculator(Calculator):
    """Convex quadratic ``E = 0.5 * x^T diag(k) x`` with distinct curvatures.

    Successive far-from-minimum RFO steps built on a *unit* model Hessian point
    in different directions, so the multi-step Hessian buffer holds genuinely
    non-parallel columns and its restoration measurably changes the next update
    (an isotropic quadratic keeps every step radial and hides the effect).
    """

    _k = np.array([4.0, 1.0, 0.25])

    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    @classmethod
    def _results(cls, coords):
        coords = np.asarray(coords, dtype=float)
        k = cls._k
        return {"energy": float(0.5 * np.sum(k * coords**2)), "forces": -k * coords}

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        results = self._results(coords)
        results["hessian"] = np.diag(self._k)
        return results


def _rfo_aniso(out_dir, start, **kwargs):
    """RFO on the anisotropic quadratic with a *unit* initial Hessian, so the
    quasi-Newton Hessian update (and its multi-step buffer) matters."""
    out_dir.mkdir(parents=True, exist_ok=True)
    geom = Geometry(["H"], np.array(start, dtype=float), coord_type="cart")
    geom.set_calculator(_AnisotropicCalculator(out_dir))
    opt = RFOptimizer(
        geom,
        hessian_init="unit",
        trust_radius=0.1,
        trust_min=1e-4,
        trust_max=0.1,
        line_search=False,
        gdiis=False,
        out_dir=out_dir,
        dump=False,
        thresh="gau_loose",
        **kwargs,
    )
    return geom, opt


def test_hessian_rejection_counters_round_trip_and_are_backward_tolerant(tmp_path) -> None:
    """The uphill-rejection counters
    ``rejections_at_floor`` / ``rejected_uphill_steps`` / the terminal
    ``uphill_rejection_stalled`` flag round-trip through the checkpoint, and a
    checkpoint written before they existed still loads (defaults retained).
    """
    geom_a, opt_a = _rfo(
        tmp_path / "a", [0.5, 0.3, 0.1], max_cycles=3,
        reject_uphill=True, max_rejections_at_floor=3,
    )
    opt_a.run()
    # Drive the rejection state away from the __init__ defaults.
    opt_a.rejections_at_floor = 2
    opt_a.rejected_uphill_steps = 5
    opt_a.uphill_rejection_stalled = True

    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)

    geom_b, opt_b = _rfo(
        tmp_path / "b", [9.9, 9.9, 9.9], max_cycles=200,
        reject_uphill=True, max_rejections_at_floor=3,
    )
    # A fresh optimizer carries only the __init__ defaults.
    assert opt_b.rejections_at_floor == 0
    assert opt_b.rejected_uphill_steps == 0
    assert opt_b.uphill_rejection_stalled is False

    checkpoint.load_and_apply(opt_b, ck)
    # The complete rejection state round-trips exactly.
    assert opt_b.rejections_at_floor == 2
    assert opt_b.rejected_uphill_steps == 5
    assert opt_b.uphill_rejection_stalled is True

    # Backward tolerance: a checkpoint written before these keys existed loads
    # and keeps the __init__ defaults (the keys are presence-guarded, not
    # required for validation).
    legacy = copy.deepcopy(checkpoint.load_payload(ck))
    for key in ("rejections_at_floor", "rejected_uphill_steps", "uphill_rejection_stalled"):
        legacy["restart_info"].pop(key, None)
    legacy_ck = tmp_path / "legacy.yaml"
    legacy_ck.write_text(yaml.safe_dump(legacy, default_flow_style=False, sort_keys=True))

    geom_c, opt_c = _rfo(
        tmp_path / "c", [9.9, 9.9, 9.9], max_cycles=200,
        reject_uphill=True, max_rejections_at_floor=3,
    )
    checkpoint.load_and_apply(opt_c, legacy_ck)  # must not raise
    assert opt_c.rejections_at_floor == 0
    assert opt_c.rejected_uphill_steps == 0
    assert opt_c.uphill_rejection_stalled is False


def test_supported_rfo_resume_terminates_at_uphill_rejection_floor_boundary(tmp_path) -> None:
    """``rejections_at_floor`` gates
    ``request_stop`` (HessianOptimizer.reject_current_uphill_step:
    ``rejections_at_floor >= max_rejections_at_floor`` -> stop).  A run
    interrupted one rejection short of that boundary, while pinned at the
    emergency trust floor, must terminate on the *next* at-floor rejection after
    resume -- exactly as the uninterrupted run would.

    """
    geom_a, opt_a = _rfo(
        tmp_path / "a", [0.5, 0.3, 0.1], max_cycles=4,
        reject_uphill=True, max_rejections_at_floor=3,
    )
    opt_a.run()
    # State of a run stuck at the emergency floor, one rejection short of stop.
    opt_a.rejections_at_floor = opt_a.max_rejections_at_floor - 1  # == 2
    opt_a.trust_radius = opt_a.rejection_trust_floor  # at the floor

    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)

    geom_b, opt_b = _rfo(
        tmp_path / "b", [9.9, 9.9, 9.9], max_cycles=200,
        reject_uphill=True, max_rejections_at_floor=3,
    )
    checkpoint.load_and_apply(opt_b, ck)
    assert opt_b.rejections_at_floor == opt_a.max_rejections_at_floor - 1
    assert opt_b.trust_radius == pytest.approx(opt_a.rejection_trust_floor)

    # A single at-floor uphill rejection reaches the stop boundary and
    # terminates the run.  (cart_coords[-2], also restored by this fix, lets the
    # rejection transaction restore the previous geometry.)
    opt_b.reject_current_uphill_step()
    assert opt_b.rejections_at_floor == opt_b.max_rejections_at_floor
    assert opt_b.uphill_rejection_stalled is True
    assert opt_b.stop_requested is True


def test_supported_rfo_round_trips_multistep_hessian_buffer(tmp_path) -> None:
    """With ``hessian_update_window >= 2`` the sliding (dx, dg)
    buffer ``_sy_buffer_S`` / ``_sy_buffer_Y`` feeds the next multi-step
    TS-BFGS Hessian update.  It round-trips through the checkpoint, and the
    resumed run reproduces the uninterrupted trajectory.

    """
    geom_full, opt_full = _rfo_aniso(
        tmp_path / "full", [9.9, 9.9, 9.9], max_cycles=30, hessian_update_window=2,
    )
    opt_full.run()

    geom_a, opt_a = _rfo_aniso(
        tmp_path / "a", [9.9, 9.9, 9.9], max_cycles=4, hessian_update_window=2,
    )
    opt_a.run()
    # The scenario actually filled the window-2 buffer.
    assert len(opt_a._sy_buffer_S) == 2 and len(opt_a._sy_buffer_Y) == 2

    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)

    geom_b, opt_b = _rfo_aniso(
        tmp_path / "b", [0.0, 0.0, 0.0], max_cycles=30, hessian_update_window=2,
    )
    # A fresh optimizer starts with an empty buffer.
    assert opt_b._sy_buffer_S == [] and opt_b._sy_buffer_Y == []
    checkpoint.load_and_apply(opt_b, ck)

    # The buffer round-trips exactly.
    assert len(opt_b._sy_buffer_S) == len(opt_a._sy_buffer_S) == 2
    for restored, original in zip(opt_b._sy_buffer_S, opt_a._sy_buffer_S):
        np.testing.assert_allclose(restored, original)
    for restored, original in zip(opt_b._sy_buffer_Y, opt_a._sy_buffer_Y):
        np.testing.assert_allclose(restored, original)

    opt_b.run()
    # Resume-equivalence.  The resume geometry (index k) is bit-identical, and
    # the first step *fed by the restored buffer* (index k+1) matches the
    # uninterrupted run to floating-point precision; an empty buffer forms a
    # rank-1 instead of rank-2 TS-BFGS update here and this step diverges.
    k = len(opt_a.coords)  # first index (re)computed after resume
    assert len(opt_b.coords) == len(opt_full.coords)
    np.testing.assert_allclose(opt_b.coords[k], opt_full.coords[k], atol=1e-12)
    np.testing.assert_allclose(opt_b.coords[k + 1], opt_full.coords[k + 1], atol=1e-7)
    # Same endpoint as the uninterrupted run.  The tolerance absorbs the FP-level
    # path difference (~1e-11 at the resume cycle) amplified over the remaining
    # cycles by the deliberately ill-conditioned (16:1) anisotropic surface;
    # two identical fresh runs are bit-identical, so this is path noise, not a
    # missing-state divergence.
    np.testing.assert_allclose(geom_b.cart_coords, geom_full.cart_coords, atol=1e-5)


def test_supported_rfo_round_trips_prev_eigvec_min_with_overlaps(tmp_path) -> None:
    """With ``rfo_overlaps=True`` the previous minimum-mode
    eigenvector ``_prev_eigvec_min`` selects the RFO root by overlap.  It
    round-trips through the checkpoint (``None`` -> the stored vector), and the
    resumed run continues to the same minimum.

    """
    geom_a, opt_a = _rfo_aniso(
        tmp_path / "a", [9.9, 9.9, 9.9], max_cycles=4, rfo_overlaps=True,
    )
    opt_a.run()
    # Overlap-based mode following actually stored a previous eigenvector.
    assert opt_a._prev_eigvec_min is not None

    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)

    geom_b, opt_b = _rfo_aniso(
        tmp_path / "b", [0.0, 0.0, 0.0], max_cycles=30, rfo_overlaps=True,
    )
    # A fresh optimizer has no previous eigenvector.
    assert opt_b._prev_eigvec_min is None
    checkpoint.load_and_apply(opt_b, ck)

    # The previous minimum-mode eigenvector round-trips exactly.
    assert opt_b._prev_eigvec_min is not None
    np.testing.assert_allclose(opt_b._prev_eigvec_min, opt_a._prev_eigvec_min)

    geom_full, opt_full = _rfo_aniso(
        tmp_path / "full", [9.9, 9.9, 9.9], max_cycles=30, rfo_overlaps=True,
    )
    opt_full.run()
    opt_b.run()
    # Resume-equivalence: with the restored root-following state the resumed run
    # reaches the same endpoint as the uninterrupted reference (the strict
    # discriminator here is the _prev_eigvec_min round-trip above; this endpoint
    # check confirms the resumed run tracks the reference, with a tolerance that
    # absorbs the FP-level path difference amplified on the ill-conditioned
    # anisotropic surface).
    assert len(opt_b.coords) == len(opt_full.coords)
    np.testing.assert_allclose(geom_b.cart_coords, geom_full.cart_coords, atol=1e-5)


def test_supported_rfo_resume_restores_cart_coords_for_uphill_rejection(tmp_path) -> None:
    """The base ``cart_coords`` history round-trips, so an uphill
    rejection *after* resume can restore the previous geometry.

    ``reject_current_trial`` restores ``cart_coords[-2]`` and requires at least
    two history entries.
    """
    geom_a, opt_a = _rfo(tmp_path / "a", [0.5, 0.3, 0.1], max_cycles=3)
    opt_a.run()
    assert len(opt_a.cart_coords) >= 2

    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)

    geom_b, opt_b = _rfo(tmp_path / "b", [9.9, 9.9, 9.9], max_cycles=200)
    # A fresh optimizer has no cart_coords history yet.
    assert opt_b.cart_coords == []
    checkpoint.load_and_apply(opt_b, ck)

    # cart_coords round-trips, aligned with the coords history (>= 2 entries).
    assert len(opt_b.cart_coords) == len(opt_b.coords) >= 2
    for restored, original in zip(opt_b.cart_coords, opt_a.cart_coords):
        np.testing.assert_allclose(restored, original)

    # An uphill rejection immediately after resume succeeds: it restores the
    # previous geometry via cart_coords[-2] instead of raising the len<2 guard.
    prev_cart = opt_b.cart_coords[-2].copy()
    n_before = len(opt_b.coords)
    opt_b.reject_current_trial()
    np.testing.assert_allclose(geom_b.cart_coords, prev_cart)
    assert len(opt_b.coords) == n_before - 1


def test_base_get_restart_info_is_backward_tolerant_without_cart_coords(tmp_path) -> None:
    """A checkpoint written before ``cart_coords`` was serialized still loads;
    the presence-guard leaves the fresh optimizer's empty cart_coords in place.
    """
    geom_a, opt_a = _rfo(tmp_path / "a", [0.5, 0.3, 0.1], max_cycles=3)
    opt_a.run()
    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)

    legacy = copy.deepcopy(checkpoint.load_payload(ck))
    legacy["restart_info"].pop("cart_coords", None)
    legacy_ck = tmp_path / "legacy.yaml"
    legacy_ck.write_text(yaml.safe_dump(legacy, default_flow_style=False, sort_keys=True))

    geom_b, opt_b = _rfo(tmp_path / "b", [9.9, 9.9, 9.9], max_cycles=200)
    checkpoint.load_and_apply(opt_b, legacy_ck)  # must not raise
    assert opt_b.cart_coords == []
    # The accepted-state coords history is still restored normally.
    assert len(opt_b.coords) >= 2


def test_legacy_checkpoint_without_cart_coords_rejects_uphill_rollback(
    tmp_path,
) -> None:
    _, opt_a = _lbfgs(
        tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=3, reject_uphill=False
    )
    opt_a.run()
    checkpoint_path = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, checkpoint_path)
    legacy = copy.deepcopy(checkpoint.load_payload(checkpoint_path))
    legacy["restart_info"].pop("cart_coords", None)
    legacy_path = tmp_path / "legacy.yaml"
    legacy_path.write_text(
        yaml.safe_dump(legacy, default_flow_style=False, sort_keys=True)
    )

    _, tolerant = _lbfgs(
        tmp_path / "tolerant",
        [9.9, 9.9, 9.9],
        max_cycles=200,
        reject_uphill=False,
    )
    checkpoint.load_and_apply(tolerant, legacy_path)
    assert tolerant.cart_coords == []

    geom_strict, strict = _lbfgs(
        tmp_path / "strict",
        [9.9, 9.9, 9.9],
        max_cycles=200,
        reject_uphill=True,
    )
    original_coords = geom_strict.cart_coords.copy()
    with pytest.raises(
        checkpoint.CheckpointValidationError,
        match="Cartesian history required for uphill-trial rollback",
    ):
        checkpoint.load_and_apply(strict, legacy_path)
    assert strict.coords == []
    np.testing.assert_allclose(geom_strict.cart_coords, original_coords)


def test_supported_lbfgs_round_trips_mu_reg(tmp_path) -> None:
    # Regularized-L-BFGS adaptive state: mu_reg is adapted per cycle and feeds
    # get_lbfgs_step, so it must round-trip through the checkpoint for a resumed
    # regularized run to preserve its trajectory.
    _, opt_a = _lbfgs(tmp_path / "mra", [0.5, 0.0, 0.0], max_cycles=3, mu_reg=0.1)
    opt_a.mu_reg = 0.037
    opt_a.tot_adapt_mu_cycles = 4
    info = opt_a._get_opt_restart_info()
    assert info["mu_reg"] == 0.037
    assert info["tot_adapt_mu_cycles"] == 4

    _, opt_b = _lbfgs(tmp_path / "mrb", [9.9, 9.9, 9.9], max_cycles=200, mu_reg=0.1)
    opt_b._set_opt_restart_info(info)
    assert opt_b.mu_reg == 0.037
    assert opt_b.tot_adapt_mu_cycles == 4


def test_lbfgs_mu_reg_backward_tolerant_without_key(tmp_path) -> None:
    _, opt_a = _lbfgs(tmp_path / "bta", [0.5, 0.0, 0.0], max_cycles=3, mu_reg=0.1)
    info = opt_a._get_opt_restart_info()
    info.pop("mu_reg", None)
    info.pop("tot_adapt_mu_cycles", None)
    _, opt_b = _lbfgs(tmp_path / "btb", [1.0, 0.0, 0.0], max_cycles=3, mu_reg=0.2)
    opt_b._set_opt_restart_info(info)
    assert opt_b.mu_reg == 0.2


def test_corrupt_nested_geom_info_is_rejected_before_any_mutation(tmp_path) -> None:
    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)
    payload = checkpoint.load_payload(ck)

    _, opt_b = _rfo(tmp_path / "b", [9.9, 9.9, 9.9], max_cycles=200)
    reference_energies = list(opt_b.energies)

    corruptions = []
    bad = copy.deepcopy(payload)
    bad["restart_info"]["geom_info"]["atoms"] = ["He"]
    corruptions.append(bad)
    bad = copy.deepcopy(payload)
    bad["restart_info"]["geom_info"]["coord_type"] = "redund"
    corruptions.append(bad)
    bad = copy.deepcopy(payload)
    bad["restart_info"]["geom_info"]["cart_coords"] = [0.0, 0.0]
    corruptions.append(bad)
    bad = copy.deepcopy(payload)
    bad["restart_info"]["geom_info"]["cart_coords"] = [0.0, float("nan"), 0.0]
    corruptions.append(bad)

    for corrupt in corruptions:
        with pytest.raises(checkpoint.CheckpointValidationError):
            checkpoint.validate_payload(corrupt, opt_b)
    assert list(opt_b.energies) == reference_energies


def test_corrupt_restart_hessian_is_rejected_before_any_mutation(tmp_path) -> None:
    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)
    payload = checkpoint.load_payload(ck)

    _, opt_b = _rfo(tmp_path / "b", [9.9, 9.9, 9.9], max_cycles=200)
    # No Hessian is set before prepare_opt, so a rejected payload must leave the
    # target without one.
    assert getattr(opt_b, "H", None) is None

    corruptions = []
    bad = copy.deepcopy(payload)
    bad["restart_info"]["H"] = [[1.0, 0.0], [0.0, 1.0]]
    corruptions.append(bad)
    bad = copy.deepcopy(payload)
    bad["restart_info"]["H"] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    corruptions.append(bad)
    bad = copy.deepcopy(payload)
    bad["restart_info"]["H"][0][0] = float("inf")
    corruptions.append(bad)

    for corrupt in corruptions:
        with pytest.raises(checkpoint.CheckpointValidationError):
            checkpoint.validate_payload(corrupt, opt_b)
    assert getattr(opt_b, "H", None) is None


def test_torch_hessian_round_trips_while_multistep_buffers_stay_numpy(tmp_path) -> None:
    import torch

    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    opt_a.H = torch.as_tensor(np.array(opt_a.H), dtype=torch.float64)
    # The multi-step update stacks these with np.column_stack, so they are and
    # must stay NumPy-only state even when H lives on another backend.
    opt_a._sy_buffer_S = [np.ones(3, dtype=np.float64)]
    opt_a._sy_buffer_Y = [np.full(3, 2.0, dtype=np.float64)]
    opt_a._prev_eigvec_min = np.zeros(3, dtype=np.float64)

    info = opt_a._get_opt_restart_info()
    assert _only_primitives(info)
    assert info["H_spec"]["backend"] == "torch"
    # Every state carries its own spec; none of them inherits H's backend.
    assert info["_sy_buffer_S_spec"]["backend"] == "numpy"
    assert info["_sy_buffer_Y_spec"]["backend"] == "numpy"
    assert info["_prev_eigvec_min_spec"]["backend"] == "numpy"

    _, opt_b = _rfo(tmp_path / "b", [0.5, 0.0, 0.0], max_cycles=2)
    opt_b._set_opt_restart_info(info)

    assert isinstance(opt_b.H, torch.Tensor)
    assert opt_b.H.dtype == torch.float64
    assert isinstance(opt_b._sy_buffer_S[0], np.ndarray)
    assert isinstance(opt_b._sy_buffer_Y[0], np.ndarray)
    assert isinstance(opt_b._prev_eigvec_min, np.ndarray)
    # The NumPy-only consumer of the buffers still works after the resume.
    assert np.column_stack(opt_b._sy_buffer_S).shape == (3, 1)
    assert np.column_stack(opt_b._sy_buffer_Y).shape == (3, 1)


def test_overlap_vector_keeps_its_own_backend_spec(tmp_path) -> None:
    import torch

    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    # H stays NumPy while the overlap eigenvector is a tensor: the two specs
    # must be independent in both directions.
    opt_a._prev_eigvec_min = torch.zeros(3, dtype=torch.float32)

    info = opt_a._get_opt_restart_info()
    assert info["H_spec"]["backend"] == "numpy"
    assert info["_prev_eigvec_min_spec"]["backend"] == "torch"
    assert info["_prev_eigvec_min_spec"]["dtype"] == "float32"
    # The full device string, including any CUDA index, is preserved verbatim.
    assert info["_prev_eigvec_min_spec"]["device"] == str(
        opt_a._prev_eigvec_min.device
    )

    _, opt_b = _rfo(tmp_path / "b", [0.5, 0.0, 0.0], max_cycles=2)
    opt_b._set_opt_restart_info(info)

    assert isinstance(opt_b.H, np.ndarray)
    assert isinstance(opt_b._prev_eigvec_min, torch.Tensor)
    assert opt_b._prev_eigvec_min.dtype == torch.float32


@pytest.mark.gpu
def test_restart_spec_records_the_full_cuda_device_string(tmp_path) -> None:
    import torch

    if not torch.cuda.is_available():
        pytest.skip("CUDA is unavailable")

    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    device = torch.device("cuda", torch.cuda.current_device())
    opt_a.H = torch.as_tensor(np.array(opt_a.H), dtype=torch.float64, device=device)

    info = opt_a._get_opt_restart_info()
    # The recorded device keeps its index, not just the device type.
    assert info["H_spec"]["device"] == str(device)
    assert ":" in info["H_spec"]["device"]
    # A CUDA tensor still serializes to plain lists.
    assert _only_primitives(info)

    _, opt_b = _rfo(tmp_path / "b", [0.5, 0.0, 0.0], max_cycles=2)
    opt_b._set_opt_restart_info(info)
    assert str(opt_b.H.device) == str(device)


def test_unavailable_checkpoint_device_falls_back_to_the_cpu(tmp_path) -> None:
    import torch

    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    info = opt_a._get_opt_restart_info()
    info["H_spec"] = {
        "backend": "torch",
        "dtype": "float64",
        "device": "cuda:31",
    }

    _, opt_b = _rfo(tmp_path / "b", [0.5, 0.0, 0.0], max_cycles=2)
    opt_b._set_opt_restart_info(info)

    assert isinstance(opt_b.H, torch.Tensor)
    assert opt_b.H.device.type == "cpu"


def test_buffer_dtype_is_preserved_without_promoting_to_float64(tmp_path) -> None:
    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    opt_a._sy_buffer_S = [np.ones(3, dtype=np.float32)]
    opt_a._sy_buffer_Y = [np.ones(3, dtype=np.float32)]

    info = opt_a._get_opt_restart_info()
    assert info["_sy_buffer_S_spec"]["dtype"] == "float32"

    _, opt_b = _rfo(tmp_path / "b", [0.5, 0.0, 0.0], max_cycles=2)
    opt_b._set_opt_restart_info(info)
    assert opt_b._sy_buffer_S[0].dtype == np.float32
    assert opt_b._sy_buffer_Y[0].dtype == np.float32


def test_numpy_hessian_state_still_restores_without_a_backend_spec(tmp_path) -> None:
    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    info = opt_a._get_opt_restart_info()
    info.pop("H_spec")

    _, opt_b = _rfo(tmp_path / "b", [0.5, 0.0, 0.0], max_cycles=2)
    opt_b._set_opt_restart_info(info)

    assert isinstance(opt_b.H, np.ndarray)


def test_restart_hessian_is_validated_against_the_resolved_active_subset(
    tmp_path,
) -> None:
    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)
    payload = checkpoint.load_payload(ck)

    _, opt_b = _rfo(tmp_path / "b", [0.5, 0.0, 0.0], max_cycles=2)
    # One active atom out of one, but only its first two DOFs are active, so the
    # full coordinate size would wrongly accept a 3x3 Hessian.
    opt_b._using_active_dofs = True
    opt_b._active_dof_indices = np.array([0, 1], dtype=int)

    with pytest.raises(checkpoint.CheckpointValidationError):
        checkpoint.validate_payload(payload, opt_b)

    active = copy.deepcopy(payload)
    active["restart_info"]["H"] = [[1.0, 0.0], [0.0, 1.0]]
    active["restart_info"]["H_spec"]["shape"] = [2, 2]
    assert checkpoint.validate_payload(active, opt_b)


def test_zero_active_dimension_round_trips_an_empty_restart_hessian(tmp_path) -> None:
    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    # A zero-dimensional active space: the (0, 0) Hessian is the correct payload.
    opt_a.H = np.empty((0, 0))
    opt_a._using_active_dofs = True
    opt_a._active_dof_indices = np.array([], dtype=int)

    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)
    payload = checkpoint.load_payload(ck)
    # Nested lists cannot carry the shape, so it travels in the spec.
    assert payload["restart_info"]["H"] == []
    assert payload["restart_info"]["H_spec"]["shape"] == [0, 0]

    _, opt_b = _rfo(tmp_path / "b", [0.5, 0.0, 0.0], max_cycles=2)
    opt_b._using_active_dofs = True
    opt_b._active_dof_indices = np.array([], dtype=int)

    # Full save -> load -> apply round trip, not validation alone.
    checkpoint.load_and_apply(opt_b, ck)
    assert opt_b.H.shape == (0, 0)

    # A three-dimensional Hessian is rejected against the zero active space.
    nonempty = copy.deepcopy(payload)
    nonempty["restart_info"]["H"] = np.eye(3).tolist()
    nonempty["restart_info"]["H_spec"]["shape"] = [3, 3]
    with pytest.raises(checkpoint.CheckpointValidationError):
        checkpoint.validate_payload(nonempty, opt_b)

    # A recorded shape that contradicts the payload size is rejected.
    inconsistent = copy.deepcopy(payload)
    inconsistent["restart_info"]["H"] = np.eye(3).tolist()
    inconsistent["restart_info"]["H_spec"]["shape"] = [0, 0]
    with pytest.raises(checkpoint.CheckpointValidationError):
        checkpoint.validate_payload(inconsistent, opt_b)


def test_legacy_empty_hessian_without_a_recorded_shape_restores_as_zero_square(
    tmp_path,
) -> None:
    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    info = opt_a._get_opt_restart_info()
    # Emulate a checkpoint written before the shape was recorded.
    info["H"] = []
    info.pop("H_spec")

    _, opt_b = _rfo(tmp_path / "b", [0.5, 0.0, 0.0], max_cycles=2)
    opt_b._set_opt_restart_info(info)
    assert opt_b.H.shape == (0, 0)


def test_restart_shapes_survive_a_full_round_trip_for_every_array(tmp_path) -> None:
    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    opt_a._sy_buffer_S = [np.zeros(3)]
    opt_a._sy_buffer_Y = [np.zeros(3)]
    opt_a._prev_eigvec_min = np.zeros(3)

    info = opt_a._get_opt_restart_info()
    assert info["H_spec"]["shape"] == [3, 3]
    assert info["_sy_buffer_S_spec"]["shape"] == [3]
    assert info["_prev_eigvec_min_spec"]["shape"] == [3]

    _, opt_b = _rfo(tmp_path / "b", [0.5, 0.0, 0.0], max_cycles=2)
    opt_b._set_opt_restart_info(info)
    assert opt_b.H.shape == (3, 3)
    assert opt_b._sy_buffer_S[0].shape == (3,)
    assert opt_b._prev_eigvec_min.shape == (3,)


def test_nonnumeric_and_two_dimensional_cart_coords_are_typed_errors(
    tmp_path,
) -> None:
    _, opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)
    payload = checkpoint.load_payload(ck)

    _, opt_b = _rfo(tmp_path / "b", [0.5, 0.0, 0.0], max_cycles=2)

    nonnumeric = copy.deepcopy(payload)
    nonnumeric["restart_info"]["geom_info"]["cart_coords"] = ["a", "b", "c"]
    with pytest.raises(checkpoint.CheckpointValidationError):
        checkpoint.validate_payload(nonnumeric, opt_b)

    nested = copy.deepcopy(payload)
    nested["restart_info"]["geom_info"]["cart_coords"] = [[0.0, 0.0, 0.0]]
    with pytest.raises(checkpoint.CheckpointValidationError):
        checkpoint.validate_payload(nested, opt_b)
