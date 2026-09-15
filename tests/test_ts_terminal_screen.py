"""Model-screen boundaries preserve exact PHVA evidence and acceptance."""

from copy import deepcopy

import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from test_ts_terminal_cadence import CountedQuadratic, run_quadratic


ELIGIBLE_STEP = np.array([0., 1e-4, 0.])
SCOPE = {"active_atoms": [3], "frozen_atoms": [0, 1, 2],
         "treatment": "constrained", "frequency_zero_cutoff_cm": 5.}


@pytest.fixture
def prior_hosp(tmp_path):
    geom, opt, calc, _, _ = run_quadratic(tmp_path, 500)
    assert opt._last_exact_validation == "higher_order"
    assert opt._last_exact_n_imaginary == 2
    assert opt._last_exact_modes is not None
    assert not opt._exact_phva_matches_current_geometry()
    assert opt._all_configured_values_met(ELIGIBLE_STEP)
    assert {key: opt._last_rigid_projection_info[key] for key in SCOPE} == SCOPE
    np.testing.assert_allclose(opt.cur_H, np.diag(calc.diagonal), atol=1e-12)
    return geom, opt, calc


def _copied_value(value):
    if hasattr(value, "detach"):
        return value.detach().cpu().numpy().copy()
    return deepcopy(value)


def screen_without_state_change(opt, calc, step=ELIGIBLE_STEP):
    """Snapshot after deliberate test mutations; screening must be read-only."""
    missing = object()
    names = [key for key in vars(opt) if key.startswith("_last_exact_")]
    names += ["_last_rigid_projection_info", "cur_H", "H",
              "exact_saddle_checks", "hessian_recalc_in"]
    before = {key: (getattr(opt, key, missing), _copied_value(getattr(opt, key, None)))
              for key in names}
    calls, diagonal = len(calc.hessian_calls), calc.diagonal.copy()
    result = opt._defer_hosp_terminal_check(step)
    for key, (original, value) in before.items():
        current = getattr(opt, key, missing)
        assert current is original, f"screen replaced {key}"
        if original is missing:  # fresh optimizer has no PHVA metadata yet
            continue
        if isinstance(value, np.ndarray):
            np.testing.assert_array_equal(_copied_value(current), value, err_msg=key)
        else:
            assert current == value, f"screen changed {key}"
    assert len(calc.hessian_calls) == calls, "screen acquired an exact Hessian"
    np.testing.assert_array_equal(calc.diagonal, diagonal)
    return result


@pytest.mark.parametrize("model, expected", [
    ([-.04, -.02, .10], True), ([-.04, .02, .10], False),
    ([.04, .02, .10], False), ([-.04, -1e-12, .10], True),
    ([-.04, np.nan, .10], False), ([.04, .02], False), (None, False),
], ids=["surplus", "fosp", "minimum", "soft-negative", "nonfinite-H",
        "unsupported-block", "missing-H"])
def test_real_model_phva_and_exact_state(prior_hosp, monkeypatch, model, expected):
    _, opt, calc = prior_hosp
    # Explicit stale approximation only: true E/g/H retain two negative roots.
    opt.cur_H = None if model is None else np.diag(model)
    original, observed = opt._mw_frequencies_and_modes, []

    def observe():
        observed.append(True)
        return original()

    monkeypatch.setattr(opt, "_mw_frequencies_and_modes", observe)
    assert screen_without_state_change(opt, calc) is expected
    if model is not None and len(model) == 3 and np.isfinite(model).all():
        assert len(observed) == 1, "ordinary model cases must exercise real PHVA"
    assert not opt._exact_terminal_candidate_matches_current_geometry()


@pytest.mark.parametrize("attribute, value", [
    ("verify_saddle", False), ("_saddle_recovery_active", True),
    ("stop_requested", True), ("flatten_enabled", True), ("hessian_xtb", True),
    ("hessian_recalc", None), ("hessian_recalc", np.inf),
    ("hessian_recalc", 0), ("hessian_recalc", -1),
    ("_last_exact_validation", "unavailable"),
    ("_last_exact_n_negative", None), ("_last_exact_n_negative", 1),
    ("_last_exact_cart_coords", None), ("_last_exact_cart_coords", np.zeros(3)),
    ("_last_exact_frequencies_cm", None),
    ("_last_exact_frequencies_cm", np.array([-100., np.nan, 50.])),
    ("_last_rigid_projection_info", None),
])
def test_ineligible_state_never_defers(prior_hosp, attribute, value):
    _, opt, calc = prior_hosp
    setattr(opt, attribute, deepcopy(value))
    assert screen_without_state_change(opt, calc) is False


@pytest.mark.parametrize("reason", ["same-coordinates", "large-step", "large-force",
                                        "energy-jump", "nonfinite-step", "never"])
def test_numerical_and_geometry_eligibility(prior_hosp, reason):
    geom, opt, calc = prior_hosp
    step = ELIGIBLE_STEP.copy()
    if reason == "same-coordinates":
        opt._last_exact_cart_coords = geom.cart_coords.copy()
    elif reason == "large-step":
        step *= 100
    elif reason == "large-force":
        opt.forces, opt.modified_forces = [np.ones(12)], []
    elif reason == "energy-jump":
        opt.energies = [0., 1.]  # nonzero step avoids Baker's zero-step exception
    elif reason == "nonfinite-step":
        step[0] = np.nan
    else:
        opt.thresh = "never"
    assert screen_without_state_change(opt, calc, step) is False


@pytest.mark.parametrize("key, mismatch", [
    ("active_atoms", [2]), ("frozen_atoms", [0, 1]),
    ("treatment", "global"), ("frequency_zero_cutoff_cm", 6.),
])
@pytest.mark.parametrize("missing", [False, True], ids=["mismatch", "missing"])
def test_each_exact_scope_key_is_required(prior_hosp, key, mismatch, missing):
    _, opt, calc = prior_hosp
    # Corrupt stored scope, not physical settings or the current PHVA kernel.
    if missing:
        del opt._last_rigid_projection_info[key]
    else:
        opt._last_rigid_projection_info[key] = mismatch
    assert screen_without_state_change(opt, calc) is False


@pytest.mark.parametrize("fault", ["none", "exception", "nan", "inf", "model-scope"])
def test_failed_model_probe_restores_exact_metadata(prior_hosp, monkeypatch, fault):
    _, opt, calc = prior_hosp
    original = opt._mw_frequencies_and_modes

    def fail_after_real_phva():
        frequencies, modes = original()  # replaces metadata as native PHVA does
        if fault == "none":
            return None
        if fault == "exception":
            raise RuntimeError("deliberate model-screen fault after metadata replacement")
        if fault == "model-scope":
            del opt._last_rigid_projection_info["active_atoms"]
        else:
            frequencies = frequencies.copy()
            frequencies[-1] = np.nan if fault == "nan" else np.inf
        return frequencies, modes

    monkeypatch.setattr(opt, "_mw_frequencies_and_modes", fail_after_real_phva)
    assert screen_without_state_change(opt, calc) is False


def test_first_stale_fosp_model_still_gets_exact_hosp_gate(tmp_path):
    """Fresh real optimizer: no fabricated prior exact result or kernel oracle."""
    geom = Geometry(["H"] * 4, [0., 0., 0., 2., 0., 0., 0., 2., 0., 0., 0., 1.],
                    coord_type="cart", freeze_atoms=[0, 1, 2])
    calc = CountedQuadratic(tmp_path, [-.04, -.02, .10])
    geom.set_calculator(calc)
    opt = RSPRFOptimizer(
        geom, hessian_init="unit", hessian_update="bofill", hessian_recalc=500,
        hessian_recalc_adapt=None, hessian_xtb=False, thresh="baker", max_cycles=1,
        verify_saddle=True, saddle_imaginary_threshold_cm=5.,
        saddle_recovery_max_cycles=0, flatten_enabled=False,
        reference_mode=np.r_[np.zeros(9), 1., 0., 0.], out_dir=tmp_path, dump=False,
    )
    calc.optimizer = opt
    opt.prepare_opt()  # unit guess; no calculator Hessian request
    opt.cur_H = np.diag([-.04, .02, .10])  # explicitly stale approximate FOSP
    physical = calc.get_forces(geom.atoms, geom.cart_coords)
    opt.forces = [physical["forces"]]
    opt.energies = [physical["energy"], physical["energy"]]
    assert not calc.hessian_calls and opt.exact_saddle_checks == 0
    assert opt._last_exact_cart_coords is None
    assert opt._all_configured_values_met(ELIGIBLE_STEP)
    assert screen_without_state_change(opt, calc) is False
    opt.validate_terminal_saddle_for_step(ELIGIBLE_STEP)  # actual exact gate
    assert len(calc.hessian_calls) == opt.exact_saddle_checks == 1
    assert opt._last_exact_n_imaginary == 2
    assert opt._last_exact_validation == "higher_order"
    assert opt._exact_phva_matches_current_geometry()
    assert not opt._exact_terminal_candidate_matches_current_geometry()
    np.testing.assert_allclose(opt.cur_H, np.diag(calc.diagonal), atol=1e-12)


@pytest.mark.parametrize("exact_values,exact_near,model_near,expected", [
    ([-100., -20.], [2.], [-2.], True),
    ([-100., 50.], [-2.], [2.], False),
])
def test_model_near_sign_uses_model_not_restored_exact_metadata(
    prior_hosp, monkeypatch, exact_values, exact_near, model_near, expected,
):
    _, opt, calc = prior_hosp
    from pysisyphus.normal_modes import frequency_partition_info

    exact_complete = np.sort(np.r_[exact_values, exact_near])
    opt._last_exact_frequencies_cm = exact_complete
    opt._last_exact_n_imaginary = int(np.count_nonzero(exact_complete < -5.))
    opt._last_exact_n_negative = 2
    opt._last_rigid_projection_info.update(frequency_partition_info(exact_complete, 5.))
    exact_info = opt._last_rigid_projection_info

    def model_packet():
        model_complete = np.sort(np.r_[-100., 50., model_near])
        opt._last_rigid_projection_info = {
            **exact_info, **frequency_partition_info(model_complete, 5.),
        }
        return model_complete, np.zeros((3, 12))

    monkeypatch.setattr(opt, "_mw_frequencies_and_modes", model_packet)
    assert screen_without_state_change(opt, calc) is expected
    assert opt._last_rigid_projection_info is exact_info


@pytest.mark.parametrize("change", ["missing-near", "incomplete", "nonfinite-near"])
def test_model_partition_failure_keeps_exact_check(prior_hosp, monkeypatch, change):
    _, opt, calc = prior_hosp
    original = opt._mw_frequencies_and_modes

    def incomplete_model():
        data = original()
        info = opt._last_rigid_projection_info
        if change == "missing-near":
            info.pop("near_zero_frequencies_cm")
        elif change == "incomplete":
            info["raw_mode_count"] += 1
        else:
            info["near_zero_frequencies_cm"] = [np.nan]
            info["raw_mode_count"] += 1
            info["near_zero_mode_count"] = 1
        return data

    monkeypatch.setattr(opt, "_mw_frequencies_and_modes", incomplete_model)
    assert screen_without_state_change(opt, calc) is False
