"""C6 — truthful scientific outcomes (mlmm_toolkit).

Every test in this file asserts the load-bearing C6 invariant: no artifact
existence and no finite fallback may promote a required nonconverged / missing
scientific leaf to success, while a genuinely converged leaf's public output is
unchanged aside from the additive outcome fields.

The falsifiers are grouped by the false-promotion path each one closes; every one
of them would have reported SUCCESS (or a wrong minimum / a cached Hessian / a
Gibbs diagram / a fabricated 0.0-substituted thermochemistry) under the pre-C6
fallback and now correctly reports nonconverged / missing / excluded.

Each production-path falsifier binds to the production helper it exercises
(imported from ``mlmm.workflows._outcomes`` / the workflow modules); the gate
logic is never re-implemented inside the test.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from mlmm.workflows._outcomes import (
    AggregateTruth,
    LeafOutcome,
    ScanPointOutcome,
    aggregate_workflow_truth,
    attach_outcomes,
    eligible_points,
    make_leaf,
    make_scan_point,
    scan_scientific_status,
    seed_eligible_mask,
)


# ---------------------------------------------------------------------------
# 1. Pure outcome types + serializer compatibility (P07)
# ---------------------------------------------------------------------------


def test_leaf_outcome_fail_closed_usability() -> None:
    # Executed + explicitly converged + finite -> usable.
    assert make_leaf("s", "a", executed=True, converged=True).usable is True
    # Explicitly not converged -> unusable, even though it executed.
    nc = make_leaf("s", "a", executed=True, converged=False)
    assert nc.usable is False and nc.reason == "not_converged"
    # Unknown convergence (None) fails closed -> unusable.
    unk = make_leaf("s", "a", executed=True, converged=None)
    assert unk.usable is False and unk.reason == "convergence_unknown"
    # Converged but energy invalid -> unusable.
    bad = make_leaf("s", "a", executed=True, converged=True, energy_valid=False)
    assert bad.usable is False and bad.reason == "energy_invalid"


def test_scan_point_seed_eligibility_fail_closed() -> None:
    ok = make_scan_point("p", executed=True, converged=True, energy=-1.0, artifact_written=True)
    assert ok.seed_eligible is True
    # Nonconverged is never a seed even with a finite energy + artifact.
    nc = make_scan_point("p", executed=True, converged=False, energy=-9.0, artifact_written=True)
    assert nc.seed_eligible is False
    # Converged but the artifact write failed -> not eligible (artifact_missing).
    na = make_scan_point("p", executed=True, converged=True, energy=-1.0, artifact_written=False)
    assert na.seed_eligible is False and na.reason == "artifact_missing"
    # NaN energy -> not eligible.
    nan = make_scan_point("p", executed=True, converged=True, energy=float("nan"), artifact_written=True)
    assert nan.seed_eligible is False and nan.reason == "energy_invalid"


def test_scan2d_failed_payload_exists_without_usable_plot_point() -> None:
    from mlmm.workflows.scan2d import _build_scan2d_result_payload

    payload = _build_scan2d_result_payload(
        records=[
            {
                "i": 0,
                "j": 0,
                "bias_converged": False,
                "energy_hartree": float("nan"),
                "artifact_written": False,
            }
        ],
        calc_cfg={"backend": "uma", "model_charge": 0, "model_mult": 1},
        pair1={"i": 1, "j": 2, "low": 1.0, "high": 2.0},
        pair2={"i": 3, "j": 4, "low": 1.0, "high": 2.0},
        files={"surface_csv": "surface.csv"},
        status="failed",
    )

    assert payload["status"] == "failed"
    assert payload["execution_status"] == "completed"
    assert payload["scientific_status"] == "failed"
    assert payload["n_points_attempted"] == 1
    assert payload["n_points_usable"] == 0
    assert payload["min_energy_hartree"] is None
    assert payload["files"] == {"surface_csv": "surface.csv"}


def test_aggregate_success_partial_failed() -> None:
    # All required usable, no missing expected -> success.
    t = aggregate_workflow_truth(
        [make_leaf("p", "seg_1", executed=True, converged=True)], ["seg_1"]
    )
    assert (t.scientific_status, t.execution_status) == ("success", "completed")
    # A usable diagnostic artifact + a missing required leaf -> partial.
    raw = LeafOutcome("p", "raw", required=True, executed=True, converged=True,
                      usable=False, reason="endpoint_hei", artifacts=("mep.pdb",))
    t2 = aggregate_workflow_truth([raw], ["seg_1"])
    assert t2.scientific_status == "partial"
    assert "missing:seg_1" in t2.status_reasons
    # Nothing usable, no artifact, execution failed -> failed.
    t3 = aggregate_workflow_truth(
        [make_leaf("p", "seg_1", required=True, executed=False, converged=None)], ["seg_1"]
    )
    assert (t3.scientific_status, t3.execution_status) == ("failed", "failed")


def test_serializer_roundtrip_and_additive_only(tmp_path: Path) -> None:
    from mlmm.core.utils import write_result_json, RESULT_JSON_SCHEMA_VERSION

    leaf = make_leaf("scan", "stage_1", executed=True, converged=True)
    truth = aggregate_workflow_truth([leaf], ["stage_1"])
    # A leaf-command result: legacy status stays "completed", outcomes are added.
    data = {"status": "completed", "min_energy_hartree": -1.5}
    attach_outcomes(data, truth=truth, stage_outcomes=[leaf])
    write_result_json(tmp_path, dict(data), command="scan")
    r = json.loads((tmp_path / "result.json").read_text())

    # Legacy contract preserved.
    assert r["status"] == "completed"
    assert r["schema_version"] == RESULT_JSON_SCHEMA_VERSION
    assert r["min_energy_hartree"] == -1.5
    # Additive outcome fields present and truthful.
    assert r["scientific_status"] == "success"
    assert r["execution_status"] == "completed"
    assert r["expected_item_ids"] == ["stage_1"]
    assert r["stage_outcomes"][0]["usable"] is True

    # An old reader that only knows `status` obtains the same type/value and is
    # unaffected by the additive fields.
    assert isinstance(r["status"], str) and r["status"] == "completed"


def test_dataclasses_are_json_safe() -> None:
    assert LeafOutcome("s", "a").to_dict()["item_id"] == "a"
    assert ScanPointOutcome("p").to_dict()["seed_eligible"] is False
    assert AggregateTruth("completed", "success").to_dict()["scientific_status"] == "success"


# ---------------------------------------------------------------------------
# 2. FALSIFIER — scan non-convergence (M50 + M48)
#    A failed point with the numerically lowest energy would have become the
#    reported minimum / baseline under the old raw-min; it must be excluded.
# ---------------------------------------------------------------------------


def test_scan_failed_low_energy_point_excluded_from_minimum() -> None:
    # These carry the seed-eligibility fields of the records scan2d/scan3d build
# (tri-state bias_converged).
    records = [
        {"i": 0, "j": 0, "energy_hartree": -1.00, "bias_converged": True, "artifact_written": True},
        {"i": 0, "j": 1, "energy_hartree": -1.20, "bias_converged": True, "artifact_written": True},
        # A FAILED point that happens to have the lowest raw energy.
        {"i": 0, "j": 2, "energy_hartree": -9.99, "bias_converged": False, "artifact_written": True},
    ]
    mask = seed_eligible_mask(records)
    assert mask == [True, True, False]

    raw_min = min(r["energy_hartree"] for r in records)
    eligible = [r for r, ok in zip(records, mask) if ok]
    eligible_min = min(r["energy_hartree"] for r in eligible)

    # OLD behavior would have reported the failed point's energy.
    assert raw_min == -9.99
    # NEW behavior: the reported minimum comes only from converged points.
    assert eligible_min == -1.20

    status, reasons = scan_scientific_status(
        [
            make_scan_point(f"{r['i']}_{r['j']}", executed=True,
                            converged=r["bias_converged"], energy=r["energy_hartree"],
                            artifact_written=r["artifact_written"])
            for r in records
        ]
    )
    assert status == "partial"
    assert reasons == ("unusable_points:1",)


def test_scan_normal_return_is_not_convergence() -> None:
    # M50: an optimizer that returns normally but reports is_converged=False (a
    # cycle-limit stop) must not be recorded as converged. seed_eligible_mask
    # reads the recorded bit; only an explicit True survives.
    normal_return_but_not_converged = {
        "energy_hartree": -1.0, "bias_converged": False,
        "artifact_written": True,
    }
    unknown_convergence = {
        "energy_hartree": -1.0, "bias_converged": None,
        "artifact_written": True,
    }
    converged = {
        "energy_hartree": -1.0, "bias_converged": True,
        "artifact_written": True,
    }
    assert seed_eligible_mask(
        [normal_return_but_not_converged, unknown_convergence, converged]
    ) == [False, False, True]


def test_scan_stage_leaf_partial_when_middle_step_fails() -> None:
    # Binds to the extracted _outcomes.combine_step_convergence helper that
    # scan.py's production stage-leaf fold calls (this test exercises the shared
    # helper, not scan.py's own code path): stage_1 all steps converged; stage_2
    # has a middle step that failed (OptimizationError) while the final step
    # converged. Legacy per-stage `converged`/`status` would still read the final
    # step, but the aggregate scientific_status must be partial because a
    # converged final step cannot hide the failure.
    from mlmm.workflows._outcomes import combine_step_convergence as _combine

    assert _combine([True, True, True]) is True
    assert _combine([True, False, True]) is False  # middle-step failure survives
    assert _combine([True, None, True]) is None     # unknown fails closed
    assert _combine([]) is None                      # no steps -> unknown

    stage1 = make_leaf("scan", "stage_1", executed=True,
                       converged=_combine([True, True, True]))
    # middle step failed, final step converged -> combined is False.
    stage2 = make_leaf("scan", "stage_2", executed=True,
                       converged=_combine([True, False, True]))
    assert stage1.usable is True and stage2.usable is False

    truth = aggregate_workflow_truth([stage1, stage2], ["stage_1", "stage_2"])
    assert truth.scientific_status == "partial"
    assert any("stage_2" in r for r in truth.status_reasons)


def test_m50_producer_records_nonconvergence_from_optimizer() -> None:
    # scan2d/scan3d record `bias_converged = optimizer_converged_bit(opt)`; scan
    # records the same bit per step.
    from mlmm.workflows._outcomes import optimizer_converged_bit

    class _FakeOpt:
        def __init__(self, conv):
            self.is_converged = conv

    assert optimizer_converged_bit(_FakeOpt(False)) is False
    assert optimizer_converged_bit(_FakeOpt(True)) is True
    # A non-boolean 1 (truthy but not convergence) collapses to unknown.
    assert optimizer_converged_bit(_FakeOpt(1)) is None

    class _NoAttr:
        pass

    assert optimizer_converged_bit(_NoAttr()) is None

    # The recorded bit drives seed eligibility: a normal-return-but-nonconverged
    # point is excluded even with a finite (lowest) energy.
    record = {"energy_hartree": -9.9, "bias_converged": optimizer_converged_bit(_FakeOpt(False))}
    assert seed_eligible_mask([record]) == [False]


def test_eligible_points_helper() -> None:
    pts = [
        make_scan_point("a", executed=True, converged=True, energy=-1.0, artifact_written=True),
        make_scan_point("b", executed=True, converged=False, energy=-2.0, artifact_written=True),
    ]
    assert [p.point_id for p in eligible_points(pts)] == ["a"]


# ---------------------------------------------------------------------------
# 3. FALSIFIER — path expected-segment miss / endpoint HEI (M29)
#    An endpoint-HEI raw diagram (segments=[]) would have been promoted to
#    success because a diagram exists; it must be partial.
# ---------------------------------------------------------------------------


def test_path_endpoint_hei_zero_segments_is_partial_not_success() -> None:
    from mlmm.workflows.path_search import _path_leaves_and_expected

    # Endpoint-HEI branch returns segments=[] but a raw R/P diagram can be drawn.
    leaves, expected = _path_leaves_and_expected(
        [], raw_artifacts=["mep.pdb", "energy_diagram_MEP.png"]
    )
    truth = aggregate_workflow_truth(leaves, expected)
    assert truth.scientific_status == "partial"  # would have been "success"
    raw = [leaf for leaf in leaves if leaf.item_id == "raw_path"][0]
    assert raw.usable is False and raw.reason == "endpoint_hei"
    assert "mep.pdb" in raw.artifacts


def test_path_engine_nonconverged_endpoint_reason_retained() -> None:
    from mlmm.workflows.path_search import _path_leaves_and_expected

    leaves, expected = _path_leaves_and_expected(
        [], raw_artifacts=["mep.pdb"], engine_converged=False
    )
    raw = [leaf for leaf in leaves if leaf.item_id == "raw_path"][0]
    assert "endpoint_hei" in raw.reason and "engine_nonconverged" in raw.reason


def test_path_summary_contract_is_versioned_and_endpoint_fail_closed(tmp_path: Path) -> None:
    from mlmm.core.utils import RESULT_JSON_SCHEMA_VERSION
    from mlmm.workflows.path_search import _enrich_path_summary_contract

    for name in ("mep.pdb", "energy_diagram_MEP.png"):
        (tmp_path / name).write_text("current\n", encoding="utf-8")
    summary = {"energy_diagrams": [{"name": "energy_diagram_MEP"}]}

    _enrich_path_summary_contract(
        summary,
        segments=[],
        out_dir=tmp_path,
        calc_cfg={"backend": "orb", "model_charge": -1, "model_mult": 1},
        command="mlmm path-search -i r.pdb -i p.pdb",
    )

    assert summary["schema_version"] == RESULT_JSON_SCHEMA_VERSION
    assert summary["status"] == "partial"
    assert summary["execution_status"] == "completed"
    assert summary["scientific_status"] == "partial"
    assert summary["stage_outcomes"][0]["item_id"] == "raw_path"
    assert summary["stage_outcomes"][0]["usable"] is False
    assert summary["charge"] == -1 and summary["spin"] == 1


def test_path_summary_log_reuses_enriched_calculator_provenance(tmp_path: Path) -> None:
    from mlmm.io.summary import write_summary_log
    from mlmm.workflows.path_search import _summary_log_provenance

    summary = {
        "mlip_backend": "orb",
        "mlip_model": "orb-v3",
        "mlip_precision": "fp64",
    }
    payload = {
        "pipeline_mode": "path-search",
        **_summary_log_provenance(summary),
    }
    assert payload == {"pipeline_mode": "path-search", **summary}

    destination = tmp_path / "summary.log"
    write_summary_log(destination, payload)
    rendered = destination.read_text(encoding="utf-8")
    assert "orb" in rendered and "orb-v3" in rendered


def test_path_summary_contract_does_not_swallow_truth_failures(
    tmp_path: Path, monkeypatch,
) -> None:
    from mlmm.workflows import path_search

    monkeypatch.setattr(
        path_search,
        "_path_leaves_and_expected",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("truth failed")),
    )
    with pytest.raises(RuntimeError, match="truth failed"):
        path_search._enrich_path_summary_contract(
            {"energy_diagrams": []},
            segments=[],
            out_dir=tmp_path,
            calc_cfg={"backend": "orb"},
            command="mlmm path-search",
        )


def test_path_single_reactive_segment_is_success() -> None:
    from mlmm.workflows.path_search import SegmentReport, _path_leaves_and_expected

    seg = SegmentReport(tag="seg_001_refine", barrier_kcal=10.0, delta_kcal=-2.0,
                        summary="", kind="seg", seg_index=1, converged=True)
    leaves, expected = _path_leaves_and_expected([seg])
    truth = aggregate_workflow_truth(leaves, expected)
    assert truth.scientific_status == "success"  # legacy behavior preserved


def test_path_nonconverged_reactive_segment_is_not_success() -> None:
    # M29: a reactive segment whose StringOptimizer hit its cycle limit
    # (is_converged=False) writes its trajectory but must NOT count toward
    # completeness. The pre-C6 producer had no converged field and would have
    # reported success; the real threading now reads the field.
    from mlmm.workflows.path_search import SegmentReport, _path_leaves_and_expected

    seg = SegmentReport(tag="seg_001_refine", barrier_kcal=10.0, delta_kcal=-2.0,
                        summary="", kind="seg", seg_index=1, converged=False)
    leaves, expected = _path_leaves_and_expected([seg])
    truth = aggregate_workflow_truth(leaves, expected)
    assert truth.scientific_status != "success"
    seg_leaf = [leaf for leaf in leaves if leaf.item_id == "segment_1"][0]
    assert seg_leaf.usable is False and seg_leaf.reason == "not_converged"

    # Convergence-unknown (None: no readable signal) also fails closed.
    seg_unk = SegmentReport(tag="seg_001_refine", barrier_kcal=10.0, delta_kcal=-2.0,
                            summary="", kind="seg", seg_index=1, converged=None)
    leaves_u, expected_u = _path_leaves_and_expected([seg_unk])
    assert aggregate_workflow_truth(leaves_u, expected_u).scientific_status != "success"


def test_path_bridge_only_is_not_success() -> None:
    from mlmm.workflows.path_search import SegmentReport, _path_leaves_and_expected

    bridge = SegmentReport(tag="bridge_01", barrier_kcal=0.0, delta_kcal=0.0,
                           summary="", kind="bridge", seg_index=1)
    leaves, expected = _path_leaves_and_expected([bridge], raw_artifacts=["mep.pdb"])
    truth = aggregate_workflow_truth(leaves, expected)
    # A path made only of a non-reactive bridge has no usable reactive segment.
    assert truth.scientific_status != "success"


# ---------------------------------------------------------------------------
# 4. FALSIFIER — IRC directional non-convergence (M42)
#    A nonconverged direction (whose trajectory + Hessian still exist) must not
#    be promoted; only the converged direction is usable.
# ---------------------------------------------------------------------------


def _irc_direction_leaves(*, forward, forward_conv, backward, backward_conv,
                          n_fwd=10, n_bwd=10):
    """Exercise _outcomes.irc_direction_leaves, the directional-leaf builder that
    irc.py calls in production (bind to the shared helper, don't copy it)."""
    from mlmm.workflows._outcomes import irc_direction_leaves
    return irc_direction_leaves(
        (
            ("forward", bool(forward), forward_conv, n_fwd,
             ["forward_irc.pdb"] if forward else []),
            ("backward", bool(backward), backward_conv, n_bwd,
             ["backward_irc.pdb"] if backward else []),
        )
    )


def test_irc_forward_converged_backward_not_is_partial() -> None:
    # Both directions requested and both trajectories/Hessians exist, but only
    # forward converged.
    leaves, expected = _irc_direction_leaves(
        forward=True, forward_conv=True, backward=True, backward_conv=False
    )
    truth = aggregate_workflow_truth(leaves, expected)
    assert truth.scientific_status == "partial"
    backward = [leaf for leaf in leaves if leaf.item_id == "backward"][0]
    assert backward.usable is False  # not promoted despite artifacts existing


def test_irc_one_sided_request_succeeds() -> None:
    # Backward explicitly disabled: its absence is optional, not a failure.
    leaves, expected = _irc_direction_leaves(
        forward=True, forward_conv=True, backward=False, backward_conv=None
    )
    truth = aggregate_workflow_truth(leaves, expected)
    assert truth.scientific_status == "success"
    assert expected == ["forward"]


def test_irc_convergence_attribute_absent_fails_closed() -> None:
    # A missing convergence attribute (None) must fail closed, not read as True.
    leaves, expected = _irc_direction_leaves(
        forward=True, forward_conv=None, backward=False, backward_conv=None
    )
    truth = aggregate_workflow_truth(leaves, expected)
    assert truth.scientific_status != "success"


def test_irc_hessian_cache_gate_condition() -> None:
    # Binds to the PRODUCTION gate irc.py uses to decide endpoint-Hessian caching
    # (irc_hessian_cache_eligible == `getattr(obj, attr, None) is True`).
    from mlmm.workflows._outcomes import irc_hessian_cache_eligible

    class _FakeEuler:
        forward_is_converged = False
        backward_is_converged = True

    e = _FakeEuler()
    # forward: nonconverged -> the Hessian key must NOT be cached.
    assert irc_hessian_cache_eligible(e, "forward_is_converged") is False
    # backward: converged -> eligible for caching.
    assert irc_hessian_cache_eligible(e, "backward_is_converged") is True
    # An absent attribute also fails closed.
    assert irc_hessian_cache_eligible(e, "never_ran_is_converged") is False
    # A non-boolean truthy value (e.g. the integer 1) is not convergence.
    e.forward_is_converged = 1  # type: ignore[assignment]
    assert irc_hessian_cache_eligible(e, "forward_is_converged") is False


# ---------------------------------------------------------------------------
# 5. FALSIFIER — FREQ/DFT name/number fallback (M28, mlmm-specific)
#    A nonzero freq exit must not be treated as success just because a
#    thermoanalysis.yaml with finite fields exists on disk, and a missing thermal
#    correction must never be replaced by a finite 0.0 / MLIP electronic energy.
# ---------------------------------------------------------------------------


def test_freq_nonzero_exit_is_not_usable_despite_finite_yaml(tmp_path: Path, monkeypatch) -> None:
    from mlmm.workflows import all as all_workflow

    fdir = tmp_path / "R"
    fdir.mkdir(parents=True, exist_ok=True)
    # A thermoanalysis.yaml with finite thermochemistry exists (e.g. a prior run
    # or a partial write) — under the old code its finite numbers would be used.
    (fdir / "thermoanalysis.yaml").write_text(
        "sum_EE_and_thermal_free_energy_ha: -123.456\n"
        "thermal_correction_free_energy_ha: 0.05\n"
    )
    pdb = tmp_path / "R.xyz"
    pdb.write_text("1\n\nH 0 0 0\n")

    # Freq child exits nonzero.
    monkeypatch.setattr(all_workflow, "_run_cli_main", lambda *a, **k: 1)
    monkeypatch.setattr(all_workflow, "_echo", lambda *a, **k: None)

    out = all_workflow._run_freq_for_state(
        pdb, 0, 1, tmp_path / "real.parm7", tmp_path / "model.pdb", False,
        fdir, None, overrides={},
    )
    # The finite YAML is NOT promoted; thermochemistry is unusable.
    assert out == {}


def test_freq_zero_exit_returns_parsed_thermo(tmp_path: Path, monkeypatch) -> None:
    from mlmm.workflows import all as all_workflow

    fdir = tmp_path / "R"
    fdir.mkdir(parents=True, exist_ok=True)
    (fdir / "thermoanalysis.yaml").write_text(
        "sum_EE_and_thermal_free_energy_ha: -123.456\n"
    )
    pdb = tmp_path / "R.xyz"
    pdb.write_text("1\n\nH 0 0 0\n")

    # Freq child succeeds (exit 0).
    monkeypatch.setattr(all_workflow, "_run_cli_main", lambda *a, **k: 0)
    monkeypatch.setattr(all_workflow, "_echo", lambda *a, **k: None)

    out = all_workflow._run_freq_for_state(
        pdb, 0, 1, tmp_path / "real.parm7", tmp_path / "model.pdb", False,
        fdir, None, overrides={},
    )
    assert out.get("sum_EE_and_thermal_free_energy_ha") == pytest.approx(-123.456)


def test_freq_no_dump_does_not_consume_stale_thermo(tmp_path: Path, monkeypatch) -> None:
    from mlmm.workflows import all as all_workflow

    fdir = tmp_path / "R"
    fdir.mkdir(parents=True)
    (fdir / "thermoanalysis.yaml").write_text(
        "sum_EE_and_thermal_free_energy_ha: -123.456\n",
        encoding="utf-8",
    )
    structure = tmp_path / "R.xyz"
    structure.write_text("1\n\nH 0 0 0\n", encoding="utf-8")
    monkeypatch.setattr(all_workflow, "_run_cli_main", lambda *a, **k: 0)

    assert all_workflow._run_freq_for_state(
        structure,
        0,
        1,
        tmp_path / "real.parm7",
        tmp_path / "model.pdb",
        False,
        fdir,
        None,
        overrides={"dump": False},
    ) == {}


@pytest.mark.parametrize("overrides, expected_count", [({}, 0), ({"symmetry_number": 3}, 1)])
def test_all_freq_forwards_symmetry_number_only_when_overridden(
    tmp_path: Path,
    monkeypatch,
    overrides: dict,
    expected_count: int,
) -> None:
    from mlmm.workflows import all as all_workflow

    structure = tmp_path / "R.xyz"
    structure.write_text("1\n\nH 0 0 0\n", encoding="utf-8")
    captured: list[str] = []

    def _capture(_name, _command, argv, **_kwargs):
        captured.extend(argv)
        return 1

    monkeypatch.setattr(all_workflow, "_run_cli_main", _capture)
    monkeypatch.setattr(all_workflow, "_echo", lambda *a, **k: None)

    all_workflow._run_freq_for_state(
        structure,
        0,
        1,
        tmp_path / "real.parm7",
        tmp_path / "model.pdb",
        False,
        tmp_path / "freq",
        None,
        overrides=overrides,
    )

    assert captured.count("--symmetry-number") == expected_count
    if expected_count:
        assert captured[captured.index("--symmetry-number") + 1] == "3"


def test_m28_thermo_gibbs_finite_gate() -> None:
    # Binds to the production finite-gates the all.py Gibbs/DFT//MLIP/MM consumers
    # use in place of the old 0.0 / MLIP substitution: a missing/nonfinite field
    # returns None (so the diagram/dict is skipped, never substituted).
    from mlmm.workflows.all import _thermo_gibbs_ha, _thermo_correction_ha

    # A usable freq payload yields the finite value.
    assert _thermo_gibbs_ha({"sum_EE_and_thermal_free_energy_ha": -123.456}) == pytest.approx(-123.456)
    assert _thermo_correction_ha({"thermal_correction_free_energy_ha": 0.05}) == pytest.approx(0.05)
    # A failed freq (empty {}) or a missing field returns None -> NOT substituted.
    assert _thermo_gibbs_ha({}) is None
    assert _thermo_correction_ha({}) is None  # would have become 0.0 pre-C6
    assert _thermo_gibbs_ha({"sum_EE_and_thermal_free_energy_ha": float("nan")}) is None
    # A nonzero freq exit returns {} from _run_freq_for_state, so the whole
    # requested R/TS/P set is incomplete and the Gibbs diagram is skipped.
    payloads = {"R": {}, "TS": {"sum_EE_and_thermal_free_energy_ha": -1.0},
                "P": {"sum_EE_and_thermal_free_energy_ha": -2.0}}
    gibbs = [_thermo_gibbs_ha(payloads[s]) for s in ("R", "TS", "P")]
    assert None in gibbs  # R failed -> the diagram must NOT be built


def test_m28_dft_energy_gates_on_dft_failed() -> None:
    # Binds to the production helper the all.py DFT / DFT//MLIP/MM consumers use in
    # place of the old `dR.get("energy",{}).get("hartree", eR)` / np.nan MLIP
    # fallback: a finite energy.hartree is NOT trusted when the DFT child failed
    # (_dft_failed) — so the DFT and DFT//MLIP/MM diagrams are skipped, never
    # substituted with an MLIP or 0.0 value.
    from mlmm.workflows.all import _dft_energy_ha, _dft_succeeded

    # A DFT payload with a finite energy but _dft_failed=True -> None (excluded).
    failed = {"energy": {"hartree": -543.21}, "_dft_failed": True}
    assert _dft_succeeded(failed) is False
    assert _dft_energy_ha(failed) is None            # would have been -543.21 pre-C6

    # A genuinely converged DFT payload returns its hartree value unchanged.
    ok = {"energy": {"hartree": -543.21}, "_dft_failed": False}
    assert _dft_succeeded(ok) is True
    assert _dft_energy_ha(ok) == pytest.approx(-543.21)

    # A missing/empty result (no _dft_failed key) fails closed to None: the
    # default treats an absent bit as failed rather than trusting a stray energy.
    assert _dft_energy_ha({}) is None
    assert _dft_energy_ha({"energy": {"hartree": -1.0}}) is None


def test_dft_mlmm_gibbs_uses_subtractive_total_not_raw_model_energy() -> None:
    from mlmm.workflows.all import _dft_mlmm_gibbs_triplet

    dft_results = {
        label: {
            "_dft_failed": False,
            "energy": {"hartree": -500.0},
            "mlmm_energy": {"E_total_ml_dft_mm_hartree": total},
        }
        for label, total in zip(("R", "TS", "P"), (-100.0, -99.9, -100.1))
    }
    thermo = {
        label: {
            "thermal_correction_free_energy_ha": correction,
            "num_imag_freq": n_imag,
        }
        for label, correction, n_imag in zip(
            ("R", "TS", "P"), (0.01, 0.02, 0.03), (0, 1, 0)
        )
    }

    assert _dft_mlmm_gibbs_triplet(dft_results, thermo) == pytest.approx(
        (-99.99, -99.88, -100.07)
    )
    missing_total = dict(dft_results)
    missing_total["TS"] = {
        "_dft_failed": False,
        "energy": {"hartree": -500.0},
    }
    assert _dft_mlmm_gibbs_triplet(missing_total, thermo) is None
    nonminimum = {label: dict(payload) for label, payload in thermo.items()}
    nonminimum["P"]["num_imag_freq"] = 1
    assert _dft_mlmm_gibbs_triplet(dft_results, nonminimum) is None


# ---------------------------------------------------------------------------
# 6. LEGACY-COMPAT — a genuinely converged leaf's public output is unchanged
#    except for the additive new fields.
# ---------------------------------------------------------------------------


def test_legacy_converged_output_is_byte_compatible(tmp_path: Path) -> None:
    from mlmm.core.utils import write_result_json

    # A representative subset of the legacy scan2d result payload for a fully
# converged run.
    legacy = {
        "status": "completed",
        "energy_reference": "bare_mlmm_pes",
        "charge": 0,
        "spin": 1,
        "min_energy_hartree": -1.2345,
        "n_grid_points": 9,
        "files": {"surface_csv": "surface.csv"},
    }

    # Post-C6, the additive outcome fields are appended; every legacy key must be
    # bit-identical.
    points = [
        make_scan_point(f"p{i}", executed=True, converged=True, energy=-1.0, artifact_written=True)
        for i in range(9)
    ]
    sci, reasons = scan_scientific_status(points)
    assert sci == "success"

    enriched = dict(legacy)
    attach_outcomes(enriched, point_outcomes=points, scientific_status=sci,
                    scientific_status_reasons=reasons)

    write_result_json(tmp_path, dict(enriched), command="scan2d")
    r = json.loads((tmp_path / "result.json").read_text())

    for key, value in legacy.items():
        assert r[key] == value, f"legacy key {key} changed"
    # Only additive keys are new; scientific_status reports the truthful success.
    assert r["scientific_status"] == "success"
    assert "scientific_status_reasons" not in r  # no reasons on a clean success
    new_keys = set(r) - set(legacy)
    # The only new keys are additive outcome / envelope fields.
    assert "point_outcomes" in new_keys and "scientific_status" in new_keys


# ---------------------------------------------------------------------------
# 7. FALSIFIER — DMF (path-opt) non-convergence via IPOPT status (M09)
#    Binds ipopt_status_to_converged + the DMF LeafOutcome. A nonconverged solve
#    (status 2) retains its trajectory but is unusable; a converged solve
#    (status 0/1) is a usable success.
# ---------------------------------------------------------------------------


def test_dmf_ipopt_status_maps_to_convergence() -> None:
    from mlmm.workflows._outcomes import ipopt_status_to_converged

    assert ipopt_status_to_converged(0) == (True, "ipopt_converged")
    assert ipopt_status_to_converged(1) == (True, "ipopt_converged")
    assert ipopt_status_to_converged(2) == (False, "ipopt_status_2")
    assert ipopt_status_to_converged(-1) == (False, "ipopt_status_-1")
    # A missing/unreadable status fails closed to unknown (never converged).
    assert ipopt_status_to_converged(None) == (None, "convergence_unknown")


def test_dmf_nonconverged_leaf_unusable_artifact_retained() -> None:
    # path_opt.cli builds `_mk_leaf("path-opt", "dmf_mep", converged=bool(dmf_res.converged),
    # artifacts=["final_geometries_trj.xyz"], reason=dmf_res.reason)`.
    from mlmm.workflows._outcomes import ipopt_status_to_converged

    conv, reason = ipopt_status_to_converged(2)
    leaf = make_leaf("path-opt", "dmf_mep", executed=True, converged=conv,
                     artifacts=["final_geometries_trj.xyz"], reason=reason)
    assert leaf.usable is False                    # not promoted by artifact
    assert "final_geometries_trj.xyz" in leaf.artifacts  # artifact retained
    # Unusable required leaf whose trajectory is retained -> partial (a reportable
    # diagnostic), never success.
    assert aggregate_workflow_truth([leaf], ["dmf_mep"]).scientific_status == "partial"


def test_dmf_converged_leaf_is_success() -> None:
    from mlmm.workflows._outcomes import ipopt_status_to_converged

    conv, reason = ipopt_status_to_converged(0)
    leaf = make_leaf("path-opt", "dmf_mep", executed=True, converged=conv,
                     artifacts=["final_geometries_trj.xyz"], reason=reason)
    assert leaf.usable is True
    assert aggregate_workflow_truth([leaf], ["dmf_mep"]).scientific_status == "success"


def test_dmf_result_legacy_contract_preserved_and_additive_leaf_added() -> None:
    # mlmm's DMF legacy contract is ALREADY convergence-aware
    # (status="converged"/"not_converged", converged=bool); the C6 change is
    # purely additive and MUST NOT flip those legacy fields. The additive
    # scientific_status leaf is fed by the SAME real IPOPT convergence bit.
    from mlmm.workflows import path_opt
    from mlmm.workflows.path_opt import DMFMepResult, _build_dmf_result_data

    conv_res = DMFMepResult(
        images=(), energies=(-1.0, -0.5, -1.2), hei_idx=1,
        converged=True, ipopt_status=0, reason="ok",
    )
    data = _build_dmf_result_data(conv_res, {"model_charge": 0, "model_mult": 1})
    # Legacy fields unchanged by C6 (still derived from the real convergence bit).
    assert data["status"] == "converged"
    assert data["converged"] is True

    nc_res = DMFMepResult(
        images=(), energies=(-1.0, -0.5, -1.2), hei_idx=1,
        converged=False, ipopt_status=2, reason="max_iter",
    )
    nc_data = _build_dmf_result_data(nc_res, {"model_charge": 0, "model_mult": 1})
    assert nc_data["status"] == "not_converged"
    assert nc_data["converged"] is False

    # The additive leaf must be fed by the ONE canonical IPOPT criterion
    # (status 0 or 1, matching path_search), not inferred from artifact
    # existence and not from the legacy status==0 bit.
    src = Path(path_opt.__file__).read_text(encoding="utf-8")
    assert "ipopt_status_to_converged(dmf_res.ipopt_status)" in src
    assert '"dmf_mep"' in src

    # Semantic falsifier for the C6 consistency fix: IPOPT status 1
    # (Solved_To_Acceptable_Level) is the case where the two axes DELIBERATELY
    # diverge. The legacy contract keeps status==0-only (converged=False), but
    # the additive scientific leaf routes through the canonical 0-or-1 criterion
    # and is therefore usable — matching path_search's DMF leaf.
    from mlmm.workflows._outcomes import (
        aggregate_workflow_truth as _agg,
        ipopt_status_to_converged as _ipopt_conv,
        make_leaf as _leaf,
    )
    acc_res = DMFMepResult(
        images=(), energies=(-1.0, -0.5, -1.2), hei_idx=1,
        converged=False, ipopt_status=1, reason="acceptable",
    )
    acc_data = _build_dmf_result_data(acc_res, {"model_charge": 0, "model_mult": 1})
    assert acc_data["status"] == "not_converged"      # legacy: status==0 only
    assert acc_data["converged"] is False
    _conv, _reason = _ipopt_conv(acc_res.ipopt_status)  # canonical: 0 or 1
    assert _conv is True
    acc_leaf = _leaf("path-opt", "dmf_mep", executed=True, converged=_conv,
                     artifacts=["final_geometries_trj.xyz"], reason=_reason)
    assert acc_leaf.usable is True
    assert _agg([acc_leaf], ["dmf_mep"]).scientific_status == "success"


# ---------------------------------------------------------------------------
# 8. FALSIFIER — the ALL-pipeline aggregate consumes truthful leaves (M42 + M29)
#    A never_stop / max-cycle IRC (trajectory present, direction nonconverged)
#    must not yield scientific_status=success in the all-pipeline aggregate,
#    while the legacy `status` string is unchanged.
# ---------------------------------------------------------------------------


def test_read_irc_outcome_gates_on_scientific_status(tmp_path: Path) -> None:
    from mlmm.workflows.all import _read_irc_outcome

    irc_dir = tmp_path / "irc"
    irc_dir.mkdir()
    # A converged both-direction IRC child result.
    (irc_dir / "result.json").write_text(json.dumps({
        "status": "completed",
        "scientific_status": "success",
        "forward_converged": True,
        "backward_converged": True,
        "files": {"finished_irc": "finished_irc_trj.xyz"},
    }))
    ok = _read_irc_outcome(irc_dir)
    assert ok["usable"] is True and ok["traj"] == "finished_irc_trj.xyz"

    # A backward direction that hit its cycle limit: trajectory still exists.
    (irc_dir / "result.json").write_text(json.dumps({
        "status": "completed",              # the IRC PROCESS ran (legacy)
        "scientific_status": "partial",     # but a direction did not converge
        "scientific_status_reasons": ["irc:backward:not_converged"],
        "forward_converged": True,
        "backward_converged": False,
        "files": {"finished_irc": "finished_irc_trj.xyz"},
    }))
    bad = _read_irc_outcome(irc_dir)
    assert bad["usable"] is False
    assert "backward" in bad["reason"]

    # A missing result.json fails closed.
    (irc_dir / "result.json").unlink()
    assert _read_irc_outcome(irc_dir)["usable"] is False


def test_all_pipeline_aggregate_excludes_nonconverged_irc() -> None:
    from mlmm.workflows.all import _pipeline_aggregate_truth

    summary = {"segments": [{
        "index": 1, "kind": "seg", "barrier_kcal": 10.0, "converged": True,
    }]}
    config = {"tsopt": True, "thermo": False, "dft": False}

    # Legacy status="success" (a trajectory exists, so _derive_pipeline_status is
    # satisfied). The IRC leaf reports backward nonconverged -> the aggregate must
    # demote scientific_status while leaving the legacy axis unchanged.
    nonconverged = [{
        "index": 1,
        "irc_traj": "finished_irc_trj.xyz",
        "irc": {"usable": False, "reason": "irc:backward:not_converged",
                "traj": "finished_irc_trj.xyz"},
        "endpoint_opt": {"reactant_converged": True, "product_converged": True},
    }]
    truth = _pipeline_aggregate_truth(
        summary, post_segments=nonconverged, config=config, legacy_status="success",
    )
    assert truth.scientific_status != "success"           # would have been success
    assert any("segment_1" in r for r in truth.status_reasons)

    # Every requested IRC direction converged + endpoints converged -> success.
    converged = [{
        "index": 1,
        "irc_traj": "finished_irc_trj.xyz",
        "irc": {"usable": True, "reason": "ok", "traj": "finished_irc_trj.xyz"},
        "endpoint_assignment": {"connectivity_validated": True},
        "endpoint_opt": {"reactant_converged": True, "product_converged": True},
    }]
    truth_ok = _pipeline_aggregate_truth(
        summary, post_segments=converged, config=config, legacy_status="success",
    )
    assert truth_ok.scientific_status == "success"


def test_all_pipeline_tsopt_only_does_not_require_mep_convergence() -> None:
    """A direct-TS segment has no MEP convergence field to gate on."""
    from mlmm.workflows.all import _pipeline_aggregate_truth

    summary = {"segments": [{"index": 1, "kind": "tsopt", "barrier_kcal": 10.0}]}
    post = [{
        "index": 1,
        "irc": {"usable": True, "reason": "ok"},
        "endpoint_assignment": {"connectivity_validated": True},
        "endpoint_opt": {"reactant_converged": True, "product_converged": True},
    }]
    truth = _pipeline_aggregate_truth(
        summary, post_segments=post, config={"tsopt": True},
        legacy_status="success",
    )
    assert truth.scientific_status == "success"
    assert "mep_convergence_unknown" not in truth.status_reasons


def test_read_opt_endpoint_converged_gates_on_status(tmp_path: Path) -> None:
    # The endpoint-opt child's convergence is read from the SAME real result.json
    # the `opt` subcommand writes (status="converged"/"not_converged"); a
    # missing/unreadable result fails closed to unknown (None) — never a silent
    # promotion.
    from mlmm.workflows.all import _read_opt_endpoint_converged

    opt_dir = tmp_path / "R"
    opt_dir.mkdir()
    # A missing result.json fails closed.
    assert _read_opt_endpoint_converged(opt_dir) is None
    # A converged endpoint opt child.
    (opt_dir / "result.json").write_text(json.dumps({"status": "converged"}))
    assert _read_opt_endpoint_converged(opt_dir) is True
    # A nonconverged endpoint opt child (max-cycle / early stop).
    (opt_dir / "result.json").write_text(json.dumps({"status": "not_converged"}))
    assert _read_opt_endpoint_converged(opt_dir) is False
    # An unreadable / status-less result fails closed to unknown.
    (opt_dir / "result.json").write_text("{ not json")
    assert _read_opt_endpoint_converged(opt_dir) is None


def test_all_pipeline_aggregate_gates_on_endpoint_opt(tmp_path: Path) -> None:
    # Production-bound: the endpoint_opt record fed to the aggregate is assembled
    # from the REAL reader (`_read_opt_endpoint_converged`) reading the opt
    # child's result.json exactly as the `all` producer does — not from a
    # hand-set literal. A nonconverged product endpoint therefore demotes the
    # aggregate below success.
    from mlmm.workflows.all import (
        _pipeline_aggregate_truth,
        _read_opt_endpoint_converged,
    )

    react_dir = tmp_path / "R"
    prod_dir = tmp_path / "P"
    react_dir.mkdir()
    prod_dir.mkdir()
    (react_dir / "result.json").write_text(json.dumps({"status": "converged"}))
    (prod_dir / "result.json").write_text(json.dumps({"status": "not_converged"}))

    summary = {"segments": [{
        "index": 1, "kind": "seg", "barrier_kcal": 10.0, "converged": True,
    }]}
    config = {"tsopt": True}
    # The producer assembles segment_log["endpoint_opt"] from the reader's bits.
    post = [{
        "index": 1,
        "irc": {"usable": True, "reason": "ok"},
        "endpoint_opt": {
            "reactant_converged": _read_opt_endpoint_converged(react_dir),
            "product_converged": _read_opt_endpoint_converged(prod_dir),
        },
    }]
    assert post[0]["endpoint_opt"]["product_converged"] is False   # from the real read
    truth = _pipeline_aggregate_truth(
        summary, post_segments=post, config=config, legacy_status="success",
    )
    assert truth.scientific_status != "success"


def test_all_producer_wires_endpoint_opt_record() -> None:
    # Regression guard for the PRODUCER: the endpoint-opt convergence gate is
    # dead unless (a) `_run_opt_for_state` emits result.json (--out-json), reads
    # the child's convergence, and returns it, and (b) BOTH the TS-only and
    # multisegment branches assemble segment_log["endpoint_opt"] from those bits.
    # A future edit that drops any of these re-opens the false-success path.
    import inspect

    from mlmm.workflows import all as _all_mod
    from mlmm.workflows.all import _run_opt_for_state

    ros_src = inspect.getsource(_run_opt_for_state)
    assert '"--out-json"' in ros_src
    assert "_read_opt_endpoint_converged(opt_dir)" in ros_src
    assert "return g_opt, final_geom_path, endpoint_converged" in ros_src

    all_src = Path(_all_mod.__file__).read_text(encoding="utf-8")
    # Both producer branches (TS-only + multisegment) assemble the record.
    assert all_src.count('segment_log["endpoint_opt"] = {') == 2
    # ...fed by the reader's returned bits, not a hard-coded literal.
    assert all_src.count("_react_opt_conv = _run_opt_for_state(") == 2
    assert all_src.count("_prod_opt_conv = _run_opt_for_state(") == 2


def test_all_pipeline_aggregate_preserves_legacy_severity() -> None:
    # The composed scientific_status is never LESS severe than the legacy axis:
    # a legacy `partial` (e.g. DFT failed) with a fully-converged IRC stays partial.
    from mlmm.workflows.all import _pipeline_aggregate_truth

    summary = {"segments": [{"index": 1, "kind": "seg", "converged": True}]}
    post = [{
        "index": 1,
        "irc": {"usable": True, "reason": "ok"},
        "endpoint_assignment": {"connectivity_validated": True},
        "endpoint_opt": {"reactant_converged": True, "product_converged": True},
    }]
    truth = _pipeline_aggregate_truth(
        summary, post_segments=post, config={"tsopt": True, "dft": True},
        legacy_status="partial", legacy_reasons=["segment 1: DFT failed (TS)"],
    )
    assert truth.scientific_status == "partial"
    assert "segment 1: DFT failed (TS)" in truth.status_reasons


def test_all_pipeline_requires_validated_mep_irc_connectivity() -> None:
    from mlmm.workflows.all import _pipeline_aggregate_truth

    summary = {
        "segments": [
            {"index": 1, "kind": "seg", "converged": True}
        ]
    }
    base = {
        "index": 1,
        "irc": {"usable": True, "reason": "ok"},
        "endpoint_opt": {
            "reactant_converged": True,
            "product_converged": True,
        },
    }

    for verdict, expected_success in ((False, False), (True, True)):
        post = {
            **base,
            "endpoint_assignment": {"connectivity_validated": verdict},
        }
        truth = _pipeline_aggregate_truth(
            summary,
            post_segments=[post],
            config={"tsopt": True},
            legacy_status="success",
        )
        assert (truth.scientific_status == "success") is expected_success


def test_all_pipeline_aggregate_post_missing_fails_closed_when_tsopt_requested() -> None:
    # Shipped-artifact fail-open: an intermediate MEP summary (post_segments not
    # yet assembled) with tsopt requested must NOT be promoted to success on the
    # MEP trajectory's existence alone. The reactive leaf fails closed
    # (post_missing) until IRC/endpoint post-processing actually runs.
    from mlmm.workflows.all import _pipeline_aggregate_truth

    summary = {"segments": [{"index": 1, "kind": "seg", "barrier_kcal": 10.0}]}
    # post_segments=[] (post ran but produced no record for this segment).
    truth = _pipeline_aggregate_truth(
        summary, post_segments=[], config={"tsopt": True}, legacy_status="success",
    )
    assert truth.scientific_status != "success"        # would have been success (fail-open)
    assert any("segment_1" in r for r in truth.status_reasons)
    # post_segments=None (the very first intermediate write, before post-processing)
    # fails closed too.
    truth_none = _pipeline_aggregate_truth(
        summary, post_segments=None, config={"tsopt": True}, legacy_status="success",
    )
    assert truth_none.scientific_status != "success"


def test_all_pipeline_aggregate_no_tsopt_uses_segment_converged() -> None:
    # A path-only final summary must gate on the segment's OWN reported
    # convergence, never a silent default-True.
    from mlmm.workflows.all import _pipeline_aggregate_truth

    cfg = {"tsopt": False}
    # A nonconverged segment must not be a success even with a barrier band.
    nc = {"segments": [{"index": 1, "kind": "seg", "barrier_kcal": 10.0, "converged": False}]}
    assert _pipeline_aggregate_truth(
        nc, post_segments=None, config=cfg, legacy_status="success",
    ).scientific_status != "success"
    # A missing convergence field fails closed (unknown), not success.
    unk = {"segments": [{"index": 1, "kind": "seg", "barrier_kcal": 10.0}]}
    assert _pipeline_aggregate_truth(
        unk, post_segments=None, config=cfg, legacy_status="success",
    ).scientific_status != "success"
    # A genuinely-converged segment IS a success (legacy behavior preserved).
    ok = {"segments": [{"index": 1, "kind": "seg", "barrier_kcal": 10.0, "converged": True}]}
    assert _pipeline_aggregate_truth(
        ok, post_segments=None, config=cfg, legacy_status="success",
    ).scientific_status == "success"


@pytest.mark.parametrize("mep_converged", [False, None])
def test_all_pipeline_post_success_cannot_promote_bad_mep(
    mep_converged: bool | None,
) -> None:
    """Successful IRC/endpoints cannot overwrite false/unknown MEP truth."""
    from mlmm.workflows.all import _pipeline_aggregate_truth

    summary = {"segments": [{
        "index": 1, "kind": "seg", "barrier_kcal": 10.0,
        "converged": mep_converged,
    }]}
    post = [{
        "index": 1,
        "irc": {"usable": True, "reason": "ok"},
        "endpoint_opt": {"reactant_converged": True, "product_converged": True},
    }]
    truth = _pipeline_aggregate_truth(
        summary, post_segments=post, config={"tsopt": True},
        legacy_status="success",
    )
    assert truth.scientific_status != "success"


def test_read_path_opt_segment_converged_is_tristate(tmp_path: Path) -> None:
    from mlmm.workflows.all import _read_path_opt_segment_converged

    assert _read_path_opt_segment_converged(tmp_path) is None
    result = tmp_path / "result.json"
    result.write_text(json.dumps({"stage_outcomes": [{"converged": True}]}))
    assert _read_path_opt_segment_converged(tmp_path) is True
    result.write_text(json.dumps({"stage_outcomes": [{"converged": False}]}))
    assert _read_path_opt_segment_converged(tmp_path) is False
    result.write_text("{ not json")
    assert _read_path_opt_segment_converged(tmp_path) is None


def test_all_path_opt_child_emits_machine_result() -> None:
    """The no-refine path-opt producer must enable its result.json contract."""
    from mlmm.workflows import all as all_workflow

    source = Path(all_workflow.__file__).read_text(encoding="utf-8")
    branch = source[source.index("else:\n        # --no-refine-path"):source.index("final_trj = path_dir", source.index("else:\n        # --no-refine-path"))]
    assert 'po_args.append("--out-json")' in branch
    assert "_read_path_opt_segment_converged(seg_out)" in branch
    assert "seg_idx = pair_pos + 1" in branch
    assert 'seg_tag = f"seg_{seg_idx:02d}"' in branch
    assert "enumerate(path_opt_segments, start=1)" in source
