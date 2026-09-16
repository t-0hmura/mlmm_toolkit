"""Completed optimizations own the workflow result; IRC stops are diagnostics."""
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from mlmm.workflows import all as workflow
from mlmm.workflows._outcomes import validate_irc_samples
from mlmm.workflows import path_search


def _fixture():
    summary = {"segments": [{"index": 1, "kind": "seg", "converged": True}],
               "energy_diagrams": [{"name": "MLIP"}]}
    post = [{"index": 1, "mlip": {}, "irc_traj": "irc.xyz",
             "tsopt": {"optimization_status": "converged", "continue_irc": True,
                       "n_negative_modes": 2, "n_imaginary_modes": 1},
             "ts_imag": {"n_imag": 1},
             "irc": {"usable": False, "reason": "integration_failed",
                     "backward_integration_stop_reason": "Predictor integration exhausted."},
             "endpoint_opt": {"reactant_converged": True, "product_converged": True,
                              "connectivity_validated": False}}]
    return summary, post


@pytest.mark.parametrize("count", [None, 0, 1, 2])
def test_actual_aggregate_uses_numerical_optimizations(count):
    summary, post = _fixture()
    post[0]["tsopt"]["n_imaginary_modes"] = count
    post[0]["ts_imag"]["n_imag"] = count
    config = {"tsopt": True}
    status, reasons = workflow._derive_pipeline_status(summary, post_segments=post, config=config)
    truth = workflow._pipeline_aggregate_truth(summary, post_segments=post, config=config,
                                              legacy_status=status, legacy_reasons=reasons)
    assert truth.scientific_status == "success"
    assert post[0]["irc"]["backward_integration_stop_reason"]
    assert post[0]["endpoint_opt"]["connectivity_validated"] is False
    post[0]["endpoint_opt"]["product_converged"] = False
    assert workflow._pipeline_aggregate_truth(summary, post_segments=post, config=config,
                                               legacy_status=status).scientific_status != "success"
    post[0]["endpoint_opt"]["product_converged"] = True
    post[0]["tsopt"]["optimization_status"] = "not_converged"
    assert workflow._pipeline_aggregate_truth(summary, post_segments=post, config=config,
                                               legacy_status=status).scientific_status != "success"


def test_irc_reader_retains_stop_reason_without_success_decision(tmp_path):
    payload = {"status": "completed", "backward_requested": True,
               "backward_integration_stop_reason": "Predictor integration exhausted.",
               "backward_integration_converged": False,
               "files": {"finished_irc": "finished_irc_trj.xyz"}}
    (tmp_path / "result.json").write_text(json.dumps(payload))
    result = workflow._read_irc_outcome(tmp_path)
    assert result["backward_integration_stop_reason"] == payload["backward_integration_stop_reason"]
    assert not {"usable", "scientific_status", "forward_status", "backward_status"} & result.keys()
    (tmp_path / "result.json").unlink()
    assert "metadata_error" in workflow._read_irc_outcome(tmp_path)


def test_irc_sample_check_accepts_finite_retained_candidates():
    validate_irc_samples([("backward", True, [-1.0, -1.1], [[0., 0., 0.], [0.1, 0., 0.]]),
                          ("forward", False, [], [])])


@pytest.mark.parametrize("energies,coords", [([], []), ([float("nan")], [[0.,0.,0.]]),
                                            ([-1.], [[float("nan"),0.,0.]]), ([-1.,-2.], [[0.,0.,0.]])])
def test_irc_sample_check_does_not_hide_unusable_data(energies, coords):
    with pytest.raises(ValueError):
        validate_irc_samples([("backward", True, energies, coords)])


def test_raw_interval_is_retained_from_path_owner_into_parent():
    segment = path_search.SegmentReport("seg_001", 1.0, 0.0, "Bond formed", seg_index=1, converged=True)
    raw = path_search._raw_path_outcome("raw_seg_002", engine_converged=False)
    leaves, expected = path_search._path_leaves_and_expected([segment], required_outcomes=[raw])
    summary = {"segments": [{"index":1,"kind":"seg","converged":True}],
               "energy_diagrams":[{}], "stage_outcomes":[leaf.to_dict() for leaf in leaves]}
    truth = workflow._pipeline_aggregate_truth(summary, post_segments=[], config={"tsopt":False},
                                              legacy_status="success")
    assert "raw_seg_002" in truth.expected_item_ids
    assert "raw_seg_002" in truth.observed_item_ids
    assert truth.scientific_status != "success"
    summary["stage_outcomes"] = [leaf.to_dict() for leaf in leaves if leaf.item_id != "raw_seg_002"]
    assert workflow._pipeline_aggregate_truth(summary, post_segments=[], config={"tsopt":False},
                                               legacy_status="success").scientific_status == "success"


def test_endpoint_failure_publishes_real_summary_log(tmp_path):
    out = tmp_path / "out"
    tsroot = out / "segments" / "seg_01"
    tsroot.mkdir(parents=True)
    manifest = workflow.InvocationManifest()
    summary = {"status":"failed", "scientific_status":"failed", "execution_status":"completed",
               "pipeline_stop":{"stage":"endpoint_opt","reason":"endpoint_execution_failed"},
               "segments":[{"index":1,"kind":"tsopt"}], "post_segments":[],
               "scientific_status_reasons":["all:segment_1:endpoint_opt:product_converged"]}
    workflow._write_endpoint_failure_summary_log(summary, out_dir=out, tsroot=tsroot,
                                                 manifest=manifest, citation_payload={})
    text = (out / "summary.log").read_text()
    assert "endpoint_opt (endpoint_execution_failed)" in text
    assert "failed" in text
    assert text == (tsroot / "summary.log").read_text()
    assert (out / "summary.log").resolve() in list(manifest.paths("output.public."))


def test_summary_relays_mode_counts_and_numerical_status(tmp_path):
    summary, post = _fixture()
    result = workflow._enrich_summary(summary, version="", pipeline_mode="path-opt",
                                      mlip_backend="uma", charge=0, spin=1, out_dir=tmp_path,
                                      post_segments=post, config={"tsopt":True})
    record = result["post_segments"][0]
    assert record["tsopt"]["n_negative_modes"] == 2
    assert record["ts_imag"]["n_negative_modes"] == 2
    assert record["ts_imag"]["n_imaginary_modes"] == 1
    assert record["ts_imag"]["optimization_status"] == "converged"


@pytest.mark.parametrize("preliminary", [False, None])
def test_final_optimizations_supersede_preliminary_convergence(preliminary):
    summary, post = _fixture()
    summary["segments"][0]["converged"] = preliminary
    summary.update(preopt_requested=True, preopt_converged=preliminary)
    summary["stage_outcomes"] = [{"stage":"path", "item_id":"preopt_endpoint_0",
        "required":True, "executed":True, "converged":preliminary, "usable":False,
        "reason":"not_converged"}]
    config = {"tsopt":True}
    def status(records):
        legacy, reasons = workflow._derive_pipeline_status(summary, post_segments=records, config=config)
        return workflow._pipeline_aggregate_truth(summary, post_segments=records, config=config,
                                                  legacy_status=legacy, legacy_reasons=reasons)
    assert status(post).scientific_status == "success"
    assert summary["preopt_converged"] is preliminary
    assert summary["stage_outcomes"][0]["converged"] is preliminary
    assert status([]).scientific_status != "success"
    post[0]["endpoint_opt"]["product_converged"] = False
    assert status(post).scientific_status != "success"
    post[0]["endpoint_opt"]["product_converged"] = True
    summary["stage_outcomes"].append(path_search._raw_path_outcome(
        "raw_unprocessed", engine_converged=preliminary).to_dict())
    incomplete = status(post)
    assert incomplete.scientific_status != "success"
    assert "raw_unprocessed" in incomplete.expected_item_ids
    assert workflow._pipeline_aggregate_truth(summary, post_segments=None, config={"tsopt":False},
                                               legacy_status="success").scientific_status != "success"


@pytest.mark.parametrize("present", [False, True])
def test_complete_result_without_ts_convergence_signal_is_not_success(present):
    summary, post = _fixture()
    if present:
        post[0]["tsopt"]["optimization_status"] = None
    else:
        del post[0]["tsopt"]["optimization_status"]
    truth = workflow._pipeline_aggregate_truth(
        summary, post_segments=post, config={"tsopt": True}, legacy_status="success",
    )
    assert truth.scientific_status != "success"
    assert any("tsopt" in reason for reason in truth.status_reasons)


@pytest.mark.parametrize("bond_changes", ["", "(no covalent changes detected)"])
def test_final_numerical_completion_is_independent_of_bond_diagnostics(bond_changes):
    summary, post = _fixture()
    summary["segments"][0]["bond_changes"] = bond_changes
    truth = workflow._pipeline_aggregate_truth(
        summary, post_segments=post, config={"tsopt": True}, legacy_status="success",
    )
    assert truth.scientific_status == "success"
    assert truth.expected_item_ids == ("segment_1",)
    post[0]["endpoint_opt"]["product_converged"] = False
    truth = workflow._pipeline_aggregate_truth(
        summary, post_segments=post, config={"tsopt": True}, legacy_status="success",
    )
    assert truth.scientific_status != "success"
