"""``all`` keeps post-processing the later segments after one stops before IRC.

Three segments run through the stubs of _all_stub_harness.py (each pair has an
interior highest image). The TS optimization of the segments in ``stopped`` does
not converge, so they stop before IRC.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import click
import pytest
from _all_stub_harness import TS_NOT_CONVERGED, run_all

from mlmm.cli.completion import completion_guard
from mlmm.core import utils
from mlmm.workflows import all as all_workflow

FRAMES = (0.0, 2.0, 0.0, 2.0)
REASON = "ts_optimization_not_converged"
REAL_RUN_CLI_MAIN = all_workflow._run_cli_main


def run_three_segments(tmp_path, monkeypatch, stopped, extra=("--tsopt",), parm=True):
    def stop_selected_segments(_context):
        child = all_workflow._run_cli_main

        def wrapped(name, cli, args, **kwargs):
            code = child(name, cli, args, **kwargs)
            if name == "tsopt":
                out = Path(args[args.index("--out-dir") + 1])
                segment = int(re.search(r"seg_(\d+)", out.as_posix()).group(1))
                if segment in stopped:
                    result = out / "result.json"
                    old = json.loads(result.read_text(encoding="utf-8"))
                    result.write_text(json.dumps({
                        **TS_NOT_CONVERGED,
                        "energy_hartree": old["energy_hartree"], "files": old["files"],
                    }), encoding="utf-8")
            return code

        monkeypatch.setattr(all_workflow, "_run_cli_main", wrapped)

    return run_all(
        tmp_path, monkeypatch, frames=FRAMES, extra=extra, parm=parm, wrap=stop_selected_segments)


@pytest.mark.parametrize("position", [1, 2, 3])
def test_a_segment_that_stops_before_irc_does_not_stop_the_others(tmp_path, monkeypatch, position):
    run = run_three_segments(tmp_path, monkeypatch, {position})
    summary = run.summary

    assert run.child_names().count("tsopt") == 3
    assert run.child_names().count("irc") == 2
    assert [item["index"] for item in summary["post_segments"]] == [1, 2, 3]
    assert (run.result.exit_code, summary["execution_status"], summary["scientific_status"]) == (
        0, "completed", "partial")
    assert summary["stopped_before_irc"] is True
    assert summary["pipeline_stop"]["segment"] == position
    assert summary["stopped_segments"] == [
        {"stage": "before_irc", "segment": position, "reason": REASON}]
    log = (run.out / "summary.log").read_text(encoding="utf-8")
    assert f"Stopped before IRC  : segment {position} ({REASON})" in log
    # A segment without an IRC would leave a gap in the one concatenated curve.
    assert not (run.out / "irc_plot_all.png").exists()
    assert not (run.out / "energy_diagram_MLIP_all.png").exists()


def test_every_segment_stopped_is_listed_and_pipeline_stop_is_the_first(tmp_path, monkeypatch):
    run = run_three_segments(tmp_path, monkeypatch, {1, 2, 3})
    summary = run.summary

    assert run.child_names().count("irc") == 0
    assert [item["index"] for item in summary["post_segments"]] == [1, 2, 3]
    assert [item["segment"] for item in summary["stopped_segments"]] == [1, 2, 3]
    assert summary["pipeline_stop"]["segment"] == 1
    assert (run.result.exit_code, summary["execution_status"], summary["scientific_status"]) == (
        0, "completed", "partial")
    assert not (run.out / "irc_plot_all.png").exists()


def test_no_stop_draws_the_irc_plot_over_all_segments(tmp_path, monkeypatch):
    run = run_three_segments(tmp_path, monkeypatch, set())
    summary = run.summary

    assert run.child_names().count("irc") == 3
    assert (run.result.exit_code, summary["scientific_status"]) == (0, "success")
    assert "pipeline_stop" not in summary and "stopped_segments" not in summary
    assert summary["stopped_before_irc"] is False
    assert (run.out / "irc_plot_all.png").is_file()


def test_resume_keeps_only_the_stops_of_the_current_run(tmp_path, monkeypatch):
    from test_all_resume_saved_parm7 import _stub_mm_parm

    _stub_mm_parm(monkeypatch, [])
    first = run_three_segments(tmp_path, monkeypatch, {3}, parm=False)
    assert first.summary["stopped_segments"][0]["segment"] == 3

    resumed = run_three_segments(
        tmp_path, monkeypatch, set(), extra=("--tsopt", "--resume-segment", "2"), parm=False)
    summary = resumed.summary

    assert resumed.child_names().count("tsopt") == 2
    assert summary["stopped_before_irc"] is False
    assert "pipeline_stop" not in summary and "stopped_segments" not in summary
    assert (resumed.result.exit_code, summary["scientific_status"]) == (0, "success")


@pytest.mark.parametrize(
    "resume_segment,stopped_now,expected",
    [
        pytest.param(3, set(), [1], id="earlier-stop-only"),
        pytest.param(2, {3}, [1, 3], id="earlier-and-current-stop"),
    ],
)
def test_resume_keeps_the_stops_of_the_segments_before_it(
    tmp_path, monkeypatch, resume_segment, stopped_now, expected
):
    from test_all_resume_saved_parm7 import _stub_mm_parm

    _stub_mm_parm(monkeypatch, [])
    first = run_three_segments(tmp_path, monkeypatch, {1}, parm=False)
    assert first.summary["stopped_segments"][0]["segment"] == 1

    resumed = run_three_segments(
        tmp_path, monkeypatch, stopped_now,
        extra=("--tsopt", "--resume-segment", str(resume_segment)), parm=False)
    summary = resumed.summary

    assert resumed.child_names().count("tsopt") == 4 - resume_segment
    assert summary["stopped_before_irc"] is True
    assert summary["pipeline_stop"]["segment"] == 1
    assert summary["stopped_segments"] == [
        {"stage": "before_irc", "segment": segment, "reason": REASON} for segment in expected]
    assert (resumed.result.exit_code, summary["scientific_status"]) == (0, "partial")
    log = (resumed.out / "summary.log").read_text(encoding="utf-8")
    stops = "; ".join(f"segment {segment} ({REASON})" for segment in expected)
    assert f"Stopped before IRC  : {stops}" in log


@pytest.mark.parametrize("failing", [{2}, {3}, {1, 2, 3}])
def test_a_failed_final_hessian_stops_the_whole_run(tmp_path, monkeypatch, failing):
    """Unlike a TS that stops before IRC, a failed final TS Hessian ends all with exit 1."""

    def fail_selected_hessians(_context):
        stub = all_workflow._run_cli_main

        def wrapped(name, cli, args, **kwargs):
            if name != "tsopt":
                return stub(name, cli, args, **kwargs)
            out = Path(args[args.index("--out-dir") + 1])
            segment = int(re.search(r"seg_(\d+)", out.as_posix()).group(1))
            if segment not in failing:
                return stub(name, cli, args, **kwargs)

            @click.command(context_settings={"ignore_unknown_options": True, "allow_extra_args": True})
            def tsopt_with_failed_hessian():
                stub(name, cli, args, **kwargs)  # writes the converged result and geometry
                old = json.loads((out / "result.json").read_text(encoding="utf-8"))
                payload = {k: v for k, v in old.items() if k not in {"execution_status", "scientific_status"}}
                payload.update(
                    hessian_status="failed", hessian_error="RuntimeError: boom",
                    saddle_validation="unavailable", n_imaginary_modes=None, n_negative_modes=None,
                    imaginary_frequencies_cm=[], reaction_mode_index=None,
                    reaction_mode_frequency_cm=None, reaction_mode_overlap=None,
                    reaction_mode_source=None,
                )
                utils.write_result_json(out, payload, command="tsopt")

            # The child verdict goes through the real dispatch of all.
            tsopt_with_failed_hessian.callback = completion_guard(tsopt_with_failed_hessian.callback)
            return REAL_RUN_CLI_MAIN(name, tsopt_with_failed_hessian, args, **kwargs)

        monkeypatch.setattr(all_workflow, "_run_cli_main", wrapped)

    run = run_all(tmp_path, monkeypatch, frames=FRAMES, extra=("--tsopt",), wrap=fail_selected_hessians)
    summary = run.summary
    first = min(failing)

    assert run.result.exit_code == 1
    assert "Error: [tsopt] tsopt exit code 1." in run.result.output
    assert run.child_names().count("tsopt") == first
    assert run.child_names().count("irc") == first - 1
    assert (summary["execution_status"], summary["scientific_status"]) == ("failed", "partial")
    assert "pipeline_stop" not in summary


def test_resume_after_a_stop_in_post_processing_keeps_the_finished_segments(tmp_path, monkeypatch):
    """A run that stops while segment 2 is post-processed can resume from segment 2."""
    from test_all_resume_saved_parm7 import _stub_mm_parm

    def stop_in_segment_2(_context):
        child = all_workflow._run_cli_main

        def wrapped(name, cli, args, **kwargs):
            out = Path(args[args.index("--out-dir") + 1])
            if name == "tsopt" and "seg_02" in out.as_posix():
                raise RuntimeError("walltime reached in segment 2")
            return child(name, cli, args, **kwargs)

        monkeypatch.setattr(all_workflow, "_run_cli_main", wrapped)

    _stub_mm_parm(monkeypatch, [])
    first = run_all(tmp_path, monkeypatch, frames=FRAMES, extra=("--tsopt",), parm=False,
                    wrap=stop_in_segment_2)
    assert first.result.exit_code == 1
    assert [item["index"] for item in first.summary["post_segments"]] == [1]
    assert first.summary["scientific_status"] != "success"

    resumed = run_three_segments(
        tmp_path, monkeypatch, set(), extra=("--tsopt", "--resume-segment", "2"), parm=False)
    summary = resumed.summary

    assert resumed.child_names().count("tsopt") == 2
    assert [item["index"] for item in summary["post_segments"]] == [1, 2, 3]
    assert (resumed.result.exit_code, summary["execution_status"], summary["scientific_status"]) == (
        0, "completed", "success")
