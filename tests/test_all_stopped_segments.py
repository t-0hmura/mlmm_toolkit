"""``all`` keeps post-processing the later segments after one stops before IRC.

Three segments run through the stubs of _all_stub_harness.py (each pair has an
interior highest image). The TS optimization of the segments in ``stopped`` does
not converge, so they stop before IRC.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest
from _all_stub_harness import TS_NOT_CONVERGED, run_all

from mlmm.workflows import all as all_workflow

FRAMES = (0.0, 2.0, 0.0, 2.0)
REASON = "ts_optimization_not_converged"


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
