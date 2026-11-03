"""The root summary.json that ``all`` keeps while post-processing runs.

Three segments run through the stubs of _all_stub_harness.py. While a segment
is post-processed, the root lists the MEP outputs at the root; a resumed run
publishes this provisional root with the kept segments before post-processing
starts, and refuses a changed post-processing precision.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest
from _all_stub_harness import Behaviour, mlmm_energy, run_all

from mlmm.workflows import all as all_workflow

FRAMES = (0.0, 2.0, 0.0, 2.0)


def _record_root_at_tsopt(monkeypatch, out: Path, seen: dict):
    """Keep the root summary.json, and which of its listed files exist, as each TS optimization starts."""

    def wrap():
        child = all_workflow._run_cli_main

        def wrapped(name, cli, args, **kwargs):
            if name == "tsopt":
                segment_dir = Path(args[args.index("--out-dir") + 1]).as_posix()
                segment = int(re.search(r"seg_(\d+)", segment_dir).group(1))
                root = json.loads((out / "summary.json").read_text(encoding="utf-8"))
                listed = list(root.get("current_output_paths", [])) + [
                    diagram["image"] for diagram in root.get("energy_diagrams", [])
                ]
                # The listed files as they exist when the TS optimization starts.
                seen[segment] = (root, [p for p in listed if (out / p).is_file()])
            return child(name, cli, args, **kwargs)

        monkeypatch.setattr(all_workflow, "_run_cli_main", wrapped)

    return wrap


def test_refine_path_run_root_lists_the_root_mep_outputs(
    tmp_path, monkeypatch
):
    seen: dict = {}
    segments = [
        {
            "index": index, "tag": f"seg_{index:03d}", "kind": "seg",
            "barrier_kcal": (mlmm_energy(1.0) - mlmm_energy(0.0)) * 627.509474,
            "delta_kcal": (mlmm_energy(2.0) - mlmm_energy(0.0)) * 627.509474,
            "bond_changes": "C1-O1 : formed", "converged": True,
        }
        for index in (1, 2, 3)
    ]
    run = run_all(
        tmp_path, monkeypatch, frames=FRAMES, extra=("--refine-path", "--tsopt"),
        behaviour=Behaviour(path_search_segments=segments),
        wrap=_record_root_at_tsopt(monkeypatch, tmp_path / "out", seen),
    )
    assert run.result.exit_code == 0, run.result.output[-3000:]
    assert sorted(seen) == [1, 2, 3]

    root, existing = seen[1]
    mep_outputs = {"energy_diagram_MEP.png", "mep_plot.png", "mep_trj.xyz", "mep_trj.pdb"}
    assert mep_outputs <= set(root["current_output_paths"])
    assert not any(path.startswith("_work") for path in root["current_output_paths"])
    images = [diagram["image"] for diagram in root["energy_diagrams"]]
    assert images == [str(run.out / "energy_diagram_MEP.png")]
    # Every listed output and the MEP image exist at the root.
    assert set(existing) == set(root["current_output_paths"]) | set(images)
    assert root["scientific_status"] != "success"


def _completed_three_segments(tmp_path, monkeypatch, extra=()):
    from test_all_resume_saved_parm7 import _stub_mm_parm

    _stub_mm_parm(monkeypatch, [])
    first = run_all(
        tmp_path, monkeypatch, frames=FRAMES, extra=("--tsopt", *extra), parm=False
    )
    assert first.result.exit_code == 0, first.result.output[-3000:]
    assert [item["index"] for item in first.summary["post_segments"]] == [1, 2, 3]
    return first


def test_resume_that_stops_before_its_first_segment_keeps_a_provisional_root(
    tmp_path, monkeypatch
):
    first = _completed_three_segments(tmp_path, monkeypatch)
    assert first.summary["scientific_status"] == "success"
    echo_section = all_workflow._echo_section

    def stop_at_post_processing(message, *args, **kwargs):
        if "Stage 4/" in str(message):
            raise RuntimeError("walltime reached before segment 2")
        return echo_section(message, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(all_workflow, "_echo_section", stop_at_post_processing)
        stopped = run_all(
            tmp_path, monkeypatch, frames=FRAMES,
            extra=("--tsopt", "--resume-segment", "2"), parm=False,
        )
    assert stopped.result.exit_code != 0
    assert "tsopt" not in stopped.child_names()
    root = stopped.summary
    assert root["scientific_status"] != "success"
    assert root["resumed_from_segment"] == 2
    assert [item["index"] for item in root["post_segments"]] == [1]

    resumed = run_all(
        tmp_path, monkeypatch, frames=FRAMES,
        extra=("--tsopt", "--resume-segment", "2"), parm=False,
    )
    summary = resumed.summary
    assert resumed.child_names().count("tsopt") == 2
    assert [item["index"] for item in summary["post_segments"]] == [1, 2, 3]
    assert (resumed.result.exit_code, summary["scientific_status"]) == (0, "success")


@pytest.mark.parametrize(
    "resume_segment,precision,refused",
    [("2", "fp32", True), ("2", "fp64", False), ("1", "fp32", False)],
)
def test_resume_refuses_a_changed_post_processing_precision(
    tmp_path, monkeypatch, resume_segment, precision, refused
):
    first = _completed_three_segments(tmp_path, monkeypatch, ("--precision-post", "fp64"))
    assert first.summary["mlip_precision_post"] == "fp64"
    before = (first.out / "summary.json").read_bytes()

    resumed = run_all(
        tmp_path, monkeypatch, frames=FRAMES, parm=False,
        extra=("--tsopt", "--resume-segment", resume_segment, "--precision-post", precision),
    )

    if refused:
        assert resumed.result.exit_code == 2
        assert resumed.calls == []
        assert "post-processing precision changed (fp64 in the saved run, fp32 now)" in (
            resumed.result.output
        )
        assert (first.out / "summary.json").read_bytes() == before
    else:
        assert resumed.result.exit_code == 0, resumed.result.output[-3000:]
        assert resumed.child_names().count("tsopt") == 4 - int(resume_segment)
        assert resumed.summary["mlip_precision_post"] == precision
