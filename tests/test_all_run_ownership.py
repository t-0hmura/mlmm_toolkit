"""Invocation-level ownership falsifiers for ``mlmm all`` (production path)."""

from __future__ import annotations

import os
from pathlib import Path

from click.testing import CliRunner
import pytest

from mlmm.workflows import all as all_workflow
from mlmm.workflows._run_session import InvocationManifest, declare_public_output
from mlmm.core.result_commit import MLMM_RUN_ID_ENV, with_current_run_id


_MIN_PDB = (
    "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
    "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C\n"
    "END\n"
)


class _FakePrepared:
    def __init__(self, path: Path) -> None:
        self.source_path = Path(path)
        self.display_path = Path(path)
        self.structure_template = None

    def cleanup(self) -> None:  # pragma: no cover - trivial
        pass


@pytest.mark.parametrize("failure", [None, RuntimeError("boom")])
@pytest.mark.parametrize("prior_run_id", [None, "prior-mcp-run"])
def test_all_binds_run_id_and_restores_process_state(
    tmp_path: Path, monkeypatch, failure, prior_run_id,
) -> None:
    from mlmm.core import utils as _utils

    inputs = [tmp_path / "left.pdb", tmp_path / "right.pdb"]
    for path in inputs:
        path.write_text(_MIN_PDB, encoding="utf-8")

    if prior_run_id is None:
        monkeypatch.delenv(MLMM_RUN_ID_ENV, raising=False)
    else:
        monkeypatch.setenv(MLMM_RUN_ID_ENV, prior_run_id)

    # A child dispatched mid-run must observe the bound run identity.
    observed: list[tuple[str | None, str | None]] = []

    monkeypatch.setattr(
        all_workflow, "prepare_input_structure", lambda p: _FakePrepared(p)
    )
    monkeypatch.setattr(all_workflow.shutil, "which", lambda _cmd: "/usr/bin/true")

    real_build = all_workflow._build_effective_args_yaml

    def build_or_fail(*args, **kwargs):
        observed.append(
            (os.environ.get(MLMM_RUN_ID_ENV), with_current_run_id({}).get("run_id"))
        )
        if failure is not None:
            raise failure
        return real_build(*args, **kwargs)

    monkeypatch.setattr(all_workflow, "_build_effective_args_yaml", build_or_fail)

    # Seed process state whose restoration must be observable. cli forces
    # pipeline_mode True internally, so seeding it to the module default False
    # makes a dropped restore detectable (restore→False vs cli's forced True);
    # echo._started defaults False, so seeding True proves its restore too.
    prior_pipeline = _utils._PIPELINE_MODE
    prior_started = all_workflow._echo_state._started
    _utils.set_pipeline_mode(False)
    all_workflow._echo_state._started = True
    try:
        result = CliRunner().invoke(
            all_workflow.cli,
            [
                "-i", str(inputs[0]), "-i", str(inputs[1]),
                "-c", "protein",
                "-q", "0",
                "--out-dir", str(tmp_path / "out"),
                "--dry-run",
            ],
        )

        if failure is None:
            assert result.exit_code == 0, result.output
        else:
            assert result.exit_code != 0

        # The run identity was bound for the in-process child boundary.
        assert observed and observed[0][0] and observed[0][1]
        assert observed[0][0] == observed[0][1]
        if prior_run_id is not None:
            assert observed[0][0] == prior_run_id

        # Exact restoration of the process-global state on every exit path.
        assert _utils._PIPELINE_MODE is False
        assert all_workflow._echo_state._started is True
        if prior_run_id is None:
            assert MLMM_RUN_ID_ENV not in os.environ
        else:
            assert os.environ[MLMM_RUN_ID_ENV] == prior_run_id
    finally:
        _utils.set_pipeline_mode(prior_pipeline)
        all_workflow._echo_state._started = prior_started


def test_key_outputs_exclude_stale_and_undeclared(tmp_path: Path) -> None:
    root = tmp_path / "out"
    current_segment = root / "segments" / "seg_01" / "ts.pdb"
    stale_segment = root / "segments" / "seg_02" / "ts.pdb"
    stale_diagram = root / "irc_plot_all.png"
    for path in (current_segment, stale_segment, stale_diagram):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("old", encoding="utf-8")

    manifest = InvocationManifest()
    summary = root / "summary.json"
    declare_public_output(manifest, root, current_segment)
    declare_public_output(manifest, root, summary)

    replacement = root / ".current"
    replacement.write_text("current", encoding="utf-8")
    replacement.replace(current_segment)
    summary.write_text("{}", encoding="utf-8")

    key_outputs = all_workflow._current_key_output_files(manifest, root)

    assert "seg_01" in key_outputs
    assert "ts.pdb" in key_outputs["seg_01"]["files"]
    assert "seg_02" not in key_outputs
    assert "irc_plot_all.png" not in key_outputs
    assert "summary.json" in key_outputs


def test_current_run_segment_diagram_declared_and_included(tmp_path: Path) -> None:
    """A per-segment energy diagram declared+claimed this run reappears in
    key_output_files.

    Guards the regression where the declared-only summary dropped the segment
    energy diagrams: the pipeline now routes them through a declare→write→claim
    helper; this exercises
    that exact sequence for a ``segments/seg_NN/`` PNG destination.
    """
    out_dir = tmp_path / "out"
    diagram = out_dir / "segments" / "seg_01" / "energy_diagram_MLIP.png"
    diagram.parent.mkdir(parents=True)

    manifest = InvocationManifest()
    # Declare before the producer writes so the pre-run baseline is captured.
    all_workflow._declare_public_output(manifest, out_dir, diagram)
    diagram.write_bytes(b"\x89PNG current-run diagram")
    all_workflow._claim_public_output(manifest, out_dir, diagram)

    key_files = all_workflow._current_key_output_files(manifest, out_dir)
    assert "seg_01" in key_files
    assert "energy_diagram_MLIP.png" in key_files["seg_01"]["files"]

    # Falsifier: an undeclared sibling diagram is still excluded.
    undeclared = out_dir / "segments" / "seg_01" / "energy_diagram_DFT.png"
    undeclared.write_bytes(b"\x89PNG undeclared diagram")
    key_files_after = all_workflow._current_key_output_files(manifest, out_dir)
    assert "energy_diagram_DFT.png" not in key_files_after["seg_01"]["files"]


def test_publish_and_seg_copy_track_only_current_run(tmp_path: Path) -> None:
    """The real producer helpers surface current-run outputs only."""

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    # Stale artifacts from a prior invocation that this run never rewrites.
    prior_seg = out_dir / "segments" / "seg_09" / "reactant.xyz"
    prior_seg.parent.mkdir(parents=True)
    prior_seg.write_text("2\nstale\nH 0 0 0\nH 0 0 0.7\n", encoding="utf-8")
    stale_plot = out_dir / "irc_plot_all.png"
    stale_plot.write_text("stale plot", encoding="utf-8")

    manifest = InvocationManifest()
    # Early root-deliverable declaration captures the pre-run baseline.
    for name in ("summary.log", "summary.json", "irc_plot_all.png"):
        all_workflow._declare_public_output(manifest, out_dir, out_dir / name)

    # Produce a current-run segment via the real copy helper.
    src = tmp_path / "reactant.xyz"
    src.write_text("2\ncurrent\nH 0 0 0\nH 0 0 0.75\n", encoding="utf-8")
    all_workflow._copy_structures_to_seg_dir(
        {"R": src}, out_dir, 1, ".xyz", manifest=manifest,
    )

    # Publish the aggregate summary via the real manifest-summary producer.
    all_workflow._publish_manifest_summary(
        out_dir / "summary.json",
        {"out_dir": str(out_dir), "status": "success"},
        manifest=manifest,
        out_dir=out_dir,
    )

    key_files = all_workflow._current_key_output_files(manifest, out_dir)

    assert "summary.json" in key_files
    assert "seg_01" in key_files
    assert "reactant.xyz" in key_files["seg_01"]["files"]
    # Stale prior-run outputs never rewritten this run are excluded.
    assert "seg_09" not in key_files
    assert "irc_plot_all.png" not in key_files

    # The persisted private manifest carries the current run identity + digest.
    internal_path = out_dir / all_workflow.WORK_DIRNAME / "_run_manifest.json"
    import json

    internal = json.loads(internal_path.read_text(encoding="utf-8"))
    assert internal["run_id"] == manifest.run_id
    assert "output.public.summary.json" in internal["produced"]
    assert "output.public.segments/seg_01/reactant.xyz" in internal["produced"]
    assert "output.public.irc_plot_all.png" not in internal["produced"]


def test_finalize_summary_includes_outputs_from_late_producers(tmp_path: Path) -> None:
    """The last summary generation owns logs and products written after JSON."""

    import json

    out_dir = tmp_path / "out"
    work_dir = out_dir / "_work" / "path_search"
    work_dir.mkdir(parents=True)
    manifest = InvocationManifest()
    late_paths = [
        out_dir / "summary.log",
        out_dir / "ml_region.pdb",
        out_dir / "energy_diagram_MLIP_all.png",
    ]
    for path in [out_dir / "summary.json", *late_paths]:
        all_workflow._declare_public_output(manifest, out_dir, path)

    summary = {"out_dir": str(out_dir), "status": "success"}
    all_workflow._publish_manifest_summary(
        out_dir / "summary.json",
        summary,
        manifest=manifest,
        out_dir=out_dir,
        mirrors=(work_dir / "summary.json",),
    )
    for path in late_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"current {path.name}\n", encoding="utf-8")

    all_workflow._finalize_current_summary(
        out_dir / "summary.json",
        summary,
        manifest=manifest,
        out_dir=out_dir,
        mirrors=(work_dir / "summary.json",),
    )

    stored = json.loads((out_dir / "summary.json").read_text(encoding="utf-8"))
    assert stored["current_output_paths"] == [
        "energy_diagram_MLIP_all.png",
        "ml_region.pdb",
        "summary.json",
        "summary.log",
    ]
    assert set(stored["key_output_files"]) == set(stored["current_output_paths"])
    assert json.loads((work_dir / "summary.json").read_text(encoding="utf-8")) == stored


def test_all_finalizes_after_each_late_summary_log_producer() -> None:
    """Every successful terminal path republishes after writing summary.log."""

    import inspect

    source = inspect.getsource(all_workflow.cli.callback)
    assert source.count("_write_pipeline_summary_log([])") == 3
    assert source.count("_write_pipeline_summary_log(post_segment_logs)") == 1
    assert source.count("_finalize_current_summary(") == 5
