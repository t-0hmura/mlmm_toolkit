"""Unit tests for the MCP runner envelope (SubcmdResult / TypedDict)."""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

from mlmm.core.result_commit import MLMM_RUN_ID_ENV
from mlmm.mcp import _runner
from mlmm.mcp._runner import (
    SubcmdResult,
    SubcmdResultDict,
    MCP_SUBCMD_RESULT_SCHEMA_VERSION,
    MCP_SUBCMD_RESULT_STATUSES,
    run_subcmd,
)


def test_subcmd_result_to_dict_carries_schema_version() -> None:
    assert MCP_SUBCMD_RESULT_SCHEMA_VERSION == "1.1"
    r = SubcmdResult(status="ok", exit_code=0, argv=["mlmm", "opt"])
    d = r.to_dict()
    assert d["schema_version"] == MCP_SUBCMD_RESULT_SCHEMA_VERSION
    assert d["status"] == "ok"
    assert d["exit_code"] == 0
    assert d["argv"] == ["mlmm", "opt"]


def test_subcmd_result_status_enum() -> None:
    assert MCP_SUBCMD_RESULT_STATUSES == (
        "ok",
        "failed",
        "summary_missing",
        "summary_parse_error",
        "summary_run_mismatch",
    )


def test_subcmd_result_dict_keys_match_to_dict() -> None:
    # SubcmdResultDict is a TypedDict; verify the runtime keys are a subset.
    typed_keys = set(SubcmdResultDict.__annotations__.keys())
    r = SubcmdResult(status="ok", exit_code=0)
    runtime_keys = set(r.to_dict().keys())
    assert runtime_keys <= typed_keys, (runtime_keys - typed_keys)


def _completed(argv) -> subprocess.CompletedProcess[str]:
    return subprocess.CompletedProcess(argv, 0, stdout="ok\n", stderr="")


@pytest.mark.parametrize(
    ("stdout", "expected_status", "expected_summary"),
    [
        ('{"changed": true}', "ok", {"changed": True}),
        ("{broken", "summary_parse_error", {}),
        ("[]", "summary_parse_error", {}),
    ],
)
def test_runner_parses_json_only_stdout(
    monkeypatch: pytest.MonkeyPatch,
    stdout: str,
    expected_status: str,
    expected_summary: dict,
) -> None:
    monkeypatch.setattr(
        "mlmm.mcp._runner.subprocess.run",
        lambda argv, **kwargs: subprocess.CompletedProcess(
            argv, 0, stdout=stdout, stderr=""
        ),
    )
    result = run_subcmd(
        ["mlmm", "bond-summary", "--json"],
        parse_stdout_json=True,
    )
    assert result.status == expected_status
    assert result.summary == expected_summary


def test_timeout_terminates_spawned_process_group(tmp_path: Path) -> None:
    if os.name != "posix":
        return

    pid_file = tmp_path / "processes.json"
    script = (
        "import json, os, pathlib, subprocess, sys, time; "
        "child = subprocess.Popen([sys.executable, '-c', "
        "'import time; time.sleep(30)']); "
        "pathlib.Path(sys.argv[1]).write_text("
        "json.dumps({'parent': os.getpid(), 'child': child.pid})); "
        "time.sleep(30)"
    )
    parent_pid = None
    try:
        result = run_subcmd(
            [sys.executable, "-c", script, str(pid_file)],
            timeout=3.0,
        )
        assert result.status == "failed"
        assert result.exit_code == 124
        assert "timeout_seconds" in (result.hint or "")
        payload = json.loads(pid_file.read_text(encoding="utf-8"))
        parent_pid = int(payload["parent"])
        deadline = time.monotonic() + 2.0
        while _runner._process_group_exists(parent_pid) and time.monotonic() < deadline:
            time.sleep(0.02)
        assert not _runner._process_group_exists(parent_pid)
    finally:
        if parent_pid is not None and _runner._process_group_exists(parent_pid):
            os.killpg(parent_pid, signal.SIGKILL)


def test_runner_binds_interpreter_source_and_current_leaf_pair(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured = {}
    source_root = str(Path(__file__).resolve().parents[1])

    def fake_run(argv, **kwargs):
        captured.update(argv=list(argv), kwargs=kwargs)
        run_id = kwargs["env"][MLMM_RUN_ID_ENV]
        content = json.dumps({"run_id": run_id, "status": "success"})
        (tmp_path / "summary.json").write_text(content)
        (tmp_path / "result.json").write_text(content)
        return _completed(argv)

    monkeypatch.setattr("mlmm.mcp._runner.subprocess.run", fake_run)
    result = run_subcmd(
        ["mlmm", "opt", "-i", "relative/input.pdb"],
        out_dir=tmp_path,
        env_overrides={"PYTHONPATH": os.pathsep.join([source_root, "/elsewhere", source_root])},
    )

    assert result.status == "ok"
    assert result.summary["run_id"] == result.run_id
    assert result.argv == captured["argv"]
    assert result.argv[:3] == [sys.executable, "-m", "mlmm"]
    assert result.argv[-2:] == ["-i", "relative/input.pdb"]
    assert "cwd" not in captured["kwargs"]
    pythonpath = captured["kwargs"]["env"]["PYTHONPATH"].split(os.pathsep)
    assert pythonpath[0] == source_root
    assert pythonpath.count(source_root) == 1


def test_runner_rejects_stale_summary_from_successful_no_write_child(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    stale = json.dumps({"run_id": "old", "status": "success"})
    (tmp_path / "summary.json").write_text(stale)
    (tmp_path / "result.json").write_text(stale)
    monkeypatch.setattr(
        "mlmm.mcp._runner.subprocess.run",
        lambda argv, **kwargs: _completed(argv),
    )

    result = run_subcmd(["mlmm", "opt"], out_dir=tmp_path)
    assert result.status == "summary_run_mismatch"
    assert result.summary == {}


def test_runner_rejects_mixed_leaf_generation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    (tmp_path / "result.json").write_text(json.dumps({"run_id": "old"}))

    def fake_run(argv, **kwargs):
        run_id = kwargs["env"][MLMM_RUN_ID_ENV]
        (tmp_path / "summary.json").write_text(json.dumps({"run_id": run_id}))
        return _completed(argv)

    monkeypatch.setattr("mlmm.mcp._runner.subprocess.run", fake_run)
    result = run_subcmd(["mlmm", "freq"], out_dir=tmp_path)
    assert result.status == "summary_run_mismatch"
    assert result.summary == {}


def test_runner_rejects_nonidentical_leaf_pair_bytes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def fake_run(argv, **kwargs):
        run_id = kwargs["env"][MLMM_RUN_ID_ENV]
        (tmp_path / "summary.json").write_text(
            json.dumps({"run_id": run_id, "status": "success"})
        )
        (tmp_path / "result.json").write_text(
            json.dumps({"status": "success", "run_id": run_id}, indent=2)
        )
        return _completed(argv)

    monkeypatch.setattr("mlmm.mcp._runner.subprocess.run", fake_run)
    result = run_subcmd(["mlmm", "opt"], out_dir=tmp_path)
    assert result.status == "summary_run_mismatch"
    assert result.summary == {}


def test_runner_accepts_current_aggregate_with_unrelated_stale_result(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    (tmp_path / "result.json").write_text(json.dumps({"run_id": "old"}))

    def fake_run(argv, **kwargs):
        run_id = kwargs["env"][MLMM_RUN_ID_ENV]
        (tmp_path / "summary.json").write_text(
            json.dumps({"run_id": run_id, "status": "success"})
        )
        return _completed(argv)

    monkeypatch.setattr("mlmm.mcp._runner.subprocess.run", fake_run)
    result = run_subcmd(["mlmm", "all"], out_dir=tmp_path)
    assert result.status == "ok"
    assert result.summary["run_id"] == result.run_id


@pytest.mark.parametrize(
    ("seed", "expected"),
    [
        (None, "summary_missing"),
        ("{broken", "summary_parse_error"),
    ],
)
def test_runner_rejects_missing_or_malformed_current_summary(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    seed: str | None,
    expected: str,
) -> None:
    if seed is not None:
        (tmp_path / "summary.json").write_text(seed)
    monkeypatch.setattr(
        "mlmm.mcp._runner.subprocess.run",
        lambda argv, **kwargs: _completed(argv),
    )
    result = run_subcmd(["mlmm", "all"], out_dir=tmp_path)
    assert result.status == expected
    assert result.summary == {}


def test_path_shadow_cannot_replace_imported_mlmm_cli(tmp_path: Path) -> None:
    shadow = tmp_path / "mlmm"
    sentinel = tmp_path / "shadow-ran"
    shadow.write_text(f"#!/bin/sh\ntouch '{sentinel}'\nexit 0\n")
    shadow.chmod(0o755)

    result = run_subcmd(
        ["mlmm", "--help"],
        env_overrides={"PATH": str(tmp_path) + os.pathsep + os.environ.get("PATH", "")},
    )
    assert result.exit_code == 0
    assert result.argv[:3] == [sys.executable, "-m", "mlmm"]
    assert not sentinel.exists()
