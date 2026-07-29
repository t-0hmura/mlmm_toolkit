"""Fault-injection tests for exact-path result publication."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from mlmm.core import result_commit
from mlmm.core.result_commit import ResultCommitError, RunIdentityError
from mlmm.core.utils import write_result_json


def _seed_pair(root: Path) -> tuple[bytes, bytes]:
    primary = b'{"run_id":"old-primary"}'
    mirror = b'{"run_id":"old-mirror"}'
    (root / "result.json").write_bytes(primary)
    (root / "summary.json").write_bytes(mirror)
    return primary, mirror


def _assert_no_staged_files(root: Path) -> None:
    assert list(root.glob(".*.tmp")) == []


def test_result_pair_serializes_once_and_has_identical_bytes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv(result_commit.MLMM_RUN_ID_ENV, "current")
    original = result_commit.serialize_json_bytes
    calls = 0

    def counted(payload):
        nonlocal calls
        calls += 1
        return original(payload)

    monkeypatch.setattr(result_commit, "serialize_json_bytes", counted)
    path = write_result_json(tmp_path, {"status": "success"}, command="opt")

    assert path == tmp_path / "result.json"
    assert calls == 1
    assert path.read_bytes() == (tmp_path / "summary.json").read_bytes()
    assert json.loads(path.read_text())["run_id"] == "current"
    _assert_no_staged_files(tmp_path)


def test_serialization_failure_preserves_old_pair(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    old_primary, old_mirror = _seed_pair(tmp_path)
    with pytest.raises(ResultCommitError, match="serialize"):
        write_result_json(
            tmp_path,
            {"status": "success", "bad": object()},
            command="opt",
        )
    assert (tmp_path / "result.json").read_bytes() == old_primary
    assert (tmp_path / "summary.json").read_bytes() == old_mirror
    _assert_no_staged_files(tmp_path)


@pytest.mark.parametrize("failed_name", ["summary.json", "result.json"])
def test_stage_failure_publishes_neither_destination(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failed_name: str
) -> None:
    old_primary, old_mirror = _seed_pair(tmp_path)
    real_stage = result_commit.stage_exact

    def fail_selected_stage(path, writer):
        if Path(path).name == failed_name:
            raise ResultCommitError("stage", Path(path), OSError("injected"))
        return real_stage(path, writer)

    monkeypatch.setattr(result_commit, "stage_exact", fail_selected_stage)
    with pytest.raises(ResultCommitError, match="stage"):
        write_result_json(tmp_path, {"status": "success"}, command="opt")
    assert (tmp_path / "result.json").read_bytes() == old_primary
    assert (tmp_path / "summary.json").read_bytes() == old_mirror
    _assert_no_staged_files(tmp_path)


def test_duplicate_mirrors_are_deduplicated_and_primary_publishes_last(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    primary = tmp_path / "result.json"
    mirror = tmp_path / "summary.json"
    real_replace = os.replace
    published = []

    def record_replace(source, destination):
        published.append(Path(destination))
        return real_replace(source, destination)

    monkeypatch.setattr(result_commit.os, "replace", record_replace)
    result_commit.commit_exact_bytes(
        primary,
        b'{"run_id":"current"}',
        mirrors=(primary, mirror, mirror),
    )
    assert published == [mirror, primary]
    assert primary.read_bytes() == mirror.read_bytes()


def test_mirror_publish_failure_preserves_old_pair(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    old_primary, old_mirror = _seed_pair(tmp_path)
    real_replace = os.replace

    def fail_mirror(source, destination):
        if Path(destination).name == "summary.json":
            raise OSError("injected mirror failure")
        return real_replace(source, destination)

    monkeypatch.setattr(result_commit.os, "replace", fail_mirror)
    with pytest.raises(ResultCommitError, match="summary.json"):
        write_result_json(tmp_path, {"status": "success"}, command="opt")
    assert (tmp_path / "result.json").read_bytes() == old_primary
    assert (tmp_path / "summary.json").read_bytes() == old_mirror
    _assert_no_staged_files(tmp_path)


def test_primary_publish_failure_leaves_two_valid_generations(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    old_primary, _ = _seed_pair(tmp_path)
    monkeypatch.setenv(result_commit.MLMM_RUN_ID_ENV, "current")
    real_replace = os.replace

    def fail_primary(source, destination):
        if Path(destination).name == "result.json":
            raise OSError("injected primary failure")
        return real_replace(source, destination)

    monkeypatch.setattr(result_commit.os, "replace", fail_primary)
    with pytest.raises(ResultCommitError, match="result.json"):
        write_result_json(tmp_path, {"status": "success"}, command="opt")
    assert (tmp_path / "result.json").read_bytes() == old_primary
    assert json.loads((tmp_path / "result.json").read_text())["run_id"] == "old-primary"
    assert json.loads((tmp_path / "summary.json").read_text())["run_id"] == "current"
    _assert_no_staged_files(tmp_path)


def test_one_file_replace_failure_preserves_old_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    destination = tmp_path / "summary.json"
    old = b'{"run_id":"old"}'
    destination.write_bytes(old)

    def fail_replace(source, target):
        raise OSError("injected")

    monkeypatch.setattr(result_commit.os, "replace", fail_replace)
    with pytest.raises(ResultCommitError, match="publish"):
        write_result_json(
            tmp_path,
            {"status": "success"},
            command="all",
            filename="summary.json",
            also_write_summary_json=False,
        )
    assert destination.read_bytes() == old
    assert not (tmp_path / "result.json").exists()
    _assert_no_staged_files(tmp_path)


def test_conflicting_caller_run_id_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv(result_commit.MLMM_RUN_ID_ENV, "current")
    with pytest.raises(RunIdentityError, match="conflicts"):
        write_result_json(
            tmp_path,
            {"status": "success", "run_id": "different"},
            command="opt",
        )
    assert list(tmp_path.iterdir()) == []


def test_atomic_write_rejects_symlinked_ancestor_without_external_write(
    tmp_path: Path,
) -> None:
    from mlmm.core.result_commit import ResultCommitError, atomic_write_exact

    root = tmp_path / "root"
    external = tmp_path / "external"
    root.mkdir()
    external.mkdir()
    (root / "segments").symlink_to(external, target_is_directory=True)
    destination = root / "segments" / "result.json"

    with pytest.raises(ResultCommitError, match="symlinked ancestor"):
        atomic_write_exact(destination, lambda stream: stream.write(b"new"))

    assert not (external / "result.json").exists()


def test_heterogeneous_primary_failure_removes_new_companion(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    primary = tmp_path / "result.pdb"
    companion = tmp_path / "result.cif"
    primary.write_bytes(b"old-pdb\n")
    real_replace = result_commit._replace_exact
    failed = False

    def fail_primary_once(staged: Path, destination: Path) -> None:
        nonlocal failed
        if destination == primary and not failed:
            failed = True
            raise OSError("injected heterogeneous primary failure")
        real_replace(staged, destination)

    monkeypatch.setattr(result_commit, "_replace_exact", fail_primary_once)
    with pytest.raises(ResultCommitError, match="heterogeneous primary failure"):
        result_commit.commit_payloads(
            primary,
            {primary: b"new-pdb\n", companion: b"new-cif\n"},
        )

    assert primary.read_bytes() == b"old-pdb\n"
    assert not companion.exists()
    _assert_no_staged_files(tmp_path)
