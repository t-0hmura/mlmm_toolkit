"""Focused ownership tests for the composite-workflow run session."""

from __future__ import annotations

import os
from pathlib import Path
import gc
import hashlib
import weakref

import pytest

from mlmm.workflows._run_session import (
    ArtifactClaimError,
    CalculatorLease,
    InvocationManifest,
    InvocationResources,
    RunSession,
    declare_public_output,
    current_output_paths,
    refresh_current_public_outputs,
)
from mlmm.core.result_commit import (
    MLMM_RUN_ID_ENV,
    with_current_run_id,
)


def test_manifest_rejects_unchanged_stale_artifact(tmp_path: Path) -> None:
    stale = tmp_path / "summary.json"
    stale.write_text("{}", encoding="utf-8")
    manifest = InvocationManifest()
    manifest.declare("path.summary", [stale])

    with pytest.raises(ArtifactClaimError, match="did not produce"):
        manifest.claim_one("path.summary")
    assert manifest.paths("path.") == []


def test_manifest_accepts_same_byte_rewrite_by_stat_identity(tmp_path: Path) -> None:
    artifact = tmp_path / "hei_seg_01.xyz"
    artifact.write_bytes(b"same bytes\n")
    manifest = InvocationManifest()
    manifest.declare("path.hei.01", [artifact])

    replacement = tmp_path / ".replacement"
    replacement.write_bytes(artifact.read_bytes())
    os.replace(replacement, artifact)

    assert manifest.claim_one("path.hei.01") == artifact.resolve()
    stamp = manifest.produced["path.hei.01"][1]
    assert stamp.sha256 is not None
    assert len(stamp.sha256) == 64


def test_manifest_does_not_hash_unchanged_stale_artifact(
    tmp_path: Path, monkeypatch,
) -> None:
    from mlmm.workflows import _run_session

    stale = tmp_path / "large-stale.bin"
    stale.write_bytes(b"stale")
    manifest = InvocationManifest()
    manifest.declare("path.stale", [stale])
    monkeypatch.setattr(
        _run_session,
        "sha256",
        lambda: (_ for _ in ()).throw(AssertionError("stale file was hashed")),
    )

    assert manifest.claim_optional("path.stale") is None


def test_manifest_reuses_digest_for_unchanged_claimed_artifact(
    tmp_path: Path, monkeypatch,
) -> None:
    from mlmm.workflows import _run_session

    artifact = tmp_path / "current.bin"
    manifest = InvocationManifest()
    manifest.declare("path.current", [artifact])
    artifact.write_bytes(b"current")
    digest_calls = 0

    def counted_sha256():
        nonlocal digest_calls
        digest_calls += 1
        return hashlib.sha256()

    monkeypatch.setattr(_run_session, "sha256", counted_sha256)

    assert manifest.claim_one("path.current") == artifact.resolve()
    assert manifest.claim_optional("path.current") == artifact.resolve()
    assert digest_calls == 1


def test_manifest_rejects_symlink_at_exact_declared_path(tmp_path: Path) -> None:
    target = tmp_path / "external.json"
    target.write_text("{}", encoding="utf-8")
    artifact = tmp_path / "summary.json"
    manifest = InvocationManifest()
    manifest.declare("path.summary", [artifact])
    artifact.symlink_to(target)

    with pytest.raises(ArtifactClaimError, match="did not produce"):
        manifest.claim_one("path.summary")

    artifact.unlink()
    artifact.write_text("{}", encoding="utf-8")
    assert manifest.claim_one("path.summary") == artifact.absolute()


def test_public_declaration_rejects_symlinked_ancestor(tmp_path: Path) -> None:
    root = tmp_path / "out"
    external = tmp_path / "external"
    root.mkdir()
    external.mkdir()
    (root / "segments").symlink_to(external, target_is_directory=True)

    with pytest.raises(ValueError, match="symlinked ancestor"):
        declare_public_output(
            InvocationManifest(),
            root,
            root / "segments" / "seg_01" / "result.json",
        )


def test_invocation_resources_reject_concurrent_output_lock(tmp_path: Path) -> None:
    first = InvocationResources()
    second = InvocationResources()
    lock_path = tmp_path / "_work" / ".run.lock"
    first.own_exclusive_lock(lock_path)

    with pytest.raises(ArtifactClaimError, match="already using output directory"):
        second.own_exclusive_lock(lock_path)

    first.close()
    second.own_exclusive_lock(lock_path)
    second.close()


def test_public_refresh_requires_producer_declaration(tmp_path: Path) -> None:
    root = tmp_path / "out"
    root.mkdir()
    undeclared = root / "concurrent.txt"
    declared = root / "summary.json"
    manifest = InvocationManifest()
    declare_public_output(manifest, root, declared)
    undeclared.write_text("unrelated", encoding="utf-8")
    declared.write_text("{}", encoding="utf-8")

    current = refresh_current_public_outputs(manifest, root)

    assert current == [declared.absolute()]
    assert undeclared.absolute() not in current


def test_current_output_paths_are_sorted_relative_and_current(tmp_path: Path) -> None:
    root = tmp_path / "out"
    current = root / "segments" / "seg_01" / "structures" / "ts.pdb"
    stale = root / "summary.log"
    for path in (current, stale):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("old", encoding="utf-8")
    manifest = InvocationManifest()
    declare_public_output(manifest, root, current)
    declare_public_output(manifest, root, stale)
    replacement = root / ".current"
    replacement.write_text("current", encoding="utf-8")
    replacement.replace(current)

    assert current_output_paths(manifest, root) == [
        "segments/seg_01/structures/ts.pdb"
    ]


def test_manifest_exposes_only_claimed_current_paths(tmp_path: Path) -> None:
    current = tmp_path / "stage_01.xyz"
    stale = tmp_path / "stage_02.xyz"
    stale.write_text("old", encoding="utf-8")
    manifest = InvocationManifest()
    manifest.declare("scan.stage.01", [current])
    manifest.declare("scan.stage.02", [stale])
    current.write_text("new", encoding="utf-8")

    manifest.claim_one("scan.stage.01")
    assert manifest.claim_optional("scan.stage.02") is None
    assert manifest.paths("scan.stage.") == [current.resolve()]


def test_same_index_stale_hei_fails_until_rewritten(tmp_path: Path) -> None:
    hei = tmp_path / "hei_seg_01.xyz"
    hei.write_text("1\nold\nH 0 0 0\n", encoding="utf-8")
    stale_run = InvocationManifest()
    stale_run.declare("path.hei.01", [hei])
    with pytest.raises(ArtifactClaimError):
        stale_run.claim_one("path.hei.01")

    current_run = InvocationManifest()
    current_run.declare("path.hei.01", [hei])
    replacement = tmp_path / "hei.rewrite"
    replacement.write_bytes(hei.read_bytes())
    os.replace(replacement, hei)
    assert current_run.claim_one("path.hei.01") == hei.resolve()


def test_repeated_invocations_have_distinct_claims_and_final_digests(
    tmp_path: Path,
) -> None:
    artifact = tmp_path / "summary.json"
    run_ids = []
    digests = []
    for payload in (b'{"run": 1}', b'{"run": 2}'):
        manifest = InvocationManifest()
        manifest.declare("path.summary", [artifact])
        replacement = tmp_path / ".summary.next"
        replacement.write_bytes(payload)
        os.replace(replacement, artifact)
        manifest.claim_one("path.summary")
        run_ids.append(manifest.run_id)
        digests.append(manifest.produced["path.summary"][1].sha256)
    assert run_ids[0] != run_ids[1]
    assert digests[0] != digests[1]


def test_unclaimed_partial_file_cannot_cross_invocations(tmp_path: Path) -> None:
    partial = tmp_path / "result.json"
    first = InvocationManifest()
    first.declare("ts.result.01", [partial])
    partial.write_text('{"partial": true}', encoding="utf-8")
    assert first.paths("ts.result.") == []

    second = InvocationManifest()
    second.declare("ts.result.01", [partial])
    with pytest.raises(ArtifactClaimError):
        second.claim_one("ts.result.01")
    assert first.run_id != second.run_id


def test_internal_manifest_is_atomic_and_carries_digest(tmp_path: Path) -> None:
    artifact = tmp_path / "result.xyz"
    manifest = InvocationManifest()
    manifest.declare("scan.stage.01", [artifact])
    artifact.write_text("1\nframe\nH 0 0 0\n", encoding="utf-8")
    manifest.claim_one("scan.stage.01")

    internal = manifest.write_internal(tmp_path / "_work" / "_run_manifest.json")
    payload = internal.read_text(encoding="utf-8")
    assert manifest.run_id in payload
    assert manifest.produced["scan.stage.01"][1].sha256 in payload


def test_invocation_resources_are_idempotent_lifo() -> None:
    calls: list[str] = []
    resources = InvocationResources()
    resources.add(lambda: calls.append("first"))
    resources.add(lambda: calls.append("second"))

    resources.close()
    resources.close()

    assert calls == ["second", "first"]


def test_run_session_binds_generated_id_and_removes_absent_environment(
    monkeypatch,
) -> None:
    monkeypatch.delenv(MLMM_RUN_ID_ENV, raising=False)
    session = RunSession()

    bound = session.own_run_id_environment()

    assert bound == session.manifest.run_id
    assert os.environ[MLMM_RUN_ID_ENV] == session.manifest.run_id
    assert with_current_run_id({})["run_id"] == session.manifest.run_id
    session.close()
    assert MLMM_RUN_ID_ENV not in os.environ


def test_run_session_retains_and_restores_existing_environment(monkeypatch) -> None:
    monkeypatch.setenv(MLMM_RUN_ID_ENV, "mcp-existing")
    session = RunSession(manifest=InvocationManifest(run_id="unused-generated"))

    bound = session.own_run_id_environment()

    assert bound == "mcp-existing"
    assert session.manifest.run_id == "mcp-existing"
    os.environ[MLMM_RUN_ID_ENV] = "child-temporary"
    session.close()
    assert os.environ[MLMM_RUN_ID_ENV] == "mcp-existing"


class _FakeGeometry:
    def __init__(self) -> None:
        self.calculator = None

    def set_calculator(self, calculator) -> None:
        self.calculator = calculator


class _ClosableCalculator:
    def __init__(self) -> None:
        self.close_count = 0

    def close(self) -> None:
        self.close_count += 1


def test_calculator_lease_attaches_exact_object_and_releases_once() -> None:
    calculator = _ClosableCalculator()
    geometries = [_FakeGeometry(), _FakeGeometry(), _FakeGeometry()]
    lease = CalculatorLease(calculator)
    for geometry in geometries:
        lease.attach(geometry)

    assert all(geometry.calculator is calculator for geometry in geometries)
    lease.release()
    lease.release()

    assert all(geometry.calculator is None for geometry in geometries)
    assert calculator.close_count == 1


def test_calculator_lease_drops_its_last_core_reference() -> None:
    geometry = _FakeGeometry()
    calculator = _ClosableCalculator()
    reference = weakref.ref(calculator)
    lease = CalculatorLease(calculator)
    lease.attach(geometry)
    del calculator

    lease.release()
    gc.collect()

    assert geometry.calculator is None
    assert reference() is None
