"""M54: the runtime cache is keyed by topology CONTENT, not the pathname alone.

Replacing a parm7's bytes at the same path (even same-size / same-mtime) must
build a new parsed runtime; byte-identical content may reuse the old one.  The
parse is fully stubbed so the test exercises the cache-key logic itself without
depending on a real Amber prmtop parse.
"""

from __future__ import annotations

import os

import pytest
import torch

from hessian_ff import workflows


class _FakeSystem:
    def __init__(self, natom: int = 3) -> None:
        self.natom = natom

    def to(self, *args, **kwargs):
        return self


class _FakeFF:
    def __init__(self, system, **kwargs) -> None:
        self.system = system


@pytest.fixture()
def counting_runtime(monkeypatch):
    counts = {"system": 0, "ff": 0}

    def fake_load_system(prmtop, device="cpu"):
        counts["system"] += 1
        return _FakeSystem()

    def fake_ff(system, **kwargs):
        counts["ff"] += 1
        return _FakeFF(system, **kwargs)

    def fake_load_coords(coords, natom, device, dtype):
        return torch.zeros((int(natom), 3), dtype=dtype)

    monkeypatch.setattr(workflows, "load_system", fake_load_system)
    monkeypatch.setattr(workflows, "ForceFieldTorch", fake_ff)
    monkeypatch.setattr(workflows, "load_coords", fake_load_coords)
    workflows.clear_runtime_cache()
    yield counts
    workflows.clear_runtime_cache()


def test_unchanged_bytes_parse_once(tmp_path, counting_runtime) -> None:
    p = tmp_path / "system.parm7"
    p.write_bytes(b"PRMTOP-CONTENT-A")

    workflows._load_runtime(p, p, "cpu", True)
    workflows._load_runtime(p, p, "cpu", True)

    assert counting_runtime["system"] == 1
    assert counting_runtime["ff"] == 1


def test_same_size_different_bytes_reparse(tmp_path, counting_runtime) -> None:
    p = tmp_path / "system.parm7"
    p.write_bytes(b"AAAA")
    workflows._load_runtime(p, p, "cpu", True)
    stat_before = p.stat()

    # Replace with different bytes of the SAME size and restore the mtime, so a
    # size/mtime-only key would (wrongly) reuse the old runtime.
    p.write_bytes(b"BBBB")
    os.utime(p, (stat_before.st_atime, stat_before.st_mtime))
    assert p.stat().st_size == stat_before.st_size

    workflows._load_runtime(p, p, "cpu", True)
    assert counting_runtime["system"] == 2


def test_byte_identical_content_reuses(tmp_path, counting_runtime) -> None:
    p = tmp_path / "system.parm7"
    p.write_bytes(b"IDENTICAL-CONTENT")
    workflows._load_runtime(p, p, "cpu", True)

    # Rewrite byte-identical content (new mtime) — the content digest is
    # unchanged, so the runtime is reused.
    p.write_bytes(b"IDENTICAL-CONTENT")
    workflows._load_runtime(p, p, "cpu", True)
    assert counting_runtime["system"] == 1


def test_stale_generation_is_evicted_on_content_change(tmp_path, counting_runtime) -> None:
    p = tmp_path / "system.parm7"
    p.write_bytes(b"GEN-1")
    workflows._load_runtime(p, p, "cpu", True)
    p.write_bytes(b"GEN-2-DIFFERENT")
    workflows._load_runtime(p, p, "cpu", True)

    # Only the current generation for this path survives in the cache.
    resolved = str(p.resolve())
    keys_for_path = [k for k in workflows._RUNTIME_CACHE if k[0] == resolved]
    assert len(keys_for_path) == 1
