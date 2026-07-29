from __future__ import annotations

import pytest
import torch


def test_failed_shim_self_check_does_not_commit_state(monkeypatch) -> None:
    from mlmm.backends import _determinism

    original = torch.Tensor.index_reduce_
    monkeypatch.setattr(_determinism, "_DONE", False)
    monkeypatch.setattr(_determinism, "_ORIG_INDEX_REDUCE", None)
    monkeypatch.setattr(torch, "allclose", lambda *_args, **_kwargs: False)

    with pytest.raises(RuntimeError, match="shim is unsafe"):
        _determinism.setup_deterministic()

    assert _determinism._ORIG_INDEX_REDUCE is None
    assert _determinism.is_deterministic_active() is False
    assert torch.Tensor.index_reduce_ is original
