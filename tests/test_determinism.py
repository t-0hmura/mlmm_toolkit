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


def test_strict_setup_seeds_python_numpy_and_torch(monkeypatch) -> None:
    import numpy as np

    from mlmm.backends import _determinism

    calls: list[tuple[str, int]] = []
    monkeypatch.setattr(_determinism, "_DONE", False)
    monkeypatch.setattr(
        _determinism,
        "_ORIG_INDEX_REDUCE",
        torch.Tensor.index_reduce_,
    )
    monkeypatch.setattr(
        _determinism.random,
        "seed",
        lambda value: calls.append(("python", value)),
    )
    monkeypatch.setattr(
        np.random,
        "seed",
        lambda value: calls.append(("numpy", value)),
    )
    monkeypatch.setattr(
        torch,
        "manual_seed",
        lambda value: calls.append(("torch", value)),
    )
    monkeypatch.setattr(torch, "use_deterministic_algorithms", lambda *_a, **_k: None)
    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)

    _determinism.setup_deterministic()

    assert calls == [("python", 0), ("numpy", 0), ("torch", 0)]
