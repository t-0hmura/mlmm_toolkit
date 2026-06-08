"""Opt-in strict-determinism setup for the MLIP backends.

Activated by the ``--deterministic`` CLI flag (present on every compute
subcommand) or the ``MLMM_STRICT_DETERMINISTIC=1`` environment variable.
Process-global and idempotent: the first activation patches the
non-deterministic ``Tensor.index_reduce_(reduce="mean")`` op (which has no
deterministic CUDA kernel in torch <= 2.8 and crashes under
``use_deterministic_algorithms``) to a ``scatter_reduce`` detour, then enables
``torch.use_deterministic_algorithms(True)`` + cuDNN determinism + fixed seeds.

Design notes:
- **Patch first.** Native strict mode is *known* to crash on
  ``index_reduce_``; probing it first would only waste a guaranteed failure.
  The shim is applied before strict mode is turned on.
- **Fail loud.** If the patch target is gone under a torch upgrade, or strict
  mode rejects an op with no deterministic kernel, raise — never silently
  degrade to ``warn_only`` (that would claim determinism while leaving
  non-deterministic kernels live).
- **Cost.** Buys bit-identical GPU output across runs, at a throughput penalty
  (deterministic scatter/reduce kernels are slower) and the fragility of a
  runtime monkey-patch. Default OFF: default GPU runs carry ~1e-7 A
  scatter/atomic drift that is chemically negligible.
"""
from __future__ import annotations

import os

_DONE = False
_ORIG_INDEX_REDUCE = None


def is_deterministic_requested() -> bool:
    """True when ``MLMM_STRICT_DETERMINISTIC=1`` (the env-var entry point used by
    CI and the direct Python API; the CLI uses ``--deterministic``)."""
    return os.environ.get("MLMM_STRICT_DETERMINISTIC") == "1"


def is_deterministic_active() -> bool:
    """True once strict-deterministic mode has actually been applied this process
    (set by ``setup_deterministic``; covers both the --deterministic flag and the
    env-var path)."""
    return _DONE


def _index_reduce_mean_deterministic(self, dim, index, source, reduce, include_self=True):
    """Drop-in for ``Tensor.index_reduce_`` routing mean reductions through
    ``scatter_reduce`` (which has a deterministic CUDA kernel)."""
    if reduce != "mean":
        return _ORIG_INDEX_REDUCE(self, dim, index, source, reduce=reduce, include_self=include_self)
    if index.dim() == 1 and source.dim() > 1:
        idx_exp = index.unsqueeze(1).expand_as(source)
    else:
        idx_exp = index
    result = self.scatter_reduce(dim, idx_exp, source, reduce="mean", include_self=include_self)
    self.copy_(result)
    return self


def setup_deterministic() -> None:
    """Enable strict-deterministic mode once (idempotent). Patch first, fail loud.

    Raises ``RuntimeError`` if determinism cannot be honoured on this torch
    build (patch target missing, shim self-check mismatch, or strict mode
    rejected) rather than degrading silently.
    """
    global _DONE, _ORIG_INDEX_REDUCE
    if _DONE:
        return
    import torch

    # Patch first (native index_reduce_ mean has no deterministic CUDA kernel).
    orig = getattr(torch.Tensor, "index_reduce_", None)
    if orig is None:
        raise RuntimeError(
            "--deterministic: torch.Tensor.index_reduce_ is missing on this "
            f"torch ({torch.__version__}); the determinism shim needs updating."
        )
    if _ORIG_INDEX_REDUCE is None:
        _ORIG_INDEX_REDUCE = orig
        # One-time equivalence self-check (catch a silent scatter_reduce vs
        # index_reduce_ mean-semantics drift across torch versions). CPU-only;
        # the native op works on CPU even when its CUDA kernel is missing.
        idx = torch.tensor([0, 0, 1])
        src = torch.tensor([1.0, 3.0, 5.0])
        ref = torch.zeros(3).index_reduce_(0, idx, src, reduce="mean", include_self=True)
        got = _index_reduce_mean_deterministic(torch.zeros(3), 0, idx, src, "mean", True)
        if not torch.allclose(ref, got):
            raise RuntimeError(
                "--deterministic: scatter_reduce detour diverges from native "
                "index_reduce_(mean); the shim is unsafe on this torch build."
            )
        torch.Tensor.index_reduce_ = _index_reduce_mean_deterministic

    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    try:
        torch.use_deterministic_algorithms(True, warn_only=False)
    except Exception as e:
        raise RuntimeError(
            f"--deterministic: torch could not enable strict deterministic "
            f"algorithms on this build ({e}). An op lacks a deterministic "
            f"kernel; run in the default (non-deterministic) mode or use a "
            f"torch build that provides it."
        )
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.manual_seed(0)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(0)
    _DONE = True
