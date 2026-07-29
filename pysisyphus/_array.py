"""Tiny array-namespace shim for the bundled pysisyphus fork.

The bundled pysisyphus carries a deliberate torch-vs-numpy dispatch in the
hot Hessian path so the same code runs both as a numpy CPU optimiser and as a
torch GPU optimiser. This shim centralizes array-backend selection and matched
primitive operations so optimizer and IRC paths do not duplicate inline
``isinstance`` branches.

Behaviour-preserving by construction: torch path and numpy path each call
the same operation (`torch.outer` / `np.outer` etc.) that the inline branch
already called, with the same dtype/device coercion via
`torch.as_tensor(b, dtype=a.dtype, device=a.device)`.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import torch

ArrayLike = Any  # numpy.ndarray | torch.Tensor (typed loosely; numpy stubs vary)


def get_xp(x: ArrayLike):
    """Return the array module backing *x* (torch for tensors, numpy otherwise).

    Use as ``xp = get_xp(H); xp.linalg.eigh(H)`` in place of
    ``if isinstance(H, torch.Tensor): torch.linalg.eigh(H) else np.linalg.eigh(H)``.
    """
    return torch if isinstance(x, torch.Tensor) else np


def to_xp(x: ArrayLike, like: ArrayLike) -> ArrayLike:
    """Cast *x* into the same backend / dtype / device as *like*, no copy if matching."""
    if isinstance(like, torch.Tensor):
        return torch.as_tensor(x, dtype=like.dtype, device=like.device)
    return np.asarray(x, dtype=getattr(like, "dtype", None))


def as_numpy(x: ArrayLike) -> np.ndarray:
    """Return a numpy view / copy of *x*, detaching torch tensors as needed."""
    if isinstance(x, torch.Tensor):
        return x.detach().cpu().numpy()
    return np.asarray(x)


def _outer(a: ArrayLike, b: ArrayLike) -> ArrayLike:
    """Outer product preserving the backend / dtype / device of *a*."""
    if isinstance(a, torch.Tensor):
        b = torch.as_tensor(b, dtype=a.dtype, device=a.device)
        return torch.outer(a, b)
    return np.outer(a, b)


def _dot(a: ArrayLike, b: ArrayLike) -> ArrayLike:
    """Dot product preserving the backend / dtype / device of *a*."""
    if isinstance(a, torch.Tensor):
        b = torch.as_tensor(b, dtype=a.dtype, device=a.device)
        return torch.dot(a, b)
    return np.dot(a, b)


def active_square(H: ArrayLike, idx, *, row_chunk_bytes: int = 2 * 1024 * 1024) -> ArrayLike:
    """Extract the active square sub-block ``H[idx][:, idx]`` with bounded peak.

    Equivalent to ``H[np.ix_(idx, idx)]`` (numpy) or
    ``H.index_select(0, idx).index_select(1, idx)`` (torch), but for torch the
    chained form first materialises a full ``(len(idx), n)`` row temporary. For
    a large Hessian that intermediate can dominate peak VRAM. This routine
    preallocates the ``(m, m)`` output and fills it in bounded row chunks, so
    the transient peak is ``output + one (chunk, n)`` block instead of
    ``output + (m, n)``. It is a pure gather (no arithmetic) and therefore
    bit-identical to the chained result. ``row_chunk_bytes`` is a per-block byte
    budget; the row count is derived from it and clamped to ``[1, m]``.

    numpy inputs keep the simple ``np.ix_`` path (numpy fancy indexing does not
    build the same oversized intermediate).
    """
    if not isinstance(H, torch.Tensor):
        idx_np = np.asarray(idx)
        return H[np.ix_(idx_np, idx_np)]

    if not isinstance(idx, torch.Tensor):
        idx = torch.as_tensor(idx, dtype=torch.long, device=H.device)
    else:
        idx = idx.to(device=H.device, dtype=torch.long)

    m = int(idx.numel())
    n = int(H.shape[1])
    out = torch.empty((m, m), dtype=H.dtype, device=H.device)
    if m == 0:
        return out
    itemsize = out.element_size()
    rows = max(1, int(row_chunk_bytes // max(1, n * itemsize)))
    rows = min(rows, m)
    for start in range(0, m, rows):
        stop = min(start + rows, m)
        block = H.index_select(0, idx[start:stop])  # (chunk, n)
        out[start:stop] = block.index_select(1, idx)  # (chunk, m)
    return out


def _eigh(H: ArrayLike):
    """Symmetric eigendecomposition preserving the backend of *H*.

    Returns ``(eigvals, eigvecs)`` in ascending order, matching both
    ``np.linalg.eigh`` and ``torch.linalg.eigh`` conventions.
    """
    if isinstance(H, torch.Tensor):
        return torch.linalg.eigh(H)
    return np.linalg.eigh(H)


__all__ = ["get_xp", "to_xp", "as_numpy", "_outer", "_dot", "_eigh", "active_square"]
