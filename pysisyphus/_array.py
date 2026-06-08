"""Tiny array-namespace shim for the bundled pysisyphus fork.

The bundled pysisyphus carries a deliberate torch-vs-numpy dispatch in the
hot Hessian path (so the same code runs both as a numpy CPU optimiser and
as a torch GPU optimiser inside our MLIP workflow). The historical pattern
was an inline `if isinstance(x, torch.Tensor): ... else: ...` repeated at
every site (~130 hits across `optimizers/`, `tsoptimizers/`, `irc/`).

This shim provides one place to ask "what backend is this array on?" and
matched primitive ops that route through the right module without changing
results. Adding more helpers here is preferred over adding more inline
`isinstance` checks.

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


def _eigh(H: ArrayLike):
    """Symmetric eigendecomposition preserving the backend of *H*.

    Returns ``(eigvals, eigvecs)`` in ascending order, matching both
    ``np.linalg.eigh`` and ``torch.linalg.eigh`` conventions.
    """
    if isinstance(H, torch.Tensor):
        return torch.linalg.eigh(H)
    return np.linalg.eigh(H)


__all__ = ["get_xp", "to_xp", "as_numpy", "_outer", "_dot", "_eigh"]
