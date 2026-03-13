"""In-process Hessian cache for the ``all`` workflow.

Stores Hessian matrices as CPU tensors (or numpy arrays) in a module-level
dict so that subsequent CLI stages executed via ``_run_cli_main()`` can
reuse them without recomputation.

Each entry may carry ``active_dofs`` — a list of DOF indices that the stored
(partial) Hessian spans.  Consumers set ``geometry.within_partial_hessian``
before assigning partial Hessians to ``geometry.cart_hessian``.
"""

import numpy as np
import torch
from typing import Any, Dict, Optional, Sequence

_cache: Dict[str, Any] = {}


def store(
    key: str,
    H,
    active_dofs: Optional[Sequence[int]] = None,
    meta: Optional[dict] = None,
) -> None:
    """Cache a (possibly partial) Hessian on CPU."""
    if isinstance(H, torch.Tensor):
        h_cpu = H.detach().cpu().clone()
    else:
        h_cpu = np.array(H, copy=True)
    _cache[key] = {
        "hessian": h_cpu,
        "active_dofs": list(active_dofs) if active_dofs is not None else None,
        "meta": meta or {},
    }


def load(key: str) -> Optional[Dict[str, Any]]:
    """Return cached entry or *None*."""
    return _cache.get(key)


def clear() -> None:
    """Drop all cached Hessians."""
    _cache.clear()
