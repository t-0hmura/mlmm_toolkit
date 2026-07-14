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
    meta_copy = dict(meta or {})
    if "cart_coords" in meta_copy:
        coords = meta_copy["cart_coords"]
        if isinstance(coords, torch.Tensor):
            coords = coords.detach().cpu().numpy()
        meta_copy["cart_coords"] = np.asarray(coords, dtype=float).copy()
    _cache[key] = {
        "hessian": h_cpu,
        "active_dofs": list(active_dofs) if active_dofs is not None else None,
        "meta": meta_copy,
    }


def load(key: str) -> Optional[Dict[str, Any]]:
    """Return cached entry or *None*."""
    return _cache.get(key)


def matches_cart_coords(
    entry: Dict[str, Any],
    cart_coords,
    *,
    atol: float = 1.0e-5,
) -> bool:
    """Whether an entry belongs to ``cart_coords``.

    Entries without a coordinate fingerprint are not safe for cross-stage
    reuse. Coordinates are in bohr and the tolerance accommodates serialized
    structure round-tripping.
    """
    cached_coords = entry.get("meta", {}).get("cart_coords")
    if cached_coords is None:
        return False
    if isinstance(cart_coords, torch.Tensor):
        cart_coords = cart_coords.detach().cpu().numpy()
    cached_arr = np.asarray(cached_coords, dtype=float).reshape(-1)
    current_arr = np.asarray(cart_coords, dtype=float).reshape(-1)
    return cached_arr.shape == current_arr.shape and bool(
        np.allclose(cached_arr, current_arr, rtol=0.0, atol=float(atol))
    )


def discard(key: str) -> None:
    """Drop one cached entry if present."""
    _cache.pop(key, None)


def clear() -> None:
    """Drop all cached Hessians."""
    _cache.clear()
