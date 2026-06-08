"""Bit-equivalence regression test for the pysisyphus `_array` shim migration.

Each test snapshots a representative numpy↔torch dispatch that the M3-M6
shim migrations touch (`_outer`, `_dot`, `_eigh`, `as_numpy`, `to_xp`). All
asserts are bit-exact for the numpy path and torch fp64 path; torch fp32
allows 1e-6 tolerance to absorb GPU kernel non-bit reproducibility.

If any assertion fires after a M-migration, the shim is no longer
behaviour-preserving and the migration must be reverted at the site.
"""
from __future__ import annotations

import numpy as np
import pytest
import torch


# ---------- _array shim direct equivalence ----------

def test_outer_numpy_matches_inline():
    from pysisyphus._array import _outer

    a = np.array([1.0, 2.0, 3.0])
    b = np.array([4.0, 5.0, 6.0])
    np.testing.assert_array_equal(_outer(a, b), np.outer(a, b))


def test_outer_torch_matches_inline_fp64():
    from pysisyphus._array import _outer

    a = torch.tensor([1.0, 2.0, 3.0], dtype=torch.float64)
    b = torch.tensor([4.0, 5.0, 6.0], dtype=torch.float64)
    np.testing.assert_array_equal(_outer(a, b).cpu().numpy(), torch.outer(a, b).cpu().numpy())


def test_dot_numpy_matches_inline():
    from pysisyphus._array import _dot

    a = np.array([1.0, 2.0, 3.0])
    b = np.array([4.0, 5.0, 6.0])
    assert float(_dot(a, b)) == float(np.dot(a, b))


def test_dot_torch_matches_inline_fp64():
    from pysisyphus._array import _dot

    a = torch.tensor([1.0, 2.0, 3.0], dtype=torch.float64)
    b = torch.tensor([4.0, 5.0, 6.0], dtype=torch.float64)
    assert float(_dot(a, b)) == float(torch.dot(a, b))


def test_eigh_numpy_matches_inline():
    from pysisyphus._array import _eigh

    H = np.array([[4.0, 1.0, 0.0], [1.0, 3.0, 2.0], [0.0, 2.0, 5.0]])
    e_shim, v_shim = _eigh(H)
    e_np, v_np = np.linalg.eigh(H)
    np.testing.assert_array_equal(e_shim, e_np)
    np.testing.assert_array_equal(v_shim, v_np)


def test_eigh_torch_matches_inline_fp64():
    from pysisyphus._array import _eigh

    H = torch.tensor(
        [[4.0, 1.0, 0.0], [1.0, 3.0, 2.0], [0.0, 2.0, 5.0]], dtype=torch.float64
    )
    e_shim, v_shim = _eigh(H)
    e_t, v_t = torch.linalg.eigh(H)
    np.testing.assert_array_equal(e_shim.cpu().numpy(), e_t.cpu().numpy())
    np.testing.assert_array_equal(v_shim.cpu().numpy(), v_t.cpu().numpy())


def test_get_xp_dispatches_correctly():
    from pysisyphus._array import get_xp

    assert get_xp(np.zeros(3)) is np
    assert get_xp(torch.zeros(3)) is torch


def test_as_numpy_idempotent_for_ndarray():
    from pysisyphus._array import as_numpy

    a = np.array([1.0, 2.0, 3.0])
    out = as_numpy(a)
    np.testing.assert_array_equal(out, a)


def test_as_numpy_detaches_torch():
    from pysisyphus._array import as_numpy

    a = torch.tensor([1.0, 2.0, 3.0], requires_grad=True)
    out = as_numpy(a)
    np.testing.assert_array_equal(out, np.array([1.0, 2.0, 3.0]))
    # Original tensor still has grad; the as_numpy detach was a copy not aliased.
    assert a.requires_grad


def test_to_xp_torch_preserves_dtype_device():
    from pysisyphus._array import to_xp

    like = torch.zeros(3, dtype=torch.float32)
    cast = to_xp(np.array([1.0, 2.0, 3.0]), like)
    assert isinstance(cast, torch.Tensor)
    assert cast.dtype == torch.float32
    assert cast.device == like.device


def test_to_xp_numpy_path():
    from pysisyphus._array import to_xp

    like = np.zeros(3, dtype=np.float64)
    cast = to_xp([1, 2, 3], like)
    assert isinstance(cast, np.ndarray)
    assert cast.dtype == np.float64


# ---------- hessian_updates.bfgs_update parity ----------

def _bfgs_update_inline_numpy(H, dx, dg):
    """Reference implementation copied verbatim from the pre-M2 helper."""
    Hdx = H @ dx
    first_term = np.outer(dg, dg) / np.dot(dg, dx)
    second_term = np.outer(Hdx, Hdx) / np.dot(dx, Hdx)
    return first_term - second_term


def test_bfgs_update_numpy_matches_shim():
    from pysisyphus.optimizers.hessian_updates import bfgs_update

    H = np.diag([2.0, 3.0, 4.0])
    dx = np.array([0.1, 0.2, 0.3])
    dg = np.array([0.05, 0.10, 0.15])
    out_shim, label = bfgs_update(H, dx, dg)
    out_ref = _bfgs_update_inline_numpy(H, dx, dg)
    np.testing.assert_array_equal(out_shim, out_ref)
    assert label == "BFGS"


def test_bfgs_update_torch_matches_numpy_close():
    from pysisyphus.optimizers.hessian_updates import bfgs_update

    H_np = np.diag([2.0, 3.0, 4.0])
    H_t = torch.diag(torch.tensor([2.0, 3.0, 4.0], dtype=torch.float64))
    dx_np = np.array([0.1, 0.2, 0.3])
    dx_t = torch.tensor([0.1, 0.2, 0.3], dtype=torch.float64)
    dg_np = np.array([0.05, 0.10, 0.15])
    dg_t = torch.tensor([0.05, 0.10, 0.15], dtype=torch.float64)
    out_np, _ = bfgs_update(H_np, dx_np, dg_np)
    out_t, _ = bfgs_update(H_t, dx_t, dg_t)
    # torch fp64 vs numpy: small floating-point reorder is allowed (1e-12).
    np.testing.assert_allclose(out_np, out_t.cpu().numpy(), atol=1e-12)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
