"""Compact/full active-space projection must be idempotent without moving DOFs."""

from types import SimpleNamespace

import numpy as np
import pytest
import torch

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer


class _Optimizer(HessianOptimizer):
    def optimize(self):
        raise NotImplementedError


@pytest.fixture(params=["numpy", "torch_cpu"])
def convert(request):
    def array(values):
        values = np.asarray(values, dtype=np.float64)
        return torch.as_tensor(values) if request.param == "torch_cpu" else values
    return array


@pytest.fixture(params=["noncontiguous_9", "h2a_1383"])
def basis(request):
    if request.param == "noncontiguous_9":
        return 9, np.array([0, 1, 2, 6, 7, 8])
    frozen = {0, 35, 75, 127, 143, 224, 227, 299, 305, 403}
    return 1383, np.array([i for i in range(1383) if i // 3 not in frozen])


def _owner(full_size, indices, *, active=True):
    optimizer = object.__new__(_Optimizer)
    optimizer._using_active_dofs = active
    optimizer._active_dof_indices = indices
    optimizer.geometry = SimpleNamespace(cart_coords=np.zeros(full_size))
    return optimizer


def _numpy(value):
    return value.detach().cpu().numpy() if isinstance(value, torch.Tensor) else value


def test_noncontiguous_compact_vector_is_unchanged(convert, basis):
    full_size, indices = basis
    optimizer = _owner(full_size, indices)
    vector = convert(np.arange(len(indices)) + 0.25)
    result = optimizer.active_from_full(vector)
    np.testing.assert_array_equal(_numpy(result), _numpy(vector))
    assert result is vector
    np.testing.assert_array_equal(optimizer.active_dof_indices, indices)


def test_repeated_projection_is_idempotent(convert, basis):
    full_size, indices = basis
    optimizer = _owner(full_size, indices)
    vector = convert(np.arange(full_size) + 0.25)
    once = optimizer.active_from_full(vector)
    twice = optimizer.active_from_full(once)
    np.testing.assert_array_equal(_numpy(twice), _numpy(once))


def test_full_vector_projects_and_lifts_in_original_order(convert, basis):
    full_size, indices = basis
    optimizer = _owner(full_size, indices)
    vector = convert(np.arange(full_size) + 0.25)
    projected = optimizer.active_from_full(vector)
    np.testing.assert_array_equal(_numpy(projected), _numpy(vector)[indices])
    lifted = optimizer.full_from_active(projected)
    expected = np.zeros(full_size)
    expected[indices] = _numpy(vector)[indices]
    np.testing.assert_array_equal(_numpy(lifted), expected)
    assert projected.dtype == vector.dtype
    if isinstance(vector, torch.Tensor):
        assert projected.device == vector.device


def test_inactive_projection_is_noop(convert):
    optimizer = _owner(9, np.array([0, 1, 2, 6, 7, 8]), active=False)
    vector = convert(np.arange(9))
    assert optimizer.active_from_full(vector) is vector


def test_missing_active_map_is_noop(convert):
    optimizer = _owner(9, None)
    vector = convert(np.arange(9))
    assert optimizer.active_from_full(vector) is vector


def test_negative_indices_are_removed_before_compact_check(convert):
    optimizer = _owner(3, np.array([-1, 0, 1]))
    vector = convert([0.25, 1.25, 2.25])
    np.testing.assert_array_equal(
        _numpy(optimizer.active_from_full(vector)), _numpy(vector)[:2]
    )
    np.testing.assert_array_equal(optimizer.active_dof_indices, [-1, 0, 1])


def test_internal_partial_hessian_keeps_internal_vector():
    optimizer = _owner(9, np.array([0, 1, 2, 6, 7, 8]))
    optimizer.H = np.diag([1.0, 2.0, 3.0, 4.0])
    optimizer.small_eigval_thresh = 1e-8
    optimizer.geometry = SimpleNamespace(
        internal=SimpleNamespace(project_hessian=lambda hessian: hessian),
        within_partial_hessian={"active_n_dof": 6}, coord_type="dlc",
        cart_coords=np.zeros(9), freeze_atoms=[1],
        hess_active_dof_indices=np.array([0, 1, 2, 6, 7, 8]),
        active_dof_indices=np.array([0, 1, 2, 6, 7, 8]), calculator=None,
    )
    vector = np.array([0.1, 0.2, 0.3, 0.4])
    projected, _, *_ = optimizer._hessian_system(vector)
    assert optimizer.using_active_dofs is False
    assert optimizer.active_from_full(vector) is vector
    np.testing.assert_array_equal(projected, vector)
