"""Regression tests for internal coordinates with a partial Cartesian Hessian."""

from types import SimpleNamespace

import numpy as np

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer


class _Internal:
    def project_hessian(self, hessian):
        return hessian


class _Optimizer(HessianOptimizer):
    def optimize(self):
        raise NotImplementedError


def test_internal_partial_hessian_does_not_apply_cartesian_active_indices() -> None:
    optimizer = object.__new__(_Optimizer)
    optimizer.H = np.diag([1.0, 2.0, 3.0, 4.0])
    optimizer.small_eigval_thresh = 1.0e-8
    optimizer.geometry = SimpleNamespace(
        internal=_Internal(),
        within_partial_hessian={"active_n_dof": 6},
        coord_type="dlc",
        cart_coords=np.zeros(9),
        freeze_atoms=[2],
        hess_active_dof_indices=np.arange(6),
        active_dof_indices=np.arange(6),
        calculator=None,
    )

    gradient = np.asarray([0.1, 0.2, 0.3, 0.4])
    projected_gradient, hessian, *_ = optimizer._hessian_system(gradient)

    assert optimizer.using_active_dofs is False
    assert np.array_equal(projected_gradient, gradient)
    assert hessian.shape == (4, 4)
    step = np.arange(4, dtype=float)
    assert np.array_equal(optimizer.full_from_active(step), step)
