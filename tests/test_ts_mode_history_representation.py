"""Small TS-mode histories follow the current eigensystem representation."""
import numpy as np
import pytest
import torch

from pysisyphus._array import as_numpy
from pysisyphus.Geometry import Geometry
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer


BACKENDS = [
    ("numpy", "numpy"), ("torch", "torch"), ("cuda", "cuda"),
    ("numpy", "torch"), ("numpy", "cuda"),
    ("torch", "numpy"), ("cuda", "numpy"),
    ("torch", "cuda"), ("cuda", "torch"),
]


def represented(values, backend, *, history=False):
    dtype = np.float32 if history else np.float64
    values = np.asarray(values, dtype=dtype)
    if backend == "numpy":
        return values.copy()
    return torch.as_tensor(values, dtype=torch.float32 if history else torch.float64,
                           device="cuda" if backend == "cuda" else "cpu")


def optimizer_for_modes(tmp_path):
    geometry = Geometry(["H", "H"], [-.5, 0., 0., .5, 0., 0.], coord_type="cart")
    optimizer = RSPRFOptimizer(
        geometry, roots=[0], hessian_init="unit", track_mode_by_overlap=True,
        verify_saddle=False, dump=False, out_dir=tmp_path,
    )
    return optimizer


@pytest.mark.parametrize("history_backend,eigen_backend", BACKENDS)
@pytest.mark.parametrize("basis", ["full", "full-to-active", "ordered-active"])
@pytest.mark.parametrize("history_owner", ["overlap", "physical"])
def test_mode_history_representation_preserves_existing_basis_and_root(
    tmp_path, history_backend, eigen_backend, basis, history_owner,
):
    if "cuda" in (history_backend, eigen_backend) and not torch.cuda.is_available():
        pytest.skip("CUDA is unavailable")
    optimizer = optimizer_for_modes(tmp_path)
    if basis == "full":
        vectors = np.eye(6)
        mode = np.eye(6)[2]
        expected_root = 2
    elif basis == "full-to-active":
        # Existing supported projection: a full Cartesian history onto atom1.
        optimizer._using_active_dofs = True
        optimizer._active_dof_indices = np.array([3, 4, 5])
        vectors = np.eye(3)
        mode = np.eye(6)[4]
        expected_root = 1
    else:
        # Already compact, in the declared same-size noncanonical atom order.
        # Representation conversion must not apply that permutation again.
        order = np.array([3, 4, 5, 0, 1, 2])
        optimizer._using_active_dofs = True
        optimizer._active_dof_indices = order.copy()
        vectors = np.eye(6)[order]
        mode = np.eye(6)[0][order]
        expected_root = 0
    values = -np.arange(vectors.shape[1], 0, -1, dtype=float)
    old_mode = represented(mode, history_backend, history=True)
    optimizer.ts_modes = old_mode.reshape(1, -1)
    optimizer.ts_mode_eigvals = represented([-1.], history_backend, history=True)
    if history_owner == "physical":
        optimizer._physical_ts_mode = old_mode
    initial_coords = optimizer.geometry.cart_coords.copy()
    eigvals = represented(values, eigen_backend)
    eigvecs = represented(vectors, eigen_backend)

    optimizer.update_ts_mode(eigvals, eigvecs)

    assert int(optimizer.roots[0]) == expected_root
    np.testing.assert_array_equal(as_numpy(optimizer.ts_modes), vectors[:, [expected_root]].T)
    np.testing.assert_array_equal(as_numpy(optimizer.ts_mode_eigvals), values[[expected_root]])
    if eigen_backend == "numpy":
        assert isinstance(optimizer.ts_modes, np.ndarray)
        assert optimizer.ts_modes.dtype == np.float64
    else:
        assert isinstance(optimizer.ts_modes, torch.Tensor)
        assert optimizer.ts_modes.dtype == eigvecs.dtype
        assert optimizer.ts_modes.device == eigvecs.device
    np.testing.assert_array_equal(optimizer.geometry.cart_coords, initial_coords)
    assert optimizer.geometry.calculator is None


@pytest.mark.parametrize("backend", ["numpy", "torch", "cuda"])
@pytest.mark.parametrize("history_owner", ["overlap", "physical"])
def test_representation_conversion_does_not_hide_incompatible_mode_width(
    tmp_path, backend, history_owner,
):
    if backend == "cuda" and not torch.cuda.is_available():
        pytest.skip("CUDA is unavailable")
    optimizer = optimizer_for_modes(tmp_path)
    mode = represented(np.ones(5), backend, history=True)
    optimizer.ts_modes = mode.reshape(1, -1)
    optimizer.ts_mode_eigvals = represented([-1.], backend, history=True)
    if history_owner == "physical":
        optimizer._physical_ts_mode = mode
    # No active-coordinate projection applies here. Preserve the existing
    # contraction failure; do not pad, truncate, or manufacture a new mode.
    with pytest.raises((ValueError, RuntimeError)):
        optimizer.update_ts_mode(represented([-6., -5., -4., -3., -2., -1.], backend),
                                 represented(np.eye(6), backend))
