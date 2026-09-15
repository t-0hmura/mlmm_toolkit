"""Only unconstrained, pure whole-system translations are numerical artifacts."""
from types import SimpleNamespace

import numpy as np
import pytest
import torch

from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer


@pytest.mark.parametrize("backend", ["numpy", "torch"])
@pytest.mark.parametrize("permuted", [False, True])
def test_translation_filter_preserves_rotations_and_mixed_modes(backend, permuted):
    opt = object.__new__(RSPRFOptimizer)
    opt.geometry = SimpleNamespace(
        coord_type="cart", freeze_atoms=[],
        cart_coords=np.array([-.7, 0., 0., .7, 0., 0.]),
    )
    permutation = np.array([4, 1, 5, 2, 3, 0]) if permuted else np.arange(6)
    opt._using_active_dofs = permuted
    opt._active_dof_indices = permutation if permuted else None
    translation = np.array([1., 0., 0., 1., 0., 0.]) / np.sqrt(2.)
    rotation = np.array([0., -1., 0., 0., 1., 0.]) / np.sqrt(2.)
    mixed = translation + .01 * rotation
    vectors = np.column_stack([translation, rotation, mixed, np.zeros(6)])
    vectors = vectors[permutation]
    if backend == "torch":
        vectors = torch.as_tensor(vectors)
    before = vectors.clone() if backend == "torch" else vectors.copy()
    mask = opt._translation_mode_mask(vectors)
    if backend == "torch":
        assert mask.dtype == torch.bool and mask.device == vectors.device
        torch.testing.assert_close(vectors, before, rtol=0., atol=0.)
        mask = mask.numpy()
    else:
        assert mask.dtype == np.bool_
        np.testing.assert_array_equal(vectors, before)
    np.testing.assert_array_equal(mask, [True, False, False, False])


@pytest.mark.parametrize("coord_type,frozen,rows", [
    ("cart", [0], 6),
    ("cart", [0], 3),
    ("cart", [], 3),
    ("dlc", [], 6),
    ("mwcartesian", [], 6),
])
def test_translation_filter_does_not_reinterpret_other_spaces(coord_type, frozen, rows):
    opt = object.__new__(RSPRFOptimizer)
    opt.geometry = SimpleNamespace(
        coord_type=coord_type, freeze_atoms=frozen, cart_coords=np.zeros(6),
    )
    opt._using_active_dofs = False
    opt._active_dof_indices = None
    np.testing.assert_array_equal(
        opt._translation_mode_mask(np.ones((rows, 1))), [False],
    )


def test_translation_filter_requires_a_complete_active_permutation():
    opt = object.__new__(RSPRFOptimizer)
    opt.geometry = SimpleNamespace(
        coord_type="cart", freeze_atoms=[], cart_coords=np.zeros(6),
    )
    opt._using_active_dofs = True
    opt._active_dof_indices = np.array([0, 1, 2, 3, 4, 4])
    np.testing.assert_array_equal(
        opt._translation_mode_mask(np.ones((6, 1))), [False],
    )
