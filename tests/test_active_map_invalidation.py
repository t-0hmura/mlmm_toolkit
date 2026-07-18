"""M51: a changed freeze mask must invalidate the memoized active atom/DOF maps.

The active-atom / active-DOF views are memoized on first access.  Before this
fix, assigning a new ``freeze_atoms`` mask left both memoized arrays stale, so a
downstream partial-Hessian / cache-identity consumer would reuse an active basis
derived from the *previous* mask.
"""

from __future__ import annotations

import numpy as np

from pysisyphus.Geometry import Geometry


def _geom() -> Geometry:
    return Geometry(
        ["H", "H", "O"],
        np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]),
        coord_type="cart",
    )


def test_freeze_mask_change_invalidates_primed_active_maps() -> None:
    geom = _geom()

    # Prime the memoized views before mutating the mask.
    assert list(geom.active_atom_indices) == [0, 1, 2]
    assert list(geom.active_dof_indices) == list(range(9))

    geom.freeze_atoms = [0]

    assert list(geom.active_atom_indices) == [1, 2]
    assert list(geom.active_dof_indices) == [3, 4, 5, 6, 7, 8]


def test_unfreeze_reinvalidates_active_maps() -> None:
    geom = _geom()
    geom.freeze_atoms = [0]
    assert list(geom.active_atom_indices) == [1, 2]

    geom.freeze_atoms = []
    assert list(geom.active_atom_indices) == [0, 1, 2]
    assert list(geom.active_dof_indices) == list(range(9))


def test_equal_normalized_mask_keeps_active_maps_equal() -> None:
    geom = _geom()
    geom.freeze_atoms = [1]
    atoms_before = list(geom.active_atom_indices)
    dofs_before = list(geom.active_dof_indices)

    # Re-assigning an equal mask must not change the resolved maps.
    geom.freeze_atoms = [1]
    assert list(geom.active_atom_indices) == atoms_before == [0, 2]
    assert list(geom.active_dof_indices) == dofs_before == [0, 1, 2, 6, 7, 8]


def test_refreshed_active_map_rejects_pre_mutation_cache_entry() -> None:
    """A refreshed active map feeds the Hessian-cache identity (M70 hook)."""

    from mlmm.io import hessian_cache

    geom = _geom()
    calc_cfg = {
        "backend": "uma",
        "model": "m",
        "precision": "fp64",
        "charge": 0,
        "spin": 1,
    }
    ident_before = hessian_cache.identity_from_context(
        geom, {**calc_cfg, "freeze_atoms": []}, role="ts"
    )
    geom.freeze_atoms = [0]
    ident_after = hessian_cache.identity_from_context(
        geom, {**calc_cfg, "freeze_atoms": [0]}, role="ts"
    )

    assert (
        ident_before["active_space"]["active_dofs"]
        != ident_after["active_space"]["active_dofs"]
    )
