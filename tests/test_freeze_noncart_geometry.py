import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.intcoords.exceptions import (
    NeedNewInternalsException,
    RebuiltInternalsException,
)


ATOMS = ["C", "C", "O", "C", "N", "H"]
COORDS = np.array(
    [
        0.0,
        0.0,
        0.0,
        1.4,
        0.0,
        0.0,
        2.1,
        1.2,
        0.0,
        3.5,
        1.2,
        0.1,
        4.2,
        0.1,
        0.2,
        5.6,
        0.1,
        0.0,
    ],
    dtype=float,
)
FREEZE = np.array([0, 2], dtype=int)


def _assert_no_frozen_primitives(geom):
    freeze = set(map(int, geom.freeze_atoms))
    for typed_prim in geom.internal.typed_prims:
        assert set(typed_prim[1:]).isdisjoint(freeze), typed_prim


@pytest.mark.parametrize("coord_type", ["redund", "dlc", "tric"])
def test_noncart_freeze_atoms_are_excluded_from_primitives_by_default(coord_type):
    geom = Geometry(ATOMS, COORDS, coord_type=coord_type, freeze_atoms=FREEZE)

    assert geom.coord_kwargs["freeze_atoms_exclude"] is True
    _assert_no_frozen_primitives(geom)


@pytest.mark.parametrize("coord_type", ["redund", "dlc", "tric"])
def test_setting_noncart_freeze_atoms_rebuilds_primitives(coord_type):
    geom = Geometry(ATOMS, COORDS, coord_type=coord_type)

    geom.freeze_atoms = FREEZE

    assert geom.coord_kwargs["freeze_atoms_exclude"] is True
    np.testing.assert_array_equal(geom.internal.freeze_atoms, FREEZE)
    _assert_no_frozen_primitives(geom)


@pytest.mark.parametrize("coord_type", ["redund", "dlc", "tric"])
def test_noncart_internal_step_preserves_frozen_cartesians(coord_type):
    geom = Geometry(ATOMS, COORDS, coord_type=coord_type, freeze_atoms=FREEZE)
    frozen_ref = geom.coords3d[FREEZE].copy()

    geom.coords = geom.coords.copy() + 1.0e-3

    np.testing.assert_allclose(geom.coords3d[FREEZE], frozen_ref, atol=1.0e-12)
    np.testing.assert_array_equal(geom.internal.freeze_atoms, FREEZE)


@pytest.mark.parametrize("coord_type", ["redund", "dlc", "tric"])
def test_reset_coords_preserves_internal_freeze_atoms(coord_type):
    geom = Geometry(ATOMS, COORDS, coord_type=coord_type, freeze_atoms=FREEZE)

    geom.reset_coords()

    np.testing.assert_array_equal(geom.internal.freeze_atoms, FREEZE)


@pytest.mark.parametrize("coord_type", ["redund", "dlc", "tric"])
def test_rebuilt_internals_reimposes_frozen_cartesians(coord_type):
    geom = Geometry(ATOMS, COORDS, coord_type=coord_type, freeze_atoms=FREEZE)
    frozen_ref = geom.coords3d[FREEZE].copy()

    def raise_rebuild(_int_step, update_constraints=False):
        rebuilt = geom.coords3d.copy()
        rebuilt[FREEZE] += 5.0
        rebuilt[1] += 0.1
        raise NeedNewInternalsException(rebuilt, invalid_inds=())

    geom.internal.transform_int_step = raise_rebuild

    with pytest.raises(RebuiltInternalsException):
        geom.coords = geom.coords.copy() + 1.0e-3

    np.testing.assert_allclose(geom.coords3d[FREEZE], frozen_ref, atol=1.0e-12)
    np.testing.assert_array_equal(geom.internal.freeze_atoms, FREEZE)
