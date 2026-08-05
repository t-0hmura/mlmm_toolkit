"""Regression tests for inline XYZ atom-row parsing."""

import numpy as np
import pytest

from pysisyphus.xyzloader import split_xyz_str


def test_inline_xyz_accepts_signed_scientific_notation() -> None:
    atoms_coords = split_xyz_str(
        "2\nscientific coordinates\nH +1e-3 -2E+2 3.0e0\nO -4E-1 +5e+0 -6E0\n"
    )

    atoms, coords = atoms_coords[0]
    assert atoms == ("H", "O")
    np.testing.assert_allclose(
        coords,
        [[1.0e-3, -2.0e2, 3.0], [-4.0e-1, 5.0, -6.0]],
        rtol=0.0,
        atol=0.0,
    )


def test_inline_xyz_rejects_extra_atom_fields() -> None:
    with pytest.raises(AssertionError):
        split_xyz_str("1\ncomment\nH 0 0 0 extra\n")
