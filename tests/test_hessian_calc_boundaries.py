"""Boundary tests for numerical Cartesian Hessians."""

from __future__ import annotations

import numpy as np
import pytest
from ase import Atoms
from ase.constraints import FixAtoms

from mlmm.io.hessian_calc import hessian_calc


class _HarmonicCalculator:
    def get_forces(self, atoms):
        return -np.asarray(atoms.get_positions(), dtype=float)


@pytest.mark.parametrize("delta", [0.0, -0.01, np.nan, np.inf])
def test_hessian_calc_rejects_invalid_delta_before_force_evaluation(delta) -> None:
    atoms = Atoms("H", positions=[[0.0, 0.0, 0.0]])

    with pytest.raises(ValueError, match="finite positive"):
        hessian_calc(atoms, object(), delta=delta)


def test_hessian_calc_all_fixed_preserves_requested_dtype() -> None:
    atoms = Atoms("H", positions=[[0.0, 0.0, 0.0]])
    atoms.set_constraint(FixAtoms(indices=[0]))

    hessian = hessian_calc(atoms, object(), dtype=np.float32)

    assert hessian.shape == (3, 3)
    assert hessian.dtype == np.float32


def test_hessian_calc_accepts_log_filename_without_directory(
    tmp_path, monkeypatch
) -> None:
    atoms = Atoms("H", positions=[[0.1, -0.2, 0.3]])
    monkeypatch.chdir(tmp_path)

    hessian = hessian_calc(
        atoms,
        _HarmonicCalculator(),
        info_path="progress.log",
    )

    np.testing.assert_allclose(hessian, np.eye(3))
    assert (tmp_path / "progress.log").is_file()
