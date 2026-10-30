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


class _CoupledQuadraticCalculator:
    """Forces of E = x.K.x/2 with a fixed dense symmetric K."""

    def __init__(self, n_atoms: int) -> None:
        a = np.random.default_rng(7).standard_normal((3 * n_atoms, 3 * n_atoms))
        self.k = a + a.T

    def get_forces(self, atoms):
        return -(self.k @ atoms.get_positions().reshape(-1)).reshape(-1, 3)


def _block_loop_expansion(h_sub, movable, n_atoms, dtype):
    """The former (3N, 3N) expansion of hessian_calc, one 3x3 block at a time."""
    h = np.zeros((3 * n_atoms, 3 * n_atoms), dtype=dtype)
    for i_local, i_atom in enumerate(movable):
        for j_local, j_atom in enumerate(movable):
            h[3 * i_atom:3 * i_atom + 3, 3 * j_atom:3 * j_atom + 3] = h_sub[
                3 * i_local:3 * i_local + 3, 3 * j_local:3 * j_local + 3
            ]
    return h


@pytest.mark.parametrize("dtype", [np.float64, np.float32])
@pytest.mark.parametrize("fixed", [[], [5, 0, 3], [1, 2, 4, 6]])
def test_hessian_calc_partial_equals_restricted_full(dtype, fixed) -> None:
    n_atoms = 7
    positions = np.random.default_rng(3).standard_normal((n_atoms, 3))
    atoms = Atoms("H7", positions=positions)
    if fixed:
        atoms.set_constraint(FixAtoms(indices=fixed))
    calc = _CoupledQuadraticCalculator(n_atoms)

    full = hessian_calc(atoms, calc, dtype=dtype)
    partial, movable = hessian_calc(atoms, calc, dtype=dtype, return_partial=True)

    np.testing.assert_array_equal(
        movable, [i for i in range(n_atoms) if i not in fixed]
    )
    assert movable.dtype == np.int64
    assert full.dtype == partial.dtype == dtype
    assert full.shape == (3 * n_atoms, 3 * n_atoms)
    idx3 = (3 * movable[:, None] + np.arange(3)).ravel()
    assert full[np.ix_(idx3, idx3)].tobytes() == partial.tobytes()
    expanded = _block_loop_expansion(partial, movable, n_atoms, dtype)
    assert full.tobytes() == expanded.tobytes()


def test_hessian_calc_all_fixed_partial_is_empty() -> None:
    atoms = Atoms("H2", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    atoms.set_constraint(FixAtoms(indices=[0, 1]))

    partial, movable = hessian_calc(
        atoms, object(), dtype=np.float32, return_partial=True
    )

    assert partial.shape == (0, 0)
    assert partial.dtype == np.float32
    assert movable.shape == (0,)
    assert movable.dtype == np.int64


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
