"""Hessian utilities and vibrational analysis helpers for the ML/MM workflow.

These routines construct Hessians, compute frequencies and handle I/O.
"""

import os
from contextlib import nullcontext
import numpy as np
import time
import torch
from typing import Sequence
from ase.data import atomic_masses
from ase.constraints import FixAtoms
from ase.atoms import Atoms
from pysisyphus.constants import AMU2AU


def clone_atoms(atoms):
    """Return a shallow copy of *atoms* keeping constraints (no PBC)."""
    new = Atoms(
        symbols=atoms.get_chemical_symbols(),
        positions=atoms.get_positions(),
        pbc=False,
    )

    if atoms.constraints:
        fix_idx = [
            i for c in atoms.constraints if isinstance(c, FixAtoms)
            for i in c.get_indices()
        ]
        if fix_idx:
            new.set_constraint(FixAtoms(indices=fix_idx))

    return new


def hessian_calc(
    atoms,
    calc,
    delta: float = 0.01,
    info_path: str | None = None,
    *,
    dtype: np.dtype = np.float64,
):
    """
    Numerically compute the full Cartesian Hessian (second derivatives).

    Finite differences are taken by displacing each *unfixed* atom by
    ``±delta`` Å along x, y, z.  Central differences give the sub-Hessian,
    which is **always expanded to the full (3N, 3N) matrix** by padding
    zeros for fixed atoms before returning.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure to differentiate.
    calc : ase.Calculator
        Calculator providing forces.
    delta : float, default 0.01 Å
        Displacement size.
    info_path : str | None
        If given, progress info is appended to this file.
    dtype : numpy dtype, default float64
        Data type of the returned Hessian.

    Returns
    -------
    ndarray
        Cartesian Hessian of shape ``(3N, 3N)`` where
        ``N == len(atoms)`` (fixed-atom rows/cols are zero).
    """
    fixed = {
        i for c in atoms.constraints if isinstance(c, FixAtoms)
        for i in c.get_indices()
    }
    movable = np.asarray([i for i in range(len(atoms)) if i not in fixed])
    m = len(movable)
    n_dof = 3 * m

    # Short-circuit: all atoms frozen → return zeros
    if m == 0:
        N = len(atoms)
        return np.zeros((3 * N, 3 * N))

    H_sub = np.empty((n_dof, n_dof), dtype=dtype)
    row = 0

    log_cm = nullcontext(None)
    if info_path is not None:
        os.makedirs(os.path.dirname(info_path), exist_ok=True)
        log_cm = open(info_path, "w", encoding="utf-8")

    with log_cm as log:
        if info_path is not None:
            log.write("Hessian calculation by numerical differentiation\n")
            log.write("------------------------------------------------\n")
            t0 = time.time()
            checkpoints = [int(m * k / 10) for k in range(1, 11)]  # every 10 %

        for count, a in enumerate(movable, start=1):
            for coord in range(3):
                plus = clone_atoms(atoms)
                plus.calc = calc
                plus.positions[a, coord] += delta
                Fp = plus.calc.get_forces(plus)

                minus = clone_atoms(atoms)
                minus.calc = calc
                minus.positions[a, coord] -= delta
                Fm = minus.calc.get_forces(minus)

                # Central difference:  −∂F/∂x = ∂²E/∂x²
                H_sub[row] = (Fm - Fp)[movable].ravel() / (2 * delta)
                row += 1

            if info_path is not None and count in checkpoints:
                pct = checkpoints.index(count) * 10 + 10
                elapsed = time.time() - t0
                total_est = elapsed * 100.0 / pct
                remaining = total_est - elapsed
                me, se = divmod(int(elapsed), 60)
                mr, sr = divmod(int(remaining), 60)
                log.write(
                    f"{pct:3d}% ({count}/{m})   Elapsed {me}m{se}s   ETA {mr}m{sr}s\n"
                )
                log.flush()

        if info_path is not None:
            elapsed = time.time() - t0
            me, se = divmod(int(elapsed), 60)
            log.write(f"Done in {me}m{se}s\n")

    N = len(atoms)
    H = np.zeros((3 * N, 3 * N), dtype=dtype)
    for i_local, i_atom in enumerate(movable):
        for j_local, j_atom in enumerate(movable):
            H[
                3 * i_atom : 3 * i_atom + 3,
                3 * j_atom : 3 * j_atom + 3,
            ] = H_sub[
                3 * i_local : 3 * i_local + 3,
                3 * j_local : 3 * j_local + 3,
            ]
    return H


def _get_masses(Z):
    """Return atomic masses (amu) for a list/array of atomic numbers."""
    return np.array([atomic_masses[z] for z in Z])


def _build_tr_basis(coords_bohr: torch.Tensor, masses_amu: torch.Tensor) -> torch.Tensor:
    """Return (3N,6) mass-weighted basis for Tx,Ty,Tz,Rx,Ry,Rz."""
    device, dtype = coords_bohr.device, coords_bohr.dtype
    N = coords_bohr.shape[0]

    m_au = masses_amu.to(dtype=dtype, device=device) * AMU2AU
    m_sqrt = torch.sqrt(m_au).reshape(-1, 1)

    com = (m_au.reshape(-1, 1) * coords_bohr).sum(0) / m_au.sum()
    x = coords_bohr - com

    eye3 = torch.eye(3, dtype=dtype, device=device)
    cols = []

    for i in range(3):                                   # Tx, Ty, Tz
        cols.append((eye3[i].repeat(N, 1) * m_sqrt).reshape(-1, 1))
    for i in range(3):                                   # Rx, Ry, Rz
        rot = torch.cross(x, eye3[i].expand_as(x), dim=1) * m_sqrt
        cols.append(rot.reshape(-1, 1))

    return torch.cat(cols, dim=1)


