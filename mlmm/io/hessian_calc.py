"""Hessian utilities and vibrational analysis helpers for the ML/MM workflow.

These routines construct Hessians, compute frequencies and handle I/O.
"""

import os
from contextlib import nullcontext
import numpy as np
import time
from ase.constraints import FixAtoms
from ase.atoms import Atoms


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
    return_partial: bool = False,
):
    """
    Numerically compute the Cartesian Hessian (second derivatives).

    Finite differences are taken by displacing each *unfixed* atom by
    ``±delta`` Å along x, y, z.  Central differences give the sub-Hessian of
    the unfixed atoms, which is expanded to the full (3N, 3N) matrix with
    zeros for fixed atoms unless ``return_partial`` is set.

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
    return_partial : bool, default False
        Return the sub-Hessian of the unfixed atoms and their indices instead
        of the full matrix.

    Returns
    -------
    ndarray
        Cartesian Hessian of shape ``(3N, 3N)`` where
        ``N == len(atoms)`` (fixed-atom rows/cols are zero).
    (ndarray, ndarray)
        With ``return_partial``: the ``(3m, 3m)`` sub-Hessian of the ``m``
        unfixed atoms and their indices (int64, ascending, the row order).
    """
    try:
        delta = float(delta)
    except (TypeError, ValueError) as exc:
        raise ValueError("delta must be a finite positive number.") from exc
    if not np.isfinite(delta) or delta <= 0.0:
        raise ValueError("delta must be a finite positive number.")

    fixed = {
        i for c in atoms.constraints if isinstance(c, FixAtoms)
        for i in c.get_indices()
    }
    movable = np.asarray(
        [i for i in range(len(atoms)) if i not in fixed], dtype=np.int64
    )
    m = len(movable)
    n_dof = 3 * m

    # Short-circuit: all atoms frozen → return zeros
    if m == 0:
        if return_partial:
            return np.zeros((0, 0), dtype=dtype), movable
        N = len(atoms)
        return np.zeros((3 * N, 3 * N), dtype=dtype)

    H_sub = np.empty((n_dof, n_dof), dtype=dtype)
    row = 0

    log_cm = nullcontext(None)
    if info_path is not None:
        info_dir = os.path.dirname(os.fspath(info_path))
        if info_dir:
            os.makedirs(info_dir, exist_ok=True)
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

    if return_partial:
        return H_sub, movable
    N = len(atoms)
    if m == N:
        return H_sub
    idx3 = (3 * movable[:, None] + np.arange(3)).ravel()
    H = np.zeros((3 * N, 3 * N), dtype=dtype)
    H[np.ix_(idx3, idx3)] = H_sub
    return H

