"""Plain NumPy ``.npy`` Hessian files for ``--dump-hess`` / ``--read-hess``.

One square float64 array: the Cartesian Hessian in Hartree/bohr^2 (not
mass-weighted), atoms in input order, for all atoms (3N x 3N), the movable
atoms only, or the atoms in the Hessian calculation.
"""

from __future__ import annotations

import zipfile
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np


def save_hessian_file(path: Path | str, hessian: Any) -> Path:
    """Atomically write ``hessian`` as a ``.npy`` array at the exact path and return it."""
    hess = np.asarray(hessian, dtype=np.float64)
    if hess.ndim != 2 or hess.shape[0] != hess.shape[1] or not np.isfinite(hess).all():
        raise ValueError(f"Hessian must be a finite square matrix, got shape {hess.shape}.")
    from mlmm.core.result_commit import atomic_write_exact

    return atomic_write_exact(Path(path), lambda stream: np.save(stream, hess))


def load_hessian_file(
    path: Path | str, *, n_atoms: int, active_dofs: Optional[Sequence[int]] = None
) -> np.ndarray:
    """Read a square, finite, symmetric ``.npy`` Hessian in the ``active_dofs`` basis.

    The file holds either all 3N DOFs (restricted here to ``active_dofs``) or
    exactly ``active_dofs``, in that order.
    """
    full = 3 * int(n_atoms)
    dofs = (
        np.arange(full, dtype=np.int64)
        if active_dofs is None
        else np.asarray(active_dofs, dtype=np.int64).reshape(-1)
    )
    hess = _read_square_hessian(path)
    n = hess.shape[0]
    if n == full and dofs.size != full:
        hess = hess[np.ix_(dofs, dofs)]
    elif n != dofs.size:
        movable = "" if dofs.size == full else f" or {dofs.size}x{dofs.size} (movable atoms only)"
        raise ValueError(
            f"Hessian in {path} is {n}x{n}; expected {full}x{full} (all atoms){movable}."
        )
    return hess


def load_hessian_file_in_bases(
    path: Path | str, *, n_atoms: int, bases: Sequence[tuple[str, Sequence[int]]]
) -> tuple[np.ndarray, list[int]]:
    """Read a ``.npy`` Hessian held in all 3N DOFs or exactly one of ``bases``.

    ``bases`` are (label, ordered DOFs) pairs of different sizes; a 3N file is
    restricted to the first one.  Returns the Hessian and the DOFs it covers.
    """
    hess = _read_square_hessian(path)
    full = 3 * int(n_atoms)
    n = hess.shape[0]
    primary = [int(d) for d in bases[0][1]]
    if n == full:
        return (hess if len(primary) == full else hess[np.ix_(primary, primary)]), primary
    for _label, dofs in bases:
        if n == len(dofs):
            return hess, [int(d) for d in dofs]
    sizes = "".join(f" or {len(d)}x{len(d)} ({label})" for label, d in bases if len(d) != full)
    raise ValueError(
        f"Hessian in {path} is {n}x{n}; expected {full}x{full} (all atoms){sizes}; "
        "write it with the same --hess-cutoff and layer settings."
    )


def _read_square_hessian(path: Path | str) -> np.ndarray:
    """Load a real, square, finite, symmetric ``.npy`` array as float64."""
    try:
        hess = np.load(path, allow_pickle=False)
    except (OSError, ValueError, EOFError, zipfile.BadZipFile) as exc:
        raise ValueError(f"Cannot read {path} as a NumPy .npy array: {exc}") from exc
    if not isinstance(hess, np.ndarray):
        hess.close()
        raise ValueError(f"{path} is not a .npy array; write the Hessian with numpy.save.")
    if hess.ndim != 2 or hess.shape[0] != hess.shape[1] or hess.dtype.kind not in "iuf":
        raise ValueError(
            f"Hessian in {path} must be a real square matrix, got {hess.dtype} {hess.shape}."
        )
    hess = np.asarray(hess, dtype=np.float64)
    if not np.isfinite(hess).all():
        raise ValueError(f"Hessian in {path} contains non-finite values.")
    if not np.allclose(hess, hess.T):
        raise ValueError(f"Hessian in {path} is not symmetric.")
    return hess
