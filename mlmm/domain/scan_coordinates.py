"""Internal-coordinate primitives shared by scan workflows."""

from __future__ import annotations

import math
from numbers import Real
from typing import Any, Iterable, Sequence, Tuple

import click
import numpy as np
from ase.geometry.geometry import (
    get_angles,
    get_angles_derivatives,
    get_dihedrals,
    get_dihedrals_derivatives,
)

from pysisyphus.constants import BOHR2ANG


Target = Tuple[Any, ...]
Range = Tuple[Any, ...]


def coordinate_kind(entry: Sequence[Any], *, is_range: bool = False) -> str:
    arity = len(entry) - (2 if is_range else 1)
    try:
        return {2: "distance", 3: "angle", 4: "dihedral"}[arity]
    except KeyError as exc:
        form = "(atoms..., low, high)" if is_range else "(atoms..., target)"
        raise click.BadParameter(f"Invalid scan coordinate {tuple(entry)!r}; expected {form}.") from exc


def coordinate_atoms(entry: Sequence[Any], *, is_range: bool = False) -> tuple[int, ...]:
    stop = -2 if is_range else -1
    return tuple(int(value) for value in entry[:stop])


def coordinate_target(entry: Sequence[Any]) -> float:
    return float(entry[-1])


def coordinate_bounds(entry: Sequence[Any]) -> tuple[float, float]:
    return float(entry[-2]), float(entry[-1])


def coordinate_unit(kind: str) -> str:
    return "angstrom" if kind == "distance" else "degree"


def coordinate_symbol(kind: str) -> str:
    return {"distance": "r", "angle": "angle", "dihedral": "dihedral"}[kind]


def canonical_axis_key(kind: str, atoms: Sequence[int]) -> tuple[str, tuple[int, ...]]:
    atom_tuple = tuple(int(i) for i in atoms)
    reverse = tuple(reversed(atom_tuple))
    return kind, min(atom_tuple, reverse)


def validate_coordinate_value(kind: str, value: float, *, context: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise click.BadParameter(f"{context} must be finite.")
    if kind == "distance" and value <= 0.0:
        raise click.BadParameter(f"{context} distance must be > 0 Å.")
    if kind == "angle" and not (0.0 < value < 180.0):
        raise click.BadParameter(f"{context} angle must be between 0 and 180 degrees.")
    return value


def parse_coordinate_entry(
    entry: Any,
    *,
    is_range: bool,
    one_based: bool,
    atom_meta: Sequence[dict[str, Any]] | None,
    context: str,
    resolve_index,
) -> tuple[Any, ...]:
    if not isinstance(entry, (list, tuple)):
        raise click.BadParameter(f"{context} must be a tuple/list.")
    kind = coordinate_kind(entry, is_range=is_range)
    n_atoms = {"distance": 2, "angle": 3, "dihedral": 4}[kind]
    n_values = 2 if is_range else 1
    if len(entry) != n_atoms + n_values or not all(
        isinstance(value, Real) for value in entry[-n_values:]
    ):
        forms = (
            "(i,j,low,high), (i,j,k,low,high), or (i,j,k,l,low,high)"
            if is_range
            else "(i,j,target), (i,j,k,target), or (i,j,k,l,target)"
        )
        raise click.BadParameter(f"{context} must be {forms}: got {entry!r}.")
    atoms = tuple(
        resolve_index(
            value,
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{context} atom {position + 1}",
        )
        for position, value in enumerate(entry[:n_atoms])
    )
    if len(set(atoms)) != len(atoms):
        raise click.BadParameter(f"{context} repeats an atom index.")
    values = tuple(
        validate_coordinate_value(kind, float(value), context=context)
        for value in entry[-n_values:]
    )
    return (*atoms, *values)


def signed_periodic_delta_deg(value: float, target: float) -> float:
    return (float(value) - float(target) + 180.0) % 360.0 - 180.0


def coordinate_delta(kind: str, value: float, target: float) -> float:
    if kind == "dihedral":
        return signed_periodic_delta_deg(value, target)
    return float(value) - float(target)


def coordinate_value(coords_bohr: np.ndarray, entry: Sequence[Any], *, is_range: bool = False) -> float:
    coords = np.asarray(coords_bohr, dtype=float).reshape(-1, 3)
    kind = coordinate_kind(entry, is_range=is_range)
    atoms = coordinate_atoms(entry, is_range=is_range)
    if kind == "distance":
        return float(np.linalg.norm(coords[atoms[0]] - coords[atoms[1]]) * BOHR2ANG)
    if kind == "angle":
        i, j, k = atoms
        return float(get_angles((coords[i] - coords[j])[None, :], (coords[k] - coords[j])[None, :])[0])
    i, j, k, l = atoms
    value = float(
        get_dihedrals(
            (coords[j] - coords[i])[None, :],
            (coords[k] - coords[j])[None, :],
            (coords[l] - coords[k])[None, :],
        )[0]
    )
    return signed_periodic_delta_deg(value, 0.0)


def coordinate_derivative(coords_bohr: np.ndarray, entry: Sequence[Any]) -> np.ndarray:
    """Return dq/dx for q in Å (distance) or radians (angles)."""
    coords = np.asarray(coords_bohr, dtype=float).reshape(-1, 3)
    deriv = np.zeros_like(coords)
    kind = coordinate_kind(entry)
    atoms = coordinate_atoms(entry)
    if kind == "distance":
        i, j = atoms
        delta = coords[i] - coords[j]
        norm = float(np.linalg.norm(delta))
        if norm < 1.0e-14:
            raise ValueError("Distance derivative is undefined for coincident atoms.")
        unit = delta / norm * BOHR2ANG
        deriv[i] = unit
        deriv[j] = -unit
        return deriv.reshape(-1)
    if kind == "angle":
        i, j, k = atoms
        local = get_angles_derivatives(
            (coords[i] - coords[j])[None, :],
            (coords[k] - coords[j])[None, :],
        )[0] * (math.pi / 180.0)
    else:
        i, j, k, l = atoms
        local = get_dihedrals_derivatives(
            (coords[j] - coords[i])[None, :],
            (coords[k] - coords[j])[None, :],
            (coords[l] - coords[k])[None, :],
        )[0] * (math.pi / 180.0)
    for atom, vector in zip(atoms, local):
        deriv[atom] += vector
    return deriv.reshape(-1)


def format_coordinate(entry: Sequence[Any], *, one_based: bool = True, is_range: bool = False) -> tuple[Any, ...]:
    atoms = coordinate_atoms(entry, is_range=is_range)
    if one_based:
        atoms = tuple(i + 1 for i in atoms)
    values = tuple(float(v) for v in entry[(-2 if is_range else -1):])
    return (*atoms, *values)


def coordinate_step_cap(kind: str, distance_step: float, angle_step: float, dihedral_step: float) -> float:
    return {
        "distance": float(distance_step),
        "angle": float(angle_step),
        "dihedral": float(dihedral_step),
    }[kind]
