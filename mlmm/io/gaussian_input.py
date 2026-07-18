"""Parser for ordinary Gaussian molecular input files.

This module intentionally does not parse Gaussian ONIOM inputs.  The latter
have a six-integer charge/layer contract and remain owned by
``mlmm.workflows.oniom_import``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import re

from ase.data import chemical_symbols


_TWO_INT_RE = re.compile(r"^\s*([+-]?\d+)\s+(\d+)\s*$")
_ELEMENT_RE = re.compile(r"^([A-Za-z]{1,2}|\d+)")


@dataclass(frozen=True)
class GaussianInput:
    """Ordered ordinary-Gaussian geometry in Angstrom."""

    elements: tuple[str, ...]
    coordinates: tuple[tuple[float, float, float], ...]
    charge: int
    multiplicity: int


def _element_from_token(token: str, *, path: Path, line_number: int) -> str:
    match = _ELEMENT_RE.match(token.strip())
    if match is None:
        raise ValueError(
            f"Invalid Gaussian element token at {path}:{line_number}: {token!r}"
        )
    raw = match.group(1)
    if raw.isdigit():
        atomic_number = int(raw)
        if not 0 < atomic_number < len(chemical_symbols):
            raise ValueError(
                f"Invalid Gaussian atomic number at {path}:{line_number}: {raw}"
            )
        return chemical_symbols[atomic_number]
    symbol = raw[0].upper() + raw[1:].lower()
    if symbol not in chemical_symbols:
        raise ValueError(
            f"Unknown Gaussian element at {path}:{line_number}: {raw!r}"
        )
    return symbol


def parse_gaussian_input(path: Path | str) -> GaussianInput:
    """Parse an ordinary Gaussian route/title/charge/Cartesian input."""

    source = Path(path)
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines()

    try:
        route_start = next(
            index for index, line in enumerate(lines) if line.lstrip().startswith("#")
        )
    except StopIteration as exc:
        raise ValueError(f"Gaussian route section not found in {source}.") from exc

    index = route_start + 1
    while index < len(lines) and lines[index].strip():
        index += 1
    while index < len(lines) and not lines[index].strip():
        index += 1
    title_start = index
    while index < len(lines) and lines[index].strip():
        index += 1
    if index == title_start:
        raise ValueError(f"Gaussian title section not found in {source}.")
    while index < len(lines) and not lines[index].strip():
        index += 1

    if index >= len(lines):
        raise ValueError(f"Gaussian charge/multiplicity line not found in {source}.")
    charge_match = _TWO_INT_RE.match(lines[index])
    if charge_match is None:
        raise ValueError(
            f"Expected an ordinary two-integer Gaussian charge/multiplicity line "
            f"at {source}:{index + 1}."
        )
    charge = int(charge_match.group(1))
    multiplicity = int(charge_match.group(2))
    if multiplicity < 1:
        raise ValueError(
            f"Gaussian multiplicity must be positive at {source}:{index + 1}."
        )

    elements: list[str] = []
    coordinates: list[tuple[float, float, float]] = []
    index += 1
    while index < len(lines) and not lines[index].strip():
        index += 1
    while index < len(lines) and lines[index].strip():
        line_number = index + 1
        fields = [part for part in re.split(r"[\s,]+", lines[index].strip()) if part]
        if len(fields) == 4:
            coordinate_fields = fields[1:4]
        elif len(fields) == 5 and re.fullmatch(r"[+-]?\d+", fields[1]):
            coordinate_fields = fields[2:5]
        else:
            raise ValueError(
                f"Unsupported ordinary Gaussian coordinate row at "
                f"{source}:{line_number}: {lines[index]!r}"
            )
        try:
            xyz = tuple(float(value) for value in coordinate_fields)
        except ValueError as exc:
            raise ValueError(
                f"Invalid Gaussian coordinates at {source}:{line_number}."
            ) from exc
        if len(xyz) != 3 or not all(math.isfinite(value) for value in xyz):
            raise ValueError(
                f"Non-finite Gaussian coordinates at {source}:{line_number}."
            )
        elements.append(
            _element_from_token(fields[0], path=source, line_number=line_number)
        )
        coordinates.append((xyz[0], xyz[1], xyz[2]))
        index += 1

    if not elements:
        raise ValueError(f"No ordinary Gaussian coordinate rows found in {source}.")
    return GaussianInput(tuple(elements), tuple(coordinates), charge, multiplicity)


__all__ = ["GaussianInput", "parse_gaussian_input"]
