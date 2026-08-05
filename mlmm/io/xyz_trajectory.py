"""Strict, lightweight XYZ-trajectory parsing shared by CLI and notebook UIs."""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Optional


_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
_NUMBER_RE = re.compile(_NUMBER)
_KEYED_PREFIX_RE = re.compile(r"\b(?:E|Energy)\s*[:=]", re.IGNORECASE)
_KEYED_VALUE_RE = re.compile(
    rf"\b(?P<key>E|Energy)\s*[:=]\s*(?P<value>{_NUMBER})",
    re.IGNORECASE,
)
_SUPPORTED_UNIT_RE = re.compile(
    r"(?:(?P<open>[\[\(\{])\s*)?"
    r"(?P<unit>Ha|Eh|hartrees?|eV|kcal(?:/mol)?)"
    r"\s*(?P<close>[\]\)\}]?)"
    r"(?=$|[\s,;])",
    re.IGNORECASE,
)
_UNIT_ASSIGN_RE = re.compile(
    r"\bunits?\s*[:=]\s*(?P<unit>[^\s,;]+)",
    re.IGNORECASE,
)
_UNIT_ASSIGN_PREFIX_RE = re.compile(r"\bunits?\s*[:=]", re.IGNORECASE)
_EXTXYZ_RE = re.compile(r"\b(?:Properties|Lattice)\s*=", re.IGNORECASE)
_PYSISYPHUS_ENERGY_RE = re.compile(
    rf"^\s*(?P<value>{_NUMBER})\s*,(?:\s.*)?$"
)
_EV_PER_HARTREE = 27.211386245988
_KCAL_PER_HARTREE = 627.5094740631


def _canonical_energy_unit(unit: str) -> str:
    lowered = unit.lower()
    if lowered in {"ha", "eh", "hartree", "hartrees"}:
        return "ha"
    if lowered.startswith("kcal"):
        return "kcal/mol"
    return lowered


def parse_xyz_energy_comment(comment: str) -> tuple[Optional[float], str]:
    """Return ``(energy_hartree, provenance)`` without guessing frame counters.

    Explicit ``E=``/``Energy=`` values may carry ``Ha``/``Eh``, ``eV``, or
    ``kcal/mol``. The bundled pysisyphus trajectory form
    ``<hartree-value> , <comment>`` is also recognized. Without either marker,
    only a comment whose entire content is one floating-point value is
    accepted. This deliberately leaves mode labels, frame counters,
    frequencies, and other numeric metadata unavailable.
    """
    keyed = list(_KEYED_VALUE_RE.finditer(comment))
    if keyed:
        converted = []
        sources = []
        is_extxyz = bool(_EXTXYZ_RE.search(comment))
        for match_index, match in enumerate(keyed):
            value = float(match.group("value"))
            if not math.isfinite(value):
                return None, "nonfinite-explicit"
            tail_stop = (
                keyed[match_index + 1].start()
                if match_index + 1 < len(keyed)
                else len(comment)
            )
            tail = comment[match.end() : tail_stop]
            whitespace = len(tail) - len(tail.lstrip())
            unit_match = _SUPPORTED_UNIT_RE.match(tail[whitespace:])
            declared_unit = unit_match.group("unit") if unit_match else None
            if unit_match is not None:
                opener = unit_match.group("open") or ""
                closer = unit_match.group("close") or ""
                matching = {"[": "]", "(": ")", "{": "}"}
                if bool(opener) != bool(closer) or (
                    opener and matching[opener] != closer
                ):
                    return None, "malformed-explicit"
            assigned_units = []
            for assignment in _UNIT_ASSIGN_RE.finditer(tail):
                token = assignment.group("unit")
                if token[:1] in "[({" or token[-1:] in "])}":
                    matching = {"[": "]", "(": ")", "{": "}"}
                    if (
                        token[:1] not in matching
                        or token[-1:] != matching[token[0]]
                    ):
                        return None, "malformed-explicit"
                    token = token[1:-1]
                if not re.fullmatch(
                    r"Ha|Eh|hartrees?|eV|kcal(?:/mol)?",
                    token,
                    re.IGNORECASE,
                ):
                    return None, "unsupported-unit"
                assigned_units.append(token)
            if len(_UNIT_ASSIGN_PREFIX_RE.findall(tail)) != len(assigned_units):
                return None, "malformed-explicit"
            normalized_units = {
                _canonical_energy_unit(unit) for unit in assigned_units
            }
            if len(normalized_units) > 1:
                return None, "conflicting-explicit"
            if assigned_units:
                assigned_unit = assigned_units[0]
                normalized_assigned = next(iter(normalized_units))
                normalized_direct = (
                    _canonical_energy_unit(declared_unit)
                    if declared_unit
                    else None
                )
                if (
                    normalized_direct is not None
                    and normalized_direct != normalized_assigned
                ):
                    return None, "conflicting-explicit"
                declared_unit = assigned_unit
            if unit_match is None:
                remainder = tail[whitespace:]
                if whitespace == 0 and remainder:
                    if remainder[0].isalnum() or remainder[0] in "._+-/":
                        return None, "malformed-explicit"
                elif remainder:
                    if remainder[0] in "[({":
                        return None, (
                            "malformed-explicit"
                            if not any(mark in remainder[1:] for mark in "])}")
                            else "unsupported-unit"
                        )
                    token = re.split(r"[\s,;]+", remainder, maxsplit=1)[0]
                    if token and token[0].isalpha() and "=" not in token:
                        return None, "unsupported-unit"
            unit = (declared_unit or "Ha").lower()
            if unit == "ev":
                value /= _EV_PER_HARTREE
                source = "explicit-eV"
            elif unit.startswith("kcal"):
                value /= _KCAL_PER_HARTREE
                source = "explicit-kcal/mol"
            elif (
                declared_unit is None
                and is_extxyz
                and match.group("key").lower() == "energy"
            ):
                value /= _EV_PER_HARTREE
                source = "explicit-extxyz-eV"
            else:
                source = "explicit-Ha" if declared_unit else "explicit-assumed-Ha"
            converted.append(value)
            sources.append(source)
        if any(abs(value - converted[0]) > 1e-12 for value in converted[1:]):
            return None, "conflicting-explicit"
        return converted[0], "+".join(dict.fromkeys(sources))
    if _KEYED_PREFIX_RE.search(comment):
        return None, "malformed-explicit"

    pysisyphus_match = _PYSISYPHUS_ENERGY_RE.fullmatch(comment)
    if pysisyphus_match is not None:
        value = float(pysisyphus_match.group("value"))
        return (
            (value, "pysisyphus-legacy-Ha")
            if math.isfinite(value)
            else (None, "nonfinite-pysisyphus")
        )

    stripped = comment.strip()
    if _NUMBER_RE.fullmatch(stripped) and (
        "." in stripped or "e" in stripped.lower()
    ):
        value = float(stripped)
        return (
            (value, "bare-assumed-Ha")
            if math.isfinite(value)
            else (None, "nonfinite-bare")
        )
    tokens = _NUMBER_RE.findall(comment)
    return None, "missing" if not tokens else "ambiguous"


def read_xyz_trajectory(
    path: Path | str,
    *,
    require_energies: bool = False,
) -> dict:
    """Parse complete XYZ blocks and return frames plus Hartree energies.

    Missing or ambiguous comments yield ``None`` energies for structure-only
    viewers. ``require_energies=True`` fails closed for energy plotting.
    """
    source = Path(path)
    lines = source.read_text(encoding="utf-8").splitlines()
    frames: list[str] = []
    energies: list[Optional[float]] = []
    provenance: list[str] = []
    atom_identity: Optional[tuple[str, ...]] = None
    cursor = 0
    while cursor < len(lines):
        while cursor < len(lines) and not lines[cursor].strip():
            cursor += 1
        if cursor >= len(lines):
            break
        frame_number = len(frames) + 1
        header = lines[cursor].strip()
        # Treat first line as atom count only if it contains exactly one integer token
        header_parts = header.split()
        if len(header_parts) != 1:
            raise ValueError(
                f"Malformed XYZ header in frame {frame_number} of {source}: {header!r}"
            )
        try:
            atom_count = int(header_parts[0])
        except ValueError as exc:
            raise ValueError(
                f"Malformed XYZ header in frame {frame_number} of {source}: {header!r}"
            ) from exc
        if atom_count <= 0:
            raise ValueError(
                f"XYZ frame {frame_number} of {source} has invalid atom count {atom_count}."
            )
        stop = cursor + atom_count + 2
        if cursor + 1 >= len(lines) or stop > len(lines):
            raise ValueError(
                f"Incomplete XYZ frame {frame_number} of {source}; "
                f"expected {atom_count} coordinate rows."
            )
        comment = lines[cursor + 1]
        frame_identity: list[str] = []
        for row_number, row in enumerate(
            lines[cursor + 2 : stop], start=1
        ):
            fields = row.split()
            if len(fields) < 4:
                raise ValueError(
                    f"Malformed coordinate row {row_number} in frame "
                    f"{frame_number} of {source}."
                )
            try:
                coordinates = tuple(float(fields[index]) for index in (1, 2, 3))
            except ValueError as exc:
                raise ValueError(
                    f"Non-numeric coordinate row {row_number} in frame "
                    f"{frame_number} of {source}."
                ) from exc
            if not all(math.isfinite(value) for value in coordinates):
                raise ValueError(
                    f"Non-finite coordinate row {row_number} in frame "
                    f"{frame_number} of {source}."
                )
            frame_identity.append(fields[0].casefold())
        identity = tuple(frame_identity)
        if atom_identity is None:
            atom_identity = identity
        elif identity != atom_identity:
            raise ValueError(
                f"XYZ atom identity/order changes in frame {frame_number} of {source}."
            )
        energy, energy_source = parse_xyz_energy_comment(comment)
        if require_energies and energy is None:
            if energy_source == "conflicting-explicit":
                reason = "Conflicting explicit XYZ energies"
            elif energy_source.startswith("nonfinite-"):
                reason = "Non-finite XYZ energy"
            elif energy_source == "unsupported-unit":
                reason = "Unsupported XYZ energy unit"
            elif energy_source == "malformed-explicit":
                reason = "Malformed explicit XYZ energy"
            elif energy_source == "ambiguous" and _NUMBER_RE.fullmatch(comment.strip()):
                reason = "Ambiguous integer-only XYZ comment"
            else:
                reason = "Missing or ambiguous XYZ energy"
            raise RuntimeError(
                f"{reason} in frame {frame_number}; write E=<value> "
                f"with an optional unit: {comment!r}"
            )
        frames.append("\n".join(lines[cursor:stop]) + "\n")
        energies.append(energy)
        provenance.append(energy_source)
        cursor = stop
    if not frames:
        raise RuntimeError(f"No XYZ frames found in {source}")
    return {
        "frames": frames,
        "energies_ha": energies,
        "energy_provenance": provenance,
        "energy_unit": "hartree",
    }
