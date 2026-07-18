"""Pure alternate-location scoring shared by MLMM structure readers.

PDB and mmCIF both allow an occupancy field to be absent.  Missing occupancy
is not the same observation as a parsed zero: any parsed mean ranks above an
all-missing conformer, and ties are resolved by first appearance.
"""

from __future__ import annotations

from collections.abc import Iterable
import math
from typing import Optional


AltlocObservation = tuple[str, Optional[float], int]


def parsed_occupancy(value: object) -> Optional[float]:
    """Return a finite occupancy value, or ``None`` when it is unavailable."""

    try:
        occupancy = float(value)
    except (TypeError, ValueError):
        return None
    return occupancy if math.isfinite(occupancy) else None


def occupancy_choice_key(
    occupancy: Optional[float],
    first_index: int,
) -> tuple[int, float, int]:
    """Ranking key where parsed zero beats missing and earliest wins ties."""

    value = parsed_occupancy(occupancy)
    return (
        int(value is not None),
        float(value) if value is not None else float("-inf"),
        -int(first_index),
    )


def choose_altloc_label(observations: Iterable[AltlocObservation]) -> str:
    """Choose one non-blank label from ordered occupancy observations.

    Means include only parsed finite values.  A label with no parsed values
    ranks below every label with a parsed value, including a mean of zero.
    """

    stats: dict[str, list[float | int]] = {}
    for raw_label, raw_occupancy, raw_index in observations:
        label = str(raw_label).strip()
        if not label:
            continue
        if label not in stats:
            stats[label] = [0.0, 0, int(raw_index)]
        value = parsed_occupancy(raw_occupancy)
        if value is not None:
            stats[label][0] = float(stats[label][0]) + value
            stats[label][1] = int(stats[label][1]) + 1

    if not stats:
        return ""

    def score(item: tuple[str, list[float | int]]) -> tuple[int, float, int]:
        _label, (occupancy_sum, occupancy_count, first_index) = item
        count = int(occupancy_count)
        mean = float(occupancy_sum) / count if count else None
        return occupancy_choice_key(mean, int(first_index))

    return max(stats.items(), key=score)[0]


__all__ = [
    "AltlocObservation",
    "choose_altloc_label",
    "occupancy_choice_key",
    "parsed_occupancy",
]
