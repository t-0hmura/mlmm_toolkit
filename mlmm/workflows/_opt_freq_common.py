"""Neutral helpers shared by the ``opt`` and ``freq`` workflows.

Both ``mlmm.workflows.opt`` and ``mlmm.workflows.freq`` need the same base
geometry / calculator keyword templates and the same 1-based -> 0-based
layer / freeze normalization.  Hosting them in this lower module — instead of
``freq`` importing them from ``opt`` — keeps both workflows independent of
each other's top-level surfaces.

This module imports only ``mlmm.core.defaults`` (foundation data) and ``click``.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from numbers import Integral
from typing import Any, Dict, List

import click

from mlmm.core.defaults import GEOM_KW_DEFAULT, MLMM_CALC_KW

# Base keyword templates (defensive copies of the defaults).  These are the same
# single shared objects that ``opt`` and ``freq`` used before the move; per-run
# callers copy them again, so behavior is unchanged.
GEOM_KW: Dict[str, Any] = dict(GEOM_KW_DEFAULT)
CALC_KW: Dict[str, Any] = dict(MLMM_CALC_KW)


def _normalize_geom_freeze(value: Any) -> List[int]:
    """Normalize YAML-provided geom.freeze_atoms to a sorted 0-based list."""
    if value is None:
        return []
    if isinstance(value, str):
        tokens = [tok.strip() for tok in value.split(",") if tok.strip()]
        try:
            items = [int(tok) for tok in tokens]
        except ValueError as exc:
            raise click.BadParameter(
                "geom.freeze_atoms must contain integers (string form)."
            ) from exc
    else:
        if (
            isinstance(value, (bytes, Mapping))
            or not isinstance(value, Iterable)
        ):
            raise click.BadParameter(
                "geom.freeze_atoms must be iterable of integers."
            )
        items = list(value)
    if any(
        isinstance(item, bool)
        or not isinstance(item, Integral)
        or int(item) < 1
        for item in items
    ):
        raise click.BadParameter(
            "geom.freeze_atoms must contain only integers >= 1."
        )
    return sorted({int(item) - 1 for item in items})


def _convert_yaml_layer_atoms_1to0(calc_cfg: dict) -> None:
    """Convert 1-based YAML layer atom indices to 0-based in-place.

    Applies to calc.hess_mm_atoms, calc.movable_mm_atoms, calc.frozen_mm_atoms.
    Only converts non-None values (None = not specified in YAML).
    """
    for key in ("hess_mm_atoms", "movable_mm_atoms", "frozen_mm_atoms"):
        val = calc_cfg.get(key)
        if val is None:
            continue
        label = f"calc.{key}"
        if (
            isinstance(val, (str, bytes, Mapping))
            or not isinstance(val, Iterable)
        ):
            raise click.BadParameter(
                f"{label} must be a sequence of 1-based integer atom indices."
            )
        items = list(val)
        if any(
            isinstance(item, bool)
            or not isinstance(item, Integral)
            or int(item) < 1
            for item in items
        ):
            raise click.BadParameter(
                f"{label} must contain only integers >= 1."
            )
        calc_cfg[key] = sorted({int(item) - 1 for item in items})
