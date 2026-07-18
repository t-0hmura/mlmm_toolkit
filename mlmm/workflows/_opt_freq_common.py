"""Neutral helpers shared by the ``opt`` and ``freq`` workflows.

Both ``mlmm.workflows.opt`` and ``mlmm.workflows.freq`` need the same base
geometry / calculator keyword templates and the same 1-based -> 0-based
layer / freeze normalization.  Hosting them in this lower module — instead of
``freq`` importing them from ``opt`` — is what removes the historical
``freq <-> opt`` import cycle (M39): both workflows import from here, and
``freq`` no longer imports ``opt``'s top-level surface.

This module imports only ``mlmm.core.defaults`` (foundation data) and ``click``.
"""

from __future__ import annotations

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
            return sorted({int(tok) - 1 for tok in tokens})
        except ValueError as exc:
            raise click.BadParameter(
                "geom.freeze_atoms must contain integers (string form)."
            ) from exc
    try:
        return sorted({int(idx) - 1 for idx in value})
    except TypeError as exc:
        raise click.BadParameter("geom.freeze_atoms must be iterable of integers.") from exc


def _convert_yaml_layer_atoms_1to0(calc_cfg: dict) -> None:
    """Convert 1-based YAML layer atom indices to 0-based in-place.

    Applies to calc.hess_mm_atoms, calc.movable_mm_atoms, calc.frozen_mm_atoms.
    Only converts non-None values (None = not specified in YAML).
    """
    for key in ("hess_mm_atoms", "movable_mm_atoms", "frozen_mm_atoms"):
        val = calc_cfg.get(key)
        if val is not None and not isinstance(val, str):
            try:
                calc_cfg[key] = sorted(int(i) - 1 for i in val)
            except (TypeError, ValueError):
                pass  # Leave as-is if not iterable of ints
