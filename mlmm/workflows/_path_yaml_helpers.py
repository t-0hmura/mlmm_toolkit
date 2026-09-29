"""Shared YAML-overlay helpers for `path-opt` and `path-search` cli().

Both subcommands accept the same single-structure optimizer YAML schema:

  opt:
    <common keys>     # OPT_BASE_KW keys overlay both LBFGS + RFO
    lbfgs: {...}      # LBFGS-only overrides
    rfo: {...}        # RFO-only overrides

The ``lbfgs`` and ``rfo`` sections may also be written at the top level or
under ``stopt:``; they configure the single-structure optimizers and never
reach StringOptimizer. Within one YAML layer, a key set to different values
in two of these places is rejected instead of one value silently winning.
"""

from __future__ import annotations

from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import click

_SINGLE_OPT_KINDS: Tuple[str, ...] = ("lbfgs", "rfo")
# The path workflows set these per segment, so differing YAML values are not a conflict.
_PER_RUN_OPT_KEYS = frozenset({"out_dir", "prefix"})


def _mapping_section(
    layer_cfg: Mapping[str, Any], path: Sequence[str]
) -> Optional[Dict[str, Any]]:
    """Return a YAML section; ``None`` when absent or empty, an error when not a mapping."""
    path = tuple(path)
    cur: Any = layer_cfg
    for depth, key in enumerate(path):
        if not isinstance(cur, Mapping):
            raise click.BadParameter(
                f"YAML section '{'.'.join(path[:depth])}' must be a mapping, "
                f"got {type(cur).__name__}."
            )
        if key not in cur:
            return None
        cur = cur[key]
        if cur is None:
            return None
    if not isinstance(cur, dict):
        raise click.BadParameter(
            f"YAML section '{'.'.join(path)}' must be a mapping, got "
            f"{type(cur).__name__}."
        )
    return cur


def _merged_kind_section(
    layer_cfg: Mapping[str, Any], kind: str
) -> Tuple[Dict[str, Any], Dict[str, str]]:
    """Merge ``<kind>``, ``opt.<kind>`` and ``stopt.<kind>`` of one layer.

    Returns the merged keys and, for each key, the section label it came from.
    """
    merged: Dict[str, Any] = {}
    origin: Dict[str, str] = {}
    for path in ((kind,), ("opt", kind), ("stopt", kind)):
        label = ".".join(path)
        for key, value in (_mapping_section(layer_cfg, path) or {}).items():
            if key in merged and merged[key] != value:
                raise click.BadParameter(
                    f"{origin[key]}.{key} and {label}.{key} conflict."
                )
            merged[key] = value
            origin.setdefault(key, label)
    return merged, origin


def apply_single_opt_yaml_layer(
    layer_cfg: Dict[str, Any],
    *,
    lbfgs_cfg: Dict[str, Any],
    rfo_cfg: Dict[str, Any],
    stopt_cfg: Dict[str, Any],
    opt_base_kw: Mapping[str, Any],
    deep_update,
) -> None:
    """Apply single-structure optimizer overrides from one YAML layer.

    Mutates ``lbfgs_cfg`` and ``rfo_cfg`` in place and removes the nested
    ``lbfgs``/``rfo`` sections that the ``stopt`` route copied into
    ``stopt_cfg``. ``deep_update`` is taken as a parameter so this module
    stays free of upward imports from ``mlmm.core.utils``.
    """
    if not isinstance(layer_cfg, dict):
        return
    opt_section = _mapping_section(layer_cfg, ("opt",)) or {}
    _mapping_section(layer_cfg, ("stopt",))
    common_updates = {k: v for k, v in opt_section.items() if k in opt_base_kw}
    for kind, target in zip(_SINGLE_OPT_KINDS, (lbfgs_cfg, rfo_cfg)):
        stopt_cfg.pop(kind, None)
        kind_updates, _ = _merged_kind_section(layer_cfg, kind)
        if common_updates:
            deep_update(target, dict(common_updates))
        if kind_updates:
            deep_update(target, kind_updates)


def check_single_opt_yaml_conflicts(
    layer_cfg: Mapping[str, Any],
    *,
    kind: str,
    opt_base_kw: Mapping[str, Any],
) -> None:
    """Reject an ``opt.<key>`` that the same layer sets differently for ``kind``.

    Only the optimizer that runs is checked, as in the ``opt`` command. CLI
    options are not counted: they overwrite the resolved values directly.
    """
    if not isinstance(layer_cfg, Mapping):
        return
    opt_section = _mapping_section(layer_cfg, ("opt",)) or {}
    kind_section, origin = _merged_kind_section(layer_cfg, kind)
    shared = (opt_section.keys() & kind_section.keys() & opt_base_kw.keys()) - _PER_RUN_OPT_KEYS
    for key in sorted(shared):
        if opt_section[key] != kind_section[key]:
            raise click.BadParameter(
                f"opt.{key} and {origin[key]}.{key} conflict."
            )


__all__ = [
    "apply_single_opt_yaml_layer",
    "check_single_opt_yaml_conflicts",
]
