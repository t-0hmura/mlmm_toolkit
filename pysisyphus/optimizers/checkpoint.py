"""Safe, atomic, explicitly bounded optimizer restart checkpoints.

The legacy ``dump_restart_info`` serialized live NumPy state with
``yaml.dump``, emitting ``!!python/object/apply:numpy...`` tags that the
constructor's ``yaml.SafeLoader`` then refuses to load — a checkpoint that
cannot be reloaded by its own loader.  It also captured the geometry *before*
the proposed step was committed, so the restored state described a different
transition phase than its own histories.

This module supplies a lower-engine checkpoint envelope that is:

* **safe** — only YAML-primitive scalars/mappings/lists reach disk;
* **atomic** — staged to a sibling, flushed, ``fsync``-ed, then ``os.replace``-d,
  so a failed write leaves any existing checkpoint byte-identical;
* **post-step** — captured after the actual transformed step is known, so a
  resumed optimizer continues the same trajectory;
* **explicitly bounded** — only classes whose complete resumable state is
  declared (LBFGS and non-TS Hessian/RF optimizers) are accepted; every other
  class fails loud with a typed error rather than writing a partial or
  approximate checkpoint;
* **validate-before-mutate** — a load validates schema, version, phase,
  optimizer id, geometry identity, required keys, finiteness, aligned history
  lengths, the nested geometry payload and any restart Hessian *completely*
  before any optimizer or geometry state changes.

The implementation stays inside the numerical engine and does not depend on a
product-layer checkpoint writer.
"""

from __future__ import annotations

import os
import tempfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np
import yaml

try:  # torch is a hard dependency of the engine, but stay defensive.
    import torch
except Exception:  # pragma: no cover - torch is always present in practice
    torch = None


CHECKPOINT_SCHEMA = "pysisyphus.optimizer-checkpoint"
CHECKPOINT_VERSION = 1
RESUMABLE_PHASE = "post_step"

# Leaf classes whose complete resumable state declaration and resume-equivalence
# controls are in place.  Default is *unsupported*.
_SUPPORTED_LEAF_CLASSES = frozenset({"RFOptimizer", "LBFGS", "HessianOptimizer"})

# Any class carrying one of these in its MRO is TS mode-following / chain-of-
# states / line-search state that is not completely declared here.
_UNSUPPORTED_BASE_CLASSES = frozenset(
    {
        "TSHessianOptimizer",
        "StringOptimizer",
        "BacktrackingOptimizer",
        "ChainOfStates",
    }
)

_REQUIRED_RESTART_KEYS = (
    "last_cycle",
    "max_cycles",
    "energies",
    "coords",
    "forces",
    "steps",
    "geom_info",
)


class CheckpointError(Exception):
    """Base class for all checkpoint failures."""


class CheckpointUnsupportedError(CheckpointError):
    """The optimizer class has no complete resumable-state declaration."""


class CheckpointSerializationError(CheckpointError):
    """A payload contained state that cannot be safely serialized."""


class CheckpointValidationError(CheckpointError):
    """A checkpoint file failed validation before any state was mutated."""


def optimizer_id(optimizer: Any) -> str:
    cls = type(optimizer)
    return f"{cls.__module__}.{cls.__qualname__}"


def _mro_names(optimizer: Any) -> set[str]:
    return {cls.__name__ for cls in type(optimizer).__mro__}


def assert_supported(optimizer: Any) -> None:
    """Raise :class:`CheckpointUnsupportedError` for any non-declared class."""

    if getattr(optimizer, "is_cos", False):
        raise CheckpointUnsupportedError(
            f"{optimizer_id(optimizer)} is a chain-of-states optimizer; its "
            "image/history state is not declared for safe checkpointing."
        )
    unsupported = _mro_names(optimizer) & _UNSUPPORTED_BASE_CLASSES
    if unsupported:
        raise CheckpointUnsupportedError(
            f"{optimizer_id(optimizer)} carries undeclared state via "
            f"{sorted(unsupported)}; refusing to write a partial checkpoint."
        )
    leaf = type(optimizer).__name__
    if leaf not in _SUPPORTED_LEAF_CLASSES:
        raise CheckpointUnsupportedError(
            f"{optimizer_id(optimizer)} is not an explicitly supported "
            "checkpoint class; complete its state declaration first."
        )


def to_safe_primitives(obj: Any) -> Any:
    """Recursively convert to YAML-primitive types; reject unsafe state."""

    # NumPy scalar types must be handled before the plain-``float``/``int``
    # branch: ``np.float64`` is a subclass of Python ``float`` and would
    # otherwise pass through unconverted and break ``yaml.safe_dump``.
    if obj is None or isinstance(obj, (bool, str)):
        return obj
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        if obj.dtype == object:
            raise CheckpointSerializationError(
                "object-dtype arrays cannot be safely serialized"
            )
        return obj.tolist()
    if isinstance(obj, (int, float)):
        return obj
    if torch is not None and isinstance(obj, torch.Tensor):
        try:
            return obj.detach().cpu().tolist()
        except Exception as exc:  # pragma: no cover - defensive
            raise CheckpointSerializationError(
                f"tensor could not be copied for serialization: {exc}"
            ) from exc
    if isinstance(obj, Mapping):
        return {str(key): to_safe_primitives(val) for key, val in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_safe_primitives(val) for val in obj]
    if callable(obj):
        raise CheckpointSerializationError(
            f"callable state {obj!r} cannot be serialized"
        )
    raise CheckpointSerializationError(
        f"unsupported checkpoint payload type: {type(obj)!r}"
    )


def _geometry_identity(geometry: Any) -> dict[str, Any]:
    return {
        "atoms": [str(atom) for atom in getattr(geometry, "atoms", [])],
        "coord_type": str(getattr(geometry, "coord_type", "cart")),
        "cart_size": int(np.asarray(geometry.cart_coords).size),
        "freeze_atoms": sorted(
            int(i) for i in getattr(geometry, "freeze_atoms", [])
        ),
    }


def build_envelope(optimizer: Any, *, phase: str = RESUMABLE_PHASE) -> dict[str, Any]:
    """Assemble a safe-primitive checkpoint envelope for a supported optimizer."""

    assert_supported(optimizer)
    restart_info = optimizer.get_restart_info()
    envelope = {
        "schema": CHECKPOINT_SCHEMA,
        "version": CHECKPOINT_VERSION,
        "phase": str(phase),
        "optimizer_id": optimizer_id(optimizer),
        "geometry_identity": _geometry_identity(optimizer.geometry),
        "cycle": {
            "completed_cycle": int(getattr(optimizer, "cur_cycle", 0)),
            "max_cycles": (
                None
                if getattr(optimizer, "max_cycles", None) is None
                else int(optimizer.max_cycles)
            ),
        },
        "restart_info": restart_info,
    }
    return to_safe_primitives(envelope)


def atomic_write_yaml(path: os.PathLike | str, payload: Mapping) -> Path:
    """Serialize *payload* and atomically replace *path*.

    Serialization, staging, flush/fsync, or replacement failure leaves any
    existing file byte-identical and removes the temporary sibling.
    """

    destination = Path(path)
    try:
        text = yaml.safe_dump(dict(payload), default_flow_style=False, sort_keys=True)
    except yaml.YAMLError as exc:
        raise CheckpointSerializationError(
            f"checkpoint payload is not YAML-safe: {exc}"
        ) from exc

    parent = destination.parent
    fd, tmp_name = tempfile.mkstemp(
        dir=str(parent), prefix=destination.name + ".", suffix=".tmp"
    )
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_path, destination)
    except BaseException:
        try:
            tmp_path.unlink()
        except OSError:
            pass
        raise
    return destination


def save_checkpoint(
    optimizer: Any,
    path: os.PathLike | str,
    *,
    phase: str = RESUMABLE_PHASE,
) -> Path:
    """Write one supported optimizer's post-step checkpoint atomically."""

    envelope = build_envelope(optimizer, phase=phase)
    return atomic_write_yaml(path, envelope)


def _validate_history_lengths(restart_info: Mapping) -> None:
    coords = restart_info["coords"]
    forces = restart_info["forces"]
    energies = restart_info["energies"]
    steps = restart_info["steps"]
    for name, value in (
        ("coords", coords),
        ("forces", forces),
        ("energies", energies),
        ("steps", steps),
    ):
        if not isinstance(value, list):
            raise CheckpointValidationError(f"history {name!r} is not a list")
    n = len(coords)
    if n == 0:
        raise CheckpointValidationError("coords history is empty")
    if len(forces) != n or len(energies) != n:
        raise CheckpointValidationError(
            "coords/forces/energies history lengths are inconsistent"
        )
    if len(steps) not in (n, n - 1):
        raise CheckpointValidationError(
            "steps history length is inconsistent with coords"
        )


def _validate_finite(restart_info: Mapping) -> None:
    for name in ("coords", "forces", "steps", "energies"):
        for entry in restart_info[name]:
            arr = np.asarray(entry, dtype=float)
            if arr.size and not np.all(np.isfinite(arr)):
                raise CheckpointValidationError(
                    f"non-finite value found in {name!r} history"
                )


def _validate_geom_info(restart_info: Mapping, optimizer: Any) -> None:
    """Validate the nested geometry payload before any state is applied.

    ``Geometry.set_restart_info`` asserts its atoms and then overwrites the
    Cartesian coordinates, so a corrupted nested payload has to be rejected
    here, together with the outer identity it must agree with.
    """

    geom_info = restart_info["geom_info"]
    if not isinstance(geom_info, Mapping):
        raise CheckpointValidationError("checkpoint geom_info is not a mapping")
    identity = _geometry_identity(optimizer.geometry)
    atoms = geom_info.get("atoms")
    if not isinstance(atoms, (list, tuple)):
        raise CheckpointValidationError("checkpoint geom_info atoms are missing")
    if [str(atom) for atom in atoms] != identity["atoms"]:
        raise CheckpointValidationError(
            "checkpoint geom_info atoms do not match the target geometry"
        )
    coord_type = geom_info.get("coord_type")
    if str(coord_type) != identity["coord_type"]:
        raise CheckpointValidationError(
            f"checkpoint geom_info coord_type {coord_type!r} does not match the "
            f"target geometry coord_type {identity['coord_type']!r}"
        )
    # Nested cart_coords validation: translate nonnumeric input to checkpoint error,
    # require 1D shape before any mutation.
    try:
        cart_coords = np.asarray(geom_info.get("cart_coords", []), dtype=float)
    except (ValueError, TypeError) as exc:
        raise CheckpointValidationError(
            f"checkpoint geom_info cart_coords contains non-numeric data: {exc}"
        ) from exc
    if cart_coords.ndim != 1:
        raise CheckpointValidationError(
            f"checkpoint geom_info cart_coords must be 1D, got shape {cart_coords.shape}"
        )
    if cart_coords.size != identity["cart_size"]:
        raise CheckpointValidationError(
            f"checkpoint geom_info has {cart_coords.size} Cartesian coordinates, "
            f"expected {identity['cart_size']}"
        )
    if not np.all(np.isfinite(cart_coords)):
        raise CheckpointValidationError(
            "non-finite value found in checkpoint geom_info cart_coords"
        )


def _restart_hessian_shape(restart_info: Mapping):
    """Shape recorded beside a restart Hessian, or ``None`` when absent."""

    spec = restart_info.get("H_spec")
    if not isinstance(spec, Mapping):
        return None
    shape = spec.get("shape")
    if not isinstance(shape, (list, tuple)):
        return None
    try:
        return tuple(int(dim) for dim in shape)
    except (TypeError, ValueError):
        raise CheckpointValidationError(
            f"checkpoint Hessian shape {shape!r} is not a sequence of integers"
        )


def _validate_hessian(restart_info: Mapping, optimizer: Any) -> None:
    """Validate a Hessian restart matrix against the target's active space."""

    if "H" not in getattr(optimizer, "required_opt_restart_keys", ()):
        return
    H = np.asarray(restart_info["H"], dtype=float)
    if H.size and not np.all(np.isfinite(H)):
        raise CheckpointValidationError("non-finite value found in checkpoint Hessian")
    # Nested lists lose the shape of an empty array, so a (0, 0) Hessian arrives
    # as []. Prefer the shape recorded with the payload; a checkpoint written
    # before that record is only accepted as (0, 0) below, where the resolved
    # active dimension is zero.
    recorded = _restart_hessian_shape(restart_info)
    if recorded is not None:
        if len(recorded) != 2 or recorded[0] != recorded[1]:
            raise CheckpointValidationError(
                f"checkpoint Hessian records shape {recorded}, expected a square "
                "matrix"
            )
        if int(recorded[0]) * int(recorded[1]) != H.size:
            raise CheckpointValidationError(
                f"checkpoint Hessian has {H.size} values, which does not match "
                f"its recorded shape {recorded}"
            )
        actual_dim = int(recorded[0])
    elif H.size == 0:
        actual_dim = 0
    else:
        if H.ndim != 2 or H.shape[0] != H.shape[1]:
            raise CheckpointValidationError(
                f"checkpoint Hessian has shape {H.shape}, expected a square matrix"
            )
        actual_dim = int(H.shape[0])
    # Use already-resolved active_dof_indices first so wrong-size active Hessians
    # cannot skip validation.
    expected = None
    if getattr(optimizer, "using_active_dofs", False):
        active_indices = getattr(optimizer, "active_dof_indices", None)
        if active_indices is not None:
            try:
                expected = int(len(np.asarray(active_indices, dtype=int)))
            except Exception:  # pragma: no cover - defensive
                pass
    # Fallback to full coordinate dimension.
    if expected is None:
        try:
            expected = int(np.asarray(optimizer.geometry.coords).size)
        except Exception:  # pragma: no cover - defensive
            expected = None
    if expected is not None and actual_dim != expected:
        raise CheckpointValidationError(
            f"checkpoint Hessian dimension {actual_dim} does not match the "
            f"resolved active dimension {expected}"
        )


def validate_payload(payload: Any, optimizer: Any) -> dict[str, Any]:
    """Validate a loaded checkpoint completely; return its ``restart_info``.

    Raises :class:`CheckpointValidationError` (or
    :class:`CheckpointUnsupportedError`) before any state is mutated.
    """

    if not isinstance(payload, Mapping):
        raise CheckpointValidationError("checkpoint payload is not a mapping")
    if payload.get("schema") != CHECKPOINT_SCHEMA:
        raise CheckpointValidationError(
            "missing or unknown checkpoint schema; refusing a legacy/unsafe "
            "restart payload"
        )
    if payload.get("version") != CHECKPOINT_VERSION:
        raise CheckpointValidationError(
            f"unsupported checkpoint version {payload.get('version')!r}"
        )
    if payload.get("phase") != RESUMABLE_PHASE:
        raise CheckpointValidationError(
            f"checkpoint phase {payload.get('phase')!r} is not resumable"
        )
    assert_supported(optimizer)
    if payload.get("optimizer_id") != optimizer_id(optimizer):
        raise CheckpointValidationError(
            f"checkpoint optimizer id {payload.get('optimizer_id')!r} does not "
            f"match target {optimizer_id(optimizer)!r}"
        )
    expected_geometry = _geometry_identity(optimizer.geometry)
    if payload.get("geometry_identity") != expected_geometry:
        raise CheckpointValidationError(
            "checkpoint geometry identity does not match the target geometry"
        )
    restart_info = payload.get("restart_info")
    if not isinstance(restart_info, Mapping):
        raise CheckpointValidationError("checkpoint is missing restart_info")
    for key in _REQUIRED_RESTART_KEYS:
        if key not in restart_info:
            raise CheckpointValidationError(
                f"checkpoint restart_info is missing required key {key!r}"
            )
    # Transactional load: also validate the subclass-specific restart keys the
    # target optimizer's ``_set_opt_restart_info`` reads unconditionally (H for
    # Hessian/RF optimizers, coord_diffs/... for LBFGS).  Without this a
    # checkpoint missing a subclass key passed validation and only failed with
    # a bare ``KeyError`` *after* ``Optimizer.set_restart_info`` had already
    # overwritten coords/energies/forces/steps -- a partial mutation with an
    # untyped error.  Validating here keeps the "reject before any state
    # change" promise for subclass state too.
    for key in getattr(optimizer, "required_opt_restart_keys", ()):
        if key not in restart_info:
            raise CheckpointValidationError(
                f"checkpoint restart_info is missing subclass-required key "
                f"{key!r} for {optimizer_id(optimizer)}"
            )
    if bool(getattr(optimizer, "reject_uphill", False)) and "cart_coords" not in restart_info:
        raise CheckpointValidationError(
            "checkpoint restart_info lacks Cartesian history required for "
            "uphill-trial rollback"
        )
    _validate_history_lengths(restart_info)
    _validate_finite(restart_info)
    _validate_geom_info(restart_info, optimizer)
    _validate_hessian(restart_info, optimizer)
    return dict(restart_info)


def load_payload(path: os.PathLike | str) -> Any:
    with open(os.fspath(path), "r") as handle:
        return yaml.safe_load(handle)


def load_and_apply(optimizer: Any, path: os.PathLike | str) -> None:
    """Load, validate, then apply a checkpoint to *optimizer*.

    Validation is complete before any mutation, so a rejected checkpoint leaves
    the optimizer and its geometry untouched.
    """

    payload = load_payload(path)
    restart_info = validate_payload(payload, optimizer)
    # ``Geometry.set_restart_info`` compares against its tuple ``atoms``; safe
    # serialization stores the ordered atoms as a list, so restore the tuple
    # before handing the payload back to the shared geometry loader.
    geom_info = restart_info.get("geom_info")
    if isinstance(geom_info, Mapping) and "atoms" in geom_info:
        geom_info = dict(geom_info)
        geom_info["atoms"] = tuple(geom_info["atoms"])
        restart_info["geom_info"] = geom_info
    optimizer.set_restart_info(restart_info)
