"""Safe, atomic, explicitly bounded optimizer restart checkpoints.

This module supplies a checkpoint envelope that is:

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
  optimizer id, geometry identity, required keys, finiteness, and aligned
  history lengths *completely* before any optimizer or geometry state changes.

It introduces no product-core dependency into the numerical engine.
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
    "cart_coords",
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


def _validate_geom_info(payload: Mapping, restart_info: Mapping) -> None:
    """Validate the nested geometry payload against the checkpoint identity.

    ``Geometry.set_restart_info`` runs only after the optimizer histories have
    already been overwritten, so a corrupted ``geom_info`` must be rejected
    here, before any mutation.
    """

    geom_info = restart_info["geom_info"]
    if not isinstance(geom_info, Mapping):
        raise CheckpointValidationError("checkpoint geom_info is not a mapping")
    for key in ("atoms", "cart_coords", "coord_type"):
        if key not in geom_info:
            raise CheckpointValidationError(
                f"checkpoint geom_info is missing required key {key!r}"
            )
    identity = payload["geometry_identity"]
    atoms = geom_info["atoms"]
    if not isinstance(atoms, (list, tuple)) or [str(a) for a in atoms] != list(
        identity["atoms"]
    ):
        raise CheckpointValidationError(
            "checkpoint geom_info atoms disagree with its geometry identity"
        )
    if str(geom_info["coord_type"]) != str(identity["coord_type"]):
        raise CheckpointValidationError(
            "checkpoint geom_info coord_type disagrees with its geometry identity"
        )
    try:
        cart_coords = np.asarray(geom_info["cart_coords"], dtype=float)
    except (TypeError, ValueError) as exc:
        raise CheckpointValidationError(
            f"checkpoint geom_info cart_coords are not numeric: {exc}"
        ) from exc
    if cart_coords.ndim != 1 or cart_coords.size != int(identity["cart_size"]):
        raise CheckpointValidationError(
            "checkpoint geom_info cart_coords do not match the geometry size "
            f"{identity['cart_size']}"
        )
    if not np.all(np.isfinite(cart_coords)):
        raise CheckpointValidationError(
            "checkpoint geom_info cart_coords contain non-finite values"
        )


def _expected_hessian_dim(optimizer: Any) -> int | None:
    geometry = getattr(optimizer, "geometry", None)
    if geometry is None:
        return None
    if getattr(optimizer, "using_active_dofs", False):
        active_dofs = getattr(optimizer, "active_dof_indices", None)
        if active_dofs is not None:
            return int(np.asarray(active_dofs).size)
        partial = getattr(geometry, "within_partial_hessian", None) or {}
        return int(partial.get("active_n_dof", 0)) or None
    coords = getattr(geometry, "coords", None)
    return None if coords is None else int(np.asarray(coords).size)


def _recorded_array_shape(spec: Any) -> tuple[int, ...] | None:
    """Return the array layout recorded beside a serialized restart array."""

    if not isinstance(spec, Mapping) or spec.get("shape") is None:
        return None
    try:
        return tuple(int(size) for size in spec["shape"])
    except (TypeError, ValueError) as exc:
        raise CheckpointValidationError(
            f"checkpoint records a non-integer array shape {spec['shape']!r}"
        ) from exc


def _validate_hessian(restart_info: Mapping, optimizer: Any) -> None:
    """Validate a restart Hessian before any optimizer state is applied."""

    if "H" not in restart_info:
        return
    try:
        hessian = np.asarray(restart_info["H"], dtype=float)
    except (TypeError, ValueError) as exc:
        raise CheckpointValidationError(
            f"checkpoint Hessian is not numeric: {exc}"
        ) from exc
    # Nested lists lose the layout of an empty array, so apply the recorded
    # shape before the squareness check; a zero active dimension is valid.
    recorded = _recorded_array_shape(restart_info.get("H_backend"))
    if recorded is not None:
        try:
            hessian = hessian.reshape(recorded)
        except ValueError as exc:
            raise CheckpointValidationError(
                f"checkpoint Hessian does not fit its recorded shape "
                f"{recorded}: {exc}"
            ) from exc
    expected = _expected_hessian_dim(optimizer)
    if hessian.ndim != 2 or hessian.shape[0] != hessian.shape[1]:
        raise CheckpointValidationError(
            f"checkpoint Hessian is not square, got shape {hessian.shape}"
        )
    if not np.all(np.isfinite(hessian)):
        raise CheckpointValidationError(
            "checkpoint Hessian contains non-finite values"
        )
    if expected is not None and hessian.shape[0] != expected:
        raise CheckpointValidationError(
            f"checkpoint Hessian dimension {hessian.shape[0]} does not match "
            f"the expected dimension {expected}"
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
    _validate_history_lengths(restart_info)
    _validate_finite(restart_info)
    _validate_geom_info(payload, restart_info)
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
