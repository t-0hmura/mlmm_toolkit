"""Small, restart-safe caches for MEP reference-mode candidates.

The file cache is authoritative across child processes and restarts.  A tiny
process-local CPU cache only avoids repeatedly decoding the same NPZ payload.
No GPU tensors or dense Hessians are stored here.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import os
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence

import numpy as np


_SCHEMA_VERSION = 1
_MAX_CANDIDATES = 3
_DIRECTION_DUPLICATE_OVERLAP = 1.0 - 1.0e-8
_RAM_CACHE: dict[str, tuple[tuple[int, int, int, int], "PathModeCache"]] = {}
_FILE_HASH_CACHE: dict[str, tuple[tuple[int, int], str]] = {}


@dataclass(frozen=True)
class PathModeCache:
    path: Path
    directions: np.ndarray
    labels: tuple[str, ...]
    metadata: dict[str, Any]

    @property
    def primary(self) -> np.ndarray:
        return np.asarray(self.directions[0], dtype=float)


def _unit(vector: np.ndarray) -> Optional[np.ndarray]:
    array = np.asarray(vector, dtype=np.float64).reshape(-1)
    norm = float(np.linalg.norm(array))
    if array.size == 0 or not np.all(np.isfinite(array)) or not np.isfinite(norm) or norm <= 0.0:
        return None
    return array / norm


def build_path_mode_candidates(
    coordinates: Sequence[np.ndarray],
    image_index: int,
    *,
    energies: Optional[Sequence[float]] = None,
    max_candidates: int = _MAX_CANDIDATES,
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Build energy-upwind, incoming, and outgoing HEI directions.

    Sign-equivalent directions are considered the same mode and removed using
    ``abs(dot)``.  The primary candidate matches the improved path tangent used
    by the path optimizer when finite energies are available.
    """
    arrays = [np.asarray(item, dtype=np.float64).reshape(-1) for item in coordinates]
    index = int(image_index)
    if len(arrays) < 2 or not 0 <= index < len(arrays):
        return np.empty((0, 0), dtype=np.float64), ()
    width = arrays[0].size
    if width == 0 or any(item.size != width or not np.all(np.isfinite(item)) for item in arrays):
        return np.empty((0, 0), dtype=np.float64), ()

    incoming_raw = arrays[index] - arrays[index - 1] if index > 0 else None
    outgoing_raw = arrays[index + 1] - arrays[index] if index + 1 < len(arrays) else None
    incoming = _unit(incoming_raw) if incoming_raw is not None else None
    outgoing = _unit(outgoing_raw) if outgoing_raw is not None else None

    primary: Optional[np.ndarray] = None
    primary_label = "secant-bisector"
    if index == 0:
        primary, primary_label = outgoing, "outgoing"
    elif index == len(arrays) - 1:
        primary, primary_label = incoming, "incoming"
    elif incoming is None:
        primary, primary_label = outgoing, "outgoing"
    elif outgoing is None:
        primary, primary_label = incoming, "incoming"
    else:
        if energies is not None and len(energies) == len(arrays):
            values = np.asarray(energies, dtype=np.float64)
            previous, current, following = (
                float(values[index - 1]),
                float(values[index]),
                float(values[index + 1]),
            )
            if np.all(np.isfinite((previous, current, following))):
                if following > current > previous:
                    tangent_raw = outgoing_raw
                elif following < current < previous:
                    tangent_raw = incoming_raw
                else:
                    next_delta = abs(following - current)
                    previous_delta = abs(previous - current)
                    delta_max = max(next_delta, previous_delta)
                    delta_min = min(next_delta, previous_delta)
                    if following >= previous:
                        tangent_raw = outgoing_raw * delta_max + incoming_raw * delta_min
                    else:
                        tangent_raw = outgoing_raw * delta_min + incoming_raw * delta_max
                primary = _unit(tangent_raw)
                primary_label = "energy-upwind"
        if primary is None:
            primary = _unit(incoming + outgoing)
        if primary is None:
            primary = _unit(arrays[index + 1] - arrays[index - 1])
            primary_label = "endpoint-chord"

    candidates: list[tuple[str, np.ndarray]] = []
    for label, candidate in (
        (primary_label, primary),
        ("incoming", incoming),
        ("outgoing", outgoing),
    ):
        if candidate is None:
            continue
        candidate = _unit(candidate)
        if candidate is None:
            continue
        if any(abs(float(np.dot(candidate, existing))) >= _DIRECTION_DUPLICATE_OVERLAP for _, existing in candidates):
            continue
        candidates.append((label, candidate))
        if len(candidates) >= max(1, int(max_candidates)):
            break

    if not candidates:
        return np.empty((0, width), dtype=np.float64), ()
    return (
        np.stack([candidate for _, candidate in candidates]).astype(np.float64, copy=False),
        tuple(label for label, _ in candidates),
    )


def _file_sha256(path: Optional[Path]) -> Optional[str]:
    if path is None:
        return None
    path = Path(path)
    if not path.is_file():
        return None
    stat = path.stat()
    key = str(path.resolve())
    signature = (int(stat.st_size), int(stat.st_mtime_ns))
    cached = _FILE_HASH_CACHE.get(key)
    if cached is not None and cached[0] == signature:
        return cached[1]
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    value = digest.hexdigest()
    _FILE_HASH_CACHE[key] = (signature, value)
    return value


def _metadata_path(cache_path: Path) -> Path:
    return Path(cache_path).with_suffix(".json")


def _atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def write_path_mode_cache(
    cache_path: Path,
    coordinates: Sequence[np.ndarray],
    image_index: int,
    *,
    energies: Optional[Sequence[float]] = None,
    trajectory_path: Optional[Path] = None,
    hei_path: Optional[Path] = None,
    atom_numbers: Optional[Iterable[int]] = None,
    primary_text_path: Optional[Path] = None,
    source: Optional[str] = None,
) -> Optional[PathModeCache]:
    """Atomically write a path-mode NPZ + JSON pair and optional legacy TXT."""
    cache_path = Path(cache_path)
    if cache_path.suffix.lower() != ".npz":
        cache_path = cache_path.with_suffix(".npz")
    directions, labels = build_path_mode_candidates(
        coordinates, image_index, energies=energies
    )
    if directions.size == 0:
        return None

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = cache_path.with_name(f".{cache_path.stem}.{os.getpid()}.tmp.npz")
    np.savez_compressed(
        temporary,
        directions=np.asarray(directions, dtype=np.float64),
        labels=np.asarray(labels, dtype="U64"),
        image_index=np.asarray([int(image_index)], dtype=np.int64),
    )
    os.replace(temporary, cache_path)

    energy_values: Optional[list[float]] = None
    if energies is not None:
        try:
            values = np.asarray(energies, dtype=np.float64)
            if values.ndim == 1 and np.all(np.isfinite(values)):
                energy_values = [float(value) for value in values]
        except (TypeError, ValueError):
            energy_values = None

    metadata: dict[str, Any] = {
        "schema_version": _SCHEMA_VERSION,
        "vector_space": "full_cartesian",
        "dtype": "float64",
        "candidate_count": int(directions.shape[0]),
        "vector_size": int(directions.shape[1]),
        "image_index": int(image_index),
        "candidate_labels": list(labels),
        "candidate_overlap_equivalence": "abs(dot)",
        "source": source,
        "trajectory": None if trajectory_path is None else str(Path(trajectory_path)),
        "trajectory_sha256": _file_sha256(trajectory_path),
        "hei": None if hei_path is None else str(Path(hei_path)),
        "hei_sha256": _file_sha256(hei_path),
        "atom_numbers": None if atom_numbers is None else [int(value) for value in atom_numbers],
        "energies": energy_values,
    }
    _atomic_write_json(_metadata_path(cache_path), metadata)

    if primary_text_path is not None:
        primary_text_path = Path(primary_text_path)
        primary_text_path.parent.mkdir(parents=True, exist_ok=True)
        temporary_text = primary_text_path.with_name(
            f".{primary_text_path.name}.{os.getpid()}.tmp"
        )
        np.savetxt(temporary_text, directions[0], fmt="%.17e")
        os.replace(temporary_text, primary_text_path)

    _RAM_CACHE.pop(str(cache_path.resolve()), None)
    return load_path_mode_cache(cache_path)


def _cache_signature(cache_path: Path, metadata_path: Path) -> tuple[int, int, int, int]:
    cache_stat = cache_path.stat()
    metadata_stat = metadata_path.stat()
    return (
        int(cache_stat.st_size),
        int(cache_stat.st_mtime_ns),
        int(metadata_stat.st_size),
        int(metadata_stat.st_mtime_ns),
    )


def load_path_mode_cache(
    cache_path: Path,
    *,
    trajectory_path: Optional[Path] = None,
    hei_path: Optional[Path] = None,
    expected_size: Optional[int] = None,
    atom_numbers: Optional[Iterable[int]] = None,
) -> PathModeCache:
    """Load and validate a path-mode cache.

    Geometry/trajectory mismatches are hard failures so callers can recompute
    rather than silently handing a stale direction to TS optimization.
    """
    cache_path = Path(cache_path)
    metadata_path = _metadata_path(cache_path)
    if not cache_path.is_file() or not metadata_path.is_file():
        raise FileNotFoundError(f"Path-mode cache pair is incomplete: {cache_path}")

    key = str(cache_path.resolve())
    signature = _cache_signature(cache_path, metadata_path)
    cached = _RAM_CACHE.get(key)
    if cached is not None and cached[0] == signature:
        result = cached[1]
    else:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if int(metadata.get("schema_version", -1)) != _SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported path-mode cache schema: {metadata.get('schema_version')!r}"
            )
        with np.load(cache_path, allow_pickle=False) as payload:
            directions = np.asarray(payload["directions"], dtype=np.float64)
            labels_raw = np.asarray(payload["labels"]).reshape(-1)
        if directions.ndim == 1:
            directions = directions.reshape(1, -1)
        if directions.ndim != 2 or directions.shape[0] < 1:
            raise ValueError("Path-mode cache contains no candidate directions")
        labels = tuple(str(value) for value in labels_raw.tolist())
        if len(labels) != directions.shape[0]:
            raise ValueError("Path-mode cache labels do not match its directions")
        normalized: list[np.ndarray] = []
        normalized_labels: list[str] = []
        for label, direction in zip(labels, directions):
            unit = _unit(direction)
            if unit is None:
                raise ValueError(f"Invalid path-mode candidate {label!r}")
            if any(abs(float(np.dot(unit, prior))) >= _DIRECTION_DUPLICATE_OVERLAP for prior in normalized):
                continue
            normalized.append(unit)
            normalized_labels.append(label)
        if not normalized:
            raise ValueError("Path-mode cache contains no usable candidate direction")
        result = PathModeCache(
            path=cache_path,
            directions=np.stack(normalized),
            labels=tuple(normalized_labels),
            metadata=metadata,
        )
        _RAM_CACHE[key] = (signature, result)

    if expected_size is not None and result.directions.shape[1] != int(expected_size):
        raise ValueError(
            "Path-mode vector size mismatch "
            f"({result.directions.shape[1]} cached, {int(expected_size)} expected)"
        )
    expected_numbers = None if atom_numbers is None else [int(value) for value in atom_numbers]
    cached_numbers = result.metadata.get("atom_numbers")
    if expected_numbers is not None and cached_numbers not in (None, expected_numbers):
        raise ValueError("Path-mode atom identities/order do not match the current structure")
    if trajectory_path is not None:
        expected_hash = result.metadata.get("trajectory_sha256")
        actual_hash = _file_sha256(Path(trajectory_path))
        if not expected_hash or actual_hash != expected_hash:
            raise ValueError("Path-mode trajectory fingerprint does not match")
    if hei_path is not None:
        expected_hash = result.metadata.get("hei_sha256")
        actual_hash = _file_sha256(Path(hei_path))
        if not expected_hash or actual_hash != expected_hash:
            raise ValueError("Path-mode HEI fingerprint does not match")
    return result


def read_reference_mode_candidates(
    path: Path,
    expected_size: int,
) -> tuple[list[np.ndarray], list[str], dict[str, Any]]:
    """Read a legacy single mode or a multi-candidate NPZ cache."""
    path = Path(path)
    metadata: dict[str, Any] = {}
    if path.suffix.lower() == ".npz":
        cache = load_path_mode_cache(path, expected_size=expected_size)
        return (
            [np.asarray(row, dtype=np.float64).copy() for row in cache.directions],
            list(cache.labels),
            dict(cache.metadata),
        )
    try:
        if path.suffix.lower() == ".npy":
            raw = np.load(path, allow_pickle=False)
        else:
            raw = np.loadtxt(path)
    except (OSError, ValueError) as exc:
        raise ValueError(f"Failed to read reference mode {path}: {exc}") from exc
    array = np.asarray(raw, dtype=np.float64)
    if array.ndim == 1:
        array = array.reshape(1, -1)
    elif array.ndim > 2:
        raise ValueError("Reference-mode file must contain one vector or a 2-D vector table")
    candidates: list[np.ndarray] = []
    labels: list[str] = []
    for index, row in enumerate(array):
        unit = _unit(row)
        if row.size != int(expected_size):
            raise ValueError(
                "--ref-mode must contain one Cartesian value per degree of "
                f"freedom ({row.size} read, {int(expected_size)} expected)"
            )
        if unit is None:
            raise ValueError("--ref-mode must be a finite, non-zero vector")
        if any(abs(float(np.dot(unit, prior))) >= _DIRECTION_DUPLICATE_OVERLAP for prior in candidates):
            continue
        candidates.append(unit)
        labels.append("reference-mode" if index == 0 else f"reference-mode-{index + 1}")
    if not candidates:
        raise ValueError("Reference-mode file contains no usable direction")
    return candidates, labels, metadata
