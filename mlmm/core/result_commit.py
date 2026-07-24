"""Atomic, exact-path commits for machine-readable MLMM artifacts."""

from __future__ import annotations

import json
import os
import tempfile
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any, BinaryIO


MLMM_RUN_ID_ENV = "MLMM_RUN_ID"


class ResultCommitError(OSError):
    """An artifact could not be staged or published at its exact path."""

    def __init__(self, phase: str, path: Path, cause: BaseException):
        self.phase = str(phase)
        self.path = Path(path)
        self.cause = cause
        super().__init__(
            getattr(cause, "errno", None),
            f"{self.phase} failed for {self.path}: {cause}",
            str(self.path),
        )


class RunIdentityError(ValueError):
    """A caller payload conflicts with the current MLMM invocation identity."""


def with_current_run_id(
    payload: Mapping[str, Any],
    *,
    run_id: str | None = None,
) -> dict[str, Any]:
    """Return a shallow payload copy carrying the current run identity.

    A caller-provided identical identity is retained.  A conflicting identity
    is rejected so an aggregate rewrite cannot relabel stale output as current.
    """

    result = dict(payload)
    current = run_id if run_id is not None else os.environ.get(MLMM_RUN_ID_ENV)
    if not current:
        return result
    current = str(current)
    supplied = result.get("run_id")
    if supplied is not None and str(supplied) != current:
        raise RunIdentityError(
            f"payload run_id {supplied!r} conflicts with current run_id {current!r}"
        )
    result.setdefault("run_id", current)
    return result


def serialize_json_bytes(payload: Mapping[str, Any]) -> bytes:
    """Serialize one immutable JSON generation exactly once."""

    try:
        return json.dumps(
            payload,
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        ).encode("utf-8")
    except Exception as exc:
        raise ResultCommitError("serialize", Path("<json-payload>"), exc) from exc


def stage_exact(path: Path, writer: Callable[[BinaryIO], None]) -> Path:
    """Stage bytes beside *path*, flushing and fsyncing before return."""

    destination = Path(path)
    temporary: Path | None = None
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        fd, raw_path = tempfile.mkstemp(
            prefix=f".{destination.name}.",
            suffix=".tmp",
            dir=destination.parent,
        )
        temporary = Path(raw_path)
        with os.fdopen(fd, "wb") as stream:
            writer(stream)
            stream.flush()
            os.fsync(stream.fileno())
        return temporary
    except Exception as exc:
        if temporary is not None:
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                pass
        if isinstance(exc, ResultCommitError):
            raise
        raise ResultCommitError("stage", destination, exc) from exc


def atomic_write_exact(path: Path, writer: Callable[[BinaryIO], None]) -> Path:
    """Publish one file atomically at exactly *path*."""

    destination = Path(path)
    staged: Path | None = None
    try:
        staged = stage_exact(destination, writer)
        try:
            _replace_exact(staged, destination)
        except Exception as exc:
            raise ResultCommitError("publish", destination, exc) from exc
        staged = None
        return destination
    finally:
        if staged is not None:
            try:
                staged.unlink(missing_ok=True)
            except OSError:
                pass


def commit_exact_bytes(
    primary: Path,
    payload: bytes,
    *,
    mirrors: Sequence[Path] = (),
) -> Path:
    """Stage every destination, then publish mirrors and primary last.

    Each directory entry is replaced atomically.  Filesystems cannot replace
    multiple names as one transaction: if primary publication fails after a
    mirror succeeds, both files remain valid but may represent different run
    IDs.  Current-run validation must reject that mixed generation.
    """

    primary_path = Path(primary)
    seen = {primary_path}
    mirror_paths: list[Path] = []
    for raw_path in mirrors:
        mirror = Path(raw_path)
        if mirror in seen:
            continue
        seen.add(mirror)
        mirror_paths.append(mirror)
    destinations = [*mirror_paths, primary_path]

    staged: dict[Path, Path] = {}
    try:
        for destination in destinations:
            staged[destination] = stage_exact(
                destination,
                lambda stream, content=payload: stream.write(content),
            )
        for destination in destinations:
            try:
                _replace_exact(staged[destination], destination)
            except Exception as exc:
                raise ResultCommitError("publish", destination, exc) from exc
            del staged[destination]
        return primary_path
    finally:
        for temporary in staged.values():
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                pass


def _replace_exact(staged: Path, destination: Path) -> None:
    """Publish one staged sibling; kept separate as a fault-injection seam."""

    os.replace(staged, destination)


def commit_payloads(
    primary: Path,
    payloads: Mapping[Path, bytes],
) -> Path:
    """Stage distinct payloads, then publish companions before *primary*.

    Every new payload and prior generation is staged before publication.  The
    primary is replaced last.  If a later publication fails, already-published
    destinations are restored atomically (or removed when newly created).
    """

    primary_path = Path(primary)
    normalized = {Path(path): bytes(payload) for path, payload in payloads.items()}
    if primary_path not in normalized:
        raise ValueError(f"primary destination {primary_path} is missing from payloads")
    destinations = [path for path in normalized if path != primary_path]
    destinations.append(primary_path)

    staged: dict[Path, Path] = {}
    prior: dict[Path, Path | None] = {}
    published: list[Path] = []
    try:
        for destination in destinations:
            staged[destination] = stage_exact(
                destination,
                lambda stream, content=normalized[destination]: stream.write(content),
            )
        for destination in destinations:
            if destination.exists():
                try:
                    prior_content = destination.read_bytes()
                    prior[destination] = stage_exact(
                        destination,
                        lambda stream, content=prior_content: stream.write(content),
                    )
                except Exception as exc:
                    if isinstance(exc, ResultCommitError):
                        raise
                    raise ResultCommitError("backup", destination, exc) from exc
            else:
                prior[destination] = None
        for destination in destinations:
            try:
                _replace_exact(staged[destination], destination)
            except Exception as exc:
                for published_destination in reversed(published):
                    backup = prior[published_destination]
                    try:
                        if backup is None:
                            published_destination.unlink(missing_ok=True)
                        else:
                            _replace_exact(backup, published_destination)
                            prior[published_destination] = None
                    except Exception as rollback_exc:
                        raise ResultCommitError(
                            "rollback", published_destination, rollback_exc
                        ) from exc
                raise ResultCommitError("replace", destination, exc) from exc
            published.append(destination)
        return primary_path
    finally:
        for temporary in (*staged.values(), *prior.values()):
            if temporary is None:
                continue
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                pass


def commit_json_exact(
    primary: Path,
    payload: Mapping[str, Any],
    *,
    mirrors: Sequence[Path] = (),
) -> Path:
    """Serialize once and atomically commit identical JSON bytes."""

    content = serialize_json_bytes(payload)
    return commit_exact_bytes(primary, content, mirrors=mirrors)
