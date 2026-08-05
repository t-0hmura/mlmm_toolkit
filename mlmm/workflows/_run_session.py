"""Invocation-local ownership for the composite ``all`` workflow.

The public output layout intentionally remains unchanged.  This module only
tracks which exact files belong to the current invocation and which in-process
resources must be released when Click closes the command context.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import MutableMapping
from dataclasses import dataclass, field
import errno
import fcntl
from hashlib import sha256
import gc
import logging
import os
from pathlib import Path
import stat as stat_module
import time
from typing import Any, Callable, Iterable, Mapping, Sequence
from uuid import uuid4

import click

from mlmm.core.result_commit import (
    MLMM_RUN_ID_ENV,
    commit_json_exact,
    symlink_ancestor,
)
from mlmm.core.defaults import SEGMENTS_DIRNAME


logger = logging.getLogger(__name__)


def _lexical_absolute(path: Path) -> Path:
    """Return an absolute path without following any filesystem links."""

    return Path(os.path.abspath(os.fspath(path)))


@dataclass(frozen=True)
class ArtifactStamp:
    """Filesystem identity captured before or after a child dispatch."""

    exists: bool
    dev: int | None = None
    ino: int | None = None
    size: int | None = None
    mtime_ns: int | None = None
    ctime_ns: int | None = None
    sha256: str | None = None

    @classmethod
    def capture(cls, path: Path, *, digest: bool = False) -> "ArtifactStamp":
        candidate = _lexical_absolute(path)
        if symlink_ancestor(candidate) is not None:
            return cls(False)
        try:
            stat = candidate.lstat()
        except (FileNotFoundError, NotADirectoryError):
            return cls(False)
        # A claimed artifact must be a regular file at the declared lexical
        # path.  Following a symlink here would let an unrelated target satisfy
        # an exact-path declaration.
        if not stat_module.S_ISREG(stat.st_mode):
            return cls(False)
        checksum = None
        if digest:
            hasher = sha256()
            flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
            try:
                descriptor = os.open(candidate, flags)
            except (FileNotFoundError, NotADirectoryError):
                return cls(False)
            except OSError as exc:
                if exc.errno == errno.ELOOP:
                    return cls(False)
                raise
            with os.fdopen(descriptor, "rb") as handle:
                opened_stat = os.fstat(handle.fileno())
                if (
                    not stat_module.S_ISREG(opened_stat.st_mode)
                    or _stat_identity(opened_stat) != _stat_identity(stat)
                ):
                    return cls(False)
                for block in iter(lambda: handle.read(1024 * 1024), b""):
                    hasher.update(block)
            checksum = hasher.hexdigest()
            try:
                final_stat = candidate.lstat()
            except (FileNotFoundError, NotADirectoryError):
                return cls(False)
            if (
                not stat_module.S_ISREG(final_stat.st_mode)
                or _stat_identity(final_stat) != _stat_identity(opened_stat)
            ):
                return cls(False)
            stat = final_stat
        return cls(
            True,
            dev=int(stat.st_dev),
            ino=int(stat.st_ino),
            size=int(stat.st_size),
            mtime_ns=int(stat.st_mtime_ns),
            ctime_ns=int(stat.st_ctime_ns),
            sha256=checksum,
        )

    def identity(self) -> tuple[int | None, ...]:
        """Return the stat identity used for current-run discrimination."""

        return (self.dev, self.ino, self.size, self.mtime_ns, self.ctime_ns)

    def as_dict(self) -> dict[str, Any]:
        return {
            "exists": self.exists,
            "dev": self.dev,
            "ino": self.ino,
            "size": self.size,
            "mtime_ns": self.mtime_ns,
            "ctime_ns": self.ctime_ns,
            "sha256": self.sha256,
        }


def _stat_identity(value: os.stat_result) -> tuple[int, ...]:
    return (
        int(value.st_dev),
        int(value.st_ino),
        int(value.st_size),
        int(value.st_mtime_ns),
        int(value.st_ctime_ns),
    )


class ArtifactClaimError(click.ClickException):
    """Raised when a declared child artifact is not current-run owned."""


@dataclass
class InvocationManifest:
    """Declare child outputs before dispatch and claim only changed files."""

    run_id: str = field(
        default_factory=lambda: os.environ.get(MLMM_RUN_ID_ENV) or str(uuid4())
    )
    started_ns: int = field(default_factory=time.time_ns)
    expected: dict[str, tuple[Path, ...]] = field(default_factory=dict)
    baseline: dict[str, dict[str, ArtifactStamp]] = field(default_factory=dict)
    produced: dict[str, tuple[Path, ArtifactStamp]] = field(default_factory=dict)

    def declare(
        self,
        key: str,
        candidates: Sequence[Path],
        *,
        snapshot: Mapping[str, ArtifactStamp] | None = None,
    ) -> None:
        """Declare ordered exact candidates and capture their prior identity."""

        logical_key = str(key)
        if logical_key in self.expected:
            raise ValueError(f"artifact key {logical_key!r} was already declared")
        normalized: list[Path] = []
        for candidate in candidates:
            path = _lexical_absolute(candidate)
            if path not in normalized:
                normalized.append(path)
        if not normalized:
            raise ValueError(f"artifact key {logical_key!r} has no candidates")
        self.expected[logical_key] = tuple(normalized)
        prior: dict[str, ArtifactStamp] = {}
        for path in normalized:
            path_key = str(path)
            if snapshot is not None:
                prior[path_key] = snapshot.get(path_key, ArtifactStamp(False))
            else:
                prior[path_key] = ArtifactStamp.capture(path)
        self.baseline[logical_key] = prior

    def _current_candidates(self, key: str) -> list[tuple[Path, ArtifactStamp]]:
        if key not in self.expected:
            raise ArtifactClaimError(f"Artifact {key!r} was not declared before dispatch.")
        current: list[tuple[Path, ArtifactStamp]] = []
        for path in self.expected[key]:
            after_stat = ArtifactStamp.capture(path)
            if not after_stat.exists:
                continue
            claimed = self.produced.get(key)
            if (
                claimed is not None
                and claimed[0] == path
                and claimed[1].identity() == after_stat.identity()
            ):
                current.append(claimed)
                continue
            before = self.baseline[key][str(path)]
            if before.exists and before.identity() == after_stat.identity():
                continue
            after = ArtifactStamp.capture(path, digest=True)
            if not after.exists:
                continue
            if (not before.exists) or before.identity() != after.identity():
                current.append((path, after))
        return current

    def claim_one(self, key: str) -> Path:
        """Claim the first ordered candidate created or changed this run."""

        logical_key = str(key)
        current = self._current_candidates(logical_key)
        if not current:
            rendered = ", ".join(str(path) for path in self.expected[logical_key])
            raise ArtifactClaimError(
                f"Current invocation did not produce {logical_key!r}; expected one of: {rendered}"
            )
        path, stamp = current[0]
        self.produced[logical_key] = (path, stamp)
        return path

    def claim_optional(self, key: str) -> Path | None:
        """Claim a current candidate when present without admitting stale files."""

        logical_key = str(key)
        current = self._current_candidates(logical_key)
        if not current:
            self.produced.pop(logical_key, None)
            return None
        path, stamp = current[0]
        self.produced[logical_key] = (path, stamp)
        return path

    def require(self, keys: Iterable[str]) -> None:
        missing = [str(key) for key in keys if str(key) not in self.produced]
        if missing:
            raise ArtifactClaimError(
                "Current invocation has not claimed required artifacts: "
                + ", ".join(missing)
            )

    def path(self, key: str) -> Path:
        self.require([key])
        return self.produced[str(key)][0]

    def paths(self, prefix: str) -> list[Path]:
        """Enumerate claimed paths only, preserving declaration order."""

        wanted = str(prefix)
        return [
            self.produced[key][0]
            for key in self.expected
            if key.startswith(wanted) and key in self.produced
        ]

    def as_dict(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "started_ns": self.started_ns,
            "expected": {
                key: [str(path) for path in paths]
                for key, paths in self.expected.items()
            },
            "baseline": {
                key: {path: stamp.as_dict() for path, stamp in stamps.items()}
                for key, stamps in self.baseline.items()
            },
            "produced": {
                key: {"path": str(path), "stamp": stamp.as_dict()}
                for key, (path, stamp) in self.produced.items()
            },
        }

    def write_internal(self, path: Path) -> Path:
        """Atomically persist the private manifest using the commit primitive."""

        return commit_json_exact(Path(path), self.as_dict())


@dataclass
class InvocationResources:
    """Idempotent LIFO cleanup stack for one Click invocation."""

    callbacks: list[Callable[[], None]] = field(default_factory=list)
    closed: bool = False

    def add(self, callback: Callable[[], None]) -> Callable[[], None]:
        if self.closed:
            raise RuntimeError("invocation resource scope is already closed")
        self.callbacks.append(callback)
        return callback

    def own_environment(
        self,
        key: str,
        value: str,
        *,
        environ: MutableMapping[str, str] | None = None,
    ) -> str:
        """Bind one environment value and restore its exact prior presence."""

        target = os.environ if environ is None else environ
        name = str(key)
        owned_value = str(value)
        prior_present = name in target
        prior_value = target[name] if prior_present else ""

        def restore() -> None:
            if prior_present:
                target[name] = prior_value
            else:
                target.pop(name, None)

        # Own restoration before mutating the process environment.
        self.add(restore)
        target[name] = owned_value
        return owned_value

    def own_cleanup(self, owner: Any) -> Any:
        cleanup = getattr(owner, "cleanup", None)
        if not callable(cleanup):
            raise TypeError("owned object must expose cleanup()")
        self.add(cleanup)
        return owner

    def own_exclusive_lock(self, path: Path) -> Path:
        """Hold one non-blocking process lock until this resource scope closes."""

        lock_path = _lexical_absolute(path)
        if symlink_ancestor(lock_path) is not None:
            raise ArtifactClaimError(
                f"Run lock has a symlinked ancestor: {lock_path}"
            )
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        flags = os.O_RDWR | os.O_CREAT | getattr(os, "O_NOFOLLOW", 0)
        try:
            descriptor = os.open(lock_path, flags, 0o600)
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            if "descriptor" in locals():
                os.close(descriptor)
            if exc.errno in {errno.EACCES, errno.EAGAIN}:
                raise ArtifactClaimError(
                    f"Another ML/MM run is already using output directory "
                    f"{lock_path.parent.parent}."
                ) from exc
            raise

        def release() -> None:
            try:
                fcntl.flock(descriptor, fcntl.LOCK_UN)
            finally:
                os.close(descriptor)

        self.add(release)
        return lock_path

    def close(self) -> None:
        if self.closed:
            return
        self.closed = True
        while self.callbacks:
            callback = self.callbacks.pop()
            try:
                callback()
            except BaseException as exc:
                logger.warning("Invocation cleanup failed: %s", exc)


@dataclass
class RunSession:
    """Field-stable owner shared by the composite workflow's facilities."""

    manifest: InvocationManifest = field(default_factory=InvocationManifest)
    resources: InvocationResources = field(default_factory=InvocationResources)

    @property
    def run_id(self) -> str:
        return self.manifest.run_id

    def own_run_id_environment(
        self,
        key: str = MLMM_RUN_ID_ENV,
        *,
        environ: MutableMapping[str, str] | None = None,
    ) -> str:
        """Expose this invocation's one run ID to every in-process child."""

        target = os.environ if environ is None else environ
        active_run_id = target.get(key)
        if active_run_id:
            self.manifest.run_id = active_run_id
        return self.resources.own_environment(
            key,
            self.manifest.run_id,
            environ=target,
        )

    def close(self) -> None:
        self.resources.close()


@dataclass
class CalculatorLease:
    """One parent-side calculator attached to a bounded geometry set."""

    calculator: Any
    geometries: list[Any] = field(default_factory=list)
    released: bool = False

    def attach(self, geometry: Any) -> Any:
        if self.released:
            raise RuntimeError("calculator lease is already released")
        geometry.set_calculator(self.calculator)
        self.geometries.append(geometry)
        return geometry

    def release(self) -> None:
        if self.released:
            return
        self.released = True
        calculator = self.calculator
        for geometry in self.geometries:
            try:
                if getattr(geometry, "calculator", None) is calculator:
                    geometry.calculator = None
            except Exception as exc:
                logger.debug("Failed to detach leased calculator: %s", exc)
        self.geometries.clear()
        close = getattr(calculator, "close", None)
        if callable(close):
            try:
                close()
            except Exception as exc:
                logger.warning("Calculator close failed: %s", exc)
        self.calculator = None
        calculator = None
        gc.collect()
        try:
            import torch

            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        except Exception as exc:
            logger.debug("CUDA cache release unavailable: %s", exc)


def refresh_current_public_outputs(
    manifest: InvocationManifest,
    root: Path,
) -> list[Path]:
    """Revalidate producer-declared public outputs only.

    Discovery is deliberately forbidden here: a changed file is public only
    after its producer declared its exact destination.
    """

    root = _lexical_absolute(root)
    for key, paths in list(manifest.expected.items()):
        if not key.startswith("output.public."):
            continue
        if len(paths) != 1:
            raise ArtifactClaimError(
                f"Public artifact {key!r} must declare exactly one destination."
            )
        expected_key = public_output_key(root, paths[0])
        if key != expected_key:
            raise ArtifactClaimError(
                f"Public artifact key {key!r} does not match {paths[0]}."
            )
        manifest.claim_optional(key)
    return manifest.paths("output.public.")


def current_key_output_files(
    manifest: InvocationManifest,
    root: Path,
) -> dict[str, Any]:
    """Build key-output metadata from claimed current-run paths only."""

    current_relative = set(current_output_paths(manifest, root))
    descriptions = {
        "summary.log": "Human-readable results summary",
        "summary.json": "Machine-readable results summary",
        "mep_trj.xyz": "Full MEP trajectory",
        "mep.xyz": "MEP trajectory",
        "mep.pdb": "Full MEP as PDB",
        "mep.cif": "Full MEP with original mmCIF identifiers",
        "ml_region.pdb": "Reusable ML-region model",
        "energy_diagram_MEP.png": "MEP energy plot",
        "mep_plot.png": "MEP energy plot (trj2fig)",
        "irc_plot_all.png": "Aggregated IRC plot",
        "energy_diagram_MLIP_all.png": "Aggregated MLIP energy plot",
        "energy_diagram_G_MLIP_all.png": "Aggregated MLIP Gibbs energy plot",
        "energy_diagram_DFT_all.png": "Aggregated DFT energy plot",
        "energy_diagram_G_DFT_plus_MLIP_all.png": (
            "Aggregated DFT//MLIP/MM Gibbs energy plot"
        ),
    }
    key_files: dict[str, Any] = {}
    segment_files: dict[str, list[str]] = defaultdict(list)
    for relative in sorted(current_relative):
        parts = Path(relative).parts
        if (
            len(parts) >= 3
            and parts[0] == SEGMENTS_DIRNAME
            and parts[1].startswith("seg_")
        ):
            segment_files[parts[1]].append(Path(*parts[2:]).as_posix())
        elif len(parts) == 1:
            key_files[relative] = descriptions.get(
                relative, "Current-run pipeline output"
            )
    for segment, files in segment_files.items():
        key_files[segment] = {
            "description": f"Per-segment results for {segment}",
            "files": files,
        }
    return key_files


def current_output_paths(
    manifest: InvocationManifest,
    root: Path,
) -> list[str]:
    """Return sorted root-relative paths claimed by the current invocation."""

    root = _lexical_absolute(root)
    current_public = refresh_current_public_outputs(manifest, root)
    return sorted({
        path.relative_to(root).as_posix()
        for path in current_public
        if path.is_relative_to(root)
    })


def public_output_key(root: Path, path: Path) -> str:
    """Return the manifest key for one exact public-layout destination."""

    public_root = _lexical_absolute(root)
    destination = _lexical_absolute(path)
    linked = symlink_ancestor(destination)
    if linked is not None:
        raise ValueError(
            f"public output {destination} has symlinked ancestor {linked}"
        )
    try:
        relative = destination.relative_to(public_root)
    except ValueError as exc:
        raise ValueError(
            f"public output {destination} is outside pipeline root {public_root}"
        ) from exc
    parts = relative.parts
    if not parts or (len(parts) > 1 and parts[0] != SEGMENTS_DIRNAME):
        raise ValueError(
            "public outputs must be root files or descendants of "
            f"{SEGMENTS_DIRNAME!r}: {destination}"
        )
    return f"output.public.{relative.as_posix()}"


def declare_public_output(
    manifest: InvocationManifest,
    root: Path,
    path: Path,
) -> str:
    """Register one producer-owned public destination before writing it."""

    destination = _lexical_absolute(path)
    key = public_output_key(root, destination)
    if key in manifest.expected:
        if manifest.expected[key] != (destination,):
            raise ValueError(f"artifact key {key!r} has a conflicting destination")
        return key
    manifest.declare(key, [destination])
    return key


def claim_public_output(
    manifest: InvocationManifest,
    root: Path,
    path: Path,
) -> Path | None:
    """Claim a previously declared public output after its producer returns."""

    key = declare_public_output(
        manifest,
        root,
        path,
    )
    return manifest.claim_optional(key)
