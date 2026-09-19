"""Verified state for resuming ``all`` post-processing at one segment."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import re
import shutil
from typing import Any, Iterable, Mapping, Sequence

import click


SCHEMA_VERSION = 1
_SEGMENT_PATH = re.compile(r"(?:^|/)segments/seg_(\d+)(?:/|$)")
_REGENERATED_ROOT_OUTPUTS = {
    "summary.json", "summary.log", "irc_plot_all.png",
    "energy_diagram_MLIP_all.png", "energy_diagram_G_MLIP_all.png",
    "energy_diagram_DFT_all.png", "energy_diagram_G_DFT_plus_MLIP_all.png",
}


def file_sha256(path: Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parm7_content_sha256(path: Path) -> str:
    """Hash Amber topology content while ignoring its generated date stamp."""

    digest = sha256()
    with Path(path).open("rb") as handle:
        first = handle.readline()
        if not first.startswith(b"%VERSION"):
            digest.update(first)
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def preferred_endpoint_trajectories(path_dir: Path) -> list[tuple[int, Path]]:
    """Return one endpoint trajectory per segment, preferring refined paths."""

    selected: dict[int, Path] = {}
    for pattern in (
        "seg_*_refine_mep/final_geometries_trj.xyz",
        "seg_*_mep/final_geometries_trj.xyz",
    ):
        for artifact in sorted(Path(path_dir).glob(pattern)):
            match = re.match(r"seg_(\d+)", artifact.parent.name)
            if match is not None:
                selected.setdefault(int(match.group(1)), artifact)
    return sorted(selected.items())


def _fingerprint(path: Path) -> dict[str, Any]:
    resolved = Path(path).resolve()
    if not resolved.is_file():
        raise click.ClickException(f"Resume identity input is missing: {resolved}")
    return {"path": str(resolved), "size": resolved.stat().st_size, "sha256": file_sha256(resolved)}


def build_resume_identity(*, inputs: Sequence[Path], ref_pdb: Path | None,
                          pipeline_mode: str, path_request: Mapping[str, Any],
                          calculator_identity: Mapping[str, Any], charge: int,
                          spin: int) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "pipeline_mode": str(pipeline_mode),
        "inputs": [_fingerprint(path) for path in inputs],
        "ref_pdb": None if ref_pdb is None else _fingerprint(ref_pdb),
        "path_request": json.loads(json.dumps(path_request, sort_keys=True)),
        "calculator_identity": json.loads(json.dumps(calculator_identity, sort_keys=True)),
        "charge": int(charge), "spin": int(spin),
    }


def validate_resume_identity(recorded: Any, current: Mapping[str, Any], *,
                             segment: int, available_segments: Iterable[int]) -> None:
    if not isinstance(recorded, Mapping):
        raise click.UsageError(
            "This output predates verified segment resume. Run the MEP once with "
            "this version before using --resume-segment."
        )
    if recorded.get("schema_version") != SCHEMA_VERSION:
        raise click.UsageError("Unsupported segment-resume metadata version.")
    available = sorted({int(value) for value in available_segments if int(value) > 0})
    if segment not in available:
        rendered = ", ".join(str(value) for value in available) or "none"
        raise click.BadParameter(
            f"segment {segment} is not present in the saved MEP; available: {rendered}",
            param_hint="--resume-segment",
        )
    for key in ("pipeline_mode", "inputs", "ref_pdb", "path_request",
                "calculator_identity", "charge", "spin"):
        if recorded.get(key) != current.get(key):
            raise click.UsageError(
                f"--resume-segment cannot reuse the saved MEP because {key} changed. "
                "Repeat the original path/extraction/calculator settings and change only "
                "post-processing options."
            )


def load_previous_manifest(path: Path, *, out_dir: Path) -> dict[str, Any]:
    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise click.ClickException(
            f"Cannot read the saved run manifest required for resume: {path}: {exc}"
        ) from exc
    produced = payload.get("produced") if isinstance(payload, dict) else None
    if not isinstance(produced, dict):
        raise click.ClickException(f"Saved run manifest is invalid: {path}")
    root = Path(out_dir).resolve()
    for key, entry in produced.items():
        if not isinstance(entry, dict) or not isinstance(entry.get("stamp"), dict):
            raise click.ClickException(f"Saved artifact entry is invalid: {key}")
        artifact = Path(str(entry.get("path", ""))).resolve()
        if not artifact.is_relative_to(root):
            raise click.ClickException(f"Saved artifact escapes the output directory: {artifact}")
        if not entry["stamp"].get("sha256"):
            raise click.ClickException(f"Saved artifact lacks a verification digest: {artifact}")
    return payload


def retained_public_artifact(path: Path, *, out_dir: Path, resume_segment: int) -> bool:
    relative = Path(path).resolve().relative_to(Path(out_dir).resolve()).as_posix()
    if relative in _REGENERATED_ROOT_OUTPUTS:
        return False
    match = _SEGMENT_PATH.search(relative)
    if match is not None:
        return int(match.group(1)) < int(resume_segment)
    return True


def invalidate_from_segment(out_dir: Path, segment: int) -> list[Path]:
    root = Path(out_dir).resolve(); removed: list[Path] = []
    segments_root = root / "segments"
    if segments_root.is_dir():
        for candidate in segments_root.glob("seg_[0-9][0-9]*"):
            try: index = int(candidate.name.split("_", 1)[1])
            except (IndexError, ValueError): continue
            if index >= int(segment):
                if candidate.is_dir() and not candidate.is_symlink(): shutil.rmtree(candidate)
                else: candidate.unlink(missing_ok=True)
                removed.append(candidate)
    for name in _REGENERATED_ROOT_OUTPUTS:
        candidate = root / name
        if candidate.is_file() or candidate.is_symlink():
            candidate.unlink(missing_ok=True); removed.append(candidate)
    return removed


def retained_post_segments(summary: Mapping[str, Any], segment: int) -> list[dict[str, Any]]:
    retained = []
    for item in summary.get("post_segments", []) or []:
        if not isinstance(item, dict): continue
        try: index = int(item.get("index", 0) or 0)
        except (TypeError, ValueError): continue
        if 0 < index < int(segment): retained.append(dict(item))
    return retained


def retained_energy_diagrams(summary: Mapping[str, Any], segment: int) -> list[dict[str, Any]]:
    retained = []
    for item in summary.get("energy_diagrams", []) or []:
        if not isinstance(item, dict): continue
        image = str(item.get("image") or item.get("diagram") or "")
        name = Path(image).name
        if name in _REGENERATED_ROOT_OUTPUTS or name.endswith("_all.png"): continue
        match = _SEGMENT_PATH.search(image.replace("\\", "/"))
        if match is not None and int(match.group(1)) >= int(segment): continue
        retained.append(dict(item))
    return retained
