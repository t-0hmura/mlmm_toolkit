#!/usr/bin/env python3
"""Validate local markdown links under docs/."""

from __future__ import annotations

import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPO_ROOT / "docs"
LINK_RE = re.compile(r"(?<!\!)\[[^\]]+\]\(([^)]+)\)")
EXTERNAL_PREFIXES = ("http://", "https://", "mailto:", "tel:")


def _normalize_target(raw: str) -> str:
    target = raw.strip()
    if target.startswith("<") and target.endswith(">"):
        target = target[1:-1].strip()
    if " " in target:
        target = target.split(" ", 1)[0]
    return target


def _iter_links(path: Path) -> list[tuple[int, str]]:
    links: list[tuple[int, str]] = []
    in_fence = False
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        stripped = line.strip()
        if stripped.startswith("```"):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        for match in LINK_RE.finditer(line):
            links.append((lineno, _normalize_target(match.group(1))))
    return links


def _iter_toctree_targets(path: Path) -> list[tuple[int, str]]:
    targets: list[tuple[int, str]] = []
    in_toctree = False
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        stripped = line.strip()
        if stripped.startswith("```"):
            marker = stripped[3:].strip().lower()
            if not in_toctree and marker.startswith("{toctree}"):
                in_toctree = True
                continue
            if in_toctree:
                in_toctree = False
            continue
        if not in_toctree:
            continue
        if not stripped or stripped.startswith(":"):
            continue
        targets.append((lineno, stripped.split()[0]))
    return targets


def _resolve_doc_target(path: Path, target: str) -> Path | None:
    if not target or target.startswith(EXTERNAL_PREFIXES) or target.startswith("#"):
        return None
    if "*" in target:
        return None
    base = (path.parent / target).resolve()
    if base.suffix:
        return base
    for candidate in (base.with_suffix(".md"), base.with_suffix(".rst"), base):
        if candidate.exists():
            return candidate
    return base.with_suffix(".md")


def _check_path(path: Path, errors: list[str]) -> None:
    for lineno, target in _iter_links(path):
        if not target or target.startswith(EXTERNAL_PREFIXES):
            continue
        if target.startswith("#"):
            continue
        target_path = target.split("#", 1)[0]
        if not target_path:
            continue
        resolved = (path.parent / target_path).resolve()
        if not resolved.exists():
            errors.append(
                f"{path.relative_to(REPO_ROOT)}:{lineno}: broken local link -> {target}"
            )

    for lineno, target in _iter_toctree_targets(path):
        resolved = _resolve_doc_target(path, target)
        if resolved is None:
            continue
        if not resolved.exists():
            errors.append(
                f"{path.relative_to(REPO_ROOT)}:{lineno}: broken toctree target -> {target}"
            )


def main() -> int:
    errors: list[str] = []
    pages = sorted(DOCS_ROOT.rglob("*.md"))
    for md in pages:
        _check_path(md, errors)

    if errors:
        print("[link-check] failed:")
        for e in errors:
            print(f"- {e}")
        return 1

    print(f"[link-check] validated local links in {len(pages)} pages.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
