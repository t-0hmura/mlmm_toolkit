#!/usr/bin/env python3
"""Validate the command-page heading template (EN/JA).

EN canonical template: a one-paragraph intro (the "when to use" guidance is
folded into the opening sentence — no separate heading), then ``## Examples``
before the ``## Workflow`` body, and none of the legacy ``When to use`` /
``Quick examples`` / ``Inputs`` headings.

JA pages still follow the legacy template (migrated in a separate pass).
"""

from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_ROOT = REPO_ROOT / "docs"

TARGETS = (
    "all",
    "scan",
    "scan2d",
    "scan3d",
    "tsopt",
    "freq",
    "opt",
    "path-opt",
    "path-search",
    "irc",
    "dft",
)

# EN canonical template.
EN_REQUIRED = ("## Examples",)
EN_BODY = "## Workflow"
EN_FORBIDDEN = (
    "## When to use",
    "## Quick examples",
    "## Common examples",
    "## Minimal example",
    "## Usage",
    "## Inputs",
)

# JA legacy template (kept until the JA pages are migrated).
JA_REQUIRED = ("## 使いどころ", "## 実行例", "## 入力")
JA_BODY = "## 処理の流れ"


def _heading_positions(path: Path) -> dict[str, int]:
    positions: dict[str, int] = {}
    for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if line.startswith("## "):
            positions.setdefault(line.strip(), i)
    return positions


def _check(
    path: Path,
    required: tuple[str, ...],
    body_heading: str,
    forbidden: tuple[str, ...],
    errors: list[str],
) -> None:
    if not path.exists():
        errors.append(f"{path}: missing file")
        return

    positions = _heading_positions(path)

    missing = [h for h in required if h not in positions]
    if missing:
        errors.append(f"{path}: missing headings: {', '.join(missing)}")

    present_forbidden = [h for h in forbidden if h in positions]
    if present_forbidden:
        errors.append(
            f"{path}: legacy headings must be removed: {', '.join(present_forbidden)}"
        )

    body_line = positions.get(body_heading)
    if body_line is not None:
        late = [h for h in required if h in positions and positions[h] > body_line]
        if late:
            errors.append(
                f"{path}: heading(s) must precede '{body_heading}': {', '.join(late)}"
            )


def main() -> int:
    errors: list[str] = []

    for name in TARGETS:
        _check(DOCS_ROOT / f"{name}.md", EN_REQUIRED, EN_BODY, EN_FORBIDDEN, errors)
        _check(DOCS_ROOT / "ja" / f"{name}.md", JA_REQUIRED, JA_BODY, (), errors)

    if errors:
        print("[intro-check] failed:")
        for e in errors:
            print(f"- {e}")
        return 1

    print(f"[intro-check] validated {len(TARGETS) * 2} pages.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
