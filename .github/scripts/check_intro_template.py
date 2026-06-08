#!/usr/bin/env python3
"""Validate intro-template headings for key docs pages (EN/JA).

Enforces the command-page template heading set: ``When to use`` /
``Quick examples`` / ``Inputs`` (EN) and ``使いどころ`` / ``実行例`` / ``入力``
(JA), in order, before the body heading.
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

# (required headings in order, body heading the intro must precede)
SETS_EN = (
    (("## When to use", "## Quick examples", "## Inputs"), "## Workflow"),
)
SETS_JA = (
    (("## 使いどころ", "## 実行例", "## 入力"), "## 処理の流れ"),
)


def _heading_positions(path: Path) -> dict[str, int]:
    positions: dict[str, int] = {}
    for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if line.startswith("## "):
            positions.setdefault(line.strip(), i)
    return positions


def _match_set(
    positions: dict[str, int], required: tuple[str, ...], body_heading: str
) -> list[str]:
    """Return [] if this set fully matches (present, ordered, before body), else issues."""
    missing = [h for h in required if h not in positions]
    if missing:
        return [f"missing headings: {', '.join(missing)}"]

    issues: list[str] = []
    req_lines = [positions[h] for h in required]
    if req_lines != sorted(req_lines):
        pairs = ", ".join(f"{h}@{positions[h]}" for h in required)
        issues.append(f"intro heading order mismatch ({pairs})")

    body_line = positions.get(body_heading)
    if body_line is not None:
        late = [h for h in required if positions[h] > body_line]
        if late:
            issues.append(
                f"intro headings must be before '{body_heading}' (late: {', '.join(late)})"
            )
    return issues


def _check_one(path: Path, sets: tuple, errors: list[str]) -> None:
    if not path.exists():
        errors.append(f"{path}: missing file")
        return

    positions = _heading_positions(path)
    # Pick the set whose required headings are all present, then validate it.
    for required, body_heading in sets:
        if all(h in positions for h in required):
            issues = _match_set(positions, required, body_heading)
            if issues:
                errors.append(f"{path}: " + "; ".join(issues))
            return

    # Neither set present: report against the new template (migration target).
    new_required = sets[0][0]
    missing = [h for h in new_required if h not in positions]
    errors.append(f"{path}: missing headings: {', '.join(missing)}")


def main() -> int:
    errors: list[str] = []

    for name in TARGETS:
        _check_one(DOCS_ROOT / f"{name}.md", SETS_EN, errors)
        _check_one(DOCS_ROOT / "ja" / f"{name}.md", SETS_JA, errors)

    if errors:
        print("[intro-check] failed:")
        for e in errors:
            print(f"- {e}")
        return 1

    print(f"[intro-check] validated {len(TARGETS) * 2} pages.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
