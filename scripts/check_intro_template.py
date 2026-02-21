#!/usr/bin/env python3
"""Validate intro-template headings for key docs pages (EN/JA)."""

from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPO_ROOT / "docs"

TARGETS = (
    "all",
    "scan",
    "tsopt",
    "freq",
    "opt",
    "path_opt",
    "path_search",
    "irc",
    "dft",
)

REQUIRED_EN = (
    "## Minimal example",
    "## Output checklist",
    "## Common examples",
)
REQUIRED_JA = (
    "## 最小例",
    "## 出力の見方",
    "## よくある例",
)


def _heading_positions(path: Path) -> dict[str, int]:
    positions: dict[str, int] = {}
    for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if line.startswith("## "):
            positions.setdefault(line.strip(), i)
    return positions


def _check_one(path: Path, required: tuple[str, ...], usage_heading: str, errors: list[str]) -> None:
    if not path.exists():
        errors.append(f"{path}: missing file")
        return

    positions = _heading_positions(path)
    missing = [h for h in required if h not in positions]
    if missing:
        errors.append(f"{path}: missing headings: {', '.join(missing)}")
        return

    req_lines = [positions[h] for h in required]
    if req_lines != sorted(req_lines):
        pairs = ", ".join(f"{h}@{positions[h]}" for h in required)
        errors.append(f"{path}: intro heading order mismatch ({pairs})")

    usage_line = positions.get(usage_heading)
    if usage_line is not None:
        late = [h for h in required if positions[h] > usage_line]
        if late:
            errors.append(
                f"{path}: intro headings must be before '{usage_heading}' "
                f"(late: {', '.join(late)})"
            )


def main() -> int:
    errors: list[str] = []

    for name in TARGETS:
        _check_one(
            DOCS_ROOT / f"{name}.md",
            REQUIRED_EN,
            "## Usage",
            errors,
        )
        _check_one(
            DOCS_ROOT / "ja" / f"{name}.md",
            REQUIRED_JA,
            "## 使用法",
            errors,
        )

    if errors:
        print("[intro-check] failed:")
        for e in errors:
            print(f"- {e}")
        return 1

    print(f"[intro-check] validated {len(TARGETS) * 2} pages.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
