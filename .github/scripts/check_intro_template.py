#!/usr/bin/env python3
"""Validate the command-page heading template (EN/JA).

Every command (one per ``docs/reference/commands/*.md``) has an EN page and a
JA page. The EN page runs ``## What it is for`` →
``## Examples`` → ``## How it works`` → (a page-specific section on reading the
result) → ``## Output files`` → ``## Main options`` → optional ``## Notes`` →
``## See also``; the JA page mirrors it.
"""

from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_ROOT = REPO_ROOT / "docs"
COMMANDS_ROOT = DOCS_ROOT / "reference" / "commands"

LAYOUT = {
    "en": {
        "ordered": (
            "## What it is for",
            "## Examples",
            "## How it works",
            "## Output files",
            "## Main options",
            "## See also",
        ),
        "notes": "## Notes",
        "replaced": (
            "## Workflow",
            "## Outputs",
            "## CLI options",
            "## YAML configuration",
            "## Exit codes",
            "## See Also",
        ),
    },
    "ja": {
        "ordered": (
            "## 主な用途",
            "## 基本的な実行例",
            "## 処理の仕組みと計算仕様",
            "## 主な出力ファイル",
            "## 主な CLI オプション",
            "## 関連ドキュメント",
        ),
        "notes": "## 使用上の注意点",
        "replaced": (
            "## 実行例",
            "## 処理の流れ",
            "## 出力",
            "## CLI オプション",
            "## YAML 設定",
            "## 終了コード",
            "## 注記",
            "## 注意事項",
            "## 関連項目",
        ),
    },
}

FORBIDDEN = {
    "en": (
        "## When to use",
        "## Quick examples",
        "## Common examples",
        "## Minimal example",
        "## Usage",
        "## Inputs",
    ),
    "ja": (
        "## 使いどころ",
        "## クイック例",
        "## 使用例",
        "## 例",
        "## 最小例",
        "## 入力",
    ),
}


def _command_names() -> list[str]:
    return sorted(
        p.stem.replace("_", "-")
        for p in COMMANDS_ROOT.glob("*.md")
        if p.stem != "index"
    )


def _heading_positions(path: Path) -> dict[str, int]:
    positions: dict[str, int] = {}
    in_fence = False
    for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if line.lstrip().startswith(("```", "~~~")):
            in_fence = not in_fence
            continue
        if not in_fence and line.startswith(("## ", "### ")):
            positions.setdefault(line.strip(), i)
    return positions


def _check_new(path: Path, pos: dict[str, int], lay: dict, errors: list[str]) -> None:
    missing = [h for h in lay["ordered"] if h not in pos]
    if missing:
        errors.append(f"{path}: missing headings: {', '.join(missing)}")
    present = [h for h in lay["ordered"] if h in pos]
    for a, b in zip(present, present[1:]):
        if pos[a] > pos[b]:
            errors.append(f"{path}: '{a}' must precede '{b}'")
    notes = lay["notes"]
    if notes in pos:
        before, after = lay["ordered"][-2], lay["ordered"][-1]
        if (before in pos and pos[notes] < pos[before]) or (after in pos and pos[notes] > pos[after]):
            errors.append(f"{path}: '{notes}' must sit between '{before}' and '{after}'")
    mixed = [h for h in lay["replaced"] if h in pos]
    if mixed:
        errors.append(f"{path}: headings of the Examples/Workflow layout remain: {', '.join(mixed)}")


def _check(path: Path, lang: str, errors: list[str]) -> None:
    if not path.exists():
        errors.append(f"{path}: missing file")
        return
    pos = _heading_positions(path)
    lay = LAYOUT[lang]

    present_forbidden = [h for h in FORBIDDEN[lang] if h in pos]
    if present_forbidden:
        errors.append(f"{path}: headings must be removed: {', '.join(present_forbidden)}")

    _check_new(path, pos, lay, errors)


def main() -> int:
    errors: list[str] = []
    names = _command_names()
    if not names:
        errors.append(f"{COMMANDS_ROOT}: no generated command pages found")

    for name in names:
        _check(DOCS_ROOT / f"{name}.md", "en", errors)
        _check(DOCS_ROOT / "ja" / f"{name}.md", "ja", errors)

    if errors:
        print("[intro-check] failed:")
        for e in errors:
            print(f"- {e}")
        return 1

    print(f"[intro-check] validated {len(names) * 2} pages.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
