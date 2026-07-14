#!/usr/bin/env python3
"""Block high-impact semantic drift in handwritten docs and skills."""

from __future__ import annotations

import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCAN_ROOTS = (REPO_ROOT / "docs", REPO_ROOT / "skills")
LEGACY_BOOL_PAGES = {
    Path("docs/cli-conventions.md"),
    Path("docs/ja/cli-conventions.md"),
    Path("docs/concepts.md"),
    Path("docs/ja/concepts.md"),
    Path("docs/glossary.md"),
    Path("docs/ja/glossary.md"),
}
VALUE_BOOL_RE = re.compile(r"--[a-z][a-z0-9-]*\s+(?:True|False|yes|no|on|off)\b", re.IGNORECASE)
STALE_PATTERNS: tuple[tuple[re.Pattern[str], str], ...] = (
    (re.compile(
        r"uma-s-1p1\s*(?:\([^)]*\))?\s*(?:is\s+)?(?:the\s+)?default\b"
        r"|\bdefault(?:\s+(?:model|is))?\s*[:=]?\s*[`'\"]?uma-s-1p1\b",
        re.I,
    ),
     "UMA's default model is uma-s-1p2"),
    (re.compile(r"\b(?:Analytical|analytical|解析[^\n]{0,10}Hessian)[^\n]{0,40}(?:UMA-only|only[^\n]{0,12}UMA|UMA[^\n]{0,12}のみ)", re.I),
     "analytical Hessians are implemented for all bundled MLIP backends"),
    (re.compile(r"(?:auto[- ]?downgrad|automatically falls? back|自動的[^\n]{0,20}FiniteDifference)", re.I),
     "an explicit analytical Hessian request must raise when unavailable"),
    (re.compile(r"energy_diagram_(?:G_)?UMA|gibbs_(?:dft_)?uma|post_segments\[i\]\.uma\b"),
     "backend-neutral output names use MLIP/mlip"),
)


def _iter_markdown() -> list[Path]:
    paths: list[Path] = []
    for root in SCAN_ROOTS:
        paths.extend(root.rglob("*.md"))
    return sorted(path for path in paths if "reference" not in path.relative_to(REPO_ROOT).parts)


def main() -> int:
    sys.path.insert(0, str(REPO_ROOT))
    from mlmm.core.defaults import DEFAULT_UMA_MODEL, MLMM_CALC_KW

    errors: list[str] = []
    if DEFAULT_UMA_MODEL != "uma-s-1p2":
        errors.append(f"source default changed: DEFAULT_UMA_MODEL={DEFAULT_UMA_MODEL!r}")
    expected = {"uma_precision": "fp32", "orb_precision": "float64", "mace_dtype": "float64"}
    for key, value in expected.items():
        if MLMM_CALC_KW.get(key) != value:
            errors.append(f"source default changed: {key}={MLMM_CALC_KW.get(key)!r}, expected {value!r}")

    for path in _iter_markdown():
        rel = path.relative_to(REPO_ROOT)
        text = path.read_text(encoding="utf-8")
        if rel not in LEGACY_BOOL_PAGES:
            for match in VALUE_BOOL_RE.finditer(text):
                line = text.count("\n", 0, match.start()) + 1
                errors.append(f"{rel}:{line}: use canonical bare/toggle bool syntax: {match.group(0)!r}")
        for pattern, message in STALE_PATTERNS:
            for match in pattern.finditer(text):
                line = text.count("\n", 0, match.start()) + 1
                errors.append(f"{rel}:{line}: {message}: {match.group(0)!r}")

    if errors:
        print(f"[docs-contract] FAIL: {len(errors)} issue(s)")
        for error in errors:
            print(error)
        return 1
    print("[docs-contract] OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
