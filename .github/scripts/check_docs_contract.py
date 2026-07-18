#!/usr/bin/env python3
"""Block high-impact semantic drift in handwritten docs and skills."""

from __future__ import annotations

import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCAN_ROOTS = (REPO_ROOT / "docs", REPO_ROOT / "skills")
sys.path.insert(0, str(Path(__file__).resolve().parent))
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

REQUIRED_SNIPPETS: dict[Path, tuple[str, ...]] = {
    Path("docs/json-output.md"): (
        "`mlip_precision`",
        "`energy_first_hartree`",
        "`never_stop_energy_bypasses`",
        "`safeguards`",
        "`hessian_npz`",
    ),
    Path("docs/ja/json-output.md"): (
        "`mlip_precision`",
        "`energy_first_hartree`",
        "`never_stop_energy_bypasses`",
        "`safeguards`",
        "`hessian_npz`",
    ),
    Path("skills/mlmm-cli/sp.md"): (
        "`mlip_backend`",
        "`mlip_model`",
        "`mlip_precision`",
    ),
    Path("skills/mlmm-cli/irc.md"): (
        "`energy_first_hartree`",
        "`never_stop_energy_bypasses`",
        "active-DOF basis",
        "inserts one underscore",
    ),
    Path("skills/mlmm-workflows-output/SKILL.md"): (
        "`mlip_model`",
        "`mlip_precision`",
    ),
    Path("docs/backends.md"): ("mlmm all", "forwards the same factory"),
    Path("docs/ja/backends.md"): ("mlmm all", "同じfactory"),
}


def _iter_markdown() -> list[Path]:
    paths: list[Path] = []
    for root in SCAN_ROOTS:
        paths.extend(root.rglob("*.md"))
    return sorted(path for path in paths if "reference" not in path.relative_to(REPO_ROOT).parts)


def main() -> int:
    sys.path.insert(0, str(REPO_ROOT))
    from mlmm.core.defaults import DEFAULT_UMA_MODEL, MLMM_CALC_KW

    from docs_command_contract import (
        bool_style_sources,
        load_root_cli,
        resolve_live_bool_options,
        validate_bool_style,
    )

    errors: list[str] = []
    if DEFAULT_UMA_MODEL != "uma-s-1p2":
        errors.append(f"source default changed: DEFAULT_UMA_MODEL={DEFAULT_UMA_MODEL!r}")
    expected = {"uma_precision": "fp32", "orb_precision": "float64", "mace_dtype": "float64"}
    for key, value in expected.items():
        if MLMM_CALC_KW.get(key) != value:
            errors.append(f"source default changed: {key}={MLMM_CALC_KW.get(key)!r}, expected {value!r}")

    # Live-derived canonical bool-style check over the full authored surface
    # (README, CONTRIBUTING, docs, skills, examples, smoke scripts). Only names
    # that resolve to a live boolean option are flagged; bool_compat runtime
    # `--flag True/False` is unaffected.
    live_bool = resolve_live_bool_options(load_root_cli())
    errors.extend(validate_bool_style(bool_style_sources(), live_bool))

    for path in _iter_markdown():
        rel = path.relative_to(REPO_ROOT)
        text = path.read_text(encoding="utf-8")
        for pattern, message in STALE_PATTERNS:
            for match in pattern.finditer(text):
                line = text.count("\n", 0, match.start()) + 1
                errors.append(f"{rel}:{line}: {message}: {match.group(0)!r}")

    for rel, snippets in REQUIRED_SNIPPETS.items():
        path = REPO_ROOT / rel
        text = path.read_text(encoding="utf-8")
        for snippet in snippets:
            if snippet not in text:
                errors.append(f"{rel}: required semantic contract missing: {snippet!r}")

    if errors:
        print(f"[docs-contract] FAIL: {len(errors)} issue(s)")
        for error in errors:
            print(error)
        return 1
    print("[docs-contract] OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
