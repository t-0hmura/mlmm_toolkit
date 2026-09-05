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
    (re.compile(r"torch_scatter|torch-scatter", re.I),
     "the current ORB extra does not require torch_scatter"),
    (re.compile(
        r"Load CUDA \*\*before\*\* activating conda"
        r"|wheel に合わせる \(cu126 ↔ 12\.6, cu129 ↔ 12\.9\)",
        re.I,
    ),
     "prebuilt PyTorch wheels do not require a matching local CUDA module"),
    (re.compile(r"result\.json\.electronic_state_verified"),
     "IRC electronic-state verification is nested under rigid_projection"),
    (re.compile(r"ref_order=(?:verified|unique-elements)"),
     "ONIOM import emits identity-verified or element-verified"),
    (re.compile(r"FD Hessians \(more VRAM\)", re.I),
     "Hessian peak memory depends on backend, system, precision, and hardware"),
    (re.compile(r"No torch / no MLIP dependency", re.I),
     "domain may use numeric torch/numpy but not MLIP runtime dependencies"),
    (re.compile(r"\bcalc\.(?:charge|spin)\b"),
     "calculator electronic-state keys are model_charge/model_mult"),
)

REQUIRED_SNIPPETS: dict[Path, tuple[str, ...]] = {
    Path("docs/json-output.md"): (
        "`mlip_precision`",
        "`energy_first_hartree`",
        "`never_stop_energy_bypasses`",
        "`safeguards`",
        "`hessian_npz`",
        "`rigid_projection.electronic_state_verified`",
        "`references`",
        "{method, citation, doi}",
        "The treatment is always\n`constrained`",
    ),
    Path("docs/ja/json-output.md"): (
        "`mlip_precision`",
        "`energy_first_hartree`",
        "`never_stop_energy_bypasses`",
        "`safeguards`",
        "`hessian_npz`",
        "`rigid_projection.electronic_state_verified`",
        "`references`",
        "{method, citation, doi}",
        "処理は常に `constrained`",
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
        'd["rigid_projection"]["electronic_state_verified"]',
    ),
    Path("skills/mlmm-workflows-output/SKILL.md"): (
        "`mlip_model`",
        "`mlip_precision`",
        "`references`",
        "{method, citation, doi}",
        "ml_region_{without,with}_linkH.xyz",
        "parm7 bonds crossing that selection",
    ),
    Path("docs/dft.md"): (
        "ml_region_without_linkH.xyz",
        "ml_region_with_linkH.xyz",
        "ml_region_without_linkH.pdb",
        "ml_region_with_linkH.pdb",
    ),
    Path("docs/ja/dft.md"): (
        "ml_region_without_linkH.xyz",
        "ml_region_with_linkH.xyz",
        "ml_region_without_linkH.pdb",
        "ml_region_with_linkH.pdb",
    ),
    Path("docs/glossary.md"): (
        "parm7 bond crossing the real-atom ML/MM selection",
        "inspection-only pocket caps",
    ),
    Path("docs/ja/glossary.md"): (
        "実在原子のML選択を横切るparm7結合",
        "抽出用リンク水素",
    ),
    Path("skills/mlmm-cli/dft.md"): (
        "ml_region_without_linkH.xyz",
        "ml_region_with_linkH.xyz",
        "ml_region_without_linkH.pdb",
        "ml_region_with_linkH.pdb",
    ),
    Path("skills/mlmm-cli/all.md"): (
        "the required child `thermoanalysis.yaml` handoff is retained even under `--no-dump`",
    ),
    Path("skills/mlmm-cli/tsopt.md"): (
        "`tsopt` always forces `reject_uphill=False`",
    ),
    Path("skills/mlmm-cli/opt.md"): (
        "final convergence check on the retained geometry",
        "convergence requires ALL of `max(|force|) <= 3e-4`",
        "deliberately tightened variant of the published",
    ),
    Path("docs/freq.md"): ("E + G_corr = G",),
    Path("docs/ja/freq.md"): ("E + G_corr = G",),
    Path("skills/mlmm-cli/freq.md"): ("E + G_corr = G",),
    Path("docs/backends.md"): ("mlmm all", "forwards the same factory"),
    Path("docs/ja/backends.md"): ("mlmm all", "同じfactory"),
    Path("docs/add-elem-info.md"): (
        "`--inplace/--no-inplace`",
        "`<input>_add_elem.pdb`",
        "Every input line is preserved byte-for-byte except columns 77–78",
    ),
    Path("docs/ja/add-elem-info.md"): (
        "`--inplace/--no-inplace`",
        "`<input>_add_elem.pdb`",
        "列 77–78 を除き、各入力行はそのまま\n保持されます",
    ),
    Path("skills/mlmm-install-backends/SKILL.md"): (
        "torch==2.13.0",
        "`cpu`, `cu126`, `cu130`, `cu132`",
    ),
    Path("skills/mlmm-install-backends/env-cuda.md"): (
        "torch==2.13.0",
        "`cu126`, `cu130`, `cu132`, and `cpu`",
        "does not require a\nmatching local CUDA toolkit",
    ),
    Path("skills/mlmm-install-backends/mace.md"): ("torch==2.13.0",),
    Path("docs/device-hpc.md"): (
        "g++ -std=c++20 -x c++ -fsyntax-only /dev/null",
        "command -v ninja",
    ),
    Path("docs/ja/device-hpc.md"): (
        "g++ -std=c++20 -x c++ -fsyntax-only /dev/null",
        "command -v ninja",
    ),
    Path("docs/oniom-import.md"): (
        "`ref_order=identity-verified`, `element-verified`, or `unverified-opt-in`",
    ),
    Path("docs/ja/oniom-import.md"): (
        "`ref_order=identity-verified` / `element-verified` / `unverified-opt-in`",
    ),
    Path("docs/irc.md"): (
        "Net charge of the ML region/model system",
        '`result.json["rigid_projection"]["electronic_state_verified"]`',
    ),
    Path("docs/ja/irc.md"): (
        "ML 領域/model system の正味電荷",
        '`result.json["rigid_projection"]["electronic_state_verified"]`',
    ),
    Path("skills/mlmm-hpc/SKILL.md"): (
        "g++ -std=c++20 -x c++ -fsyntax-only /dev/null",
        "command -v ninja",
    ),
    Path("skills/mlmm-hpc/dynamic-dispatch.md"): (
        "g++ -std=c++20 -x c++ -fsyntax-only /dev/null",
        "command -v ninja",
    ),
    Path("skills/mlmm-cli/extract.md"): (
        "`mm-parm → extract →\n  define-layer",
        "`mlmm mm-parm -i input.pdb --out-prefix system`",
        "same atom identity/order as `system.parm7`",
    ),
}

ORDERED_SNIPPETS: dict[Path, tuple[str, ...]] = {
    Path("skills/mlmm-install-backends/SKILL.md"): (
        "      - mlmm-toolkit  # install core first",
        "pip uninstall -y fairchem-core",
        "pip install mace-torch",
    ),
    Path("skills/mlmm-install-backends/mace.md"): (
        "pip install mlmm-toolkit",
        "pip uninstall -y fairchem-core",
        "pip install mace-torch",
    ),
}

YAML_ALL_SECTIONS = frozenset({
    "geom", "calc", "opt", "lbfgs", "rfo", "gs", "dmf", "irc", "freq",
    "thermo", "dft", "bias", "bond", "search", "hessian_dimer", "rsirfo",
    "stopt", "microiter",
})
YAML_ALL_SENTENCES = {
    Path("docs/yaml-reference.md"):
        "`mlmm all` consumes only the sections for its selected active stages.",
    Path("docs/ja/yaml-reference.md"):
        "`mlmm all` は、選択して有効化した stage の section だけを使用します。",
}


def _iter_markdown() -> list[Path]:
    paths: list[Path] = []
    for root in SCAN_ROOTS:
        paths.extend(root.rglob("*.md"))
    return sorted(path for path in paths if "reference" not in path.relative_to(REPO_ROOT).parts)


def main() -> int:
    sys.path.insert(0, str(REPO_ROOT))
    from mlmm.core.defaults import DEFAULT_UMA_MODEL, LBFGS_KW, MLMM_CALC_KW, RFO_KW

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

    uphill_defaults = {LBFGS_KW["uphill_tolerance"], RFO_KW["uphill_tolerance"]}
    if len(uphill_defaults) != 1:
        errors.append(f"source defaults disagree: uphill_tolerance={uphill_defaults!r}")
    uphill_default = next(iter(uphill_defaults))
    for rel in (Path("docs/yaml-reference.md"), Path("docs/ja/yaml-reference.md")):
        text = (REPO_ROOT / rel).read_text(encoding="utf-8")
        expected_literal = f"uphill_tolerance: {uphill_default:.4f}"
        if text.count(expected_literal) != 2:
            errors.append(f"{rel}: uphill_tolerance examples must match {uphill_default}")
        if YAML_ALL_SENTENCES[rel] not in text:
            errors.append(f"{rel}: conditional all-stage routing sentence is missing")
        all_sections = set()
        for line in text.splitlines():
            match = re.match(r"^\| \[`([^`]+)`\]\([^)]*\) \|.*\| ([^|]+) \|$", line)
            if match is None:
                continue
            used_by = {item.strip() for item in match.group(2).split(",")}
            if "all" in used_by:
                all_sections.add(match.group(1))
        if all_sections != YAML_ALL_SECTIONS:
            errors.append(
                f"{rel}: all-section matrix={sorted(all_sections)!r}; "
                f"expected {sorted(YAML_ALL_SECTIONS)!r}"
            )

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

    for rel, snippets in ORDERED_SNIPPETS.items():
        text = (REPO_ROOT / rel).read_text(encoding="utf-8")
        positions = [text.find(snippet) for snippet in snippets]
        if any(position < 0 for position in positions) or positions != sorted(positions):
            errors.append(
                f"{rel}: MACE environment order must be mlmm install, "
                "fairchem-core removal, then MACE install"
            )

    if errors:
        print(f"[docs-contract] FAIL: {len(errors)} issue(s)")
        for error in errors:
            print(error)
        return 1
    print("[docs-contract] OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
