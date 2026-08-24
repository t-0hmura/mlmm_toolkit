#!/usr/bin/env python3
"""Engineering-marker coverage check.

Three static gates ported from the retired ``tests/test_fingerprint.py``:

* ``# CHEMISTRY-RULE:N`` coverage — every rule this repo is responsible for
  must be annotated somewhere in the package source.
* ``# DOMAIN_PURE`` coverage — the modules that must stay backend-agnostic
  carry this marker.
* External-library import scope — MLIP-only SDKs (``fairchem`` / ``orb_models``
  / ``mace`` / ``aimnet``) may only be imported under ``backends/`` so the
  rest of the package keeps a stable, lightweight import graph.

Run from the repository root::

    python .github/scripts/check_engineering_markers.py

Exit code 0 = all gates pass, 1 = at least one violation (printed to stderr).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

PKG = "mlmm"

REQUIRED_CHEMISTRY_RULES: set[int] = {1, 2, 3, 4, 5, 6, 7, 8, 9}

REQUIRED_DOMAIN_PURE: tuple[str, ...] = (
    "mlmm/backends/mlmm_calc.py",
    "mlmm/workflows/tsopt.py",
    "mlmm/workflows/freq.py",
)

EXTERNAL_LIBS_DENY: set[str] = {"fairchem", "orb_models", "mace", "aimnet"}

REPO_ROOT = Path(__file__).resolve().parents[2]
PKG_ROOT = REPO_ROOT / PKG


def _check_chemistry_rules() -> list[str]:
    found: set[int] = set()
    for py in PKG_ROOT.rglob("*.py"):
        for m in re.finditer(r"CHEMISTRY-RULE:(\d+)", py.read_text(errors="ignore")):
            found.add(int(m.group(1)))
    missing = sorted(REQUIRED_CHEMISTRY_RULES - found)
    return [f"missing CHEMISTRY-RULE markers in {PKG}: {missing}"] if missing else []


def _check_domain_pure() -> list[str]:
    missing = [
        rel for rel in REQUIRED_DOMAIN_PURE
        if "# DOMAIN_PURE" not in (REPO_ROOT / rel).read_text(errors="ignore")
    ]
    return [f"missing # DOMAIN_PURE marker in: {missing}"] if missing else []


def _check_external_library_scope() -> list[str]:
    backends_dir = PKG_ROOT / "backends"
    violations: list[str] = []
    for py in PKG_ROOT.rglob("*.py"):
        if backends_dir in py.parents or py == backends_dir:
            continue
        text = py.read_text(errors="ignore")
        for mod in EXTERNAL_LIBS_DENY:
            if re.search(
                rf"^(?:import\s+{re.escape(mod)}\b|from\s+{re.escape(mod)})",
                text, re.M,
            ):
                violations.append(f"{py.relative_to(REPO_ROOT)}: imports {mod}")
    return ["external library import outside backends/:", *violations] if violations else []


def main() -> int:
    failures: list[str] = []
    failures += _check_chemistry_rules()
    failures += _check_domain_pure()
    failures += _check_external_library_scope()
    if failures:
        for line in failures:
            print(line, file=sys.stderr)
        return 1
    print("engineering markers OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
