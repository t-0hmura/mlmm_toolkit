#!/usr/bin/env python3
"""Run docs/reference quality checks in one place for CI and local use."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
GENERIC_TIMEOUT_ENV = "DOCS_DUMP_CASE_TIMEOUT_SEC"


def _run_step(label: str, cmd: Sequence[str], env: dict[str, str]) -> float:
    print(f"[docs-quality] {label}", flush=True)
    print(f"[docs-quality] $ {' '.join(cmd)}", flush=True)
    started = time.monotonic()
    subprocess.run(cmd, cwd=REPO_ROOT, env=env, check=True)
    elapsed = time.monotonic() - started
    print(f"[docs-quality] {label}: done in {elapsed:.2f}s", flush=True)
    return elapsed


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run docs/reference/smoke checks used by CI."
    )
    parser.add_argument(
        "--skip-dump",
        action="store_true",
        help="Skip smoke_dump_trajectories.py (faster local verification).",
    )
    parser.add_argument(
        "--dump-timeout-sec",
        type=float,
        default=None,
        help=(
            "Per-case timeout for smoke_dump_trajectories.py. "
            f"Sets {GENERIC_TIMEOUT_ENV} for this run."
        ),
    )
    args = parser.parse_args()

    env = os.environ.copy()
    if args.dump_timeout_sec is not None:
        env[GENERIC_TIMEOUT_ENV] = str(args.dump_timeout_sec)

    steps: list[tuple[str, list[str]]] = [
        ("Regenerate references and check for diff", [sys.executable, "scripts/generate_reference.py"]),
        ("Verify no diff in generated references", ["git", "diff", "--exit-code", "docs/reference/"]),
        ("Check intro template headings", [sys.executable, "scripts/check_intro_template.py"]),
        ("Check markdown local links", [sys.executable, "scripts/check_markdown_links.py"]),
        ("Check all->scan option contract", [sys.executable, "scripts/check_all_scan_contract.py"]),
        ("Smoke docs commands", [sys.executable, "scripts/smoke_docs_commands.py"]),
    ]
    if not args.skip_dump:
        steps.append(
            ("Smoke dump trajectories", [sys.executable, "scripts/smoke_dump_trajectories.py"])
        )

    total_elapsed = 0.0
    for label, cmd in steps:
        total_elapsed += _run_step(label, cmd, env)

    print(f"[docs-quality] All checks passed. total={total_elapsed:.2f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
