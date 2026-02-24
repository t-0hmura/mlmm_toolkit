#!/usr/bin/env python3
"""Smoke-test trajectory dump behavior for mlmm opt/tsopt commands."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
CLI_MODULE = "mlmm_toolkit"
TIMEOUT_ENV = "MLMM_DUMP_CASE_TIMEOUT_SEC"
GENERIC_TIMEOUT_ENV = "DOCS_DUMP_CASE_TIMEOUT_SEC"
DEFAULT_CASE_TIMEOUT_SEC = 900.0


@dataclass(frozen=True)
class Fixture:
    name: str
    input_pdb: Path
    real_parm7: Path
    charge: int
    multiplicity: int


@dataclass(frozen=True)
class Case:
    name: str
    args: list[str]
    expect_present: tuple[str, ...]
    expect_absent: tuple[str, ...] = ()


def _frame_count(xyz_path: Path) -> int:
    lines = xyz_path.read_text(encoding="utf-8").splitlines()
    i = 0
    n = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        nat = int(lines[i].strip())
        i += nat + 2
        n += 1
    return n


def _run_cli(args: list[str], timeout_sec: float | None = None) -> None:
    cmd = [sys.executable, "-m", CLI_MODULE, *args]
    try:
        proc = subprocess.run(
            cmd,
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            timeout=timeout_sec,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"Command timed out after {timeout_sec} sec: {' '.join(cmd)}"
        ) from exc
    if proc.returncode != 0:
        tail = (proc.stdout + "\n" + proc.stderr)[-5000:]
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{tail}")


def _check_runnable() -> tuple[bool, str]:
    probe = subprocess.run(
        [sys.executable, "-m", CLI_MODULE, "opt", "--help"],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
    )
    probe_out = f"{probe.stdout}\n{probe.stderr}"
    if "Command 'opt' is unavailable" in probe_out or "Missing dependency:" in probe_out:
        return False, "opt command is unavailable in this environment."
    if probe.returncode != 0:
        return False, f"probe command failed with exit code {probe.returncode}."
    return True, "ok"


def _fixture_from_dir(base: Path) -> Fixture | None:
    layered_pdb = base / "p_complex_layered.pdb"
    layered_parm7 = base / "p_complex.parm7"
    if layered_pdb.exists() and layered_parm7.exists():
        return Fixture(
            name=f"{base}/p_complex_layered",
            input_pdb=layered_pdb.resolve(),
            real_parm7=layered_parm7.resolve(),
            charge=-1,
            multiplicity=1,
        )

    complex_pdb = base / "complex.pdb"
    complex_parm7 = base / "complex.parm7"
    if complex_pdb.exists() and complex_parm7.exists():
        return Fixture(
            name=f"{base}/complex",
            input_pdb=complex_pdb.resolve(),
            real_parm7=complex_parm7.resolve(),
            charge=0,
            multiplicity=1,
        )
    return None


def _resolve_fixture() -> Fixture | None:
    candidates: list[Path] = []
    env_dir = os.environ.get("MLMM_DUMP_FIXTURE_DIR")
    if env_dir:
        candidates.append(Path(env_dir).expanduser())
    # Prefer repository-contained fixture to avoid machine-specific paths.
    candidates.append(REPO_ROOT / "hessian_ff" / "tests" / "data" / "small")

    for base in candidates:
        fixture = _fixture_from_dir(base)
        if fixture is not None:
            return fixture
    return None


def _validate_case(case: Case, base_dir: Path, timeout_sec: float | None = None) -> None:
    out_dir = base_dir / case.name
    args = [*case.args, "--out-dir", str(out_dir)]
    _run_cli(args, timeout_sec=timeout_sec)

    for rel in case.expect_present:
        p = out_dir / rel
        if not p.exists():
            raise RuntimeError(f"[{case.name}] expected '{rel}' to exist.")
        if p.suffix == ".xyz":
            n_frames = _frame_count(p)
            if n_frames <= 0:
                raise RuntimeError(f"[{case.name}] '{rel}' has no frames.")

    for rel in case.expect_absent:
        p = out_dir / rel
        if p.exists():
            raise RuntimeError(f"[{case.name}] expected '{rel}' to be absent.")


def main() -> int:
    runnable, reason = _check_runnable()
    if not runnable:
        print(f"[dump-smoke] skipped: {reason}")
        return 0

    fixture = _resolve_fixture()
    if fixture is None:
        print(
            "[dump-smoke] skipped: no fixture found. "
            "Set MLMM_DUMP_FIXTURE_DIR or provide hessian_ff/tests/data/small inputs."
        )
        return 0

    common = [
        "-i",
        str(fixture.input_pdb),
        "--real-parm7",
        str(fixture.real_parm7),
        "-q",
        str(fixture.charge),
        "-m",
        str(fixture.multiplicity),
        "--max-cycles",
        "1",
    ]

    cases = [
        Case(
            name="opt_light_dump",
            args=["opt", *common, "--opt-mode", "light", "--dump"],
            expect_present=("optimization_trj.xyz", "optimization_all_trj.xyz"),
        ),
        Case(
            name="opt_heavy_dump",
            args=["opt", *common, "--opt-mode", "heavy", "--dump"],
            expect_present=("optimization_all_trj.xyz",),
        ),
        Case(
            name="tsopt_light_dump",
            args=["tsopt", *common, "--opt-mode", "light", "--dump"],
            expect_present=("optimization_all_trj.xyz",),
        ),
        Case(
            name="tsopt_heavy_dump",
            args=["tsopt", *common, "--opt-mode", "heavy", "--dump"],
            expect_present=("optimization_all_trj.xyz",),
        ),
        Case(
            name="tsopt_heavy_nodump",
            args=["tsopt", *common, "--opt-mode", "heavy", "--no-dump"],
            expect_present=(),
            expect_absent=("optimization_trj.xyz", "optimization_all_trj.xyz"),
        ),
        Case(
            name="opt_light_dump_legacy_bool",
            args=["opt", *common, "--opt-mode", "light", "--dump", "True"],
            expect_present=("optimization_trj.xyz", "optimization_all_trj.xyz"),
        ),
    ]

    timeout_raw = ""
    timeout_env = TIMEOUT_ENV
    for env_name in (GENERIC_TIMEOUT_ENV, TIMEOUT_ENV):
        raw = os.environ.get(env_name, "").strip()
        if raw:
            timeout_raw = raw
            timeout_env = env_name
            break
    timeout_sec = float(timeout_raw) if timeout_raw else DEFAULT_CASE_TIMEOUT_SEC
    if timeout_sec <= 0:
        timeout_sec = None

    if timeout_sec is None:
        print("[dump-smoke] per-case timeout: disabled")
    else:
        print(f"[dump-smoke] per-case timeout: {timeout_sec:.1f}s (env: {timeout_env})")

    with tempfile.TemporaryDirectory(prefix="mlmm_dump_smoke_") as td:
        base_dir = Path(td)
        for idx, case in enumerate(cases, start=1):
            print(f"[dump-smoke] case {idx}/{len(cases)}: {case.name}")
            started = time.perf_counter()
            _validate_case(case, base_dir, timeout_sec=timeout_sec)
            elapsed = time.perf_counter() - started
            print(f"[dump-smoke] case ok: {case.name} ({elapsed:.1f}s)")

    print(f"[dump-smoke] validated {len(cases)} cases with fixture '{fixture.name}'.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
