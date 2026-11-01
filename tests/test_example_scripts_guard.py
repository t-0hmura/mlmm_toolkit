"""Behavior of the example run scripts, which run without ``set -e``.

A stub ``mlmm`` on PATH records the subcommands it is called with. Giving a
script an output directory that already exists must stop it before any ``mlmm``
call, and a failing step of ``run_stepwise.sh`` must stop the later steps.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
SCRIPTS = [
    "beza/run.sh",
    "comt/run.sh",
    "methyltransferase/run_all.sh",
    "methyltransferase/run_stepwise.sh",
    "toy_system/run.sh",
]


def _stub_mlmm(tmp_path: Path, fail_on: str = "") -> Path:
    """Write a stub ``mlmm`` that logs its subcommand and fails on ``fail_on``."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    calls = tmp_path / "calls.txt"
    stub = bin_dir / "mlmm"
    stub.write_text(
        "#!/bin/sh\n"
        f'echo "$1" >> "{calls}"\n'
        f'[ "$1" = "{fail_on}" ] && exit 3\n'
        "exit 0\n",
        encoding="utf-8",
    )
    stub.chmod(0o755)
    return bin_dir


def _run(script: str, work_dir: Path, bin_dir: Path) -> subprocess.CompletedProcess:
    env = dict(os.environ, PATH=f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    return subprocess.run(
        ["bash", str(EXAMPLES / script), str(work_dir)],
        cwd=work_dir.parent,
        env=env,
        capture_output=True,
        text=True,
    )


def _calls(tmp_path: Path) -> list[str]:
    calls = tmp_path / "calls.txt"
    return calls.read_text(encoding="utf-8").split() if calls.exists() else []


@pytest.mark.parametrize("script", SCRIPTS)
def test_existing_output_directory_stops_the_script(tmp_path: Path, script: str) -> None:
    bin_dir = _stub_mlmm(tmp_path)
    work_dir = tmp_path / "out"
    work_dir.mkdir()
    completed = _run(script, work_dir, bin_dir)
    assert completed.returncode != 0
    assert "File exists" in completed.stderr
    assert _calls(tmp_path) == []


def test_failing_step_stops_run_stepwise(tmp_path: Path) -> None:
    bin_dir = _stub_mlmm(tmp_path, fail_on="extract")
    completed = _run("methyltransferase/run_stepwise.sh", tmp_path / "new", bin_dir)
    assert completed.returncode != 0
    assert _calls(tmp_path) == ["mm-parm", "extract"]
