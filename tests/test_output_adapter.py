"""Tests for dependency-light private-tag output routing."""

from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

import pytest

from mlmm.core import output
from mlmm.io.trj2fig import write_csv


def test_library_output_adapter_works_before_cli_bootstrap() -> None:
    code = """
from mlmm.core.output import emit
from mlmm.core.utils import echo_run_summary, emit_optimizer_terminal_status
emit('direct-detail', detail=True)
emit('forced', force=True, raw_path=True)
echo_run_summary({'input': 'x.pdb'})
emit_optimizer_terminal_status('opt', converged=True, cycles=2, max_cycles=5)
"""
    proc = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert "direct-detail" in proc.stdout
    assert "forced" in proc.stdout
    assert "[input] x.pdb" in proc.stdout
    assert "[opt] Converged!" in proc.stdout


def test_direct_trj2fig_csv_call_uses_adapter(tmp_path: Path, capsys) -> None:
    destination = tmp_path / "energy.csv"
    write_csv(destination, [-1.0, -0.5], [0.0, 1.0], "kcal", True)
    assert destination.exists()
    assert "Saved CSV" in capsys.readouterr().out


def test_adapter_tracks_the_current_echo_marker_not_stale_bootstrap_state(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls = []

    def marked_echo(message=None, **kwargs):
        calls.append((message, kwargs))

    setattr(marked_echo, output._TAG_AWARE_MARKER, True)
    monkeypatch.setattr(output.click, "echo", marked_echo)
    output.emit(
        "marked",
        narrative=True,
        detail=True,
        force=True,
        raw_path=True,
        err=True,
    )
    assert calls[-1][1] == {
        "narrative": True,
        "detail": True,
        "force": True,
        "raw_path": True,
        "err": True,
    }

    def native_echo(message=None, **kwargs):
        calls.append((message, kwargs))

    monkeypatch.setattr(output.click, "echo", native_echo)
    output.emit(
        "native",
        narrative=True,
        detail=True,
        force=True,
        raw_path=True,
        err=True,
    )
    assert calls[-1] == ("native", {"err": True})


def test_no_native_click_echo_receives_private_or_dynamic_tags() -> None:
    private_tags = {"narrative", "detail", "force", "raw_path"}
    package_root = Path(__file__).resolve().parents[1] / "mlmm"
    violations = []
    for path in package_root.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
                continue
            if node.func.attr != "echo" or not isinstance(node.func.value, ast.Name):
                continue
            if node.func.value.id not in {"click", "_click"}:
                continue
            keywords = {kw.arg for kw in node.keywords if kw.arg is not None}
            has_dynamic_kwargs = any(kw.arg is None for kw in node.keywords)
            if (keywords & private_tags or has_dynamic_kwargs) and path.name != "output.py":
                violations.append((path.relative_to(package_root), node.lineno))
    assert violations == []


def test_lower_layers_depend_on_output_adapter_not_core_utils() -> None:
    package_root = Path(__file__).resolve().parents[1] / "mlmm"
    violations = []
    for directory in ("backends", "domain", "io"):
        for path in (package_root / directory).rglob("*.py"):
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            for node in ast.walk(tree):
                if not isinstance(node, ast.ImportFrom) or node.module != "mlmm.core.utils":
                    continue
                if any(alias.name == "emit" for alias in node.names):
                    violations.append(path.relative_to(package_root))
    assert violations == []
