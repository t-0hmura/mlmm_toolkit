"""Tests for shared CLI preflight checks."""

from __future__ import annotations

import sys
from pathlib import Path

import click
import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm requires Python >= 3.11",
)


def test_validate_existing_files_returns_paths(tmp_path: Path):
    from mlmm.preflight import validate_existing_files

    p = tmp_path / "a.txt"
    p.write_text("ok", encoding="utf-8")

    out = validate_existing_files([str(p)], option_name="-i/--input")
    assert out == [p]


def test_validate_existing_files_raises_for_missing(tmp_path: Path):
    from mlmm.preflight import validate_existing_files

    missing = tmp_path / "missing.txt"
    with pytest.raises(click.BadParameter, match="-i/--input path"):
        validate_existing_files([missing], option_name="-i/--input")


def test_ensure_commands_available_passes_for_shell_builtin_path():
    from mlmm.preflight import ensure_commands_available

    # /bin/sh should exist on Linux HPC environments.
    ensure_commands_available(["sh"], context="test")


def test_ensure_commands_available_raises_for_missing_command():
    from mlmm.preflight import ensure_commands_available

    with pytest.raises(click.ClickException, match="Missing required command"):
        ensure_commands_available(["cmd_that_should_not_exist_12345"], context="test")


def test_ensure_python_modules_available_passes_for_standard_modules():
    from mlmm.preflight import ensure_python_modules_available

    ensure_python_modules_available(["json", "pathlib"], context="test")


def test_ensure_python_modules_available_raises_for_missing_module():
    from mlmm.preflight import ensure_python_modules_available

    with pytest.raises(click.ClickException, match="Missing required Python module"):
        ensure_python_modules_available(
            {"module_that_should_not_exist_12345": "pip install module_that_should_not_exist_12345"},
            context="test",
        )
