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
    from mlmm.cli.preflight import validate_existing_files

    p = tmp_path / "a.txt"
    p.write_text("ok", encoding="utf-8")

    out = validate_existing_files([str(p)], option_name="-i/--input")
    assert out == [p]


def test_validate_existing_files_raises_for_missing(tmp_path: Path):
    from mlmm.cli.preflight import validate_existing_files

    missing = tmp_path / "missing.txt"
    with pytest.raises(click.BadParameter, match="-i/--input path"):
        validate_existing_files([missing], option_name="-i/--input")
