"""Regression for all-command multiplicity option naming."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from mlmm_toolkit.cli import cli as root_cli


def test_all_accepts_multiplicity_option() -> None:
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("single.pdb").write_text("END\n", encoding="utf-8")
        result = runner.invoke(
            root_cli,
            ["all", "-i", "single.pdb", "--tsopt", "True", "--dry-run", "True", "--multiplicity", "1"],
        )
    assert result.exit_code == 0, result.output
