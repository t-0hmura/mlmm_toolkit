"""Regression for all-command multiplicity option naming."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from mlmm.cli import cli as root_cli


def test_all_accepts_multiplicity_option() -> None:
    fixture_dir = (
        Path(__file__).resolve().parents[1]
        / "hessian_ff"
        / "tests"
        / "data"
        / "small"
    )
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(
            root_cli,
            [
                "all",
                "-i",
                str(fixture_dir / "complex.pdb"),
                "--parm",
                str(fixture_dir / "complex.parm7"),
                "--model-pdb",
                str(fixture_dir / "complex.pdb"),
                "--tsopt",
                "--dry-run",
                "--multiplicity",
                "1",
            ],
        )
    assert result.exit_code == 0, result.output
