"""Verify `mlmm --help` renders the configured command groups in order."""

from __future__ import annotations

from click.testing import CliRunner

from mlmm.cli import cli as root_cli


def test_help_renders_command_groups_in_declared_order() -> None:
    runner = CliRunner()
    result = runner.invoke(root_cli, ["--help"])
    assert result.exit_code == 0, result.output
    out = result.output

    expected_titles = ["Pipelines:", "Pipeline stages:", "Inputs & topology:", "Analysis:"]
    positions = [out.find(title) for title in expected_titles]
    for title, pos in zip(expected_titles, positions):
        assert pos >= 0, f"Missing section '{title}' in --help output:\n{out}"
    assert positions == sorted(positions), (
        f"Command sections rendered out of declared order:\n{positions}\n{out}"
    )

    assert "all" in out
    assert "tsopt" in out
