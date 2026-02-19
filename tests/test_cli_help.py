"""CLI --help smoke tests for all subcommands."""

import sys

import pytest

# The CLI modules use Python 3.10+ union syntax (str | None) which
# causes a TypeError on older interpreters at import time.
pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm_toolkit CLI requires Python >= 3.11",
)

from click.testing import CliRunner  # noqa: E402


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def cli_group():
    from mlmm_toolkit.cli import cli
    return cli


def test_main_help(runner, cli_group):
    result = runner.invoke(cli_group, ["--help"])
    assert result.exit_code == 0
    assert "all" in result.output


@pytest.mark.parametrize("subcmd", [
    "add-elem-info",
    "all",
    "define-layer",
    "dft",
    "energy-diagram",
    "extract",
    "fix-altloc",
    "freq",
    "irc",
    "mm-parm",
    "oniom-gaussian",
    "oniom-orca",
    "opt",
    "path-opt",
    "path-search",
    "scan",
    "scan2d",
    "scan3d",
    "trj2fig",
    "tsopt",
])
def test_subcommand_help(runner, cli_group, subcmd):
    result = runner.invoke(cli_group, [subcmd, "--help"])
    assert result.exit_code == 0, "{} --help failed: {}".format(subcmd, result.output)
