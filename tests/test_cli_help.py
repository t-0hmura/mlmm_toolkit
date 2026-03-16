"""CLI --help smoke tests for all subcommands."""

import sys

import pytest

# The CLI modules use Python 3.10+ union syntax (str | None) which
# causes a TypeError on older interpreters at import time.
pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm CLI requires Python >= 3.11",
)

from click.testing import CliRunner  # noqa: E402


SUBCOMMANDS = [
    "add-elem-info",
    "all",
    "define-layer",
    "dft",
    "energy-diagram",
    "extract",
    "fix-altloc",
    "freq",
    "oniom-import",
    "irc",
    "mm-parm",
    "oniom-export",
    "opt",
    "path-opt",
    "path-search",
    "scan",
    "scan2d",
    "scan3d",
    "trj2fig",
    "tsopt",
]

SCAN_SUBCOMMANDS = [
    ("scan", "--preopt"),
    ("scan2d", "--preopt"),
    ("scan3d", "--preopt"),
]

CALC_SUBCOMMANDS = [
    ("opt", "--opt-mode"),
    ("path-opt", "--mep-mode"),
    ("path-search", "--mep-mode"),
    ("tsopt", "--opt-mode"),
    ("freq", "--temperature"),
    ("irc", "--step-size"),
    ("dft", "--func-basis"),
]

UTILITY_SUBCOMMANDS = [
    ("mm-parm", "--out-prefix", "--keep-temp"),
    ("define-layer", "--model-pdb", "--radius-partial-hessian"),
    ("add-elem-info", "-o, --out", "--overwrite"),
    ("trj2fig", "--unit", "--reverse-x"),
    ("energy-diagram", "-o, --output", "--label-x"),
    ("oniom-export", "--mode", "--orcaff"),
]


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def cli_group():
    from mlmm.cli import cli
    return cli


def _has_option_header(output: str, option_prefix: str) -> bool:
    for line in output.splitlines():
        stripped = line.lstrip()
        if not stripped.startswith(option_prefix):
            continue
        tail = stripped[len(option_prefix):]
        if (not tail) or tail[0].isspace() or tail[0] in {",", "/"}:
            return True
    return False


def test_main_help(runner, cli_group):
    result = runner.invoke(cli_group, ["--help"])
    assert result.exit_code == 0
    assert "all" in result.output
    assert "opt" in result.output


def test_all_help_progressive_disclosure(runner, cli_group):
    result = runner.invoke(cli_group, ["all", "--help"])
    assert result.exit_code == 0
    assert "--help-advanced" in result.output
    assert "--scan-bias-k" not in result.output


def test_all_help_advanced_shows_hidden_options(runner, cli_group):
    result = runner.invoke(cli_group, ["all", "--help-advanced"])
    assert result.exit_code == 0
    assert "--scan-bias-k" in result.output
    assert "--freq-temperature" in result.output
    assert "--opt-mode-post" in result.output
    assert "--sopt-mode" not in result.output


def test_path_search_help_advanced_uses_opt_mode_name(runner, cli_group):
    result = runner.invoke(cli_group, ["path-search", "--help-advanced"])
    assert result.exit_code == 0
    assert "--opt-mode" in result.output
    assert "--sopt-mode" not in result.output


@pytest.mark.parametrize("subcmd,legacy_header", SCAN_SUBCOMMANDS)
def test_scan_family_help_progressive_disclosure(runner, cli_group, subcmd, legacy_header):
    result = runner.invoke(cli_group, [subcmd, "--help"])
    assert result.exit_code == 0
    assert "--help-advanced" in result.output
    assert legacy_header not in result.output


@pytest.mark.parametrize("subcmd,legacy_header", SCAN_SUBCOMMANDS)
def test_scan_family_help_advanced_shows_hidden_options(runner, cli_group, subcmd, legacy_header):
    result = runner.invoke(cli_group, [subcmd, "--help-advanced"])
    assert result.exit_code == 0


@pytest.mark.parametrize("subcmd,core_opt", CALC_SUBCOMMANDS)
def test_calc_family_help_progressive_disclosure(runner, cli_group, subcmd, core_opt):
    result = runner.invoke(cli_group, [subcmd, "--help"])
    assert result.exit_code == 0
    assert "--help-advanced" in result.output


@pytest.mark.parametrize("subcmd,_core_opt", CALC_SUBCOMMANDS)
def test_calc_family_help_advanced_shows_hidden_options(runner, cli_group, subcmd, _core_opt):
    result = runner.invoke(cli_group, [subcmd, "--help-advanced"])
    assert result.exit_code == 0


@pytest.mark.parametrize("subcmd,core_opt,hidden_opt", UTILITY_SUBCOMMANDS)
def test_utility_help_progressive_disclosure(
    runner, cli_group, subcmd, core_opt, hidden_opt
):
    result = runner.invoke(cli_group, [subcmd, "--help"])
    assert result.exit_code == 0
    assert "--help-advanced" in result.output
    assert _has_option_header(result.output, core_opt)
    assert not _has_option_header(result.output, hidden_opt)


@pytest.mark.parametrize("subcmd,_core_opt,hidden_opt", UTILITY_SUBCOMMANDS)
def test_utility_help_advanced_shows_hidden_options(
    runner, cli_group, subcmd, _core_opt, hidden_opt
):
    result = runner.invoke(cli_group, [subcmd, "--help-advanced"])
    assert result.exit_code == 0
    assert _has_option_header(result.output, hidden_opt)


def test_extract_help_progressive_disclosure(runner, cli_group):
    short = runner.invoke(cli_group, ["extract", "--help"])
    assert short.exit_code == 0
    assert "extract [OPTIONS]" in short.output
    assert _has_option_header(short.output, "-i, --input")
    assert _has_option_header(short.output, "-c, --center")
    assert _has_option_header(short.output, "--help-advanced")
    assert not _has_option_header(short.output, "--selected-resn")

    advanced = runner.invoke(cli_group, ["extract", "--help-advanced"])
    assert advanced.exit_code == 0
    assert _has_option_header(advanced.output, "--selected-resn")
    assert _has_option_header(advanced.output, "--include-H2O")


def test_fix_altloc_help_progressive_disclosure(runner, cli_group):
    short = runner.invoke(cli_group, ["fix-altloc", "--help"])
    assert short.exit_code == 0
    assert "fix-altloc [OPTIONS]" in short.output
    assert _has_option_header(short.output, "--recursive")
    assert _has_option_header(short.output, "--help-advanced")

    advanced = runner.invoke(cli_group, ["fix-altloc", "--help-advanced"])
    assert advanced.exit_code == 0
    assert _has_option_header(advanced.output, "--overwrite")
    assert _has_option_header(advanced.output, "--force")


@pytest.mark.parametrize("subcmd", SUBCOMMANDS)
def test_subcommand_help(runner, cli_group, subcmd):
    result = runner.invoke(cli_group, [subcmd, "--help"])
    assert result.exit_code == 0, "{} --help failed: {}".format(subcmd, result.output)


@pytest.mark.parametrize("subcmd", SUBCOMMANDS)
def test_subcommand_help_advanced(runner, cli_group, subcmd):
    result = runner.invoke(cli_group, [subcmd, "--help-advanced"])
    assert result.exit_code == 0, "{} --help-advanced failed: {}".format(subcmd, result.output)
