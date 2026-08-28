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
    "sp",
    "trj2fig",
    "tsopt",
    "bond-summary",
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
    ("sp", "--hess"),
]

UTILITY_SUBCOMMANDS = [
    ("mm-parm", "-o, --out-prefix", "--keep-temp"),
    ("define-layer", "--model-pdb", "--one-based"),
    ("add-elem-info", "-o, --out", "--overwrite"),
    ("trj2fig", "--unit", "--backend-model"),
    ("energy-diagram", "-o, --output", "--label-x"),
    ("oniom-export", "--mode", "--orcaff"),
    ("bond-summary", "--device", "--json"),
]

SHARED_PRIMARY_SCIENTIFIC_OPTIONS = [
    ("add-elem-info", "--inplace"),
    ("fix-altloc", "--inplace"),
    ("all", "--refine-path"),
    ("all", "--scan-max-step-size"),
    ("all", "--mep-mode"),
    ("all", "--max-nodes"),
    ("all", "-r"),
    ("all", "-m"),
    ("all", "--opt-mode"),
    ("all", "--opt-mode-post"),
    ("all", "--thresh"),
    ("all", "--parm"),
    ("all", "--model-pdb"),
    ("all", "--detect-layer"),
    ("all", "--ref-pdb"),
    ("opt", "--thresh"),
    ("opt", "--bias-k"),
    ("opt", "--dist-freeze"),
    ("opt", "--dump"),
    ("opt", "--ref-pdb"),
    ("scan", "--bias-k"),
    ("scan", "--max-step-size"),
    ("scan", "--thresh"),
    ("scan2d", "--bias-k"),
    ("scan2d", "--max-step-size"),
    ("scan2d", "--thresh"),
    ("scan2d", "--ref-pdb"),
    ("scan3d", "--bias-k"),
    ("scan3d", "--max-step-size"),
    ("scan3d", "--thresh"),
    ("scan3d", "--ref-pdb"),
    ("sp", "--ref-pdb"),
    ("tsopt", "--thresh"),
    ("tsopt", "--dump"),
    ("tsopt", "--ref-pdb"),
    ("freq", "--ref-pdb"),
    ("irc", "--ref-pdb"),
    ("dft", "--ref-pdb"),
    ("path-opt", "--ref-pdb"),
    ("path-search", "--ref-pdb"),
    ("define-layer", "--radius-freeze"),
    ("dft", "--engine"),
    ("trj2fig", "--reverse-x"),
]

SHARED_ADVANCED_SCIENTIFIC_OPTIONS = [
    ("irc", "--never-stop"),
    ("all", "--dry-run"),
    ("all", "--scan-bias-k"),
    ("all", "--scan-relax-max-cycles"),
    ("all", "--max-cycles-gsm"),
    ("all", "--gsm-param"),
    ("all", "--max-cycles-dmf"),
    ("all", "--tsopt-max-cycles"),
    ("all", "--hessian-calc-mode"),
    ("scan", "--max-cycles"),
    ("scan", "--relax-max-cycles"),
    ("scan2d", "--relax-max-cycles"),
    ("scan3d", "--relax-max-cycles"),
    ("path-opt", "--gsm-param"),
    ("path-search", "--gsm-param"),
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


@pytest.mark.parametrize(("subcommand", "option"), SHARED_PRIMARY_SCIENTIFIC_OPTIONS)
def test_shared_scientific_options_stay_in_primary_help(
    runner, cli_group, subcommand: str, option: str
) -> None:
    result = runner.invoke(cli_group, [subcommand, "--help"])
    assert result.exit_code == 0, result.output
    assert _has_option_header(result.output, option), result.output


@pytest.mark.parametrize(("subcommand", "option"), SHARED_ADVANCED_SCIENTIFIC_OPTIONS)
def test_selected_scientific_options_stay_in_advanced_help(
    runner, cli_group, subcommand: str, option: str
) -> None:
    primary = runner.invoke(cli_group, [subcommand, "--help"])
    assert primary.exit_code == 0, primary.output
    assert not _has_option_header(primary.output, option), primary.output

    advanced = runner.invoke(cli_group, [subcommand, "--help-advanced"])
    assert advanced.exit_code == 0, advanced.output
    assert _has_option_header(advanced.output, option), advanced.output


def test_main_help(runner, cli_group):
    result = runner.invoke(cli_group, ["--help"])
    assert result.exit_code == 0
    assert "all" in result.output
    assert "opt" in result.output


def test_define_layer_rejects_removed_partial_hessian_option(runner, cli_group):
    result = runner.invoke(
        cli_group, ["define-layer", "--radius-partial-hessian", "4.0"]
    )
    assert result.exit_code != 0
    assert "No such option" in result.output
    assert "--radius-partial-hessian" in result.output


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
    assert "--freeze-atoms" in result.output
    assert "--sopt-mode" not in result.output


def test_path_search_help_advanced_omits_fixed_optimizer_selector(runner, cli_group):
    result = runner.invoke(cli_group, ["path-search", "--help-advanced"])
    assert result.exit_code == 0
    assert "--opt-mode" not in result.output
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
    assert _has_option_header(result.output, legacy_header), result.output


@pytest.mark.parametrize("subcmd,core_opt", CALC_SUBCOMMANDS)
def test_calc_family_help_progressive_disclosure(runner, cli_group, subcmd, core_opt):
    result = runner.invoke(cli_group, [subcmd, "--help"])
    assert result.exit_code == 0
    assert "--help-advanced" in result.output


@pytest.mark.parametrize("subcmd,core_opt", CALC_SUBCOMMANDS)
def test_calc_family_help_advanced_shows_hidden_options(runner, cli_group, subcmd, core_opt):
    result = runner.invoke(cli_group, [subcmd, "--help-advanced"])
    assert result.exit_code == 0
    assert _has_option_header(result.output, core_opt), result.output


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
    assert _has_option_header(short.output, "--selected-resn")

    advanced = runner.invoke(cli_group, ["extract", "--help-advanced"])
    assert advanced.exit_code == 0
    assert _has_option_header(advanced.output, "--selected-resn")


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
