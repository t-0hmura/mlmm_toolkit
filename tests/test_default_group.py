"""Unit tests for shared DefaultGroup behavior."""

from __future__ import annotations

from pathlib import Path

import click
import pytest
from click.testing import CliRunner

from mlmm.cli.help_pages import (
    _configure_subcommand_help_visibility,
    _ensure_help_advanced_option,
)
from mlmm.cli.bool_compat import normalize_bool_argv as _normalize_bool_argv_impl
from mlmm.cli.default_group import DefaultGroup


def _normalize_passthrough(args, *_):
    return args, False


def _normalize_with_legacy_flag(args, *_):
    return args, True


def _make_group(
    *,
    default: str | None = None,
    lazy_subcommands: dict[str, tuple[str, str, str]] | None = None,
    normalize_bool_argv=_normalize_passthrough,
    primary_help_options: dict[str, frozenset[str]] | None = None,
    parser_wrapper_subcommands: frozenset[str] | None = None,
    parser_wrapper_bool_option_providers=None,
) -> click.Group:
    @click.group(
        cls=DefaultGroup,
        default=default,
        lazy_subcommands=lazy_subcommands or {},
        command_bool_value_options={},
        command_bool_toggle_options={},
        command_bool_toggle_negative_aliases={},
        parser_wrapper_subcommands=parser_wrapper_subcommands or frozenset(),
        parser_wrapper_bool_option_providers=parser_wrapper_bool_option_providers or {},
        subcommand_primary_help_options=primary_help_options or {},
        normalize_bool_argv=normalize_bool_argv,
        ensure_help_advanced_option=_ensure_help_advanced_option,
        configure_subcommand_help_visibility=_configure_subcommand_help_visibility,
    )
    def cli() -> None:
        """Test CLI."""

    return cli


def test_default_subcommand_is_inserted_when_no_args() -> None:
    cli = _make_group(default="all")

    @cli.command(name="all")
    def all_cmd() -> None:
        click.echo("all-called")

    runner = CliRunner()
    result = runner.invoke(cli, [])
    assert result.exit_code == 0
    assert "all-called" in result.output


def test_output_directory_command_writes_run_log(tmp_path: Path) -> None:
    cli = _make_group()

    @cli.command(name="work")
    @click.option("-o", "--out-dir", default="result_work")
    def work_cmd(out_dir: str) -> None:
        Path(out_dir).mkdir(parents=True)
        click.echo("work-complete")

    out = tmp_path / "result"
    result = CliRunner().invoke(cli, ["work", "-o", str(out)])

    assert result.exit_code == 0
    log = (out / "run.log").read_text(encoding="utf-8")
    assert "[command] mlmm work" in log
    assert "work-complete" in log


def test_help_does_not_trigger_default_subcommand() -> None:
    cli = _make_group(default="all")

    @cli.command(name="all")
    def all_cmd() -> None:
        click.echo("all-called")

    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "all-called" not in result.output
    assert "Commands:" in result.output


def test_legacy_bool_syntax_does_not_emit_deprecation_message() -> None:
    cli = _make_group(normalize_bool_argv=_normalize_with_legacy_flag)

    @cli.command()
    def ping() -> None:
        click.echo("ping-called")

    runner = CliRunner()
    result = runner.invoke(cli, ["ping"])
    assert result.exit_code == 0
    assert "ping-called" in result.output
    assert "Legacy bool syntax '--flag True/False'" not in result.output


def test_help_advanced_shows_hidden_options() -> None:
    primary = {"scan": frozenset({"--input", "--help-advanced"})}
    cli = _make_group(primary_help_options=primary)

    @cli.command(name="scan")
    @click.option("--input")
    @click.option("--advanced-opt")
    def scan_cmd(input: str | None, advanced_opt: str | None) -> None:
        _ = (input, advanced_opt)
        click.echo("scan-called")

    runner = CliRunner()
    short = runner.invoke(cli, ["scan", "--help"])
    assert short.exit_code == 0
    assert "--help-advanced" in short.output
    assert "--input" in short.output
    assert "--advanced-opt" not in short.output

    advanced = runner.invoke(cli, ["scan", "--help-advanced"])
    assert advanced.exit_code == 0
    assert "--advanced-opt" in advanced.output


def test_lazy_import_failure_is_reported_as_click_exception() -> None:
    cli = _make_group(
        lazy_subcommands={
            "broken": (
                "optional_mlmm_plugin.__this_module_should_not_exist__",
                "cli",
                "Broken command",
            )
        }
    )

    runner = CliRunner()
    result = runner.invoke(cli, ["broken"])
    assert result.exit_code != 0
    assert "Command 'broken' is unavailable because the module could not be imported." in result.output
    assert "Missing dependency:" in result.output


@pytest.mark.parametrize(
    "missing_name",
    [
        "mlmm.internal_missing",
        "pysisyphus.internal_missing",
        "thermoanalysis.internal_missing",
        "hessian_ff.internal_missing",
    ],
)
def test_lazy_internal_missing_module_is_not_masked(
    monkeypatch: pytest.MonkeyPatch,
    missing_name: str,
) -> None:
    from mlmm.cli import default_group

    cli = _make_group(
        lazy_subcommands={"broken": ("mlmm.workflows.broken", "cli", "Broken command")}
    )

    def _raise_internal_missing(*_args, **_kwargs):
        raise ModuleNotFoundError(
            f"No module named {missing_name!r}", name=missing_name
        )

    monkeypatch.setattr(default_group.importlib, "import_module", _raise_internal_missing)
    with pytest.raises(ModuleNotFoundError) as exc_info:
        cli.get_command(click.Context(cli), "broken")
    assert exc_info.value.name == missing_name


def test_lazy_plain_import_error_is_not_masked(monkeypatch: pytest.MonkeyPatch) -> None:
    from mlmm.cli import default_group

    cli = _make_group(
        lazy_subcommands={"broken": ("external_plugin", "cli", "Broken command")}
    )

    def _raise_plain_import_error(*_args, **_kwargs):
        raise ImportError("missing symbol inside imported module")

    monkeypatch.setattr(default_group.importlib, "import_module", _raise_plain_import_error)
    with pytest.raises(ImportError, match="missing symbol"):
        cli.get_command(click.Context(cli), "broken")


def test_bool_toggle_accepts_value_style_syntax_via_auto_detection() -> None:
    cli = _make_group(normalize_bool_argv=_normalize_bool_argv_impl)

    @cli.command(name="scan")
    @click.option("--detect-layer/--no-detect-layer", default=True)
    def scan_cmd(detect_layer: bool) -> None:
        click.echo(f"detect_layer={detect_layer}")

    runner = CliRunner()
    result_false = runner.invoke(cli, ["scan", "--detect-layer", "False"])
    assert result_false.exit_code == 0
    assert "detect_layer=False" in result_false.output

    result_no = runner.invoke(cli, ["scan", "--no-detect-layer"])
    assert result_no.exit_code == 0
    assert "detect_layer=False" in result_no.output


def test_single_flag_accepts_no_prefix_and_value_style_syntax() -> None:
    cli = _make_group(normalize_bool_argv=_normalize_bool_argv_impl)

    @cli.command(name="add-elem-info")
    @click.option("--overwrite", is_flag=True, default=False)
    def add_elem_info_cmd(overwrite: bool) -> None:
        click.echo(f"overwrite={overwrite}")

    runner = CliRunner()
    result_value_false = runner.invoke(cli, ["add-elem-info", "--overwrite", "False"])
    assert result_value_false.exit_code == 0
    assert "overwrite=False" in result_value_false.output

    result_no = runner.invoke(cli, ["add-elem-info", "--no-overwrite"])
    assert result_no.exit_code == 0
    assert "overwrite=False" in result_no.output


def test_toggle_with_non_no_negative_alias_accepts_no_prefix_and_values() -> None:
    cli = _make_group(normalize_bool_argv=_normalize_bool_argv_impl)

    @cli.command(name="define-layer")
    @click.option("--one-based/--zero-based", default=True)
    def define_layer_cmd(one_based: bool) -> None:
        click.echo(f"one_based={one_based}")

    runner = CliRunner()
    result_false = runner.invoke(cli, ["define-layer", "--one-based", "False"])
    assert result_false.exit_code == 0
    assert "one_based=False" in result_false.output

    result_no = runner.invoke(cli, ["define-layer", "--no-one-based"])
    assert result_no.exit_code == 0
    assert "one_based=False" in result_no.output


def test_parser_wrapper_bool_provider_enables_value_style_normalization() -> None:
    cli = _make_group(
        normalize_bool_argv=_normalize_bool_argv_impl,
        parser_wrapper_subcommands=frozenset({"extract"}),
        parser_wrapper_bool_option_providers={"extract": lambda: frozenset({"--verbose"})},
    )

    @cli.command(
        name="extract",
        context_settings={"ignore_unknown_options": True, "allow_extra_args": True},
    )
    @click.pass_context
    def extract_cmd(ctx: click.Context) -> None:
        click.echo("|".join(ctx.args))

    runner = CliRunner()
    result = runner.invoke(cli, ["extract", "--verbose", "False"])
    assert result.exit_code == 0
    assert "--no-verbose" in result.output
    assert "--verbose|False" not in result.output
