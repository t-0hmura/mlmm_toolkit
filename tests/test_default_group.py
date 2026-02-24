"""Unit tests for shared DefaultGroup behavior."""

from __future__ import annotations

import click
from click.testing import CliRunner

from mlmm_toolkit.advanced_help import (
    _configure_subcommand_help_visibility,
    _ensure_help_advanced_option,
)
from mlmm_toolkit.bool_compat import normalize_bool_argv as _normalize_bool_argv_impl
from mlmm_toolkit.default_group import DefaultGroup


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
) -> click.Group:
    @click.group(
        cls=DefaultGroup,
        default=default,
        lazy_subcommands=lazy_subcommands or {},
        command_bool_value_options={},
        command_bool_toggle_options={},
        command_bool_toggle_negative_aliases={},
        parser_wrapper_subcommands=frozenset(),
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


def test_legacy_bool_syntax_emits_deprecation_message() -> None:
    cli = _make_group(normalize_bool_argv=_normalize_with_legacy_flag)

    @cli.command()
    def ping() -> None:
        click.echo("ping-called")

    runner = CliRunner()
    result = runner.invoke(cli, ["ping"])
    assert result.exit_code == 0
    assert "ping-called" in result.output
    assert "Legacy bool syntax '--flag True/False'" in result.output


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
            "broken": (".__this_module_should_not_exist__", "cli", "Broken command")
        }
    )

    runner = CliRunner()
    result = runner.invoke(cli, ["broken"])
    assert result.exit_code != 0
    assert "Command 'broken' is unavailable because the module could not be imported." in result.output
    assert "Missing dependency:" in result.output


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
    assert "Legacy bool syntax '--flag True/False'" in result_false.output

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
    assert "Legacy bool syntax '--flag True/False'" in result_value_false.output

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
    assert "Legacy bool syntax '--flag True/False'" in result_false.output

    result_no = runner.invoke(cli, ["define-layer", "--no-one-based"])
    assert result_no.exit_code == 0
    assert "one_based=False" in result_no.output
