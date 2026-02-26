"""Helpers for progressive `--help` / `--help-advanced` behavior."""

from __future__ import annotations

import click


def _show_advanced_subcommand_help(
    ctx: click.Context, _param: click.Parameter, value: bool
) -> None:
    """Print subcommand help with advanced options and exit."""
    if not value or ctx.resilient_parsing:
        return

    if getattr(ctx.command, "_advanced_passthrough_help", False):
        try:
            ctx.command.main(args=["--help"], standalone_mode=False)
        except SystemExit as exc:
            code = getattr(exc, "code", 1)
            if code not in (None, 0):
                raise
        ctx.exit()

    hidden = getattr(ctx.command, "_advanced_hidden_options", ())
    restored: list[click.Option] = []
    for opt in hidden:
        if opt.hidden:
            opt.hidden = False
            restored.append(opt)
    try:
        click.echo(ctx.command.get_help(ctx))
    finally:
        for opt in restored:
            opt.hidden = True
    ctx.exit()


def _ensure_help_advanced_option(command: click.Command) -> click.Command:
    """Attach --help-advanced to lazily loaded subcommands when absent."""
    passthrough_help = (
        command.context_settings.get("help_option_names") == []
        if isinstance(command.context_settings, dict)
        else False
    )
    setattr(command, "_advanced_passthrough_help", passthrough_help)

    if any(
        isinstance(param, click.Option) and "--help-advanced" in param.opts
        for param in command.params
    ):
        return command

    option = click.Option(
        ["--help-advanced"],
        is_flag=True,
        is_eager=True,
        expose_value=False,
        callback=_show_advanced_subcommand_help,
        help="Show all options (including advanced settings) and exit.",
    )
    command.params.insert(0, option)
    return command


def _configure_subcommand_help_visibility(
    command_name: str,
    command: click.Command,
    primary_options_by_subcommand: dict[str, frozenset[str]],
) -> click.Command:
    """Hide advanced options from default --help for selected subcommands."""
    if hasattr(command, "_advanced_hidden_options"):
        return command

    primary_options = primary_options_by_subcommand.get(command_name)
    if not primary_options:
        return command

    hidden_options: list[click.Option] = []
    for param in command.params:
        if not isinstance(param, click.Option):
            continue
        names = set(param.opts + param.secondary_opts)
        if names & primary_options:
            continue
        if param.hidden:
            continue
        param.hidden = True
        hidden_options.append(param)

    setattr(command, "_advanced_hidden_options", tuple(hidden_options))
    return command
