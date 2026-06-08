"""CLI-level bool compatibility regression tests."""

from __future__ import annotations

import click
from click.testing import CliRunner

from mlmm.cli import cli as root_cli


def _is_unavailable_command(cmd: click.Command | None) -> bool:
    if cmd is None:
        return True
    help_text = (cmd.help or "").strip()
    return help_text.startswith("[Unavailable]")


def test_all_bool_options_accept_toggle_and_value_styles() -> None:
    runner = CliRunner()
    ctx = click.Context(root_cli)

    tested = 0
    for command_name in root_cli.list_commands(ctx):
        command = root_cli.get_command(ctx, command_name)
        if _is_unavailable_command(command):
            continue

        value_opts, toggle_opts, _aliases, single_opts = root_cli._resolve_bool_options(
            ctx, command_name
        )
        option_names = sorted(
            {
                *value_opts,
                *toggle_opts,
                *single_opts,
            }
        )
        for opt in option_names:
            if not opt.startswith("--"):
                continue

            no_opt = f"--no-{opt[2:]}"

            result_flag = runner.invoke(root_cli, [command_name, opt, "--help"])
            assert result_flag.exit_code == 0, (
                f"{command_name} should accept {opt}. Output:\n{result_flag.output}"
            )

            result_no_flag = runner.invoke(root_cli, [command_name, no_opt, "--help"])
            assert result_no_flag.exit_code == 0, (
                f"{command_name} should accept {no_opt}. Output:\n{result_no_flag.output}"
            )

            result_value_false = runner.invoke(
                root_cli, [command_name, opt, "False", "--help"]
            )
            assert result_value_false.exit_code == 0, (
                f"{command_name} should accept '{opt} False'. Output:\n"
                f"{result_value_false.output}"
            )

            # Positive value form: `--flag True` (previously only False
            # was tested, but the bool_compat helper accepts both — make sure
            # neither form regresses).
            result_value_true = runner.invoke(
                root_cli, [command_name, opt, "True", "--help"]
            )
            assert result_value_true.exit_code == 0, (
                f"{command_name} should accept '{opt} True'. Output:\n"
                f"{result_value_true.output}"
            )

            tested += 1

    assert tested > 0


def test_single_flag_negative_form_value_style() -> None:
    """`--no-FLAG True/False` must not regress on single-flag bools.

    Single-flag bools (registered in `_COMMAND_BOOL_SINGLE_FLAG_OPTIONS`)
    expose only the positive `--FLAG` name to Click; the bool_compat shim
    synthesises the `--no-FLAG` negative form. The legacy value style
    `--no-FLAG True/False` must collapse to one of:
        - `--no-FLAG False` -> `--FLAG`   (keep on)
        - `--no-FLAG True`  -> (omitted)  (turn off; matches the default)
    Emitting the literal `--no-FLAG` would crash Click ("no such option").
    """
    runner = CliRunner()
    ctx = click.Context(root_cli)

    tested = 0
    for command_name in root_cli.list_commands(ctx):
        command = root_cli.get_command(ctx, command_name)
        if _is_unavailable_command(command):
            continue
        _value, _toggle, _aliases, single_opts = root_cli._resolve_bool_options(
            ctx, command_name
        )
        for opt in single_opts:
            if not opt.startswith("--"):
                continue
            no_opt = f"--no-{opt[2:]}"
            for literal in ("True", "False"):
                result = runner.invoke(
                    root_cli, [command_name, no_opt, literal, "--help"]
                )
                assert result.exit_code == 0, (
                    f"{command_name} should accept '{no_opt} {literal}'. "
                    f"Output:\n{result.output}"
                )
                tested += 1

    # `_COMMAND_BOOL_SINGLE_FLAG_OPTIONS` may be empty in some snapshots;
    # don't require non-zero coverage so the test stays useful when the
    # registry shrinks.
    assert tested >= 0
