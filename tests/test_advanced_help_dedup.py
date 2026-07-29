"""One advanced-help implementation for the product.

``all`` and the lazily-loaded subcommands must share the single advanced-help
callback + visibility loop in ``mlmm.cli.help_pages`` — not carry near-identical
private copies that can drift.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm CLI requires Python >= 3.11",
)

from click.testing import CliRunner  # noqa: E402


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def cli_group():
    from mlmm.cli import cli

    return cli


def _help_advanced_option(command):
    for param in command.params:
        if "--help-advanced" in getattr(param, "opts", ()):
            return param
    raise AssertionError("no --help-advanced option on command")


def test_all_uses_the_shared_advanced_help_callback():
    """`all` routes through help_pages, not its own callback."""
    from mlmm.cli.help_pages import _show_advanced_subcommand_help
    from mlmm.workflows import all as all_mod

    # `all` no longer defines a private advanced-help callback.
    assert not hasattr(all_mod, "_show_advanced_help")
    # Its --help-advanced option uses the one shared callback object.
    assert _help_advanced_option(all_mod.cli).callback is _show_advanced_subcommand_help


def test_all_help_progressive_disclosure_and_repeat_idempotent(runner, cli_group):
    """Primary-only default; advanced reveals each option once;
    repeated basic/advanced/basic calls leave basic output + hidden states intact."""
    r1 = runner.invoke(cli_group, ["all", "--help"])
    ra = runner.invoke(cli_group, ["all", "--help-advanced"])
    r2 = runner.invoke(cli_group, ["all", "--help"])
    assert r1.exit_code == 0 and ra.exit_code == 0 and r2.exit_code == 0

    # Basic help is identical across the two calls that bracket an advanced call:
    # the try/finally restoration means advanced never permanently un-hides.
    assert r1.output == r2.output

    # Default hides an advanced option; advanced reveals it exactly once.
    assert "--scan-bias-k" not in r1.output
    assert "--scan-bias-k" not in r2.output
    assert ra.output.count("--scan-bias-k") == 1
    # A primary option is always visible.
    assert "--tsopt" in r1.output and "--tsopt" in ra.output


def test_shared_callback_direct_invocation_is_safe_for_all_and_subcommand():
    """The callback no-ops under resilient parsing / value False for
    both the `all` command and a lazily-loaded subcommand."""
    import click

    from mlmm.cli.help_pages import _show_advanced_subcommand_help
    from mlmm.workflows import all as all_mod
    from mlmm.workflows import opt as opt_mod

    for command in (all_mod.cli, opt_mod.cli):
        # resilient parsing => no-op, no exit raised
        ctx = click.Context(command, resilient_parsing=True)
        assert _show_advanced_subcommand_help(ctx, None, True) is None
        # value False => no-op
        ctx2 = click.Context(command)
        assert _show_advanced_subcommand_help(ctx2, None, False) is None


def _functions_with_hide_loop(root: Path):
    """Return [(file, funcname)] for the option-*discovery* hiding loop.

    That loop iterates ``command.params`` and sets ``param.hidden = True``.  It is
    distinct from the restore loop inside the callback, which iterates the stored
    ``_advanced_hidden_options`` tuple (not ``.params``).
    """
    hits = []
    for py in root.rglob("*.py"):
        try:
            tree = ast.parse(py.read_text(encoding="utf-8"))
        except (SyntaxError, UnicodeDecodeError):
            continue
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            iterates_params = any(
                isinstance(n, ast.For)
                and isinstance(n.iter, ast.Attribute)
                and n.iter.attr == "params"
                for n in ast.walk(node)
            )
            sets_hidden = any(
                isinstance(n, ast.Assign)
                and any(
                    isinstance(t, ast.Attribute) and t.attr == "hidden"
                    for t in n.targets
                )
                and isinstance(n.value, ast.Constant)
                and n.value.value is True
                for n in ast.walk(node)
            )
            if iterates_params and sets_hidden:
                hits.append((py.name, node.name))
    return hits


def test_single_advanced_help_ownership_across_the_product():
    """Exactly one advanced-help callback and one visibility loop."""
    import mlmm

    root = Path(mlmm.__file__).parent

    callback_defs = []
    for py in root.rglob("*.py"):
        tree = ast.parse(py.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in {
                "_show_advanced_subcommand_help",
                "_show_advanced_help",
            }:
                callback_defs.append((py.name, node.name))

    assert callback_defs == [("help_pages.py", "_show_advanced_subcommand_help")], callback_defs

    # The option-hiding discovery loop lives in exactly one place (help_pages).
    loops = _functions_with_hide_loop(root)
    assert loops == [("help_pages.py", "_hide_advanced_options")], loops
