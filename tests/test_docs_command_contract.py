"""Contracts for the local docs-command and live-CLI checker.

Positive controls exercise the real repository; negative controls inject a
single fault (an invented option, an unavailable lazy command, or a legacy
value-style boolean) and require the exact failure.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = REPO_ROOT / ".github" / "scripts"
for _p in (str(REPO_ROOT), str(SCRIPTS)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import docs_command_contract as dc  # noqa: E402
from mlmm.cli.default_group import build_unavailable_command  # noqa: E402


def _root_cli():
    return dc.load_root_cli()


# --------------------------------------------------------------------------- #
# Retain every command for static validation
# --------------------------------------------------------------------------- #
def test_real_docs_commands_pass_option_validation() -> None:
    commands = dc.extract_docs_commands()
    assert commands, "no docs commands extracted"
    # Bracket-bearing examples must be retained, not dropped.
    assert any(any(m in c.text for m in "[]<>") for c in commands)
    errors = dc.validate_option_names(commands, _root_cli())
    assert errors == [], "\n".join(errors)


def test_quoted_list_literal_is_data_not_synopsis() -> None:
    # A scan command with a quoted list literal stays present and validates -s.
    cmd = dc.AuthoredCommand(
        path=REPO_ROOT / "docs" / "scan.md",
        line=1,
        text='mlmm scan -i pocket.pdb --parm real.parm7 -q 0 -s "[(12,45,2.20)]"',
        executable=dc._classify_executable('-s "[(12,45,2.20)]"'),
    )
    assert dc.validate_option_names([cmd], _root_cli()) == []
    # Negative-float data list must not be mistaken for an option either.
    assert dc.looks_like_data_literal("[-205.1, -190.4, -198.7]")
    assert not dc.looks_like_data_literal("[options]")


def test_invented_option_fails_with_file_and_line() -> None:
    cmd = dc.AuthoredCommand(
        path=REPO_ROOT / "docs" / "scan.md",
        line=42,
        text='mlmm scan -i x.pdb -s "[(12,45,2.20)]" --invented-option',
        executable=False,
    )
    errors = dc.validate_option_names([cmd], _root_cli())
    assert len(errors) == 1
    assert "docs/scan.md:42" in errors[0]
    assert "--invented-option" in errors[0]


def test_unavailable_lazy_command_fails_rather_than_disappears() -> None:
    class _Ctx:
        def close(self) -> None:
            pass

    class _FakeCli:
        def make_context(self, *args, **kwargs):
            return _Ctx()

        def get_command(self, ctx, name):
            return build_unavailable_command(name, ModuleNotFoundError("x", name="x"))

    cmd = dc.AuthoredCommand(REPO_ROOT / "docs" / "all.md", 3, "mlmm all -i R.pdb", True)
    errors = dc.validate_option_names([cmd], _FakeCli())
    assert len(errors) == 1
    assert "unavailable or unknown subcommand" in errors[0]


# --------------------------------------------------------------------------- #
# Live-derived canonical boolean style
# --------------------------------------------------------------------------- #
def test_live_bool_options_nonempty_and_real_surface_clean() -> None:
    live = dc.resolve_live_bool_options(_root_cli())
    assert "--deterministic" in live
    assert "--dump" in live
    assert dc.validate_bool_style(dc.bool_style_sources(), live) == []


def test_value_style_bool_in_authored_guidance_fails(tmp_path: Path) -> None:
    live = dc.resolve_live_bool_options(_root_cli())
    bad = tmp_path / "README.md"
    bad.write_text("Run:\n\n    mlmm all -i R.pdb --deterministic True\n", encoding="utf-8")
    errors = dc.validate_bool_style([bad], live)
    assert len(errors) == 1
    assert "--flag / --no-flag" in errors[0]
    assert "--deterministic True" in errors[0]


def test_comment_line_bool_prose_is_not_flagged(tmp_path: Path) -> None:
    live = dc.resolve_live_bool_options(_root_cli())
    doc = tmp_path / "CONTRIBUTING.md"
    doc.write_text("# --dump on freq: write thermoanalysis.yaml\n", encoding="utf-8")
    assert dc.validate_bool_style([doc], live) == []

def test_every_option_with_an_effective_default_shows_it() -> None:
    """A flag whose real default lives in a config block used to render as unset
    in --help, in the generated reference and in the Colab option list, which
    derives its controls from the live CLI. Click shows a default when
    `show_default` is a non-empty string, or when it is truthy and the declared
    default is not None; the string form keeps the declared default at None so
    `cli_param_overridden` still tells an explicit value from an omission.

    Only genuinely value-less options may stay silent: required or optional
    inputs. Hidden options are checked too -- `--help-advanced` renders them,
    so a hidden knob with a real default must report it there.
    """
    import click

    from mlmm.cli import cli as root_cli

    # Matched on the user-facing flag, not the dest name: an input, an output
    # path or a YAML file has no value to report.
    ALLOWED_SILENT_OPTS = {
        "-i", "--input", "-o", "--out", "--output", "--out-dir", "--ref-pdb",
        "--ref-full-pdb", "--config", "--override", "--calc-file",
        "-c", "--center", "-q", "--charge", "-l",
        "--ligand-charge", "--freeze-atoms", "--dist-freeze", "--model-pdb",
        "--model-indices", "-s", "--scan-lists", "--label-x", "--selected-resn",
        "--modified-residue", "--ref-mode", "--hessian-ref", "--csv",
        "--parm", "--real-parm7",
    }

    ctx = click.Context(root_cli)
    missing = []
    for name in sorted(root_cli.list_commands(ctx)):
        cmd = root_cli.get_command(ctx, name)
        if cmd is None:
            continue
        for param in cmd.params:
            if not isinstance(param, click.Option) or param.expose_value is False:
                continue
            if param.required:
                continue
            if ALLOWED_SILENT_OPTS.intersection(param.opts):
                continue
            sd = param.show_default
            shown = bool(sd) if isinstance(sd, str) else (
                bool(sd) and param.default is not None
            )
            if not shown:
                missing.append(f"{name} {'/'.join(param.opts)}")
    assert missing == [], (
        "these options render no default; give them show_default=<value> "
        "(a string keeps the declared default None): " + ", ".join(missing)
    )
