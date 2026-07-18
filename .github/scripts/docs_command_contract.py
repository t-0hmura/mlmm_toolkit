#!/usr/bin/env python3
"""Local authored-command / live-CLI contract for mlmm docs and examples.

This is a field-isomorphic copy of pdb2reaction's docs-command checker shape.
mlmm keeps its own independent implementation on purpose: the two products are
independent (no cross-repo import) and only the small data shape and the
falsifier design are shared. Product-specific constants (``TOOL_NAME``, the
root CLI import, example paths) stay local to this module.

The module retains *every* authored command for static validation and
classifies execution eligibility separately. A quoted token that parses as a
Python list/tuple literal is DATA (e.g. ``-s "[(12,45,2.20)]"``), not synopsis
notation, and must not be mistaken for an option name.
"""

from __future__ import annotations

import ast
import re
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import click

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_ROOT = REPO_ROOT / "docs"
TOOL_NAME = "mlmm"
CLI_MODULE = "mlmm"

# Markdown fences whose body carries shell command examples.
_CODE_LANGS = {"", "bash", "sh", "shell", "console"}
# Shell control tokens that terminate a single command's option list.
_SHELL_BREAKS = {"|", "||", "&&", ";"}
# Synopsis notation: a command carrying any of these cannot be executed as-is
# (angle-bracket placeholders, optional-argument brackets, or data literals that
# reference specific atoms/indices the synthetic dry-run fixture cannot satisfy).
_SYNOPSIS_MARKS = ("<", ">", "[", "]")
# Legacy value-style boolean spelling in authored guidance (rejected).
_BOOL_VALUE_RE = re.compile(
    r"(?P<option>--[a-z0-9][a-z0-9-]*)(?:[ \t]+|=)"
    r"(?P<value>true|false|yes|no|on|off)\b",
    re.IGNORECASE,
)
# Landing / legacy pages that intentionally document the deprecated value-style
# boolean form for backward-compatibility guidance.
_LEGACY_BOOL_PAGES = frozenset(
    {
        Path("docs/cli-conventions.md"),
        Path("docs/ja/cli-conventions.md"),
        Path("docs/concepts.md"),
        Path("docs/ja/concepts.md"),
        Path("docs/glossary.md"),
        Path("docs/ja/glossary.md"),
    }
)


@dataclass(frozen=True)
class AuthoredCommand:
    """A single command extracted from an authored doc or example script."""

    path: Path
    line: int
    text: str
    executable: bool

    @property
    def rel(self) -> str:
        try:
            return str(self.path.relative_to(REPO_ROOT))
        except ValueError:
            return str(self.path)

    @property
    def location(self) -> str:
        return f"{self.rel}:{self.line}"


def load_root_cli() -> click.Group:
    """Import and return the live root CLI group (triggers lazy discovery)."""
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))
    from mlmm.cli import cli as root_cli  # noqa: E402

    return root_cli


# --------------------------------------------------------------------------- #
# Explicit authored/public roots
# --------------------------------------------------------------------------- #
def public_shell_examples() -> list[Path]:
    """Return the advertised example shell scripts (README's working scripts)."""
    rels = (
        "examples/toy_system/run.sh",
        "examples/methyltransferase/run_all.sh",
        "examples/methyltransferase/run_stepwise.sh",
    )
    return [REPO_ROOT / rel for rel in rels]


def bool_style_sources() -> list[Path]:
    """Authored surfaces whose command guidance must use canonical bool syntax."""
    paths: list[Path] = []
    for name in ("README.md", "CONTRIBUTING.md"):
        candidate = REPO_ROOT / name
        if candidate.exists():
            paths.append(candidate)
    for root_name in ("docs", "skills"):
        root = REPO_ROOT / root_name
        if root.exists():
            paths.extend(
                p
                for p in root.rglob("*.md")
                if "reference" not in p.relative_to(REPO_ROOT).parts
                and "_build" not in p.relative_to(REPO_ROOT).parts
            )
    examples = REPO_ROOT / "examples"
    if examples.exists():
        paths.extend(examples.rglob("*.md"))
        paths.extend(examples.rglob("*.sh"))
    smoke = REPO_ROOT / "tests" / "smoke"
    if smoke.exists():
        paths.extend(smoke.glob("*.sh"))
    return sorted(set(paths))


# --------------------------------------------------------------------------- #
# Markdown fenced-command and public-shell extraction
# --------------------------------------------------------------------------- #
def looks_like_data_literal(token: str) -> bool:
    """Return whether *token* is a Python list/tuple/number DATA literal.

    A quoted ``[(12,45,2.20)]`` or ``[-205.1, -190.4]`` is data passed to an
    option, not synopsis notation, and must never be stripped of brackets and
    mistaken for an option name.
    """
    stripped = token.strip()
    if not stripped or stripped[0] not in "[({'\"-0123456789.":
        return False
    try:
        ast.literal_eval(stripped)
    except (ValueError, SyntaxError, TypeError, MemoryError, RecursionError):
        return False
    return True


def _classify_executable(text: str) -> bool:
    """Execution eligibility for the dry-run smoke.

    A command is dry-run-eligible only when it carries no synopsis placeholder
    and no bracket notation. Quoted list/tuple DATA is retained for static
    validation but references specific atoms/indices that the synthetic dry-run
    fixture cannot satisfy, so such commands stay static-only.
    """
    return not any(mark in text for mark in _SYNOPSIS_MARKS)


def _iter_fenced_blocks(path: Path) -> Iterable[list[tuple[int, str]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    in_fence = False
    fence_lang = ""
    block: list[tuple[int, str]] = []
    for lineno, line in enumerate(lines, start=1):
        stripped = line.strip()
        if stripped.startswith("```"):
            marker = stripped[3:].strip().lower()
            if not in_fence:
                in_fence = True
                fence_lang = marker
                block = []
            else:
                if fence_lang in _CODE_LANGS:
                    yield block
                in_fence = False
                fence_lang = ""
                block = []
            continue
        if in_fence:
            block.append((lineno, line))
    # Unterminated fence: still surface any collected command lines.
    if in_fence and fence_lang in _CODE_LANGS and block:
        yield block


def _commands_from_block(block: list[tuple[int, str]]) -> list[tuple[int, str]]:
    commands: list[tuple[int, str]] = []
    current = ""
    current_line: int | None = None
    for lineno, raw in block:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if " #" in line:
            line = line.split(" #", 1)[0].rstrip()
            if not line:
                continue
        if line.startswith("$"):
            line = line[1:].strip()
        if not current:
            current_line = lineno
        current = f"{current} {line}".strip() if current else line
        if current.endswith("\\"):
            current = current[:-1].rstrip()
            continue
        commands.append((current_line if current_line is not None else lineno, current))
        current = ""
        current_line = None
    if current:
        commands.append((current_line if current_line is not None else block[-1][0], current))
    return commands


def _authored_from_pairs(path: Path, pairs: list[tuple[int, str]]) -> list[AuthoredCommand]:
    out: list[AuthoredCommand] = []
    for lineno, text in pairs:
        if not text.startswith(TOOL_NAME):
            continue
        out.append(AuthoredCommand(path, lineno, text, _classify_executable(text)))
    return out


def extract_markdown_commands(paths: Iterable[Path]) -> list[AuthoredCommand]:
    out: list[AuthoredCommand] = []
    for path in paths:
        for block in _iter_fenced_blocks(path):
            out.extend(_authored_from_pairs(path, _commands_from_block(block)))
    return out


def extract_docs_commands() -> list[AuthoredCommand]:
    return extract_markdown_commands(sorted(DOCS_ROOT.rglob("*.md")))


def extract_shell_commands(paths: Iterable[Path]) -> list[AuthoredCommand]:
    """Extract tool invocations from public shell scripts (whole file is shell)."""
    out: list[AuthoredCommand] = []
    for path in paths:
        block = list(enumerate(path.read_text(encoding="utf-8").splitlines(), start=1))
        out.extend(_authored_from_pairs(path, _commands_from_block(block)))
    return out


# --------------------------------------------------------------------------- #
# Live Click command/option discovery
# --------------------------------------------------------------------------- #
def _is_unavailable(command: click.Command | None) -> bool:
    if command is None:
        return True
    return bool((command.help or "").startswith("[Unavailable]"))


def subcommand_from_tokens(tokens: list[str]) -> str:
    if len(tokens) < 2 or tokens[1].startswith("-"):
        return "all"
    if tokens[1].startswith("<") or tokens[1].startswith("["):
        # Generic usage templates still participate in token validation; use
        # `all` as the representative command for shared help flags.
        return "all"
    return tokens[1]


def allowed_option_names(command: click.Command) -> set[str]:
    allowed = {
        opt
        for param in command.params
        if hasattr(param, "opts")
        for opt in (*param.opts, *getattr(param, "secondary_opts", ()))
    }
    allowed.update({"-h", "--help", "--help-advanced", "--version"})
    return allowed


def resolve_live_bool_options(root_cli: click.Group) -> frozenset[str]:
    """Union of every live value/toggle/single-flag boolean option name."""
    ctx = root_cli.make_context(TOOL_NAME, [], resilient_parsing=True)
    names: set[str] = set()
    try:
        for command_name in root_cli.list_commands(ctx):
            command = root_cli.get_command(ctx, command_name)
            if _is_unavailable(command):
                continue
            resolver = getattr(root_cli, "_resolve_bool_options", None)
            if resolver is None:
                continue
            value_opts, toggle_opts, _aliases, single_opts = resolver(ctx, command_name)
            names.update(value_opts)
            names.update(toggle_opts)
            names.update(single_opts)
    finally:
        ctx.close()
    return frozenset(name.lower() for name in names)


# --------------------------------------------------------------------------- #
# Authored option validation
# --------------------------------------------------------------------------- #
def validate_option_names(
    commands: Iterable[AuthoredCommand], root_cli: click.Group
) -> list[str]:
    """Validate every authored option token against the live command object.

    Executing a subcommand's ``--help`` only proves the command exists; Click's
    eager help returns before parsing invalid options. This static pass checks
    each option token against the live command and deliberately retains examples
    that contain ``<placeholder>`` or data-literal notation. An unavailable lazy
    command fails explicitly rather than disappearing from validation.
    """
    root_ctx = root_cli.make_context(TOOL_NAME, [], resilient_parsing=True)
    errors: list[str] = []
    try:
        for command in commands:
            try:
                tokens = shlex.split(command.text)
            except ValueError as exc:
                errors.append(
                    f"{command.location}: cannot parse shell words: {command.text!r}: {exc}"
                )
                continue
            if not tokens or tokens[0] != TOOL_NAME:
                continue

            subcmd = subcommand_from_tokens(tokens)
            live = root_cli.get_command(root_ctx, subcmd)
            if _is_unavailable(live):
                errors.append(
                    f"{command.location}: unavailable or unknown subcommand "
                    f"{subcmd!r}: {command.text}"
                )
                continue
            allowed = allowed_option_names(live)

            start = 2 if len(tokens) > 1 and tokens[1] == subcmd else 1
            for token in tokens[start:]:
                # A quoted list/tuple literal is DATA, not synopsis notation.
                if looks_like_data_literal(token):
                    continue
                token = token.lstrip("[").rstrip("] ,")
                if token in _SHELL_BREAKS or token == "--":
                    break
                if not token.startswith("-") or token == "-":
                    continue
                # Negative numeric option values are not option names.
                try:
                    float(token)
                    continue
                except ValueError:
                    pass
                # Synopsis blocks often compress aliases/pairs into one shell
                # word (`-b/--backend`, `--dump/--no-dump`, `--a|--b`); validate
                # each advertised spelling independently.
                names = [
                    part.split("=", 1)[0].removesuffix("...")
                    for part in re.split(r"[/|]", token)
                    if part.startswith("-")
                ]
                for name in names:
                    if name not in allowed:
                        errors.append(
                            f"{command.location}: unknown option {name!r} for "
                            f"{subcmd}: {command.text}"
                        )
    finally:
        root_ctx.close()
    return errors


# --------------------------------------------------------------------------- #
# Live-derived canonical-boolean style validation
# --------------------------------------------------------------------------- #
def validate_bool_style(
    paths: Iterable[Path], live_bool_options: frozenset[str]
) -> list[str]:
    """Reject legacy value-style booleans (``--flag True``) in authored guidance.

    Only option names that resolve to a live boolean option are flagged, so a
    value like ``--precision fp32`` is untouched. Comment lines (explanatory
    prose such as ``# --dump on freq``) are skipped; they are not command
    guidance.
    """
    errors: list[str] = []
    for path in paths:
        try:
            rel = path.relative_to(REPO_ROOT)
        except ValueError:
            rel = path
        if rel in _LEGACY_BOOL_PAGES:
            continue
        text = path.read_text(encoding="utf-8")
        for match in _BOOL_VALUE_RE.finditer(text):
            if match.group("option").lower() not in live_bool_options:
                continue
            line_start = text.rfind("\n", 0, match.start()) + 1
            line_end = text.find("\n", match.start())
            line_text = text[line_start : line_end if line_end != -1 else len(text)]
            if line_text.lstrip().startswith("#"):
                continue
            line = text.count("\n", 0, match.start()) + 1
            errors.append(
                f"{rel}:{line}: use --flag / --no-flag; found {match.group(0)!r}"
            )
    return errors
