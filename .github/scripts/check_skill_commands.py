#!/usr/bin/env python3
"""Audit fenced ``bash`` blocks in skills/**/*.md.

For every line that starts with ``mlmm <subcommand>``, parse out flag
tokens (``--foo`` and ``-f``) and verify each one is registered on that
subcommand via Click introspection. Stale flags from past edits are
reported with file:line.

This is intentionally conservative:
* shell line continuations (\\) are joined.
* placeholders / templates (``<arg>``, ``{xyz,pdb,gjf}``) are skipped.
* values for known boolean flags (``--tsopt true``, ``--no-tsopt``) are
  accepted as flag-only.
* free-form prose example fragments inside backticks (single-line
  ``mlmm extract --foo``) are also checked.

Exits non-zero on any unknown flag.
"""

from __future__ import annotations

import re
import shlex
import sys
from pathlib import Path

import click

REPO_ROOT = Path(__file__).resolve().parents[2]
SKILLS_DIR = REPO_ROOT / "skills"

sys.path.insert(0, str(REPO_ROOT))
from mlmm.cli import cli as root_cli  # noqa: E402


def _collect_subcommand_flags() -> dict[str, set[str]]:
    flags_per_cmd: dict[str, set[str]] = {}
    ctx = click.Context(root_cli)
    for name in root_cli.list_commands(ctx):
        cmd = root_cli.get_command(ctx, name)
        if cmd is None:
            continue
        flags: set[str] = set()
        for p in cmd.params:
            for opt in getattr(p, "opts", []) or []:
                flags.add(opt)
            for opt in getattr(p, "secondary_opts", []) or []:
                flags.add(opt)
        # universal flags every subcommand inherits in our setup
        flags.update({"--help", "--help-advanced"})
        flags_per_cmd[name] = flags
    return flags_per_cmd


_FENCE_RE = re.compile(r"^```(?:bash|console|sh)?\s*$")
_FENCE_END_RE = re.compile(r"^```\s*$")
_INLINE_RE = re.compile(r"`(mlmm\s+[^`]+)`")
_FLAG_RE = re.compile(r"^-{1,2}[A-Za-z][\w-]*$")


def _iter_command_lines(text: str):
    in_fence = False
    pending: list[str] = []
    pending_lineno = 0
    for lineno, line in enumerate(text.splitlines(), start=1):
        if not in_fence and _FENCE_RE.match(line):
            in_fence = True
            continue
        if in_fence and _FENCE_END_RE.match(line):
            in_fence = False
            if pending and "mlmm" in " ".join(pending):
                yield pending_lineno, " ".join(pending)
            pending = []
            continue
        if in_fence:
            stripped = line.rstrip()
            if stripped.endswith("\\"):
                if not pending:
                    pending_lineno = lineno
                pending.append(stripped[:-1].strip())
                continue
            # End of a logical command
            if pending:
                pending.append(stripped.strip())
                if "mlmm" in " ".join(pending):
                    yield pending_lineno, " ".join(pending)
                pending = []
            else:
                # standalone single-line command
                if "mlmm" in stripped:
                    yield lineno, stripped.strip()
            continue
        for m in _INLINE_RE.finditer(line):
            yield lineno, m.group(1)


def _check_command(cmd_text: str, flags_per_cmd: dict[str, set[str]]):
    try:
        tokens = shlex.split(cmd_text, posix=True)
    except ValueError:
        return []
    if len(tokens) < 2 or tokens[0] != "mlmm":
        return []
    sub = tokens[1]
    if sub not in flags_per_cmd:
        return []
    valid = flags_per_cmd[sub]
    bad: list[str] = []
    for tok in tokens[2:]:
        if tok.startswith("<") or tok.startswith("{") or tok.startswith("["):
            continue
        if not tok.startswith("-"):
            continue
        # split --foo=bar
        flag = tok.split("=", 1)[0]
        if not _FLAG_RE.match(flag):
            continue
        if flag not in valid:
            bad.append(flag)
    return bad


def main() -> int:
    flags_per_cmd = _collect_subcommand_flags()
    n_files = 0
    n_errors = 0
    for path in sorted(SKILLS_DIR.rglob("*.md")):
        text = path.read_text()
        for lineno, cmd in _iter_command_lines(text):
            bad = _check_command(cmd, flags_per_cmd)
            if bad:
                n_errors += 1
                print(f"{path.relative_to(REPO_ROOT)}:{lineno}: unknown flag(s) {bad} in: {cmd[:120]}")
        n_files += 1
    print(f"\nChecked {n_files} skill files; {n_errors} stale flag occurrences.")
    return 1 if n_errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
