#!/usr/bin/env python3
"""Audit fenced ``bash`` blocks in skills/**/*.md.

For every line that starts with ``mlmm <subcommand>``, parse out flag
tokens (``--foo`` and ``-f``), verify each one against Click introspection,
and reject examples that omit the topology needed to execute them. XYZ
examples must also supply the matching PDB topology reference.

This is intentionally conservative:
* shell line continuations (\\) are joined.
* placeholders / templates (``<arg>``, ``{xyz,pdb,gjf}``) are not treated
  as concrete XYZ paths.
* values for known boolean flags (``--tsopt true``, ``--no-tsopt``) are
  accepted as flag-only.
* free-form prose example fragments inside backticks (single-line
  ``mlmm extract --foo``) are also checked.

Exits non-zero on an unknown flag or an incomplete runnable example.
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


class CommandContract:
    def __init__(
        self,
        *,
        flags: set[str],
        parm_flags: set[str],
        charge_flags: set[str],
    ) -> None:
        self.flags = frozenset(flags)
        self.parm_flags = frozenset(parm_flags)
        self.charge_flags = frozenset(charge_flags)


def _collect_subcommand_contracts() -> dict[str, CommandContract]:
    contracts: dict[str, CommandContract] = {}
    ctx = click.Context(root_cli)
    for name in root_cli.list_commands(ctx):
        cmd = root_cli.get_command(ctx, name)
        if cmd is None:
            continue
        flags: set[str] = set()
        parm_flags: set[str] = set()
        charge_flags: set[str] = set()
        for p in cmd.params:
            opts = set(getattr(p, "opts", []) or [])
            opts.update(getattr(p, "secondary_opts", []) or [])
            for opt in opts:
                flags.add(opt)
            if getattr(p, "name", None) == "real_parm7" and p.required:
                parm_flags.update(opts)
            if getattr(p, "name", None) in {
                "charge",
                "charge_override",
                "ligand_charge",
            }:
                charge_flags.update(opts)
        # universal flags every subcommand inherits in our setup
        flags.update({"--help", "--help-advanced"})
        contracts[name] = CommandContract(
            flags=flags,
            parm_flags=parm_flags,
            charge_flags=charge_flags,
        )
    return contracts


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


def _option_values(tokens: list[str], option_names: set[str]) -> list[str]:
    """Return values following any named option until the next option token."""

    values: list[str] = []
    collecting = False
    for tok in tokens:
        flag = tok.split("=", 1)[0]
        if flag in option_names:
            collecting = True
            if "=" in tok:
                values.append(tok.split("=", 1)[1])
                collecting = False
            continue
        if tok.startswith("-"):
            collecting = False
            continue
        if collecting:
            values.append(tok)
    return values


def _is_concrete_xyz(value: str) -> bool:
    return not any(mark in value for mark in "<>{}[]") and value.lower().endswith(".xyz")


def _check_command(cmd_text: str, contracts: dict[str, CommandContract]) -> list[str]:
    try:
        tokens = shlex.split(cmd_text, posix=True)
    except ValueError:
        return []
    if len(tokens) < 2 or tokens[0] != "mlmm":
        return []
    sub = tokens[1]
    if sub not in contracts:
        return []
    contract = contracts[sub]
    valid = contract.flags
    issues: list[str] = []
    present_flags: set[str] = set()
    for tok in tokens[2:]:
        if tok.startswith("<") or tok.startswith("{") or tok.startswith("["):
            continue
        if not tok.startswith("-"):
            continue
        # split --foo=bar
        flag = tok.split("=", 1)[0]
        if not _FLAG_RE.match(flag):
            continue
        present_flags.add(flag)
        if flag not in valid:
            issues.append(f"unknown flag {flag}")

    if present_flags & {"--help", "--help-advanced"}:
        return issues

    inputs = _option_values(tokens[2:], {"-i", "--input"})
    # Bare command names and partial prose snippets are references, not runnable
    # examples. Unknown flags above are still checked in those fragments.
    if not inputs:
        return issues

    # Synopsis/template lines are flag inventories, not copy-paste commands.
    # Unknown flags are still checked above; completeness applies only to
    # concrete examples.
    template_tokens = any(
        tok.startswith("[-") or tok in {"[", "]", "..."}
        for tok in tokens[2:]
    )
    template_inputs = any(
        value == "..." or any(mark in value for mark in "<>{}")
        for value in inputs
    )
    if template_tokens or template_inputs:
        return issues

    plot_only_scan3d = sub == "scan3d" and "--csv" in present_flags
    parm_flags = set(contract.parm_flags)
    if sub in {"irc", "scan3d"} and not plot_only_scan3d:
        parm_flags.add("--parm")
    if parm_flags and not (present_flags & parm_flags):
        # ``irc`` may obtain calc.real_parm7 from YAML. Click-required topology
        # options on the other commands cannot be satisfied this way.
        if not (sub == "irc" and "--config" in present_flags):
            issues.append(f"missing topology option ({'/'.join(sorted(parm_flags))})")

    charge_commands = {
        "all", "dft", "freq", "irc", "opt", "path-opt", "path-search",
        "scan", "scan2d", "scan3d", "sp", "tsopt",
    }
    gjf_supplies_charge = bool(inputs) and all(
        not any(mark in value for mark in "<>{}[]")
        and value.lower().endswith(".gjf")
        for value in inputs
    )
    if (
        sub in charge_commands
        and not plot_only_scan3d
        and not gjf_supplies_charge
        and not (present_flags & set(contract.charge_flags))
        and "--config" not in present_flags
    ):
        expected = "/".join(sorted(contract.charge_flags)) or "-q/--charge"
        issues.append(f"missing charge option ({expected})")

    xyz_ref_commands = {
        "all",
        "dft",
        "freq",
        "irc",
        "opt",
        "path-opt",
        "path-search",
        "scan",
        "scan2d",
        "scan3d",
        "sp",
        "tsopt",
    }
    xyz_positions = (
        [i for i, value in enumerate(inputs) if _is_concrete_xyz(value)]
        if sub in xyz_ref_commands
        else []
    )
    if xyz_positions:
        refs = _option_values(tokens[2:], {"--ref-pdb"})
        required_refs = max(xyz_positions) + 1 if sub in {"path-opt", "path-search"} else 1
        if len(refs) < required_refs:
            issues.append(
                f"XYZ input requires {required_refs} corresponding --ref-pdb value(s)"
            )
    return issues


def main() -> int:
    contracts = _collect_subcommand_contracts()
    n_files = 0
    n_errors = 0
    for path in sorted(SKILLS_DIR.rglob("*.md")):
        text = path.read_text()
        for lineno, cmd in _iter_command_lines(text):
            issues = _check_command(cmd, contracts)
            if issues:
                n_errors += 1
                print(
                    f"{path.relative_to(REPO_ROOT)}:{lineno}: "
                    f"{'; '.join(issues)} in: {cmd[:120]}"
                )
        n_files += 1
    print(f"\nChecked {n_files} skill files; {n_errors} command-contract errors.")
    return 1 if n_errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
