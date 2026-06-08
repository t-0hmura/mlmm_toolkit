#!/usr/bin/env python3
"""Static drift check for the primary `--help` option registry.

Validates ``mlmm.cli.app._SUBCOMMAND_PRIMARY_HELP_OPTIONS`` against the actual
Click commands that back each subcommand:

* Every entry must include ``"--help-advanced"`` so the discoverability flag
  itself is not hidden when the registry filters non-primary options.
* Every other token in each frozenset must correspond to an actual
  ``@click.option`` opt string on the resolved Click command — long opts via
  ``Option.opts`` and short opts via ``Option.opts`` / ``Option.secondary_opts``.
  Orphan entries (dead aliases / removed flags) fail this check.

Run from the repository root::

    python .github/scripts/check_help_registry.py

Exit 0 = registry consistent, 1 = drift (printed to stderr).
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

# ``--help-advanced`` is injected by ``cli.help_pages._ensure_help_advanced_option``
# so it is NOT declared on the lazy Click commands themselves. Treat it as the
# single allow-listed token that need not resolve to a ``command.params`` opt.
_INJECTED_OPTS: frozenset[str] = frozenset({"--help-advanced"})

from mlmm.cli.app import (  # noqa: E402
    _LAZY_SUBCOMMANDS,
    _SUBCOMMAND_PRIMARY_HELP_OPTIONS,
)


def _resolve_command(spec: tuple[str, str, str]):
    module_name, attr, _help = spec
    mod = importlib.import_module(module_name)
    return getattr(mod, attr)


def _declared_opts(command) -> set[str]:
    opts: set[str] = set()
    for param in getattr(command, "params", ()):  # click.Option / click.Argument
        opts.update(getattr(param, "opts", ()) or ())
        opts.update(getattr(param, "secondary_opts", ()) or ())
    return opts


def main() -> int:
    errors: list[str] = []

    for subcmd, primary in _SUBCOMMAND_PRIMARY_HELP_OPTIONS.items():
        if "--help-advanced" not in primary:
            errors.append(
                f"[{subcmd}] missing '--help-advanced' in primary set "
                f"(would self-hide the discoverability flag)"
            )

        spec = _LAZY_SUBCOMMANDS.get(subcmd)
        if spec is None:
            errors.append(
                f"[{subcmd}] registry entry has no corresponding "
                f"_LAZY_SUBCOMMANDS spec"
            )
            continue

        try:
            command = _resolve_command(spec)
        except Exception as exc:  # pragma: no cover - import failure surfaces here
            errors.append(f"[{subcmd}] could not import {spec[0]}.{spec[1]}: {exc}")
            continue

        declared = _declared_opts(command)
        for token in sorted(primary - _INJECTED_OPTS):
            if token not in declared:
                errors.append(
                    f"[{subcmd}] primary token {token!r} is not declared by "
                    f"any @click.option on {spec[0]}.{spec[1]} "
                    f"(orphan / dead registry entry)"
                )

    if errors:
        print("help registry drift detected:", file=sys.stderr)
        for line in errors:
            print(f"  - {line}", file=sys.stderr)
        return 1

    print("help registry OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
