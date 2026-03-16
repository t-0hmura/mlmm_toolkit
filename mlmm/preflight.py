"""Shared preflight checks for CLI commands."""

from __future__ import annotations

import importlib
import shutil
from pathlib import Path
from typing import Iterable, Mapping

import click


def validate_existing_files(
    raw_paths: Iterable[str | Path],
    *,
    option_name: str,
    hint: str | None = None,
) -> list[Path]:
    """Validate that each path exists and is a file.

    Parameters
    ----------
    raw_paths
        Candidate file paths from CLI arguments.
    option_name
        Option label for error messages (for example ``-i/--input``).
    hint
        Optional extra text appended to the BadParameter message.
    """
    out: list[Path] = []
    for raw in raw_paths:
        p = Path(raw)
        if (not p.exists()) or p.is_dir():
            msg = f"{option_name} path '{raw}' not found or is a directory."
            if hint:
                msg = f"{msg} {hint}"
            raise click.BadParameter(msg)
        out.append(p)
    return out


def ensure_commands_available(commands: Iterable[str], *, context: str) -> None:
    """Raise ClickException when required executables are missing from PATH."""
    missing = [cmd for cmd in commands if shutil.which(cmd) is None]
    if not missing:
        return

    joined = ", ".join(missing)
    raise click.ClickException(
        f"[preflight] Missing required command(s) for {context}: {joined}. "
        "Please ensure they are available on PATH."
    )


def ensure_python_modules_available(
    modules: Mapping[str, str] | Iterable[str],
    *,
    context: str,
) -> None:
    """Raise ClickException when required Python modules are unavailable.

    ``modules`` can be either:
    - ``{\"module\": \"install hint\"}``
    - ``[\"module1\", \"module2\"]``
    """
    if isinstance(modules, Mapping):
        items = list(modules.items())
    else:
        items = [(name, "") for name in modules]

    missing: list[str] = []
    for module_name, _ in items:
        try:
            importlib.import_module(module_name)
        except (ImportError, ModuleNotFoundError):
            missing.append(module_name)

    if not missing:
        return

    hints = []
    for module_name, hint in items:
        if module_name in missing and hint:
            hints.append(f"{module_name}: {hint}")
    hint_text = f" Hints: {'; '.join(hints)}." if hints else ""
    raise click.ClickException(
        f"[preflight] Missing required Python module(s) for {context}: {', '.join(missing)}.{hint_text}"
    )
