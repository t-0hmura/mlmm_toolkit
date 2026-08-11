"""Shared preflight checks for CLI commands."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

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
