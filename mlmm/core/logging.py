"""Shared stdlib logging configurator for mlmm subcommands.

The per-subcommand ``-v/--verbose LEVEL`` callback enables DEBUG logging at
level 3. Direct callers may pass 0, 1, or 2+ to select WARNING, INFO, or DEBUG.
Existing ``click.echo()`` calls are unaffected.
"""
from __future__ import annotations

import logging
import sys

_DEFAULT_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
_VERBOSE_TO_LEVEL = {
    0: logging.WARNING,
    1: logging.INFO,
    2: logging.DEBUG,
}


def setup_logging(verbose: int = 0) -> None:
    """Configure stdlib logging from a numeric verbosity level. Idempotent."""
    level = _VERBOSE_TO_LEVEL.get(min(verbose, 2), logging.DEBUG)
    root = logging.getLogger()
    for h in list(root.handlers):
        root.removeHandler(h)
    handler = logging.StreamHandler(stream=sys.stderr)
    handler.setFormatter(logging.Formatter(_DEFAULT_FORMAT))
    handler.setLevel(level)
    root.addHandler(handler)
    root.setLevel(level)
