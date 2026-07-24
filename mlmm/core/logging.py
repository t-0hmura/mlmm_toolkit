"""Shared stdlib logging configurator for mlmm subcommands.

Every subcommand accepts ``-v/--verbose LEVEL`` with integer levels 0--3
(default 2). Console narration is gated separately; level 3 additionally calls
``setup_logging(3)`` so module DEBUG records are visible. Existing
``click.echo()`` calls are unaffected. Pysisyphus / fairchem handler
suppression remains owned by ``DefaultGroup``.
"""
from __future__ import annotations

import logging
import sys

_DEFAULT_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
_VERBOSE_TO_LEVEL = {
    0: logging.CRITICAL,
    1: logging.WARNING,
    2: logging.INFO,
    3: logging.DEBUG,
}


def setup_logging(verbose: int = 0) -> None:
    """Configure stdlib logging for a unified verbosity level from 0 to 3.

    Idempotent — re-calling overwrites the root handler config.
    Existing pysisyphus/fairchem handler suppression remains intact
    (those use their own logger names and are silenced separately
    via DefaultGroup._silence_pysisyphus_loggers).
    """
    normalized = max(0, min(int(verbose), 3))
    level = _VERBOSE_TO_LEVEL[normalized]
    root = logging.getLogger()
    # remove existing handlers to avoid duplicate output on re-config
    for h in list(root.handlers):
        root.removeHandler(h)
    handler = logging.StreamHandler(stream=sys.stderr)
    handler.setFormatter(logging.Formatter(_DEFAULT_FORMAT))
    handler.setLevel(level)
    root.addHandler(handler)
    root.setLevel(level)
