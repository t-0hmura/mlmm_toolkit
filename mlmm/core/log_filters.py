"""Hide two dependency log messages below ``-v 3``.

FAIR-Chem logs ``Running (<commit>) with inference settings: ...`` on the root
logger each time a UMA predict unit initializes, and choreographer (the Chrome
control used by Kaleido) logs ``Resorting to unclean kill browser.`` when Chrome
does not close within 3 s, after the PNG is already saved. Both are kept at
``-v 3``; no other record is dropped.
"""

from __future__ import annotations

import logging
import re
import sys

_FAIRCHEM_SETTINGS_RE = re.compile(r"^Running \(.*\) with inference settings: ")
_CHOREOGRAPHER_LOGGER = "choreographer.browser_async"
_CHOREOGRAPHER_KILL_MSG = "Resorting to unclean kill browser."


def cli_verbose_level() -> int:
    """Return the CLI ``-v`` level, or 0 when it was never set in this process.

    The level lives in ``mlmm.core.utils``, which imports torch. It is read
    from ``sys.modules`` so that the PNG worker process, which never imports
    it, does not pay for that import.
    """
    utils = sys.modules.get("mlmm.core.utils")
    return utils.verbose_level() if utils is not None else 0


class FairchemInferenceSettingsFilter(logging.Filter):
    """Drop FAIR-Chem's inference-settings record below ``-v 3``."""

    def filter(self, record: logging.LogRecord) -> bool:
        if cli_verbose_level() >= 3:
            return True
        return not (
            record.name == "root"
            and _FAIRCHEM_SETTINGS_RE.match(record.getMessage())
        )


class ChoreographerUncleanKillFilter(logging.Filter):
    """Drop choreographer's unclean-kill record below ``-v 3``."""

    def filter(self, record: logging.LogRecord) -> bool:
        if cli_verbose_level() >= 3:
            return True
        return record.getMessage() != _CHOREOGRAPHER_KILL_MSG


_FAIRCHEM_FILTER = FairchemInferenceSettingsFilter()
_CHOREOGRAPHER_FILTER = ChoreographerUncleanKillFilter()


def install_dependency_log_filters() -> None:
    """Attach both filters to their loggers. Idempotent."""
    root = logging.getLogger()
    if _FAIRCHEM_FILTER not in root.filters:
        root.addFilter(_FAIRCHEM_FILTER)
    # logistro sets the logger class choreographer's loggers need, so the
    # logger is created through it; without logistro there is no choreographer.
    try:
        import logistro
    except ImportError:
        return
    logger = logistro.getLogger(_CHOREOGRAPHER_LOGGER)
    if _CHOREOGRAPHER_FILTER not in logger.filters:
        logger.addFilter(_CHOREOGRAPHER_FILTER)
