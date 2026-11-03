"""Skip a test that needs a repository file the source distribution leaves out.

The sdist ships ``tests/`` but not ``.github/``, ``tests/smoke/`` or
``hessian_ff/tests/data/`` (see ``MANIFEST.in``), so a test that reads one of
them skips there and runs in a checkout.
"""

from __future__ import annotations

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


def require_repo_files(*relative: str) -> None:
    """Skip the calling test, or the module when called at import, if a file is missing.

    ``relative`` are paths from the repository root.
    """
    missing = [path for path in relative if not (REPO_ROOT / path).exists()]
    if missing:
        pytest.skip(
            f"needs {', '.join(missing)}, which this source tree does not have "
            "(the sdist does not ship it)",
            allow_module_level=True,
        )
