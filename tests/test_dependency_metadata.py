"""Version metadata of the installed package.

The dependency floors, extras and Python range in ``pyproject.toml`` are checked by
``.github/scripts/check_release_versions.py``.
"""

from __future__ import annotations


def test_version_tuple_matches_the_packaged_version() -> None:
    from packaging.version import Version
    from mlmm import _version

    assert _version.__version_tuple__ == _version.version_tuple == Version(_version.__version__).release
    assert _version.version == _version.__version__
