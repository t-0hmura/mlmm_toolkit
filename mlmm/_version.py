"""Single source of truth for the mlmm-toolkit release version."""
from __future__ import annotations

__all__ = [
    "__version__",
    "__version_tuple__",
    "version",
    "version_tuple",
    "__commit_id__",
    "commit_id",
]

version: str
__version__: str
__version_tuple__: tuple[int | str, ...]
version_tuple: tuple[int | str, ...]
commit_id: str | None
__commit_id__: str | None

__version__ = version = "0.3.7"
__version_tuple__ = version_tuple = (0, 3, 7)

__commit_id__ = commit_id = None
