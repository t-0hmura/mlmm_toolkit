# mlmm_toolkit/__init__.py

from pathlib import Path
import sys

# Prefer vendored sibling packages (e.g., pysisyphus/) over unrelated
# site-packages installations when running from an editable checkout.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

try:
    from mlmm_toolkit._version import __version__, __version_tuple__
except ImportError:
    __version__ = "0.0.0.dev0"
    __version_tuple__ = (0, 0, 0, "dev0")

__all__ = [
    "__version__",
    "__version_tuple__",
]
