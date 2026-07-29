# mlmm/__init__.py

try:
    from mlmm._version import __version__, __version_tuple__
except ImportError:
    __version__ = "0.0.0.dev0"
    __version_tuple__ = (0, 0, 0, "dev0")

__all__ = [
    "__version__",
    "__version_tuple__",
    "MLMMCore",
    "MLMMASECalculator",
    "mlmm",
    "mlmm_ase",
    "mlmm_mm_only",
]

# PEP 562 lazy resolution avoids loading torch and AmberTools wrappers when
# callers only need lightweight package metadata.
_LAZY_IMPORTS = {
    "MLMMCore": "mlmm.backends.mlmm_calc",
    "MLMMASECalculator": "mlmm.backends.mlmm_calc",
    "mlmm": "mlmm.backends.mlmm_calc",
    "mlmm_ase": "mlmm.backends.mlmm_calc",
    "mlmm_mm_only": "mlmm.backends.mlmm_calc",
}


def __getattr__(name: str):
    if name in _LAZY_IMPORTS:
        import importlib

        module = importlib.import_module(_LAZY_IMPORTS[name])
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module 'mlmm' has no attribute {name!r}")
