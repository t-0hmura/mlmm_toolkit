from __future__ import annotations

import os
import shutil
import sys
import importlib.util
from importlib import import_module
from pathlib import Path
from typing import Any, Dict, Optional

_EXT_CACHE: Dict[str, Optional[Any]] = {
    "nonbonded": None,
    "analytical_hessian": None,
    "bonded": None,
}
_EXT_ERROR: Dict[str, Optional[str]] = {
    "nonbonded": None,
    "analytical_hessian": None,
    "bonded": None,
}


def _rebuild_hint() -> str:
    return (
        "To rebuild hessian_ff native extensions in this environment:\n"
        "  conda install -c conda-forge ninja -y\n"
        "  cd $(python -c \"import hessian_ff; print(hessian_ff.__path__[0])\")/native && make clean && make"
    )


def _with_rebuild_hint(msg: str) -> str:
    txt = str(msg).strip()
    if not txt:
        txt = "native extension unavailable"
    return f"{txt}\n{_rebuild_hint()}"


def _ensure_max_jobs() -> None:
    """Set Ninja parallel compile jobs if not provided by user.

    This reduces first-build wall time for native extensions.
    """
    if "MAX_JOBS" in os.environ:
        return
    ncpu = os.cpu_count()
    if ncpu is None or ncpu < 1:
        return
    os.environ["MAX_JOBS"] = str(int(ncpu))


def try_load_native_backend(module_name: str = "hessian_ff_native") -> Optional[Any]:
    """Try importing a native extension backend module.

    Returns the module object if available, otherwise ``None``.
    """
    try:
        return import_module(module_name)
    except Exception:
        return None


def native_backend_status(module_name: str = "hessian_ff_native") -> Dict[str, str]:
    """Return a short status dict for native backend availability."""
    mod = try_load_native_backend(module_name=module_name)
    if mod is None:
        return {
            "available": "false",
            "module": module_name,
            "backend": "native_required",
            "note": "native extension module is not loaded",
        }
    return {
        "available": "true",
        "module": module_name,
        "backend": "native",
        "note": "native extension is loaded",
    }


def _cache_build_dir(build_subdir: str) -> Path:
    """Return the user-cache fallback build directory.

    Used when the package-internal directory (site-packages) is read-only.
    Location: $XDG_CACHE_HOME/mlmm-toolkit/hessian_ff/<build_subdir>
    """
    cache_root = Path(
        os.environ.get("XDG_CACHE_HOME", os.path.expanduser("~/.cache"))
    )
    return cache_root / "mlmm-toolkit" / "hessian_ff" / build_subdir


def _build_in_tree_extension(
    *,
    key: str,
    ext_name: str,
    source_files: list[str],
    build_subdir: str,
    verbose: bool,
    force_rebuild: bool,
) -> Optional[Any]:
    if _EXT_CACHE[key] is not None and not force_rebuild:
        return _EXT_CACHE[key]
    if _EXT_ERROR[key] is not None and not force_rebuild:
        return None

    here = Path(__file__).resolve().parent
    # Candidate build directories: package-internal first, then user cache.
    pkg_build_dir = here / build_subdir
    cache_build_dir = _cache_build_dir(build_subdir)
    build_dirs = [pkg_build_dir, cache_build_dir]

    def _find_prebuilt_so() -> Optional[Path]:
        """Search all candidate directories for a prebuilt .so."""
        for bd in build_dirs:
            if not bd.is_dir():
                continue
            cands = sorted(
                bd.glob(f"{ext_name}*.so"),
                key=lambda p: p.stat().st_mtime,
                reverse=True,
            )
            if cands:
                return cands[0]
        return None

    def _load_prebuilt_so(path: Path) -> Optional[Any]:
        try:
            mod_name = path.stem
            if mod_name in sys.modules:
                mod = sys.modules[mod_name]
                _EXT_CACHE[key] = mod
                _EXT_ERROR[key] = None
                return mod
            spec = importlib.util.spec_from_file_location(mod_name, str(path))
            if spec is None or spec.loader is None:
                raise ImportError(f"spec loader is unavailable for {path}")
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            _EXT_CACHE[key] = mod
            _EXT_ERROR[key] = None
            return mod
        except Exception as e:
            _EXT_ERROR[key] = _with_rebuild_hint(
                f"failed to load prebuilt extension {path}: {e}"
            )
            return None

    # Prefer prebuilt artifact when available (useful on nodes without a compiler toolchain).
    if not force_rebuild:
        prebuilt = _find_prebuilt_so()
        if prebuilt is not None:
            loaded = _load_prebuilt_so(prebuilt)
            if loaded is not None:
                return loaded

    try:
        from torch.utils.cpp_extension import load
    except Exception as e:
        _EXT_ERROR[key] = _with_rebuild_hint(
            f"torch cpp_extension import failed: {e}"
        )
        return None

    _ensure_max_jobs()

    srcs = [here / s for s in source_files]

    def _try_build_in_dir(
        build_dir: Path,
        extra_cflags: list[str],
        extra_ldflags: Optional[list[str]] = None,
    ):
        os.makedirs(build_dir, exist_ok=True)
        kwargs = dict(
            name=ext_name,
            sources=[str(s) for s in srcs],
            extra_cflags=extra_cflags,
            extra_ldflags=(extra_ldflags or []),
            build_directory=str(build_dir),
            verbose=bool(verbose),
        )
        # Some cluster conda envs do not provide ninja. Prefer distutils build
        # in that case instead of failing hard.
        if shutil.which("ninja") is None:
            kwargs["use_ninja"] = False
        try:
            return load(**kwargs)
        except TypeError:
            kwargs.pop("use_ninja", None)
            return load(**kwargs)

    # Try aggressive CPU flags first, then fall back progressively.
    build_attempts = [
        (["-Ofast", "-ffast-math", "-funroll-loops", "-march=native", "-mtune=native", "-fopenmp"], ["-fopenmp"]),
        (["-O3", "-ffast-math", "-funroll-loops", "-march=native", "-mtune=native", "-fopenmp"], ["-fopenmp"]),
        (["-O3", "-ffast-math", "-funroll-loops", "-fopenmp"], ["-fopenmp"]),
        (["-O3", "-fopenmp"], ["-fopenmp"]),
        (["-O3"], []),
    ]

    # Try each build directory: package-internal first, then user cache.
    last_err: Optional[Exception] = None
    for build_dir in build_dirs:
        try:
            os.makedirs(build_dir, exist_ok=True)
        except OSError:
            continue
        for cflags, ldflags in build_attempts:
            try:
                _EXT_CACHE[key] = _try_build_in_dir(
                    build_dir,
                    extra_cflags=cflags,
                    extra_ldflags=ldflags,
                )
                _EXT_ERROR[key] = None
                return _EXT_CACHE[key]
            except Exception as e:
                last_err = e
                continue
        # All flag combinations failed in this directory; try next.

    _EXT_CACHE[key] = None
    _EXT_ERROR[key] = _with_rebuild_hint(
        f"hessian_ff build attempts failed: {last_err}"
    )
    return None


def get_nonbonded_extension(
    *,
    verbose: bool = False,
    force_rebuild: bool = False,
) -> Optional[Any]:
    """Build/load in-tree C++ extension for nonbonded kernels.

    References:
    - torch.utils.cpp_extension.load() runtime build workflow.
    """
    return _build_in_tree_extension(
        key="nonbonded",
        ext_name="hessian_ff_nonbonded_ext",
        source_files=["nonbonded_ext.cpp"],
        build_subdir=".build_nonbonded",
        verbose=bool(verbose),
        force_rebuild=bool(force_rebuild),
    )


def nonbonded_extension_status() -> Dict[str, str]:
    ext = get_nonbonded_extension(verbose=False, force_rebuild=False)
    if ext is None:
        note = _EXT_ERROR["nonbonded"] or _with_rebuild_hint("extension unavailable")
        return {
            "available": "false",
            "backend": "native_required",
            "note": note,
        }
    return {
        "available": "true",
        "backend": "native_nonbonded_cpp",
        "note": "nonbonded extension loaded",
    }


def get_analytical_hessian_extension(
    *,
    verbose: bool = False,
    force_rebuild: bool = False,
) -> Optional[Any]:
    """Build/load in-tree C++ extension for analytical Hessian helpers."""
    return _build_in_tree_extension(
        key="analytical_hessian",
        ext_name="hessian_ff_analytical_hessian_ext",
        source_files=["analytical_hessian_ext.cpp"],
        build_subdir=".build_analytical_hessian",
        verbose=bool(verbose),
        force_rebuild=bool(force_rebuild),
    )


def analytical_hessian_extension_status() -> Dict[str, str]:
    ext = get_analytical_hessian_extension(verbose=False, force_rebuild=False)
    if ext is None:
        note = _EXT_ERROR["analytical_hessian"] or _with_rebuild_hint("extension unavailable")
        return {
            "available": "false",
            "backend": "native_analytical_hessian_optional",
            "note": note,
        }
    return {
        "available": "true",
        "backend": "native_analytical_hessian_cpp",
        "note": "analytical Hessian extension loaded",
    }


def get_bonded_extension(
    *,
    verbose: bool = False,
    force_rebuild: bool = False,
) -> Optional[Any]:
    """Build/load in-tree C++ extension for bonded energy-force kernels."""
    return _build_in_tree_extension(
        key="bonded",
        ext_name="hessian_ff_bonded_ext",
        source_files=["bonded_ext.cpp"],
        build_subdir=".build_bonded",
        verbose=bool(verbose),
        force_rebuild=bool(force_rebuild),
    )


def bonded_extension_status() -> Dict[str, str]:
    ext = get_bonded_extension(verbose=False, force_rebuild=False)
    if ext is None:
        note = _EXT_ERROR["bonded"] or _with_rebuild_hint("extension unavailable")
        return {
            "available": "false",
            "backend": "native_bonded_optional",
            "note": note,
        }
    return {
        "available": "true",
        "backend": "native_bonded_cpp",
        "note": "bonded extension loaded",
    }


def build_native_extensions(
    *,
    verbose: bool = False,
    force_rebuild: bool = False,
) -> Dict[str, str]:
    """Build/load all native extensions up front.

    This provides a practical "compile together" workflow by triggering
    all in-tree extension builds in one step before production runs.
    """
    ext_nb = get_nonbonded_extension(verbose=verbose, force_rebuild=force_rebuild)
    ext_ah = get_analytical_hessian_extension(verbose=verbose, force_rebuild=force_rebuild)
    ext_bd = get_bonded_extension(verbose=verbose, force_rebuild=force_rebuild)
    return {
        "nonbonded": "ok" if ext_nb is not None else f"error: {_EXT_ERROR['nonbonded']}",
        "analytical_hessian": "ok"
        if ext_ah is not None
        else f"error: {_EXT_ERROR['analytical_hessian']}",
        "bonded": "ok" if ext_bd is not None else f"error: {_EXT_ERROR['bonded']}",
    }
