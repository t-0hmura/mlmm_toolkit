from __future__ import annotations

import hashlib
import importlib.util
from importlib import import_module
import json
import os
import platform
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import sysconfig
import tempfile
from typing import Any, Dict, Mapping, Optional, Sequence

_FINGERPRINT_SCHEMA = 1
_FINGERPRINT_SUFFIX_LENGTH = 16
_BUILD_RECIPES = (
    (
        (
            "-Ofast",
            "-ffast-math",
            "-funroll-loops",
            "-march=native",
            "-mtune=native",
            "-fopenmp",
        ),
        ("-fopenmp",),
    ),
    (
        (
            "-O3",
            "-ffast-math",
            "-funroll-loops",
            "-march=native",
            "-mtune=native",
            "-fopenmp",
        ),
        ("-fopenmp",),
    ),
    (
        ("-O3", "-ffast-math", "-funroll-loops", "-fopenmp"),
        ("-fopenmp",),
    ),
    (("-O3", "-fopenmp"), ("-fopenmp",)),
    (("-O3",), ()),
)

# Native modules and build errors are valid only for the exact build identity.
_EXT_CACHE: Dict[tuple[str, str], Any] = {}
_EXT_ERROR: Dict[tuple[str, str], str] = {}
_LAST_FINGERPRINT: Dict[str, str] = {}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(
        payload,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_records(
    source_dir: Path, source_files: Sequence[str]
) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    for source_file in source_files:
        relative = Path(source_file)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"native source path must be relative: {source_file!r}")
        source_path = source_dir / relative
        records.append(
            {
                "name": relative.as_posix(),
                "sha256": _sha256_file(source_path),
            }
        )
    return records


def _run_compiler_probe(command: Sequence[str], *arguments: str) -> str:
    try:
        completed = subprocess.run(
            [*command, *arguments],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return f"unavailable:{type(exc).__name__}"
    return completed.stdout.strip()


def _compiler_command() -> list[str]:
    configured = os.environ.get("CXX")
    if configured:
        return shlex.split(configured) or ["c++"]
    conda_compiler = (
        Path(sys.prefix) / "bin" / f"{platform.machine()}-conda-linux-gnu-g++"
    )
    if conda_compiler.is_file() and os.access(conda_compiler, os.X_OK):
        return [str(conda_compiler)]
    return ["c++"]


def _compiler_identity() -> dict[str, Any]:
    command = _compiler_command()
    displayed_command = [Path(command[0]).name, *command[1:]]
    return {
        "command": displayed_command,
        "version": _run_compiler_probe(command, "--version"),
        "version_number": _run_compiler_probe(
            command, "-dumpfullversion", "-dumpversion"
        ),
        "target": _run_compiler_probe(command, "-dumpmachine"),
    }


def _cpu_identity() -> dict[str, Any]:
    features: set[str] = set()
    cpuinfo = Path("/proc/cpuinfo")
    try:
        contents = cpuinfo.read_text(encoding="utf-8", errors="replace")
        for line in contents.splitlines():
            label, separator, value = line.partition(":")
            if separator and label.strip().lower() in {"flags", "features"}:
                features.update(value.split())
    except OSError:
        pass

    torch_capability: Optional[str]
    try:
        import torch

        torch_capability = str(torch.backends.cpu.get_cpu_capability())
    except Exception:
        torch_capability = None
    return {
        "machine": platform.machine(),
        "processor": platform.processor(),
        "torch_capability": torch_capability,
        "features": sorted(features),
    }


def _runtime_identity() -> dict[str, Any]:
    import torch

    torch_c = getattr(torch, "_C", None)
    return {
        "python": {
            "implementation": platform.python_implementation(),
            "cache_tag": getattr(sys.implementation, "cache_tag", None),
            "soabi": sysconfig.get_config_var("SOABI"),
            "extension_suffix": sysconfig.get_config_var("EXT_SUFFIX"),
        },
        "torch": {
            "version": str(torch.__version__),
            "cuda_version": getattr(torch.version, "cuda", None),
            "hip_version": getattr(torch.version, "hip", None),
            "debug_build": bool(getattr(torch.version, "debug", False)),
            "cxx11_abi": getattr(torch_c, "_GLIBCXX_USE_CXX11_ABI", None),
            "openmp": bool(torch.backends.openmp.is_available()),
            "mkl": bool(torch.backends.mkl.is_available()),
            "mkldnn": bool(torch.backends.mkldnn.is_available()),
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "tag": sysconfig.get_platform(),
            "libc": list(platform.libc_ver()),
            "byteorder": sys.byteorder,
        },
    }


def _native_build_identity(
    source_dir: Path,
    source_files: Sequence[str],
    *,
    build_recipes: Sequence[tuple[Sequence[str], Sequence[str]]] = _BUILD_RECIPES,
    runtime_identity: Optional[Mapping[str, Any]] = None,
    compiler_identity: Optional[Mapping[str, Any]] = None,
    cpu_identity: Optional[Mapping[str, Any]] = None,
) -> tuple[dict[str, Any], str]:
    """Return the path-independent native-build payload and its SHA-256."""

    payload = {
        "schema_version": _FINGERPRINT_SCHEMA,
        "sources": _source_records(Path(source_dir), source_files),
        "build_recipes": [
            {"cflags": list(cflags), "ldflags": list(ldflags)}
            for cflags, ldflags in build_recipes
        ],
        "runtime": dict(
            _runtime_identity() if runtime_identity is None else runtime_identity
        ),
        "compiler": dict(
            _compiler_identity() if compiler_identity is None else compiler_identity
        ),
        "cpu": dict(_cpu_identity() if cpu_identity is None else cpu_identity),
    }
    fingerprint = hashlib.sha256(
        _canonical_json(payload).encode("utf-8")
    ).hexdigest()
    return payload, fingerprint


def _module_name(ext_name: str, fingerprint: str, recipe_index: int = 0) -> str:
    return (
        f"{ext_name}_r{int(recipe_index)}_"
        f"{fingerprint[:_FINGERPRINT_SUFFIX_LENGTH]}"
    )


def _sidecar_path(binary_path: Path) -> Path:
    return binary_path.with_name(f"{binary_path.name}.identity.json")


def _sidecar_payload(
    *,
    identity: Mapping[str, Any],
    fingerprint: str,
    module_name: str,
    binary_path: Path,
) -> dict[str, Any]:
    return {
        "schema_version": _FINGERPRINT_SCHEMA,
        "fingerprint": fingerprint,
        "module_name": module_name,
        "identity": dict(identity),
        "binary": {
            "name": binary_path.name,
            "sha256": _sha256_file(binary_path),
        },
    }


def _write_identity_sidecar(
    *,
    identity: Mapping[str, Any],
    fingerprint: str,
    module_name: str,
    binary_path: Path,
) -> None:
    sidecar = _sidecar_path(binary_path)
    encoded = (
        _canonical_json(
            _sidecar_payload(
                identity=identity,
                fingerprint=fingerprint,
                module_name=module_name,
                binary_path=binary_path,
            )
        )
        + "\n"
    ).encode("utf-8")
    temporary: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{sidecar.name}.",
            suffix=".tmp",
            dir=sidecar.parent,
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            handle.write(encoded)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, sidecar)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _has_valid_identity_sidecar(
    binary_path: Path,
    *,
    identity: Mapping[str, Any],
    fingerprint: str,
    module_name: str,
) -> bool:
    sidecar = _sidecar_path(binary_path)
    try:
        payload = json.loads(sidecar.read_text(encoding="utf-8"))
        expected = _sidecar_payload(
            identity=identity,
            fingerprint=fingerprint,
            module_name=module_name,
            binary_path=binary_path,
        )
    except (OSError, ValueError, TypeError):
        return False
    return payload == expected


def _expected_binary_names(module_name: str) -> tuple[str, ...]:
    ext_suffix = sysconfig.get_config_var("EXT_SUFFIX")
    names = [f"{module_name}.so", f"{module_name}.pyd"]
    if ext_suffix:
        names.append(f"{module_name}{ext_suffix}")
    return tuple(dict.fromkeys(names))


def _find_valid_prebuilt(
    build_dirs: Sequence[Path],
    *,
    identity: Mapping[str, Any],
    fingerprint: str,
    module_name: str,
) -> Optional[Path]:
    for build_dir in build_dirs:
        for name in _expected_binary_names(module_name):
            candidate = build_dir / name
            if candidate.is_file() and _has_valid_identity_sidecar(
                candidate,
                identity=identity,
                fingerprint=fingerprint,
                module_name=module_name,
            ):
                return candidate
    return None


def _build_dirs(
    here: Path, build_subdir: str, fingerprint: str
) -> tuple[Path, Path, Path]:
    local_base = Path(
        os.environ.get("TORCH_EXTENSIONS_DIR")
        or (Path(tempfile.gettempdir()) / "mlmm_hessian_ff")
    )
    package_base = here / build_subdir
    cache_base = _cache_build_dir(build_subdir)
    return (
        local_base / build_subdir / fingerprint,
        package_base / fingerprint,
        cache_base / fingerprint,
    )


def _error_for(key: str) -> Optional[str]:
    fingerprint = _LAST_FINGERPRINT.get(key)
    if fingerprint is None:
        return None
    return _EXT_ERROR.get((key, fingerprint))


def _rebuild_hint() -> str:
    return (
        "mlmm-toolkit JIT-compiles its C++ extensions on first use via "
        "torch.utils.cpp_extension. Current PyTorch releases select C++20; use "
        "a compatible compiler (validated with GCC 13.3) plus ninja "
        "(a pip dependency).\n"
        "If the build fails:\n"
        "  - check C++20 support:  g++ -std=c++20 -x c++ -fsyntax-only /dev/null\n"
        "  - install a modern compiler:  conda install -c conda-forge gxx_linux-64"
        "   (or  apt install g++ / yum install gcc-c++; on HPC, load the site's "
        "C++20-capable compiler module)\n"
        "  - a network-mounted build dir (NFS/Lustre) can hang torch's build "
        "lock; the build defaults to a local temp dir (override via "
        "TORCH_EXTENSIONS_DIR).\n"
        "To rebuild manually:\n"
        "  cd $(python -c \"import hessian_ff; print(hessian_ff.__path__[0])\")/native && make clean && make\n"
        "See docs/troubleshooting.md ('hessian_ff build / import') for details."
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


def _load_or_build_extension(
    *,
    key: str,
    ext_name: str,
    source_files: list[str],
    build_subdir: str,
    verbose: bool,
    force_rebuild: bool,
) -> Optional[Any]:
    if not force_rebuild:
        last_fingerprint = _LAST_FINGERPRINT.get(key)
        if last_fingerprint is not None:
            last_cache_key = (key, last_fingerprint)
            if last_cache_key in _EXT_CACHE:
                return _EXT_CACHE[last_cache_key]

    here = Path(__file__).resolve().parent
    identity, fingerprint = _native_build_identity(here, source_files)
    module_names = tuple(
        _module_name(ext_name, fingerprint, recipe_index)
        for recipe_index in range(len(_BUILD_RECIPES))
    )
    cache_key = (key, fingerprint)
    _LAST_FINGERPRINT[key] = fingerprint

    if cache_key in _EXT_CACHE and not force_rebuild:
        return _EXT_CACHE[cache_key]

    build_dirs = _build_dirs(here, build_subdir, fingerprint)
    prebuilt_error: Optional[Exception] = None

    # Only an exact fingerprint directory, module name, identity sidecar, and
    # binary digest can be reused. Legacy constant-name artifacts are ignored.
    if not force_rebuild:
        for module_name in module_names:
            prebuilt = _find_valid_prebuilt(
                build_dirs,
                identity=identity,
                fingerprint=fingerprint,
                module_name=module_name,
            )
            if prebuilt is None:
                continue
            try:
                spec = importlib.util.spec_from_file_location(
                    module_name, str(prebuilt)
                )
                if spec is None or spec.loader is None:
                    raise ImportError(f"spec loader is unavailable for {prebuilt}")
                loaded = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(loaded)
                sys.modules[module_name] = loaded
                _EXT_CACHE[cache_key] = loaded
                _EXT_ERROR.pop(cache_key, None)
                return loaded
            except Exception as exc:
                prebuilt_error = exc

    if cache_key in _EXT_ERROR and not force_rebuild:
        return None

    if force_rebuild:
        _EXT_CACHE.pop(cache_key, None)
        _EXT_ERROR.pop(cache_key, None)
        for module_name in module_names:
            sys.modules.pop(module_name, None)
        # Package-supplied prebuilds are immutable inputs. Only disposable,
        # out-of-tree build directories are cleared for an explicit rebuild.
        for build_dir in (build_dirs[0], build_dirs[2]):
            shutil.rmtree(build_dir, ignore_errors=True)

    try:
        from torch.utils.cpp_extension import load
    except Exception as exc:
        detail = f"torch cpp_extension import failed: {exc}"
        if prebuilt_error is not None:
            detail = f"failed to load verified prebuilt: {prebuilt_error}; {detail}"
        _EXT_ERROR[cache_key] = _with_rebuild_hint(
            detail
        )
        return None

    _ensure_max_jobs()

    srcs = [here / s for s in source_files]

    # Compile only outside the source/package tree. The package directory above
    # remains a read-only fallback for deliberately supplied verified prebuilds.
    compilation_dirs = (build_dirs[0], build_dirs[2])
    compiler_command = _compiler_command()
    restore_cxx = "CXX" not in os.environ and compiler_command != ["c++"]
    if restore_cxx:
        os.environ["CXX"] = shlex.join(compiler_command)
    original_path = os.environ.get("PATH")
    ninja = Path(sys.prefix) / "bin" / "ninja"
    restore_path = shutil.which("ninja") is None and ninja.is_file()
    if restore_path:
        prefix = str(ninja.parent)
        os.environ["PATH"] = (
            prefix if not original_path else prefix + os.pathsep + original_path
        )
    last_err: Optional[Exception] = None
    try:
        for build_dir in compilation_dirs:
            try:
                os.makedirs(build_dir, exist_ok=True)
            except OSError:
                continue
            for recipe_index, (cflags, ldflags) in enumerate(_BUILD_RECIPES):
                module_name = module_names[recipe_index]
                try:
                    loaded = load(
                        name=module_name,
                        sources=[str(source) for source in srcs],
                        extra_cflags=list(cflags),
                        extra_ldflags=list(ldflags),
                        build_directory=str(build_dir),
                        verbose=bool(verbose),
                    )
                    binary_path = Path(getattr(loaded, "__file__", "")).resolve()
                    if not binary_path.is_file():
                        raise ImportError(
                            "built extension did not expose a binary file: "
                            f"{binary_path}"
                        )
                    if binary_path.parent != build_dir.resolve():
                        raise ImportError(
                            "built extension resolved outside its fingerprinted "
                            f"directory: {binary_path}"
                        )
                    if binary_path.name not in _expected_binary_names(module_name):
                        raise ImportError(
                            "built extension has unexpected module name: "
                            f"{binary_path.name}"
                        )
                    if _source_records(here, source_files) != identity["sources"]:
                        raise RuntimeError("native source changed during compilation")
                    _write_identity_sidecar(
                        identity=identity,
                        fingerprint=fingerprint,
                        module_name=module_name,
                        binary_path=binary_path,
                    )
                    if not _has_valid_identity_sidecar(
                        binary_path,
                        identity=identity,
                        fingerprint=fingerprint,
                        module_name=module_name,
                    ):
                        raise RuntimeError(
                            "native extension identity sidecar verification failed"
                        )
                    _EXT_CACHE[cache_key] = loaded
                    _EXT_ERROR.pop(cache_key, None)
                    return loaded
                except Exception as exc:
                    last_err = exc
                    continue
    finally:
        if restore_cxx:
            os.environ.pop("CXX", None)
        if restore_path:
            if original_path is None:
                os.environ.pop("PATH", None)
            else:
                os.environ["PATH"] = original_path

    detail = f"hessian_ff build attempts failed: {last_err}"
    if prebuilt_error is not None:
        detail = f"failed to load verified prebuilt: {prebuilt_error}; {detail}"
    _EXT_ERROR[cache_key] = _with_rebuild_hint(detail)
    return None


def get_nonbonded_extension(
    *,
    verbose: bool = False,
    force_rebuild: bool = False,
) -> Optional[Any]:
    """Load or build the fingerprinted C++ extension for nonbonded kernels.

    References:
    - torch.utils.cpp_extension.load() runtime build workflow.
    """
    return _load_or_build_extension(
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
        note = _error_for("nonbonded") or _with_rebuild_hint("extension unavailable")
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
    """Load or build the fingerprinted analytical-Hessian extension."""
    return _load_or_build_extension(
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
        note = _error_for("analytical_hessian") or _with_rebuild_hint(
            "extension unavailable"
        )
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
    """Load or build the fingerprinted bonded energy-force extension."""
    return _load_or_build_extension(
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
        note = _error_for("bonded") or _with_rebuild_hint("extension unavailable")
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
    all fingerprinted extension builds in one step before production runs.
    """
    ext_nb = get_nonbonded_extension(verbose=verbose, force_rebuild=force_rebuild)
    ext_ah = get_analytical_hessian_extension(verbose=verbose, force_rebuild=force_rebuild)
    ext_bd = get_bonded_extension(verbose=verbose, force_rebuild=force_rebuild)
    return {
        "nonbonded": "ok"
        if ext_nb is not None
        else f"error: {_error_for('nonbonded')}",
        "analytical_hessian": "ok"
        if ext_ah is not None
        else f"error: {_error_for('analytical_hessian')}",
        "bonded": "ok"
        if ext_bd is not None
        else f"error: {_error_for('bonded')}",
    }
