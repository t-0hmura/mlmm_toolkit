from __future__ import annotations

from copy import deepcopy
from importlib.machinery import ModuleSpec
import json
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest

from hessian_ff.native import loader


RUNTIME_IDENTITY = {
    "python": {
        "implementation": "CPython",
        "cache_tag": "cpython-311",
        "soabi": "cpython-311-test",
        "extension_suffix": ".test.so",
    },
    "torch": {
        "version": "2.8.0",
        "cuda_version": None,
        "hip_version": None,
        "debug_build": False,
        "cxx11_abi": True,
        "openmp": True,
        "mkl": True,
        "mkldnn": True,
    },
    "platform": {
        "system": "Linux",
        "release": "test",
        "tag": "linux-test",
        "libc": ["glibc", "2.test"],
        "byteorder": "little",
    },
}
COMPILER_IDENTITY = {
    "command": ["c++"],
    "version": "test compiler 1",
    "version_number": "1",
    "target": "test-linux-gnu",
}
CPU_IDENTITY = {
    "machine": "test64",
    "processor": "test cpu",
    "torch_capability": "TEST",
    "features": ["feature_a", "feature_b"],
}
RECIPES = ((('-O3', '-march=native'), ()), (('-O3',), ()))


@pytest.fixture(autouse=True)
def _clear_loader_state():
    loader._EXT_CACHE.clear()
    loader._EXT_ERROR.clear()
    loader._LAST_FINGERPRINT.clear()
    yield
    loader._EXT_CACHE.clear()
    loader._EXT_ERROR.clear()
    loader._LAST_FINGERPRINT.clear()


def _identity(source_dir: Path, source_files=("kernel.cpp",), **overrides):
    options = {
        "build_recipes": RECIPES,
        "runtime_identity": RUNTIME_IDENTITY,
        "compiler_identity": COMPILER_IDENTITY,
        "cpu_identity": CPU_IDENTITY,
    }
    options.update(overrides)
    return loader._native_build_identity(source_dir, source_files, **options)


def test_fingerprint_is_independent_of_source_absolute_directory(tmp_path):
    source_dirs = [tmp_path / "one", tmp_path / "somewhere" / "two"]
    for source_dir in source_dirs:
        source_dir.mkdir(parents=True)
        (source_dir / "kernel.cpp").write_bytes(b"same source bytes\n")

    first_payload, first = _identity(source_dirs[0])
    second_payload, second = _identity(source_dirs[1])

    assert first == second
    assert first_payload == second_payload
    canonical = loader._canonical_json(first_payload)
    assert str(source_dirs[0]) not in canonical
    assert str(source_dirs[1]) not in canonical


def test_source_byte_change_changes_fingerprint(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    source = source_dir / "kernel.cpp"
    source.write_bytes(b"first")
    _, first = _identity(source_dir)
    source.write_bytes(b"second")
    _, second = _identity(source_dir)
    assert first != second


def test_runtime_compiler_cpu_and_recipe_changes_each_change_fingerprint(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    (source_dir / "kernel.cpp").write_bytes(b"kernel")
    _, baseline = _identity(source_dir)

    python_runtime = deepcopy(RUNTIME_IDENTITY)
    python_runtime["python"]["soabi"] = "cpython-312-test"
    torch_runtime = deepcopy(RUNTIME_IDENTITY)
    torch_runtime["torch"]["cxx11_abi"] = False
    compiler = deepcopy(COMPILER_IDENTITY)
    compiler["target"] = "different-target"
    cpu = deepcopy(CPU_IDENTITY)
    cpu["torch_capability"] = "DIFFERENT"
    recipes = ((('-O2',), ()),)

    variants = [
        {"runtime_identity": python_runtime},
        {"runtime_identity": torch_runtime},
        {"compiler_identity": compiler},
        {"cpu_identity": cpu},
        {"build_recipes": recipes},
    ]
    for variant in variants:
        _, fingerprint = _identity(source_dir, **variant)
        assert fingerprint != baseline


def test_legacy_and_wrong_fingerprint_binaries_are_not_candidates(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    (source_dir / "kernel.cpp").write_bytes(b"kernel")
    identity, fingerprint = _identity(source_dir)
    module_name = loader._module_name("example_ext", fingerprint)
    fingerprint_dir = tmp_path / fingerprint
    fingerprint_dir.mkdir()

    (fingerprint_dir / "example_ext.so").write_bytes(b"legacy")
    wrong_module = loader._module_name("example_ext", "f" * 64)
    (fingerprint_dir / f"{wrong_module}.so").write_bytes(b"wrong")

    assert loader._find_valid_prebuilt(
        [fingerprint_dir],
        identity=identity,
        fingerprint=fingerprint,
        module_name=module_name,
    ) is None


def test_prebuilt_requires_exact_sidecar_and_binary_digest(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    (source_dir / "kernel.cpp").write_bytes(b"kernel")
    identity, fingerprint = _identity(source_dir)
    module_name = loader._module_name("example_ext", fingerprint)
    fingerprint_dir = tmp_path / fingerprint
    fingerprint_dir.mkdir()
    binary = fingerprint_dir / f"{module_name}.so"
    binary.write_bytes(b"verified binary")

    def find():
        return loader._find_valid_prebuilt(
            [fingerprint_dir],
            identity=identity,
            fingerprint=fingerprint,
            module_name=module_name,
        )

    assert find() is None
    loader._write_identity_sidecar(
        identity=identity,
        fingerprint=fingerprint,
        module_name=module_name,
        binary_path=binary,
    )
    assert find() == binary

    sidecar = loader._sidecar_path(binary)
    payload = json.loads(sidecar.read_text(encoding="utf-8"))
    payload["fingerprint"] = "0" * 64
    sidecar.write_text(json.dumps(payload), encoding="utf-8")
    assert find() is None

    loader._write_identity_sidecar(
        identity=identity,
        fingerprint=fingerprint,
        module_name=module_name,
        binary_path=binary,
    )
    binary.write_bytes(b"corrupt binary")
    assert find() is None


def _mock_identity(monkeypatch, state, source_files):
    here = Path(loader.__file__).resolve().parent
    source_records = loader._source_records(here, source_files)

    def fake_identity(source_dir, selected_files):
        assert Path(source_dir) == here
        assert list(selected_files) == list(source_files)
        fingerprint = state["fingerprint"]
        return {
            "schema_version": 1,
            "sources": source_records,
            "build_recipes": [],
            "runtime": {"marker": fingerprint},
            "compiler": {},
            "cpu": {},
        }, fingerprint

    monkeypatch.setattr(loader, "_native_build_identity", fake_identity)


def _mock_cpp_load(monkeypatch, calls):
    import torch.utils.cpp_extension

    def fake_load(**kwargs):
        calls.append(kwargs)
        binary = Path(kwargs["build_directory"]) / f"{kwargs['name']}.so"
        binary.parent.mkdir(parents=True, exist_ok=True)
        binary.write_bytes(f"binary-{len(calls)}".encode("ascii"))
        return SimpleNamespace(__file__=str(binary))

    monkeypatch.setattr(torch.utils.cpp_extension, "load", fake_load)


def _build_test_extension(*, key, force_rebuild=False):
    return loader._load_or_build_extension(
        key=key,
        ext_name="test_native_ext",
        source_files=["nonbonded_ext.cpp"],
        build_subdir=".build_test_native",
        verbose=False,
        force_rebuild=force_rebuild,
    )


def test_verified_prebuild_loads_without_jit_and_mismatch_builds(
    tmp_path, monkeypatch
):
    torch_root = tmp_path / "torch"
    monkeypatch.setenv("TORCH_EXTENSIONS_DIR", str(torch_root))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))
    state = {"fingerprint": "0" * 64}
    source_files = ["nonbonded_ext.cpp"]
    _mock_identity(monkeypatch, state, source_files)

    here = Path(loader.__file__).resolve().parent
    identity, fingerprint = loader._native_build_identity(here, source_files)
    module_name = loader._module_name("test_native_ext", fingerprint)
    build_dir = loader._build_dirs(
        here, ".build_test_native", fingerprint
    )[0]
    build_dir.mkdir(parents=True)
    binary = build_dir / f"{module_name}.so"
    binary.write_bytes(b"verified prebuild")
    loader._write_identity_sidecar(
        identity=identity,
        fingerprint=fingerprint,
        module_name=module_name,
        binary_path=binary,
    )

    imported = []

    class FakeExtensionLoader:
        def __init__(self, location):
            self.location = location

        def create_module(self, spec):
            return None

        def exec_module(self, module):
            module.__file__ = self.location
            imported.append((module.__name__, self.location))

    def fake_spec_from_file_location(name, location):
        return ModuleSpec(
            name,
            FakeExtensionLoader(location),
            origin=location,
        )

    jit_calls = []
    _mock_cpp_load(monkeypatch, jit_calls)
    monkeypatch.setattr(
        loader.importlib.util,
        "spec_from_file_location",
        fake_spec_from_file_location,
    )

    prebuilt = _build_test_extension(key="test")
    assert isinstance(prebuilt, ModuleType)
    assert imported == [(module_name, str(binary))]
    assert jit_calls == []

    loader._EXT_CACHE.clear()
    loader._sidecar_path(binary).unlink()
    rebuilt = _build_test_extension(key="test")
    assert rebuilt is not None
    assert len(jit_calls) == 1
    assert jit_calls[0]["name"] == module_name


def test_error_cache_is_scoped_to_fingerprint(tmp_path, monkeypatch):
    monkeypatch.setenv("TORCH_EXTENSIONS_DIR", str(tmp_path / "torch"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))
    state = {"fingerprint": "1" * 64}
    _mock_identity(monkeypatch, state, ["nonbonded_ext.cpp"])
    calls = []
    _mock_cpp_load(monkeypatch, calls)
    loader._EXT_ERROR[("test", state["fingerprint"])] = "failure for F1"

    assert _build_test_extension(key="test") is None
    assert calls == []

    state["fingerprint"] = "2" * 64
    assert _build_test_extension(key="test") is not None
    assert len(calls) == 1


def test_error_cache_does_not_hide_new_verified_prebuild(tmp_path, monkeypatch):
    monkeypatch.setenv("TORCH_EXTENSIONS_DIR", str(tmp_path / "torch"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))
    state = {"fingerprint": "4" * 64}
    source_files = ["nonbonded_ext.cpp"]
    _mock_identity(monkeypatch, state, source_files)
    loader._EXT_ERROR[("test", state["fingerprint"])] = "initial failure"

    here = Path(loader.__file__).resolve().parent
    identity, fingerprint = loader._native_build_identity(here, source_files)
    module_name = loader._module_name("test_native_ext", fingerprint)
    build_dir = loader._build_dirs(
        here, ".build_test_native", fingerprint
    )[0]
    build_dir.mkdir(parents=True)
    binary = build_dir / f"{module_name}.so"
    binary.write_bytes(b"verified prebuild")
    loader._write_identity_sidecar(
        identity=identity,
        fingerprint=fingerprint,
        module_name=module_name,
        binary_path=binary,
    )

    class FakeExtensionLoader:
        def create_module(self, spec):
            return None

        def exec_module(self, module):
            module.__file__ = str(binary)

    monkeypatch.setattr(
        loader.importlib.util,
        "spec_from_file_location",
        lambda name, location: ModuleSpec(
            name,
            FakeExtensionLoader(),
            origin=location,
        ),
    )

    loaded = _build_test_extension(key="test")

    assert loaded is not None
    assert ("test", fingerprint) not in loader._EXT_ERROR


def test_force_rebuild_bypasses_memory_and_prebuilt_cache(tmp_path, monkeypatch):
    torch_root = tmp_path / "torch"
    monkeypatch.setenv("TORCH_EXTENSIONS_DIR", str(torch_root))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))
    state = {"fingerprint": "3" * 64}
    _mock_identity(monkeypatch, state, ["nonbonded_ext.cpp"])
    calls = []
    _mock_cpp_load(monkeypatch, calls)

    first = _build_test_extension(key="test")
    second = _build_test_extension(key="test")
    rebuilt = _build_test_extension(key="test", force_rebuild=True)

    assert first is second
    assert rebuilt is not None
    assert len(calls) == 2
    expected_name = loader._module_name("test_native_ext", state["fingerprint"])
    for call in calls:
        assert call["name"] == expected_name
        build_dir = Path(call["build_directory"])
        assert build_dir.is_relative_to(torch_root)
        assert build_dir.name == state["fingerprint"]
        assert not build_dir.is_relative_to(Path(loader.__file__).resolve().parent)


def test_successful_memory_cache_skips_recomputing_build_identity(
    tmp_path, monkeypatch,
):
    monkeypatch.setenv("TORCH_EXTENSIONS_DIR", str(tmp_path / "torch"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))
    state = {"fingerprint": "5" * 64}
    source_files = ["nonbonded_ext.cpp"]
    _mock_identity(monkeypatch, state, source_files)
    original_identity = loader._native_build_identity
    identity_calls = 0

    def counted_identity(*args, **kwargs):
        nonlocal identity_calls
        identity_calls += 1
        return original_identity(*args, **kwargs)

    monkeypatch.setattr(loader, "_native_build_identity", counted_identity)
    calls = []
    _mock_cpp_load(monkeypatch, calls)

    first = _build_test_extension(key="test")
    second = _build_test_extension(key="test")

    assert first is second
    assert identity_calls == 1
    assert len(calls) == 1
