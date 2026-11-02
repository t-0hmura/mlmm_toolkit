"""gpu4pyscf finds libcusolver in the active environment when the system cannot supply it.

CPU-only: ``gpu4pyscf`` is a stub package whose import records what
``ctypes.util.find_library`` answers, and the system lookup is a fake.
"""

from __future__ import annotations

import ctypes.util
import sys
import types
from pathlib import Path

import pytest

from mlmm.backends import pyscf_dft

SYSTEM_OTHER = "libother-from-system.so"


@pytest.fixture
def environment(tmp_path, monkeypatch):
    """A stub ``gpu4pyscf``, an environment ``lib`` with libcusolver, and no cuda-pathfinder."""
    site = tmp_path / "site"
    (site / "gpu4pyscf").mkdir(parents=True)
    stub = site / "gpu4pyscf" / "__init__.py"
    stub.write_text(
        "import ctypes.util\n"
        "cusolver_seen = ctypes.util.find_library('cusolver')\n"
        "other_seen = ctypes.util.find_library('other')\n"
    )
    monkeypatch.syspath_prepend(str(site))
    for name in [n for n in sys.modules if n == "gpu4pyscf" or n.startswith("gpu4pyscf.")]:
        monkeypatch.delitem(sys.modules, name)

    prefix = tmp_path / "prefix"
    (prefix / "lib").mkdir(parents=True)
    for name in ("libcusolver.so.11", "libcusolver.so.12", "libcusolver.so.12.3.5.3"):
        (prefix / "lib" / name).touch()
    monkeypatch.setattr(sys, "prefix", str(prefix))
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
    monkeypatch.setitem(sys.modules, "cuda.pathfinder", None)  # its import raises ImportError

    monkeypatch.setattr(pyscf_dft, "_cusolver_checked", False)
    yield types.SimpleNamespace(lib=prefix / "lib", stub=stub)
    sys.modules.pop("gpu4pyscf", None)  # the stub, before monkeypatch restores the real package


def _system_lookup(monkeypatch, cusolver):
    calls = []

    def lookup(name):
        calls.append(name)
        return cusolver if name == "cusolver" else SYSTEM_OTHER

    monkeypatch.setattr(ctypes.util, "find_library", lookup)
    return lookup, calls


def _pathfinder(monkeypatch, result):
    module = types.ModuleType("cuda.pathfinder")

    def load_nvidia_dynamic_lib(name):
        assert name == "cusolver"
        if isinstance(result, Exception):
            raise result
        return types.SimpleNamespace(abs_path=result)

    module.load_nvidia_dynamic_lib = load_nvidia_dynamic_lib
    monkeypatch.setitem(sys.modules, "cuda", types.ModuleType("cuda"))
    monkeypatch.setitem(sys.modules, "cuda.pathfinder", module)


@pytest.mark.parametrize("source", ["pathfinder", "pathfinder_fails", "pathfinder_absent"])
def test_environment_libcusolver_answers_find_library_during_the_import(
    environment, monkeypatch, source
) -> None:
    lookup, _ = _system_lookup(monkeypatch, None)
    if source == "pathfinder":
        _pathfinder(monkeypatch, "/env/lib/libcusolver.so.12")
        expected = "/env/lib/libcusolver.so.12"
    else:
        if source == "pathfinder_fails":
            _pathfinder(monkeypatch, RuntimeError("no libcusolver"))
        expected = str(environment.lib / "libcusolver.so.12")

    pyscf_dft.ensure_gpu4pyscf_cusolver()

    stub = sys.modules["gpu4pyscf"]
    assert stub.cusolver_seen == expected
    assert stub.other_seen == SYSTEM_OTHER
    assert ctypes.util.find_library is lookup


def test_working_system_lookup_changes_nothing(environment, monkeypatch) -> None:
    lookup, calls = _system_lookup(monkeypatch, "libcusolver.so.12")

    pyscf_dft.ensure_gpu4pyscf_cusolver()

    assert "gpu4pyscf" not in sys.modules
    assert ctypes.util.find_library is lookup
    assert calls == ["cusolver"]


def test_no_copy_in_the_environment_changes_nothing(environment, monkeypatch) -> None:
    for path in environment.lib.iterdir():
        path.unlink()
    lookup, _ = _system_lookup(monkeypatch, None)

    pyscf_dft.ensure_gpu4pyscf_cusolver()

    assert "gpu4pyscf" not in sys.modules
    assert ctypes.util.find_library is lookup


def test_gpu4pyscf_import_error_is_left_to_the_caller(environment, monkeypatch) -> None:
    environment.stub.write_text("raise OSError('Unable to load libcusolver')\n")
    lookup, _ = _system_lookup(monkeypatch, None)

    pyscf_dft.ensure_gpu4pyscf_cusolver()

    assert ctypes.util.find_library is lookup
    with pytest.raises(OSError, match="Unable to load libcusolver"):
        __import__("gpu4pyscf")


def test_the_lookup_runs_once_per_process(environment, monkeypatch) -> None:
    _, calls = _system_lookup(monkeypatch, "libcusolver.so.12")

    pyscf_dft.ensure_gpu4pyscf_cusolver()
    pyscf_dft.ensure_gpu4pyscf_cusolver()

    assert calls == ["cusolver"]


@pytest.mark.parametrize(("engine", "prepared"), [("gpu", True), ("cpu", False)])
def test_method_builder_prepares_cusolver_before_the_first_gpu4pyscf_import(
    monkeypatch, engine, prepared
) -> None:
    pyscf = pytest.importorskip("pyscf")
    from mlmm.core.dft_settings import resolve_dft_settings

    events = []
    monkeypatch.setattr(
        pyscf_dft, "ensure_gpu4pyscf_cusolver", lambda: events.append("prepare")
    )

    class Lowmem:
        def __init__(self, mol, xc):
            events.append("gpu4pyscf")

    package = types.ModuleType("gpu4pyscf")
    dft = types.ModuleType("gpu4pyscf.dft")
    lowmem = types.ModuleType("gpu4pyscf.dft.rks_lowmem")
    lowmem.RKS = Lowmem
    dft.rks_lowmem = lowmem
    package.dft = dft
    for module in (package, dft, lowmem):
        monkeypatch.setitem(sys.modules, module.__name__, module)

    settings = resolve_dft_settings(
        {"backend": "dft", "dft": {"func_basis": "hf/sto-3g", "engine": engine}}
    )
    mol = object() if engine == "gpu" else pyscf.M(atom="He 0 0 0", basis="sto-3g", verbose=0)
    pyscf_dft.PySCFDFTSession(settings)._make_method(mol)

    assert events == (["prepare", "gpu4pyscf"] if prepared else [])


def test_dft_command_prepares_cusolver_before_its_first_gpu4pyscf_import(
    environment, monkeypatch, tmp_path
) -> None:
    pyscf_misc = pytest.importorskip("pyscf.lib.misc")
    from ase import Atoms
    from click.testing import CliRunner

    from mlmm.cli import cli as root_cli
    from mlmm.workflows import dft

    # A real gpu4pyscf imported earlier replaces this banner with one that reads gpu4pyscf.__config__.
    monkeypatch.setattr(pyscf_misc, "format_sys_info", lambda: [])

    toy = Path(__file__).resolve().parents[1] / "examples" / "toy_system"
    if not (toy / "p_toy.parm7").exists():
        pytest.skip("toy_system example inputs not present")
    _system_lookup(monkeypatch, None)
    model = Atoms("F", positions=[(0.0, 0.0, 0.0)])
    workspace = types.SimpleNamespace(
        atoms_model=model, atoms_model_lh=model.copy(), link_pairs=[], cleanup=lambda: None,
    )
    monkeypatch.setattr(dft, "_prepare_ml_region_workspace", lambda **kwargs: workspace)

    result = CliRunner().invoke(root_cli, [
        "dft", "-i", str(toy / "p_toy.pdb"), "--parm", str(toy / "p_toy.parm7"),
        "-q", "-1", "-m", "1", "--func-basis", "hf/sto-3g", "--no-convert-files",
        "-o", str(tmp_path / "out"),
    ])

    # The stub has no gpu4pyscf.dft, so the GPU setup stops right after the import.
    assert "[gpu] GPU backend failed" in result.output, result.output
    assert sys.modules["gpu4pyscf"].cusolver_seen == str(environment.lib / "libcusolver.so.12")
