"""M56: the enforceable AST import-graph gate.

Proves the real repo is clean (no mlmm cycles, no bundled-fork -> product edge,
no core/domain -> workflows edge) AND that the gate actually fires on a known-bad
graph while NOT firing on a legitimate ``from . import child`` /
``if TYPE_CHECKING:`` graph.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / ".github" / "scripts" / "check_import_graph.py"


def _load_gate():
    spec = importlib.util.spec_from_file_location("mlmm_check_import_graph", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


GATE = _load_gate()


def _write(base: Path, rel: str, text: str = "") -> None:
    p = base / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")


# --------------------------------------------------------------------------- #
# Real repository must pass every invariant.
# --------------------------------------------------------------------------- #
def test_real_repo_is_clean():
    assert GATE.check_repo(REPO_ROOT) == []


def test_no_bundled_fork_imports_the_product():
    """Explicit falsifier: no file under pysisyphus/** (or other forks) imports mlmm."""
    mod2path, edges = GATE.build_graph(
        [(GATE.PRODUCT, str(REPO_ROOT / GATE.PRODUCT))]
        + [(f, str(REPO_ROOT / f)) for f in GATE.BUNDLED_FORKS]
    )
    assert GATE.fork_to_product_edges(mod2path, edges) == []


def test_no_product_scc_in_real_repo():
    mod2path, edges = GATE.build_graph([(GATE.PRODUCT, str(REPO_ROOT / GATE.PRODUCT))])
    assert GATE.product_multi_sccs(mod2path, edges) == []


# --------------------------------------------------------------------------- #
# Negative control: a known-bad cycle + forbidden edge MUST be detected.
# --------------------------------------------------------------------------- #
def test_gate_fires_on_known_bad_cycle(tmp_path: Path):
    root = tmp_path
    _write(root, "mlmm/__init__.py")
    _write(root, "mlmm/core/__init__.py")
    _write(root, "mlmm/core/a.py", "from mlmm.workflows import b\n")
    _write(root, "mlmm/workflows/__init__.py")
    _write(root, "mlmm/workflows/b.py", "from mlmm.core import a\n")

    mod2path, edges = GATE.build_graph([("mlmm", str(root / "mlmm"))])

    sccs = GATE.product_multi_sccs(mod2path, edges)
    assert sccs == [["mlmm.core.a", "mlmm.workflows.b"]], sccs
    # the forbidden layer edge is also reported
    assert ("mlmm.core.a", "mlmm.workflows.b") in GATE.forbidden_layer_edges(edges)


def test_gate_fires_on_fork_importing_product(tmp_path: Path):
    root = tmp_path
    _write(root, "mlmm/__init__.py")
    _write(root, "mlmm/core/__init__.py")
    _write(root, "pysisyphus/__init__.py")
    _write(root, "pysisyphus/bad.py", "from mlmm.core import thing\n")
    _write(root, "mlmm/core/thing.py")

    mod2path, edges = GATE.build_graph(
        [("mlmm", str(root / "mlmm")), ("pysisyphus", str(root / "pysisyphus"))]
    )
    assert ("pysisyphus.bad", "mlmm.core.thing") in GATE.fork_to_product_edges(mod2path, edges)


# --------------------------------------------------------------------------- #
# Positive control: package->submodule + TYPE_CHECKING must NOT be a cycle.
# --------------------------------------------------------------------------- #
def test_gate_does_not_false_positive_on_relative_and_type_checking(tmp_path: Path):
    root = tmp_path
    # ``from . import child`` (package importing its own submodule) is a normal
    # edge, not a cycle.
    _write(root, "mlmm/__init__.py", "from . import core\n")
    _write(root, "mlmm/core/__init__.py", "from . import a\n")
    _write(
        root,
        "mlmm/core/a.py",
        "from __future__ import annotations\n"
        "from typing import TYPE_CHECKING\n"
        "if TYPE_CHECKING:\n"
        "    from mlmm.workflows.b import Thing  # runtime-excluded\n",
    )
    _write(root, "mlmm/workflows/__init__.py")
    _write(root, "mlmm/workflows/b.py", "from mlmm.core import a\n")

    mod2path, edges = GATE.build_graph([("mlmm", str(root / "mlmm"))])

    # package -> submodule edges exist but form no cycle
    assert "mlmm.core" in edges["mlmm"]
    assert "mlmm.core.a" in edges["mlmm.core"]
    # TYPE_CHECKING import was excluded: core.a has no runtime edge to workflows.b
    assert "mlmm.workflows.b" not in edges.get("mlmm.core.a", set())
    assert GATE.product_multi_sccs(mod2path, edges) == []
    assert GATE.forbidden_layer_edges(edges) == []


# --------------------------------------------------------------------------- #
# The docs-quality driver must actually invoke the gate.
# --------------------------------------------------------------------------- #
def test_docs_quality_driver_runs_the_graph_gate():
    driver = (REPO_ROOT / ".github" / "scripts" / "run_docs_quality.py").read_text(
        encoding="utf-8"
    )
    assert "check_import_graph.py" in driver, (
        "run_docs_quality.py must run check_import_graph.py; removing that step "
        "disables the import-graph gate."
    )
