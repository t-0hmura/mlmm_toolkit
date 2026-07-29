#!/usr/bin/env python3
"""Enforce MLMM's product import-graph invariants.

This is the checker the architecture docs point to for the dependency direction
— unlike ``check_engineering_markers.py`` (chemistry / DOMAIN_PURE / MLIP scope),
this one actually parses the static import edges.  It asserts:

1. The ``mlmm`` package has **no strongly connected component** among its own
   modules (no import cycle).
2. No bundled fork (``pysisyphus`` / ``hessian_ff`` / ``thermoanalysis``) imports
   ``mlmm`` — the forks stay leaves that the product imports, never the reverse
   (this also rejects reverse edges into product workflows).
3. No forbidden layer edge: ``mlmm.core.* -> mlmm.workflows.*`` and
   ``mlmm.domain.* -> mlmm.workflows.*`` (``core`` / ``domain`` never import the
   application layer).

The graph builder resolves absolute imports, relative imports, package
``__init__.py``, aliases, and ``from . import child`` (a package importing its own
submodule is a normal edge, not a cycle).  ``if TYPE_CHECKING:`` blocks are
excluded — they never execute, so they cannot form a runtime cycle.
"""

from __future__ import annotations

import argparse
import ast
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

REPO_ROOT = Path(__file__).resolve().parents[2]

# Product package + bundled forks (forks must never import the product).
PRODUCT = "mlmm"
BUNDLED_FORKS = ("pysisyphus", "hessian_ff", "thermoanalysis")


# --------------------------------------------------------------------------- #
# Graph construction
# --------------------------------------------------------------------------- #
def discover_modules(roots: Iterable[Tuple[str, str]]) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Map dotted module name <-> file path for every ``.py`` under each root."""
    mod2path: Dict[str, str] = {}
    path2mod: Dict[str, str] = {}
    for _pkg_name, top in roots:
        top = os.path.abspath(top)
        if not os.path.isdir(top):
            continue
        parent = os.path.dirname(top)
        for dirpath, dirnames, filenames in os.walk(top):
            dirnames[:] = [d for d in dirnames if d != "__pycache__"]
            for fn in filenames:
                if not fn.endswith(".py"):
                    continue
                full = os.path.join(dirpath, fn)
                rel = os.path.relpath(full, parent)
                parts = rel[:-3].split(os.sep)
                if parts[-1] == "__init__":
                    parts = parts[:-1]
                dotted = ".".join(parts)
                mod2path[dotted] = full
                path2mod[full] = dotted
    return mod2path, path2mod


def _is_type_checking_test(test: ast.expr) -> bool:
    """True for ``TYPE_CHECKING`` / ``typing.TYPE_CHECKING`` guard conditions."""
    if isinstance(test, ast.Name):
        return test.id == "TYPE_CHECKING"
    if isinstance(test, ast.Attribute):
        return test.attr == "TYPE_CHECKING"
    return False


def _iter_import_nodes(node: ast.AST):
    """Yield Import / ImportFrom nodes, skipping ``if TYPE_CHECKING:`` bodies.

    The ``else`` branch of such a guard DOES run at runtime, so it is still
    traversed.
    """
    for child in ast.iter_child_nodes(node):
        if isinstance(child, (ast.Import, ast.ImportFrom)):
            yield child
        if isinstance(child, ast.If) and _is_type_checking_test(child.test):
            for sub in child.orelse:
                yield from _iter_import_nodes(sub)
            continue
        yield from _iter_import_nodes(child)


def _resolve(node: ast.AST, cur_mod: str, cur_is_pkg: bool, mod2path: Dict[str, str]) -> Set[str]:
    """Resolve one import node to the set of known target module names."""
    targets: Set[str] = set()
    if isinstance(node, ast.Import):
        for alias in node.names:
            parts = alias.name.split(".")
            for i in range(len(parts), 0, -1):
                cand = ".".join(parts[:i])
                if cand in mod2path:
                    targets.add(cand)
                    break
    elif isinstance(node, ast.ImportFrom):
        level = node.level
        module = node.module
        if level == 0:
            base = module
        else:
            cur_parts = cur_mod.split(".")
            anchor = cur_parts[:] if cur_is_pkg else cur_parts[:-1]
            up = level - 1
            if up > 0:
                anchor = anchor[:-up] if up <= len(anchor) else []
            base_parts = anchor + ([module] if module else [])
            base = ".".join(base_parts) if base_parts else None
        if base is None:
            return targets
        if base in mod2path:
            targets.add(base)
        # `from base import name` — name may itself be a submodule ``base.name``.
        for alias in node.names:
            if alias.name == "*":
                continue
            cand = f"{base}.{alias.name}"
            if cand in mod2path:
                targets.add(cand)
    return targets


def build_graph(roots: Iterable[Tuple[str, str]]) -> Tuple[Dict[str, str], Dict[str, Set[str]]]:
    """Return (mod2path, edges) for every module under the given roots."""
    roots = list(roots)
    mod2path, _path2mod = discover_modules(roots)
    edges: Dict[str, Set[str]] = defaultdict(set)
    for mod, path in mod2path.items():
        cur_is_pkg = os.path.basename(path) == "__init__.py"
        try:
            with open(path, "r", encoding="utf-8") as fh:
                tree = ast.parse(fh.read(), filename=path)
        except (SyntaxError, UnicodeDecodeError) as exc:  # pragma: no cover
            print(f"[import-graph] WARN parse {path}: {exc}", file=sys.stderr)
            continue
        for node in _iter_import_nodes(tree):
            for tgt in _resolve(node, mod, cur_is_pkg, mod2path):
                if tgt != mod:
                    edges[mod].add(tgt)
    return mod2path, edges


# --------------------------------------------------------------------------- #
# Analyses
# --------------------------------------------------------------------------- #
def tarjan_sccs(nodes: Set[str], edges: Dict[str, Set[str]]) -> List[List[str]]:
    """Tarjan's SCC over the induced subgraph on *nodes*."""
    index_counter = [0]
    stack: List[str] = []
    on_stack: Dict[str, bool] = defaultdict(bool)
    index: Dict[str, int] = {}
    low: Dict[str, int] = {}
    result: List[List[str]] = []
    sys.setrecursionlimit(1 << 20)

    def strongconnect(v: str) -> None:
        index[v] = low[v] = index_counter[0]
        index_counter[0] += 1
        stack.append(v)
        on_stack[v] = True
        for w in edges.get(v, ()):
            if w not in nodes:
                continue
            if w not in index:
                strongconnect(w)
                low[v] = min(low[v], low[w])
            elif on_stack[w]:
                low[v] = min(low[v], index[w])
        if low[v] == index[v]:
            comp = []
            while True:
                w = stack.pop()
                on_stack[w] = False
                comp.append(w)
                if w == v:
                    break
            result.append(comp)

    for v in nodes:
        if v not in index:
            strongconnect(v)
    return result


def product_multi_sccs(mod2path: Dict[str, str], edges: Dict[str, Set[str]]) -> List[List[str]]:
    """Multi-module SCCs among ``mlmm`` product modules (each sorted)."""
    nodes = {m for m in mod2path if m == PRODUCT or m.startswith(PRODUCT + ".")}
    return sorted(
        (sorted(c) for c in tarjan_sccs(nodes, edges) if len(c) > 1),
        key=lambda c: (-len(c), c),
    )


def fork_to_product_edges(mod2path: Dict[str, str], edges: Dict[str, Set[str]]) -> List[Tuple[str, str]]:
    """Edges where a bundled fork imports the product (must be none)."""
    out: List[Tuple[str, str]] = []
    for a in edges:
        top = a.split(".")[0]
        if top in BUNDLED_FORKS:
            for b in edges[a]:
                if b == PRODUCT or b.startswith(PRODUCT + "."):
                    out.append((a, b))
    return sorted(out)


def forbidden_layer_edges(edges: Dict[str, Set[str]]) -> List[Tuple[str, str]]:
    """``core -> workflows`` and ``domain -> workflows`` edges (must be none)."""
    out: List[Tuple[str, str]] = []
    for a in edges:
        for b in edges[a]:
            if (a.startswith("mlmm.core") and b.startswith("mlmm.workflows")) or (
                a.startswith("mlmm.domain") and b.startswith("mlmm.workflows")
            ):
                out.append((a, b))
    return sorted(out)


def check_repo(repo_root: Path) -> List[str]:
    """Return a list of human-readable violations (empty == clean)."""
    roots = [(PRODUCT, str(repo_root / PRODUCT))]
    roots += [(f, str(repo_root / f)) for f in BUNDLED_FORKS]
    mod2path, edges = build_graph(roots)

    problems: List[str] = []
    for scc in product_multi_sccs(mod2path, edges):
        cycle_edges = [
            f"{a} -> {b}" for a in scc for b in sorted(edges.get(a, ())) if b in scc
        ]
        problems.append(
            "import cycle among mlmm modules: " + ", ".join(scc)
            + " [" + "; ".join(cycle_edges) + "]"
        )
    for a, b in fork_to_product_edges(mod2path, edges):
        problems.append(f"bundled fork imports the product: {a} -> {b}")
    for a, b in forbidden_layer_edges(edges):
        problems.append(f"forbidden layer edge (core/domain -> workflows): {a} -> {b}")
    return problems


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=str(REPO_ROOT))
    args = parser.parse_args(argv)
    problems = check_repo(Path(args.repo_root))
    if problems:
        print("[import-graph] FAIL:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1
    print("[import-graph] OK: no mlmm cycles, no fork->product edges, "
          "no core/domain->workflows edges.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
