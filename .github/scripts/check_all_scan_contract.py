#!/usr/bin/env python3
"""Verify that options forwarded from `mlmm all` to `mlmm scan` exist on scan CLI."""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Set


REPO_ROOT = Path(__file__).resolve().parents[2]
ALL_PY = REPO_ROOT / "mlmm" / "workflows" / "all.py"
ALL_HELPERS_PY = REPO_ROOT / "mlmm" / "workflows" / "_all_helpers.py"
SCAN_PY = REPO_ROOT / "mlmm" / "workflows" / "scan.py"
# Shared option decorators (e.g. ``add_ml_layer_detection_options``) inject
# Click flags onto ``scan.cli`` from here, so they count as declared too.
COMMON_OPTIONS_PY = REPO_ROOT / "mlmm" / "cli" / "common_options.py"


def _is_click_option_call(node: ast.Call) -> bool:
    func = node.func
    return (
        isinstance(func, ast.Attribute)
        and func.attr == "option"
        and isinstance(func.value, ast.Name)
        and func.value.id == "click"
    )


def _const_flag(node: ast.AST) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        value = node.value.strip()
        if value.startswith("--"):
            return value
    return None


def _const_flags(node: ast.AST) -> Set[str]:
    """Return every literal long option contained in an expression."""
    flags: Set[str] = set()
    for child in ast.walk(node):
        flag = _const_flag(child)
        if flag is not None:
            flags.add(flag)
    return flags


def _expand_flag_spec(spec: str) -> Set[str]:
    # click toggle declarations may appear as "--a/--b".
    parts = [p.strip() for p in spec.split("/") if p.strip()]
    flags = {p for p in parts if p.startswith("--")}
    if flags:
        return flags
    return {spec}


def _collect_scan_declared_flags(scan_tree: ast.AST) -> Set[str]:
    flags: Set[str] = set()
    for node in ast.walk(scan_tree):
        if not isinstance(node, ast.Call) or not _is_click_option_call(node):
            continue
        for arg in node.args:
            flag = _const_flag(arg)
            if flag is not None:
                flags.update(_expand_flag_spec(flag))
    return flags


class _ForwardedScanFlags(ast.NodeVisitor):
    def __init__(self, function_targets: dict[str, Set[str]]) -> None:
        self.function_targets = function_targets
        self.targets: Set[str] = set()
        self.flags: Set[str] = set()

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        previous = self.targets
        self.targets = self.function_targets.get(node.name, set())
        self.generic_visit(node)
        self.targets = previous

    def visit_Assign(self, node: ast.Assign) -> None:
        if self.targets:
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in self.targets:
                    self.flags.update(_const_flags(node.value))
        self.generic_visit(node)

    def visit_Call(self, node: ast.Call) -> None:
        if self.targets:
            if isinstance(node.func, ast.Name) and node.func.id == "_append_cli_arg":
                if (
                    len(node.args) >= 2
                    and isinstance(node.args[0], ast.Name)
                    and node.args[0].id in self.targets
                ):
                    self.flags.update(_const_flags(node.args[1]))

            if (
                isinstance(node.func, ast.Attribute)
                and node.func.attr in {"append", "extend"}
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id in self.targets
            ):
                for arg in node.args:
                    self.flags.update(_const_flags(arg))
        self.generic_visit(node)


def _collect_named_function_flags(
    tree: ast.AST, function_names: Set[str],
) -> Set[str]:
    flags: Set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name in function_names:
            flags.update(_const_flags(node))
    return flags


def main() -> int:
    all_tree = ast.parse(ALL_PY.read_text(encoding="utf-8"), filename=str(ALL_PY))
    helpers_tree = ast.parse(
        ALL_HELPERS_PY.read_text(encoding="utf-8"),
        filename=str(ALL_HELPERS_PY),
    )
    scan_tree = ast.parse(SCAN_PY.read_text(encoding="utf-8"), filename=str(SCAN_PY))

    declared = _collect_scan_declared_flags(scan_tree)
    # scan.cli also pulls in flags from shared option decorators defined in
    # mlmm/cli/common_options.py (--detect-layer, --precision, ...).
    if COMMON_OPTIONS_PY.exists():
        common_tree = ast.parse(
            COMMON_OPTIONS_PY.read_text(encoding="utf-8"), filename=str(COMMON_OPTIONS_PY)
        )
        declared |= _collect_scan_declared_flags(common_tree)
    collector = _ForwardedScanFlags({"cli": {"scan_args"}})
    collector.visit(all_tree)
    helper_collector = _ForwardedScanFlags(
        {"append_backend_forwarding_args": {"args"}},
    )
    helper_collector.visit(helpers_tree)
    forwarded = collector.flags | helper_collector.flags
    forwarded |= _collect_named_function_flags(
        helpers_tree, {"build_scan_child_argv"},
    )

    missing = sorted(flag for flag in forwarded if flag not in declared)
    if missing:
        print("Detected all->scan option contract drift.")
        for flag in missing:
            print(f"  missing in scan CLI: {flag}")
        return 1

    print(f"[scan-contract] all forwards {len(forwarded)} scan options; all are declared.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
