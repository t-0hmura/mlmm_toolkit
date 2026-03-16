#!/usr/bin/env python3
"""Verify that options forwarded from `mlmm all` to `mlmm scan` exist on scan CLI."""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Set


REPO_ROOT = Path(__file__).resolve().parents[1]
ALL_PY = REPO_ROOT / "mlmm" / "all.py"
SCAN_PY = REPO_ROOT / "mlmm" / "scan.py"


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
    def __init__(self) -> None:
        self.in_cli = False
        self.flags: Set[str] = set()

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        was_in_cli = self.in_cli
        if node.name == "cli":
            self.in_cli = True
        self.generic_visit(node)
        self.in_cli = was_in_cli

    def visit_Assign(self, node: ast.Assign) -> None:
        if self.in_cli:
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "scan_args":
                    if isinstance(node.value, ast.List):
                        for elt in node.value.elts:
                            flag = _const_flag(elt)
                            if flag is not None:
                                self.flags.add(flag)
        self.generic_visit(node)

    def visit_Call(self, node: ast.Call) -> None:
        if self.in_cli:
            if isinstance(node.func, ast.Name) and node.func.id == "_append_cli_arg":
                if len(node.args) >= 2 and isinstance(node.args[0], ast.Name) and node.args[0].id == "scan_args":
                    flag = _const_flag(node.args[1])
                    if flag is not None:
                        self.flags.add(flag)

            if isinstance(node.func, ast.Attribute) and node.func.attr == "append":
                if isinstance(node.func.value, ast.Name) and node.func.value.id == "scan_args" and node.args:
                    value = node.args[0]
                    flag = _const_flag(value)
                    if flag is not None:
                        self.flags.add(flag)
                    elif isinstance(value, ast.IfExp):
                        body_flag = _const_flag(value.body)
                        else_flag = _const_flag(value.orelse)
                        if body_flag is not None:
                            self.flags.add(body_flag)
                        if else_flag is not None:
                            self.flags.add(else_flag)
        self.generic_visit(node)


def main() -> int:
    all_tree = ast.parse(ALL_PY.read_text(encoding="utf-8"), filename=str(ALL_PY))
    scan_tree = ast.parse(SCAN_PY.read_text(encoding="utf-8"), filename=str(SCAN_PY))

    declared = _collect_scan_declared_flags(scan_tree)
    collector = _ForwardedScanFlags()
    collector.visit(all_tree)
    forwarded = collector.flags

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
