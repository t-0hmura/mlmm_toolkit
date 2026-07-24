from __future__ import annotations

import ast
import importlib.util
from pathlib import Path


def _load_checker():
    path = (
        Path(__file__).parents[1]
        / ".github"
        / "scripts"
        / "check_all_scan_contract.py"
    )
    spec = importlib.util.spec_from_file_location("mlmm_scan_contract", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_checker_collects_ifexp_extend_and_named_helpers() -> None:
    module = _load_checker()
    tree = ast.parse(
        """
def cli():
    scan_args = ["--positive" if enabled else "--negative"]
    scan_args.extend(["--extended", value])

def append_backend_forwarding_args(args):
    args.extend(["--from-helper", value])

def build_scan_child_argv():
    return helper((("threshold", "--thresh", value, False),))
"""
    )
    collector = module._ForwardedScanFlags({"cli": {"scan_args"}})
    collector.visit(tree)
    helper = module._ForwardedScanFlags(
        {"append_backend_forwarding_args": {"args"}},
    )
    helper.visit(tree)
    flags = (
        collector.flags
        | helper.flags
        | module._collect_named_function_flags(tree, {"build_scan_child_argv"})
    )
    assert flags == {
        "--extended",
        "--from-helper",
        "--negative",
        "--positive",
        "--thresh",
    }
