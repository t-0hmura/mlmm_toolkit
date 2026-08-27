"""Shared completion contract for every dry-run-capable workflow."""

from __future__ import annotations

import ast
from pathlib import Path

from mlmm.core.utils import DRY_RUN_COMPLETE_MESSAGE, emit_dry_run_complete


WORKFLOWS = (
    "all", "dft", "freq", "irc", "opt", "path_opt", "path_search",
    "scan", "scan2d", "scan3d", "sp", "tsopt",
)
EXPECTED_CALLS = {name: 1 for name in WORKFLOWS} | {
    "scan": 2, "scan2d": 2, "scan3d": 2,
}
EXPECTED = "[Dry run] --dry-run completed. Input command is valid."


def _call_names(nodes: list[ast.stmt]) -> set[str]:
    names: set[str] = set()
    for root in nodes:
        for node in ast.walk(root):
            if not isinstance(node, ast.Call):
                continue
            if isinstance(node.func, ast.Name):
                names.add(node.func.id)
            elif isinstance(node.func, ast.Attribute):
                names.add(node.func.attr)
    return names


def _is_terminal(nodes: list[ast.stmt]) -> bool:
    return any(
        isinstance(node, ast.Return)
        or (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "sys"
            and node.func.attr == "exit"
        )
        for root in nodes
        for node in ast.walk(root)
    )


def test_canonical_dry_run_footer_is_exact(capsys) -> None:
    assert DRY_RUN_COMPLETE_MESSAGE == EXPECTED
    emit_dry_run_complete()
    assert capsys.readouterr().out.rstrip().splitlines()[-1] == EXPECTED


def test_every_dry_run_workflow_uses_the_canonical_terminal_footer() -> None:
    workflows = Path(__file__).resolve().parents[1] / "mlmm" / "workflows"
    for name in WORKFLOWS:
        path = workflows / f"{name}.py"
        source = path.read_text(encoding="utf-8")
        assert '"--dry-run/--no-dry-run"' in source, name
        assert "emit_dry_run_complete" in source, name
        assert source.count("emit_dry_run_complete()") == EXPECTED_CALLS[name], name
        assert "[dry-run] Validation complete" not in source, name

        tree = ast.parse(source, filename=str(path))
        terminals = [
            node for node in ast.walk(tree)
            if isinstance(node, ast.If)
            and isinstance(node.test, ast.Name)
            and node.test.id == "dry_run"
            and _is_terminal(node.body)
        ]
        assert terminals, name
        for terminal in terminals:
            calls = _call_names(terminal.body)
            assert "format_elapsed" not in calls, name
            assert calls & {"emit_dry_run_complete", "_emit_final_summary"}, name
