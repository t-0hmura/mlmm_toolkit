"""Shared completion contract for every dry-run-capable workflow."""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.core.utils import DRY_RUN_COMPLETE_MESSAGE, emit_dry_run_complete


EXPECTED = "[Dry run] --dry-run completed. Input command is valid."

_TOY = Path(__file__).resolve().parents[1] / "examples" / "toy_system"
_R, _P, _PARM = (str(_TOY / name) for name in ("r_toy.pdb", "p_toy.pdb", "p_toy.parm7"))
_SYSTEM = ["--parm", _PARM, "-q", "-1", "-m", "1"]
DRY_RUN_COMMANDS = {
    "all": ["all", "-i", _R, "-i", _P, *_SYSTEM],
    "dft": ["dft", "-i", _P, *_SYSTEM],
    "freq": ["freq", "-i", _P, *_SYSTEM],
    "irc": ["irc", "-i", _P, *_SYSTEM],
    "opt": ["opt", "-i", _R, *_SYSTEM],
    "path-opt": ["path-opt", "-i", _R, _P, *_SYSTEM],
    "path-search": ["path-search", "-i", _R, "-i", _P, *_SYSTEM],
    "scan": ["scan", "-i", _R, *_SYSTEM, "--scan-lists", "[(1,5,1.8)]"],
    "scan2d": [
        "scan2d", "-i", _R, *_SYSTEM, "--scan-lists", "[(1,5,1.8,2.2),(1,6,1.7,2.1)]",
    ],
    "scan3d": [
        "scan3d", "-i", _R, *_SYSTEM,
        "--scan-lists", "[(1,5,1.8,2.2),(1,6,1.7,2.1),(2,3,1.0,1.2)]",
    ],
    "sp": ["sp", "-i", _R, *_SYSTEM],
    "tsopt": ["tsopt", "-i", _P, *_SYSTEM],
}
WORKFLOWS_WITH_PLAN_BLOCK = {
    "dft", "freq", "irc", "opt", "path-opt", "path-search",
    "scan", "scan2d", "scan3d", "tsopt",
}
# The block is a titled key/value table; the title also appears in the echoed command line
# when this test is selected by name, so match the underlined title only.
PLAN_BLOCK = re.compile(r"^dry_run_plan\n-+$", re.MULTILINE)

needs_toy_system = pytest.mark.skipif(
    not (_TOY / "p_toy.parm7").exists(), reason="toy_system example inputs not present"
)


def test_canonical_dry_run_footer_is_exact(capsys) -> None:
    assert DRY_RUN_COMPLETE_MESSAGE == EXPECTED
    emit_dry_run_complete()
    assert capsys.readouterr().out.rstrip().splitlines()[-1] == EXPECTED


@needs_toy_system
def test_every_dry_run_workflow_uses_the_canonical_terminal_footer(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    surface = tmp_path / "surface.csv"
    surface.write_text("i,j,k,energy_kcal\n0,0,0,0.0\n", encoding="utf-8")
    commands = {**DRY_RUN_COMMANDS, "scan3d --csv": ["scan3d", "--csv", str(surface)]}
    for name, argv in commands.items():
        result = CliRunner().invoke(root_cli, [*argv, "--dry-run"])
        assert result.exit_code == 0, (name, result.output)
        assert result.output.rstrip().splitlines()[-1] == EXPECTED, name
        assert "[dry-run] Validation complete" not in result.output, name
        assert "Elapsed" not in result.output, name


@needs_toy_system
def test_dry_run_plan_block_is_shown_only_at_verbose_three(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    shown = set()
    for name, argv in DRY_RUN_COMMANDS.items():
        default = CliRunner().invoke(root_cli, [*argv, "--dry-run"])
        full = CliRunner().invoke(root_cli, [*argv, "--dry-run", "-v", "3"])
        assert default.exit_code == 0 and full.exit_code == 0, (name, full.output)
        assert not PLAN_BLOCK.search(default.output), name
        if PLAN_BLOCK.search(full.output):
            shown.add(name)
    assert shown == WORKFLOWS_WITH_PLAN_BLOCK
