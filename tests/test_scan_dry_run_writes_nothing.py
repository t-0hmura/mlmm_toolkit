"""The scan-family --dry-run leaves the output directory untouched."""

from __future__ import annotations

import importlib
import tempfile
from pathlib import Path

import pytest
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.core.utils import DRY_RUN_COMPLETE_MESSAGE


SCAN_LISTS = {
    "scan": "[(1,2,1.2)]",
    "scan2d": "[(1,2,1.0,1.2),(3,4,1.0,1.2)]",
    "scan3d": "[(1,2,1.0,1.2),(2,3,1.0,1.2),(3,4,1.0,1.2)]",
}
DRY_RUN_HELP = (
    "Resolve and validate options (input, charge/spin, --scan-lists parse) "
    "and print the planned scan, then exit without running any optimization."
)


def _write_layered_structure(path: Path) -> None:
    # Atoms 1-4 are ML (B = 0) and atoms 5-6 are movable MM (B = 10).
    lines = [
        f"HETATM{i:5d}  C{i:<2d} LIG A   1    {float(i - 1):8.3f}   0.000   0.000"
        f"  1.00{0.0 if i <= 4 else 10.0:6.2f}           C\n"
        for i in range(1, 7)
    ]
    path.write_text("".join(lines) + "END\n", encoding="utf-8")


def _invoke_dry_run(tmp_path: Path, command: str, *extra: str):
    structure = tmp_path / "system.pdb"
    _write_layered_structure(structure)
    # The dry-run never reads the parm7, so a placeholder file is enough.
    parm = tmp_path / "system.parm7"
    parm.write_text("not reached\n", encoding="utf-8")
    out_dir = tmp_path / "out"
    result = CliRunner().invoke(
        root_cli,
        [
            command, "-i", str(structure), "--parm", str(parm),
            "-q", "0", "-m", "1",
            "--scan-lists", SCAN_LISTS[command],
            "--out-dir", str(out_dir), "--dry-run", *extra,
        ],
    )
    return result, out_dir, structure


@pytest.mark.parametrize("command", sorted(SCAN_LISTS))
@pytest.mark.parametrize("layer_mode", ["bfactor", "model-indices", "model-pdb"])
def test_scan_dry_run_creates_nothing_under_out_dir(
    tmp_path: Path, monkeypatch, command: str, layer_mode: str
) -> None:
    # Import first: third-party packages create temp dirs of their own on import.
    importlib.import_module(f"mlmm.workflows.{command}")
    scratch = tmp_path / "tmp"
    scratch.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(scratch))
    extra = {
        "bfactor": (),
        "model-indices": ("--model-indices", "1-4", "--no-detect-layer"),
        "model-pdb": ("--model-pdb", str(tmp_path / "system.pdb")),
    }[layer_mode]

    result, out_dir, _ = _invoke_dry_run(tmp_path, command, *extra)

    assert result.exit_code == 0, result.output
    assert result.output.rstrip().splitlines()[-1] == DRY_RUN_COMPLETE_MESSAGE
    assert not out_dir.exists()
    # Layer files written for validation are removed again.
    assert list(scratch.iterdir()) == []


@pytest.mark.parametrize("command", sorted(SCAN_LISTS))
def test_scan_dry_run_prints_only_the_validation_summary(
    tmp_path: Path, command: str
) -> None:
    result, out_dir, _ = _invoke_dry_run(tmp_path, command)

    assert result.exit_code == 0, result.output
    assert (
        f"[{command}] --dry-run: input, charge/spin, and --scan-lists parse OK."
        in result.output
    )
    assert "scan-parsed" not in result.output
    assert result.output.rstrip().splitlines()[-1] == DRY_RUN_COMPLETE_MESSAGE
    assert not out_dir.exists()


@pytest.mark.parametrize("command", sorted(SCAN_LISTS))
def test_scan_dry_run_help_matches_the_summary_it_prints(command: str) -> None:
    module = importlib.import_module(f"mlmm.workflows.{command}")
    option = next(param for param in module.cli.params if param.name == "dry_run")

    assert option.help == DRY_RUN_HELP


@pytest.mark.parametrize("command", sorted(SCAN_LISTS))
def test_scan_dry_run_rejects_an_input_at_the_generated_model_path(
    tmp_path: Path, command: str
) -> None:
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    structure = out_dir / "model_from_bfactor.pdb"
    _write_layered_structure(structure)
    parm = tmp_path / "system.parm7"
    parm.write_text("not reached\n", encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            command, "-i", str(structure), "--parm", str(parm),
            "-q", "0", "-m", "1",
            "--scan-lists", SCAN_LISTS[command],
            "--out-dir", str(out_dir), "--dry-run",
        ],
    )

    assert result.exit_code != 0, result.output
    assert "collides with generated ML-region model path" in result.output
