"""The scan-family --dry-run prints a short input summary at the default verbosity."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.core.utils import DRY_RUN_COMPLETE_MESSAGE


def _write_scan_structure(path: Path) -> None:
    path.write_text(
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 10.00           C\n"
        "HETATM    2  C2  LIG A   1       1.000   0.000   0.000  1.00 10.00           C\n"
        "HETATM    3  C3  LIG A   1       2.000   0.000   0.000  1.00 10.00           C\n"
        "HETATM    4  C4  LIG A   1       3.000   0.000   0.000  1.00 10.00           C\n"
        "END\n",
        encoding="utf-8",
    )


@pytest.mark.parametrize(
    ("command", "scan_lists", "parsed", "flags", "not_run"),
    [
        (
            "scan",
            "[(1,2,1.2)]",
            "('[(1,2,1.2)]',) → 1 stage(s)",
            "preopt=False  endopt=False",
            "No scan / preopt was executed.",
        ),
        (
            "scan",
            "[(1,2,1.0,1.2)]",
            "('[(1,2,1.0,1.2)]',) → 2 stage(s)",
            "preopt=False  endopt=False",
            "No scan / preopt was executed.",
        ),
        (
            "scan2d",
            "[(1,2,1.0,1.2),(3,4,1.0,1.2)]",
            "[(1,2,1.0,1.2),(3,4,1.0,1.2)] → 2 axis tuples",
            "preopt=False",
            "No 2D scan was executed.",
        ),
        (
            "scan3d",
            "[(1,2,1.0,1.2),(2,3,1.0,1.2),(3,4,1.0,1.2)]",
            "[(1,2,1.0,1.2),(2,3,1.0,1.2),(3,4,1.0,1.2)] → 3 axis tuples",
            "preopt=False",
            "No 3D scan was executed.",
        ),
    ],
)
def test_scan_dry_run_prints_summary_at_default_verbosity(
    tmp_path: Path,
    command: str,
    scan_lists: str,
    parsed: str,
    flags: str,
    not_run: str,
) -> None:
    structure = tmp_path / "system.pdb"
    _write_scan_structure(structure)
    # The dry-run never reads the parm7, so a placeholder file is enough.
    parm = tmp_path / "system.parm7"
    parm.write_text("not reached\n", encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            command,
            "-i", str(structure),
            "--parm", str(parm),
            "--model-pdb", str(structure),
            "-q", "0",
            "-m", "1",
            "--scan-lists", scan_lists,
            "--out-dir", str(tmp_path / "out"),
            "--dry-run",
        ],
    )

    assert result.exit_code == 0, result.output
    tail = result.output.rstrip().splitlines()[-9:]
    tag = f"[{command}]"
    assert tail[0] == f"{tag} --dry-run: input, charge/spin, and --scan-lists parse OK."
    assert tail[1].startswith(f"{tag} input geometry  : ")
    assert tail[2] == f"{tag} resolved charge : +0"
    assert tail[3] == f"{tag} resolved spin   : 1 (multiplicity)"
    assert tail[4].startswith(f"{tag} out_dir         : ")
    assert tail[5] == f"{tag} --scan-lists    : {parsed}"
    assert tail[6] == f"{tag} {flags}"
    assert tail[7] == f"{tag} {not_run}"
    assert tail[8] == DRY_RUN_COMPLETE_MESSAGE
    # The full plan block stays behind -v 3.
    assert "dry_run_plan" not in result.output


def test_scan3d_csv_dry_run_only_checks_options(tmp_path: Path) -> None:
    csv_path = tmp_path / "surface.csv"
    csv_path.write_text(
        "i,j,k,d1_A,d2_A,d3_A,energy_hartree\n"
        + "".join(
            f"{i},{j},{k},{1.0 + i},{1.5 + j},{2.0 + k},{-(i + j + k) / 1000}\n"
            for i in (0, 1)
            for j in (0, 1)
            for k in (0, 1)
        ),
        encoding="utf-8",
    )
    out_dir = tmp_path / "scan"
    out_dir.mkdir()
    stale_plot = out_dir / "scan3d_density.html"
    stale_plot.write_text("stale\n", encoding="utf-8")
    stale_result = out_dir / "result.json"
    stale_result.write_text('{"old": true}\n', encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        ["scan3d", "--csv", str(csv_path), "--out-dir", str(out_dir), "--dry-run"],
    )

    assert result.exit_code == 0, result.output
    tail = result.output.rstrip().splitlines()[-5:]
    assert tail[0] == "[scan3d] --dry-run with --csv: option parsing OK."
    assert tail[1].startswith("[scan3d] csv input  : ")
    assert tail[2].startswith("[scan3d] out_dir    : ")
    assert tail[3] == "[scan3d] No 3D scan was executed."
    assert tail[4] == DRY_RUN_COMPLETE_MESSAGE
    # Earlier outputs in the directory are left untouched.
    assert sorted(p.name for p in out_dir.iterdir()) == ["result.json", "scan3d_density.html"]
    assert stale_plot.read_text(encoding="utf-8") == "stale\n"
    assert stale_result.read_text(encoding="utf-8") == '{"old": true}\n'
