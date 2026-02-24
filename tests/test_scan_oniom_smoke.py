"""Lightweight smoke regressions for scan2d/scan3d and ONIOM exporters."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from mlmm_toolkit.cli import cli as root_cli


_FIXTURE_DIR = (
    Path(__file__).resolve().parents[1]
    / "hessian_ff"
    / "tests"
    / "data"
    / "small"
)


def _fixture(name: str) -> Path:
    path = _FIXTURE_DIR / name
    assert path.exists(), f"missing fixture: {path}"
    return path


def test_scan2d_rejects_non_pdb_xyz_input(tmp_path: Path) -> None:
    bad_input = tmp_path / "input.txt"
    bad_input.write_text("dummy\n", encoding="utf-8")

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "scan2d",
            "-i",
            str(bad_input),
            "--real-parm7",
            str(_fixture("complex.parm7")),
            "-q",
            "0",
            "--scan-lists",
            "[(1,2,1.0,2.0),(3,4,1.0,2.0)]",
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 1
    assert "--input must be a PDB or XYZ file" in result.output


def test_scan3d_rejects_non_pdb_xyz_input(tmp_path: Path) -> None:
    bad_input = tmp_path / "input.txt"
    bad_input.write_text("dummy\n", encoding="utf-8")

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "scan3d",
            "-i",
            str(bad_input),
            "--real-parm7",
            str(_fixture("complex.parm7")),
            "-q",
            "0",
            "--scan-lists",
            "[(1,2,1.0,2.0),(3,4,1.0,2.0),(5,6,1.0,2.0)]",
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 1
    assert "--input must be a PDB or XYZ file" in result.output


def test_oniom_gaussian_export_smoke(tmp_path: Path) -> None:
    out_file = tmp_path / "model.gjf"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-gaussian",
            "--parm7",
            str(_fixture("complex.parm7")),
            "-i",
            str(_fixture("complex.pdb")),
            "--model-pdb",
            str(_fixture("complex.pdb")),
            "--no-element-check",
            "-o",
            str(out_file),
            "-q",
            "0",
            "-m",
            "1",
            "--near",
            "4.0",
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert out_file.exists()
    text = out_file.read_text(encoding="utf-8")
    assert "ONIOM" in text
    assert "0 1" in text


def test_oniom_orca_export_smoke(tmp_path: Path) -> None:
    out_file = tmp_path / "model.inp"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-orca",
            "--parm7",
            str(_fixture("complex.parm7")),
            "-i",
            str(_fixture("complex.pdb")),
            "--model-pdb",
            str(_fixture("complex.pdb")),
            "--no-element-check",
            "-o",
            str(out_file),
            "-q",
            "0",
            "-m",
            "1",
            "--near",
            "4.0",
            "--no-convert-orcaff",
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert out_file.exists()
    text = out_file.read_text(encoding="utf-8")
    assert "QMMM" in text
    assert "ORCAFFFilename" in text
