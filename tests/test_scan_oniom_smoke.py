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
            "--parm",
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
            "--parm",
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


def test_oniom_export_g16_smoke(tmp_path: Path) -> None:
    out_file = tmp_path / "model.gjf"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-export",
            "--parm",
            str(_fixture("complex.parm7")),
            "-i",
            str(_fixture("complex.pdb")),
            "--model-pdb",
            str(_fixture("complex.pdb")),
            "--no-element-check",
            "--mode",
            "g16",
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


def test_oniom_export_orca_smoke(tmp_path: Path) -> None:
    out_file = tmp_path / "model.inp"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-export",
            "--parm",
            str(_fixture("complex.parm7")),
            "-i",
            str(_fixture("complex.pdb")),
            "--model-pdb",
            str(_fixture("complex.pdb")),
            "--no-element-check",
            "--mode",
            "orca",
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


def test_oniom_export_mode_inferred_from_output_suffix_g16(tmp_path: Path) -> None:
    out_file = tmp_path / "infer_mode.com"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-export",
            "--parm",
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
    text = out_file.read_text(encoding="utf-8")
    assert "ONIOM" in text


def test_oniom_export_mode_inferred_from_output_suffix_orca(tmp_path: Path) -> None:
    out_file = tmp_path / "infer_mode.inp"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-export",
            "--parm",
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
    text = out_file.read_text(encoding="utf-8")
    assert "QMMM" in text


def test_oniom_export_mode_takes_precedence_over_output_suffix(tmp_path: Path) -> None:
    out_file = tmp_path / "forced_g16.inp"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-export",
            "--parm",
            str(_fixture("complex.parm7")),
            "-i",
            str(_fixture("complex.pdb")),
            "--model-pdb",
            str(_fixture("complex.pdb")),
            "--no-element-check",
            "--mode",
            "g16",
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
    text = out_file.read_text(encoding="utf-8")
    assert "ONIOM" in text


def test_oniom_export_errors_when_mode_and_suffix_unknown(tmp_path: Path) -> None:
    out_file = tmp_path / "unknown.ext"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-export",
            "--parm",
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

    assert result.exit_code == 1
    assert "Could not infer export mode from -o/--output" in result.output
    assert "Specify --mode (g16/orca)" in result.output


def test_oniom_export_help_shows_convert_orcaff_default_enabled() -> None:
    runner = CliRunner()
    result = runner.invoke(root_cli, ["oniom-export", "--help-advanced"], catch_exceptions=False)

    assert result.exit_code == 0, result.output
    assert "--convert-orcaff / --no-convert-orcaff" in result.output
    assert "[default: convert-orcaff]" in result.output


def test_legacy_oniom_subcommands_removed() -> None:
    runner = CliRunner()

    g16 = runner.invoke(root_cli, ["oniom-gaussian", "--help"], catch_exceptions=False)
    assert g16.exit_code != 0
    assert "No such command 'oniom-gaussian'" in g16.output

    orca = runner.invoke(root_cli, ["oniom-orca", "--help"], catch_exceptions=False)
    assert orca.exit_code != 0
    assert "No such command 'oniom-orca'" in orca.output


def test_oniom_export_default_convert_fallback_without_orca_mm(
    tmp_path: Path,
    monkeypatch,
) -> None:
    from mlmm_toolkit import oniom_export

    monkeypatch.setattr(oniom_export.shutil, "which", lambda _name: None)

    out_file = tmp_path / "model.inp"
    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-export",
            "--parm",
            str(_fixture("complex.parm7")),
            "-i",
            str(_fixture("complex.pdb")),
            "--model-pdb",
            str(_fixture("complex.pdb")),
            "--no-element-check",
            "--mode",
            "orca",
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
    assert "Run manually: cd " in result.output
    assert "orca_mm -convff -AMBER" in result.output


def test_oniom_export_orca_mm_failure_prints_manual_command(
    tmp_path: Path,
    monkeypatch,
) -> None:
    from mlmm_toolkit import oniom_export

    class _FailedProc:
        returncode = 1
        stdout = "simulated orca_mm failure"

    monkeypatch.setattr(oniom_export.shutil, "which", lambda _name: "/usr/bin/orca_mm")
    monkeypatch.setattr(oniom_export.subprocess, "run", lambda *args, **kwargs: _FailedProc())

    out_file = tmp_path / "model.inp"
    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-export",
            "--parm",
            str(_fixture("complex.parm7")),
            "-i",
            str(_fixture("complex.pdb")),
            "--model-pdb",
            str(_fixture("complex.pdb")),
            "--no-element-check",
            "--mode",
            "orca",
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

    assert result.exit_code == 1
    assert "orca_mm failed" in result.output
    assert "Run manually: cd " in result.output
    assert "orca_mm -convff -AMBER" in result.output
