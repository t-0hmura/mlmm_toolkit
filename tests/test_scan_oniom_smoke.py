"""Lightweight smoke regressions for scan2d/scan3d and ONIOM exporters."""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest
from click.testing import CliRunner

from mlmm.cli import cli as root_cli


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


@pytest.fixture
def nocmap_parm(tmp_path: Path) -> Path:
    pmd = pytest.importorskip("parmed")
    parm = pmd.load_file(str(_fixture("complex.parm7")))
    parm.cmaps[:] = []
    path = tmp_path / "complex_nocmap.parm7"
    parm.save(str(path), overwrite=True)
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
    assert "--input must be a PDB, mmCIF, or XYZ file" in result.output


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
    assert "--input must be a PDB, mmCIF, or XYZ file" in result.output


def test_multidimensional_scans_honor_configured_geometry_coordinates() -> None:
    from mlmm.workflows import scan2d, scan3d

    expected = 'geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])'
    for module in (scan2d, scan3d):
        source = Path(module.__file__).read_text(encoding="utf-8")
        assert expected in source


def test_workflows_guard_optional_companion_conversion() -> None:
    from mlmm.workflows import dft, scan, scan2d, scan3d

    conversion_calls = {
        dft: {"write_ml_region_pdb_pair"},
        scan: {"convert_xyz_to_pdb"},
        scan2d: {"convert_and_annotate_xyz_to_pdb"},
        scan3d: {"convert_and_annotate_xyz_to_pdb"},
    }
    for module, converter_names in conversion_calls.items():
        tree = ast.parse(Path(module.__file__).read_text(encoding="utf-8"))
        parents = {
            child: parent
            for parent in ast.walk(tree)
            for child in ast.iter_child_nodes(parent)
        }
        calls = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id in converter_names
        ]
        assert calls
        for call in calls:
            ancestor = parents.get(call)
            guarded = False
            while ancestor is not None:
                if isinstance(ancestor, ast.If) and any(
                    isinstance(node, ast.Name) and node.id == "convert_files"
                    for node in ast.walk(ancestor.test)
                ):
                    guarded = True
                    break
                ancestor = parents.get(ancestor)
            assert guarded, f"unguarded companion conversion in {module.__name__}:{call.lineno}"


def test_oniom_export_g16_smoke(tmp_path: Path, nocmap_parm: Path) -> None:
    out_file = tmp_path / "model.gjf"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-export",
            "--parm",
            str(nocmap_parm),
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
    assert "MLMM_REF_PDB_ORDER_V1_SHA256=" in text


def test_oniom_export_orca_smoke(tmp_path: Path, nocmap_parm: Path) -> None:
    out_file = tmp_path / "model.inp"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-export",
            "--parm",
            str(nocmap_parm),
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
    assert "# MLMM_REF_PDB_ORDER_V1_SHA256=" in text


def test_oniom_export_mode_inferred_from_output_suffix_g16(
    tmp_path: Path,
    nocmap_parm: Path,
) -> None:
    out_file = tmp_path / "infer_mode.com"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-export",
            "--parm",
            str(nocmap_parm),
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


def test_oniom_export_mode_inferred_from_output_suffix_orca(
    tmp_path: Path,
    nocmap_parm: Path,
) -> None:
    out_file = tmp_path / "infer_mode.inp"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-export",
            "--parm",
            str(nocmap_parm),
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


def test_oniom_export_mode_takes_precedence_over_output_suffix(
    tmp_path: Path,
    nocmap_parm: Path,
) -> None:
    out_file = tmp_path / "forced_g16.inp"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-export",
            "--parm",
            str(nocmap_parm),
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
    nocmap_parm: Path,
    monkeypatch,
) -> None:
    from mlmm.workflows import oniom_export

    monkeypatch.setattr(oniom_export.shutil, "which", lambda _name: None)

    out_file = tmp_path / "model.inp"
    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-export",
            "--parm",
            str(nocmap_parm),
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
    nocmap_parm: Path,
    monkeypatch,
) -> None:
    from mlmm.workflows import oniom_export

    class _FailedProc:
        returncode = 1
        stdout = "simulated orca_mm failure"

    calls = []

    def _run(cmd, **kwargs):
        calls.append((cmd, kwargs))
        return _FailedProc()

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(oniom_export.shutil, "which", lambda _name: "/usr/bin/orca_mm")
    monkeypatch.setattr(oniom_export.subprocess, "run", _run)

    out_file = tmp_path / "distinct-output" / "model.inp"
    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "oniom-export",
            "--parm",
            nocmap_parm.name,
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
    assert calls[0][0][-1] == str(nocmap_parm.resolve())
    assert Path(calls[0][1]["cwd"]) == out_file.parent.resolve()


def test_scan_forces_cartesian_even_with_coord_type_dlc(tmp_path, monkeypatch) -> None:
    """Staged scans run restrained L-BFGS with no microiteration, so DLC over the
    ML/MM system is meaningless (it crashes poly_line_search with a
    Cartesian/internal dimension mismatch). Scan must force Cartesian even when
    ``--coord-type dlc`` is requested, mirroring path-opt / path-search."""
    repo = Path(__file__).resolve().parents[1]
    pdb = repo / "examples" / "toy_system" / "r_complex_layered.pdb"
    parm = repo / "examples" / "toy_system" / "p_complex.parm7"
    if not (pdb.exists() and parm.exists()):
        pytest.skip("toy_system example inputs not present")

    # `scan` gathers --scan-lists from sys.argv (to support multiple lists), so
    # the real argv must be set for the invocation, not just the CliRunner args.
    argv = [
        "mlmm", "scan", "-i", str(pdb), "--parm", str(parm), "-q", "-1", "-m", "1",
        "--scan-lists", "[('PRE 8 C1','PRE 8 C3',2.0)]",
        "--coord-type", "dlc", "--dry-run", "-v", "3",
        "--out-dir", str(tmp_path / "scan_out"),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    result = CliRunner().invoke(root_cli, argv[1:])

    assert result.exit_code == 0, result.output
    # The resolved geometry config in the dry-run plan must be Cartesian, not DLC.
    assert "coord_type: cart" in result.output
    assert "coord_type: dlc" not in result.output
