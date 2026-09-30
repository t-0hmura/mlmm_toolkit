"""Fail-closed validation for scan axes and scientific eligibility."""

import json
import inspect
from pathlib import Path

import click
import numpy as np
import pytest
from click.testing import CliRunner

from mlmm.core.utils import (
    parse_scan_list_quads,
    parse_scan_list_triples,
    parse_scan_spec_stages,
)
from mlmm.workflows.scan2d import _rbf_support
from mlmm.workflows.scan3d import _rbf_support_3d
from mlmm.workflows.scan_common import prepare_grid_scan_output


def _write_scan_structure(path: Path) -> None:
    path.write_text(
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 10.00           C\n"
        "HETATM    2  C2  LIG A   1       1.000   0.000   0.000  1.00 10.00           C\n"
        "HETATM    3  C3  LIG A   1       2.000   0.000   0.000  1.00 10.00           C\n"
        "HETATM    4  C4  LIG A   1       3.000   0.000   0.000  1.00 10.00           C\n"
        "END\n",
        encoding="utf-8",
    )


class _ConstantScanCalculator:
    freeze_atoms = []
    analytical_2d = False

    def get_energy(self, atoms, positions, **kwargs):
        return {"energy": -1.0}

    def get_forces(self, atoms, positions, **kwargs):
        return {"energy": -1.0, "forces": np.zeros(np.asarray(positions).size)}


def _restraint_reaching_lbfgs(offset_angstrom: float = 0.0):
    """Fake optimizer: place each restrained atom at target + offset along x."""
    from types import SimpleNamespace

    from pysisyphus.constants import ANG2BOHR

    def factory(geom, *args, **kwargs):
        def run():
            restraints = getattr(geom.calculator, "_restraints", [])
            coords = np.array(geom.coords3d, dtype=float)
            for first, second, target in restraints:
                coords[second] = coords[first] + np.array(
                    [(float(target) + offset_angstrom) * ANG2BOHR, 0.0, 0.0]
                )
            geom.coords3d = coords

        return SimpleNamespace(run=run, is_converged=True)

    return factory


def _run_grid_scan(tmp_path: Path, monkeypatch, command: str, scan_lists: str, *extra: str):
    from mlmm.cli import cli as root_cli
    from mlmm.workflows import scan2d, scan3d

    module = scan2d if command == "scan2d" else scan3d
    monkeypatch.setattr(module, "mlmm", lambda **kwargs: _ConstantScanCalculator())
    monkeypatch.setattr(scan2d, "write_plotly_image", lambda *args, **kwargs: None)
    structure = tmp_path / "system.pdb"
    _write_scan_structure(structure)
    parm = tmp_path / "system.parm7"
    parm.write_text("not reached\n", encoding="utf-8")
    out_dir = tmp_path / "out"
    result = CliRunner().invoke(
        root_cli,
        [
            command, "-i", str(structure), "--parm", str(parm),
            "--model-indices", "1-4", "--no-detect-layer", "-q", "0", "-m", "1",
            "--scan-lists", scan_lists, "--max-step-size", "0.004",
            "--out-dir", str(out_dir), *extra,
        ],
    )
    return result, out_dir


@pytest.mark.parametrize(
    "raw",
    [
        "[]",
        "[(1, 1, 1.5)]",
        "[(1, 2, 1.5), (2, 1, 1.6)]",
    ],
)
def test_scan_stage_rejects_empty_self_and_duplicate_pairs(raw) -> None:
    with pytest.raises(click.BadParameter):
        parse_scan_list_triples(
            raw,
            one_based=False,
            atom_meta=None,
            option_name="--scan-lists",
        )


@pytest.mark.parametrize(
    "raw",
    [
        "[(1, 1, 1.0, 2.0), (2, 3, 1.0, 2.0)]",
        "[(1, 2, 1.0, 2.0), (2, 1, 1.0, 2.0)]",
    ],
)
def test_grid_scan_rejects_self_and_duplicate_axes(raw) -> None:
    with pytest.raises(click.BadParameter):
        parse_scan_list_quads(
            raw,
            expected_len=2,
            one_based=False,
            atom_meta=None,
            option_name="--scan-lists",
        )


def test_scan_spec_expands_bidirectional_distance_with_reset_markers(tmp_path) -> None:
    spec = tmp_path / "scan.yaml"
    spec.write_text(
        "one_based: true\n"
        "stages:\n"
        "  - [[1, 2, 1.2, 1.8]]\n",
        encoding="utf-8",
    )

    stages, one_based, snapshots, resets = parse_scan_spec_stages(
        spec,
        one_based_default=False,
        atom_meta=None,
        return_bidirectional_markers=True,
    )

    assert one_based is True
    assert stages == [[(0, 1, 1.2)], [(0, 1, 1.8)]]
    assert snapshots == frozenset({0})
    assert resets == frozenset({1})


def test_scan_spec_expands_angle_range_with_reset_markers(tmp_path) -> None:
    spec = tmp_path / "scan.yaml"
    spec.write_text(
        "one_based: true\n"
        "stages:\n"
        "  - [[1, 2, 3, 80.0, 120.0]]\n",
        encoding="utf-8",
    )
    stages, one_based, snapshots, resets = parse_scan_spec_stages(
        spec, one_based_default=False, atom_meta=None,
        return_bidirectional_markers=True,
    )
    assert one_based is True
    assert stages == [[(0, 1, 2, 80.0)], [(0, 1, 2, 120.0)]]
    assert snapshots == frozenset({0})
    assert resets == frozenset({1})


def test_scan_spec_keeps_mixed_target_and_range_order(tmp_path) -> None:
    spec = tmp_path / "scan.yaml"
    spec.write_text(
        "stages:\n"
        "  - [[1, 2, 1.4], [2, 3, 1.2, 1.8]]\n",
        encoding="utf-8",
    )
    stages, _, snapshots, resets = parse_scan_spec_stages(
        spec, one_based_default=True, atom_meta=None,
        return_bidirectional_markers=True,
    )
    assert stages == [[(0, 1, 1.4)], [(1, 2, 1.2)], [(1, 2, 1.8)]]
    assert snapshots == frozenset({1})
    assert resets == frozenset({2})


def test_all_accepts_angle_and_dihedral_scan_targets() -> None:
    from mlmm.workflows.all import _parse_scan_lists_literals

    assert _parse_scan_lists_literals(
        ("[(1,2,3,110.0),(1,2,3,4,-60.0)]",), one_based=True,
    ) == [[(1, 2, 3, 110.0), (1, 2, 3, 4, -60.0)]]


def test_all_rejects_scan_spec_file_cleanly(tmp_path: Path) -> None:
    from mlmm.workflows.all import _parse_scan_lists_literals

    spec = tmp_path / "scan.yaml"
    spec.write_text("stages:\n  - [[1, 2, 1.2]]\n", encoding="utf-8")
    with pytest.raises(click.BadParameter, match="standalone mlmm scan"):
        _parse_scan_lists_literals((str(spec),), one_based=True)


@pytest.mark.parametrize(
    ("x", "y", "expected"),
    [
        ([1.0], [2.0], (1, 0)),
        ([1.0, 2.0, 3.0], [2.0, 2.0, 2.0], (3, 1)),
        ([1.0, 2.0, 1.0], [2.0, 2.0, 3.0], (3, 2)),
    ],
)
def test_scan2d_rbf_support_requires_non_collinear_points(
    x, y, expected,
) -> None:
    assert _rbf_support(np.asarray(x), np.asarray(y)) == expected


@pytest.mark.parametrize(
    ("x", "y", "z", "expected"),
    [
        ([1.0], [2.0], [3.0], (1, 0)),
        ([0.0, 1.0, 0.0, 1.0], [0.0, 0.0, 1.0, 1.0], [2.0] * 4, (4, 2)),
        (
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            (4, 3),
        ),
    ],
)
def test_scan3d_rbf_support_requires_non_coplanar_points(
    x, y, z, expected,
) -> None:
    assert _rbf_support_3d(
        np.asarray(x),
        np.asarray(y),
        np.asarray(z),
    ) == expected


def test_grid_scan_output_replaces_only_current_generation(tmp_path: Path) -> None:
    out_dir = tmp_path / "scan"
    grid = out_dir / "grid"
    grid.mkdir(parents=True)
    (grid / "stale.xyz").write_text("stale\n", encoding="utf-8")
    (grid / "stale.pdb").write_text("stale\n", encoding="utf-8")
    fixed = [
        out_dir / "surface.csv",
        out_dir / "result.json",
        out_dir / "summary.json",
    ]
    for path in fixed:
        path.write_text("stale\n", encoding="utf-8")
    unrelated = out_dir / "notes.txt"
    unrelated.write_text("keep\n", encoding="utf-8")

    resolved, fresh_grid = prepare_grid_scan_output(
        out_dir,
        fixed_names=("surface.csv", "result.json", "summary.json"),
    )

    assert resolved == out_dir.resolve()
    assert fresh_grid == grid.resolve()
    assert list(fresh_grid.iterdir()) == []
    assert all(not path.exists() for path in fixed)
    assert unrelated.read_text(encoding="utf-8") == "keep\n"


@pytest.mark.parametrize(
    ("relative", "message"),
    [
        ("grid/spec.yaml", "reserved grid-scan output"),
        ("surface.csv", "reserved scan output"),
    ],
)
def test_grid_scan_output_rejects_input_collision(
    tmp_path: Path,
    relative: str,
    message: str,
) -> None:
    out_dir = tmp_path / "scan"
    protected = out_dir / relative
    protected.parent.mkdir(parents=True, exist_ok=True)
    protected.write_text("input\n", encoding="utf-8")

    with pytest.raises(click.UsageError, match=message):
        prepare_grid_scan_output(
            out_dir,
            fixed_names=("surface.csv",),
            protected_inputs=(protected,),
        )

    assert protected.read_text(encoding="utf-8") == "input\n"


@pytest.mark.parametrize(
    ("command", "scan_lists"),
    [
        ("scan2d", "[(1,2,1.0,1.2),(3,4,1.0,1.2)]"),
        (
            "scan3d",
            "[(1,2,1.0,1.2),(2,3,1.0,1.2),(3,4,1.0,1.2)]",
        ),
    ],
)
def test_grid_scan_collision_does_not_overwrite_config_input(
    tmp_path: Path,
    command: str,
    scan_lists: str,
) -> None:
    from mlmm.cli import cli as root_cli

    structure = tmp_path / "system.pdb"
    _write_scan_structure(structure)
    parm = tmp_path / "system.parm7"
    parm.write_text("not reached\n", encoding="utf-8")
    out_dir = tmp_path / "scan"
    out_dir.mkdir()
    config = out_dir / "result.json"
    custom = out_dir / "grid" / "custom.py"
    custom.parent.mkdir()
    custom.write_text("calculator = None\n", encoding="utf-8")
    original = json.dumps({"calc": {"calc_file": str(custom)}}) + "\n"
    config.write_text(original, encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            command,
            "-i",
            str(structure),
            "--parm",
            str(parm),
            "--model-pdb",
            str(structure),
            "-q",
            "0",
            "-m",
            "1",
            "--scan-lists",
            scan_lists,
            "--config",
            str(config),
            "--out-dir",
            str(out_dir),
        ],
    )

    assert result.exit_code == 2, result.output
    assert "reserved grid-scan output" in result.output
    assert config.read_text(encoding="utf-8") == original
    assert custom.read_text(encoding="utf-8") == "calculator = None\n"


@pytest.mark.parametrize(
    ("command", "scan_lists"),
    [
        ("scan", "[(1,2,1.2)]"),
        ("scan2d", "[(1,2,1.0,1.2),(3,4,1.0,1.2)]"),
        (
            "scan3d",
            "[(1,2,1.0,1.2),(2,3,1.0,1.2),(3,4,1.0,1.2)]",
        ),
    ],
)
def test_scan_dry_run_honors_configured_layer_detection(
    tmp_path: Path,
    command: str,
    scan_lists: str,
) -> None:
    from mlmm.cli import cli as root_cli

    structure = tmp_path / "system.pdb"
    _write_scan_structure(structure)
    parm = tmp_path / "system.parm7"
    parm.write_text("not reached\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text(
        "calc:\n  use_bfactor_layers: false\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        root_cli,
        [
            command,
            "-i",
            str(structure),
            "--parm",
            str(parm),
            "--model-pdb",
            str(structure),
            "-q",
            "0",
            "-m",
            "1",
            "--scan-lists",
            scan_lists,
            "--config",
            str(config),
            "--out-dir",
            str(tmp_path / "out"),
            "--dry-run",
            "-v",
            "3",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "detect_layer: false" in result.output


def test_scan_collision_does_not_delete_spec_input(tmp_path: Path) -> None:
    from mlmm.cli import cli as root_cli

    structure = tmp_path / "system.pdb"
    _write_scan_structure(structure)
    parm = tmp_path / "system.parm7"
    parm.write_text("not reached\n", encoding="utf-8")
    out_dir = tmp_path / "scan"
    spec = out_dir / "stage_1" / "spec.yaml"
    spec.parent.mkdir(parents=True)
    original = "one_based: true\nstages:\n  - [[1, 2, 1.2]]\n"
    spec.write_text(original, encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            "scan",
            "-i",
            str(structure),
            "--parm",
            str(parm),
            "--model-pdb",
            str(structure),
            "-q",
            "0",
            "-m",
            "1",
            "--scan-lists",
            str(spec),
            "--out-dir",
            str(out_dir),
        ],
    )

    assert result.exit_code == 2, result.output
    assert "overlaps a protected input" in result.output
    assert spec.read_text(encoding="utf-8") == original


def test_scan3d_csv_failure_invalidates_prior_plot(tmp_path: Path) -> None:
    from mlmm.cli import cli as root_cli

    csv_path = tmp_path / "invalid.csv"
    csv_path.write_text("unrelated\n1\n", encoding="utf-8")
    out_dir = tmp_path / "scan"
    out_dir.mkdir()
    stale_plot = out_dir / "scan3d_density.html"
    stale_plot.write_text("stale\n", encoding="utf-8")
    stale_result = out_dir / "result.json"
    stale_summary = out_dir / "summary.json"
    stale_result.write_text('{"old": true}\n', encoding="utf-8")
    stale_summary.write_text('{"old": true}\n', encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            "scan3d",
            "--csv",
            str(csv_path),
            "--out-dir",
            str(out_dir),
            "--no-out-json",
        ],
    )

    assert result.exit_code != 0
    assert not stale_plot.exists()
    assert '"old"' not in stale_result.read_text(encoding="utf-8")
    assert '"old"' not in stale_summary.read_text(encoding="utf-8")


def test_scan_html_plots_are_responsive_and_scan2d_keeps_native_projection() -> None:
    from mlmm.workflows import scan2d, scan3d

    scan2d_source = inspect.getsource(scan2d)
    scan3d_source = inspect.getsource(scan3d)

    assert "plane_proj = go.Surface" in scan2d_source
    assert "plane_z = z_bottom + 0.005" in scan2d_source
    assert "z=np.full_like(ZI, plane_z)" in scan2d_source
    assert "bottom_contours = go.Scatter3d" in scan2d_source
    assert "contour_z = z_bottom + 0.01" in scan2d_source
    assert "go.Figure(data=[surface3d, plane_proj, bottom_contours])" in scan2d_source
    assert 'name="2D Contour Projection (Bottom)"' in scan2d_source
    assert "Computed grid points" not in scan2d_source
    assert '"project": {"z": True}' not in scan2d_source
    for source in (scan2d_source, scan3d_source):
        assert 'config={"responsive": True, "displaylogo": False}' in source
        assert 'default_width="100%"' in source
        assert 'default_height="100%"' in source


def test_scan2d_bottom_contours_are_explicit_line_segments() -> None:
    from mlmm.workflows.scan2d import _contour_line_segments

    x_grid, y_grid = np.meshgrid(
        np.asarray([1.0, 1.5, 2.0]),
        np.asarray([1.0, 1.5, 2.0]),
    )
    z_grid = x_grid + y_grid
    line_x, line_y = _contour_line_segments(
        x_grid, y_grid, z_grid, np.asarray([2.5, 3.0])
    )

    assert line_x
    assert len(line_x) == len(line_y)
    assert any(value is None for value in line_x)
    assert all(value is None or 1.0 <= value <= 2.0 for value in line_x)
    assert all(value is None or 1.0 <= value <= 2.0 for value in line_y)


@pytest.mark.parametrize(
    ("command", "scan_lists", "expected"),
    [
        (
            "scan2d",
            "[(1,2,1.000,1.008),(3,4,1.000,1.004)]",
            [
                "point_i100_j100.xyz",
                "point_i100_j100_grid_000_001.xyz",
                "point_i100_j100_grid_001_000.xyz",
                "point_i101_j100_grid_002_001.xyz",
                "inner_path_d1_002_trj.xyz",
                "inner_path_d1_002.pdb",
            ],
        ),
        (
            "scan3d",
            "[(1,2,1.000,1.008),(2,3,1.000,1.004),(3,4,1.000,1.004)]",
            [
                "point_i100_j100_k100.xyz",
                "point_i100_j100_k100_grid_001_000_000.xyz",
                "point_i101_j100_k100_grid_002_001_001.xyz",
                "inner_path_d1_002_d2_001_trj.xyz",
                "inner_path_d1_002_d2_001.pdb",
            ],
        ),
    ],
)
def test_grid_scan_names_add_grid_indices_only_to_rounded_collisions(
    tmp_path: Path, monkeypatch, command: str, scan_lists: str, expected: list,
) -> None:
    from mlmm.workflows import scan2d, scan3d

    module = scan2d if command == "scan2d" else scan3d
    monkeypatch.setattr(module, "_make_lbfgs", _restraint_reaching_lbfgs())

    result, out_dir = _run_grid_scan(tmp_path, monkeypatch, command, scan_lists, "--dump")

    assert result.exit_code == 0, result.output
    names = {path.name for path in (out_dir / "grid").iterdir()}
    n_points = 6 if command == "scan2d" else 12
    assert len([name for name in names if name.startswith("point_") and name.endswith(".xyz")]) == n_points
    assert set(expected) <= names


@pytest.mark.parametrize(
    ("command", "scan_lists", "columns"),
    [
        (
            "scan2d",
            "[(1,2,1.000,1.008),(3,4,1.000,1.004)]",
            [
                "i", "j", "d1_A", "d2_A", "energy_hartree", "bias_converged",
                "is_preopt", "target_d1_A", "target_d2_A", "energy_kcal",
                "d1_label", "d2_label", "q1", "q1_unit", "target_q1",
                "q2", "q2_unit", "target_q2",
            ],
        ),
        (
            "scan3d",
            "[(1,2,1.000,1.008),(2,3,1.000,1.004),(3,4,1.000,1.004)]",
            [
                "i", "j", "k", "d1_A", "d2_A", "d3_A", "target_d1_A",
                "target_d2_A", "target_d3_A", "energy_hartree",
                "bias_converged", "is_preopt", "energy_kcal", "d1_label",
                "d2_label", "d3_label", "q1", "q1_unit", "target_q1",
                "q2", "q2_unit", "target_q2", "q3", "q3_unit", "target_q3",
            ],
        ),
    ],
)
def test_grid_scan_surface_records_measured_and_target_coordinates(
    tmp_path: Path, monkeypatch, command: str, scan_lists: str, columns: list,
) -> None:
    import pandas as pd

    from mlmm.workflows import scan2d, scan3d

    offset = 0.002
    module = scan2d if command == "scan2d" else scan3d
    monkeypatch.setattr(module, "_make_lbfgs", _restraint_reaching_lbfgs(offset))

    result, out_dir = _run_grid_scan(
        tmp_path, monkeypatch, command, scan_lists, "--out-json"
    )

    assert result.exit_code == 0, result.output
    surface = pd.read_csv(out_dir / "surface.csv")
    assert list(surface.columns) == columns
    grid = surface[surface["i"] >= 0]
    n_axes = 2 if command == "scan2d" else 3
    for axis in range(1, n_axes + 1):
        measured = grid[f"d{axis}_A"].to_numpy(dtype=float)
        target = grid[f"target_d{axis}_A"].to_numpy(dtype=float)
        np.testing.assert_allclose(measured, target + offset, atol=1e-9)
        np.testing.assert_allclose(grid[f"q{axis}"], measured)
        np.testing.assert_allclose(grid[f"target_q{axis}"], target)
    np.testing.assert_allclose(
        sorted(grid["target_d1_A"].unique()), [1.0, 1.004, 1.008]
    )
    payload = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    for point in payload["grid_points"]:
        np.testing.assert_allclose(
            point["coordinate_values"],
            np.asarray(point["coordinate_targets"]) + offset,
            atol=1e-9,
        )
        assert point["distances_angstrom"] == point["coordinate_values"]
        assert point["targets_angstrom"] == point["coordinate_targets"]


@pytest.mark.parametrize(
    ("command", "scan_lists"),
    [
        ("scan2d", "[(1,2,1.000,1.000),(3,4,1.000,1.004)]"),
        ("scan3d", "[(1,2,1.000,1.000),(2,3,1.000,1.000),(3,4,1.000,1.000)]"),
    ],
)
def test_grid_scan_insufficient_plot_data_skips_plots_without_error_json(
    tmp_path: Path, monkeypatch, command: str, scan_lists: str,
) -> None:
    from mlmm.workflows import scan2d, scan3d

    module = scan2d if command == "scan2d" else scan3d
    monkeypatch.setattr(module, "_make_lbfgs", _restraint_reaching_lbfgs())

    result, out_dir = _run_grid_scan(
        tmp_path, monkeypatch, command, scan_lists, "--no-out-json"
    )

    assert result.exit_code == 0, result.output
    assert "[plot] NOTE:" in result.output
    assert "[plot] ERROR:" not in result.output
    assert "Traceback" not in result.output
    assert (out_dir / "surface.csv").exists()
    assert not (out_dir / "result.json").exists()
    assert not list(out_dir.glob("*.html"))
