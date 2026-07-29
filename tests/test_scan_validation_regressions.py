"""Fail-closed validation for scan axes and scientific eligibility."""

import json
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


def test_scan_spec_expands_bidirectional_stage_with_reset_markers(tmp_path) -> None:
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
            "--no-detect-layer",
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
            "--no-detect-layer",
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
