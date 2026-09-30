"""Grid scans skip only the plots when usable points cannot span the surface."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from mlmm.workflows.scan2d import _rbf_support
from mlmm.workflows.scan3d import _rbf_support_3d


def test_rank_check_ignores_rounding_in_collinear_and_coplanar_points() -> None:
    # 0.1 + 0.2 differs from 0.3 by one ulp; the points are still on a line / plane.
    assert _rbf_support(
        np.asarray([1.0, 1.004, 1.008]), np.asarray([0.1 + 0.2, 0.3, 0.3])
    ) == (3, 1)
    assert _rbf_support_3d(
        np.asarray([1.0, 1.004, 1.0, 1.004]),
        np.asarray([1.0, 1.0, 1.004, 1.004]),
        np.asarray([0.3, 0.3, 0.3, 0.1 + 0.2]),
    ) == (4, 2)


@pytest.mark.parametrize(
    ("origin", "step"),
    [
        ((1.0, 1.0, 1.0), 0.004),
        ((179.0, -179.0, 120.0), 0.5),
        ((2.0, 3.0, 4.0), 1e-6),
    ],
)
def test_rank_check_keeps_small_and_large_valued_grids_full_rank(origin, step) -> None:
    x0, y0, z0 = origin
    x = np.asarray([x0, x0 + step, x0, x0])
    y = np.asarray([y0, y0, y0 + step, y0])
    z = np.asarray([z0, z0, z0, z0 + step])
    assert _rbf_support(x[:3], y[:3]) == (3, 2)
    assert _rbf_support_3d(x, y, z) == (4, 3)


class _ConstantScanCalculator:
    freeze_atoms = []
    analytical_2d = False

    def get_energy(self, atoms, positions, **kwargs):
        return {"energy": -1.0}

    def get_forces(self, atoms, positions, **kwargs):
        return {"energy": -1.0, "forces": np.zeros(np.asarray(positions).size)}


def _restraint_reaching_lbfgs(converged: bool = True):
    """Fake optimizer: place each restrained atom at its target along x."""
    from types import SimpleNamespace

    from pysisyphus.constants import ANG2BOHR

    def factory(geom, *args, **kwargs):
        def run():
            coords = np.array(geom.coords3d, dtype=float)
            for first, second, target in getattr(geom.calculator, "_restraints", []):
                coords[second] = coords[first] + np.array(
                    [float(target) * ANG2BOHR, 0.0, 0.0]
                )
            geom.coords3d = coords

        return SimpleNamespace(run=run, is_converged=converged)

    return factory


def _write_structure(path: Path, n_atoms: int) -> None:
    lines = [
        f"HETATM{i:5d}  C{i:<2d} LIG A   1    {2.0 * (i - 1):8.3f}   0.000   0.000"
        "  1.00 10.00           C\n"
        for i in range(1, n_atoms + 1)
    ]
    path.write_text("".join(lines) + "END\n", encoding="utf-8")


def _run_grid_scan(tmp_path, monkeypatch, command, scan_lists, *, converged=True):
    from mlmm.cli import cli as root_cli
    from mlmm.workflows import scan2d, scan3d

    module = scan2d if command == "scan2d" else scan3d
    n_atoms = 4 if command == "scan2d" else 6
    monkeypatch.setattr(module, "mlmm", lambda **kwargs: _ConstantScanCalculator())
    monkeypatch.setattr(module, "_make_lbfgs", _restraint_reaching_lbfgs(converged))
    monkeypatch.setattr(scan2d, "write_plotly_image", lambda *args, **kwargs: None)
    structure = tmp_path / "system.pdb"
    _write_structure(structure, n_atoms)
    parm = tmp_path / "system.parm7"
    parm.write_text("not reached\n", encoding="utf-8")
    out_dir = tmp_path / "out"
    result = CliRunner().invoke(
        root_cli,
        [
            command, "-i", str(structure), "--parm", str(parm),
            "--model-indices", f"1-{n_atoms}", "--no-detect-layer",
            "-q", "0", "-m", "1",
            "--scan-lists", scan_lists, "--max-step-size", "0.004",
            "--out-json", "--out-dir", str(out_dir),
        ],
    )
    return result, out_dir


def _assert_plots_skipped(result, out_dir: Path, plot_names: tuple[str, ...]) -> dict:
    assert result.exit_code == 0, result.output
    assert "[plot] NOTE:" in result.output
    assert "[plot] ERROR:" not in result.output
    assert (out_dir / "surface.csv").is_file()
    for name in plot_names:
        assert not (out_dir / name).exists()
    payload = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    assert payload["execution_status"] == "completed"
    assert "error" not in payload
    listed = [*payload["files"].values(), *payload["current_output_paths"]]
    assert "surface.csv" in listed
    assert not set(plot_names) & set(listed)
    return payload


@pytest.mark.parametrize(
    "scan_lists",
    [
        "[(1,2,1.000,1.004),(3,4,1.000,1.000)]",  # two points
        "[(1,2,1.000,1.008),(3,4,1.000,1.000)]",  # three collinear points
    ],
)
def test_scan2d_with_too_few_or_collinear_points_skips_only_plots(
    tmp_path, monkeypatch, scan_lists
) -> None:
    result, out_dir = _run_grid_scan(tmp_path, monkeypatch, "scan2d", scan_lists)

    assert "[plot] NOTE: Plots skipped" in result.output
    _assert_plots_skipped(
        result, out_dir, ("scan2d_map.png", "scan2d_landscape.html")
    )


def test_scan3d_with_coplanar_points_skips_only_the_volume_plot(
    tmp_path, monkeypatch
) -> None:
    result, out_dir = _run_grid_scan(
        tmp_path, monkeypatch, "scan3d",
        "[(1,2,1.000,1.004),(3,4,1.000,1.004),(5,6,1.000,1.000)]",
    )

    assert "[plot] NOTE: Volume plot skipped" in result.output
    payload = _assert_plots_skipped(result, out_dir, ("scan3d_density.html",))
    assert payload["n_points_usable"] == 4


def test_scan3d_csv_with_three_points_skips_the_volume_plot(tmp_path) -> None:
    from mlmm.cli import cli as root_cli

    csv_path = tmp_path / "surface.csv"
    pd.DataFrame(
        {
            "i": [0, 1, 0],
            "j": [0, 0, 1],
            "k": [0, 0, 0],
            "d1_A": [1.0, 1.5, 1.0],
            "d2_A": [2.0, 2.0, 2.5],
            "d3_A": [3.0, 3.0, 3.0],
            "energy_hartree": [-1.0, -1.001, -1.002],
        }
    ).to_csv(csv_path, index=False)
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(
        root_cli,
        ["scan3d", "--csv", str(csv_path), "--out-json", "--out-dir", str(out_dir)],
    )

    assert result.exit_code == 0, result.output
    assert "[plot] NOTE: Volume plot skipped" in result.output
    assert not (out_dir / "scan3d_density.html").exists()
    payload = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    assert payload["execution_status"] == "completed"
    assert payload["files"] == {}
    assert payload["current_output_paths"] == []


@pytest.mark.parametrize(
    ("command", "scan_lists"),
    [
        ("scan2d", "[(1,2,1.000,1.004),(3,4,1.000,1.004)]"),
        ("scan3d", "[(1,2,1.000,1.004),(3,4,1.000,1.004),(5,6,1.000,1.004)]"),
    ],
)
def test_grid_scan_without_usable_points_still_exits_1(
    tmp_path, monkeypatch, command, scan_lists
) -> None:
    result, out_dir = _run_grid_scan(
        tmp_path, monkeypatch, command, scan_lists, converged=False
    )

    assert result.exit_code == 1, result.output
    assert "[plot] No finite data for plotting." in result.output
    assert "[plot] NOTE:" not in result.output
    payload = json.loads((out_dir / "result.json").read_text())
    assert payload["execution_status"] == "completed"
    assert payload["scientific_status"] == "failed"
    assert payload["n_points_usable"] == 0
    assert payload["min_energy_hartree"] is None
    assert "status" not in payload
    assert (out_dir / "result.json").read_bytes() == (out_dir / "summary.json").read_bytes()

