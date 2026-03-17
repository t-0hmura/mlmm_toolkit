"""CLI regressions for recently added/removed options."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest
from click.testing import CliRunner

from mlmm.cli import cli as root_cli

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm CLI requires Python >= 3.11",
)


def test_freeze_links_removed_from_help_outputs() -> None:
    runner = CliRunner()
    for command_name in ("opt", "tsopt", "freq", "irc"):
        short = runner.invoke(root_cli, [command_name, "--help"])
        assert short.exit_code == 0, short.output
        assert "--freeze-links" not in short.output

        advanced = runner.invoke(root_cli, [command_name, "--help-advanced"])
        assert advanced.exit_code == 0, advanced.output
        assert "--freeze-links" not in advanced.output


def test_freeze_links_now_fails_as_unknown_option() -> None:
    runner = CliRunner()
    for command_name in ("opt", "tsopt", "freq", "irc"):
        result = runner.invoke(root_cli, [command_name, "--freeze-links"])
        assert result.exit_code != 0
        assert "No such option" in result.output


def test_path_search_help_shows_refine_mode() -> None:
    runner = CliRunner()
    result = runner.invoke(root_cli, ["path-search", "--help"])
    assert result.exit_code == 0, result.output
    assert "--refine-mode" in result.output


def test_path_opt_help_shows_fix_ends() -> None:
    runner = CliRunner()
    result = runner.invoke(root_cli, ["path-opt", "--help"])
    assert result.exit_code == 0, result.output
    assert "--fix-ends" in result.output


def test_scan3d_csv_mode_runs_without_scan_inputs(tmp_path: Path) -> None:
    csv_path = tmp_path / "surface.csv"
    out_dir = tmp_path / "scan3d_out"

    rows = []
    for i in (0.0, 1.0):
        for j in (0.0, 1.0):
            for k in (0.0, 1.0):
                rows.append(
                    {
                        "d1_A": i,
                        "d2_A": j,
                        "d3_A": k,
                        "energy_kcal": i + j + k,
                        "d1_label": "d1",
                        "d2_label": "d2",
                        "d3_label": "d3",
                    }
                )
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        ["scan3d", "--csv", str(csv_path), "--out-dir", str(out_dir)],
    )
    assert result.exit_code == 0, result.output
    assert (out_dir / "scan3d_density.html").exists()


@pytest.mark.parametrize("extra_args", [["-q", "0"], ["-m", "2"]])
def test_trj2fig_recompute_branch_triggers_on_charge_or_multiplicity(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    extra_args: list[str],
) -> None:
    from mlmm import trj2fig as trj2fig_mod

    xyz_path = tmp_path / "traj.xyz"
    xyz_path.write_text(
        "1\ncomment\nH 0.0 0.0 0.0\n1\ncomment\nH 0.0 0.0 0.1\n",
        encoding="utf-8",
    )
    out_csv = tmp_path / "energy.csv"
    called = {"value": False}

    def _fake_recompute(_traj: Path, _charge: int | None, _mult: int | None):
        called["value"] = True
        return [0.0, 0.001]

    monkeypatch.setattr(trj2fig_mod, "recompute_energies", _fake_recompute)

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        ["trj2fig", "-i", str(xyz_path), "-o", str(out_csv), *extra_args],
    )
    assert result.exit_code == 0, result.output
    assert called["value"] is True
    assert out_csv.exists()
