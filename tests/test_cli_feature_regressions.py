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
    from mlmm.io import trj2fig as trj2fig_mod

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


def test_coord_type_dlc_falls_back_to_cart_under_lbfgs() -> None:
    """`--coord-type dlc` is only meaningful with Hessian-based optimization;
    under `--opt-mode grad` (L-BFGS) it must fall back to Cartesian, while
    `--opt-mode hess` keeps DLC."""
    repo = Path(__file__).resolve().parents[1]
    in_pdb = repo / "examples" / "toy_system" / "p_complex_layered.pdb"
    parm = repo / "examples" / "toy_system" / "p_complex.parm7"
    if not (in_pdb.exists() and parm.exists()):
        pytest.skip("toy_system example inputs not present")
    runner = CliRunner()
    base = [
        "opt", "-i", str(in_pdb), "--parm", str(parm), "-q", "0",
        "--detect-layer", "--coord-type", "dlc", "--dry-run",
    ]
    grad = runner.invoke(root_cli, base + ["--opt-mode", "grad"])
    assert grad.exit_code == 0, grad.output
    assert "falling back to cart" in grad.output

    hess = runner.invoke(root_cli, base + ["--opt-mode", "hess"])
    assert hess.exit_code == 0, hess.output
    assert "falling back to cart" not in hess.output


def test_verbose_is_a_per_subcommand_option() -> None:
    """`-v/--verbose LEVEL` is injected into every subcommand (and carried by the
    parser-wrapper `extract`); it is no longer a root-group option, so a
    root-placed `-v` is rejected."""
    runner = CliRunner()
    for name in (
        "opt", "tsopt", "freq", "irc", "sp", "scan",
        "path-search", "dft", "all", "extract",
    ):
        res = runner.invoke(root_cli, [name, "--help"])
        assert res.exit_code == 0, res.output
        assert "-v, --verbose" in res.output, f"{name} --help is missing -v"

    # Root-placed `-v` no longer exists (it moved onto the subcommands).
    root = runner.invoke(root_cli, ["-v", "2", "opt", "--help"])
    assert root.exit_code != 0
    assert "No such option" in root.output

    # The level is an IntRange(0, 3); 0/1/2/3 are accepted (default 2) and
    # out-of-range values are rejected, for both the injected commands and the
    # parser-wrapper `extract`.
    for cmd in (["opt", "-v", "4"], ["extract", "-v", "4"]):
        bad = runner.invoke(root_cli, cmd)
        assert bad.exit_code != 0, cmd
        assert "is not in the range" in bad.output, cmd
