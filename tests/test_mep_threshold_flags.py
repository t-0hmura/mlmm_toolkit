"""MEP threshold surfaces: GSM preset, DMF tolerance, and parent forwarding."""

from __future__ import annotations

from pathlib import Path
from typing import Set

import click
import pytest
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.workflows._all_helpers import build_path_child_argv
from mlmm.workflows.path_opt import resolve_dmf_solve_tol

MEP_FLAGS = ("--gsm-param", "--thresh-gsm", "--thresh-dmf")


def _declared_flags(cli: click.Command) -> Set[str]:
    flags: Set[str] = set()
    for param in cli.params:
        flags.update(param.opts)
    return flags


def test_parent_and_mep_children_declare_the_same_flags() -> None:
    from mlmm.workflows.all import cli as all_cli
    from mlmm.workflows.path_opt import cli as path_opt_cli
    from mlmm.workflows.path_search import cli as path_search_cli

    for flag in MEP_FLAGS:
        assert flag in _declared_flags(all_cli)
        assert flag in _declared_flags(path_opt_cli)
        assert flag in _declared_flags(path_search_cli)


def test_path_child_argv_forwards_both_mep_thresholds() -> None:
    argv = build_path_child_argv(
        {"thresh_gsm", "thresh_dmf"},
        mep_mode="gsm",
        dmf_backend="gpu",
        max_nodes=20,
        gsm_param="equi",
        max_cycles_gsm=300,
        max_cycles_dmf=300,
        climb=True,
        dump=False,
        pre_opt=True,
        convert_files=True,
        thresh=None,
        thresh_gsm="gau",
        thresh_dmf="middle",
    )
    assert argv == ["--mep-mode", "gsm", "--thresh-gsm", "gau", "--thresh-dmf", "middle"]


def test_path_child_argv_stays_silent_without_explicit_thresholds() -> None:
    argv = build_path_child_argv(
        set(),
        mep_mode="gsm",
        dmf_backend="gpu",
        max_nodes=20,
        gsm_param="equi",
        max_cycles_gsm=300,
        max_cycles_dmf=300,
        climb=True,
        dump=False,
        pre_opt=True,
        convert_files=True,
        thresh="gau",
        thresh_gsm="gau",
        thresh_dmf="middle",
    )
    assert argv == ["--mep-mode", "gsm"]


def test_default_keeps_the_historical_tight_tolerance() -> None:
    assert resolve_dmf_solve_tol({}) == "tight"


@pytest.mark.parametrize("raw", ["tight", " Middle ", "LOOSE"])
def test_presets_are_normalized(raw: str) -> None:
    assert resolve_dmf_solve_tol({"tol": raw}) == raw.strip().lower()


def test_float_tolerance_is_passed_through() -> None:
    assert resolve_dmf_solve_tol({"tol": "0.08"}) == pytest.approx(0.08)
    assert resolve_dmf_solve_tol({"tol": 0.12}) == pytest.approx(0.12)


def test_pinned_dual_inf_tol_survives_when_no_preset_is_requested() -> None:
    cfg = {"ipopt_options": {"dual_inf_tol": 0.07}}
    # None keeps pydmf's solve() from overwriting the caller's IPOPT option.
    assert resolve_dmf_solve_tol(cfg) is None


def test_explicit_tolerance_wins_over_a_pinned_ipopt_option() -> None:
    cfg = {"tol": "loose", "ipopt_options": {"dual_inf_tol": 0.07}}
    assert resolve_dmf_solve_tol(cfg) == "loose"


def test_a_gaussian_preset_is_rejected_with_a_pointer_to_the_right_flag() -> None:
    with pytest.raises(click.ClickException) as excinfo:
        resolve_dmf_solve_tol({"tol": "gau_tight"})
    message = str(excinfo.value)
    assert "tight|middle|loose" in message
    assert "--thresh-gsm" in message


@pytest.mark.parametrize("raw", ["0", "-1", "nan", "inf", "-inf", True, False])
def test_a_nonfinite_or_non_positive_tolerance_is_rejected(raw: object) -> None:
    with pytest.raises(click.ClickException, match="positive float|finite positive"):
        resolve_dmf_solve_tol({"tol": raw})


def test_all_dry_run_rejects_an_explicit_invalid_dmf_tolerance(
    tmp_path: Path,
) -> None:
    smoke = Path(__file__).resolve().parent / "smoke"
    result = CliRunner().invoke(
        root_cli,
        [
            "all", "-i", str(smoke / "r_complex_layered.pdb"),
            str(smoke / "p_complex_layered.pdb"),
            "--parm", str(smoke / "p_complex.parm7"),
            "-q", "-1", "-m", "1", "--mep-mode", "gsm",
            "--thresh-dmf", "nan", "--dry-run",
            "--out-dir", str(tmp_path / "all"),
        ],
    )
    assert result.exit_code != 0
    assert "finite positive" in result.output


def test_all_show_config_reports_each_threshold_owner(tmp_path: Path) -> None:
    smoke = Path(__file__).resolve().parent / "smoke"
    result = CliRunner().invoke(
        root_cli,
        [
            "all", "-i", str(smoke / "r_complex_layered.pdb"),
            str(smoke / "p_complex_layered.pdb"),
            "--parm", str(smoke / "p_complex.parm7"),
            "-q", "-1", "-m", "1", "--thresh", "gau",
            "--thresh-gsm", "gau_loose", "--thresh-dmf", "middle",
            "--show-config", "--dry-run", "--out-dir", str(tmp_path / "all"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "thresh: gau" in result.output
    assert "thresh_gsm: gau_loose" in result.output
    assert "thresh_dmf: middle" in result.output
    footer = result.output.rstrip().splitlines()[-1]
    assert footer == "[Dry run] --dry-run completed. Input command is valid."
    assert "Elapsed Time for Whole Pipeline" not in result.output


@pytest.mark.parametrize(
    ("mode_args", "expected_post"),
    [
        ([], "hess"),
        (["--opt-mode", "grad"], "grad"),
        (["--opt-mode", "hess", "--opt-mode-post", "grad"], "grad"),
    ],
)
def test_all_show_config_reports_effective_optimizer_modes(
    tmp_path: Path,
    mode_args: list[str],
    expected_post: str,
) -> None:
    smoke = Path(__file__).resolve().parent / "smoke"
    result = CliRunner().invoke(
        root_cli,
        [
            "all", "-i", str(smoke / "r_complex_layered.pdb"),
            str(smoke / "p_complex_layered.pdb"),
            "--parm", str(smoke / "p_complex.parm7"),
            "-q", "-1", "-m", "1", "--show-config", "--dry-run",
            "--out-dir", str(tmp_path / "all"), *mode_args,
        ],
    )

    assert result.exit_code == 0, result.output
    assert "path_opt_mode: grad" in result.output
    assert f"post_opt_mode: {expected_post}" in result.output
    assert f"opt_mode_post: {expected_post}" in result.output


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
def test_explicit_gsm_param_overrides_yaml(tmp_path: Path, command: str) -> None:
    smoke = Path(__file__).resolve().parent / "smoke"
    config = tmp_path / "gsm.yaml"
    config.write_text("gs:\n  param: equi\n", encoding="utf-8")
    result = CliRunner().invoke(
        root_cli,
        [
            command,
            "-i",
            str(smoke / "r_complex_layered.pdb"),
            str(smoke / "p_complex_layered.pdb"),
            "--parm",
            str(smoke / "p_complex.parm7"),
            "-q",
            "-1",
            "--config",
            str(config),
            "--gsm-param",
            "energy",
            "--show-config",
            "--dry-run",
            "-v",
            "3",
            "--out-dir",
            str(tmp_path / command),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "gsm_param: energy" in result.output
