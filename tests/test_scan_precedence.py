"""Behavioral regressions for scan2d/scan3d effective optimizer settings."""

import json
from pathlib import Path

import click
import pytest
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.cli.decorators import make_is_param_explicit
from mlmm.workflows import scan as scan_workflow
from mlmm.workflows import scan2d as scan2d_workflow
from mlmm.workflows.scan_common import (
    build_scan_lbfgs_kwargs,
    resolve_scan_optimizer_configs,
)


@pytest.mark.parametrize(
    ("yaml_cfg", "explicit", "thresh", "cycles", "expected_thresh", "expected_cycles"),
    [
        ({}, set(), "baker", None, "baker", 100000),
        (
            {"opt": {"thresh": "gau_loose", "max_cycles": 77}},
            set(),
            "baker",
            10000,
            "gau_loose",
            77,
        ),
        (
            {"opt": {"thresh": "gau_loose", "max_cycles": 77}},
            {"thresh", "relax_max_cycles"},
            "baker",
            88,
            "baker",
            88,
        ),
    ],
)
def test_scan_optimizer_precedence_matrix(
    yaml_cfg,
    explicit,
    thresh,
    cycles,
    expected_thresh,
    expected_cycles,
) -> None:
    opt_cfg, lbfgs_cfg = resolve_scan_optimizer_configs(
        yaml_cfg,
        thresh=thresh,
        relax_max_cycles=cycles,
        is_param_explicit=lambda name: name in explicit,
    )

    assert opt_cfg["thresh"] == expected_thresh
    assert opt_cfg["max_cycles"] == expected_cycles
    kwargs = build_scan_lbfgs_kwargs(
        lbfgs_cfg,
        opt_cfg,
        max_step_bohr=0.2,
        out_dir=Path("result"),
        prefix="point",
    )
    assert kwargs["thresh"] == expected_thresh
    assert kwargs["max_cycles"] == expected_cycles


def test_scan_optimizer_resolution_does_not_mutate_yaml() -> None:
    yaml_cfg = {
        "opt": {"thresh": "gau_loose", "max_cycles": 17},
        "lbfgs": {"keep_last": 3},
    }
    expected = {
        "opt": {"thresh": "gau_loose", "max_cycles": 17},
        "lbfgs": {"keep_last": 3},
    }

    _, first_lbfgs = resolve_scan_optimizer_configs(
        yaml_cfg,
        thresh="baker",
        relax_max_cycles=10000,
        is_param_explicit=lambda _name: False,
    )
    first_lbfgs["keep_last"] = 99
    _, second_lbfgs = resolve_scan_optimizer_configs(
        yaml_cfg,
        thresh="baker",
        relax_max_cycles=10000,
        is_param_explicit=lambda _name: False,
    )

    assert yaml_cfg == expected
    assert second_lbfgs["keep_last"] == 3


def test_nested_opt_lbfgs_does_not_reach_scan_constructor_kwargs() -> None:
    opt_cfg, lbfgs_cfg = resolve_scan_optimizer_configs(
        {"opt": {"max_cycles": 17, "lbfgs": {"max_step": 0.12}}},
        thresh="baker",
        relax_max_cycles=10000,
        is_param_explicit=lambda _name: False,
    )

    kwargs = build_scan_lbfgs_kwargs(
        lbfgs_cfg,
        opt_cfg,
        max_step_bohr=0.2,
        out_dir=Path("result"),
        prefix="point",
    )

    assert "lbfgs" not in opt_cfg
    assert "lbfgs" not in kwargs
    assert kwargs["max_cycles"] == 17
    assert kwargs["max_step"] == 0.12


@pytest.mark.parametrize(
    ("yaml_cfg", "argv", "expected"),
    [
        ({}, [], {"thresh": "baker", "max_cycles": 100000}),
        (
            {"opt": {"thresh": "gau_loose", "max_cycles": 17}},
            [],
            {"thresh": "gau_loose", "max_cycles": 17},
        ),
        (
            {"opt": {"thresh": "gau_loose", "max_cycles": 17}},
            ["--thresh", "baker", "--relax-max-cycles", "23"],
            {"thresh": "baker", "max_cycles": 23},
        ),
    ],
)
def test_click_parameter_sources_drive_scan_precedence(
    yaml_cfg,
    argv,
    expected,
) -> None:
    @click.command()
    @click.option("--thresh", default="baker")
    @click.option("--relax-max-cycles", type=int, default=None)
    @click.pass_context
    def command(ctx, thresh, relax_max_cycles):
        opt_cfg, _ = resolve_scan_optimizer_configs(
            yaml_cfg,
            thresh=thresh,
            relax_max_cycles=relax_max_cycles,
            is_param_explicit=make_is_param_explicit(ctx),
        )
        click.echo(
            json.dumps(
                {
                    "thresh": opt_cfg["thresh"],
                    "max_cycles": opt_cfg["max_cycles"],
                }
            )
        )

    result = CliRunner().invoke(command, argv)
    assert result.exit_code == 0, result.output
    assert json.loads(result.output) == expected


def _scan_kwargs(opt_cfg, lbfgs_cfg):
    return build_scan_lbfgs_kwargs(
        lbfgs_cfg,
        opt_cfg,
        max_step_bohr=0.2,
        out_dir=Path("result"),
        prefix="point",
    )


def test_scan_optimizer_section_applies_when_opt_keeps_defaults() -> None:
    opt_cfg, lbfgs_cfg = resolve_scan_optimizer_configs(
        {"lbfgs": {"max_cycles": 50, "thresh": "gau_tight", "print_every": 7}},
        thresh="baker",
        relax_max_cycles=None,
        is_param_explicit=lambda _name: False,
    )

    kwargs = _scan_kwargs(opt_cfg, lbfgs_cfg)
    assert (kwargs["max_cycles"], kwargs["thresh"], kwargs["print_every"]) == (
        50,
        "gau_tight",
        7,
    )


@pytest.mark.parametrize(
    ("yaml_cfg", "explicit"),
    [
        ({"opt": {"max_cycles": 10}, "lbfgs": {"max_cycles": 50}}, set()),
        ({"lbfgs": {"max_cycles": 50}}, {"relax_max_cycles"}),
    ],
)
def test_scan_explicit_shared_key_conflict_is_an_error(yaml_cfg, explicit) -> None:
    with pytest.raises(
        click.BadParameter, match="opt.max_cycles and lbfgs.max_cycles conflict"
    ):
        resolve_scan_optimizer_configs(
            yaml_cfg,
            thresh="baker",
            relax_max_cycles=10,
            is_param_explicit=lambda name: name in explicit,
        )


def test_grid_scan_default_cycle_budget_matches_help() -> None:
    opt_cfg, lbfgs_cfg = resolve_scan_optimizer_configs(
        {},
        opt_defaults=scan2d_workflow.OPT_BASE_KW,
        lbfgs_defaults=scan2d_workflow.LBFGS_KW,
        thresh="baker",
        relax_max_cycles=None,
        is_param_explicit=lambda _name: False,
    )

    assert _scan_kwargs(opt_cfg, lbfgs_cfg)["max_cycles"] == 100000


def _scan_dry_run(tmp_path, monkeypatch, yaml_text, *extra):
    smoke = Path(__file__).resolve().parent / "smoke"
    config = tmp_path / "scan.yaml"
    config.write_text(yaml_text, encoding="utf-8")
    blocks = {}

    def capture(title, content, **_kwargs):
        blocks[title] = dict(content)
        return ""

    monkeypatch.setattr(scan_workflow, "pretty_block", capture)
    result = CliRunner().invoke(
        root_cli,
        [
            "scan", "-i", str(smoke / "r_complex_layered.pdb"),
            "--parm", str(smoke / "p_complex.parm7"),
            "--model-pdb", str(smoke / "pocket_r.pdb"),
            "-q", "-1", "-m", "1", "-s", "[(1,2,1.8)]",
            "--config", str(config), "--dry-run",
            "--out-dir", str(tmp_path / "out"), *extra,
        ],
    )
    return result, blocks


def test_1d_scan_rejects_conflicting_opt_and_lbfgs_cycles(tmp_path, monkeypatch):
    result, _ = _scan_dry_run(
        tmp_path, monkeypatch, "opt:\n  max_cycles: 10\nlbfgs:\n  max_cycles: 50\n"
    )

    assert result.exit_code == 2, result.output
    assert "opt.max_cycles and lbfgs.max_cycles conflict" in result.output


def test_1d_scan_keeps_yaml_opt_dump_off(tmp_path, monkeypatch):
    result, blocks = _scan_dry_run(tmp_path, monkeypatch, "opt:\n  dump: true\n")

    assert result.exit_code == 0, result.output
    assert blocks["opt"].get("dump", False) is False


def test_1d_scan_max_cycles_alias_is_hidden_and_checked(tmp_path, monkeypatch):
    help_result = CliRunner().invoke(root_cli, ["scan", "--help-advanced"])
    assert help_result.exit_code == 0, help_result.output
    assert "--max-cycles" not in help_result.output.replace("--relax-max-cycles", "")

    same, _ = _scan_dry_run(
        tmp_path, monkeypatch, "{}\n", "--max-cycles", "7", "--relax-max-cycles", "7",
    )
    assert same.exit_code == 0, same.output

    conflict, _ = _scan_dry_run(
        tmp_path, monkeypatch, "{}\n", "--max-cycles", "5", "--relax-max-cycles", "7",
    )
    assert conflict.exit_code == 2, conflict.output
    assert "--max-cycles and --relax-max-cycles conflict" in conflict.output
