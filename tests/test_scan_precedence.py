"""Behavioral regressions for scan2d/scan3d effective optimizer settings."""

import json
from pathlib import Path

import click
import pytest
from click.testing import CliRunner

from mlmm.cli.decorators import make_is_param_explicit
from mlmm.workflows.scan_common import (
    build_scan_lbfgs_kwargs,
    resolve_scan_optimizer_configs,
)


@pytest.mark.parametrize(
    ("yaml_cfg", "explicit", "thresh", "cycles", "expected_thresh", "expected_cycles"),
    [
        ({}, set(), "baker", 10000, "baker", 10000),
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


@pytest.mark.parametrize(
    ("yaml_cfg", "argv", "expected"),
    [
        ({}, [], {"thresh": "baker", "max_cycles": 10000}),
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
    @click.option("--relax-max-cycles", type=int, default=10000)
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
