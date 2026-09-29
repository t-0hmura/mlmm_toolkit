"""YAML routing of single-structure optimizer settings in path-opt and path-search."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from click.testing import CliRunner

from mlmm.cli import cli as root_cli

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm CLI requires Python >= 3.11",
)

SMOKE = Path(__file__).resolve().parent / "smoke"


def _invoke(tmp_path: Path, command: str, yaml_text: str, *extra: str):
    config = tmp_path / "config.yaml"
    config.write_text(yaml_text, encoding="utf-8")
    return CliRunner().invoke(
        root_cli,
        [
            command,
            "-i", str(SMOKE / "r_complex_layered.pdb"), str(SMOKE / "p_complex_layered.pdb"),
            "--parm", str(SMOKE / "p_complex.parm7"),
            "-q", "-1", "-m", "1",
            "--config", str(config), "--dry-run",
            "--out-dir", str(tmp_path / command), *extra,
        ],
    )


def _block(output: str, title: str) -> str:
    head = f"\n{title}\n{'-' * len(title)}\n"
    start = output.index(head) + len(head)
    end = output.find("\n\n", start)
    return output[start:] if end < 0 else output[start:end]


@pytest.mark.parametrize(
    "yaml_text", ["opt:\n  max_cycles: 7\n", "stopt:\n  lbfgs:\n    max_cycles: 7\n"]
)
def test_opt_section_configures_preoptimization(tmp_path: Path, yaml_text: str) -> None:
    result = _invoke(tmp_path, "path-opt", yaml_text, "-v", "3")

    assert result.exit_code == 0, result.output
    assert "preopt_max_cycles: 7" in _block(result.output, "dry_run_plan")


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
@pytest.mark.parametrize(
    ("yaml_text", "message"),
    [
        (
            "opt:\n  thresh: gau\nlbfgs:\n  thresh: baker\n",
            "opt.thresh and lbfgs.thresh conflict",
        ),
        (
            "opt:\n  thresh: gau\n  lbfgs:\n    thresh: baker\n",
            "opt.thresh and opt.lbfgs.thresh conflict",
        ),
    ],
)
def test_conflicting_opt_and_lbfgs_values_are_rejected(
    tmp_path: Path, command: str, yaml_text: str, message: str
) -> None:
    result = _invoke(tmp_path, command, yaml_text)

    assert result.exit_code == 2, result.output
    assert message in result.output


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
def test_conflicting_lbfgs_and_stopt_lbfgs_values_are_rejected(
    tmp_path: Path, command: str
) -> None:
    result = _invoke(
        tmp_path, command,
        "lbfgs:\n  max_cycles: 5\nstopt:\n  lbfgs:\n    max_cycles: 7\n",
    )

    assert result.exit_code == 2, result.output
    assert "lbfgs.max_cycles and stopt.lbfgs.max_cycles conflict" in result.output


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
def test_non_mapping_opt_section_is_rejected(tmp_path: Path, command: str) -> None:
    result = _invoke(tmp_path, command, "opt: 5\n")

    assert result.exit_code == 2, result.output
    assert "YAML section 'opt' must be a mapping" in result.output


def test_apply_single_opt_yaml_layer_keeps_optimizer_keys_out_of_stopt() -> None:
    from mlmm.core.utils import deep_update
    from mlmm.workflows._path_yaml_helpers import apply_single_opt_yaml_layer

    layer = {
        "opt": {"thresh": "gau_tight"},
        "stopt": {"max_cycles": 40, "lbfgs": {"max_cycles": 7}},
    }
    stopt_cfg = {"max_cycles": 40, "lbfgs": {"max_cycles": 7}}
    lbfgs_cfg = {"thresh": "gau", "max_cycles": 100}
    rfo_cfg = {"thresh": "gau", "max_cycles": 100}

    apply_single_opt_yaml_layer(
        layer,
        lbfgs_cfg=lbfgs_cfg,
        rfo_cfg=rfo_cfg,
        stopt_cfg=stopt_cfg,
        opt_base_kw={"thresh": "gau", "max_cycles": 100},
        deep_update=deep_update,
    )

    assert stopt_cfg == {"max_cycles": 40}
    assert lbfgs_cfg == {"thresh": "gau_tight", "max_cycles": 7}
    assert rfo_cfg == {"thresh": "gau_tight", "max_cycles": 100}
