"""Active contracts for the optional experimental embedcharge paths."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.core.defaults import MLMM_CALC_KW


_FIXTURE_DIR = (
    Path(__file__).resolve().parents[1] / "hessian_ff" / "tests" / "data" / "small"
)
_STRUCTURE = _FIXTURE_DIR / "p_complex_layered.pdb"
_TOPOLOGY = _FIXTURE_DIR / "p_complex.parm7"


def _base_args(command: str) -> list[str]:
    args = [
        command,
        "-i",
        str(_STRUCTURE),
    ]
    if command == "all":
        args.append(str(_STRUCTURE))
    return [
        *args,
        "--parm",
        str(_TOPOLOGY),
        "-q",
        "-1",
        "-m",
        "1",
        "--dry-run",
    ]


@pytest.mark.parametrize("command", ["all", "dft"])
def test_cli_accepts_embedcharge_in_mlip_and_dft_paths(command: str) -> None:
    result = CliRunner().invoke(
        root_cli,
        [*_base_args(command), "--embedcharge", "--embedcharge-cutoff", "8.0"],
    )
    assert result.exit_code == 0, result.output


def test_yaml_embedcharge_remains_active(tmp_path: Path) -> None:
    config = tmp_path / "embedcharge.yaml"
    config.write_text(
        yaml.safe_dump(
            {"calc": {"embedcharge": True, "embedcharge_cutoff": 8.0}}
        ),
        encoding="utf-8",
    )
    result = CliRunner().invoke(
        root_cli,
        [*_base_args("all"), "--config", str(config)],
    )
    assert result.exit_code == 0, result.output


def test_embedcharge_defaults_are_opt_in() -> None:
    assert MLMM_CALC_KW["embedcharge"] is False
    assert MLMM_CALC_KW["embedcharge_cutoff"] == pytest.approx(12.0)
