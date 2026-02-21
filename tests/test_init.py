"""Tests for the `mlmm init` template generator."""

import sys

import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm_toolkit CLI requires Python >= 3.11",
)

from click.testing import CliRunner  # noqa: E402


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def cli_group():
    from mlmm_toolkit.cli import cli

    return cli


def test_init_generates_template_file(runner, cli_group, tmp_path):
    out_path = tmp_path / "mlmm_all.config.yaml"
    result = runner.invoke(cli_group, ["init", "--out", str(out_path)])
    assert result.exit_code == 0, result.output
    assert out_path.exists()
    text = out_path.read_text(encoding="utf-8")
    assert "extract:" in text
    assert "path_search:" in text
    assert "dft:" in text

