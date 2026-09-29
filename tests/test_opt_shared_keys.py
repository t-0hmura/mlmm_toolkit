"""Precedence between the shared ``opt`` block and the lbfgs/rfo sections."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from click.testing import CliRunner

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm requires Python >= 3.11",
)


def _dry_run(tmp_path: Path, monkeypatch, yaml_text: str, *extra: str):
    opt_mod = pytest.importorskip("mlmm.workflows.opt")
    smoke = Path(__file__).resolve().parent / "smoke"
    source = smoke / "p_complex_layered.pdb"
    parm = smoke / "p_complex.parm7"
    if not source.is_file() or not parm.is_file():
        pytest.skip("smoke inputs are not present")
    config = tmp_path / "config.yaml"
    config.write_text(yaml_text, encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    return CliRunner().invoke(
        opt_mod.cli,
        [
            "-i", str(source), "--parm", str(parm), "-q", "-1", "-m", "1",
            "--config", str(config), "--dry-run", "-o", str(tmp_path / "out"),
            *extra,
        ],
    )


def test_optimizer_section_value_applies_when_opt_keeps_defaults(tmp_path, monkeypatch):
    result = _dry_run(
        tmp_path, monkeypatch, "lbfgs:\n  max_cycles: 50\n", "--opt-mode", "grad",
    )

    assert result.exit_code == 0, result.output
    assert "max_cycles=50" in result.output


def test_unused_optimizer_section_does_not_leak(tmp_path, monkeypatch):
    result = _dry_run(
        tmp_path, monkeypatch, "rfo:\n  max_cycles: 50\n", "--opt-mode", "grad",
    )

    assert result.exit_code == 0, result.output
    assert "max_cycles=100000" in result.output


@pytest.mark.parametrize(
    ("yaml_text", "extra"),
    [
        ("opt:\n  max_cycles: 10\nrfo:\n  max_cycles: 50\n", ()),
        ("rfo:\n  max_cycles: 50\n", ("--max-cycles", "10")),
    ],
)
def test_explicit_conflict_is_an_error(tmp_path, monkeypatch, yaml_text, extra):
    result = _dry_run(tmp_path, monkeypatch, yaml_text, "--opt-mode", "hess", *extra)

    assert result.exit_code == 2, result.output
    assert "opt.max_cycles and rfo.max_cycles conflict" in result.output


def test_yaml_opt_mode_selects_optimizer_unless_cli_is_explicit(tmp_path, monkeypatch):
    yaml_only = _dry_run(tmp_path, monkeypatch, "opt:\n  opt_mode: hess\n")
    explicit = _dry_run(
        tmp_path, monkeypatch, "opt:\n  opt_mode: hess\n", "--opt-mode", "grad",
    )

    assert yaml_only.exit_code == 0, yaml_only.output
    assert "RFO (hess)" in yaml_only.output
    assert explicit.exit_code == 0, explicit.output
    assert "LBFGS (grad)" in explicit.output
