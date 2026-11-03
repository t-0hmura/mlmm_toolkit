"""`opt --mm-only` names the MM backend, not the unused MLIP, in its start lines."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.core.utils import calculator_run_label

_TOY = Path(__file__).resolve().parents[1] / "examples" / "toy_system"
_R, _PARM = str(_TOY / "r_toy.pdb"), str(_TOY / "p_toy.parm7")

needs_toy_system = pytest.mark.skipif(
    not (_TOY / "p_toy.parm7").exists(), reason="toy_system example inputs not present"
)


def test_mm_only_label_names_the_mm_backend() -> None:
    assert calculator_run_label({}, mm_only=True) == "MM only (hessian_ff)"
    assert calculator_run_label({"mm_backend": "openmm"}, mm_only=True) == "MM only (openmm)"
    assert calculator_run_label({"uma_precision": "fp32"}) == "UMA (UMA-S-1.2.1 (OMol), fp32)"


@needs_toy_system
@pytest.mark.parametrize("extra", [[], ["--mm-backend", "openmm"]])
def test_opt_mm_only_dry_run_prints_the_mm_backend(tmp_path, monkeypatch, extra) -> None:
    monkeypatch.chdir(tmp_path)
    result = CliRunner().invoke(
        root_cli,
        ["opt", "-i", _R, "--parm", _PARM, "-q", "-1", "-m", "1", "--mm-only",
         "--opt-mode", "grad", "--out-dir", str(tmp_path / "out"), "--dry-run", *extra],
    )
    assert result.exit_code == 0, result.output
    backend = "openmm" if extra else "hessian_ff"
    assert f"[backend] MM only ({backend})" in result.output
    assert "[backend] UMA" not in result.output
