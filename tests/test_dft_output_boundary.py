"""DFT output-directory and stale-artifact boundary regressions."""

from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner


def test_prepare_dft_output_dir_invalidates_prior_public_results(
    tmp_path: Path,
) -> None:
    from mlmm.workflows.dft import _prepare_dft_output_dir

    out_dir = tmp_path / "dft"
    out_dir.mkdir()
    stale = [
        out_dir / "result.yaml",
        out_dir / "result.json",
        out_dir / "summary.json",
    ]
    for path in stale:
        path.write_text("stale\n", encoding="utf-8")

    resolved = _prepare_dft_output_dir(out_dir)

    assert resolved == out_dir.resolve()
    assert all(not path.exists() for path in stale)


def test_dft_error_json_uses_yaml_effective_output_dir(
    monkeypatch, tmp_path: Path,
) -> None:
    from mlmm.cli import cli as root_cli
    from mlmm.workflows import dft

    structure = tmp_path / "system.pdb"
    structure.write_text(
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 10.00           C\n"
        "END\n",
        encoding="utf-8",
    )
    parm = tmp_path / "system.parm7"
    parm.write_text("placeholder\n", encoding="utf-8")
    effective = tmp_path / "yaml-dft"
    config = tmp_path / "config.yaml"
    config.write_text(
        f"dft:\n  out_dir: {effective}\n",
        encoding="utf-8",
    )

    def fail_charge(*_args, **_kwargs):
        raise RuntimeError("probe failure")

    monkeypatch.setattr(dft, "resolve_charge_spin_or_raise", fail_charge)
    result = CliRunner().invoke(
        root_cli,
        [
            "dft", "-i", str(structure), "--parm", str(parm),
            "--model-pdb", str(structure), "-q", "0", "-m", "1",
            "--config", str(config),
        ],
    )

    assert result.exit_code != 0
    payload = json.loads((effective / "result.json").read_text())
    assert payload["status"] == "error"
    assert payload["error"] == "probe failure"
    assert not (tmp_path / "result_dft" / "result.json").exists()
