"""DFT output-directory and stale-artifact boundary regressions."""

from __future__ import annotations

import json
from pathlib import Path

import click
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
        out_dir / "ml_region_without_linkH.xyz",
        out_dir / "ml_region_with_linkH.xyz",
        out_dir / "ml_region_without_linkH.pdb",
        out_dir / "ml_region_with_linkH.pdb",
    ]
    for path in stale:
        path.write_text("stale\n", encoding="utf-8")
    unrelated = out_dir / "notes.txt"
    unrelated.write_text("keep\n", encoding="utf-8")

    resolved = _prepare_dft_output_dir(out_dir)

    assert resolved == out_dir.resolve()
    assert all(not path.exists() for path in stale)
    assert unrelated.read_text(encoding="utf-8") == "keep\n"


def test_dft_multiplicity_default_does_not_mask_yaml() -> None:
    from mlmm.workflows.dft import cli

    parameter = next(param for param in cli.params if param.name == "spin")
    assert isinstance(parameter, click.Option)
    assert parameter.default is None


def test_dft_rejects_unknown_yaml_engine_before_layer_setup(
    tmp_path: Path,
) -> None:
    from mlmm.cli import cli as root_cli

    structure = tmp_path / "system.pdb"
    structure.write_text(
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 10.00           C\n"
        "END\n",
        encoding="utf-8",
    )
    parm = tmp_path / "system.parm7"
    parm.write_text("not needed before engine validation\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text("dft:\n  engine: quantum-potato\n", encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            "dft",
            "-i",
            str(structure),
            "--parm",
            str(parm),
            "--model-pdb",
            str(structure),
            "-q",
            "0",
            "--config",
            str(config),
            "--dry-run",
        ],
    )

    assert result.exit_code != 0
    assert "dft.engine must be either 'cpu' or 'gpu'" in result.output


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
