"""Path-optimization output-boundary regressions."""

from __future__ import annotations

from pathlib import Path

import click
import pytest
from click.testing import CliRunner


def test_prepare_path_output_dir_invalidates_prior_envelopes(
    tmp_path: Path,
) -> None:
    from mlmm.workflows.path_opt import _prepare_path_output_dir

    out_dir = tmp_path / "path"
    out_dir.mkdir()
    for name in ("result.json", "summary.json"):
        (out_dir / name).write_text('{"status": "error"}\n', encoding="utf-8")
    unrelated = out_dir / "notes.txt"
    unrelated.write_text("keep\n", encoding="utf-8")

    resolved = _prepare_path_output_dir(out_dir)

    assert resolved == out_dir.resolve()
    assert not (out_dir / "result.json").exists()
    assert not (out_dir / "summary.json").exists()
    assert unrelated.read_text(encoding="utf-8") == "keep\n"


@pytest.mark.parametrize("name", ["result.json", "model_from_bfactor.pdb"])
def test_path_output_collision_preserves_reserved_input(
    tmp_path: Path,
    name: str,
) -> None:
    from mlmm.workflows.path_opt import _reject_path_output_collisions

    out_dir = tmp_path / "path"
    out_dir.mkdir()
    source = out_dir / name
    source.write_text("input\n", encoding="utf-8")

    with pytest.raises(click.UsageError, match="collides"):
        _reject_path_output_collisions(out_dir, (source,))

    assert source.read_text(encoding="utf-8") == "input\n"


def test_path_cli_preserves_config_before_early_validation(
    tmp_path: Path,
) -> None:
    from mlmm.workflows import path_opt as path_module

    repo = Path(__file__).resolve().parents[1]
    smoke = repo / "tests" / "smoke"
    out_dir = tmp_path / "path"
    out_dir.mkdir()
    config = out_dir / "result.json"
    original = b"stopt:\n  max_cycles: 0\n"
    config.write_bytes(original)

    result = CliRunner().invoke(
        path_module.cli,
        [
            "-i",
            str(smoke / "r_complex_layered.pdb"),
            str(smoke / "p_complex_layered.pdb"),
            "--parm",
            str(smoke / "p_complex.parm7"),
            "-q",
            "-1",
            "-m",
            "1",
            "--config",
            str(config),
            "--out-dir",
            str(out_dir),
        ],
    )

    assert result.exit_code == 2, result.output
    assert "collides with a reserved path-opt output" in result.output
    assert config.read_bytes() == original


def test_endpoint_identity_rejects_reordered_same_element_atoms(
    tmp_path: Path,
) -> None:
    from mlmm.core.utils import (
        prepare_input_structure,
        validate_endpoint_atom_identities,
    )

    reactant = tmp_path / "reactant.pdb"
    product = tmp_path / "product.pdb"
    first = (
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
    )
    second = (
        "HETATM    2  C2  LIG A   1       1.000   0.000   0.000  1.00  0.00           C\n"
    )
    reactant.write_text(first + second + "END\n", encoding="utf-8")
    product.write_text(second + first + "END\n", encoding="utf-8")
    prepared = [
        prepare_input_structure(reactant),
        prepare_input_structure(product),
    ]
    try:
        with pytest.raises(
            click.BadParameter,
            match="ordered atom identity mismatch at atom 1",
        ):
            validate_endpoint_atom_identities(prepared)
    finally:
        for item in prepared:
            item.cleanup()
