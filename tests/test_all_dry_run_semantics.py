from pathlib import Path
from types import SimpleNamespace

import click
import pytest
from click.testing import CliRunner

from mlmm.workflows.all import (
    _validate_all_dry_run_semantics,
    _validate_all_freeze_atom_bounds,
    cli,
)


def _pdb(path: Path) -> Path:
    path.write_text(
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "HETATM    2  C2  LIG A   1       1.500   0.000   0.000  1.00  0.00           C  \n"
        "END\n",
        encoding="utf-8",
    )
    return path


def _prepared(path: Path):
    return [SimpleNamespace(source_path=path)]


def test_dry_run_rejects_unresolved_charge_when_extraction_is_skipped(
    tmp_path: Path,
) -> None:
    source = _pdb(tmp_path / "input.pdb")

    with pytest.raises(click.ClickException, match="charge cannot be resolved"):
        _validate_all_dry_run_semantics(
            prepared_inputs=_prepared(source),
            center_spec=None,
            scan_lists_raw=(),
            scan_one_based=True,
            charge_override=None,
            yaml_model_charge=None,
            ligand_charge=None,
            spin=1,
            mm_ligand_mult=None,
        )


def test_dry_run_rejects_invalid_center_and_scan_literal(tmp_path: Path) -> None:
    source = _pdb(tmp_path / "input.pdb")
    common = dict(
        prepared_inputs=_prepared(source),
        scan_one_based=True,
        charge_override=0,
        yaml_model_charge=None,
        ligand_charge=None,
        spin=1,
        mm_ligand_mult=None,
    )

    with pytest.raises(click.BadParameter, match="Invalid -c/--center"):
        _validate_all_dry_run_semantics(
            **common,
            center_spec="MISSING",
            scan_lists_raw=(),
        )
    with pytest.raises(click.BadParameter):
        _validate_all_dry_run_semantics(
            **common,
            center_spec=None,
            scan_lists_raw=("not-a-scan",),
        )


def test_dry_run_accepts_valid_center_and_scan_syntax(tmp_path: Path) -> None:
    source = _pdb(tmp_path / "input.pdb")

    _validate_all_dry_run_semantics(
        prepared_inputs=_prepared(source),
        center_spec="LIG",
        scan_lists_raw=("[(1,2,1.4)]",),
        scan_one_based=True,
        charge_override=0,
        yaml_model_charge=None,
        ligand_charge=None,
        spin=1,
        mm_ligand_mult=None,
    )


def test_all_dry_run_does_not_print_success_for_unresolved_charge() -> None:
    root = Path(__file__).resolve().parents[1]
    source = root / "tests" / "smoke" / "p_complex_layered.pdb"
    parm = root / "tests" / "smoke" / "p_complex.parm7"
    if not source.is_file() or not parm.is_file():
        pytest.skip("smoke inputs are not present")

    result = CliRunner().invoke(
        cli,
        [
            "-i",
            str(source),
            "-i",
            str(source),
            "--parm",
            str(parm),
            "--dry-run",
        ],
    )

    assert result.exit_code != 0
    assert "charge cannot be resolved" in result.output
    assert "Dry-run validation passed" not in result.output


def test_all_freeze_atoms_are_validated_against_each_full_input(
    tmp_path: Path,
) -> None:
    source = _pdb(tmp_path / "input.pdb")
    prepared = SimpleNamespace(
        geom_path=source,
        display_path=source,
    )

    _validate_all_freeze_atom_bounds([0, 1], [prepared])
    with pytest.raises(click.BadParameter, match="atom 3.*has 2 atoms"):
        _validate_all_freeze_atom_bounds([2], [prepared])


def test_all_cli_accepts_freeze_atoms_and_reports_bounds(tmp_path: Path) -> None:
    source = _pdb(tmp_path / "input.pdb")

    result = CliRunner().invoke(
        cli,
        [
            "-i", str(source),
            "-i", str(source),
            "-q", "0",
            "--freeze-atoms", "3",
            "--dry-run",
        ],
    )

    assert result.exit_code != 0
    assert "references atom 3" in result.output
    assert "has 2 atoms" in result.output
