"""Safety contracts for ``mlmm add-elem-info`` output selection."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from mlmm.domain.add_elem_info import cli


PDB_TEXT = (
    "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00              \n"
    "TER\nEND\n"
)


def test_default_is_non_destructive_and_inplace_is_explicit(tmp_path: Path) -> None:
    source = tmp_path / "enzyme.pdb"
    source.write_text(PDB_TEXT, encoding="utf-8")
    before = source.read_bytes()

    default = CliRunner().invoke(cli, ["-i", str(source)])
    assert default.exit_code == 0, default.output
    assert source.read_bytes() == before
    assert (tmp_path / "enzyme_add_elem.pdb").is_file()

    # --overwrite controls existing element fields; it does not imply file
    # replacement. Destructive output requires the separate --inplace switch.
    field_overwrite = CliRunner().invoke(cli, ["-i", str(source), "--overwrite"])
    assert field_overwrite.exit_code == 0, field_overwrite.output
    assert source.read_bytes() == before

    inplace = CliRunner().invoke(cli, ["-i", str(source), "--inplace"])
    assert inplace.exit_code == 0, inplace.output
    assert "Wrote: " + str(source) in inplace.output


def test_explicit_output_wins_over_inplace(tmp_path: Path) -> None:
    source = tmp_path / "enzyme.pdb"
    target = tmp_path / "fixed.pdb"
    source.write_text(PDB_TEXT, encoding="utf-8")
    before = source.read_bytes()

    result = CliRunner().invoke(
        cli, ["-i", str(source), "-o", str(target), "--inplace"]
    )
    assert result.exit_code == 0, result.output
    assert target.is_file()
    assert source.read_bytes() == before
