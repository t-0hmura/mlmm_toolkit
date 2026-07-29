"""Safety contracts for ``mlmm add-elem-info`` output selection."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from mlmm.domain.add_elem_info import cli


PDB_TEXT = (
    "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00              \n"
    "TER\nEND\n"
)


def test_element_inference_disambiguates_inosine_and_numeric_water_hydrogen() -> None:
    from mlmm.domain.add_elem_info import guess_element

    assert guess_element("C1'", "I", False) == "C"
    assert guess_element("N9", "I", False) == "N"
    assert guess_element("1HW", "HOH", True) == "H"
    assert guess_element("OW", "HOH", True) == "O"
    assert guess_element(" NA ", "LIG", True) == "N"
    assert guess_element("PT  ", "LIG", True) == "Pt"
    assert guess_element(" PT ", "LIG", True) == "P"


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


def test_only_element_columns_change_and_other_records_are_preserved(
    tmp_path: Path,
) -> None:
    atom = (
        f"{'HETATM':<6}{7:>5} {'PT  ':4} {'LIG':>3} {'A':1}{1:>4}    "
        f"{0.0:8.3f}{1.0:8.3f}{2.0:8.3f}{1.0:6.2f}{10.0:6.2f}"
        f"{'':10}{'':>2}{'2+':>2}\r\n"
    )
    records = [
        "HEADER    ELEMENT FIELD TEST\r\n",
        "REMARK   1 KEEP THIS TEXT\r\n",
        "CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1\r\n",
        atom,
        "LINK         PT  LIG A   1                 O   HOH A   2\r\n",
        "CONECT    7    8\r\n",
        "END\r\n",
    ]
    source = tmp_path / "records.pdb"
    target = tmp_path / "fixed.pdb"
    source.write_bytes("".join(records).encode("ascii"))

    result = CliRunner().invoke(cli, ["-i", str(source), "-o", str(target)])

    assert result.exit_code == 0, result.output
    actual = target.read_bytes().decode("ascii").splitlines(keepends=True)
    assert actual[:3] == records[:3]
    assert actual[4:] == records[4:]
    assert actual[3][:76] == atom[:76]
    assert actual[3][76:78] == "Pt"
    assert actual[3][78:] == atom[78:]
