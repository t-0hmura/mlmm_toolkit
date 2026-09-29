"""Safety contracts for ``mlmm add-elem-info`` output selection."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner
import pytest

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
    assert guess_element(" EP ", "HOH", True) == "EP"
    assert guess_element("MW", "SOL", True) == "EP"
    assert guess_element(" SE ", "SEC", False) == "Se"
    assert guess_element("NA", "LIG", True) == "Na"
    assert guess_element(" NA ", "LIG", True) == "N"
    assert guess_element("PT  ", "LIG", True) == "Pt"
    assert guess_element(" PT ", "LIG", True) == "P"


def test_element_inference_reads_leap_halogens_and_four_character_hydrogens() -> None:
    from mlmm.domain.add_elem_info import guess_element

    assert guess_element(" CL1", "LIG", False) == "Cl"
    assert guess_element(" BR1", "LIG", False) == "Br"
    assert guess_element(" C1 ", "LIG", False) == "C"
    assert guess_element("HG11", "LIG", True) == "H"
    assert guess_element("HO2A", "LIG", True) == "H"
    assert guess_element("HG  ", "LIG", True) == "Hg"
    assert guess_element("HG1 ", "LIG", True) == "Hg"


def test_default_is_non_destructive_and_overwrite_replaces_input(tmp_path: Path) -> None:
    source = tmp_path / "enzyme.pdb"
    source.write_text(PDB_TEXT, encoding="utf-8")
    before = source.read_bytes()

    default = CliRunner().invoke(cli, ["-i", str(source)])
    assert default.exit_code == 0, default.output
    assert source.read_bytes() == before
    assert (tmp_path / "enzyme_add_elem.pdb").is_file()

    overwrite = CliRunner().invoke(cli, ["-i", str(source), "--overwrite"])
    assert overwrite.exit_code == 0, overwrite.output
    assert "Wrote: " + str(source) in overwrite.output
    assert source.read_bytes()[76:78] == b" C"


def test_inplace_remains_a_hidden_alias_of_overwrite(tmp_path: Path) -> None:
    source = tmp_path / "enzyme.pdb"
    source.write_text(PDB_TEXT, encoding="utf-8")

    inplace = CliRunner().invoke(cli, ["-i", str(source), "--inplace"])
    assert inplace.exit_code == 0, inplace.output
    assert "Wrote: " + str(source) in inplace.output

    conflict = CliRunner().invoke(
        cli, ["-i", str(source), "--overwrite", "--no-inplace"]
    )
    assert conflict.exit_code == 2
    assert "Conflicting values" in conflict.output
    assert "--inplace" not in CliRunner().invoke(cli, ["--help"]).output


@pytest.mark.parametrize("flag", ["--overwrite", "--inplace"])
def test_explicit_output_wins_over_overwrite(tmp_path: Path, flag: str) -> None:
    source = tmp_path / "enzyme.pdb"
    target = tmp_path / "fixed.pdb"
    source.write_text(PDB_TEXT, encoding="utf-8")
    before = source.read_bytes()

    result = CliRunner().invoke(
        cli, ["-i", str(source), "-o", str(target), flag]
    )
    assert result.exit_code == 0, result.output
    assert target.is_file()
    assert source.read_bytes() == before


def test_existing_element_fields_are_reinferred(tmp_path: Path) -> None:
    source = tmp_path / "enzyme.pdb"
    target = tmp_path / "fixed.pdb"
    source.write_text(PDB_TEXT.replace("              \n", "          CA  \n", 1), encoding="utf-8")

    result = CliRunner().invoke(cli, ["-i", str(source), "-o", str(target)])

    assert result.exit_code == 0, result.output
    assert target.read_text(encoding="utf-8")[76:78] == " C"
    assert "assigned/updated            : 1" in result.output


def test_output_that_aliases_the_input_requires_overwrite(tmp_path: Path) -> None:
    source = tmp_path / "enzyme.pdb"
    source.write_text(PDB_TEXT, encoding="utf-8")
    before = source.read_bytes()
    link = tmp_path / "link.pdb"
    link.symlink_to(source)

    for output in (source, link):
        refused = CliRunner().invoke(cli, ["-i", str(source), "-o", str(output)])
        assert refused.exit_code == 2
        assert "use --overwrite" in refused.output
        assert source.read_bytes() == before

    accepted = CliRunner().invoke(
        cli, ["-i", str(source), "-o", str(source), "--overwrite"]
    )
    assert accepted.exit_code == 0, accepted.output
    assert source.read_bytes()[76:78] == b" C"


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


def test_decimal_overflow_serial_shifts_the_element_field(tmp_path: Path) -> None:
    atom = (
        f"ATOM  {100000:>6}  CA  ALA A   1    "
        f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{10.0:6.2f}{'':10}  \n"
    )
    source = tmp_path / "overflow.pdb"
    target = tmp_path / "fixed.pdb"
    source.write_text(atom + "END\n", encoding="utf-8")

    result = CliRunner().invoke(cli, ["-i", str(source), "-o", str(target)])

    assert result.exit_code == 0, result.output
    line = target.read_text(encoding="utf-8").splitlines()[0]
    assert line[:77] == atom[:77]
    assert line[77:79] == " C"
