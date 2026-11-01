"""``mm-parm --add-h``: PDBFixer's His and N-terminal hydrogen names must be ones tleap accepts.

PDBFixer writes every histidine as HIS, which tleap reads as HIE, so a His carrying HD1
stops tleap with "HD1 does not have a type". It also writes the N-terminal NH3+ as
H/H2/H3, while Amber's N-terminal templates use H1/H2/H3.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mlmm.workflows import mm_parm


def _atom(serial: int, name: str, resname: str, resseq: int, element: str) -> str:
    padded = name if len(name) == 4 else f" {name:<3s}"
    return (
        f"ATOM  {serial:5d} {padded} {resname:>3s} A{resseq:4d}    "
        f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          {element:>2s}\n"
    )


def _pdb(residues) -> str:
    """Build a PDB from ``[(resname, [(atom name, element), ...]), ...]``."""
    lines, serial = [], 0
    for resseq, (resname, atoms) in enumerate(residues, start=1):
        for name, element in atoms:
            serial += 1
            lines.append(_atom(serial, name, resname, resseq, element))
    return "".join(lines) + "TER\nEND\n"


def _names(path: Path):
    return [
        (line[17:20], int(line[22:26]), line[12:16].strip())
        for line in path.read_text().splitlines()
        if line.startswith("ATOM")
    ]


HIS_RING = [("ND1", "N"), ("CE1", "C"), ("HE1", "H"), ("NE2", "N"), ("CD2", "C"), ("HD2", "H")]
BACKBONE = [("N", "N"), ("CA", "C"), ("HA", "H"), ("C", "C"), ("O", "O")]


def test_his_and_n_terminal_hydrogens_get_amber_names(tmp_path: Path) -> None:
    source = _pdb(
        [
            ("GLY", [("N", "N"), ("H", "H"), ("H2", "H"), ("H3", "H"), ("CA", "C"), ("C", "C")]),
            ("HIS", BACKBONE + [("H", "H")] + HIS_RING + [("HD1", "H")]),
            ("HIS", BACKBONE + [("H", "H")] + HIS_RING + [("HE2", "H")]),
            ("HIS", BACKBONE + [("H", "H")] + HIS_RING + [("HD1", "H"), ("HE2", "H")]),
            ("HIS", BACKBONE + [("H", "H")] + HIS_RING),
            ("ALA", BACKBONE + [("H", "H")]),
        ]
    )
    pdb = tmp_path / "fixed.pdb"
    pdb.write_text(source)

    mm_parm.rename_pdbfixer_hydrogens_for_amber(pdb)

    names = _names(pdb)
    resnames = {resseq: resname for resname, resseq, _ in names}
    assert [resnames[i] for i in range(2, 7)] == ["HID", "HIE", "HIP", "HIS", "ALA"]
    # N-terminal NH3+ is H1/H2/H3; every other backbone amide H keeps the name H.
    assert [n for r, s, n in names if s == 1 and n.startswith("H")] == ["H1", "H2", "H3"]
    for resseq in range(2, 7):
        assert ("H" in [n for _, s, n in names if s == resseq])
    assert "H1" not in [n for _, s, n in names if s > 1]
    # Nothing but the residue name and the atom name changes.
    for before, after in zip(source.splitlines(), pdb.read_text().splitlines()):
        assert before[:12] == after[:12] and before[20:] == after[20:]
        assert len(before) == len(after)


def test_non_amino_acid_residues_and_already_amber_names_are_left_alone(tmp_path: Path) -> None:
    source = _pdb(
        [
            ("LIG", [("N", "N"), ("H", "H"), ("H2", "H"), ("H3", "H"), ("C1", "C")]),
            ("HOH", [("O", "O"), ("H1", "H"), ("H2", "H")]),
            ("GLY", [("N", "N"), ("H1", "H"), ("H2", "H"), ("H3", "H"), ("CA", "C")]),
            ("HIE", BACKBONE + [("H", "H")] + HIS_RING + [("HE2", "H")]),
        ]
    )
    pdb = tmp_path / "fixed.pdb"
    pdb.write_text(source)

    mm_parm.rename_pdbfixer_hydrogens_for_amber(pdb)

    assert pdb.read_text() == source


HEAVY_TRIPEPTIDE = """\
ATOM      1  N   GLY A   1       3.326   1.548   0.000  1.00  0.00           N
ATOM      2  CA  GLY A   1       3.970   2.846   0.000  1.00  0.00           C
ATOM      3  C   GLY A   1       5.484   2.687   0.000  1.00  0.00           C
ATOM      4  O   GLY A   1       5.993   1.568   0.000  1.00  0.00           O
ATOM      5  N   HIS A   2       6.203   3.811   0.000  1.00  0.00           N
ATOM      6  CA  HIS A   2       7.652   3.794   0.000  1.00  0.00           C
ATOM      7  CB  HIS A   2       8.192   3.076  -1.232  1.00  0.00           C
ATOM      8  CG  HIS A   2       9.698   3.006  -1.321  1.00  0.00           C
ATOM      9  ND1 HIS A   2      10.345   2.385  -2.383  1.00  0.00           N
ATOM     10  CE1 HIS A   2      11.646   2.485  -2.183  1.00  0.00           C
ATOM     11  NE2 HIS A   2      11.856   3.129  -1.062  1.00  0.00           N
ATOM     12  CD2 HIS A   2      10.670   3.472  -0.492  1.00  0.00           C
ATOM     13  C   HIS A   2       8.217   5.207   0.000  1.00  0.00           C
ATOM     14  O   HIS A   2       7.466   6.179   0.000  1.00  0.00           O
ATOM     15  N   ALA A   3       9.547   5.319   0.000  1.00  0.00           N
ATOM     16  CA  ALA A   3      10.207   6.608   0.000  1.00  0.00           C
ATOM     17  CB  ALA A   3       9.824   7.421   1.232  1.00  0.00           C
ATOM     18  C   ALA A   3      11.721   6.449   0.000  1.00  0.00           C
ATOM     19  O   ALA A   3      12.231   5.331   0.000  1.00  0.00           O
ATOM     20  OXT ALA A   3      12.440   7.574   0.000  1.00  0.00           O
TER
END
"""


@pytest.mark.parametrize("ph", [5.0, 7.0])
def test_pdbfixer_output_is_renamed_for_tleap(tmp_path: Path, ph: float) -> None:
    pytest.importorskip("pdbfixer")
    source = tmp_path / "heavy.pdb"
    source.write_text(HEAVY_TRIPEPTIDE)
    fixed = tmp_path / "fixed.pdb"

    mm_parm.add_hydrogens_with_pdbfixer(source, fixed, ph)

    names = _names(fixed)
    his = {resname for resname, resseq, _ in names if resseq == 2}
    assert his <= {"HID", "HIE", "HIP"}
    assert [n for r, s, n in names if s == 1 and n.startswith("H") and not n.startswith("HA")] == [
        "H1", "H2", "H3",
    ]
    assert "H" in [n for _, s, n in names if s == 2]
