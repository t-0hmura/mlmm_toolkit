"""Inspectable ML-region structure output contracts."""

from pathlib import Path
from types import SimpleNamespace

from ase import Atoms

from mlmm.workflows.dft import (
    write_ml_region_pdb_pair,
    write_ml_region_xyz_pair,
)


def test_write_ml_region_xyz_pair_distinguishes_generated_link_h(
    tmp_path: Path,
) -> None:
    workspace = SimpleNamespace(
        atoms_model=Atoms("CN", positions=[[0.0, 0.0, 0.0], [1.3, 0.0, 0.0]]),
        atoms_model_lh=Atoms(
            "CNH",
            positions=[
                [0.0, 0.0, 0.0],
                [1.3, 0.0, 0.0],
                [-1.09, 0.0, 0.0],
            ],
        ),
        link_pairs=[(1, 3)],
    )

    without_link, with_link = write_ml_region_xyz_pair(workspace, tmp_path)

    assert without_link.name == "ml_region_without_linkH.xyz"
    assert with_link.name == "ml_region_with_linkH.xyz"
    assert without_link.read_text(encoding="utf-8").splitlines()[0] == "2"
    with_lines = with_link.read_text(encoding="utf-8").splitlines()
    assert with_lines[0] == "3"
    assert (
        with_lines[1]
        == "ML region with 1 link H; full-system ML-MM pairs: 1-3"
    )
    assert with_lines[-1].split()[0] == "H"


def test_write_ml_region_pdb_pair_preserves_model_identity_and_marks_link_h(
    tmp_path: Path,
) -> None:
    model_pdb = tmp_path / "model.pdb"
    model_pdb.write_text(
        "ATOM      7  C1  LIG A  42       0.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      9  N1  LIG A  42       1.300   0.000   0.000  1.00  0.00           N  \n"
        "END\n",
        encoding="utf-8",
    )
    workspace = SimpleNamespace(
        model_pdb=model_pdb,
        atoms_model=Atoms(
            "CN",
            positions=[[0.1, 0.2, 0.3], [1.4, 0.2, 0.3]],
        ),
        atoms_model_lh=Atoms(
            "CNH",
            positions=[
                [0.1, 0.2, 0.3],
                [1.4, 0.2, 0.3],
                [-1.0, 0.2, 0.3],
            ],
        ),
        link_pairs=[(1, 3)],
    )
    xyz_paths = write_ml_region_xyz_pair(workspace, tmp_path)

    without_pdb, with_pdb = write_ml_region_pdb_pair(
        workspace,
        tmp_path,
        xyz_paths=xyz_paths,
    )

    without_text = without_pdb.read_text(encoding="utf-8")
    with_text = with_pdb.read_text(encoding="utf-8")
    assert " C1  LIG A  42" in without_text
    assert " N1  LIG A  42" in without_text
    assert "   0.100   0.200   0.300" in without_text
    link_h_records = [
        line
        for line in with_text.splitlines()
        if line.startswith("HETATM") and line[12:16].strip() == "HL"
    ]
    assert len(link_h_records) == 1
    assert link_h_records[0][17:20] == "LKH"
    assert link_h_records[0][21:22] == "L"
    assert link_h_records[0][22:26].strip() == "1"
    assert sum(
        line.startswith(("ATOM  ", "HETATM"))
        for line in with_text.splitlines()
    ) == 3
