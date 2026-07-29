from __future__ import annotations

from pathlib import Path

from mlmm.core.utils import load_pdb_atom_metadata
from mlmm.workflows.all import _parse_scan_lists_literals


BEZA = Path(__file__).parents[1] / "examples" / "beza"
STATES = ("1.R.pdb", "2.IM.pdb", "3.P.pdb")
IDENTITY_FIELDS = (
    "is_hetatm",
    "chain",
    "resname",
    "resseq",
    "icode",
    "name",
    "altloc",
    "element",
)
SCAN_STAGES = (
    '[("CS1 SAM 320","C7 GPP 321",1.50),'
    '("CS1 SAM 320","SD SAM 320",3.30)]',
    '[("C7 GPP 321","H11 GPP 321",2.90),'
    '("OE2 GLU 186","H11 GPP 321",1.00)]',
)


def _identity(path: Path) -> list[tuple[object, ...]]:
    return [
        tuple(atom.get(field) for field in IDENTITY_FIELDS)
        for atom in load_pdb_atom_metadata(path)
    ]


def _coordinates(path: Path) -> list[str]:
    return [
        line[30:54]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]


def test_beza_states_share_one_ordered_topology() -> None:
    identities = [_identity(BEZA / name) for name in STATES]
    assert len(identities[0]) == 9215
    assert identities[1:] == [identities[0], identities[0]]

    residues = {identity[2] for identity in identities[0]}
    assert {"SAM", "GPP", "MG", "GLU"} <= residues


def test_beza_states_have_distinct_coordinates() -> None:
    coordinates = [_coordinates(BEZA / name) for name in STATES]
    assert coordinates[0] != coordinates[1]
    assert coordinates[1] != coordinates[2]


def test_beza_scan_selectors_resolve_against_reactant() -> None:
    metadata = load_pdb_atom_metadata(BEZA / "1.R.pdb")
    stages = _parse_scan_lists_literals(
        SCAN_STAGES,
        atom_meta=metadata,
        one_based=True,
    )
    assert [len(stage) for stage in stages] == [2, 2]
