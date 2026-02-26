"""Boundary and utility tests for mlmm_toolkit.extract."""

from __future__ import annotations

import io
import sys

import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm_toolkit requires Python >= 3.11",
)


def _atom_line(
    serial: int,
    atom: str,
    resname: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str,
    record: str = "ATOM",
) -> str:
    return (
        f"{record:<6}{serial:>5} {atom:<4} {resname:>3} {chain:1}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{20.00:>6.2f}          {element:>2}\n"
    )


def _sample_pdb(*, swap_first_two_atoms: bool = False, drop_last_atom: bool = False) -> str:
    lines = [
        _atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N", "ATOM"),
        _atom_line(2, "CA", "ALA", "A", 1, 1.46, 0.0, 0.0, "C", "ATOM"),
        _atom_line(3, "C1", "GPP", "A", 2, 2.50, 0.0, 0.0, "C", "HETATM"),
        _atom_line(4, "ZN", "ZN", "A", 3, 4.50, 0.0, 0.0, "ZN", "HETATM"),
    ]
    if swap_first_two_atoms:
        lines[0], lines[1] = lines[1], lines[0]
    if drop_last_atom:
        lines = lines[:-1]
    lines.extend(["TER\n", "END\n"])
    return "".join(lines)


def _parse_structure(pdb_text: str, name: str):
    from Bio import PDB

    parser = PDB.PDBParser(QUIET=True)
    return parser.get_structure(name, io.StringIO(pdb_text))


def test_parse_res_tokens_accepts_supported_formats():
    from mlmm_toolkit.extract import _parse_res_tokens

    parsed = _parse_res_tokens("123,123A,A:45,B:67C")
    assert parsed == [
        (None, 123, None),
        (None, 123, "A"),
        ("A", 45, None),
        ("B", 67, "C"),
    ]


@pytest.mark.parametrize("spec", ["", " ", "A:", "A-10", "A:1:2", "foo"])
def test_parse_res_tokens_rejects_invalid(spec):
    from mlmm_toolkit.extract import _parse_res_tokens

    with pytest.raises(ValueError):
        _parse_res_tokens(spec)


def test_parse_ligand_charge_option_numeric_and_mapping():
    from mlmm_toolkit.extract import _parse_ligand_charge_option

    total, mapping = _parse_ligand_charge_option("-3")
    assert total == -3.0
    assert mapping is None

    total, mapping = _parse_ligand_charge_option("gpp:-3, mmt:-1")
    assert total is None
    assert mapping == {"GPP": -3.0, "MMT": -1.0}


def test_parse_ligand_charge_option_rejects_invalid_mapping():
    from mlmm_toolkit.extract import _parse_ligand_charge_option

    with pytest.raises(ValueError, match="Invalid --ligand-charge token"):
        _parse_ligand_charge_option("GPP")


def test_strip_trailing_end_normalizes_newline():
    from mlmm_toolkit.extract import _strip_trailing_END

    text = "ATOM\nEND\nEND\n"
    out = _strip_trailing_END(text)
    assert out == "ATOM\n"


def test_assert_atom_ordering_identical_passes_for_identical_structures():
    from mlmm_toolkit.extract import _assert_atom_ordering_identical

    s1 = _parse_structure(_sample_pdb(), "s1")
    s2 = _parse_structure(_sample_pdb(), "s2")
    _assert_atom_ordering_identical([s1, s2])


def test_assert_atom_ordering_identical_raises_for_mismatch():
    from mlmm_toolkit.extract import _assert_atom_ordering_identical

    s1 = _parse_structure(_sample_pdb(), "s1")
    s2 = _parse_structure(_sample_pdb(swap_first_two_atoms=True), "s2")

    with pytest.raises(ValueError, match="Atom order mismatch"):
        _assert_atom_ordering_identical([s1, s2])


def test_assert_atom_ordering_identical_raises_for_atom_count_mismatch():
    from mlmm_toolkit.extract import _assert_atom_ordering_identical

    s1 = _parse_structure(_sample_pdb(), "s1")
    s2 = _parse_structure(_sample_pdb(drop_last_atom=True), "s2")

    with pytest.raises(ValueError, match="Atom count mismatch"):
        _assert_atom_ordering_identical([s1, s2])


def test_compute_charge_summary_total_charge_distribution():
    from mlmm_toolkit.extract import compute_charge_summary

    structure = _parse_structure(_sample_pdb(), "s")
    all_res = list(structure.get_residues())
    selected_ids = {r.get_full_id() for r in all_res}
    substrate_ids = {r.get_full_id() for r in all_res if r.get_resname().upper() == "GPP"}

    summary = compute_charge_summary(
        structure,
        selected_ids,
        substrate_ids,
        ligand_charge=-2.0,
    )

    assert summary["protein_charge"] == 0.0
    assert summary["ligand_total_charge"] == -2.0
    assert summary["ion_total_charge"] == 2.0
    assert summary["total_charge"] == 0.0
    assert summary["unknown_residue_charges"]["GPP"] == -2.0


def test_compute_charge_summary_mapping_mode():
    from mlmm_toolkit.extract import compute_charge_summary

    structure = _parse_structure(_sample_pdb(), "s")
    all_res = list(structure.get_residues())
    selected_ids = {r.get_full_id() for r in all_res}
    substrate_ids = {r.get_full_id() for r in all_res if r.get_resname().upper() == "GPP"}

    summary = compute_charge_summary(
        structure,
        selected_ids,
        substrate_ids,
        ligand_charge={"GPP": -3},
    )

    assert summary["ligand_total_charge"] == -3.0
    assert summary["ion_total_charge"] == 2.0
    assert summary["total_charge"] == -1.0
