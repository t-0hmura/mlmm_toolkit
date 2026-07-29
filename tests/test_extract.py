"""Boundary and utility tests for mlmm.workflows.extract."""

from __future__ import annotations

import io
import sys

import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm requires Python >= 3.11",
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
    from mlmm.workflows.extract import _parse_res_tokens

    parsed = _parse_res_tokens("123,123A,A:45,B:67C")
    assert parsed == [
        (None, 123, None),
        (None, 123, "A"),
        ("A", 45, None),
        ("B", 67, "C"),
    ]


@pytest.mark.parametrize("spec", ["", " ", "A:", "A-10", "A:1:2", "foo"])
def test_parse_res_tokens_rejects_invalid(spec):
    from mlmm.workflows.extract import _parse_res_tokens

    # ValueError is intentional — `resolve_substrate_residues` upstream
    # catches this to fall back to resname-based dispatch. CLI users still
    # see a clean error via render_cli_exception at the cli tail.
    with pytest.raises(ValueError):
        _parse_res_tokens(spec)


def test_parse_ligand_charge_option_numeric_and_mapping():
    from mlmm.workflows.extract import _parse_ligand_charge_option

    total, mapping = _parse_ligand_charge_option("-3")
    assert total == -3.0
    assert mapping is None

    total, mapping = _parse_ligand_charge_option("gpp:-3, mmt:-1")
    assert total is None
    assert mapping == {"GPP": -3.0, "MMT": -1.0}


def test_parse_ligand_charge_option_rejects_invalid_mapping():
    import click
    from mlmm.workflows.extract import _parse_ligand_charge_option

    # Invalid selectors use click.BadParameter for a clean command-line error.
    with pytest.raises(click.BadParameter, match="Invalid --ligand-charge token"):
        _parse_ligand_charge_option("GPP")


def test_strip_trailing_end_normalizes_newline():
    from mlmm.workflows.extract import _strip_trailing_END

    text = "ATOM\nEND\nEND\n"
    out = _strip_trailing_END(text)
    assert out == "ATOM\n"


def test_assert_atom_ordering_identical_passes_for_identical_structures():
    from mlmm.workflows.extract import _assert_atom_ordering_identical

    s1 = _parse_structure(_sample_pdb(), "s1")
    s2 = _parse_structure(_sample_pdb(), "s2")
    _assert_atom_ordering_identical([s1, s2])


def test_assert_atom_ordering_identical_raises_for_mismatch():
    from mlmm.workflows.extract import _assert_atom_ordering_identical

    s1 = _parse_structure(_sample_pdb(), "s1")
    s2 = _parse_structure(_sample_pdb(swap_first_two_atoms=True), "s2")

    with pytest.raises(ValueError, match="Atom order mismatch"):
        _assert_atom_ordering_identical([s1, s2])


def test_assert_atom_ordering_identical_raises_for_atom_count_mismatch():
    from mlmm.workflows.extract import _assert_atom_ordering_identical

    s1 = _parse_structure(_sample_pdb(), "s1")
    s2 = _parse_structure(_sample_pdb(drop_last_atom=True), "s2")

    with pytest.raises(ValueError, match="Atom count mismatch"):
        _assert_atom_ordering_identical([s1, s2])


def test_compute_charge_summary_total_charge_distribution():
    from mlmm.workflows.extract import compute_charge_summary

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
    from mlmm.workflows.extract import compute_charge_summary

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


def test_compute_charge_summary_does_not_override_known_ion(monkeypatch):
    from mlmm.workflows import extract as extract_module

    structure = _parse_structure(_sample_pdb(), "s")
    all_res = list(structure.get_residues())
    selected_ids = {r.get_full_id() for r in all_res}
    substrate_ids = {
        r.get_full_id() for r in all_res if r.get_resname().upper() == "GPP"
    }
    warnings = []
    monkeypatch.setattr(
        extract_module,
        "_echo_warning",
        lambda msg, *args: warnings.append(
            extract_module._format_echo_message(msg, *args)
        ),
    )

    summary = extract_module.compute_charge_summary(
        structure,
        selected_ids,
        substrate_ids,
        ligand_charge="GPP:-3,ZN:3",
    )

    assert summary["ligand_total_charge"] == -3.0
    assert summary["ion_total_charge"] == 2.0
    assert summary["total_charge"] == -1.0
    assert len(warnings) == 1
    assert "ZN" in warnings[0]
    assert "matched no unknown" in warnings[0]


def test_compute_charge_summary_accepts_exact_known_ion_mapping(monkeypatch):
    from mlmm.workflows import extract as extract_module

    structure = _parse_structure(_sample_pdb(), "s")
    all_res = list(structure.get_residues())
    selected_ids = {r.get_full_id() for r in all_res}
    substrate_ids = {
        r.get_full_id() for r in all_res if r.get_resname().upper() == "GPP"
    }
    warnings = []
    monkeypatch.setattr(
        extract_module,
        "_echo_warning",
        lambda msg, *args: warnings.append(
            extract_module._format_echo_message(msg, *args)
        ),
    )

    summary = extract_module.compute_charge_summary(
        structure,
        selected_ids,
        substrate_ids,
        ligand_charge="GPP:-3,ZN:2",
    )

    assert summary["total_charge"] == -1.0
    assert warnings == []


def test_compute_charge_summary_terminus_cap_charges():
    """A kept terminal cap carries the ionized-terminus formal charge the internal
    residue charge omits: C-terminus carboxylate (OXT) -> -1, N-terminus ammonium
    (H1/H2/H3) -> +1. No correction unless the cap is kept (keep_ncap/ccap)."""
    from mlmm.workflows.extract import compute_charge_summary

    lines = [
        # chain A: C-terminal ALA (has OXT)
        _atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        _atom_line(2, "CA", "ALA", "A", 1, 1.5, 0.0, 0.0, "C"),
        _atom_line(3, "C", "ALA", "A", 1, 2.0, 1.4, 0.0, "C"),
        _atom_line(4, "O", "ALA", "A", 1, 1.3, 2.4, 0.0, "O"),
        _atom_line(5, "OXT", "ALA", "A", 1, 3.3, 1.4, 0.0, "O"),
        _atom_line(6, "CB", "ALA", "A", 1, 2.1, -1.2, 0.0, "C"),
        "TER\n",
        # chain B: N-terminal ALA (has H1/H2/H3)
        _atom_line(7, "N", "ALA", "B", 1, 10.0, 0.0, 0.0, "N"),
        _atom_line(8, "H1", "ALA", "B", 1, 9.5, 0.8, 0.0, "H"),
        _atom_line(9, "H2", "ALA", "B", 1, 9.5, -0.8, 0.0, "H"),
        _atom_line(10, "H3", "ALA", "B", 1, 10.8, 0.0, 0.0, "H"),
        _atom_line(11, "CA", "ALA", "B", 1, 11.5, 0.0, 0.0, "C"),
        _atom_line(12, "C", "ALA", "B", 1, 12.0, 1.4, 0.0, "C"),
        _atom_line(13, "O", "ALA", "B", 1, 11.3, 2.4, 0.0, "O"),
        _atom_line(14, "CB", "ALA", "B", 1, 12.1, -1.2, 0.0, "C"),
        "TER\n", "END\n",
    ]
    structure = _parse_structure("".join(lines), "term")
    res = list(structure.get_residues())
    cterm = next(r.get_full_id() for r in res if r.get_parent().id == "A")
    nterm = next(r.get_full_id() for r in res if r.get_parent().id == "B")
    sel = {r.get_full_id() for r in res}

    # caps not kept -> both ALA neutral (regression target: previously dropped)
    assert compute_charge_summary(structure, sel, set())["protein_charge"] == 0.0
    # C-cap kept (OXT) -> -1 carboxylate
    assert compute_charge_summary(structure, sel, set(), keep_ccap_ids={cterm})["protein_charge"] == -1.0
    # N-cap kept (NH3+) -> +1 ammonium
    assert compute_charge_summary(structure, sel, set(), keep_ncap_ids={nterm})["protein_charge"] == 1.0
    # both kept -> net 0
    assert compute_charge_summary(
        structure, sel, set(), keep_ncap_ids={nterm}, keep_ccap_ids={cterm}
    )["protein_charge"] == 0.0
def test_chain_resname_selector_limits_repeated_ligands_to_requested_chain():
    import io
    from Bio import PDB
    from mlmm.workflows.extract import resolve_substrate_residues

    pdb_text = (
        "HETATM    1  C1  SAM A  10       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  C1  SAM A  20       1.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    3  C1  SAM B  10       2.000   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )
    structure = PDB.PDBParser(QUIET=True).get_structure(
        "repeated", io.StringIO(pdb_text)
    )
    assert [r.id[1] for r in resolve_substrate_residues(structure, "A:SAM")] == [10, 20]
    exact = resolve_substrate_residues(structure, "A:SAM:20")
    assert len(exact) == 1 and exact[0].get_parent().id == "A" and exact[0].id[1] == 20


def test_chain_resname_parser_accepts_ccd_and_negative_ids():
    from mlmm.workflows.extract import _parse_chain_resname_tokens

    assert _parse_chain_resname_tokens("A:1AB:-2C") == [("A", "1AB", -2, "C")]


def test_valid_selector_not_found_is_not_reinterpreted_as_resname():
    import io
    import pytest
    from Bio import PDB
    from mlmm.workflows.extract import resolve_substrate_residues

    structure = PDB.PDBParser(QUIET=True).get_structure(
        "single",
        io.StringIO(
            "HETATM    1  C1  SAM A  10       0.000   0.000   0.000  1.00  0.00           C\nEND\n"
        ),
    )
    with pytest.raises(ValueError, match="A:999"):
        resolve_substrate_residues(structure, "A:999")
    with pytest.raises(ValueError, match="A:SAM:999"):
        resolve_substrate_residues(structure, "A:SAM:999")


def test_extract_cli_executes_real_command_path(tmp_path):
    from click.testing import CliRunner
    from mlmm.cli import cli as root_cli

    input_path = tmp_path / "complex.pdb"
    output_path = tmp_path / "pocket.pdb"
    input_path.write_text(_sample_pdb(), encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            "extract",
            "-i",
            str(input_path),
            "-c",
            "GPP",
            "-l",
            "GPP:-3",
            "--no-add-linkh",
            "-o",
            str(output_path),
        ],
    )

    assert result.exit_code == 0, result.output
    assert output_path.is_file()
    assert "GPP" in output_path.read_text(encoding="utf-8")
