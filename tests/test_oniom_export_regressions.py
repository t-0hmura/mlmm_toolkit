"""Scientific-contract regressions for the ONIOM exporter and importer."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.io.oniom_identity import (
    extract_embedded_order_digest,
    format_order_marker,
    pdb_order_digest,
)
from mlmm.workflows import oniom_export
from mlmm.workflows import oniom_import


def test_topology_charge_requires_integer_proximity() -> None:
    accepted = SimpleNamespace(
        atoms=[SimpleNamespace(charge=-4.996)],
    )
    assert oniom_export._get_total_charge(accepted) == -5

    for charge in (0.49, 0.50, float("nan")):
        parm = SimpleNamespace(atoms=[SimpleNamespace(charge=charge)])
        with pytest.raises(RuntimeError, match="Topology total charge"):
            oniom_export._get_total_charge(parm)


def test_partial_model_mapping_is_fatal(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = tmp_path / "model.pdb"
    full = tmp_path / "full.pdb"
    model.write_text("MODEL\n")
    full.write_text("MODEL\n")
    model_atoms = [
        {
            "atom_name": "C1",
            "res_name": "LIG",
            "res_seq": 1,
            "idx": 0,
            "coord": np.array([0.0, 0.0, 0.0]),
            "element": "C",
        },
        {
            "atom_name": "O1",
            "res_name": "LIG",
            "res_seq": 1,
            "idx": 1,
            "coord": np.array([1.0, 0.0, 0.0]),
            "element": "O",
        },
    ]
    input_atoms = [dict(model_atoms[0])]
    monkeypatch.setattr(
        oniom_export,
        "_parse_pdb_atoms_with_meta",
        lambda path: model_atoms if Path(path) == model else input_atoms,
    )

    with pytest.raises(ValueError, match="matched 1/2"):
        oniom_export._read_qm_atoms_from_pdb(
            model,
            input_pdb=full,
            system_coords=np.array([[0.0, 0.0, 0.0]]),
            system_elements=["C"],
        )


def test_missing_gaussian_link_cap_is_fatal(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    atoms = [
        SimpleNamespace(atomic_number=6),
        SimpleNamespace(atomic_number=6),
    ]
    parm = SimpleNamespace(atoms=atoms)
    monkeypatch.setattr(
        oniom_export, "_find_qmmm_boundary_pairs", lambda *_: [(0, 1)],
    )
    monkeypatch.setattr(
        oniom_export, "_estimate_link_h_position_scaled", lambda *_a, **_k: None,
    )

    with pytest.raises(ValueError, match="requires one cap"):
        oniom_export._build_link_atom_specs(
            parm, {0}, elements=["C", "C"], strict=True,
        )


def _torsion_parm(parameter_sets: list[object]) -> SimpleNamespace:
    atom_groups = [
        [
            SimpleNamespace(
                atom_type=name,
                rmin=1.5,
                epsilon=0.1,
                idx=4 * occurrence + offset,
            )
            for offset, name in enumerate(("A", "B", "C", "D"))
        ]
        for occurrence in range(len(parameter_sets))
    ]
    dihedrals = []
    for atoms, parameter in zip(atom_groups, parameter_sets):
        dihedrals.append(
            SimpleNamespace(
                atom1=atoms[0],
                atom2=atoms[1],
                atom3=atoms[2],
                atom4=atoms[3],
                type=parameter,
                improper=False,
            )
        )
    return SimpleNamespace(
        atoms=[atom for group in atom_groups for atom in group],
        bonds=[],
        angles=[],
        dihedrals=dihedrals,
        impropers=[],
    )


def test_repeated_torsion_occurrences_do_not_scale_ambtrs() -> None:
    term = SimpleNamespace(per=3, phase=0.0, phi_k=1.4, div=2.0)
    text = oniom_export._write_gaussian_ff_params(
        _torsion_parm([term, term]),
    )
    line = next(line for line in text.splitlines() if line.startswith("AmbTrs"))
    assert " 0.700000 0.000000 1" in line


def test_conflicting_torsion_parameters_fail_closed() -> None:
    first = SimpleNamespace(per=3, phase=0.0, phi_k=1.4, div=2.0)
    second = SimpleNamespace(per=3, phase=0.0, phi_k=2.0, div=2.0)
    with pytest.raises(ValueError, match="conflicting Amber parameter sets"):
        oniom_export._write_gaussian_ff_params(
            _torsion_parm([first, second]),
        )


def test_separate_periodic_terms_for_one_atom_quartet_are_combined() -> None:
    atoms = [
        SimpleNamespace(
            atom_type=name,
            rmin=1.5,
            epsilon=0.1,
            idx=index,
        )
        for index, name in enumerate(("A", "B", "C", "D"))
    ]
    dihedrals = [
        SimpleNamespace(
            atom1=atoms[0],
            atom2=atoms[1],
            atom3=atoms[2],
            atom4=atoms[3],
            type=SimpleNamespace(per=period, phase=0.0, phi_k=magnitude, div=1.0),
            improper=False,
        )
        for period, magnitude in ((1, 0.4), (2, 0.6))
    ]
    text = oniom_export._write_gaussian_ff_params(
        SimpleNamespace(
            atoms=atoms,
            bonds=[],
            angles=[],
            dihedrals=dihedrals,
            impropers=[],
        )
    )
    line = next(line for line in text.splitlines() if line.startswith("AmbTrs"))
    assert " 0.400000 0.600000 0.000000 0.000000 1" in line


def test_oniom_import_reference_rejects_element_permutation(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference.pdb"
    reference.write_text(
        "HETATM    1  H1  MOL A   1       0.000   0.000   0.000"
        "  1.00  0.00           H\n"
        "HETATM    2  C1  MOL A   1       1.000   0.000   0.000"
        "  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    with pytest.raises(Exception, match="atom-order element mismatch"):
        oniom_import._write_layered_pdb_with_ref(
            tmp_path / "out.pdb",
            reference,
            np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            ["C", "H"],
            {0},
            {0, 1},
        )


def test_oniom_import_reference_infers_protein_ca_as_carbon(tmp_path: Path) -> None:
    reference = tmp_path / "reference.pdb"
    reference.write_text(
        "ATOM      1  CA  GLY A   1       0.000   0.000   0.000"
        "  1.00  0.00              \nEND\n",
        encoding="utf-8",
    )

    verification = oniom_import._write_layered_pdb_with_ref(
        tmp_path / "out.pdb",
        reference,
        np.array([[1.0, 2.0, 3.0]]),
        ["C"],
        {0},
        {0},
    )

    assert verification == "element-verified"


def _write_two_carbon_reference(path: Path, *, reverse: bool = False) -> None:
    rows = [
        oniom_import._format_pdb_atom_line(
            serial=1,
            atom_name="C1",
            res_name="AAA",
            chain_id="A",
            res_seq=1,
            x=0.0,
            y=0.0,
            z=0.0,
            bfac=0.0,
            element="C",
        ),
        oniom_import._format_pdb_atom_line(
            serial=2,
            atom_name="C2",
            res_name="BBB",
            chain_id="B",
            res_seq=2,
            x=1.0,
            y=0.0,
            z=0.0,
            bfac=20.0,
            element="C",
        ),
    ]
    if reverse:
        rows.reverse()
    path.write_text("".join([*rows, "END\n"]), encoding="utf-8")


def test_oniom_reference_digest_rejects_same_element_reorder_even_with_override(
    tmp_path: Path,
) -> None:
    original = tmp_path / "original.pdb"
    swapped = tmp_path / "swapped.pdb"
    _write_two_carbon_reference(original)
    _write_two_carbon_reference(swapped, reverse=True)
    digest = pdb_order_digest(original)

    with pytest.raises(Exception, match="atom identity/order"):
        oniom_import._write_layered_pdb_with_ref(
            tmp_path / "out.pdb",
            swapped,
            np.array([[0.1, 0.0, 0.0], [1.1, 0.0, 0.0]]),
            ["C", "C"],
            {0},
            {0, 1},
            expected_ref_order_digest=digest,
            allow_unverified_ref_order=True,
        )


def test_oniom_reference_digest_accepts_geometry_changes_in_same_order(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference.pdb"
    _write_two_carbon_reference(reference)
    digest = pdb_order_digest(reference)
    verification = oniom_import._write_layered_pdb_with_ref(
        tmp_path / "out.pdb",
        reference,
        np.array([[0.2, 0.3, 0.4], [1.2, 1.3, 1.4]]),
        ["C", "C"],
        {0},
        {0, 1},
        expected_ref_order_digest=digest,
    )
    assert verification == "identity-verified"


def test_markerless_repeated_elements_require_explicit_positional_opt_in(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    reference = tmp_path / "reference.pdb"
    _write_two_carbon_reference(reference)
    args = (
        tmp_path / "out.pdb",
        reference,
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        ["C", "C"],
        {0},
        {0, 1},
    )
    with pytest.raises(Exception, match="cannot be verified"):
        oniom_import._write_layered_pdb_with_ref(*args)
    verification = oniom_import._write_layered_pdb_with_ref(
        *args, allow_unverified_ref_order=True
    )
    assert verification == "unverified-opt-in"
    assert "WARNING" in capsys.readouterr().err


def test_oniom_import_cli_verifies_embedded_reference_identity(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference.pdb"
    source = tmp_path / "job.gjf"
    _write_two_carbon_reference(reference)
    source.write_text(
        "#p oniom(hf/sto-3g:amber)\n\n"
        "ONIOM molecule\n"
        f"{format_order_marker(pdb_order_digest(reference))}\n\n"
        "0 1 0 1 0 1\n"
        "C-CT-0.0 0 0.0 0.0 0.0 H\n"
        "C-CT-0.0 0 1.0 0.0 0.0 L\n\n",
        encoding="utf-8",
    )
    result = CliRunner().invoke(
        root_cli,
        [
            "oniom-import",
            "-i",
            str(source),
            "--ref-pdb",
            str(reference),
            "-o",
            str(tmp_path / "restored"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "ref_order=identity-verified" in result.output
    assert (tmp_path / "restored_layered.pdb").is_file()


def test_embedded_reference_identity_rejects_valid_and_malformed_markers(
    tmp_path: Path,
) -> None:
    source = tmp_path / "job.gjf"
    digest = "a" * 64
    source.write_text(
        f"{format_order_marker(digest)}\n"
        "MLMM_REF_PDB_ORDER_V1_SHA256=malformed\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="duplicated or conflicting"):
        extract_embedded_order_digest(source)


def test_oniom_import_public_pair_rolls_back_on_primary_failure(
    tmp_path: Path, monkeypatch
) -> None:
    from mlmm.core import result_commit

    source = tmp_path / "job.gjf"
    source.write_text(
        "#p oniom(hf/sto-3g:amber)\n\n"
        "ONIOM molecule\n\n"
        "0 1 0 1 0 1\n"
        "C-CT-0.0 0 0.0 0.0 0.0 H\n"
        "H-HC-0.0 0 1.0 0.0 0.0 L\n\n",
        encoding="utf-8",
    )
    prefix = tmp_path / "restored"
    xyz = prefix.with_suffix(".xyz")
    layered = tmp_path / "restored_layered.pdb"
    xyz.write_bytes(b"old xyz\n")
    layered.write_bytes(b"old pdb\n")
    replace = result_commit._replace_exact
    failed = False

    def fail_primary_once(staged: Path, destination: Path) -> None:
        nonlocal failed
        if destination == xyz and not failed:
            failed = True
            raise OSError("injected xyz publication failure")
        replace(staged, destination)

    monkeypatch.setattr(result_commit, "_replace_exact", fail_primary_once)
    result = CliRunner().invoke(
        root_cli,
        ["oniom-import", "-i", str(source), "-o", str(prefix)],
    )

    assert result.exit_code != 0
    assert xyz.read_bytes() == b"old xyz\n"
    assert layered.read_bytes() == b"old pdb\n"
