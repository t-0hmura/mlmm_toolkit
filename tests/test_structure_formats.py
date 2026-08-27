"""Tests for mmCIF and PDB-size-limit-safe structure I/O."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest


def _write_minimal_cif(path: Path) -> None:
    path.write_text(
        """data_input
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
HETATM 1 C C1 . SAM LONG_CHAIN 1 . ? 0.0 1.0 2.0 1.00 12.0 . 10001 SAM LONG_CHAIN C1 1
HETATM 2 O O1 . SAM LONG_CHAIN 1 . ? 1.0 1.0 2.0 1.00 13.0 . 10001 SAM LONG_CHAIN O1 1
#
""",
        encoding="utf-8",
    )


def test_cif_is_normalized_and_output_restores_auth_ids(tmp_path: Path) -> None:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from mlmm.core.utils import (
        convert_xyz_to_pdb,
        load_pdb_atom_metadata,
        prepare_input_structure,
    )

    source = tmp_path / "large-id.cif"
    _write_minimal_cif(source)
    xyz = tmp_path / "optimized.xyz"
    xyz.write_text("2\nframe\nC 3.0 4.0 5.0\nO 6.0 7.0 8.0\n", encoding="utf-8")
    out_pdb = tmp_path / "optimized.pdb"

    prepared = prepare_input_structure(source)
    try:
        assert prepared.is_cif
        assert prepared.source_path.suffix == ".pdb"
        metadata = load_pdb_atom_metadata(prepared.source_path)
        assert metadata[0]["chain"] == "LONG_CHAIN"
        assert metadata[0]["resseq"] == 10001

        convert_xyz_to_pdb(xyz, prepared.source_path, out_pdb)
        out_cif = out_pdb.with_suffix(".cif")
        assert out_cif.exists()
        data = MMCIF2Dict(str(out_cif))
        assert data["_atom_site.auth_asym_id"] == ["LONG_CHAIN", "LONG_CHAIN"]
        assert data["_atom_site.auth_seq_id"] == ["10001", "10001"]
        assert np.allclose([float(value) for value in data["_atom_site.Cartn_x"]], [3.0, 6.0])
    finally:
        prepared.cleanup()


def test_pdb_output_rejects_coordinate_field_overflow(tmp_path: Path) -> None:
    from mlmm.core.utils import convert_xyz_to_pdb

    ref = tmp_path / "ref.pdb"
    ref.write_text(
        "ATOM      1  C   MOL A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    xyz = tmp_path / "far.xyz"
    xyz.write_text("1\nframe\nC 10000.0 0.0 0.0\n", encoding="utf-8")

    with pytest.raises(ValueError, match="fixed-column PDB range"):
        convert_xyz_to_pdb(xyz, ref, tmp_path / "out.pdb")


def test_pdb_output_rejects_when_no_xyz_frame_matches_topology(tmp_path: Path) -> None:
    from mlmm.core.utils import convert_xyz_to_pdb

    ref = tmp_path / "ref.pdb"
    ref.write_text(
        "ATOM      1  C   MOL A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    xyz = tmp_path / "wrong-count.xyz"
    xyz.write_text("2\nframe\nC 0 0 0\nH 1 0 0\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Atom count mismatch"):
        convert_xyz_to_pdb(xyz, ref, tmp_path / "out.pdb")


def test_all_materializes_xyz_coordinates_on_reference_topology(
    tmp_path: Path,
) -> None:
    from mlmm.core.utils import (
        apply_ref_pdb_override,
        load_pdb_atom_metadata,
        prepare_input_structure,
    )
    from mlmm.workflows.all import _materialize_all_coordinate_inputs

    ref = tmp_path / "ref.pdb"
    ref.write_text(
        "HETATM    1  C1  LIG A   7       0.000   0.000   0.000  1.00 10.00           C\n"
        "HETATM    2  O1  LIG A   7       1.000   0.000   0.000  1.00 20.00           O\n"
        "END\n",
        encoding="utf-8",
    )
    xyz = tmp_path / "endpoint.xyz"
    xyz.write_text(
        "2\nendpoint\nC 3.125 4.250 5.375\nO 6.500 7.625 8.750\n",
        encoding="utf-8",
    )
    prepared = prepare_input_structure(xyz)
    try:
        apply_ref_pdb_override(prepared, ref)
        (materialized,) = _materialize_all_coordinate_inputs(
            [prepared],
            tmp_path / "work",
        )

        metadata = load_pdb_atom_metadata(materialized)
        atom_lines = [
            line
            for line in materialized.read_text(encoding="utf-8").splitlines()
            if line.startswith(("ATOM", "HETATM"))
        ]
        assert metadata[0]["chain"] == "A"
        assert metadata[0]["resname"] == "LIG"
        assert [float(atom_lines[0][30:38]), float(atom_lines[1][30:38])] == [
            3.125,
            6.5,
        ]
        assert [float(atom_lines[0][60:66]), float(atom_lines[1][60:66])] == [
            10.0,
            20.0,
        ]
    finally:
        prepared.cleanup()


def test_all_coordinate_materialization_preserves_colliding_reference(
    tmp_path: Path,
) -> None:
    import click

    from mlmm.core.utils import apply_ref_pdb_override, prepare_input_structure
    from mlmm.workflows.all import _materialize_all_coordinate_inputs

    work_dir = tmp_path / "work"
    ref = work_dir / "coordinate_inputs" / "endpoint_01.pdb"
    ref.parent.mkdir(parents=True)
    original = (
        "HETATM    1  C1  LIG A   7       0.000   0.000   0.000  1.00 10.00           C\n"
        "END\n"
    )
    ref.write_text(original, encoding="utf-8")
    xyz = tmp_path / "endpoint.xyz"
    xyz.write_text("1\nendpoint\nC 3.125 4.250 5.375\n", encoding="utf-8")
    prepared = prepare_input_structure(xyz)
    try:
        apply_ref_pdb_override(prepared, ref)
        with pytest.raises(click.BadParameter, match="managed all-workflow"):
            _materialize_all_coordinate_inputs([prepared], work_dir)
        assert ref.read_text(encoding="utf-8") == original
    finally:
        prepared.cleanup()


def test_all_coordinate_materialization_preserves_other_endpoint_input(
    tmp_path: Path,
) -> None:
    import click

    from mlmm.core.utils import apply_ref_pdb_override, prepare_input_structure
    from mlmm.workflows.all import _materialize_all_coordinate_inputs

    work_dir = tmp_path / "work"
    first_input = work_dir / "coordinate_inputs" / "endpoint_02.pdb"
    first_input.parent.mkdir(parents=True)
    original = (
        "HETATM    1  C1  LIG A   7       0.000   0.000   0.000  1.00 10.00           C\n"
        "END\n"
    )
    first_input.write_text(original, encoding="utf-8")
    ref = tmp_path / "ref.pdb"
    ref.write_text(original, encoding="utf-8")
    xyz = tmp_path / "endpoint.xyz"
    xyz.write_text("1\nendpoint\nC 3.125 4.250 5.375\n", encoding="utf-8")
    prepared = [
        prepare_input_structure(first_input),
        prepare_input_structure(xyz),
    ]
    try:
        apply_ref_pdb_override(prepared[1], ref)
        with pytest.raises(click.BadParameter, match="managed all-workflow"):
            _materialize_all_coordinate_inputs(prepared, work_dir)
        assert first_input.read_text(encoding="utf-8") == original
    finally:
        for item in prepared:
            item.cleanup()


@pytest.mark.parametrize(
    ("bad_frame", "message"),
    [
        pytest.param(
            "1\nbad count\nC 2 0 0\n",
            "Atom count mismatch",
            id="late-count",
        ),
        pytest.param(
            "2\nbad order\nO 2 0 0\nC 3 0 0\n",
            "Ordered elements differ",
            id="late-order",
        ),
        pytest.param(
            "2\nbad width\nC 10000 0 0\nO 3 0 0\n",
            "fixed-column PDB range",
            id="late-width",
        ),
        pytest.param(
            "2\nbad nan\nC nan 0 0\nO 3 0 0\n",
            "non-finite",
            id="late-nan",
        ),
        pytest.param(
            "2\nbad positive infinity\nC inf 0 0\nO 3 0 0\n",
            "non-finite",
            id="late-positive-inf",
        ),
        pytest.param(
            "2\nbad negative infinity\nC -inf 0 0\nO 3 0 0\n",
            "non-finite",
            id="late-negative-inf",
        ),
    ],
)
def test_xyz_overlay_prevalidates_every_frame_before_mutating_destinations(
    tmp_path: Path,
    bad_frame: str,
    message: str,
) -> None:
    from dataclasses import replace

    from mlmm.core.utils import convert_xyz_to_pdb, prepare_input_structure
    from mlmm.io.structure_formats import (
        coordinate_template_for,
        register_coordinate_template,
        unregister_coordinate_template,
    )

    source = tmp_path / "topology.cif"
    _write_minimal_cif(source)
    prepared = prepare_input_structure(source)
    xyz = tmp_path / "trajectory.xyz"
    xyz.write_text(
        "2\nvalid\nC 0 0 0\nO 1 0 0\n" + bad_frame,
        encoding="utf-8",
    )
    out_pdb = tmp_path / "result.pdb"
    out_cif = tmp_path / "result.cif"
    pdb_before = b"existing pdb generation\n"
    cif_before = b"existing cif generation\n"
    out_pdb.write_bytes(pdb_before)
    out_cif.write_bytes(cif_before)
    assert prepared.structure_template is not None
    previous_template = replace(
        prepared.structure_template,
        reason="pre-existing output generation",
    )
    register_coordinate_template(out_pdb, previous_template)

    try:
        with pytest.raises(ValueError, match=message):
            convert_xyz_to_pdb(xyz, prepared.source_path, out_pdb)

        assert out_pdb.read_bytes() == pdb_before
        assert out_cif.read_bytes() == cif_before
        assert coordinate_template_for(out_pdb) is previous_template
    finally:
        unregister_coordinate_template(out_pdb)
        prepared.cleanup()


def test_template_free_generation_removes_prior_cif_companion(tmp_path: Path) -> None:
    from mlmm.io.structure_formats import register_output_template_and_write_cif

    out_pdb = tmp_path / "current.pdb"
    out_pdb.write_text("END\n", encoding="utf-8")
    old_companion = out_pdb.with_suffix(".cif")
    old_companion.write_text("data_stale\n", encoding="utf-8")

    assert register_output_template_and_write_cif(out_pdb, None) is None
    assert not old_companion.exists()


def test_xyz_overlay_rejects_swapped_coordinate_template_before_mutation(
    tmp_path: Path,
) -> None:
    from dataclasses import replace

    from mlmm.core.utils import convert_xyz_to_pdb, prepare_input_structure
    from mlmm.io.structure_formats import (
        coordinate_template_for,
        register_coordinate_template,
        unregister_coordinate_template,
    )

    source = tmp_path / "topology.cif"
    _write_minimal_cif(source)
    xyz = tmp_path / "frame.xyz"
    xyz.write_text("2\nframe\nC 2 0 0\nO 3 0 0\n", encoding="utf-8")
    prepared = prepare_input_structure(source)
    out_pdb = tmp_path / "existing.pdb"
    out_cif = out_pdb.with_suffix(".cif")
    out_pdb.write_bytes(b"old-pdb\n")
    out_cif.write_bytes(b"old-cif\n")
    try:
        assert prepared.structure_template is not None
        swapped_template = replace(
            prepared.structure_template,
            records=tuple(reversed(prepared.structure_template.records)),
            reason="same-count swapped template",
        )
        previous_template = replace(
            prepared.structure_template,
            reason="pre-existing output generation",
        )
        register_coordinate_template(prepared.source_path, swapped_template)
        register_coordinate_template(out_pdb, previous_template)

        with pytest.raises(
            ValueError,
            match="retained coordinate template.*atom 1",
        ):
            convert_xyz_to_pdb(xyz, prepared.source_path, out_pdb)

        assert out_pdb.read_bytes() == b"old-pdb\n"
        assert out_cif.read_bytes() == b"old-cif\n"
        assert coordinate_template_for(out_pdb) is previous_template
    finally:
        unregister_coordinate_template(out_pdb)
        prepared.cleanup()


def test_xyz_overlay_valid_trajectory_emits_every_model(tmp_path: Path) -> None:
    from mlmm.core.utils import convert_xyz_to_pdb

    ref = tmp_path / "ref.pdb"
    ref.write_text(
        "ATOM      1  C   MOL A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    xyz = tmp_path / "trajectory.xyz"
    xyz.write_text(
        "1\nfirst\nC 1 0 0\n1\nsecond\nC 2 0 0\n",
        encoding="utf-8",
    )
    out = tmp_path / "trajectory.pdb"

    convert_xyz_to_pdb(xyz, ref, out)

    content = out.read_text(encoding="utf-8")
    assert content.count("MODEL") == 2
    assert content.count("ENDMDL") == 2
    assert "   1.000" in content
    assert "   2.000" in content


def test_xyz_overlay_companion_publish_failure_preserves_primary_and_registry(
    tmp_path: Path,
    monkeypatch,
) -> None:
    from mlmm.core import result_commit
    from mlmm.core.utils import convert_xyz_to_pdb, prepare_input_structure
    from mlmm.io.structure_formats import (
        coordinate_template_for,
        register_coordinate_template,
        unregister_coordinate_template,
    )

    source = tmp_path / "topology.cif"
    _write_minimal_cif(source)
    xyz = tmp_path / "frame.xyz"
    xyz.write_text("2\nframe\nC 2 0 0\nO 3 0 0\n", encoding="utf-8")
    out_pdb = tmp_path / "existing.pdb"
    out_cif = out_pdb.with_suffix(".cif")
    out_pdb.write_bytes(b"old-pdb\n")
    out_cif.write_bytes(b"old-cif\n")

    prepared = prepare_input_structure(source)
    try:
        assert prepared.structure_template is not None
        previous_template = prepared.structure_template
        register_coordinate_template(out_pdb, previous_template)
        real_replace = result_commit._replace_exact

        def fail_companion(staged: Path, destination: Path) -> None:
            if destination == out_cif:
                raise OSError("injected companion publication failure")
            real_replace(staged, destination)

        monkeypatch.setattr(result_commit, "_replace_exact", fail_companion)
        with pytest.raises(
            result_commit.ResultCommitError,
            match="companion publication failure",
        ):
            convert_xyz_to_pdb(xyz, prepared.source_path, out_pdb)

        assert out_pdb.read_bytes() == b"old-pdb\n"
        assert out_cif.read_bytes() == b"old-cif\n"
        assert coordinate_template_for(out_pdb) is previous_template
    finally:
        unregister_coordinate_template(out_pdb)
        prepared.cleanup()


def test_xyz_overlay_primary_publish_failure_rolls_back_companion_and_registry(
    tmp_path: Path,
    monkeypatch,
) -> None:
    from mlmm.core import result_commit
    from mlmm.core.utils import convert_xyz_to_pdb, prepare_input_structure
    from mlmm.io.structure_formats import (
        coordinate_template_for,
        register_coordinate_template,
        unregister_coordinate_template,
    )

    source = tmp_path / "topology.cif"
    _write_minimal_cif(source)
    xyz = tmp_path / "frame.xyz"
    xyz.write_text("2\nframe\nC 2 0 0\nO 3 0 0\n", encoding="utf-8")
    out_pdb = tmp_path / "existing.pdb"
    out_cif = out_pdb.with_suffix(".cif")
    out_pdb.write_bytes(b"old-pdb\n")
    out_cif.write_bytes(b"old-cif\n")

    prepared = prepare_input_structure(source)
    try:
        assert prepared.structure_template is not None
        previous_template = prepared.structure_template
        register_coordinate_template(out_pdb, previous_template)
        real_replace = result_commit._replace_exact
        failed = False

        def fail_primary_once(staged: Path, destination: Path) -> None:
            nonlocal failed
            if destination == out_pdb and not failed:
                failed = True
                raise OSError("injected primary publication failure")
            real_replace(staged, destination)

        monkeypatch.setattr(result_commit, "_replace_exact", fail_primary_once)
        with pytest.raises(
            result_commit.ResultCommitError,
            match="primary publication failure",
        ):
            convert_xyz_to_pdb(xyz, prepared.source_path, out_pdb)

        assert out_pdb.read_bytes() == b"old-pdb\n"
        assert out_cif.read_bytes() == b"old-cif\n"
        assert coordinate_template_for(out_pdb) is previous_template
    finally:
        unregister_coordinate_template(out_pdb)
        prepared.cleanup()


def test_pdb_with_more_than_ten_thousand_residues_uses_safe_bridge(tmp_path: Path) -> None:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from mlmm.core.utils import load_pdb_atom_metadata, prepare_input_structure
    from mlmm.io.structure_formats import write_pdb_as_mmcif
    from pysisyphus.io.pdb import parse_pdb

    source = tmp_path / "ten-thousand.pdb"
    lines = []
    for index in range(10_001):
        chain = "A"
        resseq = index + 1
        serial = index + 1
        lines.append(
            f"ATOM  {serial:5d}  CA  GLY {chain}{resseq:4d}    "
            f"{float(index % 10):8.3f}{0.0:8.3f}{0.0:8.3f}"
            f"{1.0:6.2f}{0.0:6.2f}           C\n"
        )
    lines.append("END\n")
    source.write_text("".join(lines), encoding="utf-8")

    prepared = prepare_input_structure(source)
    try:
        assert prepared.source_path != source
        assert prepared.structure_template is not None
        assert prepared.structure_template.natoms == 10_001
        atoms, _coords, fragments, _atom_map = parse_pdb(str(prepared.source_path))
        assert len(atoms) == 10_001
        assert len(fragments) == 10_001
        metadata = load_pdb_atom_metadata(prepared.source_path)
        assert metadata[0]["chain"] == "A"
        assert metadata[-1]["chain"] == "A"
        assert metadata[-1]["resseq"] == 10_001
        out_cif = tmp_path / "large-output.cif"
        write_pdb_as_mmcif(
            prepared.source_path, prepared.structure_template, out_cif
        )
        data = MMCIF2Dict(str(out_cif))
        assert data["_atom_site.auth_asym_id"][-1] == "A"
        assert data["_atom_site.auth_seq_id"][-1] == "10001"
    finally:
        prepared.cleanup()


def test_bundled_pdb_parser_distinguishes_two_letter_atoms_from_hydrogens(
    tmp_path: Path,
) -> None:
    from pysisyphus.io.pdb import parse_pdb

    records = [
        ("HG  ", "HG"),
        ("HE  ", "HE"),
        (" HG ", "H"),
        (" HE ", "H"),
        (" NH1", "NH"),
        ("ZN  ", "N"),
    ]
    lines = []
    for serial, (name, element) in enumerate(records, start=1):
        lines.append(
            f"HETATM{serial:5d} {name:4s} RES A{serial:4d}    "
            f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}"
            f"{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
        )
    path = tmp_path / "elements.pdb"
    path.write_text("".join(lines) + "END\n", encoding="utf-8")

    atoms, *_ = parse_pdb(str(path))
    assert atoms == ["Hg", "He", "H", "H", "N", "Zn"]


def test_chain_qualified_atom_selector_disambiguates_repeated_ids() -> None:
    from mlmm.core.utils import resolve_atom_spec_index

    metadata = [
        {"chain": "A", "resname": "SAM", "resseq": 12, "name": "C1"},
        {"chain": "B", "resname": "SAM", "resseq": 12, "name": "C1"},
    ]
    assert resolve_atom_spec_index("B:SAM:12:C1", metadata) == 1


def test_chain_qualified_atom_selector_keeps_duplicate_field_roles_distinct() -> None:
    from mlmm.core.utils import resolve_atom_spec_index

    metadata = [
        {"chain": "A", "resname": "SAM", "resseq": 12, "name": "A"},
        {"chain": "B", "resname": "SAM", "resseq": 12, "name": "A"},
    ]
    assert resolve_atom_spec_index("A:SAM:12:A", metadata) == 0


def test_chain_qualified_atom_selector_accepts_insertion_code() -> None:
    from mlmm.core.utils import resolve_atom_spec_index

    metadata = [
        {"chain": "A", "resname": "SAM", "resseq": 12, "icode": "A", "name": "C1"},
        {"chain": "A", "resname": "SAM", "resseq": 12, "icode": "B", "name": "C1"},
    ]
    assert resolve_atom_spec_index("A:SAM:12B:C1", metadata) == 1


def test_chain_qualified_atom_selector_keeps_chain_case_distinct() -> None:
    from mlmm.core.utils import resolve_atom_spec_index

    metadata = [
        {"chain": "A", "resname": "SAM", "resseq": 12, "name": "C1"},
        {"chain": "a", "resname": "SAM", "resseq": 12, "name": "C1"},
    ]
    assert resolve_atom_spec_index("a:SAM:12:C1", metadata) == 1


def test_atom_label_includes_chain_for_repeated_residues() -> None:
    from mlmm.core.utils import atom_label_from_meta

    metadata = [
        {"chain": "LONG_CHAIN", "resname": "SAM", "resseq": 10001, "name": "C1"}
    ]
    assert atom_label_from_meta(metadata, 0) == "LONG_CHAIN:SAM:10001:C1"


def test_atom_label_includes_insertion_code() -> None:
    from mlmm.core.utils import atom_label_from_meta

    metadata = [
        {"chain": "A", "resname": "SAM", "resseq": 12, "icode": "B", "name": "C1"}
    ]
    assert atom_label_from_meta(metadata, 0) == "A:SAM:12B:C1"


def test_hybrid36_upper_and_lower_ranges_are_decoded() -> None:
    from mlmm.io.structure_formats import _hy36decode

    assert _hy36decode(4, "A000") == 10_000
    assert _hy36decode(4, "ZZZZ") == 1_223_055
    assert _hy36decode(4, "a000") == 1_223_056
    assert _hy36decode(5, "A0000") == 100_000


def test_internal_atom_serial_wraps_without_six_digit_pdb_field() -> None:
    from mlmm.io.structure_formats import _internal_atom_serial

    assert _internal_atom_serial(0) == 1
    assert _internal_atom_serial(99_998) == 99_999
    assert _internal_atom_serial(99_999) == 1


def test_decimal_overflow_pdb_fields_are_normalized(tmp_path: Path) -> None:
    from mlmm.io.structure_formats import read_pdb_atom_sites

    source = tmp_path / "overflow.pdb"
    source.write_text(
        f"ATOM  {100000:5d}  CA  GLY A{100000:4d}    "
        f"{1.25:8.3f}{2.5:8.3f}{3.75:8.3f}{1.0:6.2f}{4.0:6.2f}           C1+\nEND\n",
        encoding="utf-8",
    )

    records, nonstandard = read_pdb_atom_sites(source)
    assert nonstandard
    assert len(records) == 1
    assert records[0].chain_id == "A"
    assert records[0].resseq == "100000"
    assert records[0].formal_charge == "1"
    assert np.allclose([records[0].x, records[0].y, records[0].z], [1.25, 2.5, 3.75])


@pytest.mark.parametrize("coordinate", ["nan", "inf", "-inf"])
def test_pdb_reader_rejects_nonfinite_coordinates(
    tmp_path: Path, coordinate: str
) -> None:
    from mlmm.io.structure_formats import read_pdb_atom_sites

    source = tmp_path / "nonfinite.pdb"
    source.write_text(
        "ATOM      1  C   MOL A   1    "
        f"{coordinate:>8}{0.0:8.3f}{0.0:8.3f}"
        f"{1.0:6.2f}{0.0:6.2f}           C\nEND\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Non-finite coordinates"):
        read_pdb_atom_sites(source)


@pytest.mark.parametrize("bfactor", [1000.0, -100.0, np.nan, np.inf])
def test_internal_pdb_rejects_bfactor_field_overflow(
    tmp_path: Path, bfactor: float
) -> None:
    from dataclasses import replace

    from mlmm.io.structure_formats import (
        read_mmcif_atom_sites,
        write_internal_pdb,
    )

    source = tmp_path / "source.cif"
    _write_minimal_cif(source)
    records = read_mmcif_atom_sites(source)
    records[0] = replace(records[0], bfactor=bfactor)

    with pytest.raises(ValueError, match="six-column PDB range"):
        write_internal_pdb(records, tmp_path / "internal.pdb")


def test_duplicate_atom_names_without_altloc_are_preserved(tmp_path: Path) -> None:
    from mlmm.io.structure_formats import (
        pdb_requires_normalization,
        read_pdb_atom_sites,
    )

    source = tmp_path / "small-molecule.pdb"
    source.write_text(
        "ATOM      1  C   MOL     1       0.000   0.000   0.000                       C\n"
        "ATOM      2  H   MOL     1       1.000   0.000   0.000                       H\n"
        "ATOM      3  H   MOL     1      -1.000   0.000   0.000                       H\n"
        "END\n",
        encoding="utf-8",
    )

    records, nonstandard = read_pdb_atom_sites(source)
    assert [record.atom_name for record in records] == ["C", "H", "H"]
    assert not nonstandard
    assert not pdb_requires_normalization(source)


def test_failed_internal_pdb_write_removes_temporary_bridge(
    tmp_path: Path, monkeypatch
) -> None:
    from mlmm.io import structure_formats

    source = tmp_path / "out-of-range.cif"
    _write_minimal_cif(source)
    source.write_text(
        source.read_text(encoding="utf-8").replace("0.0 1.0 2.0", "10000.0 1.0 2.0"),
        encoding="utf-8",
    )
    bridge_dir = tmp_path / "bridge"

    def controlled_mkdtemp(*, prefix):
        assert prefix == "mlmm_structure_"
        bridge_dir.mkdir()
        return str(bridge_dir)

    monkeypatch.setattr(structure_formats.tempfile, "mkdtemp", controlled_mkdtemp)
    with np.testing.assert_raises_regex(ValueError, "fixed-column PDB range"):
        structure_formats.normalize_structure_to_pdb(source)
    assert not bridge_dir.exists()


def test_ref_cif_atom_count_error_cleans_temporary_bridge(
    tmp_path: Path, monkeypatch
) -> None:
    from click import ClickException
    from mlmm.core import utils

    xyz = tmp_path / "one.xyz"
    xyz.write_text("1\nframe\nC 0 0 0\n", encoding="utf-8")
    ref = tmp_path / "two.cif"
    _write_minimal_cif(ref)
    bridge_dirs = []
    normalize = utils.normalize_structure_to_pdb

    def record_bridge(path):
        result = normalize(path)
        bridge_dirs.append(result[2])
        return result

    monkeypatch.setattr(utils, "normalize_structure_to_pdb", record_bridge)
    prepared = utils.prepare_input_structure(xyz)
    try:
        with pytest.raises(ClickException, match="atom count"):
            utils.apply_ref_pdb_override(prepared, ref)
    finally:
        prepared.cleanup()

    assert bridge_dirs
    assert all(not path.exists() for path in bridge_dirs)


def test_ref_pdb_override_rejects_element_order_mismatch(tmp_path: Path) -> None:
    from click import BadParameter
    from mlmm.core import utils

    xyz = tmp_path / "input.xyz"
    xyz.write_text(
        "2\ninput\nC 0 0 0\nH 1 0 0\n",
        encoding="utf-8",
    )
    ref = tmp_path / "reference.pdb"
    ref.write_text(
        "HETATM    1  H1  MOL A   1       0.000   0.000   0.000"
        "  1.00  0.00           H\n"
        "HETATM    2  C1  MOL A   1       1.000   0.000   0.000"
        "  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    prepared = utils.prepare_input_structure(xyz)
    try:
        with pytest.raises(BadParameter, match="atom-order element mismatch"):
            utils.apply_ref_pdb_override(prepared, ref)
    finally:
        prepared.cleanup()


def test_dft_dry_run_uses_one_reference_overlaid_preparation_with_explicit_charge(
    tmp_path: Path,
    monkeypatch,
) -> None:
    from types import MethodType, SimpleNamespace

    from ase import Atoms
    from click.testing import CliRunner
    from mlmm.cli import cli as root_cli
    from mlmm.workflows import dft

    xyz = tmp_path / "input.xyz"
    xyz.write_text("1\ninput\nC 1.0 2.0 3.0\n", encoding="utf-8")
    ref = tmp_path / "reference.pdb"
    ref.write_text(
        "HETATM    7  C1  SAM A   8       1.000   2.000   3.000  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    parm = tmp_path / "dummy.parm7"
    parm.write_text("dry-run placeholder\n", encoding="utf-8")

    real_prepare = dft.prepare_input_structure
    prepared_objects = []
    cleanup_calls = []
    workspace_calls = []
    workspace_cleanup_calls = []

    def recording_prepare(path):
        prepared = real_prepare(path)
        prepared_objects.append(prepared)
        original_cleanup = prepared.cleanup

        def recording_cleanup(self):
            cleanup_calls.append(self)
            original_cleanup()

        prepared.cleanup = MethodType(recording_cleanup, prepared)
        return prepared

    monkeypatch.setattr(dft, "prepare_input_structure", recording_prepare)

    def fake_workspace(**kwargs):
        workspace_calls.append(kwargs)
        return SimpleNamespace(
            atoms_model_lh=Atoms("C"),
            cleanup=lambda: workspace_cleanup_calls.append(True),
        )

    monkeypatch.setattr(dft, "_prepare_ml_region_workspace", fake_workspace)

    result = CliRunner().invoke(
        root_cli,
        [
            "dft",
            "-i",
            str(xyz),
            "--ref-pdb",
            str(ref),
            "--parm",
            str(parm),
            "-q",
            "-1",
            "-m",
            "2",
            "--model-indices",
            "1",
            "--dry-run",
            "--out-dir",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert result.output.rstrip().splitlines()[-1] == (
        "[Dry run] --dry-run completed. Input command is valid."
    )
    assert len(prepared_objects) == 1
    assert cleanup_calls == prepared_objects
    assert len(workspace_calls) == 1
    assert workspace_calls[0]["input_pdb"].resolve() == ref.resolve()
    assert workspace_calls[0]["coordinate_path"].resolve() == xyz.resolve()
    assert workspace_cleanup_calls == [True]


def test_dft_dry_run_rejects_malformed_amber_topology(tmp_path: Path) -> None:
    from click.testing import CliRunner
    from mlmm.cli import cli as root_cli

    pdb = tmp_path / "input.pdb"
    pdb.write_text(
        "HETATM    1  C1  MOL A   1       1.000   2.000   3.000"
        "  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    parm = tmp_path / "invalid.parm7"
    parm.write_text("not an Amber topology\n", encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            "dft",
            "-i",
            str(pdb),
            "--parm",
            str(parm),
            "-q",
            "0",
            "-m",
            "1",
            "--model-indices",
            "1",
            "--dry-run",
            "--out-dir",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code != 0
    assert "is not a valid Amber parm7 file" in result.output
    assert "[Dry run] --dry-run completed." not in result.output


def test_dft_workspace_keeps_xyz_coordinates_with_pdb_topology(
    tmp_path: Path,
) -> None:
    from ase import Atoms
    from ase.io import write as ase_write
    from mlmm.io.pdb_indexing import parse_pdb_ordinal_atoms
    from mlmm.workflows.dft import _prepare_ml_region_workspace

    repo = Path(__file__).resolve().parents[1]
    input_pdb = repo / "hessian_ff/tests/data/small/p_complex_layered.pdb"
    parm7 = repo / "hessian_ff/tests/data/small/p_complex.parm7"
    records = parse_pdb_ordinal_atoms(input_pdb)
    atoms = Atoms(
        symbols=[record.elem for record in records],
        positions=[record.coord for record in records],
    )
    positions = atoms.get_positions()
    positions[0] += np.array([0.00037, -0.00041, 0.00029])
    atoms.set_positions(positions)
    xyz = tmp_path / "precise.xyz"
    ase_write(str(xyz), atoms, format="xyz")

    workspace = _prepare_ml_region_workspace(
        input_pdb=input_pdb,
        coordinate_path=xyz,
        real_parm7=parm7,
        model_pdb=input_pdb,
        link_mlmm=None,
        calc_kwargs={"model_charge": 0, "model_mult": 2},
    )
    try:
        np.testing.assert_allclose(
            workspace.atoms_real.get_positions(),
            atoms.get_positions(),
            atol=1.0e-8,
        )
        np.testing.assert_allclose(
            workspace.atoms_model.get_positions(),
            atoms.get_positions(),
            atol=1.0e-8,
        )
    finally:
        workspace.cleanup()


def test_pdb_to_cif_preserves_output_occupancy_and_bfactor(tmp_path: Path) -> None:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from mlmm.io.structure_formats import (
        AtomSiteRecord,
        CoordinateTemplate,
        write_pdb_as_mmcif,
    )

    pdb = tmp_path / "marked.pdb"
    pdb.write_text(
        "HETATM    1  C1  SAM A   1       1.000   2.000   3.000  0.50 99.00           C\nEND\n",
        encoding="utf-8",
    )
    record = AtomSiteRecord(
        group_pdb="HETATM",
        element="C",
        atom_name="C1",
        altloc="",
        resname="SAM",
        chain_id="LONG_CHAIN",
        resseq="10001",
        icode="",
        occupancy=1.0,
        bfactor=12.0,
    )
    template = CoordinateTemplate((record,), pdb, "mmcif", "test")
    out = tmp_path / "marked.cif"
    write_pdb_as_mmcif(pdb, template, out)

    data = MMCIF2Dict(str(out))
    assert data["_atom_site.occupancy"] == ["0.50"]
    assert data["_atom_site.B_iso_or_equiv"] == ["99.00"]


def test_multimodel_pdb_is_written_as_multimodel_cif(tmp_path: Path) -> None:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from mlmm.io.structure_formats import (
        AtomSiteRecord,
        CoordinateTemplate,
        write_pdb_as_mmcif,
    )

    pdb = tmp_path / "trajectory.pdb"
    pdb.write_text(
        "MODEL        1\n"
        "HETATM    1  C1  SAM A   1       1.000   2.000   3.000  1.00 10.00           C\n"
        "ENDMDL\n"
        "MODEL        2\n"
        "HETATM    1  C1  SAM A   1       4.000   5.000   6.000  1.00 20.00           C\n"
        "ENDMDL\nEND\n",
        encoding="utf-8",
    )
    record = AtomSiteRecord(
        group_pdb="HETATM",
        element="C",
        atom_name="C1",
        altloc="",
        resname="SAM",
        chain_id="LONG_CHAIN",
        resseq="10001",
        icode="",
        occupancy=1.0,
        bfactor=0.0,
    )
    template = CoordinateTemplate((record,), pdb, "mmcif", "test")
    out = tmp_path / "trajectory.cif"
    write_pdb_as_mmcif(pdb, template, out)

    data = MMCIF2Dict(str(out))
    assert data["_atom_site.pdbx_PDB_model_num"] == ["1", "2"]
    assert data["_atom_site.Cartn_x"] == ["1.000000", "4.000000"]
    assert data["_atom_site.B_iso_or_equiv"] == ["10.00", "20.00"]


def test_define_layer_model_pdb_matches_one_unused_occurrence(tmp_path: Path) -> None:
    from mlmm.workflows.define_layer import (
        _get_ml_indices_from_model_pdb,
        _parse_pdb_atoms,
    )

    atom = (
        "HETATM    1  C1  SAM A   1       1.000   2.000   3.000"
        "  1.00 10.00           C\n"
    )
    full = tmp_path / "full.pdb"
    model = tmp_path / "model.pdb"
    full.write_text(atom + atom.replace(" 1.000", " 4.000", 1) + "END\n")
    model.write_text(atom + "END\n")

    assert _get_ml_indices_from_model_pdb(
        _parse_pdb_atoms(full), model, full
    ) == [0]


def test_define_layer_blank_chain_rejects_cross_chain_ambiguity(tmp_path: Path) -> None:
    from mlmm.workflows.define_layer import (
        _get_ml_indices_from_model_pdb,
        _parse_pdb_atoms,
    )

    atom_a = (
        "HETATM    1  C1  SAM A   1       1.000   2.000   3.000"
        "  1.00 10.00           C\n"
    )
    atom_b = atom_a.replace("SAM A", "SAM B").replace("    1", "    2", 1)
    full = tmp_path / "full.pdb"
    model = tmp_path / "model.pdb"
    full.write_text(atom_a + atom_b + "END\n")
    model.write_text(atom_a.replace("SAM A", "SAM  ") + "END\n")

    with pytest.raises(ValueError, match="Blank-chain.*ambiguous"):
        _get_ml_indices_from_model_pdb(_parse_pdb_atoms(full), model, full)


def test_define_layer_multimodel_warns_and_writes_only_first_model(
    tmp_path: Path, capsys
) -> None:
    from mlmm.workflows.define_layer import _define_layers_pdb

    source = tmp_path / "trajectory.pdb"
    source.write_text(
        "MODEL        1\n"
        "HETATM    1  C1  SAM A   1       1.000   2.000   3.000  1.00 10.00           C\n"
        "ENDMDL\n"
        "MODEL        2\n"
        "HETATM    1  C1  SAM A   1       4.000   5.000   6.000  1.00 20.00           C\n"
        "ENDMDL\nEND\n",
        encoding="utf-8",
    )
    output = tmp_path / "layered.pdb"

    _define_layers_pdb(source, output, model_indices=[0])

    captured = capsys.readouterr()
    text = output.read_text(encoding="utf-8")
    assert "using first model and ignoring the rest" in captured.err
    assert text.count("MODEL") == 1
    assert "   1.000" in text
    assert "   4.000" not in text


def test_cif_altloc_selection_is_coherent_per_residue(tmp_path: Path) -> None:
    from mlmm.io.structure_formats import read_mmcif_atom_sites

    source = tmp_path / "altloc.cif"
    source.write_text(
        """data_alt
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 1 ? 0 0 0 1.0 0 1 ALA A CA 1
ATOM 2 C CB A ALA A 1 1 ? 1 0 0 0.9 0 1 ALA A CB 1
ATOM 3 C CG A ALA A 1 1 ? 2 0 0 0.1 0 1 ALA A CG 1
ATOM 4 C CB B ALA A 1 1 ? 3 0 0 0.2 0 1 ALA A CB 1
ATOM 5 C CD B ALA A 1 1 ? 4 0 0 0.8 0 1 ALA A CD 1
#
""",
        encoding="utf-8",
    )

    records = read_mmcif_atom_sites(source)
    assert [record.atom_name for record in records] == ["CA", "CB", "CG"]
    assert all(record.altloc == "" for record in records)


def test_cif_altloc_parsed_zero_beats_missing_occupancy(tmp_path: Path) -> None:
    from mlmm.io.structure_formats import read_mmcif_atom_sites

    source = tmp_path / "missing-occupancy-altloc.cif"
    source.write_text(
        """data_alt
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA A ALA A 1 1 ? 1 0 0 0.0 0 1 ALA A CA 1
ATOM 2 C CA B ALA A 1 1 ? 2 0 0 ? 0 1 ALA A CA 1
#
""",
        encoding="utf-8",
    )

    records = read_mmcif_atom_sites(source)

    assert len(records) == 1
    assert records[0].x == pytest.approx(1.0)
    assert records[0].occupancy == pytest.approx(0.0)
    assert records[0].occupancy_known


def test_altloc_selection_preserves_repeated_blank_atom_names(tmp_path: Path) -> None:
    from mlmm.io.structure_formats import read_mmcif_atom_sites

    source = tmp_path / "blank-duplicates.cif"
    source.write_text(
        """data_alt
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
HETATM 1 H H . MOL A 1 1 ? 0 0 0 1.0 0 1 MOL A H 1
HETATM 2 H H . MOL A 1 1 ? 1 0 0 1.0 0 1 MOL A H 1
HETATM 3 C C1 A MOL A 1 1 ? 2 0 0 0.8 0 1 MOL A C1 1
HETATM 4 C C1 B MOL A 1 1 ? 3 0 0 0.2 0 1 MOL A C1 1
#
""",
        encoding="utf-8",
    )

    records = read_mmcif_atom_sites(source)
    assert [record.atom_name for record in records] == ["H", "H", "C1"]


def test_cif_duplicate_atom_names_survive_internal_pdb_and_structure_load(
    tmp_path: Path,
) -> None:
    from mlmm.core.utils import prepare_input_structure
    from mlmm.workflows.extract import load_structure

    source = tmp_path / "duplicate-names.cif"
    source.write_text(
        """data_duplicate
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
HETATM 1 H H . MOL A 1 1 ? 0 0 0 1.0 0 1 MOL A H 1
HETATM 2 H H . MOL A 1 1 ? 1 0 0 1.0 0 1 MOL A H 1
HETATM 3 C C1 . MOL A 1 1 ? 2 0 0 1.0 0 1 MOL A C1 1
#
""",
        encoding="utf-8",
    )

    prepared = prepare_input_structure(source)
    try:
        structure = load_structure(str(prepared.source_path), "duplicate")
    finally:
        prepared.cleanup()

    atoms = list(structure.get_atoms())
    assert len(atoms) == 3
    assert len({atom.get_name() for atom in atoms}) == 3
    assert [atom.xtra["mlmm_atom_site"].atom_name for atom in atoms] == ["H", "H", "C1"]


def test_all_segment_copy_emits_registered_cif_companion(tmp_path: Path) -> None:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from mlmm.io.structure_formats import (
        AtomSiteRecord,
        CoordinateTemplate,
        register_coordinate_template,
        unregister_coordinate_template,
    )
    from mlmm.workflows.all import _copy_structures_to_seg_dir

    xyz = tmp_path / "reactant.xyz"
    xyz.write_text("1\nreactant\nC 1 2 3\n", encoding="utf-8")
    pdb = tmp_path / "reactant.pdb"
    pdb.write_text(
        "HETATM    1  C1  MOL A   1       1.000   2.000   3.000  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    record = AtomSiteRecord(
        group_pdb="HETATM",
        element="C",
        atom_name="C1",
        altloc="",
        resname="MOL",
        chain_id="LONG_CHAIN",
        resseq="10001",
        icode="",
        occupancy=1.0,
        bfactor=0.0,
        x=1.0,
        y=2.0,
        z=3.0,
    )
    template = CoordinateTemplate((record,), tmp_path / "input.cif", "mmcif", "test")
    register_coordinate_template(pdb, template)
    try:
        seg_dir = _copy_structures_to_seg_dir(
            {"R": xyz}, tmp_path / "out", 1, ".cif"
        )
        out_pdb = seg_dir / "reactant.pdb"
        out_cif = seg_dir / "reactant.cif"
        assert out_pdb.exists()
        assert out_cif.exists()
        data = MMCIF2Dict(str(out_cif))
        assert data["_atom_site.auth_asym_id"] == ["LONG_CHAIN"]
        assert data["_atom_site.auth_seq_id"] == ["10001"]
    finally:
        unregister_coordinate_template(pdb)
        unregister_coordinate_template(tmp_path / "out" / "segments" / "seg_01" / "reactant.pdb")


def test_extract_cif_accepts_chain_resname_and_emits_cif(tmp_path: Path) -> None:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from mlmm.workflows.extract import extract_api

    source = tmp_path / "complex.cif"
    _write_minimal_cif(source)
    out_pdb = tmp_path / "new" / "nested" / "model.pdb"
    result = extract_api(
        complex_pdb=[str(source)],
        center="LONG_CHAIN:SAM",
        output=[str(out_pdb)],
        radius=0.1,
        include_h2o=False,
        add_linkh=False,
    )

    out_cif = out_pdb.with_suffix(".cif")
    assert out_pdb.exists()
    assert out_cif.exists()
    assert str(out_cif) in result["outputs"]
    data = MMCIF2Dict(str(out_cif))
    assert data["_atom_site.auth_asym_id"] == ["LONG_CHAIN", "LONG_CHAIN"]
    assert data["_atom_site.auth_seq_id"] == ["10001", "10001"]


def test_extract_multi_input_creates_nested_output_directories(tmp_path: Path) -> None:
    from mlmm.workflows.extract import extract_api

    source = tmp_path / "complex.cif"
    _write_minimal_cif(source)
    first = tmp_path / "per_file" / "a" / "model_a.pdb"
    second = tmp_path / "per_file" / "b" / "model_b.pdb"
    extract_api(
        complex_pdb=[str(source), str(source)],
        center="LONG_CHAIN:SAM",
        output=[str(first), str(second)],
        radius=0.1,
        include_h2o=False,
        add_linkh=False,
    )
    assert first.exists()
    assert second.exists()

    combined = tmp_path / "combined" / "nested" / "models.pdb"
    extract_api(
        complex_pdb=[str(source), str(source)],
        center="LONG_CHAIN:SAM",
        output=[str(combined)],
        radius=0.1,
        include_h2o=False,
        add_linkh=False,
    )
    assert combined.exists()
    assert combined.read_text().count("MODEL") == 2


def test_define_layer_matches_cif_author_identifiers_and_emits_cif(
    tmp_path: Path,
) -> None:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from mlmm.workflows.define_layer import define_layers

    full = tmp_path / "full.cif"
    model = tmp_path / "model.cif"
    _write_minimal_cif(full)
    model.write_text(
        full.read_text(encoding="utf-8").replace(
            "HETATM 2 O O1 . SAM LONG_CHAIN 1 . ? 1.0 1.0 2.0 1.00 13.0 . 10001 SAM LONG_CHAIN O1 1\n",
            "",
        ),
        encoding="utf-8",
    )
    out = tmp_path / "layered.pdb"

    layers = define_layers(full, out, model_pdb=model)

    assert layers["ml_indices"] == [0]
    data = MMCIF2Dict(str(out.with_suffix(".cif")))
    assert data["_atom_site.auth_asym_id"] == ["LONG_CHAIN", "LONG_CHAIN"]
    assert data["_atom_site.auth_seq_id"] == ["10001", "10001"]
    assert data["_atom_site.B_iso_or_equiv"] == ["0.00", "10.00"]


def test_define_layer_accepts_extracted_internal_model_pdb(tmp_path: Path) -> None:
    from mlmm.core.utils import prepare_input_structure
    from mlmm.workflows.define_layer import define_layers

    full = tmp_path / "full.cif"
    _write_minimal_cif(full)
    prepared = prepare_input_structure(full)
    try:
        first_atom = next(
            line
            for line in prepared.source_path.read_text(encoding="utf-8").splitlines()
            if line.startswith(("ATOM", "HETATM"))
        )
        model = tmp_path / "extracted-model.pdb"
        model.write_text(first_atom + "\nEND\n", encoding="utf-8")
        layers = define_layers(full, tmp_path / "layered.pdb", model_pdb=model)
    finally:
        prepared.cleanup()

    assert layers["ml_indices"] == [0]


def test_define_layer_rejects_partial_model_identity_match(tmp_path: Path) -> None:
    from mlmm.workflows.define_layer import define_layers

    full = tmp_path / "full.pdb"
    full.write_text(
        "HETATM    1  C1  SAM A   7       0.000   1.000   2.000  1.00 12.00           C\n"
        "HETATM    2  O1  SAM A   7       1.000   1.000   2.000  1.00 13.00           O\n"
        "END\n",
        encoding="utf-8",
    )
    model = tmp_path / "model.pdb"
    model.write_text(
        "HETATM    1  C1  SAM A   7       0.000   1.000   2.000  1.00 12.00           C\n"
        "HETATM    2  N1  SAM A   7       1.000   1.000   2.000  1.00 13.00           N\n"
        "END\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="absent from the full input"):
        define_layers(full, tmp_path / "layered.pdb", model_pdb=model)


@pytest.mark.parametrize("normalized_side", ["input", "model"])
def test_define_layer_matches_author_ids_when_only_one_side_is_normalized(
    tmp_path: Path, normalized_side: str,
) -> None:
    from mlmm.workflows.define_layer import define_layers

    cif = tmp_path / "author.cif"
    _write_minimal_cif(cif)
    cif.write_text(
        cif.read_text(encoding="utf-8")
        .replace("LONG_CHAIN", "A")
        .replace("10001", "7"),
        encoding="utf-8",
    )
    pdb = tmp_path / "author.pdb"
    pdb.write_text(
        "HETATM    1  C1  SAM A   7       0.000   1.000   2.000  1.00 12.00           C\n"
        "HETATM    2  O1  SAM A   7       1.000   1.000   2.000  1.00 13.00           O\n"
        "END\n",
        encoding="utf-8",
    )
    if normalized_side == "input":
        full, model = cif, pdb
    else:
        full, model = pdb, cif

    layers = define_layers(
        full, tmp_path / f"{normalized_side}.pdb", model_pdb=model
    )

    assert layers["ml_indices"] == [0, 1]


def test_define_layer_cli_reports_actual_pdb_for_plain_input_cif_request(
    tmp_path: Path,
) -> None:
    from click.testing import CliRunner

    from mlmm.cli import cli as root_cli

    full = tmp_path / "full.pdb"
    full.write_text(
        "HETATM    1  C1  SAM A   7       0.000   1.000   2.000  1.00 12.00           C\n"
        "HETATM    2  O1  SAM A   7       1.000   1.000   2.000  1.00 13.00           O\n"
        "END\n",
        encoding="utf-8",
    )
    requested = tmp_path / "layered.cif"

    result = CliRunner().invoke(
        root_cli,
        [
            "define-layer", "-i", str(full), "--model-indices", "1",
            "-o", str(requested),
        ],
    )

    assert result.exit_code == 0, result.output
    assert requested.with_suffix(".pdb").exists()
    assert not requested.exists()
    assert f"Wrote '{requested.with_suffix('.pdb')}'" in result.output
    assert f"Wrote '{requested}'" not in result.output


def test_parm_topology_order_validator_rejects_count_and_element_mismatch() -> None:
    from types import SimpleNamespace

    from mlmm.backends.mlmm_calc import validate_parmed_atom_order

    def structure(*atomic_numbers, names=None):
        atom_names = names or [f"A{i}" for i in range(len(atomic_numbers))]
        return SimpleNamespace(
            atoms=[
                SimpleNamespace(
                    atomic_number=value,
                    name=atom_names[i],
                    residue=SimpleNamespace(name="LIG", idx=0),
                )
                for i, value in enumerate(atomic_numbers)
            ]
        )

    with pytest.raises(ValueError, match="Atom-count mismatch"):
        validate_parmed_atom_order(structure(6), structure(6, 1))
    with pytest.raises(ValueError, match="Atom-order mismatch.*atom 2"):
        validate_parmed_atom_order(structure(6, 8), structure(6, 7))
    with pytest.raises(ValueError, match="Atom-order mismatch.*atom 1"):
        validate_parmed_atom_order(
            structure(6, 6, names=["C1", "C2"]),
            structure(6, 6, names=["C2", "C1"]),
        )
    validate_parmed_atom_order(
        structure(1, 1, names=["1HH1", "2HH1"]),
        structure(1, 1, names=["HH11", "HH12"]),
    )
    with pytest.raises(ValueError, match="Atom-order mismatch.*atom 1"):
        validate_parmed_atom_order(
            structure(6, names=["1C"]),
            structure(6, names=["C1"]),
        )
    validate_parmed_atom_order(structure(6, 0), structure(6, 1))


def test_dft_workspace_validates_topology_before_coordinate_assignment(
    tmp_path: Path,
    monkeypatch,
) -> None:
    from mlmm.backends import mlmm_calc
    from mlmm.workflows import dft

    repo = Path(__file__).resolve().parents[1]
    input_pdb = repo / "hessian_ff/tests/data/small/p_complex_layered.pdb"
    parm7 = repo / "hessian_ff/tests/data/small/p_complex.parm7"

    def reject(*args, **kwargs):
        raise ValueError("identity sentinel")

    monkeypatch.setattr(mlmm_calc, "validate_parmed_atom_order", reject)
    with pytest.raises(ValueError, match="identity sentinel"):
        dft._prepare_ml_region_workspace(
            input_pdb=input_pdb,
            real_parm7=parm7,
            model_pdb=input_pdb,
            link_mlmm=None,
        )


def test_annotated_conversion_emits_one_final_cif(tmp_path: Path, monkeypatch) -> None:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from mlmm.core import result_commit
    from mlmm.core import utils as core_utils
    from mlmm.io.structure_formats import (
        coordinate_template_for,
        unregister_coordinate_template,
    )

    source = tmp_path / "input.cif"
    _write_minimal_cif(source)
    xyz = tmp_path / "trajectory.xyz"
    xyz.write_text(
        "2\nframe\nC 3.12345 4.0 5.0\nO 6.54321 7.0 8.0\n",
        encoding="utf-8",
    )
    prepared = core_utils.prepare_input_structure(source)
    try:
        atom_line = next(
            line
            for line in prepared.source_path.read_text(encoding="utf-8").splitlines()
            if line.startswith(("ATOM", "HETATM"))
        )
        model = tmp_path / "model.pdb"
        model.write_text(atom_line + "\nEND\n", encoding="utf-8")
        out_pdb = tmp_path / "annotated.pdb"

        published = []
        real_replace = result_commit._replace_exact

        def record_replace(staged: Path, destination: Path) -> None:
            published.append(destination)
            real_replace(staged, destination)

        monkeypatch.setattr(result_commit, "_replace_exact", record_replace)
        core_utils.convert_and_annotate_xyz_to_pdb(
            xyz,
            prepared.source_path,
            out_pdb,
            model,
            freeze_indices_0based=[1],
        )

        assert published == [out_pdb.with_suffix(".cif"), out_pdb]
        data = MMCIF2Dict(str(out_pdb.with_suffix(".cif")))
        assert data["_atom_site.B_iso_or_equiv"] == ["0.00", "20.00"]
        assert np.allclose(
            [float(value) for value in data["_atom_site.Cartn_x"]],
            [3.12345, 6.54321],
        )
        assert coordinate_template_for(out_pdb) is prepared.structure_template
    finally:
        unregister_coordinate_template(tmp_path / "annotated.pdb")
        prepared.cleanup()


def test_annotated_conversion_rejects_swapped_template_before_mutation(
    tmp_path: Path,
) -> None:
    from dataclasses import replace

    from mlmm.core import utils as core_utils
    from mlmm.io.structure_formats import (
        coordinate_template_for,
        register_coordinate_template,
        unregister_coordinate_template,
    )

    source = tmp_path / "input.cif"
    _write_minimal_cif(source)
    xyz = tmp_path / "trajectory.xyz"
    xyz.write_text("2\nframe\nC 3 4 5\nO 6 7 8\n", encoding="utf-8")
    prepared = core_utils.prepare_input_structure(source)
    out_pdb = tmp_path / "annotated.pdb"
    out_cif = out_pdb.with_suffix(".cif")
    out_pdb.write_bytes(b"old annotated pdb\n")
    out_cif.write_bytes(b"old annotated cif\n")
    try:
        atom_line = next(
            line
            for line in prepared.source_path.read_text(encoding="utf-8").splitlines()
            if line.startswith(("ATOM", "HETATM"))
        )
        model = tmp_path / "model.pdb"
        model.write_text(atom_line + "\nEND\n", encoding="utf-8")
        assert prepared.structure_template is not None
        swapped_template = replace(
            prepared.structure_template,
            records=tuple(reversed(prepared.structure_template.records)),
            reason="same-count swapped template",
        )
        previous_template = replace(
            prepared.structure_template,
            reason="pre-existing annotated output generation",
        )
        register_coordinate_template(prepared.source_path, swapped_template)
        register_coordinate_template(out_pdb, previous_template)

        with pytest.raises(
            ValueError,
            match="retained coordinate template.*atom 1",
        ):
            core_utils.convert_and_annotate_xyz_to_pdb(
                xyz,
                prepared.source_path,
                out_pdb,
                model,
                freeze_indices_0based=[1],
            )

        assert out_pdb.read_bytes() == b"old annotated pdb\n"
        assert out_cif.read_bytes() == b"old annotated cif\n"
        assert coordinate_template_for(out_pdb) is previous_template
    finally:
        unregister_coordinate_template(out_pdb)
        prepared.cleanup()


def test_annotated_conversion_companion_failure_preserves_public_generation(
    tmp_path: Path,
    monkeypatch,
) -> None:
    from mlmm.core import result_commit
    from mlmm.core import utils as core_utils
    from mlmm.io.structure_formats import (
        coordinate_template_for,
        register_coordinate_template,
        unregister_coordinate_template,
    )

    source = tmp_path / "input.cif"
    _write_minimal_cif(source)
    xyz = tmp_path / "trajectory.xyz"
    xyz.write_text("2\nframe\nC 3 4 5\nO 6 7 8\n", encoding="utf-8")
    prepared = core_utils.prepare_input_structure(source)
    out_pdb = tmp_path / "annotated.pdb"
    out_cif = out_pdb.with_suffix(".cif")
    out_pdb.write_bytes(b"old annotated pdb\n")
    out_cif.write_bytes(b"old annotated cif\n")
    try:
        atom_line = next(
            line
            for line in prepared.source_path.read_text(encoding="utf-8").splitlines()
            if line.startswith(("ATOM", "HETATM"))
        )
        model = tmp_path / "model.pdb"
        model.write_text(atom_line + "\nEND\n", encoding="utf-8")
        assert prepared.structure_template is not None
        previous_template = prepared.structure_template
        register_coordinate_template(out_pdb, previous_template)
        real_replace = result_commit._replace_exact

        def fail_companion(staged: Path, destination: Path) -> None:
            if destination == out_cif:
                raise OSError("injected annotated companion failure")
            real_replace(staged, destination)

        monkeypatch.setattr(result_commit, "_replace_exact", fail_companion)
        with pytest.raises(
            result_commit.ResultCommitError,
            match="annotated companion failure",
        ):
            core_utils.convert_and_annotate_xyz_to_pdb(
                xyz,
                prepared.source_path,
                out_pdb,
                model,
                freeze_indices_0based=[1],
            )

        assert out_pdb.read_bytes() == b"old annotated pdb\n"
        assert out_cif.read_bytes() == b"old annotated cif\n"
        assert coordinate_template_for(out_pdb) is previous_template
    finally:
        unregister_coordinate_template(out_pdb)
        prepared.cleanup()


def test_annotated_conversion_primary_failure_rolls_back_public_generation(
    tmp_path: Path,
    monkeypatch,
) -> None:
    from mlmm.core import result_commit
    from mlmm.core import utils as core_utils
    from mlmm.io.structure_formats import (
        coordinate_template_for,
        register_coordinate_template,
        unregister_coordinate_template,
    )

    source = tmp_path / "input.cif"
    _write_minimal_cif(source)
    xyz = tmp_path / "trajectory.xyz"
    xyz.write_text("2\nframe\nC 3 4 5\nO 6 7 8\n", encoding="utf-8")
    prepared = core_utils.prepare_input_structure(source)
    out_pdb = tmp_path / "annotated.pdb"
    out_cif = out_pdb.with_suffix(".cif")
    out_pdb.write_bytes(b"old annotated pdb\n")
    out_cif.write_bytes(b"old annotated cif\n")
    try:
        atom_line = next(
            line
            for line in prepared.source_path.read_text(encoding="utf-8").splitlines()
            if line.startswith(("ATOM", "HETATM"))
        )
        model = tmp_path / "model.pdb"
        model.write_text(atom_line + "\nEND\n", encoding="utf-8")
        assert prepared.structure_template is not None
        previous_template = prepared.structure_template
        register_coordinate_template(out_pdb, previous_template)
        real_replace = result_commit._replace_exact
        failed = False

        def fail_primary_once(staged: Path, destination: Path) -> None:
            nonlocal failed
            if destination == out_pdb and not failed:
                failed = True
                raise OSError("injected annotated primary failure")
            real_replace(staged, destination)

        monkeypatch.setattr(result_commit, "_replace_exact", fail_primary_once)
        with pytest.raises(
            result_commit.ResultCommitError,
            match="annotated primary failure",
        ):
            core_utils.convert_and_annotate_xyz_to_pdb(
                xyz,
                prepared.source_path,
                out_pdb,
                model,
                freeze_indices_0based=[1],
            )

        assert out_pdb.read_bytes() == b"old annotated pdb\n"
        assert out_cif.read_bytes() == b"old annotated cif\n"
        assert coordinate_template_for(out_pdb) is previous_template
    finally:
        unregister_coordinate_template(out_pdb)
        prepared.cleanup()
