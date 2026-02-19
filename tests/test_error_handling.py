"""Error-handling regression tests (invalid files, YAML, and malformed inputs)."""

from __future__ import annotations

import sys
from pathlib import Path

import click
import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm_toolkit requires Python >= 3.11",
)


@pytest.mark.parametrize("spec", ["1-a", "0", "2-1"])
def test_parse_indices_string_rejects_invalid(spec):
    from mlmm_toolkit.utils import parse_indices_string

    with pytest.raises(click.BadParameter):
        parse_indices_string(spec)


def test_load_yaml_dict_rejects_non_mapping_root(tmp_path: Path):
    from mlmm_toolkit.utils import load_yaml_dict

    p = tmp_path / "bad_root.yaml"
    p.write_text("- a\n- b\n", encoding="utf-8")

    with pytest.raises(ValueError, match="YAML root must be a mapping"):
        load_yaml_dict(p)


def test_load_yaml_dict_rejects_malformed_yaml(tmp_path: Path):
    from mlmm_toolkit.utils import load_yaml_dict

    p = tmp_path / "malformed.yaml"
    p.write_text("geom: [1,2\n", encoding="utf-8")

    with pytest.raises(Exception):
        load_yaml_dict(p)


def test_collect_single_option_values_rejects_repeated_flags():
    from mlmm_toolkit.utils import collect_single_option_values

    argv = ["-i", "a.pdb", "-i", "b.pdb"]
    with pytest.raises(click.BadParameter, match="single -i/--input"):
        collect_single_option_values(argv, ("-i", "--input"), label="-i/--input")


def test_load_structure_rejects_missing_input_file(tmp_path: Path):
    from mlmm_toolkit.extract import load_structure

    missing = tmp_path / "no_such_input.pdb"
    with pytest.raises(FileNotFoundError):
        load_structure(str(missing), "missing")


def test_resolve_atom_spec_index_rejects_invalid_token_count():
    from mlmm_toolkit.utils import resolve_atom_spec_index

    atom_meta = [{"resname": "ALA", "resseq": 1, "name": "CA"}]
    with pytest.raises(ValueError, match="exactly 3 fields"):
        resolve_atom_spec_index("ALA-1-CA", atom_meta)


def test_resolve_atom_spec_index_rejects_ambiguous_match():
    from mlmm_toolkit.utils import resolve_atom_spec_index

    atom_meta = [
        {"resname": "ALA", "resseq": 1, "name": "CA"},
        {"resname": "ALA", "resseq": 1, "name": "CA"},
    ]
    with pytest.raises(ValueError, match="matches 2 atoms"):
        resolve_atom_spec_index("ALA 1 CA", atom_meta)


def test_load_pdb_atom_metadata_handles_missing_columns(tmp_path: Path):
    from mlmm_toolkit.utils import load_pdb_atom_metadata

    pdb_path = tmp_path / "short_cols.pdb"
    pdb_path.write_text("ATOM      1  N   ALA\nEND\n", encoding="utf-8")

    meta = load_pdb_atom_metadata(pdb_path)
    assert len(meta) == 1
    assert meta[0]["serial"] == 1
    assert meta[0]["name"] == "N"


def test_parse_ligand_charge_option_rejects_bad_mapping_token():
    from mlmm_toolkit.extract import _parse_ligand_charge_option

    with pytest.raises(ValueError, match="Invalid --ligand-charge token"):
        _parse_ligand_charge_option("GPP")


def test_parse_res_tokens_rejects_empty_specification():
    from mlmm_toolkit.extract import _parse_res_tokens

    with pytest.raises(ValueError, match="Empty -c/--center"):
        _parse_res_tokens("   ")
