"""Error-handling regression tests (invalid files, YAML, and malformed inputs)."""

from __future__ import annotations

import sys
from pathlib import Path

import click
import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm requires Python >= 3.11",
)


@pytest.mark.parametrize("spec", ["1-a", "0", "2-1"])
def test_parse_indices_string_rejects_invalid(spec):
    from mlmm.core.utils import parse_indices_string

    with pytest.raises(click.BadParameter):
        parse_indices_string(spec)


def test_load_yaml_dict_rejects_non_mapping_root(tmp_path: Path):
    from mlmm.core.utils import load_yaml_dict

    p = tmp_path / "bad_root.yaml"
    p.write_text("- a\n- b\n", encoding="utf-8")

    with pytest.raises(ValueError, match="YAML root must be a mapping"):
        load_yaml_dict(p)


def test_load_yaml_dict_rejects_malformed_yaml(tmp_path: Path):
    from mlmm.core.utils import load_yaml_dict

    p = tmp_path / "malformed.yaml"
    p.write_text("geom: [1,2\n", encoding="utf-8")

    with pytest.raises(Exception):
        load_yaml_dict(p)


def test_collect_single_option_values_rejects_repeated_flags():
    from mlmm.core.utils import collect_single_option_values

    argv = ["-i", "a.pdb", "-i", "b.pdb"]
    with pytest.raises(click.BadParameter, match="single -i/--input"):
        collect_single_option_values(argv, ("-i", "--input"), label="-i/--input")


def test_collect_option_values_accepts_grouped_and_repeated_flags():
    from mlmm.core.utils import collect_option_values

    argv = ["-i", "a.pdb", "b.cif", "--output", "x.pdb", "-i", "c.pdb"]
    assert collect_option_values(argv, ("-i", "--input")) == [
        "a.pdb",
        "b.cif",
        "c.pdb",
    ]
    assert collect_option_values(argv, ("-o", "--output")) == ["x.pdb"]


def test_collect_option_values_preserves_equals_attached_and_grouped_order():
    from mlmm.core.utils import collect_option_values

    argv = [
        "--input=first.pdb", "-i", "second.pdb", "third.pdb",
        "--output=first.out", "-osecond.out", "--ref-pdb=template.pdb",
        "--input=fourth.pdb",
    ]
    assert collect_option_values(argv, ("-i", "--input")) == [
        "first.pdb", "second.pdb", "third.pdb", "fourth.pdb",
    ]
    assert collect_option_values(argv, ("-o", "--output")) == [
        "first.out", "second.out",
    ]
    assert collect_option_values(argv, ("--ref-pdb",)) == ["template.pdb"]


def test_load_structure_rejects_missing_input_file(tmp_path: Path):
    from mlmm.workflows.extract import load_structure

    missing = tmp_path / "no_such_input.pdb"
    with pytest.raises(FileNotFoundError):
        load_structure(str(missing), "missing")


def test_resolve_atom_spec_index_rejects_invalid_token_count():
    from mlmm.core.utils import resolve_atom_spec_index

    atom_meta = [{"resname": "ALA", "resseq": 1, "name": "CA"}]
    with pytest.raises(ValueError, match="must have 3 fields.*or 4 fields"):
        resolve_atom_spec_index("ALA-1-CA", atom_meta)


def test_resolve_atom_spec_index_rejects_ambiguous_match():
    from mlmm.core.utils import resolve_atom_spec_index

    atom_meta = [
        {"resname": "ALA", "resseq": 1, "name": "CA"},
        {"resname": "ALA", "resseq": 1, "name": "CA"},
    ]
    with pytest.raises(ValueError, match="matches 2 atoms"):
        resolve_atom_spec_index("ALA 1 CA", atom_meta)


def test_load_pdb_atom_metadata_handles_missing_columns(tmp_path: Path):
    from mlmm.core.utils import load_pdb_atom_metadata

    pdb_path = tmp_path / "short_cols.pdb"
    pdb_path.write_text("ATOM      1  N   ALA\nEND\n", encoding="utf-8")

    meta = load_pdb_atom_metadata(pdb_path)
    assert len(meta) == 1
    assert meta[0]["serial"] == 1
    assert meta[0]["name"] == "N"


def test_parse_ligand_charge_option_rejects_bad_mapping_token():
    import click
    from mlmm.workflows.extract import _parse_ligand_charge_option

    # click.BadParameter provides a clean CLI error display.
    with pytest.raises(click.BadParameter, match="Invalid --ligand-charge token"):
        _parse_ligand_charge_option("GPP")


@pytest.mark.parametrize("value", ["HEM", "HEM=x", "HEM=0"])
def test_parse_ligand_mult_rejects_bad_mapping_as_click_error(value):
    import click

    from mlmm.workflows.mm_parm import parse_ligand_mult

    with pytest.raises(click.BadParameter, match="ligand-mult"):
        parse_ligand_mult(value)


def test_parse_res_tokens_rejects_empty_specification():
    from mlmm.workflows.extract import _parse_res_tokens

    # ValueError is intentional — `resolve_substrate_residues` upstream
    # catches this for the ID-vs-name dispatch fallback.
    with pytest.raises(ValueError, match="Empty -c/--center"):
        _parse_res_tokens("   ")


def test_bond_summary_reports_positional_and_malformed_input_cleanly(
    tmp_path: Path,
):
    from click.testing import CliRunner

    from mlmm.cli import cli as root_cli

    one = tmp_path / "one.xyz"
    one.write_text("1\nframe\nH 0 0 0\n", encoding="utf-8")
    too_few = CliRunner().invoke(root_cli, ["bond-summary", str(one)])
    assert too_few.exit_code != 0
    assert "Invalid value for inputs" in too_few.output
    assert "'-i'" not in too_few.output

    broken = tmp_path / "broken.xyz"
    broken.write_text("not an xyz\n", encoding="utf-8")
    malformed = CliRunner().invoke(
        root_cli, ["bond-summary", str(one), str(broken)]
    )
    assert malformed.exit_code != 0
    assert "Could not read structure" in malformed.output
    assert "Traceback" not in malformed.output
