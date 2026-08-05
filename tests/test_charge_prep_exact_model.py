from pathlib import Path

import click
import pytest

from mlmm.core.utils import PreparedInputStructure
from mlmm.workflows.charge_prep import resolve_charge_spin_or_raise


def _atom_line(
    serial, atom, resname, chain, resseq, x, y, z, element, bfactor
):
    return (
        f"ATOM  {serial:>5} {atom:<4} {resname:>3} {chain}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{bfactor:>6.2f}"
        f"          {element:>2}\n"
    )


def test_layer_charge_derivation_includes_physical_terminal_cap(tmp_path):
    pdb = tmp_path / "layered.pdb"
    pdb.write_text(
        "".join(
            [
                _atom_line(1, "N", "ALA", "A", 1, 0, 0, 0, "N", 0),
                _atom_line(2, "H1", "ALA", "A", 1, -1, 0, 0, "H", 0),
                _atom_line(3, "H2", "ALA", "A", 1, 0, -1, 0, "H", 0),
                _atom_line(4, "H3", "ALA", "A", 1, 0, 0, -1, "H", 0),
                _atom_line(5, "CA", "ALA", "A", 1, 1.4, 0, 0, "C", 0),
                _atom_line(6, "C", "ALA", "A", 1, 2.5, 0, 0, "C", 0),
                _atom_line(7, "O", "ALA", "A", 1, 3.5, 0, 0, "O", 0),
                _atom_line(8, "MG", "MG", "B", 2, 8, 0, 0, "MG", 20),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )
    prepared = PreparedInputStructure(source_path=pdb, geom_path=pdb)

    charge, spin = resolve_charge_spin_or_raise(
        prepared,
        charge=None,
        spin=None,
        ligand_charge="",
        detect_layer=True,
    )

    assert (charge, spin) == (1, 1)


def test_layer_charge_derivation_rejects_residue_split(tmp_path):
    pdb = tmp_path / "split.pdb"
    pdb.write_text(
        _atom_line(1, "N", "ALA", "A", 1, 0, 0, 0, "N", 0)
        + _atom_line(2, "CA", "ALA", "A", 1, 1, 0, 0, "C", 10)
        + "END\n",
        encoding="utf-8",
    )
    prepared = PreparedInputStructure(source_path=pdb, geom_path=pdb)

    with pytest.raises(click.ClickException, match="splits residue.*-q/--charge"):
        resolve_charge_spin_or_raise(
            prepared,
            charge=None,
            spin=1,
            ligand_charge="ALA:0",
            detect_layer=True,
        )


def test_model_indices_require_explicit_charge_for_exact_selection():
    prepared = PreparedInputStructure(
        source_path=Path("dummy.pdb"),
        geom_path=Path("dummy.pdb"),
    )

    with pytest.raises(Exception, match="--model-indices"):
        resolve_charge_spin_or_raise(
            prepared,
            charge=None,
            spin=1,
            ligand_charge="LIG:-1",
            model_indices_spec="1-3",
        )


def test_explicit_charge_bypasses_automatic_model_selection():
    prepared = PreparedInputStructure(
        source_path=Path("dummy.pdb"),
        geom_path=Path("dummy.pdb"),
    )

    assert resolve_charge_spin_or_raise(
        prepared,
        charge=-2,
        spin=1,
        ligand_charge="LIG:-1",
        model_indices_spec="1-3",
    ) == (-2, 1)


@pytest.mark.parametrize(
    ("yaml_cfg", "expected"),
    [
        ({"calc": {"model_charge": -2, "model_mult": 3}}, (-2, 3)),
        ({"mlmm": {"model_charge": -1, "model_mult": 2}}, (-1, 2)),
        (
            {
                "calc": {"model_charge": -3, "model_mult": 4},
                "mlmm": {"model_charge": 1, "model_mult": 5},
            },
            (-3, 4),
        ),
    ],
)
def test_yaml_charge_spin_precedence(yaml_cfg, expected):
    prepared = PreparedInputStructure(
        source_path=Path("dummy.pdb"),
        geom_path=Path("dummy.pdb"),
    )

    assert resolve_charge_spin_or_raise(
        prepared,
        charge=None,
        spin=None,
        yaml_cfg=yaml_cfg,
    ) == expected


def test_canonical_calculator_section_shadows_legacy_per_key_fallback():
    prepared = PreparedInputStructure(
        source_path=Path("dummy.pdb"),
        geom_path=Path("dummy.pdb"),
    )

    with pytest.raises(click.ClickException, match="charge is unresolved"):
        resolve_charge_spin_or_raise(
            prepared,
            charge=None,
            spin=None,
            yaml_cfg={
                "calc": {"model_charge": None, "model_mult": None},
                "mlmm": {"model_charge": -1, "model_mult": 2},
            },
        )


@pytest.mark.parametrize("invalid", [True, False, 1.5, "1"])
def test_yaml_charge_rejects_non_integer_types(invalid):
    prepared = PreparedInputStructure(
        source_path=Path("dummy.pdb"),
        geom_path=Path("dummy.pdb"),
    )

    with pytest.raises(Exception, match="must be an integer"):
        resolve_charge_spin_or_raise(
            prepared,
            charge=None,
            spin=None,
            yaml_cfg={"calc": {"model_charge": invalid, "model_mult": 1}},
        )


def test_explicit_cli_charge_spin_override_yaml():
    prepared = PreparedInputStructure(
        source_path=Path("dummy.pdb"),
        geom_path=Path("dummy.pdb"),
    )

    assert resolve_charge_spin_or_raise(
        prepared,
        charge=2,
        spin=5,
        yaml_cfg={"calc": {"model_charge": -2, "model_mult": 3}},
    ) == (2, 5)


@pytest.mark.parametrize(
    ("yaml_cfg", "expected_name"),
    [
        ({"calc": {"model_pdb": "canonical.pdb"}}, "canonical.pdb"),
        ({"mlmm": {"model_pdb": "legacy.pdb"}}, "legacy.pdb"),
        (
            {
                "calc": {"model_pdb": "canonical.pdb"},
                "mlmm": {"model_pdb": "legacy.pdb"},
            },
            "canonical.pdb",
        ),
    ],
)
def test_yaml_model_pdb_is_the_charge_derivation_source(
    tmp_path,
    monkeypatch,
    yaml_cfg,
    expected_name,
):
    from mlmm.workflows import charge_prep

    prepared = PreparedInputStructure(
        source_path=tmp_path / "layered-full.pdb",
        geom_path=tmp_path / "layered-full.pdb",
    )
    seen = {}

    def derive(path, ligand_charge, *, select_bfactor_layer, prefix=""):
        seen.update(
            path=Path(path),
            select_bfactor_layer=select_bfactor_layer,
        )
        return -3

    monkeypatch.setattr(
        charge_prep,
        "_derive_charge_from_ligand_charge",
        derive,
    )

    assert resolve_charge_spin_or_raise(
        prepared,
        charge=None,
        spin=1,
        ligand_charge="LIG:-1",
        yaml_cfg=yaml_cfg,
    ) == (-3, 1)
    assert seen["path"] == Path(expected_name)
    assert seen["select_bfactor_layer"] is False


def test_explicit_model_pdb_overrides_yaml_for_charge_derivation(
    tmp_path,
    monkeypatch,
):
    from mlmm.workflows import charge_prep

    prepared = PreparedInputStructure(
        source_path=tmp_path / "layered-full.pdb",
        geom_path=tmp_path / "layered-full.pdb",
    )
    explicit = tmp_path / "explicit.pdb"
    seen = {}

    def derive(path, ligand_charge, *, select_bfactor_layer, prefix=""):
        seen["path"] = Path(path)
        return 2

    monkeypatch.setattr(
        charge_prep,
        "_derive_charge_from_ligand_charge",
        derive,
    )

    assert resolve_charge_spin_or_raise(
        prepared,
        charge=None,
        spin=1,
        ligand_charge="LIG:-1",
        model_pdb=explicit,
        yaml_cfg={"calc": {"model_pdb": "yaml.pdb"}},
    ) == (2, 1)
    assert seen["path"] == explicit


def test_invalid_yaml_multiplicity_is_rejected():
    prepared = PreparedInputStructure(
        source_path=Path("dummy.pdb"),
        geom_path=Path("dummy.pdb"),
    )

    with pytest.raises(Exception, match="multiplicity"):
        resolve_charge_spin_or_raise(
            prepared,
            charge=None,
            spin=None,
            yaml_cfg={"calc": {"model_charge": 0, "model_mult": 0}},
        )
