"""Contracts for canonical ML-region and boundary-link identity resolution."""

from __future__ import annotations

from pathlib import Path

import pytest

from mlmm.io.pdb_indexing import resolve_mlmm_atoms


def _atom_line(
    serial: int,
    name: str,
    resname: str,
    chain: str,
    resseq: int,
    x: float,
    *,
    icode: str = "",
    altloc: str = "",
    element: str = "C",
) -> str:
    return (
        f"ATOM  {serial:5d} {name:^4s}{altloc:1s}{resname:>3s} {chain:1s}"
        f"{resseq:4d}{icode:1s}   {x:8.3f}{0.0:8.3f}{0.0:8.3f}"
        f"  1.00  0.00          {element:>2s}\n"
    )


def _write(path: Path, lines: list[str]) -> None:
    path.write_text("".join(lines) + "END\n", encoding="utf-8")


def test_model_identity_retains_chain_and_insertion_code(tmp_path: Path) -> None:
    full = tmp_path / "full.pdb"
    model = tmp_path / "model.pdb"
    lines = [
        _atom_line(1, "C1", "LIG", "A", 1, 0.0, icode="A"),
        _atom_line(2, "C1", "LIG", "A", 1, 5.0, icode="B"),
        _atom_line(3, "C1", "LIG", "B", 1, 10.0, icode="B"),
    ]
    _write(full, lines)
    _write(model, [lines[1]])

    resolved = resolve_mlmm_atoms(full, model)

    assert resolved.model_indices == (2,)


def test_manual_pairs_preserve_requested_pairing_and_order(tmp_path: Path) -> None:
    full = tmp_path / "full.pdb"
    model = tmp_path / "model.pdb"
    lines = [
        _atom_line(1, "C1", "LIG", "A", 1, 0.0),
        _atom_line(2, "CA", "ALA", "A", 2, 1.4),
        _atom_line(3, "C2", "LIG", "A", 1, 3.0),
        _atom_line(4, "CA", "GLY", "A", 3, 4.4),
    ]
    _write(full, lines)
    _write(model, [lines[0], lines[2]])

    resolved = resolve_mlmm_atoms(
        full,
        model,
        [
            ("A:LIG:1:C2", "A:GLY:3:CA"),
            ("A:LIG:1:C1", "A:ALA:2:CA"),
        ],
    )

    assert resolved.model_indices == (1, 3)
    assert resolved.link_pairs == ((3, 4), (1, 2))


def test_manual_pairs_reject_ambiguous_and_wrong_side_selectors(
    tmp_path: Path,
) -> None:
    full = tmp_path / "full.pdb"
    model = tmp_path / "model.pdb"
    lines = [
        _atom_line(1, "C1", "LIG", "A", 1, 0.0),
        _atom_line(2, "C1", "LIG", "B", 1, 4.0),
        _atom_line(3, "CA", "ALA", "A", 2, 1.4),
    ]
    _write(full, lines)
    _write(model, [lines[0]])

    with pytest.raises(ValueError, match="matches 2 atoms"):
        resolve_mlmm_atoms(full, model, [("LIG 1 C1", "A:ALA:2:CA")])
    with pytest.raises(ValueError, match="outside model_pdb"):
        resolve_mlmm_atoms(
            full,
            model,
            [("B:LIG:1:C1", "A:ALA:2:CA")],
        )


def test_model_mapping_rejects_reordering_and_raw_altloc(tmp_path: Path) -> None:
    full = tmp_path / "full.pdb"
    model = tmp_path / "model.pdb"
    lines = [
        _atom_line(1, "C1", "LIG", "A", 1, 0.0),
        _atom_line(2, "C2", "LIG", "A", 1, 1.5),
    ]
    _write(full, lines)
    _write(model, [lines[1], lines[0]])
    with pytest.raises(ValueError, match="atom order differs"):
        resolve_mlmm_atoms(full, model)

    altloc_lines = [
        _atom_line(1, "C1", "LIG", "A", 1, 0.0, altloc="A"),
    ]
    _write(full, altloc_lines)
    _write(model, altloc_lines)
    with pytest.raises(ValueError, match="nonblank altLoc"):
        resolve_mlmm_atoms(full, model)


def test_model_mapping_rejects_element_mismatch(tmp_path: Path) -> None:
    full = tmp_path / "full.pdb"
    model = tmp_path / "model.pdb"
    _write(
        full,
        [_atom_line(1, "CA", "LIG", "A", 1, 0.0, element="C")],
    )
    _write(
        model,
        [_atom_line(1, "CA", "LIG", "A", 1, 0.0, element="N")],
    )

    with pytest.raises(ValueError, match="element N.*element C"):
        resolve_mlmm_atoms(full, model)


def test_explicit_empty_manual_links_disable_auto_detection(tmp_path: Path) -> None:
    full = tmp_path / "full.pdb"
    model = tmp_path / "model.pdb"
    lines = [
        _atom_line(1, "C1", "LIG", "A", 1, 0.0),
        _atom_line(2, "CA", "ALA", "A", 2, 1.4),
    ]
    _write(full, lines)
    _write(model, [lines[0]])

    assert resolve_mlmm_atoms(full, model, None).link_pairs == ((1, 2),)
    assert resolve_mlmm_atoms(full, model, []).link_pairs == ()
