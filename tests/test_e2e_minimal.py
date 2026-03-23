"""Minimal end-to-end smoke tests for lightweight workflow paths."""

from __future__ import annotations

import sys
from pathlib import Path

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


def _write_pdb(path: Path, *, drop_last_atom: bool = False) -> Path:
    lines = [
        _atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N", "ATOM"),
        _atom_line(2, "CA", "ALA", "A", 1, 1.46, 0.0, 0.0, "C", "ATOM"),
        _atom_line(3, "C1", "GPP", "A", 2, 2.50, 0.0, 0.0, "C", "HETATM"),
        _atom_line(4, "O1", "GPP", "A", 2, 3.20, 0.3, 0.0, "O", "HETATM"),
    ]
    if drop_last_atom:
        lines = lines[:-1]
    lines.extend(["TER\n", "END\n"])
    path.write_text("".join(lines), encoding="utf-8")
    return path


def test_extract_api_single_structure_smoke(tmp_path: Path):
    from mlmm.extract import extract_api

    inp = _write_pdb(tmp_path / "input.pdb")
    out = tmp_path / "pocket.pdb"

    result = extract_api(
        complex_pdb=[str(inp)],
        center="GPP",
        output=[str(out)],
        radius=2.6,
        radius_het2het=0.0,
        include_h2o=False,
        exclude_backbone=False,
        add_linkh=False,
        ligand_charge="-2",
        verbose=False,
    )

    assert out.exists()
    assert result["outputs"] == [str(out)]
    assert len(result["counts"]) == 1
    assert result["counts"][0]["kept_atoms"] <= result["counts"][0]["raw_atoms"]
    assert "charge_summary" in result


def test_extract_api_multi_structure_rejects_atom_count_mismatch(tmp_path: Path):
    from mlmm.extract import extract_api

    p1 = _write_pdb(tmp_path / "s1.pdb")
    p2 = _write_pdb(tmp_path / "s2.pdb", drop_last_atom=True)

    with pytest.raises(ValueError, match="Atom count mismatch"):
        extract_api(
            complex_pdb=[str(p1), str(p2)],
            center="GPP",
            output=[str(tmp_path / "multi.pdb")],
            add_linkh=False,
            verbose=False,
        )


def test_extract_api_rejects_unknown_center(tmp_path: Path):
    from mlmm.extract import extract_api

    inp = _write_pdb(tmp_path / "input_unknown_center.pdb")

    with pytest.raises(ValueError, match="Residue name 'NOPE' not found"):
        extract_api(
            complex_pdb=[str(inp)],
            center="NOPE",
            output=[str(tmp_path / "pocket_unknown.pdb")],
            add_linkh=False,
            verbose=False,
        )
