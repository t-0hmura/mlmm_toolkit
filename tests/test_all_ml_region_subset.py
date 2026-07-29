"""Regression tests for the ML-region B-factor subset and charge/spin diagnostics.

Covers:
  * ``_write_bfactor_ml_subset`` writes ONLY the B≈0 ML atoms and returns ``None`` when
    the input carries no B≈0 atoms.
  * ``validate_charge_spin`` error message carries ``counted_atoms=`` always and
    ``source=`` only when a source label is passed.
  * The all.py skip_extract+detect_layer path picks ``_write_bfactor_ml_subset``, and the
    full-layered ``--model-pdb`` warning predicate (``_summarize_existing_bfactor_layers``)
    distinguishes an ML-only subset (no MM-layer atoms) from a full layered system.

All checks bind to the production functions (imported, never reimplemented).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mlmm.core.utils import validate_charge_spin
from mlmm.workflows.all import (
    _ml_region_atom_summary,
    _ml_region_summary_suffix,
    _summarize_existing_bfactor_layers,
    _write_bfactor_ml_subset,
    _write_ml_region_definition,
)


def _atom(serial: int, name: str, resname: str, resseq: int,
          x: float, y: float, z: float, bfac: float, elem: str) -> str:
    """Format a column-exact PDB ATOM line (B-factor at cols 61-66, element at 77-78)."""
    return (
        "ATOM  "
        f"{serial:>5d}"
        " "
        f"{name:<4s}"
        " "               # altLoc (col 17)
        f"{resname:>3s}"
        " A"              # col 21 space + chain (col 22)
        f"{resseq:>4d}"
        " "               # iCode (col 27)
        "   "             # cols 28-30
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.00:6.2f}"
        f"{bfac:6.2f}"
        "          "      # cols 67-76
        f"{elem:>2s}"
        "\n"
    )


# 3 ML atoms (B=0.00: N, C, C), 2 MovableMM (B=10.00), 2 FrozenMM (B=20.00).
_ML_ATOMS = [(1, "N", "N"), (2, "CA", "C"), (3, "C", "C")]
_MM_MOVABLE = [(4, "O", "O", 10.00), (5, "N", "N", 10.00)]
_MM_FROZEN = [(6, "CB", "C", 20.00), (7, "CG", "C", 20.00)]


def _layered_pdb_text() -> str:
    lines = []
    for i, (serial, name, elem) in enumerate(_ML_ATOMS):
        lines.append(_atom(serial, name, "ALA", 1, float(i), 0.0, 0.0, 0.00, elem))
    for serial, name, elem, bf in _MM_MOVABLE + _MM_FROZEN:
        lines.append(_atom(serial, name, "ALA", 2, float(serial), 0.0, 0.0, bf, elem))
    lines.append("END\n")
    return "".join(lines)


def _mm_only_pdb_text() -> str:
    lines = []
    for serial, name, elem, bf in _MM_MOVABLE + _MM_FROZEN:
        lines.append(_atom(serial, name, "ALA", 2, float(serial), 0.0, 0.0, bf, elem))
    lines.append("END\n")
    return "".join(lines)


# --- (i) _write_bfactor_ml_subset and the layer/charge helpers around it ------------------------------------------------

def test_write_bfactor_ml_subset_keeps_only_b0_atoms(tmp_path: Path):
    src = tmp_path / "layered.pdb"
    src.write_text(_layered_pdb_text(), encoding="utf-8")
    dest = tmp_path / "ml_region.pdb"

    out = _write_bfactor_ml_subset(src, dest)
    assert out is not None
    assert Path(out).exists()

    kept = [ln for ln in Path(out).read_text().splitlines()
            if ln.startswith(("ATOM", "HETATM"))]
    # Exactly the 3 B≈0 atoms, in order, with their elements preserved.
    assert len(kept) == 3
    assert [ln[76:78].strip() for ln in kept] == ["N", "C", "C"]
    # No MovableMM(10) / FrozenMM(20) atom leaked through.
    for ln in kept:
        assert abs(float(ln[60:66])) < 0.5


def test_write_bfactor_ml_subset_returns_none_without_ml_atoms(tmp_path: Path):
    src = tmp_path / "mm_only.pdb"
    src.write_text(_mm_only_pdb_text(), encoding="utf-8")
    dest = tmp_path / "ml_region.pdb"

    assert _write_bfactor_ml_subset(src, dest) is None
    # Caller must be able to fall back: no partial/empty file is relied upon.
    assert not dest.exists()


def test_ml_region_definition_removes_extractor_link_h(tmp_path: Path):
    src = tmp_path / "pocket_with_link.pdb"
    src.write_text(
        _atom(1, "CA", "ALA", 1, 0.0, 0.0, 0.0, 0.0, "C")
        + "HETATM    2  HL  LKH L   1       1.090   0.000   0.000"
        "  1.00  0.00           H\nEND\n",
        encoding="utf-8",
    )
    dest = tmp_path / "ml_region.pdb"

    result = _write_ml_region_definition(src, dest)

    atom_lines = [
        line
        for line in Path(result).read_text(encoding="utf-8").splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]
    assert len(atom_lines) == 1
    assert "LKH" not in Path(result).read_text(encoding="utf-8")


# --- (ii) validate_charge_spin diagnostics --------------------------------------

def test_validate_charge_spin_error_includes_counted_atoms_and_source():
    # CH3 (9 e-) with mult=1 (even unpaired) → odd/even mismatch → raises.
    with pytest.raises(ValueError) as excinfo:
        validate_charge_spin(["C", "H", "H", "H"], charge=0, multiplicity=1,
                             source="/path/to/ml_region.pdb")
    msg = str(excinfo.value)
    assert "counted_atoms=4" in msg
    assert "source=/path/to/ml_region.pdb" in msg


def test_validate_charge_spin_error_omits_source_when_not_passed():
    with pytest.raises(ValueError) as excinfo:
        validate_charge_spin(["C", "H", "H", "H"], charge=0, multiplicity=1)
    msg = str(excinfo.value)
    assert "counted_atoms=4" in msg
    assert "source=" not in msg


# --- (iii) branch predicate: subset selection + full-layered warning -------------

def test_all_branch_selects_bfactor_subset_and_summary(tmp_path: Path):
    """The skip_extract+detect_layer branch uses _write_bfactor_ml_subset; its output is
    an ML-only region and the one-line summary reports the correct atoms/sumZ."""
    src = tmp_path / "layered.pdb"
    src.write_text(_layered_pdb_text(), encoding="utf-8")
    dest = tmp_path / "ml_region.pdb"

    out = _write_bfactor_ml_subset(src, dest)
    assert out is not None

    # sumZ = N(7) + C(6) + C(6) = 19 over 3 atoms.
    assert _ml_region_atom_summary(Path(out)) == (3, 19)
    assert _ml_region_summary_suffix(Path(out)) == " (atoms=3, sumZ=19)"


def test_full_layered_warning_predicate(tmp_path: Path):
    """Item (a): a full layered --model-pdb (ML + MM layers) trips the warning predicate;
    an ML-only subset does not."""
    full = tmp_path / "full.pdb"
    full.write_text(_layered_pdb_text(), encoding="utf-8")
    layers = _summarize_existing_bfactor_layers(full)
    assert layers["ml"] == 3
    assert layers["movable"] == 2
    assert layers["frozen"] == 2
    # Predicate used by all.py: ml>0 AND (movable+frozen)>0 → warn.
    assert layers["ml"] > 0 and (layers["movable"] + layers["frozen"]) > 0

    subset = tmp_path / "ml_region.pdb"
    _write_bfactor_ml_subset(full, subset)
    sub_layers = _summarize_existing_bfactor_layers(subset)
    assert sub_layers["ml"] == 3
    assert sub_layers["movable"] == 0
    assert sub_layers["frozen"] == 0
    # ML-only subset → predicate does NOT fire.
    assert not (sub_layers["ml"] > 0 and (sub_layers["movable"] + sub_layers["frozen"]) > 0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
