"""M26: one immutable MLMM residue catalog.

Charge inference (``mlmm.workflows.extract``) and element inference
(``mlmm.domain.add_elem_info``) must read the SAME canonical residue tables
(``mlmm.core.residue_data``), and ``--modified-residue`` must never mutate the
canonical table.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm requires Python >= 3.11",
)


def _atom_line(serial, atom, resname, chain, resseq, x, y, z, element, record="ATOM"):
    return (
        f"{record:<6}{serial:>5} {atom:<4} {resname:>3} {chain:1}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{20.00:>6.2f}          {element:>2}\n"
    )


def _write_pdb(path: Path) -> Path:
    lines = [
        _atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N", "ATOM"),
        _atom_line(2, "CA", "ALA", "A", 1, 1.46, 0.0, 0.0, "C", "ATOM"),
        _atom_line(3, "C1", "GPP", "A", 2, 2.50, 0.0, 0.0, "C", "HETATM"),
        _atom_line(4, "O1", "GPP", "A", 2, 3.20, 0.3, 0.0, "O", "HETATM"),
        "TER\n",
        "END\n",
    ]
    path.write_text("".join(lines), encoding="utf-8")
    return path


# Representative key/charge pairs for the terminal, ion and D-amino partitions.
_TERMINAL_REPS = {"CGLU": -2, "NLYS": +2, "CTER": -1, "NTER": +1}
_ION_REPS = {"NA": +1, "ZN": +2, "FE": +3, "CL": -1}
_D_AMINO_REPS = {"DAL": 0, "DAR": +1, "DAS": -1}


def test_extract_and_element_inference_share_the_same_canonical_tables():
    """Falsifier 1: both inference paths expose identical key/charge pairs."""
    import mlmm.core.residue_data as rd
    import mlmm.workflows.extract as extract
    import mlmm.domain.add_elem_info as add_elem

    canonical = dict(rd.AMINO_ACIDS)

    # extract's working copy is a faithful copy of the canonical table.
    assert dict(extract.AMINO_ACIDS) == canonical
    # element inference derives its protein-residue set from the same table.
    assert set(add_elem.PROTEIN_RES) == set(canonical.keys())

    # Representative charge pairs agree across the exposed table.
    for rn, q in {**_TERMINAL_REPS, **_D_AMINO_REPS}.items():
        assert canonical[rn] == q
        assert extract.AMINO_ACIDS[rn] == q
    for rn, q in _ION_REPS.items():
        assert rd.ION[rn] == q
        assert extract.ION[rn] == q
    for water in ("HOH", "WAT", "SOL"):
        assert water in rd.WATER_RES and water in extract.WATER_RES

    # Terminal-name sets are shared frozensets, not local drift.
    assert extract.C_TERMINAL_RESNAMES is rd.C_TERMINAL_RESNAMES
    assert extract.N_TERMINAL_RESNAMES is rd.N_TERMINAL_RESNAMES


def test_canonical_tables_are_read_only():
    """Falsifier 3: direct mutation of the canonical catalog is a typed error."""
    import mlmm.core.residue_data as rd

    with pytest.raises(TypeError):
        rd.AMINO_ACIDS["ZZZ"] = 7  # type: ignore[index]
    with pytest.raises(TypeError):
        rd.ION["ZZZ"] = 7  # type: ignore[index]
    with pytest.raises(AttributeError):
        rd.WATER_RES.add("ZZZ")  # type: ignore[attr-defined]


def test_two_modified_residue_requests_leave_canonical_identical(tmp_path: Path):
    """Falsifier 2: modified-residue requests never corrupt the canonical table.

    Runs the REAL ``extract_api`` code path twice with ``--modified-residue`` and
    asserts the canonical mapping is byte/value-identical to a fresh snapshot and
    the working copy is restored after every request.
    """
    import mlmm.core.residue_data as rd
    import mlmm.workflows.extract as extract
    from mlmm.workflows.extract import extract_api

    fresh_canonical = dict(rd.AMINO_ACIDS)
    fresh_working = dict(extract.AMINO_ACIDS)
    assert "ZZ9" not in fresh_canonical and "ZZ8" not in fresh_canonical

    inp = _write_pdb(tmp_path / "input.pdb")

    for i, spec in enumerate(("ZZ9:-2", "ZZ8:1")):
        out = tmp_path / f"pocket_{i}.pdb"
        extract_api(
            complex_pdb=[str(inp)],
            center="GPP",
            output=[str(out)],
            radius=2.6,
            radius_het2het=0.0,
            modified_residue=spec,
            verbose=False,
        )
        # Canonical table is immutable and thus byte/value-identical to fresh.
        assert dict(rd.AMINO_ACIDS) == fresh_canonical
        # Working copy is restored (the request-local addition is gone).
        assert dict(extract.AMINO_ACIDS) == fresh_working
        assert "ZZ9" not in extract.AMINO_ACIDS
        assert "ZZ8" not in extract.AMINO_ACIDS
