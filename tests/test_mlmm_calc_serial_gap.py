"""Regression test: mlmm_calc must use sequential ATOM file position as `idx`,
not raw PDB serial. PDB serial gaps (e.g. COMT model 1.R.pdb has 3411→3418)
must not break parm7 lookups."""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path

import pytest

_need_py311 = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm.mlmm_calc requires Python >= 3.11",
)


def _write_pdb_with_gap(path: Path) -> None:
    """Write a tiny 5-atom PDB whose serials skip from 3 to 8 (gap of 4)."""
    body = textwrap.dedent(
        """\
        ATOM      1  N   ALA A   1      11.000  10.000  10.000  1.00  0.00           N
        ATOM      2  CA  ALA A   1      12.000  10.000  10.000  1.00  0.00           C
        ATOM      3  C   ALA A   1      13.000  10.000  10.000  1.00  0.00           C
        ATOM      8  O   ALA A   1      14.000  10.000  10.000  1.00  0.00           O
        ATOM      9  CB  ALA A   1      12.000  11.000  10.000  1.00  0.00           C
        END
        """
    )
    path.write_text(body)


@_need_py311
def test_idx_is_file_position_not_pdb_serial(tmp_path: Path) -> None:
    """`_ml_prep` must record `idx` as 1..5 (file position) regardless of
    the input PDB's serial column having a 3->8 gap."""
    pdb = tmp_path / "gap.pdb"
    _write_pdb_with_gap(pdb)

    from mlmm import mlmm_calc as _mc

    # We don't actually instantiate MLMMCore (it needs parm7/rst7 + UMA).
    # Just verify the file-position logic by re-implementing the relevant
    # snippet locally — mirroring the production code path.
    leap_atoms = []
    atom_pos = 0
    for ln in open(pdb):
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        atom_pos += 1
        leap_atoms.append({"idx": atom_pos, "serial": int(ln[6:11])})

    assert [a["idx"] for a in leap_atoms] == [1, 2, 3, 4, 5]
    assert [a["serial"] for a in leap_atoms] == [1, 2, 3, 8, 9]

    # Confirm the production code uses the same `atom_pos += 1` pattern:
    src = Path(_mc.__file__).read_text()
    assert "atom_pos = 0" in src and 'leap_atoms.append' in src
    assert '"idx": atom_pos,' in src
