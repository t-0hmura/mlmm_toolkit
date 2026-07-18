"""PDB serial gaps must not break sequential parm7 atom lookups."""

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

    from mlmm.backends.mlmm_calc import MLMMCore
    from mlmm.io.pdb_indexing import parse_pdb_ordinal_atoms

    atoms = parse_pdb_ordinal_atoms(pdb)
    assert [atom.idx for atom in atoms] == [1, 2, 3, 4, 5]
    assert [atom.serial for atom in atoms] == [1, 2, 3, 8, 9]

    # Exercise production `_ml_prep` without running the heavyweight
    # constructor.  Taking the whole input as the model avoids link atoms.
    core = MLMMCore.__new__(MLMMCore)
    core.input_pdb = str(pdb)
    core.model_pdb = str(pdb)
    core.link_mlmm = None
    ml_ids, links, element_pairs = core._ml_prep()

    assert ml_ids == ["1", "2", "3", "4", "5"]
    assert links == []
    assert element_pairs == []


@_need_py311
def test_dft_loader_delegates_to_ordinal_parser(tmp_path: Path, monkeypatch) -> None:
    from mlmm.io.pdb_indexing import parse_pdb_ordinal_atoms
    from mlmm.workflows import dft

    pdb = tmp_path / "gap.pdb"
    _write_pdb_with_gap(pdb)
    expected = parse_pdb_ordinal_atoms(pdb)
    calls = []

    def recording_parser(path):
        calls.append(Path(path))
        return expected

    monkeypatch.setattr(dft, "parse_pdb_ordinal_atoms", recording_parser)

    assert dft._load_input_atoms(pdb) is expected
    assert calls == [pdb]


@_need_py311
def test_ordinal_parser_accepts_hybrid36_and_opaque_diagnostic_serials(
    tmp_path: Path,
) -> None:
    from mlmm.io.pdb_indexing import parse_pdb_ordinal_atoms

    source = tmp_path / "hybrid36.pdb"
    base = (
        "ATOM      1  C   MOL A   1       0.000   0.000   0.000  1.00  0.00           C\n"
    )
    source.write_text(
        base[:6] + "A0000" + base[11:]
        + base[:6] + "A0002" + base[11:]
        + base[:6] + "*****" + base[11:]
        + "END\n",
        encoding="utf-8",
    )

    atoms = parse_pdb_ordinal_atoms(source)

    assert [atom.idx for atom in atoms] == [1, 2, 3]
    assert [atom.serial for atom in atoms] == [100_000, 100_002, "*****"]
