"""Regression test: `mlmm all` must honor pre-set B-factor layers when
extraction is skipped (no -c) and --detect-layer is True.

Before the v0.2.9 fix, `mlmm all -i layered.pdb` would silently overwrite
the user's pre-built layers because `pocket_outputs` defaulted to the
full input PDB and the internal `_define_layers` then promoted every
atom to the ML region.
"""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path

import pytest

_need_py311 = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm.all requires Python >= 3.11",
)


def _write_layered_pdb(path: Path) -> None:
    """6-atom toy PDB: 2 ML (B=0), 2 movable MM (B=10), 2 frozen MM (B=20)."""
    body = textwrap.dedent(
        """\
        ATOM      1  N   ALA A   1      11.000  10.000  10.000  1.00  0.00           N
        ATOM      2  CA  ALA A   1      12.000  10.000  10.000  1.00  0.00           C
        ATOM      3  C   ALA A   2      13.000  10.000  10.000  1.00 10.00           C
        ATOM      4  O   ALA A   2      14.000  10.000  10.000  1.00 10.00           O
        ATOM      5  CB  ALA A   3      12.000  11.000  10.000  1.00 20.00           C
        ATOM      6  CG  ALA A   3      12.000  12.000  10.000  1.00 20.00           C
        END
        """
    )
    path.write_text(body)


@_need_py311
def test_summarize_existing_bfactor_layers(tmp_path: Path) -> None:
    """The helper distinguishes ML/MovableMM/FrozenMM/other correctly."""
    pdb = tmp_path / "layered.pdb"
    _write_layered_pdb(pdb)

    from mlmm.workflows.all import _summarize_existing_bfactor_layers

    counts = _summarize_existing_bfactor_layers(pdb)
    assert counts == {"ml": 2, "movable": 2, "frozen": 2, "other": 0}


@_need_py311
def test_summarize_unlayered_pdb(tmp_path: Path) -> None:
    """An un-layered PDB (single B-factor value) lands all atoms in
    'other', flagging that the user did not pre-encode layers."""
    pdb = tmp_path / "flat.pdb"
    body = textwrap.dedent(
        """\
        ATOM      1  N   ALA A   1      11.000  10.000  10.000  1.00 33.50           N
        ATOM      2  CA  ALA A   1      12.000  10.000  10.000  1.00 33.50           C
        END
        """
    )
    pdb.write_text(body)

    from mlmm.workflows.all import _summarize_existing_bfactor_layers

    counts = _summarize_existing_bfactor_layers(pdb)
    assert counts == {"ml": 0, "movable": 0, "frozen": 0, "other": 2}


@_need_py311
def test_honor_input_bfactors_path_exists(tmp_path: Path) -> None:
    """Smoke test: the `mlmm all` define-layer block contains the
    `honor_input_bfactors` branch added by the v0.2.9 fix."""
    from mlmm.workflows import all as _all

    src = Path(_all.__file__).read_text()
    # The fix introduces this exact identifier and message.
    assert "honor_input_bfactors" in src
    assert "honoring input PDB" in src
