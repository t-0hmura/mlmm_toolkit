"""Regression tests for the ``mlmm dft`` model-parm7 builder (v0.3.3).

``dft._prepare_ml_region_workspace`` and ``mlmm_calc._mk_model_parm7`` build the same
artefact -- the sliced ML-region ``model.parm7`` -- and had silently diverged: the dft
builder applied neither the CMAP drop nor the LJ-table normalization that the ML/MM
builder has always applied.

Covers:
  * ParmEd leaves ``LENNARD_JONES_*COEF`` at the *parent's* length whenever the selection
    uses fewer atom types than the full system, producing a ``model.parm7`` that mlmm's
    own MM backend refuses to load; ``_normalize_prmtop_lj_tables`` repairs it.
  * Why the smoke suite never caught this: its ML region happens to span every atom type
    in the fixture, so the stale table is accidentally the right length. A test that only
    exercised that region would pass against the broken builder.
  * Both builders stay in step (anti-divergence guard).

All checks bind to the production functions (imported, never reimplemented).
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

parmed = pytest.importorskip("parmed")

from mlmm.backends.mlmm_calc import MLMMCore, _normalize_prmtop_lj_tables
from mlmm.core.utils import parse_layer_indices_from_bfactors
from mlmm.workflows.dft import _prepare_ml_region_workspace

FIXTURE = Path(__file__).parent / "smoke" / "p_complex.parm7"
LAYERED_PDB = Path(__file__).parent / "smoke" / "r_complex_layered.pdb"
for _fixture in (FIXTURE, LAYERED_PDB):
    assert _fixture.exists(), f"committed smoke fixture is missing: {_fixture}"

# Residues 3-5 are a contiguous backbone stretch: it uses 10 of the fixture's 18 atom
# types, and retains exactly one complete CMAP term.
BACKBONE_RESIDUES = {3, 4, 5}


def _slice(residues):
    top = parmed.load_file(str(FIXTURE))
    selection = [a.idx for r in top.residues if r.idx in residues for a in r.atoms]
    model = top[selection]
    model.box = None
    return top, model


def _written_lj_lengths(path):
    raw = parmed.amber.AmberFormat(str(path))
    ntypes = raw.parm_data["POINTERS"][1]
    return ntypes, len(raw.parm_data["LENNARD_JONES_ACOEF"])


def test_sliced_model_parm7_needs_lj_normalization(tmp_path):
    """A fewer-atom-types selection writes an internally inconsistent parm7."""
    _, model = _slice(BACKBONE_RESIDUES)
    out = tmp_path / "model.parm7"
    model.save(str(out), overwrite=True)

    ntypes, acoef = _written_lj_lengths(out)
    assert acoef != ntypes * (ntypes + 1) // 2, (
        "fixture no longer reproduces the ParmEd slicing defect; "
        f"NTYPES={ntypes} ACOEF={acoef}"
    )
    with pytest.raises(Exception):
        parmed.load_file(str(out))

    _normalize_prmtop_lj_tables(str(out))

    ntypes, acoef = _written_lj_lengths(out)
    assert acoef == ntypes * (ntypes + 1) // 2
    assert len(parmed.load_file(str(out)).atoms) == len(model.atoms)


def test_normalization_is_a_noop_for_a_full_atom_type_selection(tmp_path):
    """The smoke's ML region spans every atom type, which is why it never caught this.

    Pinning the near-miss keeps a future reviewer from concluding the smoke lane covers
    the defect.
    """
    top = parmed.load_file(str(FIXTURE))
    full_ntypes = top.parm_data["POINTERS"][1]
    bfactors = [
        float(line[60:66])
        for line in LAYERED_PDB.read_text().splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]
    selection = sorted(parse_layer_indices_from_bfactors(bfactors)["ml_indices"])
    model = top[selection]
    model.box = None
    out = tmp_path / "model.parm7"
    model.save(str(out), overwrite=True)

    ntypes, acoef = _written_lj_lengths(out)
    assert ntypes == full_ntypes, (
        "the near-miss no longer holds: this slice must span every atom type for the "
        "test to document why the smoke lane stayed green against the broken builder"
    )
    assert acoef == ntypes * (ntypes + 1) // 2
    _normalize_prmtop_lj_tables(str(out))
    assert _written_lj_lengths(out) == (ntypes, acoef)


def test_both_model_parm7_builders_apply_the_same_normalizations():
    """dft and ML/MM build the same artefact and must not diverge again."""
    for builder in (_prepare_ml_region_workspace, MLMMCore._mk_model_parm7):
        src = inspect.getsource(builder)
        assert "model.cmaps[:] = []" in src, f"{builder.__name__} lost the CMAP drop"
        assert "_normalize_prmtop_lj_tables(" in src, (
            f"{builder.__name__} lost the LJ-table normalization"
        )


def test_dft_builder_honours_use_cmap():
    assert "use_cmap" in inspect.signature(_prepare_ml_region_workspace).parameters
    src = inspect.getsource(_prepare_ml_region_workspace)
    assert "if not use_cmap:" in src
