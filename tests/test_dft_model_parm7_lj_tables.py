"""Regression tests for the shared ML-region model-parm7 builder (v0.3.3).

Every ONIOM ``E_model_low`` -- ML/MM through ``MLMMCore`` and QM/MM through ``mlmm
dft`` -- must describe the model region identically, or the DFT correction (a
difference against the ML/MM energy) carries whatever the two paths disagree on.
``dft`` used to slice and save its own model and had drifted on both the CMAP drop
and the LJ tables; ``write_model_parm7`` is now the single owner and both call it.

Covers:
  * ParmEd leaves ``LENNARD_JONES_*COEF`` at the *parent's* length whenever the selection
    uses fewer atom types than the full system, producing a ``model.parm7`` that mlmm's
    own MM backend refuses to load; ``_normalize_prmtop_lj_tables`` repairs it.
  * Why the smoke suite never caught this: its ML region happens to span every atom type
    in the fixture, so the stale table is accidentally the right length. A test that only
    exercised that region would pass against the broken builder.
  * A second builder cannot reappear (anti-divergence guard).

All checks bind to the production functions (imported, never reimplemented).
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

parmed = pytest.importorskip("parmed")

from mlmm.backends.mlmm_calc import (
    MLMMCore,
    _normalize_prmtop_lj_tables,
    apply_cmap_policy,
    write_model_parm7,
)
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


def test_exactly_one_owner_of_the_cmap_and_lj_rules():
    """Every ONIOM layer must go through one implementation of each rule.

    dft used to slice and save its own model, which is how it drifted away from
    the ML/MM path on both CMAP and the LJ tables. Pin the single owner so a
    second builder cannot reappear.
    """
    root = Path(__file__).resolve().parents[1] / "mlmm"
    writers = []
    for path in sorted(root.rglob("*.py")):
        text = path.read_text(encoding="utf-8")
        if "cmaps[:] = []" in text or "_normalize_prmtop_lj_tables(" in text:
            writers.append(str(path.relative_to(root)))
    assert writers == ["backends/mlmm_calc.py"], writers


def test_both_oniom_layers_share_one_cmap_policy():
    """`E_real_low` and `E_model_low` must come from the same MM Hamiltonian.

    A CMAP term lying entirely inside the ML region only cancels in the
    subtractive total when both layers treat CMAP the same way. Stripping the
    model alone left it uncancelled -- an empirical backbone potential stacked
    on the high-level description.
    """
    core = inspect.getsource(MLMMCore.__init__)
    assert "apply_cmap_policy(real_top, use_cmap)" in core, (
        "the ML/MM real layer no longer applies the shared CMAP policy"
    )
    dft_src = inspect.getsource(_prepare_ml_region_workspace)
    assert "apply_cmap_policy(real_top, use_cmap)" in dft_src, (
        "the dft real layer no longer applies the shared CMAP policy"
    )
    assert "apply_cmap_policy(model, use_cmap)" in inspect.getsource(write_model_parm7)


def test_cmap_policy_default_is_off_and_strips_both_layers():
    top_real = parmed.load_file(str(FIXTURE))
    assert len(top_real.cmaps) > 0, "fixture must carry CMAP terms to be meaningful"
    apply_cmap_policy(top_real, False)
    assert len(top_real.cmaps) == 0

    kept = parmed.load_file(str(FIXTURE))
    apply_cmap_policy(kept, True)
    assert len(kept.cmaps) > 0


def test_dft_and_mlmm_both_delegate_to_the_shared_builder():
    for consumer in (_prepare_ml_region_workspace, MLMMCore._mk_model_parm7):
        src = inspect.getsource(consumer)
        assert "write_model_parm7(" in src, (
            f"{consumer.__name__} no longer uses the shared model-parm7 builder"
        )


def test_shared_builder_honours_use_cmap():
    assert "use_cmap" in inspect.signature(write_model_parm7).parameters
    assert "use_cmap" in inspect.signature(_prepare_ml_region_workspace).parameters
    assert "use_cmap" in inspect.signature(apply_cmap_policy).parameters
    # The rule itself lives in apply_cmap_policy; the builder must delegate.
    assert "if not use_cmap:" in inspect.getsource(apply_cmap_policy)
