"""Activation contract for the multi-structure ``all`` pre-alignment step.

Regression guard for C4b / M27. The ``all`` pre-alignment block (global
coordinate continuity across path-opt segments, ``mlmm/workflows/all.py``
"Global pre-alignment for coordinate continuity across segments") builds its
shared ML/MM calculator through the unified request-local template
(``_stage_calc_kwargs``).

Before C4b the block hand-built its calc kwargs WITHOUT ``model_pdb``, so
``mlmm(**kw)`` reached ``MLMMCore.__init__`` where
``shutil.copy(model_pdb, ...)`` received ``None`` and raised ``TypeError``.
The surrounding ``except Exception`` then silently skipped pre-alignment, so
the step had been dead since v0.3.0. C4b routes the pre-align site through
``_stage_calc_kwargs`` like every other calculator site, which always threads
``model_pdb`` -- activating the step.

These tests pin the fixed contract at the calculator-kwargs boundary (a fast,
GPU-free CI guard). End-to-end runtime activation is validated separately by a
real multi-structure ``all`` path-opt smoke.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import yaml

from mlmm.backends.mlmm_calc import mlmm
from mlmm.workflows.all import _resolve_calculator_template, _stage_calc_kwargs


def _prealign_stage_kwargs(tmp_path: Path) -> dict:
    """Mirror the all.py pre-align call site exactly (use_bfactor_layers=True)."""
    cfg = tmp_path / "args.yaml"
    cfg.write_text(yaml.safe_dump({"calc": {"backend": "uma"}}))
    template = _resolve_calculator_template(
        cfg,
        backend=None,
        embedcharge=False,
        embedcharge_explicit=False,
        embedcharge_cutoff=None,
        link_atom_method=None,
        mm_backend=None,
        use_cmap=None,
    )
    return _stage_calc_kwargs(
        template,
        input_pdb=tmp_path / "pocket_000.pdb",
        real_parm7=tmp_path / "real.parm7",
        model_pdb=tmp_path / "ml_region.pdb",
        charge=0,
        spin=1,
        use_bfactor_layers=True,
    )


def test_prealign_calc_kwargs_thread_model_pdb(tmp_path: Path) -> None:
    """The exact regression: the pre-align calc must carry a real model_pdb.

    Base omitted model_pdb -> None -> MLMMCore shutil.copy(None) TypeError ->
    dead block. The fixed contract threads the ML-region path.
    """
    kw = _prealign_stage_kwargs(tmp_path)
    assert "model_pdb" in kw
    assert kw["model_pdb"] == str(tmp_path / "ml_region.pdb")
    assert kw["model_pdb"] not in (None, "None")
    # state identity for the pre-align site
    assert kw["input_pdb"] == str(tmp_path / "pocket_000.pdb")
    assert kw["real_parm7"] == str(tmp_path / "real.parm7")
    assert kw["model_charge"] == 0
    assert kw["model_mult"] == 1
    assert kw["use_bfactor_layers"] is True


def test_prealign_calc_kwargs_do_not_reintroduce_the_dead_construction(
    tmp_path: Path,
) -> None:
    """The pre-align kwargs must be constructible by ``mlmm`` without the base
    failure mode: ``model_pdb`` is a real ``mlmm`` parameter and is non-None, so
    MLMMCore's unconditional ``shutil.copy(model_pdb, ...)`` receives a path.
    """
    kw = _prealign_stage_kwargs(tmp_path)
    params = inspect.signature(mlmm.__init__).parameters
    accepts_var_kw = any(
        p.kind is inspect.Parameter.VAR_KEYWORD for p in params.values()
    )
    unknown = set(kw) - set(params)
    assert accepts_var_kw or not unknown, f"pre-align kwargs mlmm rejects: {unknown}"
    # model_pdb must map to a declared parameter (routed to MLMMCore), non-None.
    assert "model_pdb" in params
    assert kw["model_pdb"] is not None
