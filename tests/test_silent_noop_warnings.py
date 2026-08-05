"""Flags and settings that do nothing must say so.

Requested settings that cannot be honored must emit a diagnostic instead of
appearing to have been applied. Output-only guards do not change behavior,
thresholds, or defaults.
"""

from pathlib import Path

import mlmm.workflows.freq as freq_mod
import mlmm.workflows.mm_parm as mm_parm_mod
import mlmm.workflows.tsopt as tsopt_mod


def _src(mod) -> str:
    return Path(mod.__file__).read_text(encoding="utf-8")


def test_unresolvable_layer_sets_are_reported() -> None:
    # Empty layer sets collapse every --active-dof-mode to ALL atoms, and the
    # only previous signal was a debug-level log line.
    src = _src(freq_mod)
    assert "could not resolve ML/MM layer sets" in src
    assert "falls back to ALL atoms regardless of" in src


def test_mode_pdb_fallback_is_reported() -> None:
    # The file still appears, but without the reference topology its atom names
    # and residues are not the input's.
    src = _src(freq_mod)
    assert "mode PDB fell back to plain ASE output" in src
    assert "mode PDB fell back to plain ASE output" in _src(tsopt_mod)


def test_unknown_path_mode_diagnostic_is_not_reported_as_lost() -> None:
    src = _src(tsopt_mod)
    marker = src.index("Path-correlated mode was lost")
    assert "is False" in src[marker - 400 : marker]


def test_unused_ligand_charge_entries_are_reported() -> None:
    # A --ligand-charge entry for a residue tleap already knows is never read;
    # a mistyped residue name is exactly how the wrong net charge gets used.
    src = _src(mm_parm_mod)
    assert "matched no residue needing parameters" in src


def test_the_layer_fallback_keeps_its_return_contract() -> None:
    """The warning must not change what callers receive."""
    src = _src(freq_mod)
    i = src.index("could not resolve ML/MM layer sets")
    # `return empty` still follows the echo, so callers are unaffected.
    assert "return empty" in src[i : i + 600]


def test_dft_slices_the_model_with_the_same_cmap_rule_as_mlmm() -> None:
    """The CMAP policy must reach the ordinary ML/MM core used by dft.

    `dft.py` used to build its own model workspace and never read the setting,
    so with a backbone-containing ML region the ML/MM and QM/MM paths computed
    E_model_low on different models -- and `result.yaml` recorded
    `use_cmap: False` either way. dft now delegates to the shared builder, so
    the rule lives in one place and the flag has to reach it.
    """
    import inspect

    import mlmm.workflows.dft as dft_mod
    from mlmm.backends.mlmm_calc import apply_cmap_policy, write_model_parm7

    assert "if not use_cmap:" in inspect.getsource(apply_cmap_policy)
    assert "cmaps[:] = []" in inspect.getsource(apply_cmap_policy)
    assert "apply_cmap_policy(model, use_cmap)" in inspect.getsource(write_model_parm7)

    src = _src(dft_mod)
    assert "MLMMCalculator(**calculator_kwargs)" in src
    assert "workspace.core.compute(" in src
    # Threaded, not read from a global: both call sites must pass it.
    assert src.count('use_cmap=bool(calc_kw.get("use_cmap", True))') == 2
    assert "use_cmap: bool = True" in src


def test_dft_cmap_default_matches_the_product_default() -> None:
    from mlmm.core.defaults import MLMM_CALC_KW

    assert MLMM_CALC_KW["use_cmap"] is True
