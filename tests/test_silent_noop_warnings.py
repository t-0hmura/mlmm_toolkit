"""Flags and settings that do nothing must say so.

Each case below was a surface where the user asked for something, the code
could not honour it, and the run continued in silence -- so the result read as
if the request had been applied. Output-only guards: no behaviour, threshold or
default changes with them.
"""

from pathlib import Path

import mlmm.workflows.freq as freq_mod
import mlmm.workflows.mm_parm as mm_parm_mod


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
