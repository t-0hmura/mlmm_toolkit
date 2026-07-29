"""The IRC active basis preserves ML + MovableMM and is never
silently cropped to ML + link-parent DOFs, and the IRC Hessian device policy is
GPU-first with an explicit CUDA request that never silently falls back to CPU.

The ordered-basis crop was the defect: reducing the seeded Hessian to the
ML-macro (ML + link-parent) sub-block drops every MovableMM coordinate, so the
bundled integrator's ``_full`` expansion writes zero displacement into those
atoms.  These controls exercise the real pysisyphus ``IRC._full`` scatter and the
real device-policy helper.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11), reason="mlmm requires Python >= 3.11"
)

from pysisyphus.irc.IRC import IRC

from mlmm.workflows._microiteration import resolve_hessian_device

IRC_SRC = Path(__file__).resolve().parents[1] / "mlmm" / "workflows" / "irc.py"


class _CoordsGeom:
    """Geometry stub exposing only the Cartesian coord vector ``_full`` needs."""

    def __init__(self, n_atoms):
        self._coords = np.zeros(3 * n_atoms, dtype=float)

    @property
    def coords(self):
        return self._coords


def _irc_with_basis(n_atoms, act_dofs):
    irc = IRC.__new__(IRC)
    irc.geometry = _CoordsGeom(n_atoms)
    irc._act_dofs = np.asarray(act_dofs, dtype=int)
    return irc


def test_full_writes_displacement_into_movable_mm_when_basis_preserved():
    # atoms: 0 = ML, 1 = link parent, 2 = MovableMM (3 atoms, 9 DOFs).
    # A PRESERVED ML+MovableMM basis includes atom 2's DOFs.
    preserved = [0, 1, 2, 3, 4, 5, 6, 7, 8]  # all three atoms active
    irc = _irc_with_basis(3, preserved)
    full = irc._full(np.ones(len(preserved)))
    movable_dofs = [6, 7, 8]
    assert np.all(full[movable_dofs] != 0.0), "MovableMM atom must receive displacement"


def test_full_zeros_movable_mm_under_the_old_ml_macro_crop():
    # The old crop keeps only ML + link-parent DOFs (atoms 0,1) and drops the
    # MovableMM atom (2). This control pins WHY the crop was wrong: _full writes
    # zero into the dropped MovableMM coordinates.
    cropped = [0, 1, 2, 3, 4, 5]  # atoms 0,1 only
    irc = _irc_with_basis(3, cropped)
    full = irc._full(np.ones(len(cropped)))
    movable_dofs = [6, 7, 8]
    assert np.all(full[movable_dofs] == 0.0)


def test_irc_source_has_no_ml_macro_crop_path():
    # The ML-macro crop must remain absent.
    text = IRC_SRC.read_text(encoding="utf-8")
    for banned in ("_macro_atoms", "_try_reduce", "ML macro sub-block"):
        assert banned not in text, banned


# --------------------------------------------------------------------------
# Device policy — GPU-first; explicit CUDA never silently selects CPU
# --------------------------------------------------------------------------


def test_explicit_cuda_stays_cuda_when_available():
    dev, reason = resolve_hessian_device("cuda", cuda_available=True)
    assert dev == "cuda" and reason == "explicit_cuda"


def test_explicit_cuda_errors_when_unavailable_never_silent_cpu():
    with pytest.raises(ValueError):
        resolve_hessian_device("cuda", cuda_available=False)


def test_auto_is_gpu_first():
    assert resolve_hessian_device("auto", cuda_available=True) == ("cuda", "auto_gpu_first")
    assert resolve_hessian_device("auto", cuda_available=False) == ("cpu", "auto_no_cuda")


def test_explicit_cpu_is_a_calibrated_choice():
    assert resolve_hessian_device("cpu", cuda_available=True) == ("cpu", "explicit_cpu")
