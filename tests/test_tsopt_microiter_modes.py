"""F1 regression: the microiter macro loop must build for every Hessian TS optimizer.

Microiteration was historically RS-I-RFO-only. It now also drives RS-P-RFO and TRIM
macro steps. These tests pin the two invariants that make that safe without running a
full MLIP optimization:

  * ``TSOPT_CLASS_MAP`` exposes all three Hessian TS optimizers, and
  * ``_build_rsirfo_kwargs(mode=...)`` disables the RS-I-RFO-only torch line search for
    TRIM / RS-P-RFO (which reject it) while leaving it intact for RS-I-RFO, and each
    resulting kwarg set constructs its optimizer.
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm requires Python >= 3.11",
)


@pytest.fixture(scope="module")
def tsopt_mod():
    return pytest.importorskip("mlmm.workflows.tsopt")


def _tiny_geom():
    """A minimal 3-atom geometry (bohr) — enough to construct an optimizer."""
    from pysisyphus.Geometry import Geometry

    ang2bohr = 1.0 / 0.529177210903
    coords = np.array([0, 0, 0, 0.96, 0, 0, -0.24, 0.93, 0], dtype=float) * ang2bohr
    return Geometry(atoms=("O", "H", "H"), coords=coords)


def test_tsopt_class_map_covers_all_hessian_modes(tsopt_mod):
    assert set(tsopt_mod.TSOPT_CLASS_MAP) == {"rsirfo", "rsprfo", "trim"}


@pytest.mark.parametrize("mode", ["rsirfo", "rsprfo", "trim"])
def test_microiter_macro_optimizer_builds(tsopt_mod, mode, tmp_path):
    from mlmm.core.defaults import RSIRFO_KW

    kw = tsopt_mod._build_rsirfo_kwargs(
        dict(RSIRFO_KW),
        max_cycles=1,
        out_dir=tmp_path,
        macro_thresh="baker",
        mode=mode,
    )
    if mode == "rsirfo":
        # For RS-I-RFO the builder leaves min_line_search alone.
        assert kw.get("min_line_search") is not False
    else:
        # For TRIM / RS-P-RFO the macro kwargs pin the line search off.
        assert kw.get("min_line_search") is False
        assert kw.get("max_line_search") is False

    opt = tsopt_mod.TSOPT_CLASS_MAP[mode](_tiny_geom(), **kw)
    assert type(opt).__name__ in ("RSIRFOptimizer", "RSPRFOptimizer", "TRIM")
