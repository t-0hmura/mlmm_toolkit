"""Hessian TS optimizer construction in the microiteration path."""

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
        assert "min_line_search" not in kw
        assert "max_line_search" not in kw
    elif mode == "rsprfo":
        assert kw["min_line_search"] is False
        assert kw["max_line_search"] is False
    else:
        assert "min_line_search" not in kw
        assert "max_line_search" not in kw

    opt = tsopt_mod.TSOPT_CLASS_MAP[mode](_tiny_geom(), **kw)
    assert type(opt).__name__ in ("RSIRFOptimizer", "RSPRFOptimizer", "TRIM")


def test_rsprfo_honors_explicit_line_search_values(tsopt_mod, tmp_path):
    kw = tsopt_mod._build_rsirfo_kwargs(
        {"min_line_search": True, "max_line_search": True},
        max_cycles=1,
        out_dir=tmp_path,
        mode="rsprfo",
    )

    assert kw["min_line_search"] is True
    assert kw["max_line_search"] is True


@pytest.mark.parametrize("mode", ["rsirfo", "trim"])
def test_unused_line_search_values_are_removed(tsopt_mod, mode, tmp_path):
    kw = tsopt_mod._build_rsirfo_kwargs(
        {"min_line_search": True, "max_line_search": True},
        max_cycles=1,
        out_dir=tmp_path,
        mode=mode,
    )

    assert "min_line_search" not in kw
    assert "max_line_search" not in kw


def test_hessian_ts_kwargs_require_one_root(tsopt_mod, tmp_path):
    with pytest.raises(tsopt_mod.click.BadParameter, match="exactly one root"):
        tsopt_mod._build_rsirfo_kwargs(
            {"roots": [0, 1]},
            max_cycles=1,
            out_dir=tmp_path,
            mode="rsprfo",
        )


def test_ts_kwargs_drop_ordinary_rfo_overlap_tracking(tsopt_mod, tmp_path):
    kw = tsopt_mod._build_rsirfo_kwargs(
        {"rfo_overlaps": True},
        max_cycles=1,
        out_dir=tmp_path,
        mode="rsprfo",
    )

    assert "rfo_overlaps" not in kw
