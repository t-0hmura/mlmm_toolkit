"""Pin the default rigid-mode treatment.

``tr_projection`` decides whether frozen-atom freq / IRC / Dimer-seed / TS
validation project rigid modes out of the active block. The selected treatment
can change frequencies, imaginary-mode counts, and derived thermochemistry, so
the product and bundled engine must share one explicit default.

The default is ``constrained``: it projects only an actual null space of the
constrained problem.  ``legacy-active`` -- projecting the ACTIVE FRAGMENT's own
TR out of a Hessian whose frozen block makes those directions non-null -- is a
Rayleigh-Ritz compression biased toward a SMALLER n_imag, so it can hide a real
imaginary mode.  It is retained only as a named opt-in for reproducing output
made with the earlier treatment.
"""

from __future__ import annotations

import pytest

from pysisyphus.tr_projection import (
    DEFAULT_TR_PROJECTION,
    TR_PROJECTION_MODES,
    normalize_tr_projection_mode,
)


def test_default_is_the_correct_projector() -> None:
    """The default must project only an actual null space, never the active fragment's own TR."""
    assert DEFAULT_TR_PROJECTION == "constrained"


def test_unset_resolves_to_the_default() -> None:
    """``--tr-projection`` unset must resolve to the declared default."""
    assert normalize_tr_projection_mode(None) == DEFAULT_TR_PROJECTION


def test_geometry_defaults_agree_with_the_single_source_of_truth() -> None:
    """The product's GEOM defaults must not fork from the engine's default."""
    from mlmm.core.defaults import GEOM_KW_DEFAULT

    assert GEOM_KW_DEFAULT["tr_projection"] == DEFAULT_TR_PROJECTION


def test_engine_fallbacks_agree_with_the_single_source_of_truth() -> None:
    """IRC and the TS optimizer must not carry their own literal fallback.

    Both read ``getattr(self.geometry, "tr_projection", DEFAULT_TR_PROJECTION)``.
    A geometry without the attribute must land on the same default as everything
    else, not on a stale literal.
    """
    import inspect

    from pysisyphus.irc import IRC
    from pysisyphus.tsoptimizers import TSHessianOptimizer

    for module in (IRC, TSHessianOptimizer):
        src = inspect.getsource(module)
        assert '"tr_projection", "constrained"' not in src, (
            f"{module.__name__} carries a literal default that can drift from "
            f"DEFAULT_TR_PROJECTION"
        )


def test_legacy_remains_available_as_an_explicit_opt_in() -> None:
    """Legacy projection output remains reproducible through explicit opt-in."""
    assert "legacy-active" in TR_PROJECTION_MODES
    assert normalize_tr_projection_mode("legacy-active") == "legacy-active"
    assert normalize_tr_projection_mode("LEGACY-ACTIVE") == "legacy-active"


def test_unknown_mode_is_rejected_loudly() -> None:
    with pytest.raises(ValueError, match="Unknown TR projection mode"):
        normalize_tr_projection_mode("no-such-mode")
