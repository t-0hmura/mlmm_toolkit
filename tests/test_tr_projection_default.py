"""Pin the default rigid-mode treatment.

``tr_projection`` decides whether frozen-atom freq / IRC / Dimer-seed / TS
validation project rigid modes out of the active block. The selected treatment
can change frequencies, imaginary-mode counts, and derived thermochemistry, so
the product and bundled engine must share one explicit default.

The only treatment is ``constrained``: it projects only an actual null space of
the constrained problem. The superseded ``legacy-active`` treatment -- projecting
the active fragment's own TR out of a Hessian whose frozen block makes those
directions non-null -- was a Rayleigh-Ritz compression biased toward a smaller
n_imag, so it could hide a real imaginary mode. It has been removed.
"""

from __future__ import annotations

import pytest

from pysisyphus.tr_projection import (
    DEFAULT_TR_PROJECTION,
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


def test_unknown_mode_is_rejected_loudly() -> None:
    with pytest.raises(ValueError, match="Unknown TR projection mode"):
        normalize_tr_projection_mode("no-such-mode")
