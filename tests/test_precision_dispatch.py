"""Unit tests for mlmm.backends.apply_precision_to_calc_cfg.

Pins the --precision routing contract: per-backend kwarg translation,
the fp64 ⇒ H_double coupling, and the warning emitted when a config
deliberately requests the inconsistent fp64 + H_double=False combination.
"""

from __future__ import annotations

import warnings

import pytest

from mlmm.backends import apply_precision_to_calc_cfg


def test_fp64_routes_uma_precision_kwarg() -> None:
    cfg = {"backend": "uma"}
    apply_precision_to_calc_cfg(cfg, "fp64")
    assert cfg["uma_precision"] == "fp64"


def test_fp64_default_forces_h_double_silently() -> None:
    cfg = {"backend": "uma"}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        apply_precision_to_calc_cfg(cfg, "fp64")
    assert cfg["H_double"] is True
    assert caught == []


def test_fp64_overrides_explicit_false_with_warning() -> None:
    cfg = {"backend": "uma", "H_double": False}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        apply_precision_to_calc_cfg(cfg, "fp64")
    assert cfg["H_double"] is True
    assert len(caught) == 1
    assert "fp64" in str(caught[0].message)


def test_fp32_leaves_h_double_untouched() -> None:
    cfg = {"backend": "uma", "H_double": False}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        apply_precision_to_calc_cfg(cfg, "fp32")
    assert cfg["H_double"] is False
    assert caught == []


def test_fp64_rejected_for_aimnet2() -> None:
    with pytest.raises(ValueError):
        apply_precision_to_calc_cfg({"backend": "aimnet2"}, "fp64")


def test_invalid_precision_rejected() -> None:
    with pytest.raises(ValueError):
        apply_precision_to_calc_cfg({"backend": "uma"}, "fp16")
