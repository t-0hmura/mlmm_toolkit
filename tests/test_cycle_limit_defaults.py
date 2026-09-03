"""Regression coverage for finite product cycle defaults."""

from __future__ import annotations

import click
import pytest

from mlmm.core.defaults import (
    DMF_KW,
    DFT_KW,
    IRC_KW,
    MICROITER_KW,
    OPT_BASE_KW,
    STOPT_KW,
)
from mlmm.core.utils import optional_positive_int


def test_internal_cycle_defaults_keep_finite_backstops() -> None:
    assert OPT_BASE_KW["max_cycles"] == 100000
    assert DMF_KW["max_cycles"] == 3000
    assert STOPT_KW["max_cycles"] == 300
    assert STOPT_KW["stop_in_when_full"] == 300
    assert IRC_KW["max_cycles"] == 125
    assert DFT_KW["max_cycle"] == 100
    assert MICROITER_KW["micro_max_cycles"] == 100000


def test_bundled_engine_defaults_remain_finite() -> None:
    from inspect import signature

    from pysisyphus.irc.IRC import IRC
    from pysisyphus.optimizers.Optimizer import Optimizer

    assert signature(Optimizer).parameters["max_cycles"].default == 150
    assert signature(IRC).parameters["max_cycles"].default == 125


def test_none_is_uncapped_and_zero_is_invalid() -> None:
    assert optional_positive_int(None, "cycles") is None
    assert optional_positive_int(5, "cycles") == 5
    with pytest.raises(click.BadParameter):
        optional_positive_int(0, "cycles")
