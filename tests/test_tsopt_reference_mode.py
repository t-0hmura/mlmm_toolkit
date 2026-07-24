"""Validation tests for the path-derived TS reference mode."""

from __future__ import annotations

import numpy as np
import pytest
import click

from mlmm.core.defaults import TSOPT_MODE_ALIASES
from mlmm.core.utils import normalize_choice
from mlmm.workflows.tsopt import (
    _load_reference_mode,
    _validate_reference_mode_optimizer,
)


def test_reference_mode_is_normalized(tmp_path) -> None:
    path = tmp_path / "mode.txt"
    np.savetxt(path, [3.0, 4.0, 0.0])

    np.testing.assert_allclose(_load_reference_mode(path, 3), [0.6, 0.8, 0.0])


def test_reference_mode_rejects_wrong_size(tmp_path) -> None:
    path = tmp_path / "mode.npy"
    np.save(path, np.ones(4))

    with pytest.raises(ValueError, match="4 read, 3 expected"):
        _load_reference_mode(path, 3)


def test_reference_mode_rejects_zero_vector(tmp_path) -> None:
    path = tmp_path / "mode.txt"
    np.savetxt(path, np.zeros(3))

    with pytest.raises(ValueError, match="finite, non-zero"):
        _load_reference_mode(path, 3)


@pytest.mark.parametrize("mode", ["dimer", "grad", "light"])
def test_reference_mode_rejects_dimer_aliases(mode, tmp_path) -> None:
    path = tmp_path / "mode.txt"
    path.write_text("1 0 0\n", encoding="utf-8")
    normalized = normalize_choice(
        mode,
        param="--opt-mode",
        alias_groups=TSOPT_MODE_ALIASES,
        allowed_hint="grad|hess|dimer|rsirfo|trim|rsprfo",
    )

    with pytest.raises(click.BadParameter, match="requires a Hessian TS optimizer"):
        _validate_reference_mode_optimizer(normalized, path)


@pytest.mark.parametrize("mode", ["hess", "rsirfo", "rsprfo", "trim"])
def test_reference_mode_accepts_hessian_optimizers(mode, tmp_path) -> None:
    path = tmp_path / "mode.txt"
    normalized = normalize_choice(
        mode,
        param="--opt-mode",
        alias_groups=TSOPT_MODE_ALIASES,
        allowed_hint="grad|hess|dimer|rsirfo|trim|rsprfo",
    )

    _validate_reference_mode_optimizer(normalized, path)
