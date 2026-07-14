"""Validation tests for the path-derived TS reference mode."""

from __future__ import annotations

import numpy as np
import pytest

from mlmm.workflows.tsopt import _load_reference_mode


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
