"""Shared configurable frequency-zero cutoff contracts."""

from __future__ import annotations

import numpy as np
import pytest
import torch

from mlmm.core.defaults import FREQ_KW
from pysisyphus.normal_modes import (
    filter_resolved_modes,
    normalize_frequency_zero_cutoff_cm,
    resolved_imaginary_mask,
)


def test_default_cutoff_is_five_cm_inverse() -> None:
    assert FREQ_KW["zero_cutoff_cm"] == 5.0


def test_cutoff_is_configurable_and_boundary_is_removed() -> None:
    frequencies = np.array([-7.01, -7.0, -6.99, 0.0, 7.0, 7.01])
    modes = torch.eye(len(frequencies))

    filtered, filtered_modes = filter_resolved_modes(frequencies, modes, 7.0)

    np.testing.assert_allclose(filtered, [-7.01, 7.01])
    assert tuple(filtered_modes.shape) == (2, 6)
    assert resolved_imaginary_mask(frequencies, 7.0).tolist() == [
        True, False, False, False, False, False
    ]


@pytest.mark.parametrize("value", [-1.0, float("nan"), float("inf")])
def test_cutoff_rejects_invalid_values(value) -> None:
    with pytest.raises(ValueError):
        normalize_frequency_zero_cutoff_cm(value)
