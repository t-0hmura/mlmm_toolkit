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
    filter_info = {}

    filtered, filtered_modes = filter_resolved_modes(
        frequencies,
        modes,
        7.0,
        filter_info=filter_info,
    )

    np.testing.assert_allclose(filtered, [-7.01, 7.01])
    assert tuple(filtered_modes.shape) == (2, 6)
    assert filter_info == {
        "frequency_zero_cutoff_cm": 7.0,
        "raw_mode_count": 6,
        "resolved_mode_count": 2,
        "near_zero_mode_count": 4,
        "near_zero_frequencies_cm": [-7.0, -6.99, 0.0, 7.0],
    }
    assert resolved_imaginary_mask(frequencies, 7.0).tolist() == [
        True, False, False, False, False, False
    ]


def test_mode_accounting_is_one_compact_line() -> None:
    from mlmm.workflows.freq import _format_mode_accounting

    line = _format_mode_accounting(
        410,
        414,
        {
            "effective_rank": 3,
            "frequency_zero_cutoff_cm": 5.0,
            "near_zero_mode_count": 1,
            "near_zero_frequencies_cm": [-4.2],
        },
    )

    assert line == (
        "410 modes = 414 active DOF - 3 rigid - 1 near-zero "
        "(|ν|≤5.0 cm⁻¹)"
    )
    assert "\n" not in line
    assert "raw" not in line


@pytest.mark.parametrize("value", [-1.0, float("nan"), float("inf")])
def test_cutoff_rejects_invalid_values(value) -> None:
    with pytest.raises(ValueError):
        normalize_frequency_zero_cutoff_cm(value)
