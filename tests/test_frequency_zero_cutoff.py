"""Shared configurable frequency-zero cutoff contracts."""

from __future__ import annotations

import numpy as np
import pytest
import torch

from mlmm.core.defaults import FREQ_KW
from pysisyphus.normal_modes import (
    DEFAULT_FREQUENCY_ZERO_CUTOFF_CM,
    DEFAULT_IMAGINARY_EIGENVALUE_THRESHOLD,
    frequency_criterion_info,
    warn_legacy_frequency_cutoff,
    filter_resolved_modes,
    normalize_frequency_zero_cutoff_cm,
    resolved_imaginary_mask,
)


def test_default_cutoff_derives_from_original_mass_weighted_eigenvalue() -> None:
    from pysisyphus.helpers_pure import eigval_to_wavenumber

    assert DEFAULT_IMAGINARY_EIGENVALUE_THRESHOLD == 1e-6
    assert FREQ_KW["zero_cutoff_cm"] == eigval_to_wavenumber(1e-6)
    assert 5.1 < DEFAULT_FREQUENCY_ZERO_CUTOFF_CM < 5.2
    eigenvalues = np.array([-1.01e-6, -1e-6, -.99e-6, -1e-12, 0., 1e-12])
    frequencies = eigval_to_wavenumber(eigenvalues)
    np.testing.assert_array_equal(resolved_imaginary_mask(frequencies), eigenvalues < -1e-6)
    # The old rounded 5 cm^-1 rule gives a different decision near the bound.
    assert resolved_imaginary_mask(frequencies, 5.0)[2]
    info = frequency_criterion_info()
    assert info["imaginary_mode_criterion"] == "mass_weighted_eigenvalue"
    assert info["imaginary_eigenvalue_threshold"] == 1e-6
    assert info["imaginary_eigenvalue_units"] == "hartree/(bohr^2*amu)"


def test_explicit_legacy_cutoff_is_visible_and_warns():
    with pytest.warns(FutureWarning, match="overrides the original"):
        warn_legacy_frequency_cutoff(5.0)
    assert frequency_criterion_info(5.0)["imaginary_mode_criterion"] == "legacy_frequency_cutoff"


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
        411,
        414,
        {
            "effective_rank": 3,
            "frequency_zero_cutoff_cm": 5.0,
            "near_zero_mode_count": 1,
            "near_zero_frequencies_cm": [-4.2],
        },
    )

    assert line == (
        "411 modes = 414 active DOF - 3 rigid; 1 near-zero retained "
        "(|ν|≤5.0 cm⁻¹)"
    )
    assert "\n" not in line
    assert "raw" not in line


@pytest.mark.parametrize("value", [-1.0, float("nan"), float("inf")])
def test_cutoff_rejects_invalid_values(value) -> None:
    with pytest.raises(ValueError):
        normalize_frequency_zero_cutoff_cm(value)


@pytest.mark.parametrize("cutoff", [0.0, 5.0])
def test_strict_count_keeps_negative_zero_window_and_signed_zero(cutoff):
    from pysisyphus.normal_modes import _strict_negative_count

    values = np.array([-10., -5., -1e-9, -0., 0., 1e-9, 5., 10.])
    before = values.copy()
    info = {}
    resolved, _ = filter_resolved_modes(values, np.eye(len(values)), cutoff, filter_info=info)
    assert _strict_negative_count(resolved, info) == 3
    assert np.count_nonzero(resolved_imaginary_mask(resolved, cutoff)) == (3 if cutoff == 0 else 1)
    np.testing.assert_array_equal(values, before)


@pytest.mark.parametrize("frequencies,info", [
    ([20.], None), ([20.], {}), ([20.], {"raw_mode_count": 1}),
    ([20.], {"raw_mode_count": 2, "near_zero_frequencies_cm": []}),
    ([20.], {"raw_mode_count": True, "near_zero_frequencies_cm": []}),
    ([20.], {"raw_mode_count": 1, "near_zero_frequencies_cm": [], "resolved_mode_count": 2}),
    ([20.], {"raw_mode_count": 2, "near_zero_frequencies_cm": [np.nan]}),
    ([np.inf], {"raw_mode_count": 1, "near_zero_frequencies_cm": []}),
    ([np.nan], {"raw_mode_count": 1, "near_zero_frequencies_cm": []}),
])
def test_incomplete_or_nonfinite_partition_cannot_certify_strict_zero(frequencies, info):
    from pysisyphus.normal_modes import _strict_negative_count

    assert _strict_negative_count(frequencies, info) is None
