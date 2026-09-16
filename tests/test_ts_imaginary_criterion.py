"""Physical classification uses the original eigenvalue criterion, not raw signs."""
from types import SimpleNamespace

import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import eigval_to_wavenumber
from pysisyphus.normal_modes import (
    DEFAULT_FREQUENCY_ZERO_CUTOFF_CM, frequency_partition_info,
    resolved_imaginary_mask,
)
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer


def test_direct_optimizer_constructor_uses_geometry_imaginary_default(tmp_path):
    geometry = Geometry(["H"], [0., 0., 0.])
    optimizer = RSPRFOptimizer(geometry, hessian_init="unit", dump=False,
                              out_dir=tmp_path)
    assert optimizer.saddle_imaginary_threshold_cm == eigval_to_wavenumber(1e-6)
    assert optimizer.small_eigval_thresh == 1e-8
    values = np.array([-1.01e-6, -1e-6, -.99e-6, 1e-9])
    # Geometry.get_imag_frequencies owns the original strict eigenvalue test.
    frequencies = eigval_to_wavenumber(values)
    geometry.get_normal_modes = lambda _hessian: (frequencies, values, None, None)
    np.testing.assert_array_equal(
        frequencies[resolved_imaginary_mask(frequencies, optimizer.saddle_imaginary_threshold_cm)],
        geometry.get_imag_frequencies(),
    )


@pytest.mark.parametrize("second,selected,verified", [(-.99e-6, 1, True), (-1.01e-6, 2, False)])
def test_exact_phva_keeps_raw_two_and_classifies_selected_count(second, selected, verified):
    coordinates = np.zeros(3)
    frequencies = eigval_to_wavenumber(np.array([-1e-3, second, 1e-6]))
    modes = np.eye(3)
    optimizer = RSPRFOptimizer.__new__(RSPRFOptimizer)
    optimizer.geometry = SimpleNamespace(cart_coords=coordinates.copy())
    optimizer.saddle_imaginary_threshold_cm = DEFAULT_FREQUENCY_ZERO_CUTOFF_CM
    optimizer.roots = [0]
    optimizer.reference_mode = None
    optimizer.cur_cycle = 7
    printed, stopped = [], []
    optimizer.table = SimpleNamespace(print=printed.append)
    optimizer.request_stop = stopped.append
    optimizer._last_rigid_projection_info = frequency_partition_info(frequencies)
    optimizer._mw_frequencies_and_modes = lambda: (frequencies, modes)
    optimizer._recovery_mode_from_mw = lambda _modes, index: _modes[index].copy()
    optimizer._record_exact_saddle_candidate = lambda: None
    has_modes, _, performed = optimizer._verify_exact_vibrational_structure(None, None)
    assert performed and has_modes
    assert optimizer._last_exact_n_negative == 2
    assert optimizer._last_exact_n_imaginary == selected
    assert optimizer._last_exact_saddle_verified is verified
    assert optimizer._last_exact_validation == ("first_order" if verified else "higher_order")
    assert optimizer.higher_order_saddle_checks == int(not verified)
    assert optimizer._exact_saddle_matches_current_geometry() is verified
    assert not stopped
    assert any(f"n_imag={selected}, n_negative=2" in line for line in printed)
    np.testing.assert_array_equal(optimizer.geometry.cart_coords, coordinates)
    optimizer.geometry.cart_coords[0] = 1e-3
    assert not optimizer._exact_saddle_matches_current_geometry()


def test_incomplete_phva_still_clears_saddle_diagnostics():
    optimizer = RSPRFOptimizer.__new__(RSPRFOptimizer)
    optimizer.geometry = SimpleNamespace(cart_coords=np.zeros(3))
    frequencies = eigval_to_wavenumber(np.array([-1e-3, -.99e-6, 1e-6]))
    optimizer._mw_frequencies_and_modes = lambda: (frequencies, np.eye(3))
    optimizer._last_rigid_projection_info = frequency_partition_info(frequencies)
    optimizer._last_rigid_projection_info["raw_mode_count"] += 1
    optimizer._last_exact_saddle_verified = True
    optimizer.table = SimpleNamespace(print=lambda *_args: None)
    stopped = []
    optimizer.request_stop = stopped.append
    optimizer._verify_exact_vibrational_structure(None, None)
    assert optimizer._last_exact_validation == "unavailable"
    assert optimizer._last_exact_saddle_verified is False
    assert optimizer._last_exact_n_negative is None
    assert optimizer._last_exact_n_imaginary is None
    assert len(stopped) == 1
