"""Local peak selection and the preceding basin of an MEP energy profile."""

import numpy as np
from scipy.signal import find_peaks

from pysisyphus.constants import AU2KCALPERMOL


def select_mep_peak(energies):
    """Return (peak, preceding basin); a monotone profile selects an endpoint."""
    values = np.asarray(energies, dtype=float)
    if values.size == 0:
        raise ValueError("Cannot select an HEI from an empty energy profile.")
    if not np.all(np.isfinite(values)):
        raise ValueError("Cannot select an HEI from non-finite energies.")
    peaks, _ = find_peaks(values)
    if not peaks.size:
        peak = len(values) - 1 if values[-1] > values[0] else 0
        return peak, 0
    peak = int(peaks[np.argmax(values[peaks])])
    minima, _ = find_peaks(-values[:peak + 1])
    basin = int(minima[-1]) if minima.size else 0
    return peak, basin


def mep_barrier_kcal(energies):
    """Return the selected local peak minus its preceding basin, in kcal/mol."""
    peak, basin = select_mep_peak(energies)
    if not 0 < peak < len(energies) - 1:
        return None
    return float(energies[peak] - energies[basin]) * AU2KCALPERMOL
