"""Physical spectra and their vectors must survive reporting thresholds."""
import copy

import numpy as np
import pytest
import torch

from pysisyphus.normal_modes import (
    _frequencies_cm_and_modes, _strict_negative_count,
    frequency_partition_info, resolved_imaginary_mask,
)
from pysisyphus.tr_projection import active_tr_basis

COORDS = np.array([[0., 0., 0.], [2., 0., 0.], [0., 2., 0.],
                   [0., 0., 2.], [2., 2., 2.]])


@pytest.mark.parametrize('partial', [False, True])
def test_complete_physical_pairs_are_independent_of_reporting_cutoff(partial):
    # Three anchored atoms leave six unconstrained physical Cartesian DOF.
    # Small curvatures produce soft signed modes and an exact zero eigenpair.
    diagonal = np.array([-1e-3, -1e-9, 0., 1e-9, 1e-8, 1e-3])
    active_h = torch.diag(torch.tensor(diagonal, dtype=torch.float64))
    full_h = torch.eye(15, dtype=torch.float64)
    full_h[9:, 9:] = active_h
    source = active_h if partial else full_h
    outputs = []
    for cutoff in (0., 5., 7.):
        info = {}
        freq, modes = _frequencies_cm_and_modes(
            source.clone(), [1] * 5, COORDS, torch.device('cpu'),
            freeze_idx=[0, 1, 2], projection_info=info,
            frequency_zero_cutoff_cm=cutoff,
        )
        assert freq.shape == (6,) and modes.shape == (6, 15)
        assert info['frequency_representation'] == 'complete'
        assert info['effective_rank'] == 0 and info['raw_mode_count'] == 6
        assert _strict_negative_count(freq, info) == 2
        assert np.count_nonzero(resolved_imaginary_mask(freq, cutoff)) == (2 if cutoff == 0. else 1)
        assert freq[1] < 0. < freq[3] < freq[4] < 5.
        assert freq[2] == 0.
        assert torch.count_nonzero(modes[:, :9]) == 0
        torch.testing.assert_close(modes @ modes.T, torch.eye(6, dtype=torch.float64))
        # Equal masses make the unweighted eigensystem an independent residual.
        torch.testing.assert_close(active_h @ modes[:, 9:].T,
                                   modes[:, 9:].T * torch.tensor(diagonal))
        outputs.append((freq, modes))
    for freq, modes in outputs[1:]:
        np.testing.assert_array_equal(freq, outputs[0][0])
        torch.testing.assert_close(modes, outputs[0][1], rtol=0., atol=0.)


@pytest.mark.parametrize('cutoff', [0., 5., 7.])
def test_complete_partition_counts_each_negative_once(cutoff):
    values = np.array([-20., -5., -.2, -0., 0., .2, 5., 20.])
    info = frequency_partition_info(values, cutoff)
    assert _strict_negative_count(values, info) == 3
    assert info['raw_mode_count'] == 8
    assert info['resolved_mode_count'] + info['near_zero_mode_count'] == 8
    for key, value in [('raw_mode_count', 9), ('near_zero_mode_count', True),
                       ('near_zero_frequencies_cm', []),
                       ('frequency_representation', 'unknown')]:
        damaged = copy.deepcopy(info)
        damaged[key] = value
        assert _strict_negative_count(values, damaged) is None


def test_tall_frozen_projection_avoids_quadratic_unused_allocation(monkeypatch):
    rng = np.random.default_rng(42)
    coords = torch.tensor(rng.normal(size=(602, 3)), dtype=torch.float64)
    real_svd = torch.linalg.svd
    observed = []

    def allocation_guard(matrix, *args, **kwargs):
        if matrix.shape[0] > 1000:
            assert not kwargs.get('full_matrices', True), 'unused quadratic U allocation'
            observed.append(tuple(matrix.shape))
        return real_svd(matrix, *args, **kwargs)

    monkeypatch.setattr(torch.linalg, 'svd', allocation_guard)
    basis, info = active_tr_basis(coords, torch.ones(602, dtype=torch.float64), [600, 601])
    assert (1800, 6) in observed  # frozen rows; global basis SVD is also reduced
    assert basis.shape == (6, 0) and info.effective_rank == 0


def test_soft_positive_modes_reach_thermochemistry_without_cutoff_dependence():
    from thermoanalysis.QCData import QCData
    from thermoanalysis.config import WORKFLOW_THERMO_POLICY
    from thermoanalysis.thermo import thermochemistry

    spectrum = np.array([-20., -.2, 0., .2, 3.1595608552, 20.])
    results = []
    for cutoff in (0., 5., 7.):
        frequency_partition_info(spectrum, cutoff)
        qc = QCData({'coords3d': COORDS, 'masses': np.ones(5),
                     'wavenumbers': spectrum.copy(), 'scf_energy': -1., 'mult': 1},
                    point_group='c1', mult=1)
        tr = thermochemistry(qc, 298.15, pressure=101325.,
                             **WORKFLOW_THERMO_POLICY.thermochemistry_kwargs())
        np.testing.assert_array_equal(tr.wavenumbers, spectrum[spectrum > 0.])
        assert np.isfinite([tr.ZPE, tr.G, tr.S_tot, tr.c_tot]).all()
        results.append([tr.ZPE, tr.G, tr.S_tot, tr.c_tot])
    np.testing.assert_array_equal(results, [results[0]] * 3)
