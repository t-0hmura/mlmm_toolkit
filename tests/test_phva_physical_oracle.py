"""Independent distance-potential and geometric boundary oracles for both bundles."""
from itertools import combinations

import numpy as np
import pytest
import torch
from ase.data import atomic_masses
from scipy.constants import c, physical_constants

from pysisyphus.normal_modes import _frequencies_cm_and_modes
from pysisyphus.tr_projection import active_tr_basis


REFERENCE = np.array([[0., 0., 0.], [1.2, 0., 0.],
                      [0.1, 1.1, 0.], [0.2, 0.3, 1.3]])
NUMBERS = [6, 1, 8, 14]


def distance_model(coords, *, signed=False):
    """Exact derivatives of sum k_ij (r_ij-r0_ij)^2 / 2; no product helpers."""
    energy, gradient, hessian = 0., np.zeros(12), np.zeros((12, 12))
    for i, j in combinations(range(4), 2):
        delta = coords[i] - coords[j]
        length = np.linalg.norm(delta)
        reference = np.linalg.norm(REFERENCE[i] - REFERENCE[j])
        unit = delta / length
        stiffness = -0.3 if signed and (i, j) == (2, 3) else 0.2
        stretch = length - reference
        energy += 0.5 * stiffness * stretch ** 2
        pair_gradient = stiffness * stretch * unit
        pair_hessian = stiffness * (
            np.outer(unit, unit)
            + stretch / length * (np.eye(3) - np.outer(unit, unit))
        )
        ia, ja = slice(3*i, 3*i+3), slice(3*j, 3*j+3)
        gradient[ia] += pair_gradient
        gradient[ja] -= pair_gradient
        hessian[ia, ia] += pair_hessian
        hessian[ja, ja] += pair_hessian
        hessian[ia, ja] -= pair_hessian
        hessian[ja, ia] -= pair_hessian
    return energy, gradient, hessian


def independent_frequency_factor():
    hartree = physical_constants['Hartree energy'][0]
    bohr = physical_constants['Bohr radius'][0]
    amu = physical_constants['atomic mass constant'][0]
    return np.sqrt(hartree / (bohr ** 2 * amu)) / (2 * np.pi * c * 100)


@pytest.mark.parametrize('frozen, rigid_rank', [([], 6), ([0], 3), ([0, 1], 1), ([0, 1, 2], 0)])
@pytest.mark.parametrize('signed', [False, True])
def test_stationary_distance_potential_phva(frozen, rigid_rank, signed):
    energy, gradient, full_hessian = distance_model(REFERENCE, signed=signed)
    assert energy == 0.
    np.testing.assert_array_equal(gradient, 0.)
    active = [i for i in range(4) if i not in frozen]
    dofs = [3*i+j for i in active for j in range(3)]
    mass = np.repeat(atomic_masses[NUMBERS][active], 3)
    active_hessian = full_hessian[np.ix_(dofs, dofs)]
    weighted = active_hessian / np.sqrt(mass[:, None] * mass[None, :])
    eigenvalues = np.linalg.eigvalsh(weighted)
    # This exact stationary distance model has only the independently known
    # compatible rigid null space. All other roots are well separated from zero.
    zero = np.abs(eigenvalues) < 1.e-10
    assert np.count_nonzero(zero) == rigid_rank
    resolved = eigenvalues[~zero]
    assert np.min(np.abs(resolved)) > 1.e-4
    assert np.count_nonzero(resolved < 0) == int(signed)
    expected = np.sign(resolved) * np.sqrt(np.abs(resolved)) * independent_frequency_factor()
    for raw in (full_hessian, active_hessian):
        info = {}
        frequencies, modes = _frequencies_cm_and_modes(
            torch.tensor(raw, dtype=torch.float64), NUMBERS, REFERENCE,
            torch.device('cpu'), freeze_idx=frozen, projection_info=info,
        )
        assert info['effective_rank'] == rigid_rank
        # Independent CODATA and ASE constants differ slightly by release.
        np.testing.assert_allclose(frequencies, expected, rtol=2.e-7, atol=1.e-5)
        vectors = modes.numpy()[:, dofs].T
        np.testing.assert_allclose(vectors.T @ vectors, np.eye(len(resolved)), atol=1.e-11)
        np.testing.assert_allclose(weighted @ vectors, vectors * resolved,
                                   atol=1.e-11, rtol=1.e-9)
        if frozen:
            frozen_dofs = [3*i+j for i in frozen for j in range(3)]
            np.testing.assert_array_equal(modes.numpy()[:, frozen_dofs], 0.)


def test_nonstationary_rotation_has_gradient_term():
    displaced = REFERENCE.copy()
    displaced[3] += [0.08, -0.04, 0.05]
    _, gradient, hessian = distance_model(displaced)
    axis = np.array([0.2, -0.3, 0.7])
    rotation = np.cross(axis, displaced).ravel()
    rotated_gradient = np.cross(axis, gradient.reshape(-1, 3)).ravel()
    assert np.linalg.norm(rotated_gradient) > 1.e-3
    np.testing.assert_allclose(hessian @ rotation, rotated_gradient, atol=1.e-14)
    translation = np.tile([0.2, 0.3, -0.1], 4)
    np.testing.assert_allclose(hessian @ translation, 0., atol=1.e-14)


@pytest.mark.parametrize('dtype', [torch.float32, torch.float64])
def test_resolved_near_collinear_rigid_ranks(dtype):
    for height, full_rank, constrained_rank in [(0., 5, 1), (1.e-5, 6, 0)]:
        triangle = np.array([[-1., 0., 0.], [0., height, 0.], [1., 0., 0.]])
        _, info = active_tr_basis(torch.tensor(triangle, dtype=dtype),
                                  torch.ones(3), [0, 1, 2])
        assert info.full_rigid_rank == info.effective_rank == full_rank
        anchors = np.vstack((triangle, [0., 0., 1.], [0.3, 1., 0.2]))
        basis, info = active_tr_basis(torch.tensor(anchors, dtype=dtype),
                                      torch.ones(5), [3, 4])
        assert info.effective_rank == constrained_rank
        if constrained_rank:
            # The only compatible motion is rotation around the anchor x-axis.
            expected = np.cross([1., 0., 0.], anchors[3:]).ravel()
            expected /= np.linalg.norm(expected)
            np.testing.assert_allclose(basis.numpy() @ basis.numpy().T,
                                       np.outer(expected, expected), atol=1.e-8)
