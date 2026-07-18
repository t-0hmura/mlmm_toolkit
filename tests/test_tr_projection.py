import pytest
import numpy as np
import torch

from pysisyphus.Geometry import Geometry
from pysisyphus.constants import AMU2AU
from pysisyphus.tr_projection import active_tr_basis, project_hessian_inplace


COORDS = torch.tensor(
    [
        [0.0, 0.0, 0.0],
        [1.2, 0.0, 0.0],
        [0.1, 1.1, 0.0],
        [0.2, 0.3, 1.3],
        [1.0, 0.8, 0.7],
    ],
    dtype=torch.float64,
)
MASSES = torch.tensor([12.0, 1.0, 16.0, 14.0, 32.0], dtype=torch.float64)


@pytest.mark.parametrize(
    ("frozen", "expected_rank"),
    [([], 6), ([0], 3), ([0, 1], 1), ([0, 1, 2], 0)],
)
def test_constrained_rigid_null_rank(frozen, expected_rank):
    active = [i for i in range(len(COORDS)) if i not in frozen]
    basis, info = active_tr_basis(COORDS, MASSES, active, mode="constrained")

    assert info.effective_rank == expected_rank
    assert info.full_rigid_rank == 6
    assert basis.shape == (3 * len(active), expected_rank)
    if expected_rank:
        torch.testing.assert_close(
            basis.T @ basis,
            torch.eye(expected_rank, dtype=basis.dtype),
            atol=1.0e-12,
            rtol=1.0e-12,
        )


def test_legacy_active_uses_isolated_fragment_projection():
    active = [2, 3, 4]
    constrained, constrained_info = active_tr_basis(
        COORDS, MASSES, active, mode="constrained"
    )
    legacy, legacy_info = active_tr_basis(
        COORDS, MASSES, active, mode="legacy-active"
    )

    assert constrained_info.effective_rank == 1
    assert legacy_info.effective_rank == 6
    assert constrained.shape[1] < legacy.shape[1]


def test_no_freeze_modes_have_the_same_projector():
    active = list(range(len(COORDS)))
    constrained, _ = active_tr_basis(COORDS, MASSES, active, mode="constrained")
    legacy, _ = active_tr_basis(COORDS, MASSES, active, mode="legacy-active")
    torch.testing.assert_close(
        constrained @ constrained.T,
        legacy @ legacy.T,
        atol=1.0e-12,
        rtol=1.0e-12,
    )


def test_linear_system_has_five_full_rigid_modes():
    coords = torch.tensor(
        [[-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        dtype=torch.float64,
    )
    masses = torch.ones(3, dtype=torch.float64)
    basis, info = active_tr_basis(coords, masses, [0, 1, 2])
    assert info.full_rigid_rank == 5
    assert info.effective_rank == 5
    assert basis.shape == (9, 5)


def test_rank_zero_projection_is_an_exact_noop():
    active = [3, 4]
    basis, info = active_tr_basis(COORDS, MASSES, active, mode="constrained")
    assert info.effective_rank == 0
    hessian = torch.arange(36, dtype=torch.float64).reshape(6, 6)
    original = hessian.clone()
    returned = project_hessian_inplace(hessian, basis)
    assert returned.data_ptr() == hessian.data_ptr()
    torch.testing.assert_close(hessian, original)


def test_full_and_active_hessian_phva_are_equivalent():
    from mlmm.workflows.freq import _frequencies_cm_and_modes

    generator = torch.Generator().manual_seed(41)
    raw = torch.randn((15, 15), generator=generator, dtype=torch.float64)
    full_hessian = raw.T @ raw + 0.2 * torch.eye(15, dtype=torch.float64)
    active_atoms = [2, 3, 4]
    active_dofs = [3 * atom + axis for atom in active_atoms for axis in range(3)]
    active_hessian = full_hessian[active_dofs][:, active_dofs].clone()
    atomic_numbers = [6, 1, 8, 7, 16]
    device = torch.device("cpu")
    full_info = {}
    active_info = {}

    full_freqs, _ = _frequencies_cm_and_modes(
        full_hessian.clone(),
        atomic_numbers,
        COORDS.numpy(),
        device,
        freeze_idx=[0, 1],
        projection_info=full_info,
    )
    active_freqs, _ = _frequencies_cm_and_modes(
        active_hessian,
        atomic_numbers,
        COORDS.numpy(),
        device,
        freeze_idx=[0, 1],
        projection_info=active_info,
    )

    np.testing.assert_allclose(full_freqs, active_freqs, rtol=1.0e-10, atol=1.0e-8)
    assert full_info["treatment"] == active_info["treatment"] == "constrained"
    assert full_info["effective_rank"] == active_info["effective_rank"] == 1


def test_constrained_rank_zero_does_not_drop_active_modes():
    from mlmm.workflows.freq import _frequencies_cm_and_modes

    atomic_numbers = [6, 1, 8, 7, 16]
    active_hessian = torch.eye(6, dtype=torch.float64)
    constrained, _ = _frequencies_cm_and_modes(
        active_hessian.clone(),
        atomic_numbers,
        COORDS.numpy(),
        torch.device("cpu"),
        freeze_idx=[0, 1, 2],
        tr_projection="constrained",
    )
    legacy, _ = _frequencies_cm_and_modes(
        active_hessian.clone(),
        atomic_numbers,
        COORDS.numpy(),
        torch.device("cpu"),
        freeze_idx=[0, 1, 2],
        tr_projection="legacy-active",
    )

    assert len(constrained) == 6
    assert len(legacy) < len(constrained)


@pytest.mark.parametrize("dtype", [torch.float32, torch.float64])
@pytest.mark.parametrize(
    ("frozen", "expected_rank"),
    [([], 6), ([0], 3), ([0, 1], 1), ([0, 1, 2], 0)],
)
def test_ranks_do_not_depend_on_floating_point_dtype(dtype, frozen, expected_rank):
    active = [i for i in range(len(COORDS)) if i not in frozen]
    _, info = active_tr_basis(
        COORDS.to(dtype=dtype), MASSES.to(dtype=dtype), active
    )
    assert info.effective_rank == expected_rank


def test_active_atom_order_is_preserved():
    sorted_active = [2, 3, 4]
    permuted_active = [4, 2, 3]
    sorted_basis, _ = active_tr_basis(COORDS, MASSES, sorted_active)
    permuted_basis, info = active_tr_basis(COORDS, MASSES, permuted_active)
    row_order = [
        3 * sorted_active.index(atom) + axis
        for atom in permuted_active
        for axis in range(3)
    ]
    expected = (sorted_basis @ sorted_basis.T)[row_order][:, row_order]
    torch.testing.assert_close(
        permuted_basis @ permuted_basis.T,
        expected,
        atol=1.0e-12,
        rtol=1.0e-12,
    )
    assert info.active_atoms == tuple(permuted_active)


def test_projection_is_correct_for_asymmetric_hessian():
    basis, _ = active_tr_basis(COORDS, MASSES, list(range(len(COORDS))))
    generator = torch.Generator().manual_seed(87)
    hessian = torch.randn((15, 15), generator=generator, dtype=torch.float64)
    projector = torch.eye(15, dtype=torch.float64) - basis @ basis.T
    expected = projector @ hessian @ projector
    project_hessian_inplace(hessian, basis)
    torch.testing.assert_close(hessian, expected, atol=1.0e-12, rtol=1.0e-12)


@pytest.mark.parametrize("frozen", [[0], [0, 1], [0, 1, 2]])
@pytest.mark.parametrize("mode", ["constrained", "legacy-active"])
@pytest.mark.parametrize("container", ["numpy", "torch"])
def test_geometry_partial_normal_modes_use_input_representation(
    frozen, mode, container
):
    atoms = ["C", "H", "O", "N", "S"]
    geom = Geometry(
        atoms,
        COORDS.numpy().reshape(-1),
        freeze_atoms=frozen,
        tr_projection=mode,
    )
    active = [i for i in range(len(atoms)) if i not in frozen]
    active_dofs = [3 * atom + axis for atom in active for axis in range(3)]
    geom.within_partial_hessian = {
        "active_n_dof": len(active_dofs),
        "full_n_dof": 3 * len(atoms),
        "active_dofs": active_dofs,
        "active_atoms": active,
    }
    hessian = np.eye(len(active_dofs), dtype=float)
    if container == "torch":
        hessian = torch.as_tensor(hessian, dtype=torch.float64)

    _, eigvals, mw_modes, cart_modes = geom.get_normal_modes(hessian)
    _, info = active_tr_basis(COORDS, torch.as_tensor(geom.masses), active, mode=mode)
    expected_modes = len(active_dofs) - info.effective_rank
    cart_array = (
        cart_modes.detach().cpu().numpy()
        if isinstance(cart_modes, torch.Tensor)
        else cart_modes
    )
    frozen_dofs = [3 * atom + axis for atom in frozen for axis in range(3)]

    assert len(eigvals) == expected_modes
    assert mw_modes.shape == (len(active_dofs), expected_modes)
    assert cart_modes.shape == (3 * len(atoms), expected_modes)
    assert np.all(np.isfinite(cart_array))
    assert np.allclose(cart_array[frozen_dofs], 0.0)
    np.testing.assert_allclose(np.linalg.norm(cart_array, axis=0), 1.0)


def test_geometry_full_projector_has_numpy_torch_parity_without_mutation():
    atoms = ["C", "H", "O", "N", "S"]
    geom = Geometry(atoms, COORDS.numpy().reshape(-1), freeze_atoms=[0])
    active = [1, 2, 3, 4]
    active_dofs = [3 * atom + axis for atom in active for axis in range(3)]
    geom.within_partial_hessian = {
        "active_n_dof": len(active_dofs),
        "full_n_dof": 15,
        "active_dofs": active_dofs,
        "active_atoms": active,
    }
    hessian_np = np.diag(np.linspace(1.0, 2.0, len(active_dofs)))
    hessian_torch = torch.as_tensor(hessian_np.copy(), dtype=torch.float64)
    original_np = hessian_np.copy()
    original_torch = hessian_torch.clone()

    projected_np, projector_np = geom.eckart_projection(
        hessian_np, return_P=True, full=True
    )
    projected_torch, projector_torch = geom.eckart_projection(
        hessian_torch, return_P=True, full=True
    )

    np.testing.assert_allclose(hessian_np, original_np)
    torch.testing.assert_close(hessian_torch, original_torch)
    np.testing.assert_allclose(projector_torch.numpy(), projector_np, atol=1.0e-12)
    np.testing.assert_allclose(projected_torch.numpy(), projected_np, atol=1.0e-12)
    np.testing.assert_allclose(projector_np, projector_np.T, atol=1.0e-12)
    np.testing.assert_allclose(projector_np @ projector_np, projector_np, atol=1.0e-12)


@pytest.mark.parametrize("frozen", [[0], [0, 1], [0, 1, 2]])
@pytest.mark.parametrize("mode", ["constrained", "legacy-active"])
@pytest.mark.parametrize("container", ["numpy", "torch"])
def test_geometry_full_and_active_hessians_have_identical_phva_spectra(
    frozen, mode, container
):
    atoms = ["C", "H", "O", "N", "S"]
    active = [atom for atom in range(len(atoms)) if atom not in frozen]
    active_dofs = [3 * atom + axis for atom in active for axis in range(3)]
    metadata = {
        "active_n_dof": len(active_dofs),
        "full_n_dof": 3 * len(atoms),
        "active_dofs": active_dofs,
        "active_atoms": active,
    }
    generator = np.random.default_rng(43)
    raw = generator.normal(size=(15, 15))
    full_hessian = raw.T @ raw + np.diag(np.linspace(0.2, 0.5, 15))
    active_hessian = full_hessian[np.ix_(active_dofs, active_dofs)].copy()
    if container == "torch":
        full_hessian = torch.as_tensor(full_hessian, dtype=torch.float64)
        active_hessian = torch.as_tensor(active_hessian, dtype=torch.float64)
    full_original = full_hessian.clone() if container == "torch" else full_hessian.copy()
    active_original = active_hessian.clone() if container == "torch" else active_hessian.copy()

    full_geom = Geometry(
        atoms, COORDS.numpy().reshape(-1), freeze_atoms=frozen, tr_projection=mode
    )
    active_geom = Geometry(
        atoms, COORDS.numpy().reshape(-1), freeze_atoms=frozen, tr_projection=mode
    )
    full_geom.within_partial_hessian = dict(metadata)
    active_geom.within_partial_hessian = dict(metadata)

    full_freqs, full_eigvals, full_mw_modes, full_cart_modes = full_geom.get_normal_modes(
        full_hessian
    )
    active_freqs, active_eigvals, active_mw_modes, active_cart_modes = active_geom.get_normal_modes(
        active_hessian
    )
    full_again = full_geom.get_normal_modes(full_hessian)
    active_again = active_geom.get_normal_modes(active_hessian)
    full_cart_modes = np.asarray(
        full_cart_modes.detach().cpu() if isinstance(full_cart_modes, torch.Tensor) else full_cart_modes
    )
    active_cart_modes = np.asarray(
        active_cart_modes.detach().cpu() if isinstance(active_cart_modes, torch.Tensor) else active_cart_modes
    )
    frozen_dofs = [3 * atom + axis for atom in frozen for axis in range(3)]

    np.testing.assert_allclose(full_eigvals, active_eigvals, rtol=1.0e-11, atol=1.0e-12)
    np.testing.assert_allclose(full_freqs, active_freqs, rtol=1.0e-11, atol=1.0e-8)
    np.testing.assert_allclose(full_freqs, full_again[0], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(active_freqs, active_again[0], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(full_eigvals, full_again[1], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(active_eigvals, active_again[1], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(np.abs(np.asarray(full_mw_modes)), np.abs(np.asarray(full_again[2])))
    np.testing.assert_allclose(np.abs(np.asarray(active_mw_modes)), np.abs(np.asarray(active_again[2])))
    np.testing.assert_allclose(full_cart_modes[frozen_dofs], 0.0, atol=0.0)
    np.testing.assert_allclose(active_cart_modes[frozen_dofs], 0.0, atol=0.0)
    if container == "torch":
        torch.testing.assert_close(full_hessian, full_original, rtol=0.0, atol=0.0)
        torch.testing.assert_close(active_hessian, active_original, rtol=0.0, atol=0.0)
    else:
        np.testing.assert_array_equal(full_hessian, full_original)
        np.testing.assert_array_equal(active_hessian, active_original)


@pytest.mark.parametrize("container", ["numpy", "torch"])
def test_geometry_cached_partial_hessian_survives_modes_and_thermochemistry(container):
    atoms = ["C", "H", "O", "N", "S"]
    active = [2, 3, 4]
    active_dofs = [3 * atom + axis for atom in active for axis in range(3)]
    geom = Geometry(atoms, COORDS.numpy().reshape(-1), freeze_atoms=[0, 1])
    geom.within_partial_hessian = {
        "active_n_dof": len(active_dofs),
        "full_n_dof": 15,
        "active_dofs": active_dofs,
        "active_atoms": active,
    }
    hessian = np.diag(np.linspace(1.0, 2.0, len(active_dofs)))
    if container == "torch":
        hessian = torch.as_tensor(hessian, dtype=torch.float64)
    geom.cart_hessian = hessian
    geom.energy = 0.0
    original = hessian.clone() if container == "torch" else hessian.copy()

    first = geom.get_normal_modes()
    second = geom.get_normal_modes()
    geom.get_thermoanalysis()

    np.testing.assert_allclose(first[0], second[0], rtol=0.0, atol=0.0)
    if container == "torch":
        torch.testing.assert_close(geom.cart_hessian, original, rtol=0.0, atol=0.0)
    else:
        np.testing.assert_array_equal(geom.cart_hessian, original)


def test_geometry_rejects_partial_hessian_dofs_in_a_different_order():
    geom = Geometry(["C", "H", "O"], COORDS[:3].numpy().reshape(-1))
    geom.within_partial_hessian = {
        "active_n_dof": 6,
        "full_n_dof": 9,
        "active_atoms": [1, 2],
        "active_dofs": [6, 7, 8, 3, 4, 5],
    }

    with pytest.raises(ValueError, match="complete xyz triplets"):
        geom.get_normal_modes(np.eye(9))


def test_geometry_preexisting_positional_signature_remains_valid():
    geom = Geometry(
        ["H"],
        np.zeros(3),
        None,
        "cart",
        None,
        None,
        None,
        "old-comment",
        "old-name",
    )
    assert geom.comment == "old-comment"
    assert geom.name == "old-name"
    assert geom.tr_projection == "constrained"


def test_collinear_and_coincident_anchor_ranks():
    collinear = COORDS.clone()
    collinear[:3] = torch.tensor(
        [[-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        dtype=torch.float64,
    )
    _, collinear_info = active_tr_basis(collinear, MASSES, [3, 4])
    assert collinear_info.effective_rank == 1

    coincident = COORDS.clone()
    coincident[:3] = 0.0
    _, coincident_info = active_tr_basis(coincident, MASSES, [3, 4])
    assert coincident_info.effective_rank == 3


def test_invalid_rank_tolerance_is_rejected():
    with pytest.raises(ValueError, match="rtol"):
        active_tr_basis(COORDS, MASSES, range(len(COORDS)), rtol=0.0)


def test_all_frozen_phva_is_an_explicit_error():
    from mlmm.workflows.freq import _frequencies_cm_and_modes

    with pytest.raises(ValueError, match="at least one active atom"):
        _frequencies_cm_and_modes(
            torch.eye(15, dtype=torch.float64),
            [6, 1, 8, 7, 16],
            COORDS.numpy(),
            torch.device("cpu"),
            freeze_idx=[0, 1, 2, 3, 4],
        )


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA is unavailable")
def test_cuda_cpu_projector_parity_for_float32_two_anchor_case():
    active = [2, 3, 4]
    cpu_basis, cpu_info = active_tr_basis(
        COORDS.float(), MASSES.float(), active
    )
    gpu_basis, gpu_info = active_tr_basis(
        COORDS.float().cuda(), MASSES.float().cuda(), active
    )
    assert cpu_info.effective_rank == gpu_info.effective_rank == 1
    torch.testing.assert_close(
        cpu_basis @ cpu_basis.T,
        (gpu_basis @ gpu_basis.T).cpu(),
        atol=1.0e-12,
        rtol=1.0e-12,
    )


@pytest.mark.parametrize("frozen", [[], [0], [0, 1]])
@pytest.mark.parametrize("active_block", [False, True])
def test_dimer_root_selection_excludes_projector_zero_modes(frozen, active_block):
    from mlmm.workflows.tsopt import (
        _mode_direction_by_root,
        _mode_direction_by_root_from_Hact,
    )

    n_atoms = len(COORDS)
    active = [atom for atom in range(n_atoms) if atom not in frozen]
    masses_au = MASSES * AMU2AU
    full_hessian = torch.diag(torch.repeat_interleave(MASSES, 3))

    if active_block:
        active_dofs = [3 * atom + axis for atom in active for axis in range(3)]
        mode = _mode_direction_by_root_from_Hact(
            full_hessian[active_dofs][:, active_dofs],
            COORDS.numpy(),
            [6, 1, 8, 7, 16],
            masses_au,
            active,
            torch.device("cpu"),
        )
    else:
        mode = _mode_direction_by_root(
            full_hessian,
            COORDS,
            masses_au,
            freeze_idx=frozen,
        )

    mode_active = torch.as_tensor(mode[active].reshape(-1), dtype=torch.float64)
    mode_mw = mode_active * torch.sqrt(torch.repeat_interleave(MASSES[active], 3))
    mode_mw /= torch.linalg.norm(mode_mw)
    basis, info = active_tr_basis(COORDS, masses_au, active)

    assert info.effective_rank > 0
    assert np.all(np.isfinite(mode))
    torch.testing.assert_close(
        basis.T @ mode_mw,
        torch.zeros(info.effective_rank, dtype=torch.float64),
        atol=2.0e-10,
        rtol=0.0,
    )


@pytest.mark.parametrize("frozen", [[], [0], [0, 1]])
def test_dimer_root_order_is_preserved_in_compact_complement(frozen):
    from mlmm.workflows.tsopt import (
        _mode_direction_by_root,
        _mode_direction_by_root_from_Hact,
    )

    n_atoms = len(COORDS)
    active = [atom for atom in range(n_atoms) if atom not in frozen]
    active_dofs = [3 * atom + axis for atom in active for axis in range(3)]
    masses_au = MASSES * AMU2AU
    basis, info = active_tr_basis(COORDS, masses_au, active)
    complete, _ = torch.linalg.qr(basis, mode="complete")
    complement = complete[:, info.effective_rank :]
    spectrum = torch.arange(1, complement.shape[1] + 1, dtype=torch.float64)
    spectrum[:2] = torch.tensor([-2.0, -1.0], dtype=torch.float64)
    hessian_mw = complement @ torch.diag(spectrum) @ complement.T
    sqrt_m = torch.sqrt(torch.repeat_interleave(MASSES[active], 3))
    hessian_active = sqrt_m[:, None] * hessian_mw * sqrt_m[None, :]
    hessian_full = torch.diag(torch.repeat_interleave(MASSES, 3))
    hessian_full[np.ix_(active_dofs, active_dofs)] = hessian_active

    modes_full = []
    modes_active = []
    for root in (0, 1):
        mode_full = _mode_direction_by_root(
            hessian_full.clone(),
            COORDS,
            masses_au,
            root=root,
            freeze_idx=frozen,
        )
        mode_active = _mode_direction_by_root_from_Hact(
            hessian_active.clone(),
            COORDS.numpy(),
            [6, 1, 8, 7, 16],
            masses_au,
            active,
            torch.device("cpu"),
            root=root,
        )
        for mode, collected in ((mode_full, modes_full), (mode_active, modes_active)):
            mw = torch.as_tensor(mode[active].reshape(-1), dtype=torch.float64) * sqrt_m
            mw /= torch.linalg.norm(mw)
            collected.append(mw)
            assert abs(float(torch.dot(mw, complement[:, root]))) > 1.0 - 1.0e-10

    # root=0 uses iterative LOBPCG, whose returned orthogonality is typically
    # at the 1e-8 level for this deliberately degenerate positive complement.
    assert abs(float(torch.dot(modes_full[0], modes_full[1]))) < 1.0e-7
    for full, partial in zip(modes_full, modes_active):
        assert abs(float(torch.dot(full, partial))) > 1.0 - 1.0e-10
