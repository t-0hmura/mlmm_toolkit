"""Focused regressions for bundled pysisyphus correctness fixes."""

from __future__ import annotations

import json
import inspect
from types import SimpleNamespace
import weakref

import numpy as np
import pytest
import torch

from pysisyphus.intcoords.BondedFragment import BondedFragment
from pysisyphus.intcoords.Cartesian import CartesianX, CartesianY
from pysisyphus.intcoords.LinearDisplacement import LinearDisplacement
from pysisyphus.intcoords.RedundantCoords import RedundantCoords
from pysisyphus.intcoords.Torsion import Torsion
from pysisyphus.intcoords.exceptions import NeedNewInternalsException
from pysisyphus.intcoords.update import transform_int_step, update_internals
from pysisyphus.io.qcschema import geom_from_qcschema
from pysisyphus.Geometry import Geometry
from pysisyphus.irc.EulerPC import EulerPC
from pysisyphus.irc.DWI import DWI
from pysisyphus.irc.IRC import IRC
from pysisyphus.irc.Instanton import Instanton
from pysisyphus.linalg import quaternion_to_rot_mat
from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.cos.GrowingChainOfStates import GrowingChainOfStates
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.hessian_updates import bofill_update
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from thermoanalysis.constants import C
from thermoanalysis.thermo import qrrho_vibrational_part_func


def _qcschema_payload() -> dict:
    return {
        "molecule": {
            "symbols": ["H", "H"],
            "geometry": [0.0, 0.0, 0.0, 0.0, 0.0, 1.4],
        },
        "comment": "x" * 300,
    }


def test_geometry_hessian_ownership_transfer_keeps_other_results() -> None:
    geom = Geometry(("H",), np.zeros(3), coord_type="cart")
    hessian = torch.eye(3, dtype=torch.float64)
    geom.cart_hessian = hessian
    geom.results = {"energy": -1.0, "hessian": hessian}

    taken = geom.take_hessian()

    assert taken is hessian
    assert geom._hessian is None
    assert geom.results == {"energy": -1.0}


def test_generic_hessian_recalculation_releases_all_old_dense_aliases() -> None:
    old_hessian = np.eye(3)
    old_ref = weakref.ref(old_hessian)

    class ReplacementGeometry:
        @property
        def hessian(self):
            assert old_ref() is None
            return 2.0 * np.eye(3)

    optimizer = SimpleNamespace(
        forces=[np.ones(3)],
        adapt_norm=None,
        hessian_recalc_adapt=None,
        hessian_recalc_in=0,
        hessian_recalc=5,
        H=old_hessian,
        cur_H=old_hessian,
        hessian_xtb=False,
        geometry=ReplacementGeometry(),
        using_active_dofs=False,
        cur_cycle=1,
        log=lambda _message: None,
    )
    del old_hessian

    HessianOptimizer.update_hessian(optimizer)

    np.testing.assert_allclose(optimizer.H, 2.0 * np.eye(3))
    assert optimizer.cur_H is None


def test_qcschema_accepts_json_text_and_file_paths(tmp_path) -> None:
    payload = _qcschema_payload()
    text = json.dumps(payload)
    from_text = geom_from_qcschema(text)
    path = tmp_path / "molecule.json"
    path.write_text(text)
    from_path = geom_from_qcschema(path)

    assert from_text.atoms == from_path.atoms == ("H", "H")
    np.testing.assert_allclose(from_text.coords, from_path.coords)


def test_bonded_fragment_gradient_is_local_to_bond_endpoints() -> None:
    coords = np.array([
        [4.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ])
    value, gradient = BondedFragment._calculate(
        coords, (0, 1), gradient=True, bond_indices=(1, 2),
    )

    assert value == pytest.approx(2.0)
    np.testing.assert_allclose(
        gradient.reshape(-1, 3),
        [[0.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0]],
    )


@pytest.mark.parametrize("z", [-0.5, 0.0, 0.5])
def test_torsion_jacobian_matches_finite_difference(z) -> None:
    coords = np.array([
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, z],
    ])
    indices = [0, 1, 2, 3]
    analytic = Torsion._jacobian(coords, indices).reshape(coords.size, coords.size)
    numeric = np.zeros_like(analytic)
    step = 1.0e-5
    for column in range(coords.size):
        plus = coords.copy().reshape(-1)
        minus = coords.copy().reshape(-1)
        plus[column] += step
        minus[column] -= step
        plus_gradient = Torsion._calculate(
            plus.reshape(-1, 3), indices, gradient=True,
        )[1]
        minus_gradient = Torsion._calculate(
            minus.reshape(-1, 3), indices, gradient=True,
        )[1]
        numeric[:, column] = (plus_gradient - minus_gradient) / (2 * step)
    np.testing.assert_allclose(analytic, numeric, atol=2.0e-6, rtol=2.0e-6)


def test_complement_linear_displacement_derivatives_match_scalar() -> None:
    coords = np.array(
        [
            [-1.1, 0.2, 0.1],
            [0.0, 0.0, 0.0],
            [1.3, 0.4, -0.2],
        ],
        dtype=float,
    )
    indices = [0, 1, 2]
    cross_vec = LinearDisplacement._get_cross_vec(coords, indices)
    _, analytic_gradient = LinearDisplacement._calculate(
        coords,
        indices,
        gradient=True,
        complement=True,
        cross_vec=cross_vec.copy(),
    )
    numeric_gradient = np.zeros(9, dtype=float)
    scalar_step = 1.0e-5
    for dof in range(coords.size):
        plus = coords.copy().reshape(-1)
        minus = coords.copy().reshape(-1)
        plus[dof] += scalar_step
        minus[dof] -= scalar_step
        value_plus = LinearDisplacement._calculate(
            plus.reshape(-1, 3),
            indices,
            complement=True,
            cross_vec=cross_vec.copy(),
        )
        value_minus = LinearDisplacement._calculate(
            minus.reshape(-1, 3),
            indices,
            complement=True,
            cross_vec=cross_vec.copy(),
        )
        numeric_gradient[dof] = (value_plus - value_minus) / (2 * scalar_step)
    np.testing.assert_allclose(
        analytic_gradient, numeric_gradient, atol=2.0e-8, rtol=2.0e-8,
    )

    analytic = LinearDisplacement._jacobian(
        coords,
        indices,
        complement=True,
        cross_vec=cross_vec.copy(),
    ).reshape(9, 9)
    numeric = np.zeros((9, 9), dtype=float)
    step = 1.0e-5
    for column in range(coords.size):
        plus = coords.copy().reshape(-1)
        minus = coords.copy().reshape(-1)
        plus[column] += step
        minus[column] -= step
        plus_gradient = LinearDisplacement._calculate(
            plus.reshape(-1, 3),
            indices,
            gradient=True,
            complement=True,
            cross_vec=cross_vec.copy(),
        )[1]
        minus_gradient = LinearDisplacement._calculate(
            minus.reshape(-1, 3),
            indices,
            gradient=True,
            complement=True,
            cross_vec=cross_vec.copy(),
        )[1]
        numeric[:, column] = (plus_gradient - minus_gradient) / (2 * step)
    np.testing.assert_allclose(analytic, numeric, atol=3.0e-6, rtol=3.0e-6)


def test_invalid_dihedral_and_bend_indices_are_unioned(monkeypatch) -> None:
    class FakePrimitive:
        def __init__(self, indices):
            self.indices = list(indices)

        def calculate(self, coords3d, gradient=False):
            return 0.0, np.zeros(coords3d.size)

    primitives = [
        FakePrimitive((0, 1, 2)),
        FakePrimitive((0,)),
        FakePrimitive((0, 1, 2, 3)),
    ]
    monkeypatch.setattr(
        "pysisyphus.intcoords.update.dihedral_valid", lambda *_: False,
    )
    monkeypatch.setattr(
        "pysisyphus.intcoords.update.bend_valid", lambda *_: True,
    )

    with pytest.raises(NeedNewInternalsException) as caught:
        update_internals(
            np.zeros((4, 3)),
            np.zeros(3),
            primitives,
            dihedral_inds=[2],
            rotation_inds=[],
            bend_inds=[0],
            check_dihedrals=True,
            check_bends=True,
        )
    assert caught.value.invalid_inds == [2]


def test_all_constraints_must_be_satisfied_after_backtransform() -> None:
    primitives = [CartesianX([0]), CartesianY([0])]
    _, cart_step, failed = transform_int_step(
        int_step=np.array([0.0, 1.0]),
        old_cart_coords=np.zeros(3),
        cur_internals=np.zeros(2),
        Bt_inv_prim=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        primitives=primitives,
        dihedral_inds=[],
        rotation_inds=[],
        bend_inds=[],
        constrained_inds=[0, 1],
        update_constraints=False,
    )
    assert not failed
    np.testing.assert_allclose(cart_step, 0.0, atol=1.0e-12)


def test_redundant_coords_string_uses_existing_index_properties() -> None:
    coords = RedundantCoords.__new__(RedundantCoords)
    coords._bond_inds = [0]
    coords._bend_inds = [1, 2]
    coords._dihedral_inds = [3]
    assert str(coords) == "RedundantCoords(1 bonds, 2 bends, 1 dihedrals)"


def test_instanton_hessian_and_analytical_flag_follow_images() -> None:
    instanton = Instanton.__new__(Instanton)
    instanton.images = [
        SimpleNamespace(cart_hessian=np.eye(3), is_analytical_2d=True),
        SimpleNamespace(cart_hessian=2.0 * np.eye(3), is_analytical_2d=True),
    ]
    np.testing.assert_allclose(
        instanton.cart_hessian,
        np.diag([1.0, 1.0, 1.0, 2.0, 2.0, 2.0]),
    )
    assert instanton.is_analytical_2d is True


def test_quaternion_rotation_matrix_is_orthogonal() -> None:
    quaternion = np.array([0.5, -0.5, 0.5, 0.5])
    rotation = quaternion_to_rot_mat(quaternion)
    np.testing.assert_allclose(rotation.T @ rotation, np.eye(3), atol=1.0e-12)
    assert np.linalg.det(rotation) == pytest.approx(1.0)


def test_qrrho_partition_function_uses_frequency_in_effective_inertia() -> None:
    frequencies = np.array([20.0, 50.0, 1600.0]) * 100.0 * C
    q_vib, q_vib_v0 = qrrho_vibrational_part_func(
        298.15, frequencies, I_mean=10.0, cutoff=100.0, alpha=4,
    )
    assert q_vib == pytest.approx(0.04492554859960663, rel=1.0e-13)
    assert q_vib_v0 == pytest.approx(2.148803489705469, rel=1.0e-13)


@pytest.mark.parametrize(
    ("kind", "eigvals", "gradient"),
    [
        ("max", np.array([-0.2]), np.array([0.1])),
        ("min", np.array([0.3]), np.array([0.12])),
    ],
)
def test_rsprfo_partition_derivative_matches_secular_finite_difference(
    kind, eigvals, gradient,
) -> None:
    opt = RSPRFOptimizer.__new__(RSPRFOptimizer)
    opt.rfo_dict = {"max": (None, "max"), "min": (None, "min")}
    opt.log = lambda *_: None
    alpha = 10.0
    eps = 1.0e-5

    minus = opt.solve_rfo_secular(
        eigvals, gradient, alpha - eps, kind=kind,
    )
    center = opt.solve_rfo_secular(eigvals, gradient, alpha, kind=kind)
    plus = opt.solve_rfo_secular(
        eigvals, gradient, alpha + eps, kind=kind,
    )
    assert minus is not None and center is not None and plus is not None
    numeric = (
        np.dot(plus[0], plus[0]) - np.dot(minus[0], minus[0])
    ) / (2.0 * eps)
    analytic = opt._partition_dstep2_dalpha(
        alpha, center[1], center[0], eigvals, gradient,
    )
    assert analytic == pytest.approx(numeric, rel=2.0e-7, abs=1.0e-10)


def test_rsprfo_fallback_enforces_full_trust_radius() -> None:
    opt = RSPRFOptimizer.__new__(RSPRFOptimizer)
    opt.trust_radius = 0.1
    opt.log = lambda *_: None
    step = opt._restrict_final_step(np.array([0.12, 0.16]))
    assert np.linalg.norm(step) == pytest.approx(0.1)


def test_rfoptimizer_rejects_oversized_accelerated_displacement() -> None:
    opt = RFOptimizer.__new__(RFOptimizer)
    opt.trust_radius = 0.1
    opt.log = lambda *_: None
    ref_step = np.array([0.08, 0.0])
    result = opt._accept_accelerated_step(
        np.array([0.08, 0.0]), np.array([0.08, 0.0]), ref_step,
    )
    np.testing.assert_allclose(result, ref_step)

    accepted = opt._accept_accelerated_step(
        np.array([0.04, 0.0]), np.array([0.03, 0.0]), ref_step,
    )
    np.testing.assert_allclose(accepted, [0.07, 0.0])


def test_rfoptimizer_accelerated_step_keeps_reference_tensor_representation() -> None:
    opt = RFOptimizer.__new__(RFOptimizer)
    opt.trust_radius = 0.1
    opt.log = lambda *_: None
    ref_step = torch.tensor([0.08, 0.0], dtype=torch.float64)

    accepted = opt._accept_accelerated_step(
        np.array([0.04, 0.0]), torch.tensor([0.03, 0.0]), ref_step,
    )

    assert isinstance(accepted, torch.Tensor)
    assert accepted.dtype == ref_step.dtype
    assert accepted.device == ref_step.device
    torch.testing.assert_close(accepted, torch.tensor([0.07, 0.0], dtype=torch.float64))


def test_full_string_budget_requests_nonconverged_stop(monkeypatch) -> None:
    monkeypatch.setattr(
        Optimizer, "check_convergence", lambda *_args, **_kwargs: (False, "no"),
    )
    opt = StringOptimizer.__new__(StringOptimizer)
    opt.geometry = SimpleNamespace(fully_grown=True)
    opt.stop_in = 1
    opt.stop_in_when_full = 1
    opt.stop_requested = False
    opt.stop_reason = ""
    opt.log = lambda *_: None

    converged, _ = opt.check_convergence()
    assert converged is False
    assert opt.stop_requested is True
    assert opt.stop_reason == "full-string cycle budget exhausted"


def test_automatic_climbing_never_selects_fixed_endpoint() -> None:
    cos = SimpleNamespace(
        get_hei_index=lambda: 2,
        moving_indices=np.array([1]),
        climb="one",
        started_climbing=True,
        fixed_climb_indices=None,
        log=lambda *_: None,
    )
    assert ChainOfStates.get_climbing_indices(cos) == ()


def test_growing_string_reparametrization_guards_zero_density() -> None:
    assert GrowingString._reparam_step_fraction(0.0, 0.0, 1.0e-3) is None
    with pytest.raises(ValueError, match="coincident parameter densities"):
        GrowingString._reparam_step_fraction(0.1, 0.0, 1.0e-3)


def test_growing_string_rejects_fewer_than_two_nodes() -> None:
    with pytest.raises(ValueError, match="at least 2"):
        GrowingChainOfStates([], lambda: None, max_nodes=1)


def test_optimizer_convergence_vector_excludes_frozen_cartesian_dofs() -> None:
    geometry = SimpleNamespace(
        coord_type="cart",
        active_dof_indices=np.array([0, 1, 2]),
        cart_coords=np.zeros(12),
    )
    opt = RFOptimizer.__new__(RFOptimizer)
    opt.is_cos = False
    opt.geometry = geometry
    active = opt._active_convergence_vector(
        np.array([4.0e-4] * 3 + [0.0] * 9),
    )
    assert np.sqrt(np.mean(active**2)) == pytest.approx(4.0e-4)


def test_cos_convergence_vector_uses_moving_active_dofs_only() -> None:
    image = SimpleNamespace(
        coord_type="cart", active_dof_indices=np.array([0, 1, 2]),
    )
    opt = RFOptimizer.__new__(RFOptimizer)
    opt.is_cos = True
    opt.geometry = SimpleNamespace(
        coords_length=6,
        moving_indices=np.array([1]),
        images=[image, image, image],
    )
    vector = np.arange(18.0)
    np.testing.assert_allclose(
        opt._active_convergence_vector(vector), vector[[6, 7, 8]],
    )


def test_equal_energy_upwinding_tangent_falls_back_to_path_geometry() -> None:
    class Image:
        def __init__(self, coords):
            self.coords = np.asarray(coords, dtype=float)
            self.energy = 0.0

        def __sub__(self, other):
            return self.coords - other.coords

    cos = ChainOfStates.__new__(ChainOfStates)
    cos.images = [
        Image([0.0, 0.0]),
        Image([1.0, 0.5]),
        Image([2.0, 0.0]),
    ]
    cos.started_climbing_lanczos = False

    tangent = cos.get_tangent(1, kind="upwinding")

    np.testing.assert_allclose(tangent, [1.0, 0.0])


def test_cos_public_exports_are_bound() -> None:
    import pysisyphus.cos as cos

    assert cos.ChainOfStates is ChainOfStates
    assert cos.GrowingChainOfStates is GrowingChainOfStates
    assert cos.GrowingString.__name__ == "GrowingString"


def test_irc_rms_gradient_uses_integration_basis() -> None:
    irc = IRC.__new__(IRC)
    irc._act_dofs = np.array([0, 1, 2])
    full = np.array([1.2e-3] * 3 + [0.0] * 9)
    assert irc.active_rms_gradient(full) == pytest.approx(1.2e-3)


def test_irc_releases_dense_direction_state_before_next_branch(
    monkeypatch,
) -> None:
    irc = IRC.__new__(IRC)
    irc.mw_hessian = object()
    irc.dwi = SimpleNamespace(hessians=[object(), object()])
    emptied = []
    monkeypatch.setattr("torch.cuda.is_available", lambda: True)
    monkeypatch.setattr("torch.cuda.empty_cache", lambda: emptied.append(True))

    irc._release_direction_hessian_state()

    assert irc.mw_hessian is None
    assert irc.dwi is None
    assert emptied == [True]


def test_irc_mass_weights_hessian_without_mutating_seed() -> None:
    irc = IRC.__new__(IRC)
    irc.mm_inv2 = np.array([0.5, 1.5, 2.0])
    seed = torch.arange(9, dtype=torch.float64).reshape(3, 3)
    before = seed.clone()

    weighted = irc._mw_hessian_active(seed)

    expected = irc.mm_inv2[:, None] * before.numpy() * irc.mm_inv2[None, :]
    np.testing.assert_allclose(weighted.numpy(), expected)
    assert torch.equal(seed, before)
    assert weighted.data_ptr() != seed.data_ptr()


def test_irc_releases_initial_and_finished_hessian_owners(
    monkeypatch,
) -> None:
    monkeypatch.setattr("torch.cuda.is_available", lambda: False)
    irc = IRC.__new__(IRC)
    irc.init_hessian = torch.eye(6)
    irc._release_initial_hessian()
    assert irc.init_hessian is None
    assert irc.init_hessian_shape == (6, 6)

    irc.dwi = SimpleNamespace(hessians=[torch.eye(2), torch.eye(2)])
    irc.mw_hessian = torch.eye(2)
    irc.backward = False
    irc._release_finished_interpolation_state()
    assert irc.dwi is None
    assert irc.mw_hessian is None

    retained = torch.eye(2)
    irc.dwi = SimpleNamespace(hessians=[torch.eye(2), torch.eye(2)])
    irc.mw_hessian = retained
    irc.backward = True
    irc._release_finished_interpolation_state()
    assert irc.dwi is None
    assert irc.mw_hessian is retained


def test_euler_bofill_update_is_in_place_and_matches_dense_formula() -> None:
    irc = EulerPC.__new__(EulerPC)
    original = torch.tensor(
        [[2.0, 0.2, 0.1], [0.2, 1.5, -0.3], [0.1, -0.3, 1.1]],
        dtype=torch.float64,
    )
    dx = torch.tensor([0.4, -0.2, 0.3], dtype=torch.float64)
    dg = torch.tensor([0.1, 0.5, -0.4], dtype=torch.float64)
    dense_update, _ = bofill_update(original.clone(), dx, dg)
    expected = original + dense_update
    irc.mw_hessian = original.clone()
    irc.hessian_update_func = bofill_update
    storage = irc.mw_hessian.data_ptr()

    key = irc._apply_hessian_update(dx, dg)

    assert key == "Bofill"
    assert irc.mw_hessian.data_ptr() == storage
    torch.testing.assert_close(irc.mw_hessian, expected)


@pytest.mark.parametrize("roots", ([1], [0, 1]))
def test_ts_fixed_root_mode_preserves_configured_roots(roots) -> None:
    optimizer = RSPRFOptimizer.__new__(RSPRFOptimizer)
    optimizer.small_eigval_thresh = 1.0e-8
    optimizer.log_negative_eigenvalues = lambda *_args: None
    optimizer._physical_ts_mode = None
    optimizer.track_mode_by_overlap = False
    optimizer.roots = list(roots)

    eigvals = np.array([-2.0, -1.0, 0.5])
    eigvecs = np.eye(3)
    optimizer.update_ts_mode(eigvals, eigvecs)

    np.testing.assert_array_equal(optimizer.roots, roots)
    np.testing.assert_allclose(optimizer.ts_modes, eigvecs[:, roots].T)


@pytest.mark.parametrize("roots", ([0, 0], [3]))
def test_ts_fixed_root_mode_rejects_invalid_roots(roots) -> None:
    optimizer = RSPRFOptimizer.__new__(RSPRFOptimizer)
    optimizer.small_eigval_thresh = 1.0e-8
    optimizer.log_negative_eigenvalues = lambda *_args: None
    optimizer._physical_ts_mode = None
    optimizer.track_mode_by_overlap = False
    optimizer.roots = list(roots)

    with pytest.raises(ValueError):
        optimizer.update_ts_mode(np.array([-2.0, -1.0, 0.5]), np.eye(3))


def test_dwi_evicts_old_hessian_before_copying_replacement() -> None:
    dwi = DWI(maxlen=2)
    matrices = [torch.eye(3) * value for value in (1.0, 2.0, 3.0)]
    for value, matrix in enumerate(matrices):
        dwi.update(
            np.array([value]),
            float(value),
            np.array([value]),
            matrix,
            copy_hessian=True,
        )

    assert len(dwi.hessians) == 2
    torch.testing.assert_close(dwi.hessians[0], matrices[1])
    torch.testing.assert_close(dwi.hessians[1], matrices[2])
    assert dwi.hessians[1].data_ptr() != matrices[2].data_ptr()
    source = inspect.getsource(DWI.update)
    assert source.index("self.hessians.popleft()") < source.index(
        "hessian.detach().clone()"
    )


def test_euler_corrector_reaches_unweighted_target_for_heavy_mass() -> None:
    class ConstantDWI:
        @staticmethod
        def interpolate(coords, gradient=True):
            return 0.0, np.ones_like(coords)

    irc = EulerPC.__new__(EulerPC)
    irc._m_sqrt = np.array([4.0])  # oxygen-like mass sqrt
    irc._act_dofs = np.array([0])
    irc.log = lambda *_: None
    start = np.zeros(1)
    corrected = irc.corrector_step(start, 0.1, ConstantDWI())
    unweighted_length = np.linalg.norm((corrected - start) / irc._m_sqrt)
    assert unweighted_length == pytest.approx(0.1, abs=1.0e-4)


def test_euler_corrector_rejects_incomplete_zero_gradient() -> None:
    class ZeroDWI:
        @staticmethod
        def interpolate(coords, gradient=True):
            return 0.0, np.zeros_like(coords)

    irc = EulerPC.__new__(EulerPC)
    irc._m_sqrt = np.array([1.0])
    irc._act_dofs = np.array([0])
    irc.log = lambda *_: None
    with pytest.raises(RuntimeError, match="zero or non-finite gradient"):
        irc.corrector_step(np.zeros(1), 0.1, ZeroDWI())


def test_directional_irc_trajectory_contains_terminal_frame(tmp_path) -> None:
    irc = IRC.__new__(IRC)
    irc.atoms = ("H",)
    irc._m_sqrt = np.ones(3)
    irc.get_path_for_fn = lambda filename: str(tmp_path / filename)
    # Forward data has already been reversed by IRC.irc() for stitched-path
    # order: endpoint -> TS-adjacent.
    irc.irc_coords = [
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 0.0]),
    ]
    irc.irc_gradients = [np.zeros(3), np.zeros(3)]
    irc.irc_mw_coords = [coords.copy() for coords in irc.irc_coords]
    irc.irc_mw_gradients = [np.zeros(3), np.zeros(3)]
    irc.irc_energies = [-1.1, -1.0]
    irc.all_coords = []
    irc.all_gradients = []
    irc.all_mw_coords = []
    irc.all_mw_gradients = []
    irc.all_energies = []
    irc.converged = True
    irc.integration_stop_reason = ""
    irc.energy_increased = False
    irc.energy_converged = True
    irc.never_stop = False
    irc.cur_cycle = 1

    irc.set_data("forward")

    atom_lines = [
        line
        for line in (tmp_path / "forward_irc_trj.xyz").read_text().splitlines()
        if line.strip().startswith("H ")
    ]
    assert len(atom_lines) == len(irc.forward_energies) == 2
    # Directional file is chronological TS -> endpoint and includes endpoint.
    assert float(atom_lines[0].split()[1]) == pytest.approx(0.0)
    assert float(atom_lines[-1].split()[1]) == pytest.approx(
        0.529177, rel=1e-6
    )


def test_full_hessian_normal_modes_honor_geometry_freezes() -> None:
    geometry = Geometry(
        ["C"] * 5,
        np.array([
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
            [2.0, 2.0, 2.0],
        ]).reshape(-1),
        coord_type="cart",
        freeze_atoms=[0, 1, 2],
    )
    _, _, _, cart_modes = geometry.get_normal_modes(np.eye(15), full=True)

    assert cart_modes.shape == (15, 6)
    np.testing.assert_allclose(cart_modes[:9], 0.0)
    assert geometry.within_partial_hessian is None


def test_final_hessian_thermochemistry_uses_keyword_contract(monkeypatch) -> None:
    import pysisyphus.helpers as helpers

    calls = []

    class FakeGeometry:
        cart_hessian = np.eye(3)

        @staticmethod
        def mass_weigh_hessian(hessian):
            return hessian

        @staticmethod
        def eckart_projection(hessian):
            return hessian

        @staticmethod
        def get_thermoanalysis(*, T, p):
            calls.append((T, p))
            return "thermo"

    monkeypatch.setattr(helpers, "report_isotopes", lambda *_args: None)
    result = helpers.do_final_hessian(FakeGeometry(), T=310.0, p=98_000.0)

    assert calls == [(310.0, 98_000.0)]
    assert result.thermo == "thermo"


@pytest.mark.parametrize(
    "update_name",
    ["ts_bfgs_update", "ts_bfgs_update_org", "ts_bfgs_update_revised"],
)
@pytest.mark.parametrize("device", ["cpu", "cuda"])
def test_ts_bfgs_updates_preserve_tensor_device_and_numpy_values(
    update_name: str, device: str
) -> None:
    if device == "cuda" and not torch.cuda.is_available():
        pytest.skip("CUDA is unavailable")
    from pysisyphus.optimizers import hessian_updates

    update = getattr(hessian_updates, update_name)
    hessian = np.array(
        [[-1.0, 0.1, 0.0], [0.1, 2.0, 0.2], [0.0, 0.2, 3.0]],
        dtype=np.float64,
    )
    step = np.array([0.2, -0.1, 0.3], dtype=np.float64)
    gradient_delta = hessian @ step + np.array([0.1, 0.05, -0.02])
    expected, _ = update(hessian, step, gradient_delta)

    actual, _ = update(
        torch.as_tensor(hessian, device=device),
        torch.as_tensor(step, device=device),
        torch.as_tensor(gradient_delta, device=device),
    )

    assert isinstance(actual, torch.Tensor)
    assert actual.device.type == device
    assert actual.dtype == torch.float64
    np.testing.assert_allclose(actual.detach().cpu().numpy(), expected)
