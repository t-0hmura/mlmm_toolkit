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
import pysisyphus.intcoords.augment_bonds as augment_bonds_module
from pysisyphus.intcoords.DummyTorsion import DummyTorsion
from pysisyphus.intcoords.PrimTypes import PrimTypes
from pysisyphus.intcoords.exceptions import (
    NeedNewInternalsException,
    PrimitiveNotDefinedException,
)
from pysisyphus.intcoords.setup_fast import find_bonds, find_bonds_for_geom
from pysisyphus.intcoords.update import transform_int_step, update_internals
from pysisyphus.io.qcschema import geom_from_qcschema
from pysisyphus.Geometry import Geometry, get_trans_rot_vectors
from pysisyphus.calculators.Dimer import Dimer
import pysisyphus.normal_modes as nm
from pysisyphus.irc.EulerPC import EulerPC
from pysisyphus.irc.DWI import DWI
from pysisyphus.irc.IRC import IRC
from pysisyphus.irc.Instanton import Instanton
from pysisyphus.helpers import geom_loader
from pysisyphus.linalg import quaternion_to_rot_mat
from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.cos.GrowingChainOfStates import GrowingChainOfStates
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.closures import bfgs_multiply
import pysisyphus.optimizers.gdiis as gdiis_module
from pysisyphus.optimizers.hessian_updates import (
    bofill_update,
    damped_bfgs_update,
    flowchart_update,
)
from pysisyphus.optimizers.HessianOptimizer import dummy_hessian_update
from pysisyphus.optimizers.poly_fit import quartic_fit
from pysisyphus.modefollow.lanczos import lanczos
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from pysisyphus.tsoptimizers.TSHessianOptimizer import TSHessianOptimizer
from thermoanalysis.constants import AMU2KG, C, KB, PLANCK, R
from thermoanalysis.thermo import (
    chai_head_gordon_weights,
    qrrho_vibrational_part_func,
    vibrational_heat_capacity,
    vibrational_part_funcs,
)


def test_dimer_bias_force_norm_axis_error_keeps_translation_step(
    monkeypatch,
) -> None:
    dimer = object.__new__(Dimer)
    dimer.freeze_atoms = np.empty(0, dtype=int)
    dimer.rigid_basis = None
    dimer.rigid_basis_getter = None
    dimer.rotation_remove_trans = False
    dimer._coords0 = None
    dimer._energy0 = None
    dimer._f0 = None
    dimer._f1 = None
    dimer.force_evals = 0
    dimer.N = np.array([1.0, 0.0, 0.0])
    dimer.rotation_disable = True
    dimer.rotation_disable_pos_curv = True
    dimer.write_orientations = False
    dimer.length = 0.1
    dimer.bias_translation = True
    dimer.curvature = lambda *_args: 1.0
    dimer.calculator = SimpleNamespace(
        get_energy=lambda _atoms, _coords: {"energy": 0.0},
        get_forces=lambda _atoms, _coords: {
            "energy": 0.0,
            "forces": np.array([1.0, 0.0, 0.0]),
        },
    )
    dimer.gaussians = []
    dimer.get_gaussian_energies = lambda _coords: 0.0
    dimer.get_gaussian_forces = lambda _coords, sum_: np.zeros(3)
    dimer.trans_force_f_perp = True
    dimer.calc_counter = 0
    messages = []
    dimer.log = lambda message="": messages.append(message)
    dimer.make_fn = lambda name: name
    monkeypatch.setattr(np, "savetxt", lambda *args, **kwargs: None)

    result = dimer.get_forces(["H"], np.zeros(3))

    np.testing.assert_allclose(result["forces"], [-1.0, 0.0, 0.0])
    assert "Skipping calculation of norm(bias_forces)" in messages


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
    temperature = 298.15
    frequencies = np.array([20.0, 50.0, 1600.0]) * 100.0 * C
    q_vib, q_vib_v0 = qrrho_vibrational_part_func(
        temperature, frequencies, I_mean=10.0, cutoff=100.0, alpha=4,
    )
    harmonic, harmonic_v0 = vibrational_part_funcs(temperature, frequencies)
    weights = chai_head_gordon_weights(frequencies, 100.0, 4)
    mu = PLANCK / (8.0 * np.pi**2 * frequencies)
    inertia = 10.0 * 1.0e-20 * AMU2KG
    effective_inertia = mu * inertia / (mu + inertia)
    free_rotor = np.sqrt(
        8.0 * np.pi**3 * effective_inertia * KB * temperature / PLANCK**2
    )
    expected = np.exp(
        np.sum(weights * np.log(harmonic) + (1.0 - weights) * np.log(free_rotor))
    )
    expected_v0 = np.exp(
        np.sum(
            weights * np.log(harmonic_v0)
            + (1.0 - weights) * np.log(free_rotor)
        )
    )
    assert q_vib == pytest.approx(expected, rel=1.0e-13)
    assert q_vib_v0 == pytest.approx(expected_v0, rel=1.0e-13)


def test_vibrational_heat_capacity_is_stable_in_both_limits() -> None:
    high_frequency = np.array([3000.0]) * 100.0 * C
    assert vibrational_heat_capacity(1.0, high_frequency) == pytest.approx(0.0)
    low_frequency = np.array([1.0e-9])
    assert vibrational_heat_capacity(298.15, low_frequency) == pytest.approx(R)


def test_regularized_lbfgs_scales_the_secant_correction() -> None:
    result = bfgs_multiply(
        [np.array([1.0, -1.0])],
        [np.array([-2.0, 2.0])],
        np.array([0.3, -0.4]),
        gamma_mult=False,
        mu_reg=0.1,
    )
    np.testing.assert_allclose(result, [3.45, -3.55])


def test_gediis_returns_none_when_the_inner_solve_fails(monkeypatch) -> None:
    monkeypatch.setattr(
        gdiis_module,
        "minimize",
        lambda *_args, **_kwargs: SimpleNamespace(success=False),
    )
    result = gdiis_module.gediis(
        np.array([[0.0, 0.0], [0.1, 0.0]]),
        np.array([0.0, 0.1]),
        np.array([[0.1, 0.0], [0.2, 0.0]]),
    )
    assert result is None


def test_geom_loader_dispatches_trajectory_suffix_before_xyz(tmp_path) -> None:
    trajectory = tmp_path / "two_trj.xyz"
    trajectory.write_text(
        "1\nfirst\nH 0 0 0\n1\nsecond\nH 0 0 1\n",
        encoding="utf-8",
    )
    geometries = geom_loader(trajectory, iterable=True)
    assert len(geometries) == 2
    assert geom_loader(f"{trajectory}[1]").comment == "second"


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

    class Image:
        coords = np.zeros(1)

        def copy(self, **_kwargs):
            return self

        def __sub__(self, _other):
            return np.zeros(1)

    string = SimpleNamespace(
        images=[Image(), Image()],
        lf_ind=0,
        sk=0.1,
        get_cur_param_density=lambda: np.zeros(2),
        reset_geometries=lambda _image: None,
    )
    with pytest.raises(ValueError, match="zero path density"):
        GrowingString.get_new_image(string, 0)


def test_max_line_search_projects_endpoint_gradients(monkeypatch) -> None:
    captured = {}

    def capture_fit(**kwargs):
        captured.update(kwargs)
        return None

    monkeypatch.setattr(
        "pysisyphus.tsoptimizers.TSHessianOptimizer.poly_fit.quartic_fit",
        capture_fit,
    )
    step = np.array([0.5, -0.25])
    g0 = np.array([2.0, 4.0])
    g1 = np.array([-3.0, 1.0])
    optimizer = SimpleNamespace(
        max_line_search=True,
        min_line_search=False,
        cur_cycle=1,
        energies=[0.0, 0.2],
        forces=[-g0, -g1],
        steps=[step],
        logger=None,
        do_line_search=TSHessianOptimizer.do_line_search,
    )

    TSHessianOptimizer.step_and_grad_from_line_search(
        optimizer,
        0.2,
        g1,
        np.eye(2),
        np.array([], dtype=int),
        np.array([0, 1]),
    )

    assert captured["g0"] == pytest.approx(step.dot(g0))
    assert captured["g1"] == pytest.approx(step.dot(g1))
    assert captured["maximize"] is True


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


@pytest.mark.parametrize("act_dofs", [np.array([0, 1, 2]), [0, 1, 2]])
def test_irc_rms_gradient_uses_integration_basis(act_dofs) -> None:
    irc = IRC.__new__(IRC)
    irc._act_dofs = act_dofs
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


def test_euler_corrector_degrades_instead_of_aborting_on_oscillation(capsys) -> None:
    """An oscillating DWI must cost one corrector, not the whole IRC.

    The corrector descends the two-point DWI *interpolation*, not the real PES,
    so a reversal there is an interpolation artefact. Raising instead of
    returning the last non-oscillating point aborted complete ``mlmm all`` runs
    from inside a healthy IRC (smoke test73), because this branch is also the
    integration loop's escape hatch.
    """

    class OscillatingDWI:
        """1-D well at the origin; fixed-length descent overshoots and flips."""

        @staticmethod
        def interpolate(coords, gradient=True):
            return 0.0, coords.copy()

    irc = EulerPC.__new__(EulerPC)
    irc._m_sqrt = np.array([1.0])
    irc._act_dofs = np.array([0])
    irc.log = lambda *_: None
    start = np.array([0.02])

    corrected = irc.corrector_step(start, 0.1, OscillatingDWI())

    assert np.all(np.isfinite(corrected))
    # Degraded, not aborted: short of the requested 0.1 but still advancing, so
    # the caller gets a usable geometry and the IRC keeps going.
    advance = float(np.linalg.norm(corrected - start))
    assert 0.0 < advance < 0.1
    assert "oscillated" in capsys.readouterr().out


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


def test_find_bonds_reports_empty_bonds_as_index_pairs() -> None:
    bonds = find_bonds(["Ne"], np.zeros((1, 3)))

    assert bonds.shape == (0, 2)
    assert np.issubdtype(bonds.dtype, np.integer)
    # The default PDB serialization sorts along axis 1.
    assert np.sort(bonds, axis=1).shape == (0, 2)


def test_find_bonds_for_geom_forwards_the_requested_bond_factor() -> None:
    geom = Geometry(["C", "C"], np.array([0.0, 0.0, 0.0, 0.0, 0.0, 5.4]))

    assert find_bonds_for_geom(geom).shape == (0, 2)
    assert find_bonds_for_geom(geom, bond_factor=2.0).shape == (1, 2)


def test_defined_primitives_use_original_indices_with_leading_frozen_atoms() -> None:
    atoms = ["C"] * 8
    coords3d = np.zeros((8, 3))
    coords3d[:, 0] = np.arange(8) * 2.9

    red = RedundantCoords(
        atoms,
        coords3d,
        freeze_atoms=[0, 1],
        freeze_atoms_exclude=True,
        define_prims=[(PrimTypes.BOND, 2, 5)],
    )

    assert (PrimTypes.BOND, 2, 5) in red.typed_prims
    # Original indices must not be mapped a second time onto mobile atoms.
    assert (PrimTypes.BOND, 4, 7) not in red.typed_prims


def test_defined_primitives_on_excluded_frozen_atoms_are_rejected() -> None:
    atoms = ["C"] * 8
    coords3d = np.zeros((8, 3))
    coords3d[:, 0] = np.arange(8) * 2.9

    with pytest.raises(PrimitiveNotDefinedException):
        RedundantCoords(
            atoms,
            coords3d,
            freeze_atoms=[0, 1],
            freeze_atoms_exclude=True,
            define_prims=[(PrimTypes.BOND, 0, 5)],
        )


def test_dummy_torsion_gradient_is_finite_without_a_global_x_component() -> None:
    # The inner bond lies in the yz-plane, so a zero sign would collapse the
    # synthetic point onto the central atom.
    coords3d = np.array([
        [1.0, 0.4, -0.3],
        [0.0, 0.0, 0.0],
        [0.0, 2.1, 0.0],
    ])
    value, gradient = DummyTorsion._calculate(coords3d, [0, 1, 2], gradient=True)

    assert np.isfinite(value)
    assert np.isfinite(gradient).all()


def test_augment_bonds_preserves_frozen_atoms_and_isotopes(monkeypatch) -> None:
    atoms = ["C", "C", "H", "H"]
    coords3d = np.array([
        [0.0, 0.0, 0.0],
        [2.8, 0.0, 0.0],
        [-2.0, 0.0, 0.0],
        [4.8, 0.0, 0.0],
    ])
    geom = Geometry(
        atoms,
        coords3d.flatten(),
        coord_type="redund",
        freeze_atoms=[0],
        isotopes=((2, 2.0),),
    )
    geom.cart_hessian = np.zeros((12, 12))

    monkeypatch.setattr(
        augment_bonds_module, "find_missing_strong_bonds", lambda *a, **kw: [(0, 3)]
    )
    new_geom = augment_bonds_module.augment_bonds(geom)

    assert new_geom is not geom
    np.testing.assert_array_equal(new_geom.freeze_atoms, geom.freeze_atoms)
    assert new_geom.isotopes == geom.isotopes
    assert (PrimTypes.AUX_BOND, 0, 3) in new_geom.internal.typed_prims


def test_augment_bonds_keeps_the_remaining_coordinate_options(monkeypatch) -> None:
    atoms = ["C", "C", "H", "H"]
    coords3d = np.array([
        [0.0, 0.0, 0.0],
        [2.8, 0.0, 0.0],
        [-2.0, 0.0, 0.0],
        [4.8, 0.0, 0.0],
    ])
    geom = Geometry(
        atoms,
        coords3d.flatten(),
        coord_type="redund",
        coord_kwargs={
            "define_prims": [(PrimTypes.BOND, 1, 2)],
            "bonds_only": False,
        },
    )
    geom.cart_hessian = np.zeros((12, 12))

    monkeypatch.setattr(
        augment_bonds_module, "find_missing_strong_bonds", lambda *a, **kw: [(0, 3)]
    )
    new_geom = augment_bonds_module.augment_bonds(geom)

    assert new_geom.coord_kwargs["bonds_only"] is False
    assert new_geom.coord_kwargs["define_prims"] == [
        (PrimTypes.BOND, 1, 2),
        (PrimTypes.AUX_BOND, 0, 3),
    ]


def test_damped_bfgs_and_flowchart_updates_preserve_the_hessian_backend() -> None:
    dx = np.array([0.1, -0.2, 0.05])
    dg = np.array([0.4, 0.1, -0.3])
    H_np = np.array([
        [1.0, 0.2, 0.0],
        [0.2, 0.9, 0.1],
        [0.0, 0.1, 1.3],
    ])
    H_torch = torch.tensor(H_np, dtype=torch.float64)

    for update in (damped_bfgs_update, flowchart_update):
        expected, expected_key = update(H_np, dx, dg)
        actual, actual_key = update(H_torch, dx, dg)

        assert actual_key == expected_key
        assert isinstance(actual, torch.Tensor)
        assert actual.dtype == H_torch.dtype
        assert actual.device == H_torch.device
        np.testing.assert_allclose(actual.cpu().numpy(), expected, atol=1e-12)

    zeros, key = dummy_hessian_update(H_torch, dx, dg)
    assert key == "no"
    assert isinstance(zeros, torch.Tensor)
    assert zeros.dtype == H_torch.dtype
    assert bool((zeros == 0.0).all())


def test_gediis_weights_use_the_quadratic_form_of_the_inverse_hessian() -> None:
    coords = np.array([[0.0, 0.0], [0.2, 0.0], [0.1, 0.3]])
    energies = np.array([0.3, 0.1, 0.2])
    forces = np.array([[0.2, 0.1], [0.05, -0.1], [-0.1, 0.15]])
    # Strongly off-diagonal, so a row-summed contraction differs from f^T H^-1 f.
    hessian = np.array([[1.0, 0.8], [0.8, 1.2]])
    hessian_inv = np.linalg.pinv(hessian, rcond=1e-6)

    R = coords[::-1]
    f = forces[::-1]
    Rifi = np.einsum("ik,ik->i", R, f)
    Rjfi = np.einsum("jk,ik->ji", R, f)
    quadratic_form = np.einsum("ki,ij,kj->k", f, hessian_inv, f)
    row_summed = np.einsum("ki,ji,ki->k", f, hessian_inv, f)
    assert quadratic_form[0] != pytest.approx(row_summed[0])

    captured = {}
    original = gdiis_module.minimize

    def spy(fun, *args, **kwargs):
        captured["value"] = fun(np.array([1.0, 0.0, 0.0]))
        return original(fun, *args, **kwargs)

    gdiis_module.minimize = spy
    try:
        gdiis_module.gediis(coords, energies, forces, hessian=hessian)
    finally:
        gdiis_module.minimize = original

    # Eq. (5) of the reference at the first vertex.
    expected = 0.5 * quadratic_form[0] - Rjfi[0, 0] + Rifi[0]
    assert captured["value"] == pytest.approx(expected)


def test_quartic_fit_returns_none_for_a_degenerate_quadratic_coefficient() -> None:
    # Equal endpoint energies with vanishing projected gradients.
    assert quartic_fit(-1.0, -1.0, 0.0, 0.0) is None
    # A well conditioned fit still interpolates.
    assert quartic_fit(0.371, 0.301, 0.377, -0.222) is not None


def test_lanczos_stops_at_an_exact_residual_breakdown() -> None:
    # Quadratic PES, so the gradient difference is exactly H @ dx.
    hessian = np.diag([-0.5, 1.0, 2.0])

    def grad_getter(coords):
        return hessian @ coords

    w_min, mode = lanczos(
        np.zeros(3),
        grad_getter,
        guess=np.array([1.0, 0.0, 0.0]),
        max_cycles=10,
    )

    assert w_min == pytest.approx(-0.5, abs=1e-6)
    assert np.isfinite(mode).all()

    with pytest.raises(ValueError):
        lanczos(np.zeros(3), grad_getter, guess=np.zeros(3))


def test_normal_modes_retain_every_low_complement_root() -> None:
    # Rank-zero constrained PHVA: the whole active block is the complement, so
    # every root, including deliberately tiny ones, must survive.
    atomic_numbers = [6, 1, 1, 1]
    coords_bohr = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.0, 0.0, 2.0],
    ])
    # One active atom, so the constrained rigid space has rank zero.
    active_hessian = torch.diag(torch.tensor([3.0e-9, -4.0e-9, 0.5], dtype=torch.float64))

    freqs, modes = nm._frequencies_cm_and_modes(
        active_hessian.clone(),
        atomic_numbers,
        coords_bohr,
        torch.device("cpu"),
        freeze_idx=[1, 2, 3],
        frequency_zero_cutoff_cm=0.0,
    )

    assert len(freqs) == 3
    assert modes.shape == (3, 12)
    assert (freqs < 0.0).sum() == 1
    # The tiny roots are far below the historical 5.14 cm^-1 magnitude floor.
    assert abs(freqs[freqs < 0.0][0]) < 1.0


def test_imaginary_frequencies_exclude_small_positive_eigenvalues() -> None:
    geom = Geometry(["H", "H"], np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.4]))

    def fake_normal_modes(hessian=None):
        return np.array([3.0, -7.0]), np.array([1.0e-9, -2.0e-6]), None

    geom.get_normal_modes = fake_normal_modes
    imag = geom.get_imag_frequencies()

    np.testing.assert_allclose(imag, [-7.0])


def test_trans_rot_vectors_are_translation_invariant_for_a_linear_molecule() -> None:
    coords = np.array([0.0, 0.0, -1.4, 0.0, 0.0, 0.0, 0.0, 0.0, 1.4])
    masses = np.array([12.0, 12.0, 12.0])

    here = get_trans_rot_vectors(coords, masses)
    shifted = get_trans_rot_vectors(coords + np.tile([7.0, -3.0, 5.0], 3), masses)

    # A linear molecule keeps rigid rank five under any translation.
    assert here.shape[0] == 5
    assert shifted.shape[0] == here.shape[0]


def test_exact_phva_order_ignores_subthreshold_negative_root() -> None:
    # A strong reaction mode plus a numerically soft negative root below the
    # configured 5 cm^-1 saddle/export/recovery threshold.
    freqs_cm = np.array([-450.0, -3.2, 12.0])
    modes = torch.eye(3, dtype=torch.float64)

    # RSPRFOptimizer is the concrete TSHessianOptimizer used by the product.
    opt = RSPRFOptimizer.__new__(RSPRFOptimizer)
    opt._mw_frequencies_and_modes = lambda: (freqs_cm, modes)
    opt._recovery_mode_from_mw = lambda _modes, index: np.eye(3)[index]
    opt.geometry = SimpleNamespace(cart_coords=np.zeros(3))
    opt.reference_mode = None
    opt.roots = [0]
    opt.saddle_imaginary_threshold_cm = 5.0
    opt.higher_order_saddle_checks = 0
    opt.max_higher_order_checks = 99
    opt.cur_cycle = 7
    opt.table = SimpleNamespace(print=lambda *_a, **_kw: None)
    opt.request_stop = lambda *_a: None
    opt._record_exact_saddle_candidate = lambda: None
    opt._last_exact_target_mode_reanchored = False

    has_saddle_modes, physical_mode, verified = opt._verify_exact_vibrational_structure(
        None, None
    )

    # The soft -3.2 cm^-1 root does not change the certified saddle order.
    assert opt._last_exact_n_imaginary == 1
    assert opt._last_exact_saddle_verified is True
    assert opt._last_exact_saddle_cycle == 7
    assert has_saddle_modes is True
    assert verified is True
    assert physical_mode is not None


def test_bonded_fragment_jacobian_embeds_the_bond_second_derivative() -> None:
    from pysisyphus.intcoords.derivatives import d2q_b

    coords3d = np.array([
        [4.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ])
    fragment = BondedFragment([0, 1], bond_indices=[1, 2])

    expected_block = d2q_b(*coords3d[1], *coords3d[2]).reshape(6, 6)
    np.testing.assert_allclose(
        fragment.jacobian(coords3d).reshape(6, 6), expected_block
    )

    internal = RedundantCoords.__new__(RedundantCoords)
    internal.coords3d = coords3d
    internal.primitives = [fragment]
    K = internal.get_K_matrix([1.0])

    expected = np.zeros((9, 9))
    endpoint_coords = [3, 4, 5, 6, 7, 8]
    expected[np.ix_(endpoint_coords, endpoint_coords)] = expected_block
    np.testing.assert_allclose(K, expected)


def test_collapsed_path_returns_zero_coord_diffs_instead_of_nans() -> None:
    from pysisyphus.helpers import get_coords_diffs

    collapsed = get_coords_diffs(np.zeros((4, 6)))
    np.testing.assert_array_equal(collapsed, np.zeros(4))
    assert np.isfinite(collapsed).all()

    # A nondegenerate path is still normalized to one.
    spread = get_coords_diffs(np.array([[0.0], [1.0], [3.0]]))
    np.testing.assert_allclose(spread, [0.0, 1.0 / 3.0, 1.0])


def test_parallel_cos_restores_pal_on_the_images_it_evaluated() -> None:
    class _Calc:
        def __init__(self, pal):
            self.pal = pal

    images = [
        SimpleNamespace(calculator=_Calc(pal)) for pal in (8, 6, 4, 2)
    ]

    class _Client:
        @staticmethod
        def scheduler_info():
            return {"workers": {"w0": {}, "w1": {}}}

        @staticmethod
        def map(_func, items):
            return list(items)

        @staticmethod
        def gather(futures):
            return list(futures)

    cos = ChainOfStates.__new__(ChainOfStates)
    cos.images = images
    cos.log = lambda *_a, **_kw: None
    cos.get_dask_client = lambda: _Client()

    # Cached endpoints: only the two moving images are evaluated, so restoring
    # by position would leave image 3 with a reduced pal.
    image_indices = [1, 3]
    cos.concurrent_force_calcs([images[1], images[3]], image_indices)

    assert [image.calculator.pal for image in images] == [8, 6, 4, 2]


def test_euler_corrector_keeps_the_last_advancing_point_on_immediate_reversal(
    capsys,
) -> None:
    """An immediate DWI reversal must not return the zero-advance start point."""

    class FlippingDWI:
        """Gradient sign flips on every call, so the second step reverses."""

        def __init__(self):
            self.calls = 0

        def interpolate(self, coords, gradient=True):
            self.calls += 1
            sign = 1.0 if self.calls % 2 else -1.0
            return 0.0, np.array([sign])

    irc = EulerPC.__new__(EulerPC)
    irc._m_sqrt = np.array([1.0])
    irc._act_dofs = np.array([0])
    irc.log = lambda *_: None
    start = np.zeros(1)

    corrected = irc.corrector_step(start, 0.1, FlippingDWI())

    assert np.all(np.isfinite(corrected))
    # The last finite advancing microstep, not the original point.
    assert float(np.linalg.norm(corrected - start)) > 0.0
    assert "oscillated" in capsys.readouterr().out


def test_irc_integration_failure_outranks_a_small_gradient() -> None:
    irc = IRC.__new__(IRC)
    irc.never_stop = False
    irc.past_inflection = True
    irc.rms_grad_thresh = 1.0
    irc.hard_rms_grad_thresh = None
    irc.energy_increase_thresh = 0.0
    irc.energy_thresh = 1.0e-6
    irc.require_pos_def_hessian = False

    # Integration failure in the same macrostep as a converged gradient.
    irc.integration_stop_requested = True
    irc.integration_stop_reason = "Predictor integration exhausted."
    irc.converged = False
    irc.energy_increased = False
    irc.energy_converged = False
    assert irc._gradient_converged(0.0) is True
    # The numerical failure has unconditional priority, so this direction is
    # never published as converged.
    assert irc.integration_stop_requested and not irc.converged


def test_irc_ordinary_energy_rise_outranks_physical_convergence() -> None:
    irc = IRC.__new__(IRC)
    irc.never_stop = False
    irc.past_inflection = True
    irc.rms_grad_thresh = 1.0
    irc.energy_increase_thresh = 1.0e-6

    # An ordinary-mode energy rise is reported instead of convergence.
    assert irc._energy_increase_exceeds_tolerance(-10.0, -9.0) is True
    irc.energy_increased = True
    irc.energy_converged = False
    assert irc._energy_stop_message() == "Energy increased!"

    # never_stop still bypasses the physical energy stop only.
    irc.never_stop = True
    assert irc._energy_stop_message() == ""
