"""Energy/force/Hessian consistency for harmonic distance restraints."""

from __future__ import annotations

import numpy as np
import pytest
import torch

from pysisyphus.constants import ANG2BOHR
from pysisyphus.Geometry import Geometry as PySiGeometry

from mlmm.workflows.freq import _calc_full_hessian_torch
from mlmm.workflows.opt import (
    _flatten_all_imag_modes_for_geom,
    _seed_rfo_initial_hessian,
)
from mlmm.workflows.restraints import (
    HarmonicFixAtoms,
    HarmonicBiasCalculator,
    harmonic_pair_energy_forces_hessian,
)


@pytest.mark.parametrize("k_fix", [0.0, -1.0, np.nan, np.inf])
def test_harmonic_fix_atoms_rejects_invalid_force_constant(k_fix) -> None:
    with pytest.raises(ValueError, match="finite positive"):
        HarmonicFixAtoms(
            indices=[0],
            ref_positions=np.zeros((1, 3)),
            k_fix=k_fix,
        )


def _pair_system():
    coords = np.array(
        [[0.0, 0.0, 0.0], [1.35 * ANG2BOHR, 0.21 * ANG2BOHR, 0.0]],
        dtype=float,
    )
    return coords, [(0, 1, 1.0)]


def test_harmonic_pair_hessian_matches_force_finite_difference() -> None:
    coords, pairs = _pair_system()
    k = 0.37
    _, _, analytical = harmonic_pair_energy_forces_hessian(
        coords, k, pairs, need_hessian=True
    )
    assert analytical is not None

    eps = 1.0e-5
    numerical = np.zeros_like(analytical)
    for column in range(coords.size):
        plus = coords.copy().reshape(-1)
        minus = coords.copy().reshape(-1)
        plus[column] += eps
        minus[column] -= eps
        _, force_plus, _ = harmonic_pair_energy_forces_hessian(
            plus, k, pairs, need_hessian=False
        )
        _, force_minus, _ = harmonic_pair_energy_forces_hessian(
            minus, k, pairs, need_hessian=False
        )
        numerical[:, column] = -(force_plus - force_minus) / (2.0 * eps)

    assert np.allclose(analytical, numerical, atol=2.0e-10, rtol=2.0e-9)
    assert np.allclose(analytical, analytical.T, atol=0.0, rtol=0.0)


def test_harmonic_pair_force_is_negative_energy_gradient() -> None:
    coords, pairs = _pair_system()
    k = 0.37
    _, force, _ = harmonic_pair_energy_forces_hessian(
        coords, k, pairs, need_hessian=False
    )

    eps = 1.0e-5
    gradient = np.zeros(coords.size, dtype=float)
    for column in range(coords.size):
        plus = coords.copy().reshape(-1)
        minus = coords.copy().reshape(-1)
        plus[column] += eps
        minus[column] -= eps
        e_plus, _, _ = harmonic_pair_energy_forces_hessian(
            plus, k, pairs, need_hessian=False
        )
        e_minus, _, _ = harmonic_pair_energy_forces_hessian(
            minus, k, pairs, need_hessian=False
        )
        gradient[column] = (e_plus - e_minus) / (2.0 * eps)

    assert np.allclose(force, -gradient, atol=2.0e-10, rtol=2.0e-9)


@pytest.mark.parametrize(
    ("coords", "k", "pairs", "message"),
    [
        (
            np.zeros((2, 3)),
            1.0,
            [(0, 0, 1.0)],
            "two distinct atoms",
        ),
        (
            np.zeros((2, 3)),
            1.0,
            [(0, 2, 1.0)],
            "outside the valid range",
        ),
        (
            np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            1.0,
            [(0, 1, 0.0)],
            "greater than zero",
        ),
        (
            np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            -1.0,
            [(0, 1, 1.0)],
            "non-negative",
        ),
        (
            np.zeros((2, 3)),
            1.0,
            [(0, 1, 1.0)],
            "coincident atoms",
        ),
    ],
)
def test_invalid_harmonic_pair_restraints_fail_closed(
    coords: np.ndarray,
    k: float,
    pairs,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        harmonic_pair_energy_forces_hessian(
            coords,
            k,
            pairs,
            need_hessian=True,
        )


class _StaticBase:
    def __init__(self, result, *, freeze_atoms=()):
        self.result = result
        self.freeze_atoms = list(freeze_atoms)
        self.hessian_calls = 0

    def get_energy(self, elem, coords):
        return self.result

    def get_forces(self, elem, coords):
        return self.result

    def get_hessian(self, elem, coords):
        self.hessian_calls += 1
        return self.result


@pytest.mark.parametrize("use_torch", [False, True])
def test_restrained_hessian_preserves_full_representation_and_final_mask(
    use_torch: bool,
) -> None:
    coords, pairs = _pair_system()
    hessian = np.eye(6, dtype=np.float32)
    if use_torch:
        hessian = torch.as_tensor(hessian)
    base_result = {
        "energy": 1.0,
        "forces": np.ones(6, dtype=np.float32),
        "hessian": hessian,
    }
    base = _StaticBase(base_result, freeze_atoms=[0])
    wrapper = HarmonicBiasCalculator(base, k=2.0, pairs=pairs)

    result = wrapper.get_hessian(["H", "H"], coords)

    assert isinstance(result["hessian"], torch.Tensor) is use_torch
    assert result["hessian"].dtype == hessian.dtype
    if use_torch:
        assert result["hessian"].device == hessian.device
    hessian_np = (
        result["hessian"].detach().cpu().numpy()
        if use_torch
        else result["hessian"]
    )
    assert np.array_equal(hessian_np[:3, :], np.zeros((3, 6)))
    assert np.array_equal(hessian_np[:, :3], np.zeros((6, 3)))
    assert np.array_equal(result["forces"][:3], np.zeros(3))
    assert base_result["energy"] == 1.0
    assert np.array_equal(base_result["forces"], np.ones(6, dtype=np.float32))
    assert np.array_equal(
        hessian.detach().cpu().numpy() if use_torch else hessian,
        np.eye(6, dtype=np.float32),
    )


def test_restrained_hessian_maps_through_partial_metadata() -> None:
    coords, pairs = _pair_system()
    metadata = {
        "active_atoms": np.array([1]),
        "active_dofs": np.array([3, 4, 5]),
        "active_n_dof": 3,
        "full_n_dof": 6,
    }
    base_result = {
        "energy": 1.0,
        "forces": np.ones(6, dtype=np.float64),
        "hessian": torch.eye(3, dtype=torch.float64),
        "within_partial_hessian": metadata,
    }
    wrapper = HarmonicBiasCalculator(
        _StaticBase(base_result, freeze_atoms=[0]), k=2.0, pairs=pairs
    )

    result = wrapper.get_hessian(["H", "H"], coords)

    assert result["hessian"].shape == (3, 3)
    assert result["hessian"].dtype == torch.float64
    assert result["within_partial_hessian"] is metadata
    assert np.array_equal(result["forces"][:3], np.zeros(3))


def test_restrained_full_hessian_preserves_core_hessian_only_mask() -> None:
    coords, pairs = _pair_system()
    base_hessian = np.eye(6, dtype=np.float64)
    base_hessian[3:, :] = 0.0
    base_hessian[:, 3:] = 0.0

    class RefreshingBase(_StaticBase):
        def __init__(self):
            super().__init__(
                {
                    "energy": 0.0,
                    "forces": np.zeros(6, dtype=np.float64),
                    "hessian": base_hessian,
                },
                freeze_atoms=[],
            )
            self.core = type(
                "Core",
                (),
                {"effective_hess_freeze_atoms": []},
            )()

        def get_hessian(self, elem, coords):
            self.core.effective_hess_freeze_atoms = [1]
            return super().get_hessian(elem, coords)

    wrapper = HarmonicBiasCalculator(RefreshingBase(), k=2.0, pairs=pairs)

    result = wrapper.get_hessian(["H", "H"], coords)

    assert np.array_equal(result["hessian"][3:, :], np.zeros((3, 6)))
    assert np.array_equal(result["hessian"][:, 3:], np.zeros((6, 3)))
    assert np.any(np.abs(result["forces"][3:]) > 0.0)


def test_empty_restraints_match_base_values_and_clone_buffers() -> None:
    coords, _ = _pair_system()
    base_result = {
        "energy": 1.25,
        "forces": np.arange(6, dtype=float),
        "hessian": torch.arange(36, dtype=torch.float64).reshape(6, 6),
    }
    wrapper = HarmonicBiasCalculator(_StaticBase(base_result), k=2.0, pairs=[])

    result = wrapper.get_hessian(["H", "H"], coords)

    assert result["energy"] == base_result["energy"]
    assert np.array_equal(result["forces"], base_result["forces"])
    assert torch.equal(result["hessian"], base_result["hessian"])
    assert result["forces"] is not base_result["forces"]
    assert result["hessian"] is not base_result["hessian"]


class _Geometry:
    atoms = ["H", "H"]
    cart_coords = np.zeros(6, dtype=float)
    within_partial_hessian = None


def test_explicit_calculator_evaluates_the_exact_active_pes() -> None:
    geometry = _Geometry()
    expected = torch.eye(6, dtype=torch.float64) * 2.0
    evaluator = _StaticBase(
        {
            "energy": 3.5,
            "forces": np.zeros(6),
            "hessian": expected,
        }
    )

    result, energy = _calc_full_hessian_torch(
        geometry,
        {"input_pdb": "must-not-be-constructed"},
        torch.device("cpu"),
        calculator=evaluator,
    )

    assert evaluator.hessian_calls == 1
    assert torch.equal(result, expected)
    assert energy == 3.5


def test_restrained_rfo_seed_uses_exact_wrapper_and_never_reads_irc_cache(
    monkeypatch,
) -> None:
    class Geometry:
        cart_coords = np.zeros(6, dtype=float)
        cart_hessian = None

    geometry = Geometry()
    evaluator = object()
    expected = torch.eye(6, dtype=torch.float64) * 7.0
    seen = []

    def fake_hessian(geom, cfg, device, *, calculator=None, **kwargs):
        seen.append((geom, cfg, device, calculator, kwargs))
        return expected, 0.0

    def forbidden_cache(*_args, **_kwargs):
        raise AssertionError("restrained RFO must not read the IRC cache")

    monkeypatch.setattr(
        "mlmm.workflows.freq._calc_full_hessian_torch", fake_hessian
    )
    monkeypatch.setattr("mlmm.io.hessian_cache.load", forbidden_cache)

    source = _seed_rfo_initial_hessian(
        geometry,
        {"ml_device": "cpu", "backend": "sentinel"},
        evaluator,
        restraints_active=True,
    )

    assert source == "restrained"
    assert geometry.cart_hessian is expected
    assert seen[0][0] is geometry
    assert seen[0][1] == {"ml_device": "cpu", "backend": "sentinel"}
    assert seen[0][2] == torch.device("cpu")
    assert seen[0][3] is evaluator
    assert seen[0][4] == {"refresh_geom_meta": True}


def test_flatten_energy_probes_use_exact_active_evaluator() -> None:
    class Geometry:
        atoms = ["H", "H"]

        def __init__(self):
            self._cart_coords = np.zeros(6, dtype=float)
            self._internal_coords = np.array([91.0, 92.0], dtype=float)

        @property
        def cart_coords(self):
            return self._cart_coords

        @cart_coords.setter
        def cart_coords(self, value):
            self._cart_coords = np.asarray(value, dtype=float).copy()

        @property
        def coords(self):
            return self._internal_coords

        @coords.setter
        def coords(self, value):
            raise AssertionError("Cartesian trials must not use the internal setter")

    class Evaluator:
        def __init__(self):
            self.calls = []

        def get_energy(self, atoms, coords):
            self.calls.append((atoms, np.asarray(coords).copy()))
            return {"energy": float(np.sum(np.asarray(coords) ** 2))}

    geometry = Geometry()
    evaluator = Evaluator()

    did_flatten = _flatten_all_imag_modes_for_geom(
        geometry,
        np.array([1.0, 1.0]),
        {"backend": "sentinel"},
        np.array([-100.0]),
        torch.tensor([[1.0, 0.0, 0.0, -1.0, 0.0, 0.0]]),
        5.0,
        0.1,
        calculator=evaluator,
    )

    assert did_flatten is True
    assert len(evaluator.calls) == 3
    assert all(atoms is geometry.atoms for atoms, _coords in evaluator.calls)
    assert all(coords.shape == (6,) for _atoms, coords in evaluator.calls)
    assert np.array_equal(geometry.coords, np.array([91.0, 92.0]))


def test_flatten_real_internal_geometry_probes_cartesian_coordinates() -> None:
    atoms = ["C", "C", "O", "C", "N", "H"]
    coords = np.array(
        [
            0.0, 0.0, 0.0,
            1.4, 0.0, 0.0,
            2.1, 1.2, 0.0,
            3.5, 1.2, 0.1,
            4.2, 0.1, 0.2,
            5.6, 0.1, 0.0,
        ],
        dtype=float,
    )
    geometry = PySiGeometry(atoms, coords, coord_type="redund")
    assert geometry.coords.size != geometry.cart_coords.size

    class Evaluator:
        def __init__(self):
            self.coords = []

        def get_energy(self, _atoms, cart_coords):
            snapshot = np.asarray(cart_coords, dtype=float).copy()
            self.coords.append(snapshot)
            return {"energy": float(np.sum(snapshot ** 2))}

    evaluator = Evaluator()
    mode = torch.zeros((1, coords.size), dtype=torch.float64)
    mode[0, 0] = 1.0
    mode[0, 3] = -1.0

    did_flatten = _flatten_all_imag_modes_for_geom(
        geometry,
        np.ones(len(atoms), dtype=float),
        {"backend": "sentinel"},
        np.array([-100.0]),
        mode,
        5.0,
        0.1,
        calculator=evaluator,
    )

    assert did_flatten is True
    assert len(evaluator.coords) == 3
    assert all(snapshot.shape == (coords.size,) for snapshot in evaluator.coords)
    np.testing.assert_allclose(
        0.5 * (evaluator.coords[1] + evaluator.coords[2]),
        evaluator.coords[0],
        atol=1.0e-12,
    )
    selected = min(evaluator.coords[1:], key=lambda point: float(np.sum(point ** 2)))
    np.testing.assert_allclose(geometry.cart_coords, selected, atol=1.0e-12)
