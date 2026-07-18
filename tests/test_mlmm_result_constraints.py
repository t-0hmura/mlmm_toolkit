"""Final ML/MM force/Hessian masking after link and embed composition."""

from __future__ import annotations

from types import MethodType, SimpleNamespace

import numpy as np
import torch
from ase import Atoms

from mlmm.backends.mlmm_calc import (
    MLMMASECalculator,
    MLMMCore,
    _MLHighOut,
    _MMLowOut,
)


def _force_only_core(*, with_link: bool) -> MLMMCore:
    core = object.__new__(MLMMCore)
    core._n_real = 2
    core.freeze_atoms = [0]
    core.hess_freeze_atoms = []
    core.selection_indices = [0]
    core._idx_map_real_to_model = {0: 0}
    core.print_vram = False
    core.print_timing = False
    core.ml_device = torch.device("cpu")
    core.calc_real_low = SimpleNamespace(device="cpu")
    core.embedcharge = False
    core._embed_correction = None
    core.link_atom_method = "scaled"

    atoms_real = Atoms("HH", positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    atoms_model = Atoms("H", positions=[[0.0, 0.0, 0.0]])
    atoms_model_lh = (
        Atoms("HH", positions=[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
        if with_link
        else atoms_model.copy()
    )
    core._atoms_model_tpl = atoms_model.copy()

    def prep(self, _coords):
        links = [(1, 0, 1, 0.5)] if with_link else []
        return (
            atoms_real.copy(),
            atoms_model.copy(),
            atoms_model_lh.copy(),
            links,
            [0],
        )

    def eval_high(self, _atoms, _freeze_model, *, return_hessian):
        forces = (
            np.array([[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]])
            if with_link
            else np.zeros((1, 3), dtype=float)
        )
        return _MLHighOut(E=0.0, F=forces, H=None, timing={})

    def eval_low(self, _atoms_real, _atoms_model, *, return_hessian):
        return _MMLowOut(
            E_real=0.0,
            F_real=np.zeros((2, 3), dtype=float),
            E_model=0.0,
            F_model=np.zeros((1, 3), dtype=float),
            H_real=None,
            H_model=None,
            active_atoms_from_fd=None,
            timing={},
        )

    core._prep_3_layer_atoms = MethodType(prep, core)
    core._eval_ml_high = MethodType(eval_high, core)
    core._eval_mm_low = MethodType(eval_low, core)
    return core


def test_final_mask_runs_after_nonzero_link_force_redistribution() -> None:
    core = _force_only_core(with_link=True)

    result = core.compute(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        return_forces=True,
    )

    assert np.array_equal(result["forces"][0], np.zeros(3))
    assert np.array_equal(result["forces"][1], np.array([2.0, 0.0, 0.0]))


def test_partial_link_hessian_retains_the_unconstrained_endpoint_block() -> None:
    core = _force_only_core(with_link=True)
    core.mm_fd = True
    core.H_dtype = torch.float64
    core.H_np_dtype = np.float64
    core.symmetrize_hessian = True
    core.return_partial_hessian = True

    def eval_high(self, _atoms, _freeze_model, *, return_hessian):
        hessian = torch.zeros((2, 3, 2, 3), dtype=torch.float64)
        hessian[1, :, 1, :] = torch.eye(3, dtype=torch.float64)
        return _MLHighOut(
            E=0.0,
            F=np.zeros((2, 3), dtype=float),
            H=hessian,
            timing={},
        )

    def eval_low(self, _atoms_real, _atoms_model, *, return_hessian):
        return _MMLowOut(
            E_real=0.0,
            F_real=np.zeros((2, 3), dtype=float),
            E_model=0.0,
            F_model=np.zeros((1, 3), dtype=float),
            H_real=np.zeros((3, 3), dtype=float),
            H_model=np.zeros((3, 3), dtype=float),
            active_atoms_from_fd=np.array([1], dtype=int),
            timing={},
        )

    core._eval_ml_high = MethodType(eval_high, core)
    core._eval_mm_low = MethodType(eval_low, core)

    result = core.compute(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        return_hessian=True,
    )

    assert result["hessian"].shape == (1, 3, 1, 3)
    assert torch.equal(
        result["hessian"].reshape(3, 3),
        torch.eye(3, dtype=torch.float64).mul(0.25),
    )
    assert np.array_equal(
        result["within_partial_hessian"]["active_atoms"], np.array([1])
    )


def test_final_mask_runs_after_nonzero_embed_force_correction() -> None:
    core = _force_only_core(with_link=False)
    core.embedcharge = True
    core.embedcharge_cutoff = None
    core.model_charge = 0
    core.model_mult = 1

    class EmbedCorrection:
        def compute_correction(self, **_kwargs):
            return 0.0, np.array([[7.0, 0.0, 0.0]]), None

    core._embed_correction = EmbedCorrection()
    core._get_mm_charges = MethodType(
        lambda self, indices: np.zeros(len(indices), dtype=float), core
    )

    result = core.compute(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        return_forces=True,
    )

    assert np.array_equal(result["forces"][0], np.zeros(3))


def _mapping_core() -> MLMMCore:
    core = object.__new__(MLMMCore)
    core._n_real = 3
    core.freeze_atoms = [0]
    core.hess_freeze_atoms = [2]
    core.selection_indices = [0, 1]
    core._idx_map_real_to_model = {0: 0, 1: 1}
    core._update_active_dof_mappings()
    return core


def test_partial_hessian_metadata_excludes_all_effective_constraints() -> None:
    core = _mapping_core()

    assert core.effective_hess_freeze_atoms == [0, 2]
    assert core.hess_active_atoms == [1]
    metadata = core._build_within_partial_hessian()
    assert np.array_equal(metadata["active_atoms"], np.array([1]))
    assert np.array_equal(metadata["active_dofs"], np.array([3, 4, 5]))

    hessian = torch.eye(3, dtype=torch.float64).reshape(1, 3, 1, 3)
    result = core._finalize_result_constraints(
        {
            "energy": 0.0,
            "forces": np.ones((3, 3), dtype=float),
            "hessian": hessian,
            "within_partial_hessian": metadata,
        }
    )
    assert result["hessian"] is hessian
    assert np.array_equal(result["forces"][0], np.zeros(3))
    assert np.array_equal(result["forces"][2], np.ones(3))


def test_full_hessian_final_mask_is_exact_and_preserves_symmetry() -> None:
    core = _mapping_core()
    raw = torch.arange(81, dtype=torch.float64).reshape(9, 9)
    symmetric = (raw + raw.T).mul(0.5)
    active_block_before = symmetric[3:6, 3:6].clone()
    hessian = symmetric.reshape(3, 3, 3, 3)
    result = core._finalize_result_constraints(
        {
            "energy": 0.0,
            "forces": np.ones((3, 3), dtype=float),
            "hessian": hessian,
        }
    )

    square = result["hessian"].reshape(9, 9)
    constrained = [0, 1, 2, 6, 7, 8]
    assert torch.equal(square[constrained, :], torch.zeros((6, 9)))
    assert torch.equal(square[:, constrained], torch.zeros((9, 6)))
    assert torch.equal(square[3:6, 3:6], active_block_before)
    assert torch.equal(square, square.T)
    assert np.array_equal(result["forces"][0], np.zeros(3))
    assert np.array_equal(result["forces"][2], np.ones(3))


def test_ase_adapter_exposes_finalized_core_forces_without_recomposition() -> None:
    expected = np.array([[0.0, 0.0, 0.0], [2.0, -1.0, 3.0]])

    class Core:
        def compute(self, coords, *, return_forces, return_hessian):
            assert return_forces is True
            assert return_hessian is False
            return {"energy": 1.5, "forces": expected.copy()}

    atoms = Atoms("HH", positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    atoms.calc = MLMMASECalculator(core=Core())

    assert np.array_equal(atoms.get_forces(), expected)
    assert atoms.get_potential_energy() == 1.5
