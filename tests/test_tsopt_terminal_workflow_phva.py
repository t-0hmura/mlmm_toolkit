"""Terminal PHVA scope and raw-Hessian ordering through the real CLI workflow."""

import json
from types import SimpleNamespace

import numpy as np
import pytest
import torch
from click.testing import CliRunner

from pysisyphus.Geometry import Geometry
from pysisyphus.normal_modes import DEFAULT_FREQUENCY_ZERO_CUTOFF_CM, resolved_imaginary_mask
from mlmm.core import calc_eval
from mlmm.io import hessian_cache
from mlmm.workflows import freq, tsopt


@pytest.mark.parametrize("widen", [False, True], ids=["reuse", "wider-final"])
@pytest.mark.parametrize("soft_added_root", [False, True], ids=["resolved-extra", "soft-extra"])
@pytest.mark.parametrize("final_energy", ["finite", "raises", "nonfinite"])
def test_real_terminal_workflow_phva_scope_and_raw_order(tmp_path, monkeypatch, widen, soft_added_root, final_energy):
    # Three noncollinear fixed anchors remove all compatible rigid motions.
    # E = 1/2 (x-x0)^T D (x-x0), with negative x curvature on atoms3 and4.
    # The macro block contains atoms3,5; restoring atom4 adds exactly one root.
    x0 = np.array([[0., 0., 0.], [2., 0., 0.], [0., 2., 0.],
                   [0., 0., 2.], [2., 2., 2.], [3., 1., 2.]]).ravel()
    diagonal = np.arange(1., 19.) * 0.01
    diagonal[[9, 12]] *= -1
    if soft_added_root:
        diagonal[12] = -1e-8
    full_h = torch.diag(torch.as_tensor(diagonal, dtype=torch.float64))
    macro_frozen = [0, 1, 2, 4]
    final_frozen = [0, 1, 2] if widen else macro_frozen
    final_atoms = [3, 4, 5] if widen else [3, 5]
    macro_order = [5, 3]  # optimizer-owned rows differ from calculator metadata
    calls, exported, observed = [], [], {}

    def dofs(atoms):
        return [3 * atom + axis for atom in atoms for axis in range(3)]

    class QuadraticCalculator:
        def __init__(self, **kwargs):
            self.freeze_atoms = list(kwargs.get("freeze_atoms", []))
            released = 4 not in self.freeze_atoms
            self.order = [4, 5, 3] if released else [3, 5]
            self.core = SimpleNamespace(
                ml_indices=[3, 5], hess_mm_indices=[4] if released else [],
                movable_mm_indices=[4] if released else [],
                frozen_layer_indices=self.freeze_atoms,
                hess_active_atoms=self.order,
            )

        def get_energy(self, atoms, coords):
            displacement = np.asarray(coords) - x0
            return {"energy": float(0.5 * np.dot(diagonal * displacement, displacement))}

        def get_forces(self, atoms, coords):
            return {**self.get_energy(atoms, coords),
                    "forces": -diagonal * (np.asarray(coords) - x0)}

        def get_hessian(self, atoms, coords):
            idx = dofs(self.order)
            raw = full_h[idx][:, idx].clone()
            calls.append({"order": list(self.order), "raw": raw.clone()})
            return {**self.get_forces(atoms, coords), "hessian": raw,
                    "within_partial_hessian": {
                        "active_atoms": list(self.order), "active_dofs": idx,
                        "active_n_dof": len(idx), "full_n_dof": x0.size,
                    }}

    def load_geometry(_path, **kwargs):
        geom = Geometry(["H"] * 6, x0.copy(), **kwargs)
        observed["geometry"] = geom
        return geom

    def completed_macro(geom, calc_cfg, *args, **kwargs):
        assert list(geom.freeze_atoms) == final_frozen
        entry_frozen = list(geom.freeze_atoms)
        geom.freeze_atoms = macro_frozen
        macro_cfg = {**calc_cfg, "freeze_atoms": macro_frozen}
        raw, _ = tsopt._freq_calc_full_hessian_torch(
            geom, macro_cfg, torch.device("cpu"), refresh_geom_meta=True,
        )
        h_analysis, active, _, _ = tsopt._reconcile_hessian_analysis_basis(
            raw, geom, [3, 5],
        )
        projection = {}
        frequencies, modes = tsopt._modes_from_Hact_embedded(
            h_analysis, geom.atomic_numbers, geom.cart_coords.reshape(-1, 3),
            active, torch.device("cpu"), tr_projection=geom.tr_projection,
            projection_info=projection, frequency_zero_cutoff_cm=DEFAULT_FREQUENCY_ZERO_CUTOFF_CM,
        )
        # The optimizer's terminal record includes its configured analysis cutoff.
        projection["frequency_zero_cutoff_cm"] = DEFAULT_FREQUENCY_ZERO_CUTOFF_CM
        permutation = [3, 4, 5, 0, 1, 2]
        optimizer = SimpleNamespace(
            is_converged=True, is_stalled=False, cur_cycle=0, stop_reason="",
            cur_H=raw[permutation][:, permutation].clone(),
            active_dof_indices=dofs(macro_order),
            _last_exact_cart_coords=geom.cart_coords.copy(),
            _last_exact_frequencies_cm=frequencies.copy(),
            _last_exact_modes=modes.clone(), _last_rigid_projection_info=projection,
            _last_exact_target_mode_index=0, _last_exact_target_mode_overlap=0.91,
        )
        assert np.count_nonzero(resolved_imaginary_mask(frequencies)) == 1
        assert projection["effective_rank"] == 0
        observed["optimizer"] = optimizer
        # Match the existing driver's restoration boundary; cached optimizer data
        # retains its macro basis while geometry/calculator regain the entry mask.
        geom.freeze_atoms = entry_frozen
        geom.set_calculator(QuadraticCalculator(**calc_cfg))
        return {"optimizer": optimizer, "converged": True, "cycles": 1,
                "safeguards": {}, "micro_cycles": 0, "outcome": None}

    def export_modes(_geom, frequencies, modes, *_args, **_kwargs):
        exported.append((frequencies.copy(), modes.clone()))
        return 0  # file-writing boundary only; do not manufacture a verdict

    monkeypatch.chdir(tmp_path)
    # Isolate the in-process raw store; its real store/load implementations run.
    monkeypatch.setattr(hessian_cache, "_cache", {})
    for module in (tsopt, freq, calc_eval):
        monkeypatch.setattr(module, "mlmm", QuadraticCalculator)
    source, parm, config = (tmp_path / name for name in ("input.pdb", "input.parm7", "cpu.yaml"))
    source.write_text("Input parsing is replaced at the prepared-structure boundary.\n")
    parm.write_text("No Amber parser or MM engine is constructed.\n")
    config.write_text("calc:\n  ml_device: cpu\nhessian_dimer:\n  device: cpu\n")
    prepared = SimpleNamespace(original_path=source, source_path=source,
                               geom_path=source, cleanup=lambda: None)
    monkeypatch.setattr(tsopt, "prepare_input_structure", lambda *_a: prepared)
    monkeypatch.setattr(tsopt, "resolve_charge_spin_or_raise", lambda *_a, **_k: (0, 1))
    monkeypatch.setattr(tsopt, "resolve_ml_layer_assignment", lambda **_k: (source, None))
    monkeypatch.setattr(tsopt, "geom_loader", load_geometry)
    monkeypatch.setattr(tsopt, "_run_microiter_tsopt", completed_macro)
    monkeypatch.setattr(tsopt, "_write_all_imag_modes", export_modes)
    energy_calls = []
    actual_energy = tsopt._calc_energy

    def terminal_energy(*args, **kwargs):
        energy_calls.append(True)
        if final_energy == "raises":
            raise RuntimeError("injected final energy failure")
        if final_energy == "nonfinite":
            return float("nan")
        return actual_energy(*args, **kwargs)

    monkeypatch.setattr(tsopt, "_calc_energy", terminal_energy)
    output = tmp_path / "terminal"
    result = CliRunner().invoke(tsopt.cli, [
        "-i", str(source), "--parm", str(parm), "-q", "0", "-m", "1",
        "--config", str(config), "--out-dir", str(output), "--opt-mode", "hess",
        "--microiter", "--active-dof-mode", "all", "--no-flatten",
        "--no-convert-files", "--no-dump", "--out-json",
        "--freeze-atoms", ",".join(str(atom + 1) for atom in final_frozen),
    ])
    assert result.exit_code == 0, result.output
    report = json.loads((output / "result.json").read_text())
    expected_count = 2 if widen and not soft_added_root else 1
    expected_strict = 2 if widen else 1  # independently fixed by diagonal[9,12]
    assert report["hessian_status"] == "completed", (report, result.output)
    assert report["hessian_error"] is None
    assert report["n_imaginary_modes"] == expected_count
    assert report["n_negative_modes"] == expected_strict
    assert report["saddle_validation"] == ("higher_order" if expected_count > 1 else "first_order")
    assert report["saddle_order_verified"] is (expected_count == 1)
    assert report["optimization_status"] == "converged"
    assert report["status"] == ("converged" if final_energy == "finite" else "energy_missing")
    assert report["energy_hartree"] == (0.0 if final_energy == "finite" else None)
    assert len(energy_calls) == 1
    assert report["reaction_mode_index"] == 0
    assert report["reaction_mode_overlap"] == (None if widen else 0.91)
    assert report["reaction_mode_source"] == ("lowest-imaginary" if widen else "mep-reference-overlap")
    assert report["imaginary_mode_criterion"] == "mass_weighted_eigenvalue"
    assert report["imaginary_eigenvalue_threshold"] == 1e-6
    projection = report["rigid_projection"]
    assert projection["active_atoms"] == final_atoms
    assert projection["frozen_atoms"] == final_frozen
    assert projection["effective_rank"] == 0
    expected_source = "tsopt_exact" if widen else "optimizer_terminal_exact_phva"
    assert projection["source"] == expected_source
    assert bool(projection.get("reused_without_hessian_recalculation", False)) is (not widen)
    assert len(calls) == 1 + int(widen)  # one macro H; zero/one terminal H
    assert len(exported) == 1
    frequencies, modes = exported[0]
    assert frequencies.size == 3 * len(final_atoms)
    assert np.count_nonzero(resolved_imaginary_mask(frequencies)) == expected_count
    assert torch.count_nonzero(modes[:, dofs(final_frozen)]).item() == 0
    geom = observed["geometry"]
    assert list(geom.freeze_atoms) == final_frozen
    np.testing.assert_array_equal(geom.cart_coords, x0)
    raw_cache = hessian_cache.load("ts")
    expected_order = [4, 5, 3] if widen else macro_order
    expected_dofs = dofs(expected_order)
    assert raw_cache["active_dofs"] == expected_dofs
    assert raw_cache["meta"]["source"] == expected_source
    assert torch.equal(raw_cache["hessian"], full_h[expected_dofs][:, expected_dofs])
    np.testing.assert_array_equal(raw_cache["meta"]["cart_coords"], x0)
    # Raw ownership is checked after the actual analysis, not against its output.
    macro_dofs = dofs(macro_order)
    assert torch.equal(observed["optimizer"].cur_H, full_h[macro_dofs][:, macro_dofs])
    if widen:
        assert calls[-1]["order"] == expected_order
        assert torch.equal(raw_cache["hessian"], calls[-1]["raw"])
