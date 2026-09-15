"""Run the real frequency CLI through eigenpairs, files and QRRHO."""
import importlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import torch
import yaml
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.workflows import freq as workflow
from thermoanalysis.config import WORKFLOW_THERMO_POLICY


@pytest.mark.parametrize("partial", [False, True])
def test_frequency_cli_retains_soft_modes_and_thermal_input(tmp_path, monkeypatch, partial):
    # Evaluator boundary only is replaced: real structure preparation, active
    # projection, diagonalization, mode writers, JSON/YAML and thermal code run.
    coords = np.array([[0., 0., 0.], [2., 0., 0.], [0., 2., 0.],
                       [0., 0., 2.], [2., 2., 2.]])
    source = tmp_path / "five_carbons.pdb"
    source.write_text("".join(
        f"HETATM{i:5d}  C{i:<2d} MOL A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  \n"
        for i, (x, y, z) in enumerate(coords, 1)
    ) + "END\n")
    # Carbon masses are equal; each of the six free Cartesian eigenpairs is
    # distinct. Three non-collinear frozen anchors give rigid rank zero.
    diagonal = torch.tensor([-1e-3, -1e-9, 0., 1e-9, 1e-8, 1e-3], dtype=torch.float64)
    active_h = torch.diag(diagonal)
    full_h = torch.eye(15, dtype=torch.float64)
    full_h[9:, 9:] = active_h
    evaluations = []

    def evaluate(geometry, config, device, **kwargs):
        assert list(geometry.freeze_atoms) == [0, 1, 2]
        evaluations.append(geometry.cart_coords.copy())
        geometry.within_partial_hessian = ({
            "active_n_dof": 6, "full_n_dof": 15,
            "active_atoms": [3, 4], "active_dofs": list(range(9, 15)),
        } if partial else None)
        hessian = (active_h if partial else full_h).clone()
        return hessian, -1.

    monkeypatch.setattr(workflow, "_calc_full_hessian_torch", evaluate)
    monkeypatch.setattr(workflow, "mlmm", lambda **kwargs: SimpleNamespace(
        ml_indices=list(range(5)), hess_mm_indices=[], movable_mm_indices=[],
        frozen_layer_indices=[],
    ))
    cache = importlib.import_module("mlmm.io.hessian_cache")
    monkeypatch.setattr(cache, "load_matching", lambda *args, **kwargs: None)
    thermal_module = importlib.import_module("thermoanalysis.thermo")
    original_thermo = thermal_module.thermochemistry
    observed = []

    def observe(qc, *args, **kwargs):
        result = original_thermo(qc, *args, **kwargs)
        observed.append((qc.wavenumbers.copy(), np.asarray(result.wavenumbers).copy(),
                         np.asarray(qc.coords3d).copy()))
        return result

    monkeypatch.setattr(thermal_module, "thermochemistry", observe)
    payloads = []
    for cutoff, sort in [(0., "value"), (5., "abs"), (7., "value")]:
        output = tmp_path / f"freq_{cutoff:g}"
        config = tmp_path / f"config_{cutoff:g}.yaml"
        config.write_text(yaml.safe_dump({"freq": {"zero_cutoff_cm": cutoff},
                                         "thermo": {"symmetry_number": 1}}))
        extra = ["--parm", str(source), "--model-indices", "1-5", "--no-detect-layer", "--active-dof-mode", "all", "--hess-device", "cpu"]
        result = CliRunner().invoke(root_cli, [
            "freq", "-i", str(source), "-q", "0", "-m", "1",
            "--freeze-atoms", "1,2,3", "--max-write", "6", "--n-frames", "4",
            "--sort", sort, "--no-convert-files", "--dump", "--out-json",
            "--config", str(config), "-o", str(output), *extra,
        ])
        assert result.exit_code == 0, result.output + repr(result.exception)
        payload = json.loads((output / "result.json").read_text())
        summary = yaml.safe_load((output / "thermoanalysis.yaml").read_text())
        raw = np.asarray(payload["frequencies_cm"])
        assert raw.shape == (6,)
        assert raw[1] < 0. < raw[3] < raw[4] < 5. and raw[2] == 0.
        assert payload["frequency_representation"] == "complete"
        assert payload["n_modes"] == 6 and payload["n_negative_modes"] == 2
        assert payload["n_imaginary"] == (2 if cutoff == 0. else 1)
        assert summary["num_imag_freq"] == payload["n_imaginary"]
        assert summary["n_negative_modes"] == 2
        assert summary["thermo_policy"] == WORKFLOW_THERMO_POLICY.as_dict()
        assert payload["rigid_projection"]["effective_rank"] == 0
        np.testing.assert_array_equal(observed[-1][0], raw)
        np.testing.assert_array_equal(observed[-1][1], raw[raw > 0.])
        # QCData centers/rotates to principal axes; pair distances retain the
        # Angstrom unit and atom identity without assuming that orientation.
        thermal_coords = observed[-1][2]
        np.testing.assert_allclose(
            np.linalg.norm(thermal_coords[:, None] - thermal_coords[None, :], axis=-1),
            np.linalg.norm(coords[:, None] - coords[None, :], axis=-1), atol=1e-6,
        )
        table = np.loadtxt(output / "frequencies_cm-1.txt")
        np.testing.assert_allclose(np.sort(table[:, 1]), raw, atol=5.1e-5)
        trajectories = [output / p for p in payload["files"]["mode_files"] if p.endswith("_trj.xyz")]
        assert len(trajectories) == 6 and all(p.is_file() for p in trajectories)
        # The exported low-positive mode must move an active atom, preserving
        # the frequency/vector pairing; the anchored coordinates stay fixed.
        small_positive = next(p for p in trajectories if f"{raw[3]:+.2f}cm-1" in p.name)
        from ase.io import read
        frames = read(small_positive, index=":")
        assert len(frames) == 4
        for frame in frames:
            np.testing.assert_allclose(frame.positions[:3], coords[:3], atol=1e-5)
        assert max(np.linalg.norm(frame.positions[3:] - coords[3:]) for frame in frames) > .01
        payloads.append(payload)
    assert len(evaluations) == len(observed) == 3
    for payload in payloads[1:]:
        np.testing.assert_array_equal(payload["frequencies_cm"], payloads[0]["frequencies_cm"])
        assert payload["thermochemistry"] == payloads[0]["thermochemistry"]
