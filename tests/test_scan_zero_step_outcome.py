"""Real scan CLI producer; ML/MM calculator and optimizer execution are replaced.

All cases consume the same finite physical input and request its current bond
distance. The scheduling, result files and outcome aggregation stay real.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pytest
from click.testing import CliRunner
from pysisyphus.constants import BOHR2ANG
from pysisyphus.helpers import geom_loader

from mlmm.cli import cli as root_cli
from mlmm.core.utils import prepare_input_structure
from mlmm.workflows import scan as workflow


@pytest.mark.parametrize(
    "endopt_requested,endopt_converged,finite_energy,expected_usable",
    [
        pytest.param(False, None, True, True, id="no-optimizer-finite"),
        pytest.param(False, None, False, False, id="no-optimizer-nan-energy"),
        pytest.param(True, True, True, True, id="endopt-converged"),
        pytest.param(True, False, True, False, id="endopt-not-converged"),
        pytest.param(True, None, True, False, id="endopt-convergence-unknown"),
        pytest.param(True, True, False, False, id="endopt-converged-nan-energy"),
    ],
)
def test_scan_zero_step_cli_outcome(
    tmp_path, monkeypatch, endopt_requested, endopt_converged,
    finite_energy, expected_usable,
):
    coords = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.]])
    source = tmp_path / "three_carbons.pdb"
    source.write_text("".join(
        f"HETATM{i:5d}  C{i:<2d} MOL A   1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  \n"
        for i, (x, y, z) in enumerate(coords, 1)
    ) + "END\n")
    original_bytes = source.read_bytes()
    # Use the same real PDB -> Cartesian -> Angstrom conversion as scan. This
    # prevents a rounding-sized requested displacement from becoming one step.
    with prepare_input_structure(source) as prepared:
        initial_geom = geom_loader(prepared.geom_path, coord_type="cart")
        initial_bohr = np.asarray(initial_geom.cart_coords).reshape(-1, 3).copy()
    initial_angstrom = initial_bohr * BOHR2ANG
    target = float(np.linalg.norm(initial_angstrom[0] - initial_angstrom[1]))
    energy = -1.0 if finite_energy else float("nan")
    evaluations = []
    optimizer_calls = []
    optimizer_runs = []

    class ConstantCalculator:
        freeze_atoms = []
        analytical_2d = False

        def get_energy(self, atoms, positions, **kwargs):
            evaluations.append(np.asarray(positions).reshape(-1, 3).copy())
            return {"energy": energy}

        def get_forces(self, atoms, positions, **kwargs):
            return {"energy": energy, "forces": np.zeros(np.asarray(positions).size)}

    def optimizer(geometry, *args, **kwargs):
        optimizer_calls.append(kwargs.get("prefix"))
        np.testing.assert_allclose(geometry.cart_coords.reshape(-1, 3), initial_bohr)
        return SimpleNamespace(
            run=lambda: optimizer_runs.append(True),
            is_converged=endopt_converged,
            is_stalled=False,
            stop_reason=("" if endopt_converged is True else
                         "max_cycles" if endopt_converged is False else "unknown"),
        )

    monkeypatch.setattr(workflow, "mlmm", lambda **kwargs: ConstantCalculator())
    monkeypatch.setattr(workflow, "LBFGS", optimizer)
    monkeypatch.setattr(workflow, "_has_bond_change", lambda *args, **kwargs: (False, ""))
    output = tmp_path / "scan"
    # Explicit model indices prevent layer inference. The calculator boundary
    # is replaced, so this existing path is never parsed as an Amber topology.
    result = CliRunner().invoke(root_cli, [
        "scan", "-i", str(source), "-q", "0", "-m", "1",
        "--parm", str(source), "--model-indices", "1-3", "--no-detect-layer",
        "--no-preopt", "--endopt" if endopt_requested else "--no-endopt",
        "--no-convert-files", "--out-json", "--one-based",
        "--scan-lists", f"[(1,2,{target!r})]", "-o", str(output),
    ])
    assert result.exit_code == 0, result.output + repr(result.exception)
    assert source.read_bytes() == original_bytes
    assert optimizer_calls == (["endopt"] if endopt_requested else [])
    assert len(optimizer_runs) == int(endopt_requested)
    assert evaluations
    for actual in evaluations:
        np.testing.assert_allclose(actual, initial_bohr, rtol=0., atol=1e-12)

    final_xyz = output / "stage_01/result.xyz"
    assert final_xyz.is_file()
    final_geom = geom_loader(final_xyz, coord_type="cart")
    np.testing.assert_allclose(final_geom.cart_coords.reshape(-1, 3), initial_bohr, atol=1e-6)
    payload = json.loads((output / "result.json").read_text())
    assert payload["status"] == "completed"
    assert payload["n_stages"] == 1
    stage, = payload["stages"]
    leaf, = payload["stage_outcomes"]
    assert stage["n_steps"] == 0
    assert stage["energies_hartree"] == []
    assert stage["converged"] is (endopt_converged if endopt_requested else None)
    assert stage["final_energy_hartree"] == (-1.0 if finite_energy else None)
    assert leaf["item_id"] == "stage_1"
    assert leaf["executed"] is True
    assert leaf["converged"] is (endopt_converged if endopt_requested else None)
    assert leaf["usable"] is expected_usable
    assert (payload["scientific_status"] == "success") is expected_usable
    if not finite_energy:
        assert leaf["reason"] == "energy_invalid"
    elif not endopt_requested:
        assert leaf["reason"] == "no_optimization_requested"
    elif endopt_converged is not True:
        assert leaf["reason"] == ("not_converged" if endopt_converged is False
                                   else "convergence_unknown")
