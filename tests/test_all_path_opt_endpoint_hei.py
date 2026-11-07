"""``all`` selects an internal MEP peak and measures its preceding-basin barrier.

Profiles without an internal peak must produce no barrier or TS dispatch.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from ase import Atoms
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.core import utils
from mlmm.core.result_commit import MLMM_RUN_ID_ENV, with_current_run_id
from mlmm.workflows import all as all_workflow
from mlmm.workflows import dft
from test_all_tsopt_bond_diagnostic import ReachedTSDispatch, _coords, _pdb, _write


def _frame(index, energy):
    atoms = "".join(
        f"{element} {x:.3f} {y:.3f} {z:.3f}\n"
        for element, (x, y, z) in zip(("C", "O"), _coords(index))
    )
    return f"2\nE={energy:.3f} unit=hartree frame{index}\n{atoms}"


@pytest.mark.parametrize(
    "energies,hei,barrier_au",
    [
        pytest.param((0.0, 0.1, 0.02), 1, 0.1, id="interior-hei-control"),
        pytest.param((0.1, 0.05, 0.0), 0, None, id="reactant-end-hei"),
        pytest.param((0.0, 0.05, 0.1), 2, None, id="product-end-hei"),
        pytest.param((0.0, 0.1, 0.02, 0.15), 1, 0.1, id="peak-below-product"),
        pytest.param((0.04, 0.0, 0.1, 0.02, 0.15), 2, 0.1, id="preceding-basin"),
    ],
)
def test_path_opt_branch_builds_no_segment_for_an_endpoint_hei(tmp_path, monkeypatch, energies, hei, barrier_au):
    monkeypatch.delenv(MLMM_RUN_ID_ENV, raising=False)
    monkeypatch.setattr(utils, "_CONVERT_FILES_ENABLED", True)
    out = tmp_path / "out"
    inputs = [_write(tmp_path / f"input{i}.pdb", _pdb(i)) for i in (0, 2)]
    # Only initial ML/MM workspace construction is replaced. This fixture is
    # not a force field and must never reach a calculator or topology reader.
    parm = _write(tmp_path / "unused.parm7", "dispatch fixture; not a force field\n")
    model = Atoms("CO", positions=_coords(0))
    workspace = SimpleNamespace(
        atoms_model=model, atoms_model_lh=model.copy(), model_pdb=inputs[0],
        link_pairs=[], cleanup=lambda: None,
    )
    monkeypatch.setattr(dft, "_prepare_ml_region_workspace", lambda **kwargs: workspace)
    ts_calls = []

    def child(name, _cli, args, **kwargs):
        assert name == "path_opt"
        child_out = Path(args[args.index("--out-dir") + 1])
        _write(child_out / "final_geometries_trj.xyz", "".join(_frame(i, e) for i, e in enumerate(energies)))
        _write(child_out / "hei.xyz", _frame(hei, energies[hei]))
        _write(child_out / "hei.pdb", _pdb(hei))
        payload = {
            "optimization_status": "converged", "converged": True,
            "preopt_requested": False, "preopt_converged": None,
            "path_optimizers": [],
            "stage_outcomes": [{"stage": "path-opt", "item_id": "gsm_mep",
                                "executed": True, "converged": True,
                                "usable": True, "required": True}],
        }
        _write(child_out / "result.json", json.dumps(with_current_run_id(payload)))

    def stop_at_ts(hei_path, *args, **kwargs):
        ts_calls.append(Path(hei_path))
        raise ReachedTSDispatch

    def diagram(prefix, *, labels, energies_eh, **kwargs):
        return {"name": prefix.stem, "labels": labels, "energies_au": list(energies_eh)}

    def no_calculator(*args, **kwargs):
        raise AssertionError("The TS entry sentinel must precede calculator construction")

    monkeypatch.setattr(all_workflow, "_run_cli_main", child)
    monkeypatch.setattr(all_workflow, "_mlmm_calc", no_calculator)
    monkeypatch.setattr(all_workflow, "run_trj2fig", lambda *args, **kwargs: None)
    monkeypatch.setattr(all_workflow, "close_matplotlib_figures", lambda: None)
    monkeypatch.setattr(all_workflow, "_write_segment_energy_diagram", diagram)
    monkeypatch.setattr(all_workflow._path_search, "_has_bond_change", lambda *a, **k: (True, "Bond formed"))
    monkeypatch.setattr(all_workflow, "_run_tsopt_on_hei", stop_at_ts)
    args = ["all", *[arg for path in inputs for arg in ("-i", str(path))]]
    args += ["-q", "0", "-m", "1", "--parm", str(parm), "--model-pdb", str(inputs[0]),
             "--out-dir", str(out), "--no-preopt", "--convert-files",
             "--no-tsopt-from-mep-tan", "--tsopt"]

    if 0 < hei < len(energies) - 1:
        with pytest.raises(ReachedTSDispatch):
            CliRunner().invoke(root_cli, args, catch_exceptions=False)
        assert len(ts_calls) == 1 and ts_calls[0].stem == "hei_seg_01"
        summary = json.loads((out / "_work/path_opt/summary.json").read_text())
        segment, = summary["segments"]
        from pysisyphus.constants import AU2KCALPERMOL
        assert segment["barrier_kcal"] == pytest.approx(barrier_au * AU2KCALPERMOL)
        assert segment["index"] == 1 and "stage_outcomes" not in summary
        diagram, = summary["energy_diagrams"]
        assert diagram["labels"] == ["R", "TS1", "P"]
        assert diagram["energies_au"] == pytest.approx([energies[0], energies[hei], energies[-1]])
        return

    # The interval is not a reaction step, so the run ends as a failed result.
    result = CliRunner().invoke(root_cli, args)
    assert result.exit_code == 1, result.output + repr(result.exception)
    assert ts_calls == []
    assert f"[seg_01] WARNING: HEI is at an endpoint (idx={hei}). Returning the raw path." in result.output
    summary = json.loads((out / "summary.json").read_text())
    assert summary["segments"] == [] and summary["n_segments"] == 0 and summary["n_images"] == 3
    assert summary["execution_status"] == "completed"
    assert summary["scientific_status"] == "failed"
    raw, = [leaf for leaf in summary["stage_outcomes"] if leaf["item_id"] == "raw_seg_01"]
    assert raw["stage"] == "path" and raw["usable"] is False and raw["reason"] == "endpoint_hei"
    assert "rate_limiting_step" not in summary
    # The MEP diagram shows no TS for the interval.
    diagram, = summary["energy_diagrams"]
    assert diagram["labels"] == ["R", "P"]
    assert diagram["energies_au"] == pytest.approx([energies[0], energies[-1]])


@pytest.mark.parametrize(
    "second,labels",
    [
        pytest.param((0.02, 0.08, 0.01), ["R", "TS1", "IM1_1", "IM1_2", "TS2", "P"], id="two-interior-heis"),
        pytest.param((0.02, 0.01, 0.0), ["R", "TS1", "IM1_1", "IM1_2", "P"], id="second-endpoint-hei"),
    ],
)
def test_path_opt_mep_diagram_covers_every_interval(tmp_path, monkeypatch, second, labels):
    monkeypatch.delenv(MLMM_RUN_ID_ENV, raising=False)
    monkeypatch.setattr(utils, "_CONVERT_FILES_ENABLED", True)
    out = tmp_path / "out"
    inputs = [_write(tmp_path / f"input{i}.pdb", _pdb(i)) for i in (0, 2, 4)]
    parm = _write(tmp_path / "unused.parm7", "dispatch fixture; not a force field\n")
    model = Atoms("CO", positions=_coords(0))
    workspace = SimpleNamespace(
        atoms_model=model, atoms_model_lh=model.copy(), model_pdb=inputs[0],
        link_pairs=[], cleanup=lambda: None,
    )
    monkeypatch.setattr(dft, "_prepare_ml_region_workspace", lambda **kwargs: workspace)
    profiles = {"seg_01_mep": (0.0, 0.1, 0.02), "seg_02_mep": second}

    def child(name, _cli, args, **kwargs):
        assert name == "path_opt"
        child_out = Path(args[args.index("--out-dir") + 1])
        energies = profiles[child_out.name]
        start = 2 * (int(child_out.name[4:6]) - 1)
        _write(child_out / "final_geometries_trj.xyz",
               "".join(_frame(start + i, e) for i, e in enumerate(energies)))
        _write(child_out / "hei.xyz", _frame(start + 1, energies[1]))
        _write(child_out / "hei.pdb", _pdb(start + 1))
        payload = {
            "optimization_status": "converged", "converged": True,
            "preopt_requested": False, "preopt_converged": None,
            "path_optimizers": [],
            "stage_outcomes": [{"stage": "path-opt", "item_id": "gsm_mep",
                                "executed": True, "converged": True,
                                "usable": True, "required": True}],
        }
        _write(child_out / "result.json", json.dumps(with_current_run_id(payload)))

    def stop_at_ts(*args, **kwargs):
        raise ReachedTSDispatch

    def diagram(prefix, *, labels, energies_eh, **kwargs):
        return {"name": prefix.stem, "labels": labels, "energies_au": list(energies_eh)}

    monkeypatch.setattr(all_workflow, "_run_cli_main", child)
    monkeypatch.setattr(all_workflow, "run_trj2fig", lambda *args, **kwargs: None)
    monkeypatch.setattr(all_workflow, "close_matplotlib_figures", lambda: None)
    monkeypatch.setattr(all_workflow, "_write_segment_energy_diagram", diagram)
    monkeypatch.setattr(all_workflow._path_search, "_has_bond_change", lambda *a, **k: (True, "Bond formed"))
    monkeypatch.setattr(all_workflow, "_run_tsopt_on_hei", stop_at_ts)
    args = ["all", *[arg for path in inputs for arg in ("-i", str(path))]]
    args += ["-q", "0", "-m", "1", "--parm", str(parm), "--model-pdb", str(inputs[0]),
             "--out-dir", str(out), "--no-preopt", "--convert-files",
             "--no-tsopt-from-mep-tan", "--tsopt"]

    with pytest.raises(ReachedTSDispatch):
        CliRunner().invoke(root_cli, args, catch_exceptions=False)
    summary = json.loads((out / "_work/path_opt/summary.json").read_text())
    diagram, = summary["energy_diagrams"]
    assert diagram["labels"] == labels
    expected = [0.0, 0.1, 0.02, *(second if len(labels) == 6 else (second[0], second[-1]))]
    assert diagram["energies_au"] == pytest.approx(expected)
