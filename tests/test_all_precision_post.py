"""``all --precision-post`` sets the precision of the stages after the MEP.

TSOPT, IRC, the post-IRC endpoint optimizations and freq (children and their
parent-side calculators) use ``--precision-post`` when it is given and
``--precision`` when it is not; scan, path-opt and path-search always use
``--precision``. Children are stubbed at the dispatcher and the parent
calculator is a recording stub, so no backend runs.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml
from ase import Atoms
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.core import utils
from mlmm.core.result_commit import MLMM_RUN_ID_ENV, with_current_run_id
from mlmm.workflows import all as all_workflow
from mlmm.workflows import dft
from test_all_tsopt_bond_diagnostic import _coords, _frame, _pdb, _write


class ReachedChild(BaseException):
    """Escape the parent pipeline at a child dispatch."""


def _child_precision(args):
    """Precision a child resolves: its --precision, else the UMA key of --config."""
    if "--precision" in args:
        return args[args.index("--precision") + 1]
    config = yaml.safe_load(
        Path(args[args.index("--config") + 1]).read_text(encoding="utf-8")
    )
    return config["calc"]["uma_precision"]


_TS_RESULT = {
    "optimization_status": "converged", "converged": True,
    "execution_status": "completed", "scientific_status": "success",
    "hessian_status": "completed", "saddle_validation": "first_order",
    "n_imaginary_modes": 1, "imaginary_frequencies_cm": [-500.0],
    "reaction_mode_index": 0, "reaction_mode_frequency_cm": -500.0,
    "energy_hartree": 0.1,
}


class _EnergyOnly:
    """Parent-side calculator stub: a finite energy, nothing else."""

    def get_energy(self, _atoms, _coords):
        return {"energy": -1.0}

    def close(self):
        pass


def _invoke(tmp_path, monkeypatch, n_inputs, extra, precision_post, child):
    monkeypatch.delenv(MLMM_RUN_ID_ENV, raising=False)
    monkeypatch.setattr(utils, "_CONVERT_FILES_ENABLED", True)
    inputs = [_write(tmp_path / f"input{i}.pdb", _pdb(i)) for i in (0, 2)[:n_inputs]]
    # Only the ML/MM workspace construction is replaced; this parm7 fixture is
    # not a force field and must never reach a calculator or topology reader.
    parm = _write(tmp_path / "unused.parm7", "dispatch fixture; not a force field\n")
    model = Atoms("CO", positions=_coords(0))
    workspace = SimpleNamespace(
        atoms_model=model, atoms_model_lh=model.copy(), model_pdb=inputs[0],
        link_pairs=[], cleanup=lambda: None,
    )
    monkeypatch.setattr(dft, "_prepare_ml_region_workspace", lambda **kwargs: workspace)
    monkeypatch.setattr(all_workflow, "_run_cli_main", child)
    args = ["all", *[arg for path in inputs for arg in ("-i", str(path))]]
    args += ["-q", "0", "-m", "1", "--parm", str(parm), "--model-pdb", str(inputs[0]),
             "--out-dir", str(tmp_path / "out"), "--no-preopt",
             "--precision", "tf32", *extra]
    if precision_post is not None:
        args += ["--precision-post", precision_post]
    with pytest.raises(ReachedChild):
        CliRunner().invoke(root_cli, args, catch_exceptions=False)


def _run(tmp_path, monkeypatch, n_inputs, extra, precision_post):
    calls, calculators = [], []

    def child(name, _cli, args, **kwargs):
        args = list(args)
        # Read now: the child YAML is removed when the run unwinds.
        calls.append((name, args, _child_precision(args)))
        out = Path(args[args.index("--out-dir") + 1])
        if name == "path_opt":
            _write(out / "final_geometries_trj.xyz", "".join(_frame(i) for i in range(3)))
            _write(out / "hei.xyz", _frame(1))
            _write(out / "hei.pdb", _pdb(1))
            _write(out / "result.json", json.dumps(with_current_run_id({
                "optimization_status": "converged", "converged": True,
                "preopt_requested": False, "preopt_converged": None,
                "path_optimizers": [],
                "stage_outcomes": [{"stage": "path-opt", "item_id": "gsm_mep",
                                    "executed": True, "converged": True,
                                    "usable": True, "required": True}],
            })))
        elif name == "tsopt":
            _write(out / "final_geometry.xyz", _frame(1))
            _write(out / "final_geometry.pdb", _pdb(1))
            _write(out / "result.json", json.dumps(with_current_run_id(_TS_RESULT)))
        elif name == "irc":
            trajectory = _write(out / "finished_irc_trj.xyz", _frame(0) + _frame(2))
            _write(out / "result.json", json.dumps(with_current_run_id(
                {"files": {"finished_irc": str(trajectory)}}
            )))
        elif name == "opt":
            _write(out / "final_geometry.xyz", _frame(0))
            _write(out / "final_geometry.pdb", _pdb(0))
            _write(out / "result.json", json.dumps(with_current_run_id(
                {"optimization_status": "converged", "converged": True}
            )))
        else:
            raise ReachedChild
        return 0

    def calculator(**cfg):
        calculators.append(dict(cfg))
        return _EnergyOnly()

    monkeypatch.setattr(all_workflow, "_mlmm_calc", calculator)
    monkeypatch.setattr(all_workflow, "run_trj2fig", lambda *args, **kwargs: None)
    monkeypatch.setattr(all_workflow, "close_matplotlib_figures", lambda: None)
    monkeypatch.setattr(
        all_workflow, "_write_segment_energy_diagram",
        lambda prefix, *, labels, energies_eh, **kwargs: {
            "name": prefix.stem, "labels": labels, "energies_au": list(energies_eh),
        },
    )
    monkeypatch.setattr(
        all_workflow._path_search, "_has_bond_change",
        lambda *args, **kwargs: (True, "Bond formed"),
    )
    _invoke(tmp_path, monkeypatch, n_inputs,
            ["--no-tsopt-from-mep-tan", *extra], precision_post, child)
    return calls, calculators


_PRE_MEP = {
    "path_opt": (2, []),
    "path_search": (2, ["--refine-path"]),
    "scan": (1, ["--scan-lists", "[(1,2,1.50)]"]),
}


@pytest.mark.parametrize("precision_post", [None, "fp64"], ids=["post-unset", "post-fp64"])
@pytest.mark.parametrize("child_name", list(_PRE_MEP))
def test_pre_mep_children_keep_precision(tmp_path, monkeypatch, child_name, precision_post):
    n_inputs, extra = _PRE_MEP[child_name]
    captured = []

    def first_child(name, _cli, args, **kwargs):
        captured.append((name, list(args), _child_precision(list(args))))
        raise ReachedChild

    def no_calculator(*args, **kwargs):
        raise AssertionError("No calculator may be built before the child dispatch")

    monkeypatch.setattr(all_workflow, "_mlmm_calc", no_calculator)
    _invoke(tmp_path, monkeypatch, n_inputs, extra, precision_post, first_child)

    (name, child_args, resolved), = captured
    assert name == child_name
    assert "--precision" not in child_args
    assert resolved == "tf32"


@pytest.mark.parametrize("precision_post", [None, "fp64"], ids=["post-unset", "post-fp64"])
def test_tsopt_only_post_stages_follow_precision_post(tmp_path, monkeypatch, precision_post):
    expected = precision_post or "tf32"
    calls, calculators = _run(
        tmp_path, monkeypatch, 1, ["--tsopt", "--thermo"], precision_post
    )

    assert [name for name, _, _ in calls] == ["tsopt", "irc", "opt", "opt", "freq"]
    for _name, child_args, resolved in calls:
        assert resolved == expected
        assert ("--precision" in child_args) is (precision_post is not None)
    # TS probe, IRC endpoint calculator, and the two endpoint probes.
    assert [cfg["uma_precision"] for cfg in calculators] == [expected] * 4


@pytest.mark.parametrize("precision_post", [None, "fp64"], ids=["post-unset", "post-fp64"])
def test_mep_post_stages_follow_precision_post(tmp_path, monkeypatch, precision_post):
    expected = precision_post or "tf32"
    calls, calculators = _run(
        tmp_path, monkeypatch, 2, ["--tsopt", "--thermo"], precision_post
    )

    assert [name for name, _, _ in calls] == ["path_opt", "tsopt", "irc", "opt", "opt", "freq"]
    (_, mep_args, mep_resolved), *post_calls = calls
    assert "--precision" not in mep_args
    assert mep_resolved == "tf32"
    for _name, child_args, resolved in post_calls:
        assert resolved == expected
        assert ("--precision" in child_args) is (precision_post is not None)
    assert [cfg["uma_precision"] for cfg in calculators] == [expected] * 4


def test_summary_records_post_precision(tmp_path):
    from mlmm.io.summary import write_summary_log

    for post, recorded in ((None, "tf32"), ("fp64", "fp64")):
        summary = {"segments": []}
        all_workflow._enrich_summary(
            summary, pipeline_mode="path-opt", out_dir=tmp_path,
            mlip_backend="uma", mlip_model="uma-s-1p2p1", mlip_precision="tf32",
            mlip_precision_post=post, charge=0, spin=1,
            config={"tsopt": True, "thermo": False, "dft": False},
        )
        assert summary["mlip_precision"] == "tf32"
        assert summary["mlip_precision_post"] == recorded

        dest = tmp_path / f"summary_{recorded}.log"
        write_summary_log(dest, {"root_out_dir": str(tmp_path), **summary})
        text = dest.read_text(encoding="utf-8")
        assert "MLIP precision      : tf32" in text
        assert ("MLIP post precision : fp64" in text) is (post is not None)
