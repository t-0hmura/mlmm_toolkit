"""Endpoint hard failures stop real all-branch control flow before consumers.

The endpoint-boundary tests run the whole ``mlmm all`` callback (TS-only and
MEP routes) with the child commands stubbed at the dispatcher
(``_all_stub_harness.py``); exception handlers, retention, the stage gate and
cleanup execute unchanged.
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from _all_stub_harness import Behaviour, coords, run_all, write
from mlmm.workflows import all as workflow
from mlmm.io import hessian_cache

def test_failed_ts_only_summary_does_not_invent_a_barrier(tmp_path):
    stop = {"stage": "endpoint_opt", "reason": "endpoint_execution_failed"}
    summary = {
        "segments": [{"index": 1, "tag": "seg_01", "kind": "tsopt"}],
        "pipeline_stop": stop, "energy_diagrams": [],
    }
    post = [{
        "index": 1, "kind": "tsopt", "pipeline_stop": stop,
        "tsopt": {"continue_irc": True, "saddle_validation": "first_order",
                  "optimization_status": "converged", "hessian_status": "completed",
                  "n_imaginary_modes": 1, "energy_valid": True,
                  "structure_valid": True},
        "irc": {"usable": True, "reason": "ok"}, "irc_traj": "irc.xyz",
        "endpoint_opt": {
            "reactant_converged": None,
            "product_converged": True,
            "failures": {"reactant": {"error_type": "ValueError", "error": "boom"}},
        },
    }]
    workflow._enrich_summary(
        summary, pipeline_mode="tsopt-only", out_dir=tmp_path,
        mlip_backend="orb", mlip_model="orb_v3_conservative_omol",
        charge=0, spin=1, post_segments=post,
        config={"tsopt": True, "thermo": True, "dft": True},
    )
    assert summary["scientific_status"] == "partial"
    assert summary["execution_status"] == "failed"
    assert any(
        "reactant_execution_failed" in reason
        for reason in summary["scientific_status_reasons"]
    )
    assert "rate_limiting_step" not in summary


def test_valid_ts1_with_one_nonconverged_endpoint_is_partial_not_failed(tmp_path):
    summary = {
        "segments": [{"index": 1, "kind": "tsopt"}],
        "energy_diagrams": [],
    }
    post = [{
        "index": 1,
        "kind": "tsopt",
        "tsopt": {
            "optimization_status": "converged",
            "hessian_status": "completed",
            "n_imaginary_modes": 1,
            "energy_valid": True,
            "structure_valid": True,
        },
        "endpoint_opt": {
            "reactant_converged": True,
            "product_converged": False,
        },
    }]

    truth = workflow._pipeline_aggregate_truth(
        summary,
        post_segments=post,
        config={"tsopt": True},
        legacy_status="partial",
    )

    assert truth.scientific_status == "partial"
    assert truth.execution_status == "completed"
    assert "all:segment_1:endpoint_opt:product_not_converged" in truth.status_reasons


class Geometry:
    def __init__(self, coordinate):
        self.cart_coords = np.array([coordinate, 0.0, 0.0])
        self.energy = -1.0
        self.calculator = None
        self.freeze_atoms = []

    def set_calculator(self, calculator):
        self.calculator = calculator

    def as_xyz(self):
        return f"1\nobserved\nH {self.cart_coords[0]} 0 0\n"


_ROLE = {0: "R", 1: "P"}
_STATE = {"R": "reactant", "P": "product"}
_MODE = {
    "hard_error": "error", "nonfinite_error": "nonfinite",
    "converged": "converged", "not_converged": "not_converged",
}


def _xyz_values(path: Path):
    rows = Path(path).read_text(encoding="utf-8").splitlines()[2:4]
    return [tuple(float(v) for v in row.split()[1:4]) for row in rows]


def _record_irc_leases(monkeypatch, leases):
    real = workflow._irc_and_match

    def irc_and_match(*args, **kwargs):
        result = real(*args, **kwargs)
        leases.append(result["calculator_lease"])
        return result

    return lambda _hooks: monkeypatch.setattr(workflow, "_irc_and_match", irc_and_match)


@pytest.mark.parametrize("branch", ["ts", "seg"])
@pytest.mark.parametrize("failed_endpoint", [0, 1])
@pytest.mark.parametrize("dump", [False, True])
@pytest.mark.parametrize("outcome", ["hard_error", "nonfinite_error", "converged", "not_converged"])
def test_endpoint_boundary_retains_provenance_and_stops_consumers(
    monkeypatch, tmp_path, branch, failed_endpoint, dump, outcome,
):
    role = _ROLE[failed_endpoint]
    leases = []
    run = run_all(
        tmp_path, monkeypatch,
        frames=(1.0,) if branch == "ts" else (0.0, 2.0),
        extra=["--tsopt", "--thermo", *(["--dump"] if dump else [])],
        behaviour=Behaviour(opt={role: _MODE[outcome]}),
        wrap=_record_irc_leases(monkeypatch, leases),
    )
    failed = outcome in {"hard_error", "nonfinite_error"}
    assert run.result.exit_code == (1 if failed else 0), run.result.output[-3000:]
    segment = run.out / "segments" / "seg_01"
    structures = segment / "structures"
    endpoint_dir = segment / "endpoint_opt"

    # Both endpoints are attempted even when the first one fails.
    assert [Path(args[args.index("--out-dir") + 1]).name for args in run.child_args("opt")] == ["R", "P"]
    assert (structures / "ts.xyz").is_file()
    # TS-only presents the higher-energy IRC endpoint (frame 2) as the reactant;
    # the MEP route keeps the MEP orientation (frame 0 is the reactant).
    irc_frames = {"reactant": 2.0, "product": 0.0} if branch == "ts" else {"reactant": 0.0, "product": 2.0}
    for state, frame in irc_frames.items():
        np.testing.assert_allclose(_xyz_values(structures / f"{state}_irc.xyz"), coords(frame), atol=1e-6)
    (lease,) = leases
    assert lease.released
    summary = run.summary
    (post,) = summary["post_segments"]

    if failed:
        state = _STATE[role]
        assert not (structures / "reactant.xyz").exists()
        assert not (structures / "product.xyz").exists()
        failure = json.loads((endpoint_dir / "failure.json").read_text())
        assert set(failure["failures"]) == {state}
        info = failure["failures"][state]
        if outcome == "hard_error":
            assert info["error_type"] == "RuntimeError"
            assert info["error"] == f"endpoint numerical failure ({role})"
        else:
            assert info["error_type"] == "ClickException"
            assert "Nonfinite final coordinates" in info["error"]
        assert info["geometry_role"] == "parent_irc_input_not_child_iterate"
        # The retained structure is the parent IRC input, not the child iterate.
        np.testing.assert_allclose(_xyz_values(info["irc_structure"]), coords(irc_frames[state]), atol=1e-6)
        assert (endpoint_dir / role / "iterate.xyz").is_file()
        # No refined consumer runs after the stage gate.
        assert "freq" not in run.child_names()
        assert not (segment / "energy_diagram_MLIP.png").exists()
        assert post["pipeline_stop"]["stage"] == "endpoint_opt"
        assert post["pipeline_stop"]["reason"] == "endpoint_execution_failed"
        assert summary["execution_status"] == "failed"
    else:
        assert (structures / "reactant.xyz").is_file()
        assert (structures / "product.xyz").is_file()
        assert endpoint_dir.exists() is (dump or outcome == "not_converged")
        assert run.child_names().count("freq") == 3
        assert "pipeline_stop" not in post
        converged = {state: True for state in ("reactant", "product")}
        converged[_STATE[role]] = outcome == "converged"
        assert post["endpoint_opt"]["reactant_converged"] is converged["reactant"]
        assert post["endpoint_opt"]["product_converged"] is converged["product"]
    for opt_role, state in _STATE.items():
        if not (failed and opt_role == role):
            assert post["endpoint_opt"][state]["n_opt_cycles"] == {"R": 3, "P": 4}[opt_role]


def test_ts_only_higher_energy_irc_endpoint_is_reactant(monkeypatch, tmp_path):
    loads, checkpoints = [], []
    real_load = hessian_cache.load
    real_checkpoint = workflow._save_role_scf_checkpoint
    real_irc = workflow._irc_and_match
    working = write(tmp_path / "ts.chk", "checkpoint fixture\n")

    def record_load(key):
        loads.append(key)
        return real_load(key)

    def irc_with_checkpoint(*args, **kwargs):
        result = real_irc(*args, **kwargs)
        result["scf_checkpoint"] = working
        return result

    def record_checkpoint(calculator, geometry_path, destination, *, manifest, key):
        checkpoints.append((key, Path(geometry_path).name))
        return real_checkpoint(calculator, geometry_path, destination, manifest=manifest, key=key)

    def hooks(_hooks):
        monkeypatch.setattr(hessian_cache, "load", record_load)
        monkeypatch.setattr(workflow, "_irc_and_match", irc_with_checkpoint)
        monkeypatch.setattr(workflow, "_save_role_scf_checkpoint", record_checkpoint)

    # IRC frame 0 (left, E=-1.00) is lower than frame 2 (right, E=-0.99).
    run = run_all(tmp_path, monkeypatch, frames=(1.0,), extra=["--tsopt"],
                  behaviour=Behaviour(opt={"P": "error"}), wrap=hooks)

    assert run.result.exit_code == 1, run.result.output[-3000:]
    structures = run.out / "segments" / "seg_01" / "structures"
    np.testing.assert_allclose(_xyz_values(structures / "reactant_irc.xyz"), coords(2.0), atol=1e-6)
    np.testing.assert_allclose(_xyz_values(structures / "product_irc.xyz"), coords(0.0), atol=1e-6)
    assert [key for key in loads if key.startswith("irc_")] == ["irc_right", "irc_left"]
    opt_inputs = [Path(args[args.index("-i") + 1]) for args in run.child_args("opt")]
    assert [Path(args[args.index("--out-dir") + 1]).name for args in run.child_args("opt")] == ["R", "P"]
    np.testing.assert_allclose(_xyz_values(opt_inputs[0]), coords(2.0), atol=1e-6)
    assert checkpoints == [
        ("post.01.scf_checkpoint.R", "reactant_irc.xyz"),
        ("post.01.scf_checkpoint.P", "product_irc.xyz"),
    ]
    failure = json.loads((run.out / "segments/seg_01/endpoint_opt/failure.json").read_text())
    assert set(failure["failures"]) == {"product"}
    (post,) = run.summary["post_segments"]
    assignment = post["endpoint_assignment"]
    assert assignment["policy"] == "higher_energy_endpoint_as_reactant"
    assert (assignment["left_role"], assignment["right_role"]) == ("product", "reactant")
    assert post["endpoint_opt"]["reactant_converged"] is True
    assert post["endpoint_opt"]["product_converged"] is None


@pytest.mark.parametrize(
    "outcome",
    ["converged", "not_converged", "stalled", "missing_result", "malformed_result",
     "failed", "missing_geometry", "stale", "nonfinite_coords", "nonfinite_energy", "energy_error"],
)
def test_endpoint_child_requires_current_finite_output(monkeypatch, tmp_path, outcome):
    import click
    from mlmm.workflows._run_session import ArtifactClaimError

    source = tmp_path / "input.pdb"
    source.write_text("END\n")
    output = tmp_path / "opt"
    output.mkdir()
    result_path = output / "result.json"
    final_path = output / "final_geometry.xyz"
    if outcome == "stale":
        result_path.write_text('{"optimization_status": "converged"}')
        final_path.write_text("prior invocation")
    prepared = SimpleNamespace(
        geom_path=source, source_path=source, cleanup=lambda: None,
    )
    monkeypatch.setattr(workflow, "prepare_input_structure", lambda *_a: prepared)

    def child(*_args, **_kwargs):
        if outcome == "stale":
            return
        if outcome != "missing_result":
            status = outcome if outcome in {"converged", "not_converged", "stalled", "failed"} else "converged"
            result_path.write_text("{" if outcome == "malformed_result" else json.dumps({
                "optimization_status": status, "stop_reason": "plateau" if outcome == "stalled" else None,
                "n_opt_cycles": 7, "max_cycles": 12,
            }))
        if outcome != "missing_geometry":
            final_path.write_text("terminal coordinates")

    terminal = Geometry(np.nan if outcome == "nonfinite_coords" else 2)
    if outcome == "nonfinite_energy":
        terminal.energy = np.inf
    if outcome == "energy_error":
        def fail_energy(_self):
            raise RuntimeError("endpoint energy evaluation failed")
        monkeypatch.setattr(Geometry, "energy", property(fail_energy), raising=False)
    monkeypatch.setattr(workflow, "_run_cli_main", child)
    monkeypatch.setattr(workflow, "geom_loader", lambda *_a, **_kw: terminal)
    monkeypatch.setattr(workflow, "_stage_calc_kwargs", lambda *_a, **_kw: {})
    closed = []
    monkeypatch.setattr(workflow, "_mlmm_calc", lambda **_kw: SimpleNamespace(close=lambda: closed.append(True)))

    endpoint_outcome = {}

    def run():
        return workflow._run_opt_for_state(
            source, 0, 1, tmp_path / "system.parm7", source, True,
            output, None, "grad", resolved_calc_template=object(),
            outcome=endpoint_outcome,
        )

    if outcome in {"converged", "not_converged", "stalled"}:
        result, path, converged = run()
        assert result is terminal and path == final_path
        assert converged is (outcome == "converged")
        assert endpoint_outcome["optimization_status"] == outcome
        assert endpoint_outcome["n_opt_cycles"] == 7
        assert endpoint_outcome["max_cycles"] == 12
        if outcome == "stalled":
            assert json.loads(result_path.read_text())["stop_reason"] == "plateau"
    elif outcome == "energy_error":
        with pytest.raises(RuntimeError, match="energy evaluation failed"):
            run()
        assert terminal.calculator is None and closed == [True]
    else:
        with pytest.raises((click.ClickException, ArtifactClaimError)):
            run()
