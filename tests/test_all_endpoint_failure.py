"""Endpoint hard failures stop real all-branch control flow before consumers.

The contiguous production statements from IRC snapshot publication through the
optimized-structure handoff are compiled without CLI setup. Optimizers and I/O
boundaries are injected; exception handlers, retention, return/continue, and
cleanup execute unchanged.
"""
from __future__ import annotations

import ast
import builtins
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

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
            "endpoint_1_converged": None,
            "endpoint_2_converged": True,
            "failures": {"endpoint_1": {"error_type": "ValueError", "error": "boom"}},
        },
    }]
    workflow._enrich_summary(
        summary, version="", pipeline_mode="tsopt-only", out_dir=tmp_path,
        mlip_backend="orb", mlip_model="orb_v3_conservative_omol",
        charge=0, spin=1, post_segments=post,
        config={"tsopt": True, "thermo": True, "dft": True},
    )
    assert summary["scientific_status"] == "partial"
    assert summary["execution_status"] == "failed"
    assert any(
        "endpoint_1_execution_failed" in reason
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
            "endpoint_1_converged": True,
            "endpoint_2_converged": False,
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
    assert "all:segment_1:endpoint_opt:endpoint_2_not_converged" in truth.status_reasons


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


def endpoint_statements(branch):
    tree = ast.parse(Path(workflow.__file__).read_text())
    owner = "tsroot" if branch == "ts" else "seg_dir"
    save_name = "_save_single_geom_for_tools"
    for node in ast.walk(tree):
        body = getattr(node, "body", None)
        if not isinstance(body, list):
            continue
        endpoint_index = next((
            index for index, statement in enumerate(body)
            if isinstance(statement, ast.Assign)
            and any(isinstance(target, ast.Name) and target.id == "endpoint_opt_dir"
                    for target in statement.targets)
            and isinstance(statement.value, ast.BinOp)
            and isinstance(statement.value.left, ast.Name)
            and statement.value.left.id == owner
        ), None)
        if endpoint_index is None:
            continue
        saves = [
            index for index, statement in enumerate(body)
            if any(isinstance(call, ast.Call) and isinstance(call.func, ast.Name)
                   and call.func.id == save_name for call in ast.walk(statement))
        ]
        previous = [index for index in saves if index < endpoint_index]
        following = [index for index in saves if index > endpoint_index]
        start, end = previous[-3], following[1]
        guard = next(statement for statement in body[endpoint_index:end]
                     if isinstance(statement, ast.If)
                     and isinstance(statement.test, ast.Name)
                     and statement.test.id == "_endpoint_failures")
        # All real refined consumers follow this stage gate. P segment
        # frequency/DFT are in the outer loop, so their line numbers are checked
        # in that loop as well.
        scope = node
        if branch == "seg":
            scope = next(parent for parent in ast.walk(tree)
                         if isinstance(parent, ast.For)
                         and any(child is node for child in ast.walk(parent)))
        consumers = [
            call for call in ast.walk(scope)
            if isinstance(call, ast.Call) and isinstance(call.func, ast.Name)
            and call.func.id in {
                "_run_freq_for_state", "_run_dft_for_state", "_run_dft_sequence",
                "_write_public_energy_diagram", "_write_public_segment_diagram",
                "_copy_structures_to_seg_dir",
            }
            and call.lineno > body[endpoint_index].lineno
        ]
        assert consumers
        assert all(call.lineno > guard.end_lineno for call in consumers)
        return copy.deepcopy(body[start:end + 1])
    raise AssertionError(f"Endpoint branch not found: {branch}")


@pytest.mark.parametrize("branch", ["ts", "seg"])
@pytest.mark.parametrize("failed_endpoint", [0, 1])
@pytest.mark.parametrize("dump", [False, True])
@pytest.mark.parametrize("outcome", ["hard_error", "nonfinite_error", "converged", "not_converged"])
def test_endpoint_boundary_retains_provenance_and_stops_consumers(
    monkeypatch, tmp_path, branch, failed_endpoint, dump, outcome,
):
    statements = endpoint_statements(branch)
    monkeypatch.setattr(hessian_cache, "load", lambda *_a: None)
    monkeypatch.setattr(hessian_cache, "discard", lambda *_a: None)
    events, records = [], []
    root = tmp_path / "segment"
    structures = root / "structures"
    structures.mkdir(parents=True)
    geometries = [Geometry(1), Geometry(2)]
    lease = SimpleNamespace(release=lambda: events.append("release"))

    def save(geom, _ref, directory, name):
        events.append(name)
        path = directory / f"{name}.xyz"
        path.write_text(geom.as_xyz())
        return (path, path)

    opt_calls = []

    def optimize(*args, **kwargs):
        index = len(opt_calls)
        opt_calls.append(index)
        child = args[6]
        geom = Geometry(index + 1)
        child.mkdir(parents=True, exist_ok=True)
        geom.cart_coords[0] = 10 + index
        (child / "iterate.xyz").write_text(geom.as_xyz())
        if index == failed_endpoint and outcome in {"hard_error", "nonfinite_error"}:
            if outcome == "nonfinite_error":
                geom.cart_coords[0] = np.nan
            raise ValueError("endpoint numerical failure")
        return geom, child / "iterate.xyz", outcome != "not_converged"

    def write_record(path, payload, **_kwargs):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload))
        records.append(payload)
        return path

    def enrich(summary, **kwargs):
        summary["post_segments"] = kwargs["post_segments"]
        summary["scientific_status"] = "failed"
        records.append(summary)

    names = {node.id for statement in statements for node in ast.walk(statement)
             if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load)}
    namespace = {name: None for name in names if not hasattr(builtins, name)}
    namespace.update(vars(workflow))
    namespace.update({
        "gL": geometries[0], "gR": geometries[1], "gT": Geometry(3),
        "g_react_irc": geometries[0], "g_prod_irc": geometries[1],
        "g_react": geometries[0], "g_prod": geometries[1],
        "tsroot": root, "seg_dir": root, "struct_dir": structures,
        "out_dir": tmp_path, "dump": dump, "eL": 1.0, "eR_raw": 0.0,
        "reverse_irc": False, "irc_res": {"calculator_lease": lease},
        "seg_idx": 1, "seg_tag": "seg_01", "s": {}, "segment_log": {},
        "post_segment_logs": [], "endpoint_assignment": {},
        "_tsopt_record": {}, "_tsopt_decision": {},
        "_tsopt_payload": {}, "calculator_lease": lease,
        "resolved_calc_template": SimpleNamespace(materialize=lambda: {}),
        "_save_single_geom_for_tools": save, "_run_opt_for_state": optimize,
        "_validate_optimized_endpoint_pair": lambda *_a, **_kw: {},
        "commit_json": write_record, "commit_json_exact": write_record,
        "_enrich_summary": enrich,
        "_publish_manifest_summary": write_record,
        "_finalize_current_summary": write_record,
        "_copy_public_logged": lambda *_a, **_kw: None,
        "_persist_run_manifest": lambda *_a: None,
        "_write_endpoint_failure_summary_log": lambda *_a, **_kw: None,
        "_all_method_citation_payload": lambda: {},
        "_freeze_atoms_for_log": lambda: [],
        "_emit_final_summary": lambda *_a, **_kw: None,
    })
    namespace["__builtins__"] = __builtins__
    function = ast.parse("def exercise():\n    pass\n").body[0]
    if branch == "seg":
        wrapper = ast.parse("for iteration in [0]:\n    pass\n").body[0]
        wrapper.body = statements
        function.body = [wrapper]
    else:
        function.body = statements
    inputs = ("gL", "gR", "gT", "g_react", "g_prod")
    function.args.args = [ast.arg(arg=name) for name in inputs]
    exec(compile(ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[])),
                 str(workflow.__file__), "exec"), namespace)
    namespace["exercise"](*(namespace[name] for name in inputs))

    states = ["reactant", "product"] if branch == "seg" else ["endpoint_1", "endpoint_2"]
    assert opt_calls == [0, 1]
    assert "ts" in events
    for index, state in enumerate(states):
        assert f"H {index + 1}.0 0 0" in (structures / f"{state}_irc.xyz").read_text()
    failed = outcome in {"hard_error", "nonfinite_error"}
    if failed:
        assert not any(state in events for state in states)
        assert (root / "endpoint_opt" / "failure.json").is_file()
        failure = json.loads((root / "endpoint_opt" / "failure.json").read_text())
        info = failure["failures"][states[failed_endpoint]]
        assert info["error_type"] == "ValueError"
        assert info["error"] == "endpoint numerical failure"
        assert Path(info["irc_structure"]).is_file()
        assert "release" in events
        assert list((root / "endpoint_opt").rglob("iterate.xyz"))
        assert info["geometry_role"] == "parent_irc_input_not_child_iterate"
        assert geometries[failed_endpoint].cart_coords[0] == failed_endpoint + 1
        if branch == "ts":
            assert any(record.get("scientific_status") == "failed" for record in records)
        else:
            assert namespace["post_segment_logs"][0]["pipeline_stop"]["stage"] == "endpoint_opt"
    else:
        assert all(state in events for state in states)
        assert (root / "endpoint_opt").exists() is (dump or outcome == "not_converged")


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
        result_path.write_text('{"status": "converged"}')
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
                "status": status, "stop_reason": "plateau" if outcome == "stalled" else None,
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

    def run():
        return workflow._run_opt_for_state(
            source, 0, 1, tmp_path / "system.parm7", source, True,
            output, None, "grad", resolved_calc_template=object(),
        )

    if outcome in {"converged", "not_converged", "stalled"}:
        result, path, converged = run()
        assert result is terminal and path == final_path
        assert converged is (outcome == "converged")
        if outcome == "stalled":
            assert json.loads(result_path.read_text())["stop_reason"] == "plateau"
    elif outcome == "energy_error":
        with pytest.raises(RuntimeError, match="energy evaluation failed"):
            run()
        assert terminal.calculator is None and closed == [True]
    else:
        with pytest.raises((click.ClickException, ArtifactClaimError)):
            run()
