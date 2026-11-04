"""Endpoint hard failures stop real all-branch control flow before consumers.

The endpoint-boundary tests run the whole ``mlmm all`` callback (TS-only and
MEP routes) with the child commands stubbed at the dispatcher
(``_all_stub_harness.py``); exception handlers, retention, the stage gate and
cleanup execute unchanged.
"""
from __future__ import annotations

import contextlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import click
import numpy as np
import pytest
import torch
from click.testing import CliRunner

from pysisyphus.constants import ANG2BOHR
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError
from _all_stub_harness import Behaviour, coords, pdb_text, read_values, run_all, write, xyz_text
from mlmm.cli import cli as root_cli
from mlmm.cli.completion import completion_guard, record_completion
from mlmm.core import calc_eval
from mlmm.workflows import all as workflow
from mlmm.workflows import freq
from mlmm.workflows import opt as opt_workflow
from mlmm.io import hessian_cache

REAL_RUN_CLI_MAIN = workflow._run_cli_main

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


def _record_irc_leases(monkeypatch, leases):
    real = workflow._irc_and_match

    def irc_and_match(*args, **kwargs):
        result = real(*args, **kwargs)
        leases.append(result["calculator_lease"])
        return result

    return lambda: monkeypatch.setattr(workflow, "_irc_and_match", irc_and_match)


@pytest.mark.parametrize("branch", ["ts", "seg"])
@pytest.mark.parametrize("dump", [False, True])
# A converged run does not depend on the endpoint; every other outcome is set on R and on P.
@pytest.mark.parametrize(("failed_endpoint", "outcome"), [
    (0, "hard_error"), (1, "hard_error"), (0, "nonfinite_error"), (1, "nonfinite_error"),
    (0, "converged"), (0, "not_converged"), (1, "not_converged"),
])
def test_endpoint_boundary_retains_provenance_and_stops_consumers(
    monkeypatch, tmp_path, branch, failed_endpoint, dump, outcome,
):
    role = _ROLE[failed_endpoint]
    leases = []
    run = run_all(
        tmp_path, monkeypatch,
        frames=(1.0,) if branch == "ts" else (0.0, 2.0),
        extra=["--tsopt", "--thermo", "--dft", *(["--dump"] if dump else [])],
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
        np.testing.assert_allclose(read_values(structures / f"{state}_irc.xyz"), coords(frame), atol=1e-6)
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
        np.testing.assert_allclose(read_values(info["irc_structure"]), coords(irc_frames[state]), atol=1e-6)
        assert (endpoint_dir / role / "iterate.xyz").is_file()
        # No refined consumer runs after the stage gate.
        assert "freq" not in run.child_names()
        assert "dft" not in run.child_names()
        assert not (segment / "energy_diagram_MLIP.png").exists()
        assert post["pipeline_stop"]["stage"] == "endpoint_opt"
        assert post["pipeline_stop"]["reason"] == "endpoint_execution_failed"
        assert summary["execution_status"] == "failed"
    else:
        assert (structures / "reactant.xyz").is_file()
        assert (structures / "product.xyz").is_file()
        assert endpoint_dir.exists() is (dump or outcome == "not_converged")
        if outcome == "converged":
            assert run.child_names().count("freq") == 3
            assert run.child_names().count("dft") == 3
            assert "pipeline_stop" not in post
        else:
            # An endpoint did not converge, so freq and DFT are skipped.
            assert "freq" not in run.child_names()
            assert "dft" not in run.child_names()
            assert post["pipeline_stop"] == {
                "stage": "endpoint_opt", "segment": 1, "reason": "endpoint_not_converged",
            }
        converged = {state: True for state in ("reactant", "product")}
        converged[_STATE[role]] = outcome == "converged"
        assert post["endpoint_opt"]["reactant_converged"] is converged["reactant"]
        assert post["endpoint_opt"]["product_converged"] is converged["product"]
    for opt_role, state in _STATE.items():
        if not (failed and opt_role == role):
            assert post["endpoint_opt"][state]["n_opt_cycles"] == {"R": 3, "P": 4}[opt_role]


@pytest.mark.parametrize("outcome", ["not_converged", "hard_error"])
def test_endpoint_opt_exit_code_goes_through_the_real_dispatch(monkeypatch, tmp_path, outcome):
    """Like the opt command, the endpoint child exits 1 when it does not converge
    (a completed run) and when it fails with an error; only the error fails the run,
    and an endpoint that does not converge skips freq."""

    def opt_with_its_exit_code():
        stub = workflow._run_cli_main

        def dispatch(name, cli, args, **kwargs):
            if name != "opt":
                return stub(name, cli, args, **kwargs)
            out = Path(args[args.index("--out-dir") + 1])

            @click.command(context_settings={"ignore_unknown_options": True, "allow_extra_args": True})
            def opt():
                try:
                    stub(name, cli, args, **kwargs)  # writes the result and geometry
                except RuntimeError:
                    sys.exit(1)
                status = json.loads((out / "result.json").read_text(encoding="utf-8"))["optimization_status"]
                record_completion({"status": status}, command="opt")

            opt.callback = completion_guard(opt.callback)
            return REAL_RUN_CLI_MAIN(name, opt, args, **kwargs)

        monkeypatch.setattr(workflow, "_run_cli_main", dispatch)

    run = run_all(
        tmp_path, monkeypatch, frames=(1.0,), extra=["--tsopt", "--thermo"],
        behaviour=Behaviour(opt={"P": _MODE[outcome]}), wrap=opt_with_its_exit_code,
    )
    summary = run.summary
    (post,) = summary["post_segments"]

    if outcome == "not_converged":
        assert run.result.exit_code == 0, run.result.output[-3000:]
        assert (summary["execution_status"], summary["scientific_status"]) == ("completed", "partial")
        assert post["endpoint_opt"]["product_converged"] is False
        assert post["pipeline_stop"]["reason"] == "endpoint_not_converged"
        assert "freq" not in run.child_names()
    else:
        assert run.result.exit_code == 1
        assert summary["execution_status"] == "failed"
        assert post["pipeline_stop"]["reason"] == "endpoint_execution_failed"
        assert "freq" not in run.child_names()


def _product_status_in_segment_1(monkeypatch, status):
    """The product endpoint opt of segment 1 ends with ``status``; the others converge."""

    def wrap():
        stub = workflow._run_cli_main

        def dispatch(name, cli, args, **kwargs):
            code = stub(name, cli, args, **kwargs)
            out = Path(args[args.index("--out-dir") + 1])
            # The endpoint child writes to segments/seg_NN/endpoint_opt/<R|P>.
            if name == "opt" and out.name == "P" and out.parent.parent.name == "seg_01":
                result = out / "result.json"
                payload = json.loads(result.read_text(encoding="utf-8"))
                payload.update(optimization_status=status, converged=False)
                result.write_text(json.dumps(payload), encoding="utf-8")
            return code

        monkeypatch.setattr(workflow, "_run_cli_main", dispatch)

    return wrap


@pytest.mark.parametrize("status", ["not_converged", "stalled"])
@pytest.mark.parametrize("branch", ["ts", "seg"])
def test_unconverged_endpoint_skips_freq_and_dft_of_its_segment(monkeypatch, tmp_path, branch, status):
    """freq and DFT are not run for a segment whose endpoint did not converge; the
    next segment is still processed and the result is partial with exit 0."""
    run = run_all(
        tmp_path, monkeypatch,
        frames=(1.0,) if branch == "ts" else (0.0, 2.0, 0.0),
        extra=["--tsopt", "--thermo", "--dft"],
        wrap=_product_status_in_segment_1(monkeypatch, status),
    )
    summary = run.summary
    first = summary["post_segments"][0]
    segment_dir = run.out / "segments" / "seg_01"

    def segments_of(child):
        return [
            Path(args[args.index("--out-dir") + 1]).parent.parent.name
            for args in run.child_args(child)
        ]

    assert run.result.exit_code == 0, run.result.output[-3000:]
    assert (summary["execution_status"], summary["scientific_status"]) == ("completed", "partial")
    assert summary["scientific_status_reasons"] == ["all:segment_1:endpoint_opt:product_not_converged"]
    for stage in ("thermo", "dft"):
        assert (
            f"[{stage}] seg_01: skipped because the product endpoint optimization did not converge"
            in run.result.output
        )
    assert first["pipeline_stop"] == {
        "stage": "endpoint_opt", "segment": 1, "reason": "endpoint_not_converged",
    }
    assert first["endpoint_opt"]["product"]["optimization_status"] == status
    assert isinstance(first["mlip"], dict)
    assert not {"gibbs_mlip", "dft", "gibbs_dft_mlip", "thermo_symmetry"} & set(first)
    assert (segment_dir / "energy_diagram_MLIP.png").is_file()
    assert not (segment_dir / "freq").exists() and not (segment_dir / "dft").exists()
    # A G or DFT diagram over all segments needs every segment.
    assert (run.out / "energy_diagram_MLIP_all.png").is_file()
    for stem in ("G_MLIP", "DFT", "G_DFT_plus_MLIP"):
        assert not (run.out / f"energy_diagram_{stem}_all.png").exists()
    assert summary["rate_limiting_step"]["method"] == "MLIP"
    if branch == "ts":
        assert segments_of("freq") == [] and segments_of("dft") == []
        assert "dft_status" not in summary["config"]
    else:
        assert segments_of("freq") == ["seg_02"] * 3 and segments_of("dft") == ["seg_02"] * 3
        second = summary["post_segments"][1]
        assert "pipeline_stop" not in second
        assert {"gibbs_mlip", "dft", "gibbs_dft_mlip"} <= set(second)


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

    def hooks():
        monkeypatch.setattr(hessian_cache, "load", record_load)
        monkeypatch.setattr(workflow, "_irc_and_match", irc_with_checkpoint)
        monkeypatch.setattr(workflow, "_save_role_scf_checkpoint", record_checkpoint)

    # IRC frame 0 (left, E=-1.00) is lower than frame 2 (right, E=-0.99).
    run = run_all(tmp_path, monkeypatch, frames=(1.0,), extra=["--tsopt"],
                  behaviour=Behaviour(opt={"P": "error"}), wrap=hooks)

    assert run.result.exit_code == 1, run.result.output[-3000:]
    structures = run.out / "segments" / "seg_01" / "structures"
    np.testing.assert_allclose(read_values(structures / "reactant_irc.xyz"), coords(2.0), atol=1e-6)
    np.testing.assert_allclose(read_values(structures / "product_irc.xyz"), coords(0.0), atol=1e-6)
    assert [key for key in loads if key.startswith("irc_")] == ["irc_right", "irc_left"]
    opt_inputs = [Path(args[args.index("-i") + 1]) for args in run.child_args("opt")]
    assert [Path(args[args.index("--out-dir") + 1]).name for args in run.child_args("opt")] == ["R", "P"]
    np.testing.assert_allclose(read_values(opt_inputs[0]), coords(2.0), atol=1e-6)
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


# An endpoint optimization that stops early: the real opt command runs on the
# two-atom fixture with an ML/MM stand-in (atom 0 is ML, atom 1 movable MM).
_K = 0.01  # Hartree/Bohr^2


class _Well:
    """E = 1/2 K |x - minimum|^2 (Hartree, Bohr); ``minimum`` is set per test."""

    minimum = np.zeros(6)

    def __init__(self, **kwargs):
        self.freeze_atoms = list(kwargs.get("freeze_atoms", []))
        self.core = SimpleNamespace(
            ml_indices=[0], hess_mm_indices=[1], movable_mm_indices=[1],
            frozen_layer_indices=[], hess_active_atoms=[0, 1], mlmm_links=[],
            _ml_backend=None,
        )

    def get_energy(self, atoms, coords):
        d = np.asarray(coords, dtype=float).reshape(-1) - self.minimum
        return {"energy": float(0.5 * _K * d @ d)}

    def get_forces(self, atoms, coords):
        d = np.asarray(coords, dtype=float).reshape(-1) - self.minimum
        return {**self.get_energy(atoms, coords), "forces": -_K * d}

    def get_hessian(self, atoms, coords):
        return {**self.get_forces(atoms, coords), "hessian": _K * torch.eye(6, dtype=torch.float64)}


class _ZeroStepLBFGS(LBFGS):
    """LBFGS whose second step has zero length, so pysisyphus raises ZeroStepLength."""

    def optimize(self):
        step = super().optimize()
        return 0.0 * step if self.cur_cycle >= 1 else step


class _RejectedRFO(RFOptimizer):
    """RFO whose second cycle fails, as on a rejected trust-radius update."""

    def optimize(self):
        if self.cur_cycle >= 1:
            raise OptimizationError("Actual energy change too high!")
        return super().optimize()


# path: (opt arguments, patched optimizer, its exception text, the standalone error line)
_EARLY_STOPS = {
    "lbfgs": (["--opt-mode", "grad"], "LBFGS", _ZeroStepLBFGS, "",
              "ERROR: Step length fell below the minimum allowed (ZeroStepLength)."),
    "rfo": (["--opt-mode", "hess", "--no-microiter"], "RFOptimizer", _RejectedRFO,
            "Actual energy change too high!",
            "ERROR: Optimization failed - Actual energy change too high!"),
    # The macro (ONIOM) step and the MM relaxation of the microiteration.
    "macro": (["--opt-mode", "hess", "--microiter"], "RFOptimizer", _RejectedRFO,
              "Actual energy change too high!",
              "ERROR: Optimization failed - Actual energy change too high!"),
    "micro": (["--opt-mode", "hess", "--microiter"], "LBFGS", _ZeroStepLBFGS, "",
              "ERROR: Step length fell below the minimum allowed (ZeroStepLength)."),
}


def _bohr(values):
    return np.ravel(values) * ANG2BOHR


def _install_endpoint_opt(monkeypatch, path, minimum):
    """Let the real opt command run on the fixture; its optimizer stops early on ``path``."""
    _args, name, optimizer, _text, _error = _EARLY_STOPS[path]
    monkeypatch.setattr(opt_workflow, name, optimizer)
    monkeypatch.setattr(_Well, "minimum", np.asarray(minimum, dtype=float))
    monkeypatch.setattr(hessian_cache, "_cache", {})
    for module in (opt_workflow, freq, calc_eval):
        monkeypatch.setattr(module, "mlmm", _Well)
    monkeypatch.setattr(opt_workflow, "mlmm_mm_only", lambda _core, **kw: _Well(**kw))
    monkeypatch.setattr(opt_workflow, "resolve_charge_spin_or_raise", lambda *_a, **_k: (0, 1))
    monkeypatch.setattr(opt_workflow, "load_pdb_atom_metadata", lambda *_a: [])

    def layers(*, source_path, calc_cfg, **_kwargs):
        calc_cfg["model_pdb"] = str(source_path)
        return source_path, None

    monkeypatch.setattr(opt_workflow, "resolve_ml_layer_assignment", layers)


def _well_minimum(frame, path):
    """The C atom is pulled away from the start; O only on the MM-relaxation path."""
    offset = np.zeros((2, 3))
    offset[0] = (-1.0, 0.0, -1.0)
    if path == "micro":
        offset[1] = (1.0, 0.0, 1.0)
    return _bohr(coords(frame)) + offset.ravel()


@pytest.mark.parametrize("endpoint", [True, False])
@pytest.mark.parametrize("path", sorted(_EARLY_STOPS))
def test_opt_that_stops_early_is_not_converged_only_for_an_endpoint(
    monkeypatch, tmp_path, path, endpoint,
):
    """For an endpoint of ``all``, an opt that stops early keeps its last geometry
    and is not converged, as in pdb2reaction ``all``; the standalone opt command
    reports the error, as pdb2reaction ``opt`` does."""
    args, _name, _optimizer, text, error_line = _EARLY_STOPS[path]
    _install_endpoint_opt(monkeypatch, path, _well_minimum(2.0, path))
    pdb = write(tmp_path / "product.pdb", pdb_text(2.0))
    xyz = write(tmp_path / "product.xyz", xyz_text(2.0))
    parm = write(tmp_path / "fixture.parm7", "Not read: the calculator is replaced.\n")
    config = write(tmp_path / "cpu.yaml", "calc:\n  ml_device: cpu\n")
    out = tmp_path / "opt"
    block = (
        opt_workflow.keep_last_geometry_on_early_stop("seg_01_product")
        if endpoint else contextlib.nullcontext()
    )
    with block:
        result = CliRunner().invoke(root_cli, [
            "opt", "-i", str(xyz), "--ref-pdb", str(pdb), "--parm", str(parm),
            "--model-pdb", str(pdb), "-q", "0", "-m", "1", "--config", str(config),
            "--out-dir", str(out), "--out-json", "--no-convert-files", *args,
        ])
    report = json.loads((out / "result.json").read_text(encoding="utf-8"))

    assert result.exit_code == 1, result.output[-3000:]
    if endpoint:
        assert (
            f"[endpoint-opt] WARNING: optimization for 'seg_01_product' terminated early "
            f"({text}); using last geometry."
        ) in result.output
        assert error_line not in result.output
        assert (
            report["execution_status"], report["scientific_status"], report["optimization_status"],
        ) == ("completed", "failed", "not_converged")
        assert "error_type" not in report
        assert report["n_opt_cycles"] == {"lbfgs": 2, "rfo": 2, "macro": 1, "micro": 0}[path]
        # The result is the last geometry, not the input.
        final = _bohr(read_values(out / "final_geometry.xyz"))
        assert not np.allclose(final, _bohr(coords(2.0)), atol=1e-6)
        expected = 0.5 * _K * np.sum((final - _Well.minimum) ** 2)
        assert report["energy_hartree"] == pytest.approx(expected, rel=1e-6)
    else:
        assert error_line in result.output
        assert (report["execution_status"], report["scientific_status"]) == ("failed", "failed")
        assert report["error_type"] == ("OptimizationError" if text else "ZeroStepLength")


@pytest.mark.parametrize("path", ["macro", "lbfgs"])
@pytest.mark.parametrize("branch", ["ts", "seg"])
def test_endpoint_that_stops_early_skips_freq_and_dft_of_its_segment(
    monkeypatch, tmp_path, branch, path,
):
    """An endpoint opt that stops early (OptimizationError on the default
    microiteration path, ZeroStepLength with LBFGS) is not converged, as in
    pdb2reaction: the segment keeps the last geometry, freq and DFT of the
    segment are skipped, and the run is partial with exit 0."""
    _args, _name, _optimizer, text, _error = _EARLY_STOPS[path]
    # The product IRC endpoint: frame 0 on the TS-only route, frame 2 on the MEP route.
    frame = 0.0 if branch == "ts" else 2.0

    def wrap():
        _install_endpoint_opt(monkeypatch, path, _well_minimum(frame, path))
        stub = workflow._run_cli_main

        def dispatch(name, cli, args, **kwargs):
            # The product endpoint runs the real opt command.
            if name == "opt" and Path(args[args.index("--out-dir") + 1]).name == "P":
                return REAL_RUN_CLI_MAIN(name, cli, args, **kwargs)
            return stub(name, cli, args, **kwargs)

        monkeypatch.setattr(workflow, "_run_cli_main", dispatch)

    run = run_all(
        tmp_path, monkeypatch,
        frames=(1.0,) if branch == "ts" else (0.0, 2.0),
        extra=["--tsopt", "--thermo", "--dft",
               *(["--opt-mode-post", "grad"] if path == "lbfgs" else [])],
        wrap=wrap,
    )
    summary = run.summary
    (post,) = summary["post_segments"]
    tag = "product" if branch == "ts" else "seg_01_product"
    segment = run.out / "segments" / "seg_01"
    child = segment / "endpoint_opt" / "P"

    assert run.result.exit_code == 0, run.result.output[-3000:]
    assert (summary["execution_status"], summary["scientific_status"]) == ("completed", "partial")
    assert summary["scientific_status_reasons"] == ["all:segment_1:endpoint_opt:product_not_converged"]
    assert (
        f"[endpoint-opt] WARNING: optimization for '{tag}' terminated early ({text}); "
        "using last geometry."
    ) in run.result.output
    for stage in ("thermo", "dft"):
        assert (
            f"[{stage}] seg_01: skipped because the product endpoint optimization did not converge"
            in run.result.output
        )
    assert post["pipeline_stop"] == {
        "stage": "endpoint_opt", "segment": 1, "reason": "endpoint_not_converged",
    }
    assert post["endpoint_opt"]["product_converged"] is False
    assert post["endpoint_opt"]["product"]["optimization_status"] == "not_converged"
    assert "failures" not in post["endpoint_opt"]
    assert "freq" not in run.child_names() and "dft" not in run.child_names()
    child_result = json.loads((child / "result.json").read_text(encoding="utf-8"))
    assert (child_result["execution_status"], child_result["optimization_status"]) == (
        "completed", "not_converged",
    )
    # The segment keeps the last geometry of the stopped optimization.
    product = read_values(segment / "structures" / "product.xyz")
    np.testing.assert_allclose(product, read_values(child / "final_geometry.xyz"), atol=1e-5)
    assert not np.allclose(product, coords(frame), atol=1e-5)
