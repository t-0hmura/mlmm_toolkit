"""Exercise real OPT/TSOPT config merging and dispatch without model evaluation.

Input/calculator boundaries follow tests/test_tsopt_terminal_workflow_phva.py.
Dispatch cases stop at the microdriver or ordinary optimizer.run, not at dry-run.
Current commands support --config plus explicit CLI overrides, not --override;
these tests deliberately do not inject the dormant override-YAML resolver.
Two completion cases use an explicitly unconverged optimizer stub and the real
JSON writer. Their initial Hessian and final energy are mocked; no physical
optimization, terminal curvature, or scientific convergence is validated.
"""

from copy import deepcopy
import importlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from click.testing import CliRunner


_SENTINEL_EXIT = 86
_MOCK_FINAL_ENERGY = -1.25  # Serialization input only, not a calculated energy.
_MOCK_STOP_REASON = "mock optimizer did not converge"
_FROZEN = (0, 1, 2)
_COORDS = np.array(
    [[0., 0., 0.], [2., 0., 0.], [0., 2., 0.],
     [0., 0., 2.], [2., 2., 2.]]
).ravel()


class _DispatchReached(SystemExit):
    def __init__(self):
        super().__init__(_SENTINEL_EXIT)


@pytest.fixture
def dispatch_cli(tmp_path, monkeypatch):
    from pysisyphus.Geometry import Geometry

    monkeypatch.chdir(tmp_path)
    source, parm = tmp_path / "input.pdb", tmp_path / "input.parm7"
    source.write_text("REMARK Prepared-structure boundary is mocked.\nEND\n")
    parm.write_text("Topology and model constructors are mocked.\n")

    def invoke(command, mode, config_embed, cli_embed, microiter=True, *, complete=False):
        module = importlib.import_module(f"mlmm.workflows.{command}")
        events, calculations, seeds = [], [], []
        energy_calls = []
        prepared = SimpleNamespace(
            original_path=source, source_path=source, geom_path=source,
            cleanup=lambda: None,
        )

        class Calculator:
            def __init__(self, **kwargs):
                self.config = deepcopy(kwargs)
                self.freeze_atoms = list(kwargs.get("freeze_atoms", []))
                self.core = SimpleNamespace(
                    ml_indices=[3], hess_mm_indices=[4], movable_mm_indices=[],
                    frozen_layer_indices=list(_FROZEN), mlmm_links=[],
                    hess_active_atoms=[3, 4],
                    embedcharge=bool(kwargs.get("embedcharge", False)),
                )
                calculations.append(self)

        def stop(branch, geom, cfg, selected, optimizer_kwargs=None, *, terminate=True):
            events.append({
                "branch": branch, "mode": selected, "config": deepcopy(cfg),
                "freeze": tuple(map(int, geom.freeze_atoms)),
                "coords": geom.cart_coords.copy(), "calculator": geom.calculator,
                "optimizer_kwargs": deepcopy(optimizer_kwargs),
            })
            if terminate:
                raise _DispatchReached()

        def optimizer_class(selected):
            class Optimizer:
                def __init__(self, geom, **kwargs):
                    self.geom, self.kwargs = geom, kwargs

                def run(self):
                    stop("ordinary", self.geom, self.geom.calculator.config,
                         selected, self.kwargs, terminate=not complete)
                    # Completion is a serialization probe, not a fake success.
                    self.is_converged = False
                    self.is_stalled = False
                    self.cur_cycle = -1  # No numerical optimization cycle ran.
                    self.stop_reason = _MOCK_STOP_REASON
                    self.final_fn = self.get_path_for_fn("final_geometry.xyz")
                    self.final_fn.write_text(self.geom.as_xyz(), encoding="utf-8")

                def get_path_for_fn(self, name):
                    return Path(self.kwargs["out_dir"]) / name
            return Optimizer

        def microdriver(geom, *args, **kwargs):
            # OPT also passes base_calc before calc_cfg; TSOPT does not.
            cfg = args[1] if command == "opt" else args[0]
            if command == "opt":
                partition = kwargs["partition"]
                assert partition.macro_active_atoms == (3,)
                assert partition.micro_active_atoms == (4,)
            stop("micro", geom, cfg, "rfo" if command == "opt" else kwargs["mode"])

        def resolve_layer(**kwargs):
            kwargs["calc_cfg"]["model_pdb"] = str(source)
            return source, None

        def seed(geom, cfg, *args, **kwargs):
            seeds.append((tuple(map(int, geom.freeze_atoms)), deepcopy(cfg)))
            # TSOPT stores this constant matrix before reaching optimizer.run.
            return np.eye(geom.cart_coords.size)

        def final_energy(geom, calc_or_cfg):
            cfg = calc_or_cfg.config if command == "opt" else calc_or_cfg
            assert cfg["embedcharge"] is True
            assert cfg["freeze_atoms"] == list(_FROZEN)
            assert tuple(map(int, geom.freeze_atoms)) == _FROZEN
            np.testing.assert_array_equal(geom.cart_coords, _COORDS)
            energy_calls.append(_MOCK_FINAL_ENERGY)
            return _MOCK_FINAL_ENERGY

        monkeypatch.setattr(module, "prepare_input_structure", lambda *_a: prepared)
        monkeypatch.setattr(module, "resolve_charge_spin_or_raise", lambda *_a, **_k: (0, 1))
        monkeypatch.setattr(module, "resolve_ml_layer_assignment", resolve_layer)
        monkeypatch.setattr(
            module, "geom_loader",
            lambda _path, **kwargs: Geometry(["H"] * 5, _COORDS.copy(), **kwargs),
        )
        monkeypatch.setattr(module, "mlmm", Calculator)
        monkeypatch.setattr(module, f"_run_microiter_{command}", microdriver)
        if command == "opt":
            monkeypatch.setattr(module, "RFOptimizer", optimizer_class("rfo"))
            monkeypatch.setattr(module, "_seed_rfo_initial_hessian", seed)
        else:
            monkeypatch.setattr(module, "_calc_full_hessian_torch", seed)
            for selected in ("rsirfo", "rsprfo", "trim"):
                monkeypatch.setitem(module.TSOPT_CLASS_MAP, selected, optimizer_class(selected))
        if complete:
            energy_name = "unbiased_energy_hartree" if command == "opt" else "_calc_energy"
            monkeypatch.setattr(module, energy_name, final_energy)
            if command == "tsopt":
                # Final active-atom resolution constructs a temporary calculator,
                # even though this unconverged run must skip terminal PHVA.
                freq = importlib.import_module("mlmm.workflows.freq")
                monkeypatch.setattr(freq, "mlmm", Calculator)

        config = {
            "geom": {"coord_type": "cart", "freeze_atoms": [1, 2]},
            "calc": {"ml_device": "cpu", "embedcharge_cutoff": 7.5},
        }
        if command == "tsopt":
            config["hessian_dimer"] = {"device": "cpu"}
        if config_embed is not None:
            config["calc"]["embedcharge"] = config_embed
        config_path = tmp_path / "config.yaml"
        config_path.write_text(json.dumps(config))  # JSON is valid YAML; keeps bool types.
        args = [
            "-i", str(source), "--parm", str(parm), "-q", "0", "-m", "1",
            "--config", str(config_path), "--out-dir", str(tmp_path / "output"),
            "--opt-mode", mode, "--freeze-atoms", "3", "--max-cycles", "7",
            "--microiter" if microiter else "--no-microiter",
            "--no-flatten", "--no-dump", "--no-convert-files",
        ]
        if cli_embed is not None:
            args.append("--embedcharge" if cli_embed else "--no-embedcharge")
        if complete:
            args.append("--out-json")
        result = CliRunner().invoke(module.cli, args)
        assert result.exit_code == (0 if complete else _SENTINEL_EXIT), result.output
        assert len(events) == 1, (events, result.output)
        if complete:
            assert energy_calls == [_MOCK_FINAL_ENERGY]
            output = tmp_path / "output"
            report = json.loads((output / "result.json").read_text())
            assert json.loads((output / "summary.json").read_text()) == report
            events[0]["result_json"] = report
        return events[0], calculations, seeds

    return invoke


def _assert_retained(event, calculations, seeds, expected_branch, expected_embed, expected_mode):
    assert event["branch"] == expected_branch
    assert event["mode"] == expected_mode
    assert event["freeze"] == _FROZEN
    np.testing.assert_array_equal(event["coords"], _COORDS)
    assert event["config"]["embedcharge"] is expected_embed
    assert event["config"]["freeze_atoms"] == list(_FROZEN)
    assert event["config"]["embedcharge_cutoff"] == 7.5
    for calculator in calculations:
        assert calculator.config["embedcharge"] is expected_embed
        assert tuple(map(int, calculator.freeze_atoms)) == _FROZEN
    if expected_branch == "ordinary":
        assert event["calculator"] in calculations
        assert event["calculator"].core.embedcharge is expected_embed
        assert event["optimizer_kwargs"]["max_cycles"] == 7
        assert len(seeds) == 1
        assert seeds[0][0] == _FROZEN
        assert seeds[0][1]["embedcharge"] is expected_embed
    else:
        assert seeds == []  # Stop before any macro or micro optimization begins.


@pytest.mark.parametrize("command,expected_mode", [("opt", "rfo"), ("tsopt", "rsprfo")])
@pytest.mark.parametrize(
    "config_embed,cli_embed,microiter,expected_embed,expected_branch",
    [
        pytest.param(None, True, True, True, "ordinary", id="cli-enabled"),
        pytest.param(True, None, True, True, "ordinary", id="config-enabled"),
        pytest.param(False, True, True, True, "ordinary", id="cli-overrides-config-to-true"),
        pytest.param(True, False, True, False, "micro", id="cli-overrides-config-to-false"),
        pytest.param(False, None, True, False, "micro", id="config-disabled"),
        pytest.param(None, None, True, False, "micro", id="default-unembedded"),
        pytest.param(True, None, False, True, "ordinary", id="explicit-no-microiter"),
    ],
)
def test_embed_micro_dispatch_uses_final_config(
    dispatch_cli, command, expected_mode, config_embed, cli_embed,
    microiter, expected_embed, expected_branch,
):
    observed = dispatch_cli(command, "hess", config_embed, cli_embed, microiter)
    _assert_retained(*observed, expected_branch, expected_embed, expected_mode)


@pytest.mark.parametrize("mode", ["rsirfo", "rsprfo", "trim"])
def test_embedded_fallback_preserves_explicit_ts_optimizer(dispatch_cli, mode):
    observed = dispatch_cli("tsopt", mode, True, None)
    _assert_retained(*observed, "ordinary", True, mode)


@pytest.mark.parametrize("command,expected_mode", [("opt", "rfo"), ("tsopt", "rsprfo")])
def test_embedded_fallback_serializes_nonconverged_completion(
    dispatch_cli, command, expected_mode,
):
    observed = dispatch_cli(command, "hess", True, None, complete=True)
    _assert_retained(*observed, "ordinary", True, expected_mode)
    assert observed[0]["optimizer_kwargs"]["flatten_enabled"] is False
    report = observed[0]["result_json"]
    assert report["command"] == command
    assert report["status"] == "not_converged"
    assert report["stop_reason"] == _MOCK_STOP_REASON
    assert report["energy_hartree"] == _MOCK_FINAL_ENERGY
    assert report["n_opt_cycles"] == 0
    assert report["n_freeze_atoms"] == len(_FROZEN)
    assert report["microiteration"] == {
        "requested": True, "used": False, "fallback_reason": "embedcharge",
    }
    assert "n_micro_cycles" not in report
    if command == "tsopt":
        assert report["optimizer"] == expected_mode
        assert report["optimization_status"] == "not_converged"
        assert report["flatten_enabled"] is False
        assert report["hessian_status"] == "skipped"
        assert report["saddle_validation"] == "unavailable"
        assert report["saddle_order_verified"] is False
        assert report["n_imaginary_modes"] is None
        assert report["imaginary_frequencies_cm"] is None
