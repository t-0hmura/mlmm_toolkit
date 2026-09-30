"""tsopt reports why --flatten did not run and which nested YAML sections it ignores."""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pytest
import torch
from click.testing import CliRunner
from pysisyphus.Geometry import Geometry

from mlmm.core.utils import unused_nested_yaml_sections
from mlmm.workflows import tsopt

TWO_IMAGINARY = (-100.0, -50.0) + (30.0,) * 7
BUDGET_BEFORE = "max-cycles budget exhausted before flattening"
BUDGET_DURING = "max-cycles budget exhausted during flattening"
FINAL_FREQ_SKIPPED = "final Hessian skipped (--skip-final-freq)"


def _load_geometry(_path, **kwargs):
    return Geometry(("O", "H", "H"), np.array([0., 0., 0., 1.8, 0., 0., -.4, 1.7, 0.]),
                    coord_type=kwargs.get("coord_type", "cart"))


def _patch_inputs(monkeypatch, tmp_path):
    source, parm = tmp_path / "input.pdb", tmp_path / "input.parm7"
    source.write_text("Prepared-structure parsing is replaced at its boundary.\n")
    parm.write_text("No MM engine or topology parser is constructed.\n")
    prepared = SimpleNamespace(original_path=source, source_path=source,
                               geom_path=source, cleanup=lambda: None)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(tsopt, "prepare_input_structure", lambda *_a, **_k: prepared)
    monkeypatch.setattr(tsopt, "resolve_charge_spin_or_raise", lambda *_a, **_k: (0, 1))
    monkeypatch.setattr(tsopt, "resolve_ml_layer_assignment", lambda **_k: (source, None))
    monkeypatch.setattr(tsopt, "geom_loader", _load_geometry)
    monkeypatch.setattr(tsopt, "mlmm", lambda **_k:
                        SimpleNamespace(core=SimpleNamespace(hess_active_atoms=[0, 1, 2])))
    monkeypatch.setattr(tsopt, "_calc_full_hessian_torch", lambda *_a, **_k:
                        torch.eye(9, dtype=torch.float64))
    monkeypatch.setattr(tsopt, "_torch_device", lambda *_a: torch.device("cpu"))
    return source, parm


def _run_microiter(monkeypatch, tmp_path, *, runs, max_cycles, extra_args):
    """``runs[i] = (cycles, converged)`` for the i-th microiteration TS run;
    ``cycles=None`` spends the whole budget given to that run."""
    source, parm = _patch_inputs(monkeypatch, tmp_path)
    calls = []
    phva_calls = []

    def fake_microiter(geometry, calc_cfg, rsirfo_cfg, lbfgs_cfg, opt_cfg, *args, **kwargs):
        cycles, converged = runs[len(calls)]
        calls.append(opt_cfg.get("max_cycles"))
        used = int(opt_cfg["max_cycles"]) if cycles is None else cycles
        optimizer = SimpleNamespace(
            is_converged=converged, is_stalled=False, cur_cycle=used - 1,
            stop_reason="", forces=[],
        )
        return {"optimizer": optimizer, "converged": converged, "cycles": used,
                "safeguards": {}, "micro_cycles": 0, "outcome": None}

    def exact_frequencies(optimizer, geometry, **_kwargs):
        phva_calls.append(1)
        return np.array(TWO_IMAGINARY), torch.eye(9, dtype=torch.float64), {}, None

    monkeypatch.setattr(tsopt, "_run_microiter_tsopt", fake_microiter)
    monkeypatch.setattr(tsopt, "_optimizer_exact_frequency_data", exact_frequencies)
    monkeypatch.setattr(tsopt, "_flatten_once_with_modes_for_geom", lambda *_a, **_k: True)
    monkeypatch.setattr(tsopt, "_write_all_imag_modes", lambda *_a, **_k: 0)
    monkeypatch.setattr(tsopt, "_calc_energy", lambda *_a, **_k: -1.0)

    out_dir = tmp_path / "out"
    result = CliRunner().invoke(tsopt.cli, [
        "-i", str(source), "--parm", str(parm), "-q", "0", "-m", "1",
        "-o", str(out_dir), "--opt-mode", "rsprfo", "--microiter",
        "--active-dof-mode", "all", "--max-cycles", str(max_cycles),
        "--no-dump", "--no-convert-files", "--out-json", *extra_args,
    ])
    result_json = out_dir / "result.json"
    assert result_json.is_file(), result.output
    return json.loads(result_json.read_text(encoding="utf-8")), result, phva_calls


def test_unconverged_max_cycles_stop_records_budget_reason(monkeypatch, tmp_path):
    payload, result, phva_calls = _run_microiter(
        monkeypatch, tmp_path, runs=[(None, False)], max_cycles=3,
        extra_args=["--flatten"],
    )

    assert payload["flatten_requested"] is True, result.output
    assert payload["flatten_skip_reason"] == BUDGET_BEFORE, result.output
    assert "Reached --max-cycles budget; skipping flatten loop." in result.output
    # S2: an unconverged max-cycles stop computes no Hessian.
    assert payload["hessian_status"] == "skipped"
    assert phva_calls == []


def test_skip_final_freq_records_skipped_hessian_reason(monkeypatch, tmp_path):
    payload, result, phva_calls = _run_microiter(
        monkeypatch, tmp_path, runs=[(1, True)], max_cycles=10,
        extra_args=["--flatten", "--skip-final-freq"],
    )

    assert payload["optimization_status"] == "converged", result.output
    assert payload["flatten_skip_reason"] == FINAL_FREQ_SKIPPED
    assert phva_calls == []


def test_unconverged_flatten_retry_records_budget_reason(monkeypatch, tmp_path):
    payload, result, phva_calls = _run_microiter(
        monkeypatch, tmp_path, runs=[(1, True), (None, False)], max_cycles=3,
        extra_args=["--flatten"],
    )

    assert payload["flatten_skip_reason"] == BUDGET_DURING, result.output
    assert len(phva_calls) == 1


@pytest.mark.parametrize(
    "runs, extra_args",
    [([(None, False)], []), ([(1, True)], ["--skip-final-freq"])],
    ids=["max_cycles", "skip_final_freq"],
)
def test_no_reason_without_flatten(monkeypatch, tmp_path, runs, extra_args):
    payload, result, _ = _run_microiter(
        monkeypatch, tmp_path, runs=runs, max_cycles=3,
        extra_args=["--no-flatten", *extra_args],
    )

    assert payload["flatten_requested"] is False, result.output
    assert payload["flatten_skip_reason"] is None


# ------------------------------------------------------------ YAML NOTE

def test_unused_nested_sections_are_listed():
    yaml_cfg = {
        "opt": {"max_cycles": 5, "lbfgs": {"max_step": 0.1}, "rfo": {"trust_radius": 0.2}},
        "freq": {"thermo": {"temperature": 310.0}},
    }
    assert unused_nested_yaml_sections(yaml_cfg) == ["opt.lbfgs", "opt.rfo", "freq.thermo"]
    # tsopt reads opt.lbfgs for the microiteration MM relaxation.
    assert unused_nested_yaml_sections(yaml_cfg, read=(("opt", "lbfgs"),)) == [
        "opt.rfo", "freq.thermo",
    ]


def test_read_or_empty_nested_sections_are_not_listed():
    yaml_cfg = {"opt": {"lbfgs": {"max_step": 0.1}, "rfo": None}, "freq": {"thermo": {}}}
    assert unused_nested_yaml_sections(yaml_cfg, read=(("opt", "lbfgs"),)) == []
    assert unused_nested_yaml_sections({"geom": {"coord_type": "cart"}}) == []


def _tsopt_until_constructor(monkeypatch, tmp_path, config):
    class ConstructorReached(RuntimeError):
        pass

    class CaptureOptimizer:
        def __init__(self, geometry, **kwargs):
            raise ConstructorReached("ConstructorReached: no model evaluation requested.")

    source, parm = _patch_inputs(monkeypatch, tmp_path)
    monkeypatch.setitem(tsopt.TSOPT_CLASS_MAP, "rsprfo", CaptureOptimizer)
    path = tmp_path / "config.yaml"
    path.write_text(json.dumps(config))
    result = CliRunner().invoke(tsopt.cli, [
        "-i", str(source), "--parm", str(parm), "-q", "0", "-m", "1",
        "-o", str(tmp_path / "out"), "--no-flatten", "--no-dump", "--no-convert-files",
        "--no-microiter", "--config", str(path),
    ])
    assert result.exit_code == 1 and "ConstructorReached" in result.output, result.output
    return result


def test_tsopt_notes_ignored_nested_sections_once(monkeypatch, tmp_path):
    result = _tsopt_until_constructor(monkeypatch, tmp_path, {
        "opt": {"max_cycles": 7, "lbfgs": {"max_step": 0.1}, "rfo": {"trust_radius": 0.2}},
        "freq": {"thermo": {"temperature": 310.0}},
    })

    notes = [line for line in result.output.splitlines() if "NOTE: Ignoring YAML" in line]
    assert notes == [
        "[tsopt] NOTE: Ignoring YAML sections that tsopt does not use: "
        "opt.rfo, freq.thermo."
    ]


def test_tsopt_prints_no_note_without_unused_sections(monkeypatch, tmp_path):
    result = _tsopt_until_constructor(monkeypatch, tmp_path, {
        "opt": {"max_cycles": 7, "lbfgs": {"max_step": 0.1}},
    })

    assert "NOTE: Ignoring YAML" not in result.output
