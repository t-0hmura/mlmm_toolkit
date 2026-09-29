"""Nested YAML spellings: opt.lbfgs = lbfgs, opt.rfo = rfo, freq.thermo = thermo."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import click
import numpy as np
import pytest
from click.testing import CliRunner

from mlmm.core.utils import apply_yaml_overrides

SMOKE = Path(__file__).resolve().parent / "smoke"


def test_both_spellings_are_combined_and_conflicts_rejected():
    lbfgs_cfg = {"max_step": 0.3, "memory": 5}
    yaml_cfg = {"lbfgs": {"memory": 7}, "opt": {"lbfgs": {"max_step": 0.2, "memory": 7}}}
    apply_yaml_overrides(yaml_cfg, [(lbfgs_cfg, (("lbfgs",), ("opt", "lbfgs")))])
    assert lbfgs_cfg == {"max_step": 0.2, "memory": 7}

    yaml_cfg["lbfgs"]["memory"] = 9
    with pytest.raises(click.BadParameter, match="lbfgs.memory and opt.lbfgs.memory conflict"):
        apply_yaml_overrides(yaml_cfg, [({}, (("lbfgs",), ("opt", "lbfgs")))])

    with pytest.raises(
        click.BadParameter, match="thermo.temperature and freq.thermo.temperature conflict"
    ):
        apply_yaml_overrides(
            {"thermo": {"temperature": 300.0}, "freq": {"thermo": {"temperature": 310.0}}},
            [({}, (("thermo",), ("freq", "thermo")))],
        )


def test_nested_sections_never_reach_the_parent_target():
    opt_cfg, freq_cfg, thermo_cfg = {"max_cycles": 1}, {"max_write": 1}, {"temperature": 298.15}
    yaml_cfg = {
        "opt": {"max_cycles": 5, "lbfgs": {"max_step": 0.1}, "rfo": {"trust_radius": 0.2}},
        "freq": {"max_write": 3, "thermo": {"temperature": 310.0}},
    }
    # No rfo target is registered, as in tsopt.
    apply_yaml_overrides(
        yaml_cfg,
        [
            (opt_cfg, (("opt",),)),
            (freq_cfg, (("freq",),)),
            (thermo_cfg, (("thermo",), ("freq", "thermo"))),
        ],
    )
    assert opt_cfg == {"max_cycles": 5}
    assert freq_cfg == {"max_write": 3}
    assert thermo_cfg == {"temperature": 310.0}


def test_dimer_line_search_section_is_not_a_nested_spelling():
    simple_cfg, lbfgs_cfg = {}, {"max_step": 0.3}
    apply_yaml_overrides(
        {"hessian_dimer": {"lbfgs": {"max_step": 0.1}}},
        [(simple_cfg, (("hessian_dimer",),)), (lbfgs_cfg, (("lbfgs",), ("opt", "lbfgs")))],
    )
    assert simple_cfg == {"lbfgs": {"max_step": 0.1}}
    assert lbfgs_cfg == {"max_step": 0.3}


def test_calc_and_mlmm_sections_are_still_one_section():
    calc_cfg = {}
    apply_yaml_overrides(
        {"calc": {"a": 1}, "mlmm": {"b": 2, "model_indices": [1, 2]}},
        [(calc_cfg, (("calc",), ("mlmm",)))],
    )
    assert calc_cfg == {"a": 1, "b": 2}
    with pytest.raises(click.BadParameter, match="Conflicting YAML values for calc.a and mlmm.a"):
        apply_yaml_overrides(
            {"calc": {"a": 1}, "mlmm": {"a": 2}}, [({}, (("calc",), ("mlmm",)))],
        )


def _scan_configs(yaml_cfg, kind):
    from mlmm.workflows.scan_common import resolve_scan_optimizer_configs

    return resolve_scan_optimizer_configs(
        yaml_cfg, kind=kind, thresh="baker", relax_max_cycles=10000,
        is_param_explicit=lambda name: False,
    )


@pytest.mark.parametrize("kind", ["lbfgs", "rfo"])
def test_scan_reads_nested_optimizer_section(kind):
    opt_cfg, sopt_cfg = _scan_configs({"opt": {kind: {"max_step": 0.05}}}, kind)
    assert sopt_cfg["max_step"] == 0.05
    assert "lbfgs" not in opt_cfg and "rfo" not in opt_cfg


def test_scan_rejects_conflicting_spellings():
    with pytest.raises(click.BadParameter, match="rfo.max_step and opt.rfo.max_step conflict"):
        _scan_configs({"rfo": {"max_step": 0.05}, "opt": {"rfo": {"max_step": 0.06}}}, "rfo")


def test_scan_nested_value_counts_as_explicit():
    with pytest.raises(click.BadParameter, match="opt.max_cycles and rfo.max_cycles conflict"):
        _scan_configs({"opt": {"max_cycles": 10, "rfo": {"max_cycles": 50}}}, "rfo")


def _opt_dry_run(tmp_path: Path, monkeypatch, config: dict, *extra: str):
    from mlmm.workflows import opt as opt_workflow

    path = tmp_path / "config.yaml"
    path.write_text(json.dumps(config))
    monkeypatch.chdir(tmp_path)
    return CliRunner().invoke(opt_workflow.cli, [
        "-i", str(SMOKE / "p_complex_layered.pdb"), "--parm", str(SMOKE / "p_complex.parm7"),
        "-q", "-1", "-m", "1", "--config", str(path), "--dry-run",
        "-o", str(tmp_path / "out"), *extra,
    ])


def test_opt_rejects_conflicting_spellings(tmp_path, monkeypatch):
    result = _opt_dry_run(
        tmp_path, monkeypatch,
        {"lbfgs": {"max_cycles": 50}, "opt": {"lbfgs": {"max_cycles": 60}}},
        "--opt-mode", "grad",
    )
    assert result.exit_code == 2, result.output
    assert "lbfgs.max_cycles and opt.lbfgs.max_cycles conflict" in result.output


def test_tsopt_constructor_gets_no_nested_optimizer_section(monkeypatch, tmp_path):
    import torch
    from pysisyphus.Geometry import Geometry

    from mlmm.workflows import tsopt

    class ConstructorReached(RuntimeError):
        pass

    captured = {}

    class CaptureOptimizer:
        def __init__(self, geometry, **kwargs):
            captured.update(kwargs)
            raise ConstructorReached("ConstructorReached: no model evaluation requested.")

    def load_geometry(_path, **kwargs):
        return Geometry(("O", "H", "H"), np.array([0., 0., 0., 1.8, 0., 0., -.4, 1.7, 0.]),
                        coord_type=kwargs.get("coord_type", "cart"))

    source, parm = tmp_path / "input.pdb", tmp_path / "input.parm7"
    source.write_text("Prepared-structure parsing is replaced for constructor capture.\n")
    parm.write_text("No MM engine or topology parser is constructed.\n")
    prepared = SimpleNamespace(original_path=source, source_path=source,
                               geom_path=source, cleanup=lambda: None)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(tsopt, "prepare_input_structure", lambda *_a, **_k: prepared)
    monkeypatch.setattr(tsopt, "resolve_charge_spin_or_raise", lambda *_a, **_k: (0, 1))
    monkeypatch.setattr(tsopt, "resolve_ml_layer_assignment", lambda **_k: (source, None))
    monkeypatch.setattr(tsopt, "geom_loader", load_geometry)
    monkeypatch.setattr(tsopt, "mlmm", lambda **_k:
                        SimpleNamespace(core=SimpleNamespace(hess_active_atoms=[0, 1, 2])))
    monkeypatch.setattr(tsopt, "_calc_full_hessian_torch", lambda *_a, **_k:
                        torch.eye(9, dtype=torch.float64))
    monkeypatch.setattr(tsopt, "_torch_device", lambda *_a: torch.device("cpu"))
    monkeypatch.setitem(tsopt.TSOPT_CLASS_MAP, "rsprfo", CaptureOptimizer)
    config = tmp_path / "config.yaml"
    config.write_text(json.dumps(
        {"opt": {"max_cycles": 7, "lbfgs": {"max_step": 0.1}, "rfo": {"trust_radius": 0.2}}}
    ))
    result = CliRunner().invoke(tsopt.cli, [
        "-i", str(source), "--parm", str(parm), "-q", "0", "-m", "1",
        "-o", str(tmp_path / "out"), "--no-flatten", "--no-dump", "--no-convert-files",
        "--no-microiter", "--config", str(config),
    ])
    assert result.exit_code == 1 and "ConstructorReached" in result.output, result.output
    assert captured["max_cycles"] == 7
    assert "lbfgs" not in captured and "rfo" not in captured
