"""F2: `--workers` / `--workers-per-node` wiring.

Covers the calc_cfg router (``apply_workers_to_calc_cfg``): it sets the keys,
leaves them alone when no CLI value is given, and rejects an explicit
``Analytical`` Hessian request when workers > 1 (the parallel MLIP predictor
exposes no autograd model).
"""

from __future__ import annotations

import sys
import warnings

import pytest
import yaml

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm requires Python >= 3.11",
)


def _router():
    return pytest.importorskip("mlmm.backends").apply_workers_to_calc_cfg


def test_workers_keys_set_from_cli():
    cfg: dict = {}
    _router()(cfg, 4, 2)
    assert cfg["workers"] == 4
    assert cfg["workers_per_node"] == 2


def test_workers_none_is_noop():
    cfg = {"workers": 1, "workers_per_node": 1}
    _router()(cfg, None, None)
    assert cfg["workers"] == 1
    assert cfg["workers_per_node"] == 1


def test_workers_gt1_rejects_analytical_hessian():
    cfg = {"hessian_calc_mode": "Analytical"}
    with pytest.raises(ValueError, match="Analytical Hessian.*workers>1"):
        _router()(cfg, 2, 1)


def test_workers_eq1_keeps_analytical_hessian():
    cfg = {"hessian_calc_mode": "Analytical"}
    _router()(cfg, 1, 1)
    assert cfg["hessian_calc_mode"] == "Analytical"


def test_non_uma_worker_parallelism_is_reported_and_removed():
    cfg = {"backend": "orb"}
    with pytest.warns(UserWarning, match="does not use UMA worker parallelism"):
        _router()(cfg, 4, 2)
    assert "workers" not in cfg
    assert "workers_per_node" not in cfg


def test_non_uma_default_worker_values_are_silent():
    cfg = {"backend": "orb"}
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _router()(cfg, 1, 1)
    assert "workers" not in cfg
    assert "workers_per_node" not in cfg


def test_all_injects_workers_into_finite_difference_child_config(tmp_path):
    from mlmm.workflows.all import _inject_coord_type_into_args_yaml

    source = tmp_path / "all.yaml"
    source.write_text(
        "calc:\n  hessian_calc_mode: FiniteDifference\n", encoding="utf-8"
    )
    effective = _inject_coord_type_into_args_yaml(
        source, None, workers=2, workers_per_node=1
    )
    assert effective is not None and effective != source
    payload = yaml.safe_load(effective.read_text(encoding="utf-8"))
    assert payload["calc"]["workers"] == 2
    assert payload["calc"]["workers_per_node"] == 1


def test_all_rejects_parallel_analytical_child_config(tmp_path):
    from mlmm.workflows.all import _inject_coord_type_into_args_yaml

    source = tmp_path / "all.yaml"
    source.write_text("calc:\n  hessian_calc_mode: Analytical\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Analytical Hessian.*workers>1"):
        _inject_coord_type_into_args_yaml(
            source, None, workers=2, workers_per_node=1
        )


def test_direct_core_rejects_analytical_parallelism_before_io():
    """The Python API must fail before copying inputs or loading a model."""
    core_cls = pytest.importorskip("mlmm.backends.mlmm_calc").MLMMCore
    with pytest.raises(ValueError, match="Analytical Hessian.*workers>1"):
        core_cls(
            input_pdb="does-not-exist.pdb",
            real_parm7="does-not-exist.parm7",
            model_pdb="does-not-exist-model.pdb",
            workers=2,
            hessian_calc_mode="Analytical",
        )


def test_opt_registers_workers_options():
    """The --workers options must be attached to a real subcommand (opt)."""
    runner = pytest.importorskip("click.testing").CliRunner()
    root = pytest.importorskip("mlmm.cli").cli
    res = runner.invoke(root, ["opt", "--help-advanced"])
    assert "--workers" in res.output
    assert "--workers-per-node" in res.output
