"""F2: `--workers` / `--workers-per-node` wiring.

Covers the calc_cfg router (``apply_workers_to_calc_cfg``): it sets the keys,
leaves them alone when no CLI value is given, and downgrades an explicit
``Analytical`` Hessian request to ``FiniteDifference`` when workers > 1 (the
parallel MLIP predictor exposes no autograd model).
"""

from __future__ import annotations

import sys

import pytest

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


def test_workers_gt1_downgrades_analytical_hessian():
    cfg = {"hessian_calc_mode": "Analytical"}
    _router()(cfg, 2, 1)
    assert cfg["workers"] == 2
    assert cfg["hessian_calc_mode"] == "FiniteDifference"


def test_workers_eq1_keeps_analytical_hessian():
    cfg = {"hessian_calc_mode": "Analytical"}
    _router()(cfg, 1, 1)
    assert cfg["hessian_calc_mode"] == "Analytical"


def test_opt_registers_workers_options():
    """The --workers options must be attached to a real subcommand (opt)."""
    runner = pytest.importorskip("click.testing").CliRunner()
    root = pytest.importorskip("mlmm.cli").cli
    res = runner.invoke(root, ["opt", "--help-advanced"])
    assert "--workers" in res.output
    assert "--workers-per-node" in res.output
