"""Contracts for uphill rejection in minimum and transition-state searches.

mlmm's `all` runs post-IRC endpoint re-optimization by invoking the `opt` child
CLI, so the toggle is threaded into `_run_opt_for_state` and forwarded to that
child as `--reject-uphill` / `--no-reject-uphill`. These tests pin:

1. both `opt` and `all` expose the shipped default-off toggle;
2. an omitted parent toggle remains omitted for the child optimizer;
3. `_run_opt_for_state` faithfully forwards its resolved toggle, and the
   parent `all` workflow preserves omission as ``None``;
4. TS optimizers force uphill rejection off even when YAML-like input tries to
   re-enable it.
"""

from __future__ import annotations

from pathlib import Path

import click
import pytest

import mlmm.workflows.all as allmod
from mlmm.cli import cli as root_cli
from mlmm.core.defaults import RFO_KW
from mlmm.workflows.tsopt import (
    _build_rsirfo_kwargs,
    _force_ts_reject_uphill_off,
)


def test_shipped_default_is_reject_uphill_off() -> None:
    assert RFO_KW["reject_uphill"] is False
    assert RFO_KW["uphill_tolerance"] == 1e-4


def test_ts_rfo_forces_reject_uphill_off(tmp_path: Path) -> None:
    kwargs = _build_rsirfo_kwargs(
        {"reject_uphill": True},
        max_cycles=10,
        out_dir=tmp_path,
    )
    assert kwargs["reject_uphill"] is False


def test_ts_dimer_forces_reject_uphill_off() -> None:
    hostile_yaml_cfg = {"reject_uphill": True}
    effective = _force_ts_reject_uphill_off(hostile_yaml_cfg)
    assert effective["reject_uphill"] is False
    assert hostile_yaml_cfg["reject_uphill"] is True


@pytest.mark.parametrize("command", ["opt", "all"])
def test_toggle_is_exposed_with_default_off(command: str) -> None:
    ctx = click.Context(root_cli)
    cmd = root_cli.get_command(ctx, command)
    assert cmd is not None
    param = next(
        (p for p in cmd.params if isinstance(p, click.Option) and p.name == "reject_uphill"),
        None,
    )
    assert param is not None, f"{command} is missing --reject-uphill"
    assert param.opts == ["--reject-uphill"]
    assert param.secondary_opts == ["--no-reject-uphill"]
    assert param.is_bool_flag is True
    assert param.default is False


class _StopBeforeChild(Exception):
    """Raised by the stub opt runner to capture forwarded args without running."""


class _PreparedStub:
    def __init__(self, path: Path) -> None:
        self.geom_path = path
        self.source_path = path

    def cleanup(self) -> None:  # noqa: D401 - test stub
        pass


@pytest.mark.parametrize(
    "reject_uphill, expected_token",
    [(None, None), (True, "--reject-uphill"), (False, "--no-reject-uphill")],
    ids=["default-none", "explicit-on", "explicit-off"],
)
def test_endpoint_child_forwards_reject_uphill(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, reject_uphill, expected_token
) -> None:
    captured: dict[str, list[str]] = {}

    monkeypatch.setattr(
        allmod, "prepare_input_structure", lambda p: _PreparedStub(Path(p))
    )

    def _capture(_name, _cli, args, **_kwargs):
        captured["args"] = list(args)
        raise _StopBeforeChild()

    monkeypatch.setattr(allmod, "_run_cli_main", _capture)

    pdb = tmp_path / "endpoint.pdb"
    pdb.write_text("END\n", encoding="utf-8")

    with pytest.raises(_StopBeforeChild):
        allmod._run_opt_for_state(
            pdb,
            0,
            1,
            tmp_path / "real.parm7",
            tmp_path / "model.pdb",
            False,
            tmp_path / "out",
            None,
            "hess",
            resolved_calc_template=object(),
            reject_uphill=reject_uphill,
        )

    args = captured["args"]
    if expected_token is None:
        assert "--reject-uphill" not in args
        assert "--no-reject-uphill" not in args
    else:
        assert expected_token in args


@pytest.mark.parametrize(
    "extra, expected_eff",
    [(["--no-reject-uphill"], False), ([], None), (["--reject-uphill"], True)],
    ids=["off-arm", "default-unchanged", "on-arm"],
)
def test_all_gate_resolves_endpoint_reject_uphill(tmp_path: Path, extra, expected_eff) -> None:
    """Positive control for the benchmark's on/off arms.

    Parse `all` the way the CLI does and reproduce the exact resolution line
    ``_reject_uphill_eff = bool(reject_uphill) if _is_param_explicit(...) else None``.
    The default arm therefore forwards no explicit override.
    """
    from mlmm.workflows.all import cli as all_cli
    from mlmm.cli.decorators import make_is_param_explicit

    pdb = tmp_path / "x.pdb"
    pdb.write_text("END\n", encoding="utf-8")
    ctx = all_cli.make_context("all", ["-i", str(pdb), "--tsopt", "True", *extra])
    is_explicit = make_is_param_explicit(ctx)
    eff = bool(ctx.params["reject_uphill"]) if is_explicit("reject_uphill") else None
    assert eff is expected_eff
