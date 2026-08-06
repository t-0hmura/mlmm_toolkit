"""Contracts for the opt-in energy-plateau stop (``--stop-plateau``).

A plateau stop reports ``stalled``; it never reports convergence and it never
replaces ``--max-cycles`` as the real bound.  Because a run that stops there
leaves the geometry unconverged, it is opt-in.  These tests pin:

1. the shipped default is off in every optimizer config block AND in the
   bundled ``Optimizer`` signature, so a construction site that passes explicit
   kwargs instead of a config block cannot silently enable it;
2. ``opt`` / ``tsopt`` / ``all`` all expose the toggle and its two value
   options, with the toggle defaulting to off;
3. every flag ``all`` forwards to a child exists on that child;
4. the ``all`` gate resolves to ``None`` unless the toggle is explicit, and an
   explicit toggle reaches the endpoint ``opt`` child and the ``tsopt``
   overrides.

The MM micro relaxation's exemption is pinned end-to-end on both microiteration
drivers in ``test_microiteration_cycle_carry.py``.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import click
import pytest

import mlmm.core.defaults as defaults_mod
import mlmm.workflows.all as allmod
from mlmm.cli import cli as root_cli
from mlmm.workflows._all_helpers import build_tsopt_overrides

PLATEAU_KEYS = (
    "energy_plateau",
    "energy_plateau_thresh",
    "energy_plateau_window",
)


def test_shipped_default_is_plateau_off() -> None:
    """Sweep every config block, not just ``OPT_BASE_KW``: a derived block that
    re-declared the key would reintroduce the default-on behavior for exactly
    one optimizer."""
    seen = 0
    for name, value in vars(defaults_mod).items():
        if not isinstance(value, dict):
            continue
        for cfg_name, cfg in [(name, value)] + [
            (f"{name}[{k!r}]", v) for k, v in value.items() if isinstance(v, dict)
        ]:
            if "energy_plateau" not in cfg:
                continue
            seen += 1
            assert cfg["energy_plateau"] is False, cfg_name
    assert seen >= 2, "no optimizer config block declares energy_plateau"
    assert defaults_mod.OPT_BASE_KW["energy_plateau"] is False
    assert defaults_mod.OPT_BASE_KW["energy_plateau_thresh"] == 1e-4
    assert defaults_mod.OPT_BASE_KW["energy_plateau_window"] == 50


def test_bundled_optimizer_default_is_plateau_off() -> None:
    """``align_freeze`` builds LBFGS from explicit kwargs with no config block,
    so the library default is the effective default there."""
    from pysisyphus.optimizers.Optimizer import Optimizer

    default = inspect.signature(Optimizer.__init__).parameters["energy_plateau"].default
    assert default is False


@pytest.mark.parametrize("command", ["opt", "tsopt", "all"])
def test_toggle_is_exposed_with_default_off(command: str) -> None:
    ctx = click.Context(root_cli)
    cmd = root_cli.get_command(ctx, command)
    params = {p.name: p for p in cmd.params if isinstance(p, click.Option)}

    toggle = params.get("stop_plateau")
    assert toggle is not None, f"{command} is missing --stop-plateau"
    assert toggle.opts == ["--stop-plateau"]
    assert toggle.secondary_opts == ["--no-stop-plateau"]
    assert toggle.is_bool_flag is True
    assert toggle.default is False

    for name, type_name in (
        ("stop_plateau_thresh", "float"),
        ("stop_plateau_window", "integer"),
    ):
        option = params.get(name)
        assert option is not None, f"{command} is missing --{name.replace('_', '-')}"
        assert option.default is None
        assert option.type.name == type_name


@pytest.mark.parametrize("child", ["opt", "tsopt"])
def test_forwarded_flags_exist_on_the_child(child: str) -> None:
    """``all`` forwards these by literal string; a typo would be silent."""
    ctx = click.Context(root_cli)
    cmd = root_cli.get_command(ctx, child)
    flags = {
        opt
        for param in cmd.params
        for opt in list(param.opts) + list(param.secondary_opts)
    }
    assert {
        "--stop-plateau",
        "--no-stop-plateau",
        "--stop-plateau-thresh",
        "--stop-plateau-window",
    } <= flags


@pytest.mark.parametrize(
    "extra, expected_eff",
    [([], None), (["--stop-plateau"], True), (["--no-stop-plateau"], False)],
    ids=["default-unchanged", "on-arm", "off-arm"],
)
def test_all_gate_resolves_stop_plateau(tmp_path: Path, extra, expected_eff) -> None:
    """Reproduce the resolution line ``_stop_plateau_eff = bool(stop_plateau)
    if explicit else None``: the default arm forwards no override at all, so
    the children keep their own default-off configuration."""
    from mlmm.workflows.all import cli as all_cli
    from mlmm.cli.decorators import make_is_param_explicit

    pdb = tmp_path / "x.pdb"
    pdb.write_text("END\n", encoding="utf-8")
    ctx = all_cli.make_context("all", ["-i", str(pdb), "--tsopt", "True", *extra])
    is_explicit = make_is_param_explicit(ctx)
    eff = bool(ctx.params["stop_plateau"]) if is_explicit("stop_plateau") else None
    assert eff is expected_eff


@pytest.mark.parametrize(
    "stop_plateau, expected_token",
    [(None, None), (True, "--stop-plateau"), (False, "--no-stop-plateau")],
    ids=["default-none", "explicit-on", "explicit-off"],
)
def test_endpoint_child_forwards_stop_plateau(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, stop_plateau, expected_token
) -> None:
    captured: dict[str, list[str]] = {}

    class _StopBeforeChild(Exception):
        pass

    class _PreparedStub:
        def __init__(self, path: Path) -> None:
            self.geom_path = path
            self.source_path = path

        def cleanup(self) -> None:  # noqa: D401 - test stub
            pass

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
            stop_plateau=stop_plateau,
            stop_plateau_window=None if stop_plateau is None else 7,
        )

    args = captured["args"]
    if expected_token is None:
        assert "--stop-plateau" not in args
        assert "--no-stop-plateau" not in args
        assert "--stop-plateau-window" not in args
    else:
        assert expected_token in args
        assert args[args.index("--stop-plateau-window") + 1] == "7"


@pytest.mark.parametrize(
    "stop_plateau, expected",
    [(None, {}), (True, {"stop_plateau": True}), (False, {"stop_plateau": False})],
    ids=["default-none", "explicit-on", "explicit-off"],
)
def test_tsopt_overrides_carry_stop_plateau(stop_plateau, expected) -> None:
    overrides = build_tsopt_overrides(
        tsopt_max_cycles=None,
        dump=False,
        dump_override_requested=False,
        tsopt_out_dir=None,
        hessian_calc_mode=None,
        opt_mode_post_norm=None,
        opt_mode_post_set=False,
        opt_mode_set=False,
        tsopt_opt_mode_default=None,
        convert_files=False,
        convert_files_explicit=False,
        thresh_post_forward=None,
        flatten_explicit=False,
        flatten=None,
        skip_final_freq=False,
        skip_final_freq_explicit=False,
        stop_plateau=stop_plateau,
        stop_plateau_thresh=None,
        stop_plateau_window=None,
    )
    assert overrides == expected
