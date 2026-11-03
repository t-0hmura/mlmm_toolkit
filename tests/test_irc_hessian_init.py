"""The IRC initial Hessian follows ``irc.hessian_init``."""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pytest
import torch
from click.testing import CliRunner

from pysisyphus.Geometry import Geometry
from _repo_files import require_repo_files


class _EulerPCReached(RuntimeError):
    pass


def _run_irc_until_eulerpc(tmp_path, monkeypatch, irc_section=None, *, finished=None, extra_args=()):
    """Run the irc CLI with stub inputs up to the EulerPC constructor.

    With ``finished`` (an object standing for an EulerPC that has run), the
    constructor returns it and the command continues to its outputs.
    """
    from mlmm.workflows import irc as irc_module

    source, parm = tmp_path / "input.pdb", tmp_path / "input.parm7"
    source.write_text("Structure parsing is replaced for this test.\n")
    parm.write_text("No MM topology is read.\n")
    prepared = SimpleNamespace(
        original_path=source, source_path=source, geom_path=source,
        cleanup=lambda: None,
    )
    hessian_calls = []
    captured = {}

    def load_geometry(_path, **_kwargs):
        return Geometry(
            ("O", "H", "H"),
            np.array([0.0, 0.0, 0.0, 1.8, 0.0, 0.0, -0.4, 1.7, 0.0]),
            coord_type="cart",
        )

    def fake_hessian(geometry, _calc_cfg, _device, **_kwargs):
        hessian_calls.append(geometry)
        return torch.eye(geometry.cart_coords.size, dtype=torch.float64), None

    def fake_eulerpc(geometry, **kwargs):
        captured["hessian_init"] = kwargs.get("hessian_init")
        captured["seeded"] = geometry._hessian is not None
        if finished is not None:
            return finished
        raise _EulerPCReached("EulerPCReached: no IRC integration requested.")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(irc_module, "prepare_input_structure", lambda *_a, **_k: prepared)
    monkeypatch.setattr(irc_module, "resolve_charge_spin_or_raise", lambda *_a, **_k: (0, 1))
    monkeypatch.setattr(irc_module, "resolve_ml_layer_assignment", lambda **_k: (source, None))
    monkeypatch.setattr(irc_module, "geom_loader", load_geometry)
    monkeypatch.setattr(irc_module, "mlmm", lambda **_k: SimpleNamespace(core=None))
    monkeypatch.setattr(irc_module, "_calc_full_hessian_torch", fake_hessian)
    monkeypatch.setattr(irc_module, "_torch_device", lambda *_a: torch.device("cpu"))
    monkeypatch.setattr(irc_module, "EulerPC", fake_eulerpc)

    args = ["-i", str(source), "--parm", str(parm), "-q", "0", "-m", "1",
            "--out-dir", str(tmp_path / "out"), *extra_args]
    if irc_section is not None:
        config = tmp_path / "config.yaml"
        config.write_text(json.dumps({"irc": irc_section}), encoding="utf-8")
        args += ["--config", str(config)]
    result = CliRunner().invoke(irc_module.cli, args)
    if finished is None:
        assert "EulerPCReached" in result.output, result.output
    else:
        captured["result"] = result
    captured["device_line"] = "[device] IRC Hessian device" in result.output
    return hessian_calls, captured


class _FinishedEulerPC:
    """What the irc command reads from an EulerPC that has run both directions."""

    forward = backward = True
    forward_is_converged = backward_is_converged = False
    forward_mw_hessian = mw_hessian = None
    never_stop = False
    max_cycles = 100
    init_hessian_shape = (9, 9)
    _act_atoms = (0, 1, 2)
    forward_integration_converged = True
    backward_integration_converged = False
    forward_integration_stop_reason = "energy_minimum"
    backward_integration_stop_reason = "max_cycles"
    forward_downhill_departure_valid = True
    backward_downhill_departure_valid = False
    forward_energy_increased = False
    backward_energy_increased = True

    def __init__(self, out_dir):
        self.out_dir = out_dir
        frames = [np.zeros(9) for _ in range(5)]
        self.forward_energies = self.backward_energies = [-1.0, -1.1, -1.2, -1.3, -1.4]
        self.forward_coords = self.backward_coords = frames
        self.all_energies = [-1.4, -1.2, -1.0, -1.2, -1.4]
        self.ts_energy = -1.0

    def run(self):
        pass

    def get_path_for_fn(self, filename):
        return str(self.out_dir / filename)


def test_irc_result_json_carries_every_direction_field_that_all_reads(
    tmp_path, monkeypatch
) -> None:
    from mlmm.workflows.all import _read_irc_outcome

    out_dir = tmp_path / "out"
    finished = _FinishedEulerPC(out_dir)
    _hessian_calls, captured = _run_irc_until_eulerpc(
        tmp_path, monkeypatch, finished=finished, extra_args=("--out-json",)
    )

    result = captured["result"]
    assert result.exit_code == 0, result.output
    outcome = _read_irc_outcome(out_dir)
    for direction in ("forward", "backward"):
        assert outcome[f"{direction}_requested"] is True
        for suffix in ("integration_converged", "integration_stop_reason",
                       "downhill_departure_valid", "energy_increased"):
            assert outcome[f"{direction}_{suffix}"] == getattr(finished, f"{direction}_{suffix}")
        assert outcome[f"{direction}_short_branch"] is False
        assert outcome[f"n_frames_{direction}"] == 5


@pytest.mark.parametrize(("child", "prefix"), [(False, "--"), (True, "--irc-")])
def test_few_frames_warning_names_the_options_of_the_running_command(
    tmp_path, monkeypatch, child, prefix
) -> None:
    from mlmm.core import utils

    finished = _FinishedEulerPC(tmp_path / "out")
    finished.forward_energies = [-1.0, -1.1]
    finished.forward_coords = [np.zeros(9) for _ in range(2)]
    # Under `all` the parent's --irc-* options are the ones a user can pass.
    monkeypatch.setattr(utils, "_CHILD_MODE", child)
    _hessian_calls, captured = _run_irc_until_eulerpc(tmp_path, monkeypatch, finished=finished)

    output = captured["result"].output
    assert "IRC stopped after only a few frames in forward." in output
    assert f"for example {prefix}step-size 0.05." in output
    assert f"consider {prefix}never-stop; it is opt-in." in output


def test_irc_endpoint_bond_changes_use_the_yaml_bond_section(tmp_path, monkeypatch) -> None:
    from mlmm.domain import bond_changes

    seen = []
    real = bond_changes.compare_structures

    def record(g1, g2, **kwargs):
        seen.append(kwargs)
        return real(g1, g2, **kwargs)

    class WithEndpoints(_FinishedEulerPC):
        def run(self):
            for name in ("finished_first.xyz", "finished_last.xyz"):
                (self.out_dir / name).write_text(
                    "3\n\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n", encoding="utf-8")

    monkeypatch.setattr(bond_changes, "compare_structures", record)
    config = tmp_path / "config.yaml"
    config.write_text(json.dumps({"bond": {"bond_factor": 1.3, "metal_delta_fraction": 0.25}}),
                      encoding="utf-8")
    _run_irc_until_eulerpc(
        tmp_path, monkeypatch, finished=WithEndpoints(tmp_path / "out"),
        extra_args=("--out-json", "--config", str(config)),
    )
    (call,) = seen
    assert (call["bond_factor"], call["metal_delta_fraction"]) == (1.3, 0.25)


def test_cached_ts_hessian_seeds_the_irc_and_is_logged_at_the_detail_tier(
    tmp_path, monkeypatch
) -> None:
    from mlmm.io import hessian_cache
    from mlmm.workflows import irc as irc_module

    monkeypatch.setattr(hessian_cache, "identity_from_context", lambda *_a, **_k: {})
    monkeypatch.setattr(
        hessian_cache, "load_matching",
        lambda *_a, **_k: {"hessian": torch.eye(9, dtype=torch.float64), "active_dofs": None},
    )
    messages = []
    monkeypatch.setattr(irc_module, "emit", lambda text="", **tags: messages.append((text, tags)))

    hessian_calls, captured = _run_irc_until_eulerpc(tmp_path, monkeypatch)

    assert hessian_calls == []
    assert captured["seeded"] is True
    reuse = [tags for text, tags in messages if "Reusing cached TS Hessian from tsopt" in text]
    assert reuse == [{"detail": True}]


def test_cached_ts_hessian_of_another_basis_is_dropped_and_logged_at_the_detail_tier(
    tmp_path, monkeypatch
) -> None:
    from mlmm.io import hessian_cache
    from mlmm.workflows import irc as irc_module

    monkeypatch.setattr(hessian_cache, "identity_from_context", lambda *_a, **_k: {})
    monkeypatch.setattr(
        hessian_cache, "load_matching",
        lambda *_a, **_k: {"hessian": torch.eye(3, dtype=torch.float64), "active_dofs": [0, 1, 2]},
    )
    messages = []
    monkeypatch.setattr(irc_module, "emit", lambda text="", **tags: messages.append((text, tags)))

    hessian_calls, _captured = _run_irc_until_eulerpc(tmp_path, monkeypatch)

    assert len(hessian_calls) == 1
    mismatch = [tags for text, tags in messages if "active-DOF basis does not match" in text]
    assert mismatch == [{"detail": True}]


def test_default_hessian_init_computes_and_seeds_the_initial_hessian(
    tmp_path, monkeypatch
) -> None:
    hessian_calls, captured = _run_irc_until_eulerpc(tmp_path, monkeypatch)

    assert len(hessian_calls) == 1
    assert captured == {"hessian_init": "calc", "seeded": True, "device_line": True}


@pytest.mark.parametrize("hessian_init", ["unit", "ts_hessian.h5"])
def test_other_hessian_init_is_left_to_eulerpc(
    tmp_path, monkeypatch, hessian_init
) -> None:
    hessian_calls, captured = _run_irc_until_eulerpc(
        tmp_path, monkeypatch, {"hessian_init": hessian_init}
    )

    assert hessian_calls == []
    assert captured == {
        "hessian_init": hessian_init, "seeded": False, "device_line": False,
    }


@pytest.mark.parametrize("command", ["irc", "freq"])
def test_hess_device_cuda_without_cuda_is_rejected_before_dry_run(
    tmp_path, monkeypatch, command
) -> None:
    require_repo_files("tests/smoke")
    import importlib
    from pathlib import Path

    module = importlib.import_module(f"mlmm.workflows.{command}")
    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    smoke = Path(__file__).resolve().parent / "smoke"

    result = CliRunner().invoke(module.cli, [
        "-i", str(smoke / "p_complex_layered.pdb"),
        "--parm", str(smoke / "p_complex.parm7"),
        "-q", "-1", "-m", "1", "--out-dir", str(tmp_path / "out"),
        "--hess-device", "cuda", "--dry-run",
    ])

    assert result.exit_code == 1, result.output
    assert "--hess-device cuda was requested but no CUDA device is available" in result.output
