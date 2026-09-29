"""The IRC initial Hessian follows ``irc.hessian_init``."""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pytest
import torch
from click.testing import CliRunner

from pysisyphus.Geometry import Geometry


class _EulerPCReached(RuntimeError):
    pass


def _run_irc_until_eulerpc(tmp_path, monkeypatch, irc_section=None):
    """Run the irc CLI with stub inputs up to the EulerPC constructor."""
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
            "--out-dir", str(tmp_path / "out")]
    if irc_section is not None:
        config = tmp_path / "config.yaml"
        config.write_text(json.dumps({"irc": irc_section}), encoding="utf-8")
        args += ["--config", str(config)]
    result = CliRunner().invoke(irc_module.cli, args)
    assert "EulerPCReached" in result.output, result.output
    captured["device_line"] = "[device] IRC Hessian device" in result.output
    return hessian_calls, captured


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
