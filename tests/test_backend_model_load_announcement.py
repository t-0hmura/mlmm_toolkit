"""Pins the announcement that brackets the first load of each MLIP model."""

from __future__ import annotations

import pytest
import torch

from mlmm.backends import mlmm_calc


class _StubBackend:
    def __init__(self, **kwargs):
        self.kwargs = kwargs


@pytest.fixture(autouse=True)
def _fresh_announcement_state(monkeypatch):
    monkeypatch.setattr(mlmm_calc, "_ANNOUNCED_MODEL_LOADS", set())


def _create(backend: str = "uma", **kwargs):
    return mlmm_calc._create_ml_backend(
        backend, ml_device=torch.device("cpu"), **kwargs
    )


@pytest.mark.parametrize(
    "backend, stub_name, model_kwarg, model",
    [
        ("uma", "_UMABackend", "uma_model", "uma-s-1p2"),
        ("mace", "_MACEBackend", "mace_model", "MACE-OMOL-0"),
    ],
)
def test_first_load_is_bracketed_then_silent(
    backend, stub_name, model_kwarg, model, capsys, monkeypatch
) -> None:
    monkeypatch.setattr(mlmm_calc, stub_name, _StubBackend)
    created = _create(backend, **{model_kwarg: model})
    assert isinstance(created, _StubBackend)
    out = capsys.readouterr().out
    assert f"[backend] Preparing MLIP model ({backend} / {model})..." in out
    assert "[backend] Done." in out
    assert out.index("Preparing MLIP model") < out.index("[backend] Done.")
    lines = out.splitlines()
    preparing = next(i for i, line in enumerate(lines) if "Preparing MLIP model" in line)
    done = lines.index("[backend] Done.")
    assert preparing == 0 or lines[preparing - 1] == ""
    assert done + 1 < len(lines) and lines[done + 1] == ""

    _create(backend, **{model_kwarg: model})
    assert "Preparing MLIP model" not in capsys.readouterr().out


def test_each_model_is_announced_once(capsys, monkeypatch) -> None:
    monkeypatch.setattr(mlmm_calc, "_UMABackend", _StubBackend)
    _create("uma", uma_model="uma-s-1p2")
    _create("uma", uma_model="uma-m-1p1")
    out = capsys.readouterr().out
    assert out.count("Preparing MLIP model") == 2
    assert "(uma / uma-s-1p2)" in out and "(uma / uma-m-1p1)" in out


def test_failed_load_is_announced_again(capsys, monkeypatch) -> None:
    class _Boom:
        def __init__(self, **kwargs):
            raise RuntimeError("checkpoint download failed")

    monkeypatch.setattr(mlmm_calc, "_UMABackend", _Boom)
    with pytest.raises(RuntimeError):
        _create("uma", uma_model="uma-s-1p2")
    out = capsys.readouterr().out
    assert "Preparing MLIP model" in out
    assert "[backend] Done." not in out

    monkeypatch.setattr(mlmm_calc, "_UMABackend", _StubBackend)
    _create("uma", uma_model="uma-s-1p2")
    retry = capsys.readouterr().out
    assert "Preparing MLIP model" in retry and "[backend] Done." in retry
