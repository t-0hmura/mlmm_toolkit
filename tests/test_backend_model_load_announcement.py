"""Pins the announcement that brackets the first load of each MLIP model."""

from __future__ import annotations

from pathlib import Path

import pytest
import torch

from mlmm.backends import mlmm_calc
from mlmm.core.output import mlip_model_label


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
    expected = {
        "uma": "UMA / UMA-S-1.2 (OMol)",
        "mace": "MACE / MACE-OMOL-0",
    }[backend]
    assert f"[backend] Preparing MLIP model ({expected})..." in out
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
    assert "(UMA / UMA-S-1.2 (OMol))" in out
    assert "(UMA / UMA-M-1.1 (OMol))" in out


def test_factory_forwards_analytical_hessian_requirement(monkeypatch) -> None:
    monkeypatch.setattr(mlmm_calc, "_UMABackend", _StubBackend)

    created = _create("uma", analytical_hessian=True)

    assert created.kwargs["analytical_hessian"] is True


@pytest.mark.parametrize(
    ("mode", "expected"),
    [("Analytical", True), ("FiniteDifference", False)],
)
def test_mlmmcore_forwards_effective_hessian_mode_to_backend(
    monkeypatch, mode, expected
) -> None:
    captured = {}

    def fake_create(_backend, **kwargs):
        captured.update(kwargs)
        return _StubBackend(**kwargs)

    monkeypatch.setattr(mlmm_calc, "_create_ml_backend", fake_create)
    monkeypatch.setattr(
        mlmm_calc, "hessianffCalculator", lambda **_kwargs: object()
    )
    example = Path(__file__).resolve().parents[1] / "examples" / "toy_system"
    core = mlmm_calc.MLMMCore(
        input_pdb=str(example / "r_complex.pdb"),
        real_parm7=str(example / "p_toy.parm7"),
        model_pdb=str(example / "ml_region_r.pdb"),
        model_charge=-1,
        ml_device="cpu",
        hessian_calc_mode=mode,
    )
    try:
        assert captured["analytical_hessian"] is expected
    finally:
        core.cleanup()


def test_uma_task_is_part_of_the_model_announcement(capsys, monkeypatch) -> None:
    monkeypatch.setattr(mlmm_calc, "_UMABackend", _StubBackend)
    _create("uma", uma_model="uma-s-1p2", uma_task_name="omol")
    _create("uma", uma_model="uma-s-1p2", uma_task_name="omat")
    out = capsys.readouterr().out
    assert "UMA-S-1.2 (OMol)" in out
    assert "UMA-S-1.2 (OMat)" in out


def test_uma_default_and_explicit_omol_share_one_announcement(
    capsys, monkeypatch
) -> None:
    monkeypatch.setattr(mlmm_calc, "_UMABackend", _StubBackend)
    _create("uma", uma_model="uma-s-1p2")
    _create("uma", uma_model="uma-s-1p2", uma_task_name="omol")
    assert capsys.readouterr().out.count("Preparing MLIP model") == 1


def test_public_model_labels_canonicalize_supported_aliases() -> None:
    assert mlip_model_label("mace", "off:small") == "MACE-OFF23-small"
    assert (
        mlip_model_label("orb", "orb-v3-conservative-omol")
        == "ORB-v3-conservative-OMol"
    )


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
