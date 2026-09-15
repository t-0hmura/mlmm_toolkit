"""FAIR-Chem owns its lazy-initialization device transfer."""

from __future__ import annotations

import numpy as np
import pytest
import torch
from ase import Atoms

from mlmm.backends import mlmm_calc
from mlmm.backends.mlmm_calc import _UMABackend

# UMA analytical Hessians need the EFS force graph, not backbone training.
import torch


def _scope_model(wrapped=True):
    inner = torch.nn.Module()
    inner.backbone = torch.nn.Sequential(torch.nn.Linear(3, 3), torch.nn.Dropout(0.2))
    inner.backbone.composition_dropout = 0.1
    head = torch.nn.Sequential(torch.nn.Linear(3, 1), torch.nn.Dropout(0.3))
    entry = torch.nn.Module()
    entry.head = head
    inner.output_heads = torch.nn.ModuleDict({
        "energyandforcehead": entry,
        "other": torch.nn.Linear(3, 1),
    })
    model = torch.nn.Module() if wrapped else inner
    if wrapped:
        model.module = inner
    model.eval()
    return model, inner, head


def _scope_state(model):
    return (
        tuple((name, module.training, getattr(module, "p", None))
              for name, module in model.named_modules()),
        tuple((name, parameter.requires_grad)
              for name, parameter in model.named_parameters()),
    )


@pytest.mark.parametrize("wrapped", [False, True])
@pytest.mark.parametrize("fail", [False, True])
def test_uma_head_scope_restores_mixed_state(wrapped, fail):
    model, inner, head = _scope_model(wrapped)
    # A saved training backbone is made eval during derivatives, then restored.
    model.train()
    inner.backbone[0].eval()
    head[0].eval()
    next(model.parameters()).requires_grad_(False)
    before = _scope_state(model)
    try:
        with _uma_analytical_head_scope(model):
            assert not model.training
            assert all(not module.training for module in inner.backbone.modules())
            assert not inner.output_heads["other"].training
            assert head.training and head[0].training
            assert not head[1].training and head[1].p == 0.0
            assert inner.backbone.composition_dropout == 0.1
            assert inner.backbone[1].p == 0.2
            assert not any(parameter.requires_grad for parameter in model.parameters())
            if fail:
                raise ValueError("injected derivative failure")
    except ValueError as exc:
        assert fail and str(exc) == "injected derivative failure"
    assert _scope_state(model) == before
    assert inner.backbone.composition_dropout == 0.1


@pytest.mark.parametrize("layout", ["missing", "not_module", "backbone_alias"])
def test_uma_head_scope_rejects_unknown_layout_without_mutation(layout):
    model, inner, head = _scope_model()
    if layout == "missing":
        del inner.output_heads["energyandforcehead"]
    elif layout == "not_module":
        del inner.output_heads["energyandforcehead"].head
        inner.output_heads["energyandforcehead"].head = object()
    else:
        inner.output_heads["energyandforcehead"].head = inner.backbone
    before = _scope_state(model)
    with pytest.raises(RuntimeError, match="UMA analytical Hessian requires"):
        with _uma_analytical_head_scope(model):
            pytest.fail("Unsupported UMA layout must not fall back to global training.")
    assert _scope_state(model) == before


def test_uma_head_scope_restores_after_partial_preparation(monkeypatch):
    model, _, head = _scope_model()
    before = _scope_state(model)
    original_train = head.train

    def interrupted_train(mode=True):
        original_train(mode)
        if mode:
            raise ValueError("injected head preparation failure")
        return head

    monkeypatch.setattr(head, "train", interrupted_train)
    with pytest.raises(ValueError, match="head preparation failure"):
        with _uma_analytical_head_scope(model):
            pytest.fail("Preparation did not fail.")
    assert _scope_state(model) == before


class _ScopePredictor:
    def __init__(self, model, inner, head, fail=False):
        self.model, self.inner, self.head = model, inner, head
        self.fail = fail
        self.derivative_calls = 0

    def predict(self, batch):
        if self.head.training:
            self.derivative_calls += 1
            assert not self.model.training
            assert all(not module.training for module in self.inner.backbone.modules())
            assert not self.inner.output_heads["other"].training
            assert self.inner.backbone.composition_dropout == 0.1
            assert not any(p.requires_grad for p in self.model.parameters())
            if self.fail:
                raise ValueError("injected predictor failure")
        # Deterministic stand-in for UMA's training-only functional routing.
        coefficient = 11.0 if self.inner.backbone.training else 1.0
        energy = coefficient * (batch.pos ** 2).sum()
        # Like UMA, the ordinary force evaluation consumes the energy graph;
        # EFS-head training must retain it for the subsequent Hessian.
        forces = -torch.autograd.grad(
            energy, batch.pos,
            create_graph=self.head.training, retain_graph=self.head.training,
        )[0]
        return {"energy": energy.reshape(1), "forces": forces}


@pytest.mark.parametrize("fail", [False, True])
def test_uma_analytical_owner_uses_head_scope_and_restores(fail):
    from types import SimpleNamespace

    model, inner, head = _scope_model()
    next(model.parameters()).requires_grad_(False)
    before = _scope_state(model)
    predictor = _ScopePredictor(model, inner, head, fail=fail)
    batch = SimpleNamespace(pos=torch.tensor([[0.2, 0.3, 0.4]], dtype=torch.float64,
                                            requires_grad=True))
    compute = _scope_owner_call(predictor, batch)
    if fail:
        with pytest.raises(ValueError, match="injected predictor failure"):
            compute()
    else:
        hessian = compute()
        torch.testing.assert_close(hessian.reshape(3, 3), 2.0 * torch.eye(3, dtype=torch.float64))
    assert predictor.derivative_calls == 1
    assert _scope_state(model) == before

from mlmm.backends.mlmm_calc import _uma_analytical_head_scope


def _scope_owner_call(predictor, batch):
    backend = object.__new__(_UMABackend)
    backend.predictor = predictor
    backend._device = torch.device("cpu")
    return lambda: backend.hessian_analytical(batch, 1, dtype=torch.float64)



@pytest.mark.parametrize(
    ("precision", "expected_dtype"),
    [("fp32", "float32"), ("fp64", "float64")],
)
def test_analytical_mode_uses_noncompiled_precision_matched_settings(
    monkeypatch, precision, expected_dtype
) -> None:
    captured = {}

    class FakeSettings:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

    class FakePretrained:
        @staticmethod
        def get_predict_unit(_model, **kwargs):
            captured.update(kwargs)
            return object()

    monkeypatch.setattr(mlmm_calc, "HAS_FAIRCHEM", True)
    monkeypatch.setattr(mlmm_calc, "_UMAInferenceSettings", FakeSettings)
    monkeypatch.setattr(
        mlmm_calc, "pretrained_mlip", FakePretrained(), raising=False
    )
    monkeypatch.setattr(mlmm_calc, "AtomicData", object, raising=False)
    monkeypatch.setattr(mlmm_calc, "data_list_collater", object(), raising=False)

    _UMABackend(
        ml_device=torch.device("cpu"),
        precision=precision,
        analytical_hessian=True,
    )

    settings = captured["inference_settings"]
    assert settings.kwargs == {
        "compile": False,
        "base_precision_dtype": expected_dtype,
    }


def test_serial_uma_prediction_supplies_a_cpu_batch() -> None:
    class FakeData:
        dataset = None

    class FakeAtomicData:
        @staticmethod
        def from_ase(*_args, **_kwargs):
            return FakeData()

    class FakeBatch:
        def __init__(self):
            self.pos = torch.zeros((1, 3), dtype=torch.float32)

        def to(self, _device):
            pytest.fail("MLMM must leave FAIR-Chem input device transfer to FAIR-Chem")

    class FakePredictor:
        @staticmethod
        def predict(batch):
            assert batch.pos.device.type == "cpu"
            return {
                "energy": torch.zeros(1),
                "forces": torch.zeros((1, 3)),
            }

    backend = object.__new__(_UMABackend)
    backend._AtomicData = FakeAtomicData
    backend._data_list_collater = lambda *_args, **_kwargs: FakeBatch()
    backend.predictor = FakePredictor()
    backend._device = torch.device("cuda")
    backend.precision = "fp32"
    backend.uma_task_name = "omol"
    backend.model_charge = 0
    backend.model_mult = 1
    backend.parallel_predict = False
    backend._uma_max_neigh = None
    backend._uma_radius = None

    energy, forces, batch = backend.eval(
        Atoms("H", positions=np.zeros((1, 3))), need_grad=True
    )

    assert energy == 0.0
    assert forces.shape == (1, 3)
    assert batch.pos.requires_grad
