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


def _lazy_preparation_backend(monkeypatch, *, device, precision, workers=1):
    """Model the PR298 preparation order with CPU tensors and logical devices."""
    from contextlib import nullcontext
    from types import SimpleNamespace

    dtype = torch.float64 if precision == "fp64" else torch.float32
    events = []

    class Batch:
        def __init__(self, data):
            self.pos = data.pos
            self.charge, self.spin = data.charge, data.spin
            self.device = self.charge_device = self.spin_device = "cpu"

        def to(self, target):
            # No CUDA allocation: these labels track the ordering contract.
            self.device = self.charge_device = self.spin_device = str(target)
            events.append(("batch_to", str(target)))
            return self

    class LazyPredictor:
        def __init__(self):
            if workers == 1:
                self.model = torch.nn.Module()
                self.model.backbone = torch.nn.Module()
            self.model_device = "cpu"
            self.initialized = False

        def predict(self, batch):
            if not self.initialized:
                assert (batch.device == batch.charge_device == batch.spin_device
                        == self.model_device), "lazy preparation device mismatch"
                events.append(("prepare", self.model_device))
                self.model_device = device
                self.initialized = True
            assert batch.device == "cpu", "the predictor must receive a host batch"
            assert (batch.charge, batch.spin) == (-1, 2)
            assert batch.pos.dtype == dtype and batch.pos.device.type == "cpu"
            batch.to(device)  # FAIR-Chem's transfer follows its lazy preparation.
            assert batch.device == self.model_device
            events.append(("predict", batch.device))
            return {"energy": (batch.pos ** 2).sum().reshape(1),
                    "forces": -2 * batch.pos}

    predictor = LazyPredictor()  # Deliberately exposes no move_to_device hook.

    class AtomicData:
        @staticmethod
        def from_ase(atoms, **kwargs):
            keys = kwargs["r_data_keys"]
            return SimpleNamespace(
                pos=torch.tensor(atoms.positions, dtype=kwargs["target_dtype"]),
                charge=atoms.info.get("charge", 0) if "charge" in keys else 0,
                spin=atoms.info.get("spin", 0) if "spin" in keys else 0,
                dataset=None,
            )

    def collate(data, **kwargs):
        assert len(data) == 1 and data[0].dataset == "omol"
        assert kwargs == {"otf_graph": True}
        return Batch(data[0])

    def serial_factory(model, **kwargs):
        assert model == "uma-s-1p2" and kwargs["device"] == device
        events.append(("construct", "serial"))
        return predictor

    def parallel_factory(**kwargs):
        assert kwargs["device"] == device and kwargs["num_workers"] == workers
        events.append(("construct", "parallel"))
        return predictor

    pretrained = SimpleNamespace(
        get_predict_unit=serial_factory,
        pretrained_checkpoint_path_from_name=lambda _model: "/unused-checkpoint",
        get_reference_energies=lambda *_args, **_kwargs: {},
    )
    monkeypatch.setattr(mlmm_calc, "HAS_FAIRCHEM", True)
    monkeypatch.setattr(mlmm_calc, "pretrained_mlip", pretrained, raising=False)
    monkeypatch.setattr(mlmm_calc, "AtomicData", AtomicData, raising=False)
    monkeypatch.setattr(mlmm_calc, "data_list_collater", collate, raising=False)
    monkeypatch.setattr(mlmm_calc, "_UMAInferenceSettings", SimpleNamespace)
    monkeypatch.setattr(mlmm_calc, "ParallelMLIPPredictUnit", parallel_factory)
    monkeypatch.setattr(mlmm_calc, "guess_inference_settings", lambda name: name)
    monkeypatch.setattr(mlmm_calc.torch.cuda, "device", lambda _device: nullcontext())
    backend = _UMABackend(uma_model="uma-s-1p2", model_charge=-1, model_mult=2,
                          ml_device=torch.device(device), precision=precision, workers=workers)
    assert events == [("construct", "parallel" if workers > 1 else "serial")]
    assert not predictor.initialized
    return backend, predictor, events


@pytest.mark.parametrize("entry", ["eval", "energy", "forces_tensor"])
@pytest.mark.parametrize("device", ["cpu", "cuda"])
@pytest.mark.parametrize("precision", ["fp32", "fp64"])
def test_fresh_uma_entries_prepare_before_predictor_owned_transfer(
    monkeypatch, entry, device, precision,
):
    backend, predictor, events = _lazy_preparation_backend(
        monkeypatch, device=device, precision=precision,
    )
    positions = np.array([[0., 0., 0.], [.757, .586, 0.], [-.757, .586, 0.]])
    for displacement in (0., .001):
        current = positions.copy()
        current[1, 0] += displacement
        atoms = Atoms("OHH", positions=current)
        if entry == "eval":
            energy, force, batch = backend.eval(atoms, need_grad=True)
            assert batch.pos.requires_grad
            assert energy == pytest.approx((current ** 2).sum(), rel=1e-6)
            np.testing.assert_allclose(force, -2 * current, rtol=1e-6)
        elif entry == "energy":
            assert backend.energy(atoms) == pytest.approx((current ** 2).sum(), rel=1e-6)
        else:
            force = backend.forces_tensor(atoms)
            assert not force.requires_grad
            assert force.dtype == (torch.float64 if precision == "fp64" else torch.float32)
            np.testing.assert_allclose(force.numpy(), -2 * current, rtol=1e-6)
    assert predictor.initialized
    assert events[1:] == [("prepare", "cpu"), ("batch_to", device), ("predict", device),
                          ("batch_to", device), ("predict", device)]


def test_parallel_uma_constructor_keeps_host_handoff(monkeypatch):
    backend, _, events = _lazy_preparation_backend(
        monkeypatch, device="cuda", precision="fp32", workers=2,
    )
    assert backend.parallel_predict and not backend.supports_analytical_hessian
    _, force, _ = backend.eval(Atoms("OHH", positions=np.ones((3, 3))))
    np.testing.assert_array_equal(force, -2 * np.ones((3, 3)))
    assert events == [("construct", "parallel"), ("prepare", "cpu"),
                      ("batch_to", "cuda"), ("predict", "cuda")]


def test_lazy_preparation_control_detects_premature_device_transfer(monkeypatch):
    backend, predictor, _ = _lazy_preparation_backend(
        monkeypatch, device="cuda", precision="fp32",
    )
    original_collate = backend._data_list_collater

    def premature_collate(*args, **kwargs):
        return original_collate(*args, **kwargs).to("cuda")

    monkeypatch.setattr(backend, "_data_list_collater", premature_collate)
    with pytest.raises(AssertionError, match="lazy preparation device mismatch"):
        backend.energy(Atoms("OHH", positions=np.zeros((3, 3))))
    assert not predictor.initialized
