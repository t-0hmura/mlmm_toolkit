"""FAIR-Chem owns its lazy-initialization device transfer."""

from __future__ import annotations

import numpy as np
import pytest
import torch
from ase import Atoms

from mlmm.backends import mlmm_calc
from mlmm.backends.mlmm_calc import _UMABackend


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
