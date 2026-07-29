from __future__ import annotations

from types import ModuleType
import sys

import pytest
import torch

from mlmm.backends import mlmm_calc


@pytest.mark.parametrize("dtype", ["float32", "float64"])
def test_mace_anicc_uses_supported_factory_signature(
    monkeypatch,
    dtype: str,
) -> None:
    captured = {}
    calculators = ModuleType("mace.calculators")

    def mace_anicc(**kwargs):
        captured["factory"] = kwargs
        return "raw-model" if kwargs.get("return_raw_model") else "calculator"

    class FakeMACECalculator:
        def __init__(self, **kwargs):
            captured["calculator"] = kwargs

    calculators.mace_anicc = mace_anicc
    calculators.mace_mp = lambda **_kwargs: object()
    calculators.mace_off = lambda **_kwargs: object()
    calculators.mace_omol = lambda **_kwargs: object()
    calculators.MACECalculator = FakeMACECalculator
    mace_package = ModuleType("mace")
    mace_package.__path__ = []

    monkeypatch.setattr(mlmm_calc, "HAS_MACE", True)
    monkeypatch.setitem(sys.modules, "mace", mace_package)
    monkeypatch.setitem(sys.modules, "mace.calculators", calculators)

    mlmm_calc._MACEBackend(
        mace_model="MACE-ANICC",
        mace_dtype=dtype,
        ml_device=torch.device("cpu"),
    )

    if dtype == "float64":
        assert captured == {"factory": {"device": "cpu"}}
    else:
        assert captured == {
            "factory": {"device": "cpu", "return_raw_model": True},
            "calculator": {
                "models": "raw-model",
                "device": "cpu",
                "default_dtype": "float32",
            },
        }
