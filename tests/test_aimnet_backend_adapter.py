"""AIMNet2 backend follows the public AIMNet 0.2 ASE adapter contract."""

from __future__ import annotations

import sys
from types import ModuleType

import numpy as np
import torch
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from mlmm.backends import mlmm_calc


class _BaseCalculator:
    def __init__(self, *, model: str, device: str):
        self.model = model
        self.device = device


class _ASECalculator(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, *, base_calc, charge: int, mult: int):
        super().__init__()
        self.base_calc = base_calc
        self.charge = charge
        self.mult = mult

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.results = {
            "energy": -1.0,
            "forces": np.zeros((len(self.atoms), 3), dtype=np.float64),
        }

    def get_hessian(self, atoms=None):
        return np.eye(3 * len(atoms), dtype=np.float64)


def test_aimnet_backend_uses_public_ase_adapter(monkeypatch) -> None:
    aimnet = ModuleType("aimnet")
    aimnet.__path__ = []
    calculators = ModuleType("aimnet.calculators")
    calculators.AIMNet2Calculator = _BaseCalculator
    calculators.AIMNet2ASE = _ASECalculator
    aimnet.calculators = calculators
    monkeypatch.setitem(sys.modules, "aimnet", aimnet)
    monkeypatch.setitem(sys.modules, "aimnet.calculators", calculators)
    monkeypatch.setattr(mlmm_calc, "HAS_AIMNET2", True)

    backend = mlmm_calc._AIMNet2Backend(
        aimnet2_model="aimnet2",
        model_charge=-1,
        model_mult=1,
        ml_device=torch.device("cpu"),
    )
    atoms = Atoms("H2", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
    energy, forces, opaque = backend.eval(atoms)
    hessian = backend.hessian_analytical(opaque, 2, dtype=torch.float64)

    assert isinstance(backend._ase_calc, _ASECalculator)
    assert backend._ase_calc.base_calc is backend._aimnet_base_calc
    assert backend._ase_calc.charge == -1
    assert energy == -1.0
    assert forces.shape == (2, 3)
    assert torch.equal(hessian.reshape(6, 6), torch.eye(6, dtype=torch.float64))
