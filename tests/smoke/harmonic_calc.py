"""Minimal energy/forces-only ASE calculator for release smoke tests."""

import numpy as np
from ase.calculators.calculator import Calculator, all_changes


class HarmonicCalculator(Calculator):
    implemented_properties = ["energy", "forces"]

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        if self.atoms.info.get("charge") != -1:
            raise RuntimeError("custom calculator did not receive charge=-1 in atoms.info")
        if self.atoms.info.get("mult") != 1 or self.atoms.info.get("spin") != 1:
            raise RuntimeError("custom calculator did not receive singlet multiplicity in atoms.info")
        positions = np.asarray(self.atoms.get_positions(), dtype=float)
        displacement = positions - positions.mean(axis=0, keepdims=True)
        self.results = {
            "energy": 0.5 * float(np.sum(displacement * displacement)),
            "forces": -displacement,
        }


def get_calculator(charge=0, spin=1, device="auto", **kwargs):
    if charge != -1 or spin != 1:
        raise RuntimeError(f"unexpected factory state: charge={charge}, spin={spin}")
    if device not in {"cpu", "cuda"}:
        raise RuntimeError(f"unexpected factory device: {device}")
    return HarmonicCalculator()
