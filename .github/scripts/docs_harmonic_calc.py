"""Dependency-light ASE calculator used by documentation smoke tests."""

import numpy as np
from ase.calculators.calculator import Calculator, all_changes


class CenteredHarmonicCalculator(Calculator):
    implemented_properties = ["energy", "forces"]

    def calculate(
        self,
        atoms=None,
        properties=("energy", "forces"),
        system_changes=all_changes,
    ):
        super().calculate(atoms, properties, system_changes)
        positions = np.asarray(self.atoms.get_positions(), dtype=float)
        displacement = positions - positions.mean(axis=0, keepdims=True)
        self.results = {
            "energy": 0.5 * float(np.sum(displacement * displacement)),
            "forces": -displacement,
        }


def get_calculator(**kwargs):
    return CenteredHarmonicCalculator()
