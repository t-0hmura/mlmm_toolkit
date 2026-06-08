"""Harmonic restraint calculator wrappers (position fix + bond-length bias).

Two pysisyphus-style Calculator wrappers consumed by multiple workflow stages:

- ``HarmonicFixAtoms`` — harmonic *position* restraint on a subset of atoms.
  Used by ``path_opt`` and ``path_search`` (incl. the DMF path optimizer) to pin
  pre-selected atoms with a quadratic well around their reference coordinates.

- ``HarmonicBiasCalculator`` — harmonic *distance* restraint on selected atom
  pairs (1-based / 0-based indexing per caller). Used by ``scan`` / ``scan2d`` /
  ``scan3d`` for bond-length staged scans, and by ``opt`` for ad-hoc distance
  biasing. Wraps a base UMA-style calculator and adds the bias E/F per evaluation.

Both classes are pure-Python (numpy only) and do not import any MLIP SDK, so they
belong with the workflow orchestration layer rather than ``io`` or ``backends``.
"""
from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np
from ase.calculators.calculator import Calculator
from pysisyphus.constants import ANG2BOHR, AU2EV


# eV/Å² → Hartree/Bohr² conversion (= k_evAA * H_EVAA_2_AU)
EV2AU = 1.0 / AU2EV
H_EVAA_2_AU = EV2AU / ANG2BOHR / ANG2BOHR


class HarmonicFixAtoms(Calculator):
    """Harmonic position restraint on a subset of atoms (ASE Calculator).

    Energy = 1/2 * k_fix * Σ_i |r_i − r_i^ref|² (sum over the fixed indices).
    Used in path_opt / path_search (incl. DMF) to pin atoms with a soft well.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, indices, ref_positions, k_fix=300.0):
        super().__init__()
        idx = np.asarray(indices, dtype=int).ravel()
        if idx.size == 0:
            raise ValueError("HarmonicFixAtoms requires at least one index.")
        ref_pos = np.asarray(ref_positions, dtype=float)
        if ref_pos.shape != (idx.size, 3):
            raise ValueError(
                f"ref_positions must have shape ({idx.size}, 3), got {ref_pos.shape}"
            )
        self.indices = idx
        self.ref_positions = ref_pos
        self.k_fix = float(k_fix)

    def calculate(self, atoms, properties, system_changes):
        super().calculate(atoms, properties, system_changes)
        pos = atoms.get_positions().astype(float)
        disp = pos[self.indices] - self.ref_positions
        energy = 0.5 * self.k_fix * np.sum(disp ** 2)
        forces = np.zeros_like(pos, dtype=float)
        forces[self.indices] = -self.k_fix * disp
        self.results = {
            "energy": float(energy),
            "forces": forces,
        }


class HarmonicBiasCalculator:
    """Wrap a base UMA-style calculator with harmonic distance restraints.

    Per-pair bias: Energy = 1/2 * k * (r_ij − target)² for each (i, j, target) tuple.
    Forces are added to the base calculator's force output. Indices are 0-based
    Cartesian atom indices.

    Used by scan / scan2d / scan3d for bond-length staged scans, and by opt for
    ad-hoc distance biasing. Reusable for any future DMF distance-restraint case
    (same pattern: wrap base calc + add per-pair harmonic E/F).
    """

    def __init__(self, base_calc, k: float = 10.0, pairs: Optional[List[Tuple[int, int, float]]] = None):
        self.base = base_calc
        self.k_evAA = float(k)
        self.k_au_bohr2 = self.k_evAA * H_EVAA_2_AU
        self._pairs: List[Tuple[int, int, float]] = list(pairs or [])

    def set_pairs(self, pairs: List[Tuple[int, int, float]]) -> None:
        self._pairs = [(int(i), int(j), float(t)) for (i, j, t) in pairs]

    def _bias_energy_forces_bohr(self, coords_bohr: np.ndarray) -> Tuple[float, np.ndarray]:
        coords = np.array(coords_bohr, dtype=float).reshape(-1, 3)
        n = coords.shape[0]
        E_bias = 0.0
        F_bias = np.zeros((n, 3), dtype=float)
        k = self.k_au_bohr2
        for (i, j, target_ang) in self._pairs:
            if not (0 <= i < n and 0 <= j < n):
                continue
            rij_vec = coords[i] - coords[j]
            rij = float(np.linalg.norm(rij_vec))
            if rij < 1e-14:
                continue
            target_bohr = float(target_ang) * ANG2BOHR
            diff_bohr = rij - target_bohr
            E_bias += 0.5 * k * diff_bohr * diff_bohr
            u = rij_vec / max(rij, 1e-14)
            Fi = -k * diff_bohr * u
            F_bias[i] += Fi
            F_bias[j] -= Fi
        return E_bias, F_bias.reshape(-1)

    def get_forces(self, elem, coords):
        coords_bohr = np.asarray(coords, dtype=float).reshape(-1, 3)
        base = self.base.get_forces(elem, coords_bohr)
        E0 = float(base["energy"])
        F0 = np.asarray(base["forces"], dtype=float).reshape(-1)
        Ebias, Fbias = self._bias_energy_forces_bohr(coords_bohr)
        return {"energy": E0 + Ebias, "forces": F0 + Fbias}

    def get_energy(self, elem, coords):
        coords_bohr = np.asarray(coords, dtype=float).reshape(-1, 3)
        E0 = float(self.base.get_energy(elem, coords_bohr)["energy"])
        Ebias, _ = self._bias_energy_forces_bohr(coords_bohr)
        return {"energy": E0 + Ebias}

    def get_energy_and_forces(self, elem, coords):
        res = self.get_forces(elem, coords)
        return res["energy"], res["forces"]

    def get_energy_and_gradient(self, elem, coords):
        res = self.get_forces(elem, coords)
        return res["energy"], -np.asarray(res["forces"], dtype=float).reshape(-1)

    def __getattr__(self, name: str):
        return getattr(self.base, name)
