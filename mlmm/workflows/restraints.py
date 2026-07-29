"""Harmonic restraint calculator wrappers (position fix + bond-length bias).

Two pysisyphus-style Calculator wrappers consumed by multiple workflow stages:

- ``HarmonicFixAtoms`` — harmonic *position* restraint on a subset of atoms.
  Used by ``path_opt`` and ``path_search`` (incl. the DMF path optimizer) to pin
  pre-selected atoms with a quadratic well around their reference coordinates.

- ``HarmonicBiasCalculator`` — harmonic *distance* restraint on selected atom
  pairs (1-based / 0-based indexing per caller). Used by ``scan`` / ``scan2d`` /
  ``scan3d`` for bond-length staged scans, and by ``opt`` for ad-hoc distance
  biasing. Wraps a base ML/MM calculator and adds the bias E/F/H on one PES.

Both classes are pure-Python (numpy only) and do not import any MLIP SDK, so they
belong with the workflow orchestration layer rather than ``io`` or ``backends``.
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np
from ase.calculators.calculator import Calculator
from pysisyphus.constants import ANG2BOHR, AU2EV

from mlmm.core.pes_composition import (
    clone_pes_result,
    compose_additive_pes_result,
)


# eV/Å² → Hartree/Bohr² conversion (= k_evAA * H_EVAA_2_AU)
EV2AU = 1.0 / AU2EV
H_EVAA_2_AU = EV2AU / ANG2BOHR / ANG2BOHR


def _calculator_chain(root):
    """Yield a calculator and its ``base``/``core`` delegates once each."""

    pending = [root]
    seen = set()
    while pending:
        calculator = pending.pop(0)
        if calculator is None or id(calculator) in seen:
            continue
        seen.add(id(calculator))
        yield calculator
        for attribute in ("base", "core"):
            try:
                delegate = getattr(calculator, attribute)
            except AttributeError:
                continue
            if delegate is not calculator:
                pending.append(delegate)


def harmonic_pair_energy_forces_hessian(
    coords_bohr: np.ndarray,
    k_au_bohr2: float,
    pairs: Sequence[Tuple[int, int, float]],
    *,
    need_hessian: bool = True,
) -> Tuple[float, np.ndarray, Optional[np.ndarray]]:
    """Evaluate harmonic pair energy, force, and exact Cartesian Hessian.

    Coordinates and returned quantities use atomic units; pair targets remain
    in Angstrom to match the public restraint configuration.
    """

    coords = np.asarray(coords_bohr, dtype=float).reshape(-1, 3)
    n_atoms = coords.shape[0]
    force = np.zeros((n_atoms, 3), dtype=float)
    hessian = (
        np.zeros((3 * n_atoms, 3 * n_atoms), dtype=float)
        if need_hessian
        else None
    )
    energy = 0.0
    identity = np.eye(3, dtype=float)
    k = float(k_au_bohr2)
    if not np.isfinite(coords).all():
        raise ValueError("Harmonic restraint coordinates must be finite.")
    if not np.isfinite(k) or k < 0.0:
        raise ValueError(
            "Harmonic restraint force constant must be finite and non-negative."
        )

    for pair_index, (i_raw, j_raw, target_ang) in enumerate(pairs, start=1):
        i, j = int(i_raw), int(j_raw)
        if not (0 <= i < n_atoms and 0 <= j < n_atoms):
            raise ValueError(
                f"Harmonic restraint pair {pair_index} uses atom index "
                f"({i}, {j}) outside the valid range 0..{n_atoms - 1}."
            )
        if i == j:
            raise ValueError(
                f"Harmonic restraint pair {pair_index} must use two distinct atoms."
            )
        target = float(target_ang)
        if not np.isfinite(target) or target <= 0.0:
            raise ValueError(
                f"Harmonic restraint pair {pair_index} target must be finite "
                "and greater than zero."
            )
        delta = coords[i] - coords[j]
        distance = float(np.linalg.norm(delta))
        if distance < 1.0e-14:
            raise ValueError(
                f"Harmonic restraint pair {pair_index} has coincident atoms; "
                "its direction is undefined."
            )
        target_bohr = target * ANG2BOHR
        displacement = distance - target_bohr
        unit = delta / distance

        energy += 0.5 * k * displacement * displacement
        pair_force = -k * displacement * unit
        force[i] += pair_force
        force[j] -= pair_force

        if hessian is not None:
            outer = np.outer(unit, unit)
            block = k * (
                outer + (displacement / distance) * (identity - outer)
            )
            i_slice = slice(3 * i, 3 * i + 3)
            j_slice = slice(3 * j, 3 * j + 3)
            hessian[i_slice, i_slice] += block
            hessian[j_slice, j_slice] += block
            hessian[i_slice, j_slice] -= block
            hessian[j_slice, i_slice] -= block

    return float(energy), force.reshape(-1), hessian


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
        try:
            k_fix = float(k_fix)
        except (TypeError, ValueError) as exc:
            raise ValueError("k_fix must be a finite positive number.") from exc
        if not np.isfinite(k_fix) or k_fix <= 0.0:
            raise ValueError("k_fix must be a finite positive number.")
        self.indices = idx
        self.ref_positions = ref_pos
        self.k_fix = k_fix

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
    Forces and exact Cartesian Hessians are added to the base calculator's
    output. Indices are 0-based Cartesian atom indices.

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
        energy, forces, _ = harmonic_pair_energy_forces_hessian(
            coords_bohr,
            self.k_au_bohr2,
            self._pairs,
            need_hessian=False,
        )
        return energy, forces

    @property
    def _constrained_atoms(self) -> Sequence[int]:
        atoms = getattr(self.base, "freeze_atoms", ())
        return tuple(() if atoms is None else atoms)

    @property
    def _hessian_constrained_atoms(self) -> Sequence[int]:
        """Resolve ML/MM's stricter Hessian mask through wrapper delegates.

        ``MLMMCore.compute`` refreshes ``effective_hess_freeze_atoms`` at its
        evaluation boundary, so callers intentionally read this property only
        after the base Hessian evaluation.  The fallback union also supports
        calculator implementations that expose the two source sets but not
        the derived attribute.
        """

        calculators = tuple(_calculator_chain(self.base))
        effective_atoms = set()
        found_effective = False
        for calculator in calculators:
            try:
                effective = getattr(calculator, "effective_hess_freeze_atoms")
            except AttributeError:
                continue
            if effective is not None:
                found_effective = True
                effective_atoms.update(int(atom) for atom in effective)
        if found_effective:
            return tuple(sorted(effective_atoms))

        constrained = set()
        for calculator in calculators:
            for attribute in ("freeze_atoms", "hess_freeze_atoms"):
                try:
                    atoms = getattr(calculator, attribute)
                except AttributeError:
                    continue
                if atoms is not None:
                    constrained.update(int(atom) for atom in atoms)
        return tuple(sorted(constrained))

    def get_forces(self, elem, coords):
        coords_bohr = np.asarray(coords, dtype=float).reshape(-1, 3)
        base = self.base.get_forces(elem, coords_bohr)
        Ebias, Fbias = self._bias_energy_forces_bohr(coords_bohr)
        return compose_additive_pes_result(
            base,
            n_atoms=coords_bohr.shape[0],
            energy_delta=Ebias,
            force_delta_full=Fbias,
            constrained_atoms=self._constrained_atoms,
        )

    def get_energy(self, elem, coords):
        coords_bohr = np.asarray(coords, dtype=float).reshape(-1, 3)
        base = self.base.get_energy(elem, coords_bohr)
        if not self._pairs:
            return clone_pes_result(base)
        Ebias, _ = self._bias_energy_forces_bohr(coords_bohr)
        return compose_additive_pes_result(
            base,
            n_atoms=coords_bohr.shape[0],
            energy_delta=Ebias,
        )

    def get_hessian(self, elem, coords):
        coords_bohr = np.asarray(coords, dtype=float).reshape(-1, 3)
        base = self.base.get_hessian(elem, coords_bohr)
        if not self._pairs:
            return clone_pes_result(base)
        energy, forces, hessian = harmonic_pair_energy_forces_hessian(
            coords_bohr,
            self.k_au_bohr2,
            self._pairs,
            need_hessian=True,
        )
        assert hessian is not None
        return compose_additive_pes_result(
            base,
            n_atoms=coords_bohr.shape[0],
            energy_delta=energy,
            force_delta_full=forces,
            hessian_delta_full=hessian,
            constrained_atoms=self._constrained_atoms,
            hessian_constrained_atoms=self._hessian_constrained_atoms,
        )

    def get_energy_and_forces(self, elem, coords):
        res = self.get_forces(elem, coords)
        return res["energy"], res["forces"]

    def get_energy_and_gradient(self, elem, coords):
        res = self.get_forces(elem, coords)
        return res["energy"], -np.asarray(res["forces"], dtype=float).reshape(-1)

    def __getattr__(self, name: str):
        return getattr(self.base, name)
