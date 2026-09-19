"""Harmonic restraint calculator wrappers (position and internal coordinates).

Two pysisyphus-style Calculator wrappers consumed by multiple workflow stages:

- ``HarmonicFixAtoms`` — harmonic *position* restraint on a subset of atoms.
  Used by ``path_opt`` and ``path_search`` (incl. the DMF path optimizer) to pin
  pre-selected atoms with a quadratic well around their reference coordinates.

- ``HarmonicBiasCalculator`` — harmonic distance, angle or dihedral restraints.
  Used by ``scan`` / ``scan2d`` / ``scan3d`` and by ``opt`` for distance biasing.
  It wraps a base ML/MM calculator and adds the restraint derivatives on one PES.

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
from mlmm.domain.scan_coordinates import (
    coordinate_delta,
    coordinate_derivative,
    coordinate_kind,
    coordinate_target,
    coordinate_value,
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


def harmonic_internal_energy_forces_hessian(
    coords_bohr: np.ndarray,
    k_ev: float,
    restraints: Sequence[Tuple],
    *,
    need_hessian: bool = True,
) -> Tuple[float, np.ndarray, Optional[np.ndarray]]:
    """Evaluate distance/angle/dihedral harmonic restraints."""
    coords = np.asarray(coords_bohr, dtype=float).reshape(-1, 3)
    if restraints and all(coordinate_kind(item) == "distance" for item in restraints):
        return harmonic_pair_energy_forces_hessian(
            coords, float(k_ev) * H_EVAA_2_AU, restraints, need_hessian=need_hessian
        )
    if not np.isfinite(coords).all():
        raise ValueError("Harmonic restraint coordinates must be finite.")
    k_hartree = float(k_ev) * EV2AU
    if not np.isfinite(k_hartree) or k_hartree < 0.0:
        raise ValueError("Harmonic restraint force constant must be finite and non-negative.")

    def energy_force(at: np.ndarray) -> Tuple[float, np.ndarray]:
        energy = 0.0
        force = np.zeros(at.size, dtype=float)
        for restraint in restraints:
            kind = coordinate_kind(restraint)
            value = coordinate_value(at, restraint)
            delta_native = coordinate_delta(kind, value, coordinate_target(restraint))
            delta = delta_native if kind == "distance" else np.deg2rad(delta_native)
            derivative = coordinate_derivative(at, restraint)
            energy += 0.5 * k_hartree * delta * delta
            force -= k_hartree * delta * derivative
        return float(energy), force

    energy, force = energy_force(coords)
    hessian = None
    if need_hessian:
        size = coords.size
        hessian = np.zeros((size, size), dtype=float)
        affected = sorted({3 * int(atom) + axis
                           for item in restraints for atom in item[:-1]
                           for axis in range(3)})
        step = 1.0e-4
        flat = coords.reshape(-1)
        for column in affected:
            plus = flat.copy(); plus[column] += step
            minus = flat.copy(); minus[column] -= step
            _, force_plus = energy_force(plus.reshape(-1, 3))
            _, force_minus = energy_force(minus.reshape(-1, 3))
            hessian[:, column] = -(force_plus - force_minus) / (2.0 * step)
        hessian = 0.5 * (hessian + hessian.T)
    return energy, force, hessian


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
    """Add harmonic distance, angle or dihedral restraints to a calculator.

    Atom indices are 0-based. Distances and distance force constants use Å and
    eV/Å²; angular targets and force constants use degrees and eV/rad².
    """

    def __init__(self, base_calc, k: float = 10.0, pairs: Optional[List[Tuple[int, int, float]]] = None):
        self.base = base_calc
        self.k_evAA = float(k)
        self.k_au_bohr2 = self.k_evAA * H_EVAA_2_AU
        self._restraints: List[Tuple] = list(pairs or [])

    @property
    def _pairs(self):
        return self._restraints

    def set_pairs(self, pairs: List[Tuple[int, int, float]]) -> None:
        self.set_restraints(pairs)

    def set_restraints(self, restraints: Sequence[Tuple]) -> None:
        self._restraints = [tuple(item) for item in restraints]

    def _bias_energy_forces_bohr(self, coords_bohr: np.ndarray) -> Tuple[float, np.ndarray]:
        energy, forces, _ = harmonic_internal_energy_forces_hessian(
            coords_bohr,
            self.k_evAA,
            self._restraints,
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
        if not self._restraints:
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
        if not self._restraints:
            return clone_pes_result(base)
        energy, forces, hessian = harmonic_internal_energy_forces_hessian(
            coords_bohr,
            self.k_evAA,
            self._restraints,
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
