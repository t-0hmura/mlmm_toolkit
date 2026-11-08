"""Native MACE-POLAR embedding using graph_longrange's external-source operators."""

from __future__ import annotations

from contextlib import contextmanager
import threading

import numpy as np
import torch
from ase import Atoms


class MACEPolarEmbedding:
    """Return the embedded-minus-vacuum correction on ML and MM coordinates."""

    def __init__(self, calculator, *, hessian_step=1e-3):
        try:
            from graph_longrange.external_source_energy import GTOElectrostaticExternalSourceEnergy
            from graph_longrange.external_source_features import GTOElectrostaticExternalSourceFeatures
        except ImportError as exc:
            raise ImportError(
                "Native MACE-POLAR embedding requires graph-longrange >= 0.4.4: "
                "pip install 'mace-torch[polar]>=0.3.17' 'graph-longrange>=0.4.4'"
            ) from exc

        models = calculator.models
        if len(models) != 1 or type(models[0]).__name__ != "PolarMACE":
            raise ValueError("Native embedding requires a single PolarMACE model")
        if getattr(calculator, "use_compile", False):
            raise ValueError("Native MACE-POLAR embedding requires an uncompiled calculator")
        self.calculator, self.model = calculator, models[0]
        self.hessian_step = float(hessian_step)
        if not np.isfinite(self.hessian_step) or self.hessian_step <= 0:
            raise ValueError("Embedding Hessian step must be positive and finite")
        self._lock = threading.RLock()
        self._energy = GTOElectrostaticExternalSourceEnergy.from_energy(self.model.coulomb_energy)
        self._features = GTOElectrostaticExternalSourceFeatures.from_features(
            self.model.electric_potential_descriptor, external_scale=0.5,
        )

    @contextmanager
    def _connected(self):
        energy, features = self.model.coulomb_energy, self.model.electric_potential_descriptor
        try:
            self.model.coulomb_energy = self._energy
            self.model.electric_potential_descriptor = self._features
            yield
        finally:
            self._energy.clear_external_sources()
            self._features.clear_external_sources()
            self.model.coulomb_energy = energy
            self.model.electric_potential_descriptor = features

    def _evaluate(self, symbols, coordinates, mm_coordinates, charges, charge, multiplicity, forces):
        # Translate both regions together; PDB offsets must not enlarge the MACE box.
        origin = 0.5 * (coordinates.min(axis=0) + coordinates.max(axis=0))
        atoms = Atoms(symbols, positions=coordinates - origin)
        atoms.info.update(charge=int(charge), spin=int(multiplicity))
        parameter = next(self.model.parameters())
        pos = torch.tensor(coordinates - origin, dtype=parameter.dtype, device=parameter.device,
                           requires_grad=forces)
        mm_pos = torch.tensor(mm_coordinates - origin, dtype=parameter.dtype, device=parameter.device,
                              requires_grad=forces)
        env = pos.new_zeros((len(charges), (self.model.atomic_multipoles_max_l + 1) ** 2))
        env[:, 0] = torch.as_tensor(charges, dtype=pos.dtype, device=pos.device)
        batch = torch.zeros(len(charges), dtype=torch.long, device=pos.device)
        graph = self.calculator._atoms_to_batch(atoms).to_dict()
        graph["positions"] = pos
        with self._lock, torch.enable_grad(), self._connected():
            for block in (self._energy, self._features):
                block.set_external_sources(env, mm_pos, batch)
            embedded = self.model(dict(graph), compute_force=False, training=False,
                                  compute_stress=False, compute_virials=False)["energy"].sum()
            for block in (self._energy, self._features):
                block.clear_external_sources()
            vacuum = self.model(dict(graph), compute_force=False, training=False,
                                compute_stress=False, compute_virials=False)["energy"].sum()
            delta = embedded - vacuum
            derivative = None
            if forces:
                ml_grad, mm_grad = torch.autograd.grad(delta, (pos, mm_pos))
                derivative = -torch.cat((ml_grad, mm_grad)).detach().cpu().numpy()
            result = float(delta.detach())
        if not np.isfinite(result) or (derivative is not None and not np.isfinite(derivative).all()):
            raise ValueError("MACE-POLAR embedding returned non-finite energy or forces")
        return result, derivative

    def compute_correction(self, symbols, coords_ml_ang, mm_coords_ang, mm_charges,
                           charge, multiplicity, *, need_forces=False, need_hessian=False,
                           hessian_mm_rows=None):
        ml = np.asarray(coords_ml_ang, dtype=float).reshape(len(symbols), 3)
        mm = np.asarray(mm_coords_ang, dtype=float).reshape(-1, 3)
        q = np.asarray(mm_charges, dtype=float).reshape(-1)
        if len(q) != len(mm) or not all(np.isfinite(x).all() for x in (ml, mm, q)):
            raise ValueError("Embedding coordinates and charges must be finite with matching lengths")
        rows = list(range(len(mm))) if hessian_mm_rows is None else list(hessian_mm_rows)
        if any(not isinstance(i, (int, np.integer)) or i < 0 or i >= len(mm) for i in rows):
            raise ValueError("Embedding Hessian MM rows must be valid atom indices")
        if len(set(rows)) != len(rows):
            raise ValueError("Embedding Hessian MM rows must be unique")
        active = list(range(len(ml))) + [len(ml) + i for i in rows]
        n = len(ml)
        if not len(q):
            return (0., np.zeros((n, 3)) if need_forces or need_hessian else None,
                    np.zeros((3*n, 3*n)) if need_hessian else None)
        e, f = self._evaluate(symbols, ml, mm, q, charge, multiplicity, need_forces or need_hessian)
        h = None
        if need_hessian:
            full = np.concatenate((ml, mm))
            h = np.empty((3*len(active), 3*len(active)))
            for column, component in enumerate(3*i+j for i in active for j in range(3)):
                plus, minus = full.copy(), full.copy()
                plus.flat[component] += self.hessian_step
                minus.flat[component] -= self.hessian_step
                _, fp = self._evaluate(symbols, plus[:n], plus[n:], q, charge, multiplicity, True)
                _, fm = self._evaluate(symbols, minus[:n], minus[n:], q, charge, multiplicity, True)
                h[:, column] = -(fp[active] - fm[active]).reshape(-1) / (2*self.hessian_step)
            h = 0.5*(h+h.T)
        return e, f, h
