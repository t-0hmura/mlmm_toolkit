"""Stateful PySCF/GPU4PySCF high-level backend for ML/MM workflows."""

from __future__ import annotations

import hashlib
import json
import os
import uuid
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import numpy as np
import torch

from pysisyphus.constants import ANG2BOHR, AU2EV

from mlmm.core.dft_settings import DFTSettings, resolve_dft_settings


def _to_numpy(value):
    getter = getattr(value, "get", None)
    if callable(getter):
        value = getter()
    return None if value is None else np.asarray(value)


def _safe_value(value: Any) -> bool:
    if value is None or isinstance(value, (str, int, float, bool)):
        return True
    if isinstance(value, (list, tuple)):
        return all(_safe_value(item) for item in value)
    if isinstance(value, Mapping):
        return all(isinstance(key, str) and _safe_value(item) for key, item in value.items())
    return False


def _apply_attributes(obj: Any, values: Mapping[str, Any], path: str) -> None:
    for name, value in values.items():
        if not _safe_value(value):
            raise ValueError(f"{path}.{name} is not a YAML-safe value.")
        if not hasattr(obj, name) or callable(getattr(obj, name)):
            raise ValueError(f"Unknown or callable PySCF attribute {path}.{name}.")
        setattr(obj, name, value)


def _write_checkpoint_generation(path: Path, generation: str) -> None:
    import h5py

    with h5py.File(path, "a") as handle:
        handle.attrs["mlmm_record_generation"] = generation


def _read_checkpoint_generation(path: Path) -> str:
    import h5py

    with h5py.File(path, "r") as handle:
        value = handle.attrs["mlmm_record_generation"]
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    return str(value)


def _validated_checkpoint_record(
    loaded_mol: Any,
    record: Mapping[str, Any],
    symbols: Sequence[str],
    coords_ang: np.ndarray,
) -> Dict[str, Any]:
    energy = float(record["e_tot"])
    mo_coeff = _to_numpy(record["mo_coeff"])
    mo_occ = _to_numpy(record["mo_occ"])
    mo_energy = _to_numpy(record["mo_energy"])
    if not np.isfinite(energy) or any(
        value is None or not np.all(np.isfinite(value))
        for value in (mo_coeff, mo_occ, mo_energy)
    ):
        raise ValueError("checkpoint contains non-finite SCF data")
    if mo_coeff.ndim not in (2, 3) or mo_occ.shape != mo_energy.shape:
        raise ValueError("checkpoint orbital arrays have inconsistent ranks")
    if mo_coeff.shape[-1] != mo_occ.shape[-1]:
        raise ValueError("checkpoint orbital dimensions are inconsistent")
    if mo_coeff.ndim == 3 and mo_coeff.shape[0] != mo_occ.shape[0]:
        raise ValueError("checkpoint spin dimensions are inconsistent")
    if int(loaded_mol.nao_nr()) != int(mo_coeff.shape[-2]):
        raise ValueError("checkpoint AO dimension is inconsistent")
    loaded_symbols = tuple(
        str(loaded_mol.atom_symbol(index)) for index in range(int(loaded_mol.natm))
    )
    if loaded_symbols != tuple(map(str, symbols)):
        raise ValueError("checkpoint Mole symbols do not match metadata")
    loaded_coords = np.asarray(
        loaded_mol.atom_coords(unit="Angstrom"), dtype=float
    ).reshape(-1, 3)
    if loaded_coords.shape != coords_ang.shape or not np.allclose(
        loaded_coords, coords_ang, atol=1.0e-8, rtol=0.0
    ):
        raise ValueError("checkpoint Mole coordinates do not match metadata")
    return {
        "energy": energy,
        "mol": loaded_mol,
        "mo_coeff": mo_coeff.copy(),
        "mo_occ": mo_occ.copy(),
        "mo_energy": mo_energy.copy(),
    }


def _subtract_packed_mm_hcore(h1e, mm_hcore, pack_tril):
    """Subtract a square GPU MM potential from a packed host Hcore."""
    mm_hcore_packed = pack_tril(mm_hcore)
    to_host = getattr(mm_hcore_packed, "get", None)
    if callable(to_host):
        mm_hcore_packed = to_host()
    mm_hcore_packed = np.asarray(mm_hcore_packed)
    if h1e.shape != mm_hcore_packed.shape:
        raise RuntimeError(
            "GPU low-memory QMMM produced inconsistent packed Hcore shapes: "
            f"{h1e.shape} != {mm_hcore_packed.shape}."
        )
    h1e -= mm_hcore_packed
    return h1e


def _attach_gpu_lowmem_point_charges(
    mf,
    coords: np.ndarray,
    charges: np.ndarray,
    *,
    unit: str,
):
    """Attach GPU point charges while retaining rks_lowmem's packed Hcore."""
    from pyscf import lib as pyscf_lib
    from pyscf.qmmm import mm_mole
    from gpu4pyscf.gto.int3c1e import int1e_grids
    from gpu4pyscf.lib.cupy_helper import pack_tril
    from gpu4pyscf.qmmm import itrf as gpu_qmmm_itrf

    class _PackedHcoreQMMMSCF(gpu_qmmm_itrf.QMMMSCF):
        def get_hcore(self, mol=None):
            if mol is None:
                mol = self.mol
            h1e = super(gpu_qmmm_itrf.QMMMSCF, self).get_hcore(mol)
            if np.ndim(h1e) != 1:
                raise RuntimeError(
                    "GPU low-memory QMMM expected a packed one-dimensional Hcore."
                )
            mm_mol = self.mm_mol
            expnts = (
                mm_mol.get_zetas()
                if mm_mol.charge_model == "gaussian"
                else None
            )
            mm_hcore = int1e_grids(
                mol,
                mm_mol.atom_coords(),
                charges=mm_mol.atom_charges(),
                charge_exponents=expnts,
            )
            return _subtract_packed_mm_hcore(h1e, mm_hcore, pack_tril)

    mm_mol = mm_mole.create_mm_mol(coords, charges, unit=unit)
    mixin = _PackedHcoreQMMMSCF
    return pyscf_lib.set_class(mixin(mf, mm_mol), (mixin, mf.__class__))


class PySCFDFTSession:
    """Own one live SCF method, last-good orbitals, and an exact-geometry cache."""

    def __init__(self, settings: DFTSettings):
        self.settings = settings
        self._scanner = None
        self._using_rks_lowmem = False
        self._last_key: Optional[str] = None
        self._cache: Dict[str, Any] = {}
        self._last_good: Optional[Dict[str, Any]] = None
        self._mm_coords: Optional[np.ndarray] = None
        self._mm_charges: Optional[np.ndarray] = None
        self._last_symbols: Optional[Tuple[str, ...]] = None
        self._last_coords: Optional[np.ndarray] = None
        self._checkpoint_attempted = False
        self._pending_checkpoint = False
        self.checkpoint_status: Optional[Dict[str, Any]] = None
        self.metrics: list[Dict[str, Any]] = []

    @staticmethod
    def _key(symbols, coords, mm_coords, mm_charges) -> str:
        digest = hashlib.sha256("\0".join(map(str, symbols)).encode())
        for value in (coords, mm_coords, mm_charges):
            if value is not None:
                digest.update(np.ascontiguousarray(value, dtype=np.float64).tobytes())
        return digest.hexdigest()

    def set_embedding(
        self, mm_coords_ang: Optional[np.ndarray], mm_charges: Optional[np.ndarray]
    ) -> None:
        self._mm_coords = None if mm_coords_ang is None else np.asarray(mm_coords_ang, dtype=float).reshape(-1, 3)
        self._mm_charges = None if mm_charges is None else np.asarray(mm_charges, dtype=float).reshape(-1)

    def _make_mol(self, symbols: Sequence[str], coords_ang: np.ndarray):
        try:
            from pyscf import gto, lib
        except ImportError as exc:
            raise ImportError("PySCF is required for backend='dft'.") from exc
        lib.num_threads(int(self.settings.nprocs))
        mol = gto.Mole()
        mol.atom = [(str(s), tuple(map(float, xyz))) for s, xyz in zip(symbols, coords_ang)]
        mol.unit = "Angstrom"
        mol.basis = self.settings.basis
        mol.charge = self.settings.charge
        mol.spin = self.settings.multiplicity - 1
        mol.verbose = self.settings.verbose
        if self.settings.memory_mb is not None:
            mol.max_memory = int(self.settings.memory_mb)
        _apply_attributes(mol, self.settings.pyscf.get("mol", {}), "pyscf.mol")
        if (
            not mol.ecp
            and self.settings.basis.casefold().startswith("def2")
            and any(int(gto.charge(str(symbol))) >= 37 for symbol in symbols)
        ):
            mol.ecp = self.settings.basis
        mol.build()
        return mol

    def _build_scanner(self, mol) -> None:
        from pyscf import dft, qmmm, scf

        unrestricted = self.settings.multiplicity != 1
        self._using_rks_lowmem = self.settings.use_rks_lowmem
        if self._using_rks_lowmem:
            try:
                from gpu4pyscf.dft import rks_lowmem
            except ImportError as exc:
                raise RuntimeError(
                    "gpu4pyscf.dft.rks_lowmem is required by the default "
                    "closed-shell GPU low-memory path. Install a compatible "
                    "GPU4PySCF release or explicitly use --no-lowmem."
                ) from exc
            xc = "HF" if self.settings.is_hf else self.settings.functional
            mf = rks_lowmem.RKS(mol, xc=xc)
        elif self.settings.is_hf:
            mf = scf.UHF(mol) if unrestricted else scf.RHF(mol)
        else:
            mf = dft.UKS(mol) if unrestricted else dft.RKS(mol)
            mf.xc = self.settings.functional
        mf.conv_tol = self.settings.conv_tol
        mf.max_cycle = self.settings.max_cycle
        mf.chkfile = None
        _apply_attributes(mf, self.settings.pyscf.get("mf", {}), "pyscf.mf")
        if hasattr(mf, "grids"):
            mf.grids.level = self.settings.grid_level
            _apply_attributes(mf.grids, self.settings.pyscf.get("grids", {}), "pyscf.grids")
        if self.settings.density_fit:
            kwargs = {}
            if self.settings.auxbasis:
                kwargs["auxbasis"] = self.settings.auxbasis
            kwargs.update({
                key: value
                for key, value in self.settings.pyscf.get("density_fit", {}).items()
                if key != "enabled"
            })
            mf = mf.density_fit(**kwargs)
            _apply_attributes(mf.with_df, self.settings.pyscf.get("with_df", {}), "pyscf.with_df")
        if self.settings.engine == "gpu" and not self._using_rks_lowmem:
            try:
                mf = mf.to_gpu()
            except Exception as exc:
                raise RuntimeError(f"Could not initialize GPU4PySCF: {exc}") from exc
        if self.settings.embedcharge:
            if self._mm_coords is None or self._mm_charges is None:
                raise ValueError("DFT electrostatic embedding requires MM coordinates and charges.")
            if self.settings.engine == "gpu":
                if self._using_rks_lowmem:
                    mf = _attach_gpu_lowmem_point_charges(
                        mf, self._mm_coords, self._mm_charges, unit="Angstrom"
                    )
                else:
                    from gpu4pyscf import qmmm as gpu_qmmm

                    mf = gpu_qmmm.mm_charge(
                        mf, self._mm_coords, self._mm_charges, unit="Angstrom"
                    )
            else:
                mf = qmmm.mm_charge(mf, self._mm_coords, self._mm_charges, unit="Angstrom")
        mf.chkfile = None
        self._scanner = mf if self._using_rks_lowmem else mf.as_scanner()
        if self._last_good is not None:
            def engine_array(value):
                if self.settings.engine == "gpu":
                    import cupy

                    return cupy.asarray(value)
                return np.asarray(value).copy()

            self._scanner.mo_coeff = engine_array(self._last_good["mo_coeff"])
            self._scanner.mo_occ = engine_array(self._last_good["mo_occ"])
            self._scanner.mo_energy = engine_array(self._last_good["mo_energy"])
            self._scanner._last_mol_fp = mol.ao_loc.copy()

    def _update_mm_mol(self) -> None:
        if not self.settings.embedcharge or self._scanner is None:
            return
        if self.settings.engine == "gpu":
            from gpu4pyscf.qmmm import mm_mole
        else:
            from pyscf.qmmm import mm_mole
        self._scanner.mm_mol = mm_mole.create_mm_mol(
            self._mm_coords, self._mm_charges, unit="Angstrom"
        )

    def _run_scf(self, mol) -> Tuple[float, str]:
        had_scanner = self._scanner is not None
        guess = (
            "checkpoint"
            if self._pending_checkpoint
            else "previous_density"
            if had_scanner
            else "fresh"
        )
        if self._scanner is None:
            self._build_scanner(mol)
            self._update_mm_mol()
            if self._using_rks_lowmem:
                dm0 = (
                    self._scanner.make_rdm1()
                    if self._last_good is not None
                    else None
                )
                energy = float(self._scanner.kernel(dm0=dm0))
            else:
                energy = float(self._scanner(mol))
        elif self._using_rks_lowmem:
            dm0 = self._scanner.make_rdm1()
            self._scanner = None
            self._build_scanner(mol)
            self._update_mm_mol()
            energy = float(self._scanner.kernel(dm0=dm0))
        else:
            self._update_mm_mol()
            energy = float(self._scanner(mol))
        if not bool(getattr(self._scanner, "converged", False)):
            self._scanner = None
            saved = self._last_good
            self._last_good = None
            try:
                self._build_scanner(mol)
            finally:
                self._last_good = saved
            self._update_mm_mol()
            energy = float(
                self._scanner.kernel()
                if self._using_rks_lowmem
                else self._scanner(mol)
            )
            guess = "fresh_retry"
            if not bool(getattr(self._scanner, "converged", False)):
                self._scanner = None
                advice = (
                    " If sufficient GPU and host memory are available, retry with "
                    "--no-lowmem; standard density-fitted SCF may converge more robustly."
                    if self.settings.lowmem
                    else ""
                )
                raise RuntimeError(
                    "PySCF SCF did not converge from either reused or fresh density."
                    + advice
                )
        self._pending_checkpoint = False
        return energy, guess

    def _metadata(
        self,
        symbols,
        coords,
        *,
        generation: str,
        mm_coords: Optional[np.ndarray] = None,
        mm_charges: Optional[np.ndarray] = None,
    ) -> Dict[str, Any]:
        identity = json.dumps(
            self.settings.scientific_identity(), sort_keys=True, separators=(",", ":")
        )
        return {
            "schema": 2,
            "generation": generation,
            "symbols": list(map(str, symbols)),
            "coordinates_angstrom": np.asarray(coords, dtype=float).tolist(),
            "mm_coordinates_angstrom": (
                None if mm_coords is None else np.asarray(mm_coords).tolist()
            ),
            "mm_charges": (
                None if mm_charges is None else np.asarray(mm_charges).tolist()
            ),
            "scientific_identity_sha256": hashlib.sha256(identity.encode()).hexdigest(),
        }

    def load_scf_checkpoint(self, path, symbols, coords) -> bool:
        destination = Path(path)
        metadata_path = destination.with_suffix(destination.suffix + ".json")
        self._checkpoint_attempted = True
        if not destination.is_file() or not metadata_path.is_file():
            self.checkpoint_status = {
                "loaded": False,
                "reason": "checkpoint_or_metadata_missing",
                "path": str(destination),
            }
            return False
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            if metadata.get("schema") != 2:
                self.checkpoint_status = {
                    "loaded": False,
                    "reason": "schema_mismatch",
                    "path": str(destination),
                }
                return False
            expected = self._metadata(
                symbols,
                coords,
                generation=str(metadata.get("generation", "")),
                mm_coords=self._mm_coords,
                mm_charges=self._mm_charges,
            )
            generation = str(metadata.get("generation", ""))
            if not generation or _read_checkpoint_generation(destination) != generation:
                self.checkpoint_status = {
                    "loaded": False,
                    "reason": "generation_mismatch",
                    "path": str(destination),
                }
                return False
            if metadata.get("symbols") != expected["symbols"]:
                self.checkpoint_status = {
                    "loaded": False, "reason": "symbols_mismatch", "path": str(destination)
                }
                return False
            if metadata.get("scientific_identity_sha256") != expected["scientific_identity_sha256"]:
                self.checkpoint_status = {
                    "loaded": False, "reason": "settings_mismatch", "path": str(destination)
                }
                return False
            for key in (
                "coordinates_angstrom", "mm_coordinates_angstrom", "mm_charges"
            ):
                saved_value = metadata.get(key)
                expected_value = expected.get(key)
                if saved_value is None or expected_value is None:
                    if saved_value is not expected_value:
                        self.checkpoint_status = {
                            "loaded": False, "reason": f"{key}_mismatch", "path": str(destination)
                        }
                        return False
                    continue
                saved = np.asarray(saved_value, dtype=float)
                target = np.asarray(expected_value, dtype=float)
                if saved.shape != target.shape or not np.allclose(saved, target, atol=1.0e-8, rtol=0.0):
                    self.checkpoint_status = {
                        "loaded": False, "reason": f"{key}_mismatch", "path": str(destination)
                    }
                    return False
            from pyscf.scf import chkfile

            loaded_mol, record = chkfile.load_scf(str(destination))
            if int(loaded_mol.charge) != int(self.settings.charge) or int(
                loaded_mol.spin
            ) != int(self.settings.multiplicity) - 1:
                raise ValueError("checkpoint Mole charge/spin mismatch")
            self._last_good = _validated_checkpoint_record(
                loaded_mol, record, tuple(map(str, symbols)), np.asarray(coords, dtype=float)
            )
            self._last_good.update(
                {
                    "symbols": tuple(map(str, symbols)),
                    "coordinates_angstrom": np.asarray(coords, dtype=float).copy(),
                    "mm_coordinates_angstrom": (
                        None if self._mm_coords is None else self._mm_coords.copy()
                    ),
                    "mm_charges": (
                        None if self._mm_charges is None else self._mm_charges.copy()
                    ),
                }
            )
            self._scanner = None
            self._pending_checkpoint = True
            self._last_key = None
            self._cache = {}
            self.checkpoint_status = {
                "loaded": True, "reason": "matched", "path": str(destination)
            }
            return True
        except (OSError, ValueError, KeyError, TypeError):
            self.checkpoint_status = {
                "loaded": False, "reason": "invalid_checkpoint", "path": str(destination)
            }
            return False

    def save_last_checkpoint(self, path=None) -> Optional[Path]:
        if not self.settings.save_scf_checkpoint:
            raise RuntimeError("SCF checkpoint saving is disabled.")
        destination_value = path or self.settings.checkpoint_path
        if destination_value is None or self._last_symbols is None or self._last_coords is None:
            return None
        from pyscf.scf import chkfile

        destination = Path(destination_value)
        destination.parent.mkdir(parents=True, exist_ok=True)
        generation = uuid.uuid4().hex
        tmp = destination.with_name(f"{destination.name}.{generation}.tmp")
        metadata_path = destination.with_suffix(destination.suffix + ".json")
        metadata_tmp = metadata_path.with_name(
            f"{metadata_path.name}.{generation}.tmp"
        )
        record = self._last_good
        if record is None:
            self._commit_last_good()
            record = self._last_good
        try:
            chkfile.dump_scf(
                record["mol"],
                str(tmp),
                float(record["energy"]),
                record["mo_energy"],
                record["mo_coeff"],
                record["mo_occ"],
            )
            _write_checkpoint_generation(tmp, generation)
            metadata_tmp.write_text(
                json.dumps(
                    self._metadata(
                        record["symbols"],
                        record["coordinates_angstrom"],
                        generation=generation,
                        mm_coords=record["mm_coordinates_angstrom"],
                        mm_charges=record["mm_charges"],
                    ),
                    indent=2,
                    sort_keys=True,
                ) + "\n",
                encoding="utf-8",
            )
            os.replace(tmp, destination)
            os.replace(metadata_tmp, metadata_path)
        finally:
            tmp.unlink(missing_ok=True)
            metadata_tmp.unlink(missing_ok=True)
        return destination

    def evaluate(self, atoms, *, need_forces: bool = True) -> Dict[str, Any]:
        symbols = tuple(atoms.get_chemical_symbols())
        coords = np.asarray(atoms.get_positions(), dtype=float)
        if not self._checkpoint_attempted:
            self._checkpoint_attempted = True
            if self.settings.checkpoint_path:
                self.load_scf_checkpoint(
                    self.settings.checkpoint_path, symbols, coords
                )
        key = self._key(symbols, coords, self._mm_coords, self._mm_charges)
        if key == self._last_key and (not need_forces or self._cache.get("forces_au") is not None):
            result = dict(self._cache)
            result["cache_hit"] = True
            return result
        if key != self._last_key:
            mol = self._make_mol(symbols, coords)
            energy, guess = self._run_scf(mol)
            self._last_key = key
            self._last_symbols = symbols
            self._last_coords = coords.copy()
            self._cache = {
                "energy_au": energy, "forces_au": None, "mm_forces_au": None,
                "guess_source": guess, "cache_hit": False,
            }
            if self.settings.save_scf_checkpoint:
                self._commit_last_good()
            else:
                # The live SCF method supplies the next-geometry guess; do not
                # duplicate its full MO arrays on host when persistence is disabled.
                self._last_good = None
            self.metrics.append({
                "guess_source": guess,
                "cycles": int(getattr(self._scanner, "cycles", -1)),
                "converged": True,
            })
        if need_forces and self._cache.get("forces_au") is None:
            grad = self._scanner.nuc_grad_method()
            if hasattr(grad, "auxbasis_response") and self.settings.density_fit:
                grad.auxbasis_response = True
            self._cache["forces_au"] = -_to_numpy(grad.kernel()).reshape(-1, 3)
            if self.settings.embedcharge:
                dm = self._scanner.make_rdm1()
                if getattr(dm, "ndim", 0) == 3:
                    dm = dm.sum(axis=0)
                mm_gradient = _to_numpy(grad.grad_hcore_mm(dm)) + _to_numpy(grad.grad_nuc_mm())
                self._cache["mm_forces_au"] = -mm_gradient.reshape(-1, 3)
        return dict(self._cache)

    def _commit_last_good(self) -> None:
        self._last_good = {
            "energy": float(self._scanner.e_tot),
            "mol": self._scanner.mol,
            "mo_coeff": _to_numpy(self._scanner.mo_coeff).copy(),
            "mo_occ": _to_numpy(self._scanner.mo_occ).copy(),
            "mo_energy": _to_numpy(self._scanner.mo_energy).copy(),
            "symbols": self._last_symbols,
            "coordinates_angstrom": (
                None if self._last_coords is None else self._last_coords.copy()
            ),
            "mm_coordinates_angstrom": (
                None if self._mm_coords is None else self._mm_coords.copy()
            ),
            "mm_charges": (
                None if self._mm_charges is None else self._mm_charges.copy()
            ),
        }

    def close(self) -> None:
        """Release all live electronic and embedding state owned by this session."""

        self._scanner = None
        self._last_good = None
        self._cache.clear()
        self._last_key = None
        self._last_symbols = None
        self._last_coords = None
        self._mm_coords = None
        self._mm_charges = None
        self._pending_checkpoint = False

    def hessian(self, atoms) -> np.ndarray:
        if self.settings.embedcharge:
            raise RuntimeError("Embedded DFT uses the complete ML/MM full-force finite-difference Hessian.")
        self.evaluate(atoms, need_forces=True)
        hess = self._scanner.Hessian()
        if hasattr(hess, "auxbasis_response") and self.settings.density_fit:
            hess.auxbasis_response = 2
        raw = _to_numpy(hess.kernel())
        n_atoms = len(atoms)
        return raw.transpose(0, 2, 1, 3).reshape(3 * n_atoms, 3 * n_atoms)


class DFTBackend:
    """Adapter implementing the existing MLMM high-level backend contract."""

    name = "dft"

    def __init__(self, settings: DFTSettings):
        self.settings = settings
        self.session = PySCFDFTSession(settings)
        self._device = torch.device("cuda" if settings.engine == "gpu" else "cpu")
        self._mm_indices: list[int] = []
        self._last_mm_forces: Optional[np.ndarray] = None
        self._closed = False

    @property
    def device(self) -> torch.device:
        return self._device

    @property
    def supports_analytical_hessian(self) -> bool:
        return not self.settings.embedcharge

    def set_embedding(self, coords, charges, indices) -> None:
        self._mm_indices = list(map(int, indices))
        self.session.set_embedding(coords, charges)

    def eval(self, atoms, need_grad: bool = True):
        result = self.session.evaluate(atoms, need_forces=need_grad)
        forces = (
            np.zeros((len(atoms), 3), dtype=float)
            if result["forces_au"] is None
            else np.asarray(result["forces_au"]) * AU2EV * ANG2BOHR
        )
        self._last_mm_forces = (
            None if result["mm_forces_au"] is None
            else np.asarray(result["mm_forces_au"]) * AU2EV * ANG2BOHR
        )
        return float(result["energy_au"]) * AU2EV, forces, atoms.copy()

    def energy(self, atoms) -> float:
        return float(self.session.evaluate(atoms, need_forces=False)["energy_au"]) * AU2EV

    def mm_forces(self) -> Tuple[list[int], Optional[np.ndarray]]:
        return list(self._mm_indices), self._last_mm_forces

    def save_scf_checkpoint(self, path=None, atoms=None):
        if atoms is not None:
            self.session.evaluate(atoms, need_forces=False)
        return self.session.save_last_checkpoint(path)

    def load_scf_checkpoint(self, path, atoms=None) -> bool:
        if atoms is not None:
            return self.session.load_scf_checkpoint(
                path,
                tuple(atoms.get_chemical_symbols()),
                np.asarray(atoms.get_positions(), dtype=float),
            )
        self.session.settings = DFTSettings(
            **{
                **self.session.settings.__dict__,
                "checkpoint_path": str(path),
            }
        )
        self.session._checkpoint_attempted = False
        return Path(path).is_file()

    def close(self) -> None:
        if self._closed:
            return
        try:
            if self.settings.save_scf_checkpoint:
                self.session.save_last_checkpoint()
        finally:
            self.session.close()
            self._last_mm_forces = None
            self._mm_indices.clear()
            self._closed = True

    def hessian_analytical(self, opaque, n_atoms: int, *, dtype: torch.dtype):
        hessian = self.session.hessian(opaque) * AU2EV * ANG2BOHR * ANG2BOHR
        return torch.as_tensor(hessian, dtype=dtype, device=self._device).reshape(n_atoms, 3, n_atoms, 3)

    def hessian_fd(self, atoms, freeze_model, *, eps_ang=1.0e-3, dtype=torch.float32, device=torch.device("cpu")):
        n_atoms = len(atoms)
        H = torch.zeros((3 * n_atoms, 3 * n_atoms), dtype=dtype, device=device)
        coord0 = atoms.get_positions().copy()
        frozen = set(map(int, freeze_model))
        for k in range(3 * n_atoms):
            if k // 3 in frozen:
                continue
            plus = atoms.copy()
            minus = atoms.copy()
            plus.positions = coord0.copy()
            minus.positions = coord0.copy()
            plus.positions.reshape(-1)[k] += eps_ang
            minus.positions.reshape(-1)[k] -= eps_ang
            _, fp, _ = self.eval(plus, need_grad=True)
            _, fm, _ = self.eval(minus, need_grad=True)
            H[:, k] = torch.as_tensor(-(fp - fm).reshape(-1) / (2 * eps_ang), dtype=dtype, device=device)
        self.eval(atoms, need_grad=True)
        return H.reshape(n_atoms, 3, n_atoms, 3)


def create_dft_backend(
    values: Mapping[str, Any], *, model_charge: Optional[int] = None,
    model_mult: Optional[int] = None,
) -> DFTBackend:
    charge = values.get("charge", 0) if model_charge is None else model_charge
    multiplicity = (
        values.get("multiplicity", 1) if model_mult is None else model_mult
    )
    settings = resolve_dft_settings({
        "backend": "dft",
        "model_charge": int(charge),
        "model_mult": int(multiplicity),
        "embedcharge": values.get("embedcharge", False),
        "embedcharge_cutoff": values.get("embedcharge_cutoff", 12.0),
        "dft_settings": dict(values),
    })
    return DFTBackend(settings)


__all__ = ["DFTBackend", "PySCFDFTSession", "create_dft_backend"]
