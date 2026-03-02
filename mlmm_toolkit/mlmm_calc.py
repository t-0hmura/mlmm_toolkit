"""
ONIOM-like ML/MM calculator coupling FAIR-Chem UMA (ML) and hessian_ff (MM).

Example:
    calc = mlmm(input_pdb="input.pdb", real_parm7="real.parm7", model_pdb="model.pdb", charge=0)

For detailed documentation, see: docs/mlmm_calc.md
"""

from __future__ import annotations

import os
import shutil
import tempfile
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple
from concurrent.futures import ThreadPoolExecutor

import click
import numpy as np
import torch
import torch.nn as nn

from ase import Atoms
from ase.io import read
from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixAtoms

import parmed as pmd
from hessian_ff import ForceFieldTorch, load_coords, load_system
from hessian_ff.analytical_hessian import build_analytical_hessian

# Optional OpenMM import
try:
    import openmm as mm
    from openmm import app, unit, Platform
    from openmm.unit import ScaledUnit, joule
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False

try:
    from fairchem.core import pretrained_mlip
    from fairchem.core.datasets.atomic_data import AtomicData
    from fairchem.core.datasets import data_list_collater
except ImportError as e:
    raise ImportError(
        "fairchem-core is required. Install with `pip install fairchem-core` "
        "and ensure Hugging Face authentication is configured."
    ) from e

# ---------- PySisyphus unit constants ----------
from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AU2EV, AU2KCALPERMOL
EV2AU = 1.0 / AU2EV  # eV → Hartree
KCALMOL2EV = AU2EV / AU2KCALPERMOL  # kcal/mol -> eV


# ======================================================================
#                           Utilities
# ======================================================================

def _fixed_indices_from_constraints(atoms: Atoms) -> set[int]:
    fixed: set[int] = set()
    for c in atoms.constraints or []:
        if isinstance(c, FixAtoms):
            fixed.update(int(i) for i in c.get_indices())
    return fixed


def _normalize_prmtop_lj_tables(parm7_path: str) -> None:
    """Normalize LJ table lengths in parm7 files generated from sliced structures.

    ParmEd slicing can leave ``LENNARD_JONES_*COEF`` longer than the ``POINTERS``
    ``NTYPES`` expectation. Trim only the trailing unused tail when detected.
    """
    from parmed.amber import AmberFormat, AmberParm

    try:
        AmberParm(parm7_path)
        return
    except Exception as exc:
        msg = str(exc)
        if (
            "FLAG LENNARD_JONES_ACOEF" not in msg
            and "FLAG LENNARD_JONES_BCOEF" not in msg
        ):
            raise

    af = AmberFormat(parm7_path)
    pointers = list(af.parm_data.get("POINTERS", []))
    if len(pointers) < 2:
        raise ValueError(f"Invalid POINTERS section in parm7: {parm7_path}")
    ntypes = int(pointers[1])
    expected = ntypes * (ntypes + 1) // 2

    changed = False
    for key in ("LENNARD_JONES_ACOEF", "LENNARD_JONES_BCOEF"):
        values = list(af.parm_data.get(key, []))
        if len(values) == expected:
            continue
        if len(values) < expected:
            raise ValueError(
                f"{key} has {len(values)} entries but expected at least {expected} "
                f"from NTYPES={ntypes} in {parm7_path}."
            )
        af.parm_data[key] = values[:expected]
        changed = True

    if changed:
        af.write_parm(parm7_path)

    # Validate normalized topology immediately.
    AmberParm(parm7_path)


# ======================================================================
#                    hessian_ff (MM) -> ASE calculator
# ======================================================================

def _expand_partial_hessian(
    h_sub: np.ndarray,
    active_atoms: np.ndarray,
    n_atoms: int,
    *,
    dtype: np.dtype,
) -> np.ndarray:
    h_full = np.zeros((3 * n_atoms, 3 * n_atoms), dtype=dtype)
    for i_local, i_atom in enumerate(active_atoms):
        i0 = 3 * int(i_atom)
        for j_local, j_atom in enumerate(active_atoms):
            j0 = 3 * int(j_atom)
            h_full[i0:i0 + 3, j0:j0 + 3] = h_sub[
                3 * i_local:3 * i_local + 3,
                3 * j_local:3 * j_local + 3,
            ]
    return h_full


class hessianffCalculator(Calculator):
    """Calculator for MM. hessian_ff-backed."""

    implemented_properties = ["energy", "forces"]

    def __init__(
        self,
        parm7: str,
        rst7: Optional[str] = None,
        *,
        device: str = "auto",
        cuda_idx: int = 0,
        threads: int = 16,
        **kwargs,
    ):
        super().__init__(**kwargs)

        requested = str(device).lower()
        if requested not in {"auto", "cpu"}:
            raise ValueError(
                "MM backend 'hessian_ff' is CPU-only. "
                f"Got device={device!r}. Use mm_device='cpu' or 'auto'."
            )

        self.device = "cpu"
        self.cuda_idx = int(cuda_idx)
        self.threads = int(threads)
        if self.threads > 0:
            torch.set_num_threads(self.threads)

        self.system = load_system(parm7, device="cpu").to(dtype=torch.float64)
        self.ff = ForceFieldTorch(self.system)
        self.natom = int(self.system.natom)
        self._coords_dtype = torch.float64
        self._coords_device = torch.device("cpu")
        self._coord_buf = torch.empty((self.natom, 3), dtype=self._coords_dtype, device=self._coords_device)

        if rst7 is not None:
            xyz = load_coords(rst7, natom=self.natom, device=self._coords_device, dtype=self._coords_dtype)
            self._coord_buf.copy_(xyz)

    def _positions_to_tensor(self, positions_ang: np.ndarray) -> torch.Tensor:
        arr = np.asarray(positions_ang, dtype=np.float64)
        if arr.shape != (self.natom, 3):
            raise ValueError(
                f"Coordinate shape mismatch for '{type(self).__name__}': "
                f"got {arr.shape}, expected ({self.natom}, 3)."
            )
        self._coord_buf.copy_(torch.as_tensor(arr, dtype=self._coords_dtype, device=self._coords_device))
        return self._coord_buf

    def _energy_forces_from_positions(self, positions_ang: np.ndarray) -> Tuple[float, np.ndarray]:
        xyz = self._positions_to_tensor(positions_ang)
        out, force = self.ff.energy_force(xyz, force_calc_mode="Analytical")
        energy_ev = float(out["E_total"].detach().cpu()) * KCALMOL2EV
        forces_ev = force.detach().cpu().numpy().astype(np.float64, copy=False) * KCALMOL2EV
        return energy_ev, forces_ev

    def calculate(self, atoms: Atoms = None, properties=None, system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        if atoms is None:
            raise ValueError("ASE Atoms is required for MM evaluation.")
        energy_ev, forces_ev = self._energy_forces_from_positions(atoms.get_positions())
        self.results = {"energy": energy_ev, "forces": forces_ev}

    def analytical_hessian(
        self,
        atoms: Atoms,
        *,
        info_path: Optional[str] = None,
        dtype: np.dtype = np.float64,
        return_partial_hessian: bool = False,
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        fixed = _fixed_indices_from_constraints(atoms)
        active_atoms = np.asarray([i for i in range(len(atoms)) if i not in fixed], dtype=int)

        if active_atoms.size == 0:
            if return_partial_hessian:
                return np.zeros((0, 0), dtype=dtype), active_atoms
            return np.zeros((3 * len(atoms), 3 * len(atoms)), dtype=dtype), None

        if info_path is not None:
            dir_ = os.path.dirname(info_path)
            if dir_:
                os.makedirs(dir_, exist_ok=True)
            with open(info_path, "w", encoding="utf-8") as log:
                log.write("Analytical Hessian (hessian_ff)\n")
                log.write("--------------------------------\n")
                log.write(f"n_active_atoms = {active_atoms.size}\n")
                log.flush()

        xyz = self._positions_to_tensor(atoms.get_positions())
        h_local, _ = build_analytical_hessian(
            system=self.system,
            coords=xyz,
            active_atoms=active_atoms.tolist(),
        )
        h_sub = h_local.detach().cpu().numpy().astype(np.float64, copy=False) * KCALMOL2EV
        h_sub = np.asarray(h_sub, dtype=dtype)

        if return_partial_hessian:
            return h_sub, active_atoms

        h_full = _expand_partial_hessian(h_sub, active_atoms, len(atoms), dtype=dtype)
        return h_full, None

    def finite_difference_hessian(
        self,
        atoms: Atoms,
        *,
        delta: float = 1e-3,
        info_path: Optional[str] = None,
        dtype: np.dtype = np.float64,
        return_partial_hessian: bool = False,
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        _ = float(delta)  # Kept for backward-compatible signature.
        return self.analytical_hessian(
            atoms,
            info_path=info_path,
            dtype=dtype,
            return_partial_hessian=return_partial_hessian,
        )



# ======================================================================
#                           OpenMM Calculator
# ======================================================================
class OpenMMCalculator(Calculator):
    """
    ASE Calculator wrapper for OpenMM backend (finite-difference Hessian).

    This calculator uses OpenMM for MM force field evaluation and supports
    CUDA/CPU platforms. Unlike hessianffCalculator, it computes Hessians
    via numerical finite differences.

    Parameters
    ----------
    parm7 : str
        Path to Amber parm7 topology file.
    rst7 : str
        Path to Amber rst7 coordinate file.
    device : str, default "auto"
        Platform selection: "auto", "cuda", or "cpu".
    cuda_idx : int, default 0
        CUDA device index when device="cuda".
    threads : int, default 16
        Number of CPU threads when device="cpu".
    """

    implemented_properties = ["energy", "forces"]

    def __init__(
        self,
        parm7: str,
        rst7: str,
        *,
        device: str = "auto",
        cuda_idx: int = 0,
        threads: int = 16,
        **kwargs,
    ):
        super().__init__(**kwargs)

        if not HAS_OPENMM:
            raise ImportError(
                "OpenMM is required for OpenMMCalculator. "
                "Install with: conda install -c conda-forge openmm"
            )

        # Auto-detect device
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"

        # Platform selection
        if device == "cuda":
            platform = Platform.getPlatformByName("CUDA")
            properties = {
                "CudaDeviceIndex": str(cuda_idx),
                "CudaPrecision": "double",
                "DeterministicForces": "true",
                "CudaUseBlockingSync": "true",
            }
        else:
            platform = Platform.getPlatformByName("CPU")
            properties = {"Threads": str(threads)}

        # Load Amber topology and coordinates
        self.prmtop = app.AmberPrmtopFile(parm7)
        inpcrd = app.AmberInpcrdFile(rst7)

        # Create OpenMM system and context
        self.system = self.prmtop.createSystem(
            nonbondedMethod=app.NoCutoff,
            rigidWater=False
        )
        self.integrator = mm.VerletIntegrator(0 * unit.femtoseconds)
        self.context = mm.Context(self.system, self.integrator, platform, properties)
        self.context.setPositions(inpcrd.positions)

    def calculate(self, atoms: Atoms = None, properties=None, system_changes=all_changes):
        """Compute energy and forces for the given atoms."""
        super().calculate(atoms, properties, system_changes)

        # Define eV unit for OpenMM
        ev_base_unit = ScaledUnit(1.602176634e-19, joule, "electron volt", "eV")
        eV = unit.Unit({ev_base_unit: 1.0})

        # Update positions and get state
        self.context.setPositions(atoms.get_positions() * unit.angstrom)
        state = self.context.getState(getEnergy=True, getForces=True)

        # Extract energy and forces in eV units
        energy = state.getPotentialEnergy().value_in_unit(eV / unit.item)
        forces = state.getForces(asNumpy=True).value_in_unit(eV / unit.angstrom / unit.item)

        self.results = {"energy": energy, "forces": forces}

    def finite_difference_hessian(
        self,
        atoms: Atoms,
        *,
        delta: float = 0.01,
        info_path: Optional[str] = None,
        dtype: np.dtype = np.float64,
        return_partial_hessian: bool = False,
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """
        Compute Hessian via finite differences using hessian_calc utility.

        Parameters
        ----------
        atoms : Atoms
            Structure to differentiate.
        delta : float, default 0.01
            Displacement size in Angstrom.
        info_path : str | None
            Progress log file path.
        dtype : numpy dtype, default float64
            Data type for the Hessian matrix.
        return_partial_hessian : bool, default False
            If True, return only the active sub-Hessian and active atom indices.

        Returns
        -------
        H_full : ndarray
            Full (3N, 3N) Hessian matrix in eV/Å².
        active_atoms : ndarray | None
            Active atom indices (only if return_partial_hessian=True).
        """
        from .hessian_calc import hessian_calc

        H_full = hessian_calc(atoms, self, delta=delta, info_path=info_path, dtype=dtype)

        if return_partial_hessian:
            fixed = _fixed_indices_from_constraints(atoms)
            active_atoms = np.asarray([i for i in range(len(atoms)) if i not in fixed])
            return H_full, active_atoms

        return H_full, None


# ======================================================================
#                               ML/MM Core (UMA)
# ======================================================================

@dataclass(frozen=True)
class _MLHighOut:
    E: float
    F: np.ndarray
    H: Optional[torch.Tensor]
    timing: Dict[str, float | str]


@dataclass(frozen=True)
class _MMLowOut:
    E_real: float
    F_real: np.ndarray
    E_model: float
    F_model: np.ndarray
    H_real: Optional[np.ndarray]
    H_model: Optional[np.ndarray]
    active_atoms_from_fd: Optional[np.ndarray]
    timing: Dict[str, float | str]


class MLMMCore:
    """UMA‑only ONIOM-like ML/MM engine."""

    def __init__(
        self,
        *,
        input_pdb: str,
        real_parm7: str,
        model_pdb: str,
        model_charge: Optional[int] = 0,
        model_mult: int = 1,
        link_mlmm: List[Tuple[str, str]] | None = None,
        uma_model: str = "uma-s-1p1",
        uma_task_name: str = "omol",
        mm_fd: bool = True,
        mm_fd_dir: Optional[str] = None,
        mm_fd_delta: float = 1e-3,
        symmetrize_hessian: bool = True,
        print_timing: bool = True,
        print_vram: bool = True,
        H_double: bool = False,
        ml_device: str = "auto",
        ml_cuda_idx: int = 0,
        mm_backend: str = "hessian_ff",
        mm_device: str = "cpu",
        mm_cuda_idx: int = 0,
        mm_threads: int = 16,
        freeze_atoms: List[int] | None = None,
        ml_hessian_mode: str = "Analytical",
        hessian_calc_mode: Optional[str] = None,
        return_partial_hessian: bool = True,
        hess_cutoff: Optional[float] = None,
        movable_cutoff: Optional[float] = None,
        use_bfactor_layers: bool = False,
        hess_mm_atoms: Optional[List[int]] = None,
        movable_mm_atoms: Optional[List[int]] = None,
        frozen_mm_atoms: Optional[List[int]] = None,
    ):
        self._tmpdir_obj = tempfile.TemporaryDirectory()
        self.tmpdir: str = self._tmpdir_obj.name
        for src, dst in [(input_pdb, "input.pdb"), (real_parm7, "real.parm7"), (model_pdb, "model.pdb")]:
            shutil.copy(src, os.path.join(self.tmpdir, dst))

        self.input_pdb = os.path.join(self.tmpdir, "input.pdb")
        self.real_parm7 = os.path.join(self.tmpdir, "real.parm7")
        self.real_rst7 = os.path.join(self.tmpdir, "real.rst7")
        self.model_pdb = os.path.join(self.tmpdir, "model.pdb")
        self.model_parm7 = os.path.join(self.tmpdir, "model.parm7")
        self.model_rst7 = os.path.join(self.tmpdir, "model.rst7")

        real_top = pmd.load_file(self.real_parm7)
        start_struct = pmd.load_file(self.input_pdb)
        real_n_atoms = int(len(real_top.atoms))
        start_n_atoms = int(len(start_struct.atoms))
        if start_n_atoms != real_n_atoms:
            raise ValueError(
                "Atom-count mismatch between input structure and real topology: "
                f"input_pdb='{input_pdb}' has {start_n_atoms} atoms, "
                f"real_parm7='{real_parm7}' expects {real_n_atoms} atoms. "
                "Provide a full-system input structure consistent with the parm7."
            )
        real_top.coordinates = start_struct.coordinates
        real_top.box = None
        real_top.save(self.real_parm7, overwrite=True)
        real_top.save(self.real_rst7, overwrite=True)

        self.link_mlmm = link_mlmm
        self.ml_ID, self.mlmm_links = self._ml_prep()
        self.selection_indices = self._mk_model_parm7()

        self.hess_cutoff = hess_cutoff
        self.movable_cutoff = movable_cutoff
        self.use_bfactor_layers = use_bfactor_layers
        self._original_input_pdb = input_pdb
        self._explicit_hess_mm_atoms = hess_mm_atoms
        self._explicit_movable_mm_atoms = movable_mm_atoms
        self._explicit_frozen_mm_atoms = frozen_mm_atoms
        self._compute_layer_indices(real_top.coordinates)

        self.freeze_atoms = [] if freeze_atoms is None else list(freeze_atoms)
        if self.frozen_layer_indices:
            self.freeze_atoms = sorted(set(self.freeze_atoms) | set(self.frozen_layer_indices))

        hess_set = set(self.hess_indices)
        all_atoms = set(range(len(real_top.atoms)))
        self.hess_freeze_atoms = sorted(all_atoms - hess_set)

        self.return_partial_hessian = bool(return_partial_hessian)

        self._n_real = len(real_top.atoms)
        self._idx_map_real_to_model = {idx: pos for pos, idx in enumerate(self.selection_indices)}
        self._update_active_dof_mappings()

        self.H_double = bool(H_double)
        self.H_dtype = torch.float64 if self.H_double else torch.float32
        self.H_np_dtype = np.float64 if self.H_double else np.float32

        self.mm_fd = mm_fd
        self.mm_fd_dir = mm_fd_dir
        self.mm_fd_delta = mm_fd_delta
        self.symmetrize_hessian = symmetrize_hessian
        self.print_timing = bool(print_timing)
        self.print_vram = bool(print_vram)
        if self.mm_fd_dir and not os.path.exists(self.mm_fd_dir):
            os.makedirs(self.mm_fd_dir, exist_ok=True)

        if ml_device == "auto":
            ml_device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device_str = ml_device
        self.ml_device = torch.device(f"cuda:{ml_cuda_idx}" if ml_device == "cuda" else "cpu")

        self._AtomicData = AtomicData
        self._data_list_collater = data_list_collater
        self.predictor = pretrained_mlip.get_predict_unit(uma_model, device=self.device_str)
        self.predictor.model.eval()
        for m in self.predictor.model.modules():
            if isinstance(m, nn.Dropout):
                m.p = 0.0

        self.uma_task_name = uma_task_name
        self.model_charge = int(0 if model_charge is None else model_charge)
        self.model_mult = int(model_mult)

        backbone = getattr(self.predictor.model, "module", self.predictor.model).backbone
        self._uma_max_neigh = getattr(backbone, "max_neighbors", None)
        self._uma_radius = getattr(backbone, "cutoff", None)

        # MM backend selection: hessian_ff or openmm
        self.mm_backend = str(mm_backend).strip().lower()
        if self.mm_backend == "openmm":
            self.calc_real_low = OpenMMCalculator(
                parm7=self.real_parm7, rst7=self.real_rst7,
                device=mm_device, cuda_idx=mm_cuda_idx, threads=mm_threads
            )
            self.calc_model_low = OpenMMCalculator(
                parm7=self.model_parm7, rst7=self.model_rst7,
                device=mm_device, cuda_idx=mm_cuda_idx, threads=mm_threads
            )
        elif self.mm_backend == "hessian_ff":
            self.calc_real_low = hessianffCalculator(
                parm7=self.real_parm7, rst7=None,
                device=mm_device, cuda_idx=mm_cuda_idx, threads=mm_threads
            )
            self.calc_model_low = hessianffCalculator(
                parm7=self.model_parm7, rst7=None,
                device=mm_device, cuda_idx=mm_cuda_idx, threads=mm_threads
            )
        else:
            raise ValueError(
                f"Unknown mm_backend '{mm_backend}'. Choose 'hessian_ff' or 'openmm'."
            )

        mode_in = hessian_calc_mode if hessian_calc_mode is not None else ml_hessian_mode
        mode = (mode_in or "Analytical").strip().lower()
        self._ml_hessian_mode = "analytical" if mode.startswith("analyt") else "fd"

        self._atoms_real_tpl = read(self.input_pdb)
        self._atoms_model_tpl = read(self.model_pdb)
        tmp = self._atoms_model_tpl.copy()
        for _ in self.mlmm_links:
            tmp += Atoms("H", positions=[[0.0, 0.0, 0.0]])
        self._atoms_model_LH_tpl = tmp

    @staticmethod
    def _pdb_atom_key(line: str) -> str:
        return f"{line[12:16].strip()} {line[17:20].strip()} {line[22:26].strip()}"

    def _ml_prep(self) -> Tuple[List[str], List[Tuple[int, int]]]:
        ml_region = set()
        with open(self.model_pdb) as fh:
            for ln in fh:
                if ln.startswith(("ATOM", "HETATM")):
                    ml_region.add(self._pdb_atom_key(ln))

        leap_atoms: List[Dict] = []
        with open(self.input_pdb) as fh:
            for ln in fh:
                if not ln.startswith(("ATOM", "HETATM")):
                    continue
                leap_atoms.append(
                    {
                        "idx": int(ln[6:11]),
                        "id": self._pdb_atom_key(ln),
                        "elem": ln[76:78].strip(),
                        "coord": np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]),
                    }
                )

        ml_ID = [str(a["idx"]) for a in leap_atoms if a["id"] in ml_region]

        if self.link_mlmm:
            processed = [(" ".join(q.split()[:3]), " ".join(m.split()[:3])) for q, m in self.link_mlmm]

            ml_indices: List[int] = []
            mm_indices: List[int] = []
            for a in leap_atoms:
                for qnm, mnm in processed:
                    if a["id"] == qnm:
                        ml_indices.append(a["idx"])
                    elif a["id"] == mnm:
                        mm_indices.append(a["idx"])

            if len(set(ml_indices)) != len(ml_indices) or len(set(mm_indices)) != len(mm_indices):
                raise ValueError("Duplicated ML or MM indices in link specification.")
            mlmm_links = list(zip(ml_indices, mm_indices))
        else:
            threshold = 1.7
            ml_set = {a["idx"] for a in leap_atoms if a["id"] in ml_region}
            coords = {a["idx"]: a["coord"] for a in leap_atoms}
            elem = {a["idx"]: a["elem"] for a in leap_atoms}

            ml_indices: List[int] = []
            mm_indices: List[int] = []
            for qidx in ml_set:
                for a in leap_atoms:
                    midx = a["idx"]
                    if midx in ml_set:
                        continue
                    if (
                        np.linalg.norm(coords[midx] - coords[qidx]) < threshold
                        and (
                            (elem[midx] == "C" and elem[qidx] == "C")
                            or (elem[midx] == "N" and elem[qidx] == "C")
                            or (elem[midx] == "C" and elem[qidx] == "N")
                        )
                    ):
                        ml_indices.append(qidx)
                        mm_indices.append(midx)

            if len(set(ml_indices)) != len(ml_indices) or len(set(mm_indices)) != len(mm_indices):
                raise ValueError(
                    "Automatic link detection produced duplicate pairs. Specify 'link_mlmm' manually."
                )
            mlmm_links = list(zip(ml_indices, mm_indices))

        return ml_ID, mlmm_links

    def _mk_model_parm7(self) -> List[int]:
        real = pmd.load_file(self.real_parm7, self.real_rst7)
        real.box = None
        ml_atoms = [real.atoms[int(i) - 1] for i in self.ml_ID]
        selection = [a.idx for a in ml_atoms]

        if len(selection) == len(real.atoms):
            shutil.copy(self.real_parm7, self.model_parm7)
            shutil.copy(self.real_rst7, self.model_rst7)
            return selection

        model = real[selection]
        model.box = None
        model.save(self.model_parm7, overwrite=True)
        _normalize_prmtop_lj_tables(self.model_parm7)
        model.save(self.model_rst7, overwrite=True)
        return selection

    def _compute_layer_indices(self, coords: np.ndarray) -> None:
        self.ml_indices = sorted(self.selection_indices)

        n_atoms = int(coords.shape[0])
        all_indices = set(range(n_atoms))
        mm_indices = all_indices - set(self.ml_indices)

        has_explicit = (
            self._explicit_hess_mm_atoms is not None
            or self._explicit_movable_mm_atoms is not None
            or self._explicit_frozen_mm_atoms is not None
        )
        if has_explicit:
            explicit_hess = set(self._explicit_hess_mm_atoms or [])
            explicit_movable = set(self._explicit_movable_mm_atoms or [])
            explicit_frozen = set(self._explicit_frozen_mm_atoms or [])

            for idx_set, name in [
                (explicit_hess, "hess_mm_atoms"),
                (explicit_movable, "movable_mm_atoms"),
                (explicit_frozen, "frozen_mm_atoms"),
            ]:
                for idx in idx_set:
                    if idx < 0 or idx >= n_atoms:
                        raise ValueError(f"Invalid atom index {idx} in {name}: must be 0 <= idx < {n_atoms}")
                    if idx in self.ml_indices:
                        raise ValueError(f"Atom index {idx} in {name} is also in ML region (model_pdb)")

            self.hess_mm_indices = sorted(explicit_hess & mm_indices)
            self.movable_mm_indices = sorted(explicit_movable & mm_indices)
            self.frozen_layer_indices = sorted(explicit_frozen & mm_indices)

            assigned_mm = explicit_hess | explicit_movable | explicit_frozen
            unassigned_mm = mm_indices - assigned_mm
            self.movable_mm_indices = sorted(set(self.movable_mm_indices) | unassigned_mm)

            self.hess_indices = sorted(self.ml_indices + self.hess_mm_indices)
            self.movable_indices = sorted(self.ml_indices + self.hess_mm_indices + self.movable_mm_indices)
            return

        if self.use_bfactor_layers:
            from .utils import read_bfactors_from_pdb, parse_layer_indices_from_bfactors, has_valid_layer_bfactors
            from pathlib import Path

            bfactors = read_bfactors_from_pdb(Path(self._original_input_pdb))
            if has_valid_layer_bfactors(bfactors):
                layer_info = parse_layer_indices_from_bfactors(bfactors)

                movable_from_layer = set(layer_info["movable_mm_indices"]) & mm_indices
                frozen_from_layer = set(layer_info["frozen_indices"]) & mm_indices
                hess_from_layer = set(layer_info["hess_mm_indices"]) & mm_indices

                # Unassigned MM atoms default to movable.
                assigned_mm = movable_from_layer | frozen_from_layer | hess_from_layer
                unassigned_mm = mm_indices - assigned_mm
                movable_pool = set(movable_from_layer) | set(unassigned_mm)

                # Hessian-target MM selection:
                #   1) If hess_cutoff is set, use distance-to-ML over movable MM pool.
                #   2) Otherwise, keep any Layer-2 assignments (if present).
                hess_mm: set[int]
                if self.hess_cutoff is not None:
                    ml_coords = coords[self.ml_indices]

                    def min_dist_to_ml(atom_idx: int) -> float:
                        atom_coord = coords[atom_idx]
                        dists = np.linalg.norm(ml_coords - atom_coord, axis=1)
                        return float(np.min(dists))

                    hess_cut = float(self.hess_cutoff)
                    hess_mm = {idx for idx in movable_pool if min_dist_to_ml(idx) <= hess_cut}
                else:
                    hess_mm = set(hess_from_layer)

                movable_mm = movable_pool - hess_mm

                self.hess_mm_indices = sorted(hess_mm)
                self.movable_mm_indices = sorted(movable_mm)
                self.frozen_layer_indices = sorted(frozen_from_layer)

                self.hess_indices = sorted(self.ml_indices + self.hess_mm_indices)
                self.movable_indices = sorted(self.ml_indices + self.hess_mm_indices + self.movable_mm_indices)
                return

        if self.hess_cutoff is None and self.movable_cutoff is None:
            self.hess_mm_indices = sorted(mm_indices)
            self.movable_mm_indices = []
            self.frozen_layer_indices = []
            self.hess_indices = sorted(self.ml_indices + self.hess_mm_indices)
            self.movable_indices = sorted(self.ml_indices + self.hess_mm_indices)
            return

        ml_coords = coords[self.ml_indices]

        def min_dist_to_ml(atom_idx: int) -> float:
            atom_coord = coords[atom_idx]
            dists = np.linalg.norm(ml_coords - atom_coord, axis=1)
            return float(np.min(dists))

        hess_mm: List[int] = []
        movable_mm: List[int] = []
        frozen_mm: List[int] = []

        hess_cut = self.hess_cutoff if self.hess_cutoff is not None else float("inf")
        mov_cut = self.movable_cutoff if self.movable_cutoff is not None else float("inf")

        for idx in mm_indices:
            d = min_dist_to_ml(idx)
            if d <= hess_cut:
                hess_mm.append(idx)
            elif d <= mov_cut:
                movable_mm.append(idx)
            else:
                frozen_mm.append(idx)

        self.hess_mm_indices = sorted(hess_mm)
        self.movable_mm_indices = sorted(movable_mm)
        self.frozen_layer_indices = sorted(frozen_mm)
        self.hess_indices = sorted(self.ml_indices + self.hess_mm_indices)
        self.movable_indices = sorted(self.ml_indices + self.hess_mm_indices + self.movable_mm_indices)

    def _update_active_dof_mappings(self) -> None:
        freeze_set = set(self.freeze_atoms)
        self.active_atoms_real = [i for i in range(self._n_real) if i not in freeze_set]
        self.n_active_real = len(self.active_atoms_real)
        self.full_to_active_real = {a: i for i, a in enumerate(self.active_atoms_real)}
        self.active_to_full_real = {i: a for i, a in enumerate(self.active_atoms_real)}

        hess_freeze_set = set(self.hess_freeze_atoms)
        self.hess_active_atoms = [i for i in range(self._n_real) if i not in hess_freeze_set]
        self.n_hess_active = len(self.hess_active_atoms)
        self.full_to_hess_active = {a: i for i, a in enumerate(self.hess_active_atoms)}
        self.hess_active_to_full = {i: a for i, a in enumerate(self.hess_active_atoms)}

        self.ml_hess_active_indices = [
            self.full_to_hess_active[i] for i in self.selection_indices if i in self.full_to_hess_active
        ]

        self.freeze_model = [
            self._idx_map_real_to_model[i] for i in self.freeze_atoms if i in self._idx_map_real_to_model
        ]

    def _build_within_partial_hessian(self) -> Dict[str, np.ndarray | int | str]:
        """Build metadata for a partial (Hessian-target-only) Hessian."""
        n_real = int(self._n_real)
        active_atoms = np.asarray(self.hess_active_atoms, dtype=int)
        active_n_atoms = int(active_atoms.size)

        active_dofs = np.empty(active_n_atoms * 3, dtype=int)
        for i, a in enumerate(active_atoms):
            base = 3 * int(a)
            active_dofs[3 * i:3 * i + 3] = (base, base + 1, base + 2)

        full_to_active = -np.ones(n_real, dtype=int)
        if active_n_atoms:
            full_to_active[active_atoms] = np.arange(active_n_atoms, dtype=int)

        return {
            "kind": "hess-target-only",
            "active_atoms": active_atoms,
            "active_dofs": active_dofs,
            "active_to_full": active_atoms.copy(),
            "full_to_active": full_to_active,
            "full_n_atoms": n_real,
            "full_n_dof": int(3 * n_real),
            "active_n_atoms": active_n_atoms,
            "active_n_dof": int(3 * active_n_atoms),
        }

    def _prep_3_layer_atoms(self, real_coord: np.ndarray):
        atoms_real = self._atoms_real_tpl.copy()
        atoms_real.set_positions(real_coord)

        atoms_model = self._atoms_model_tpl.copy()
        atoms_model_LH = self._atoms_model_LH_tpl.copy()

        for i, ridx in enumerate(self.ml_ID):
            pos = atoms_real[int(ridx) - 1].position
            atoms_model[i].position = pos
            atoms_model_LH[i].position = pos

        added_link_atoms = []
        base_model_len = len(self._atoms_model_tpl)
        for k, (ml_idx, mm_idx) in enumerate(self.mlmm_links):
            ml_i = ml_idx - 1
            mm_i = mm_idx - 1
            ml_elem = atoms_real[ml_i].symbol
            if ml_elem == "C":
                dist = 1.09
            elif ml_elem == "N":
                dist = 1.01
            else:
                raise ValueError(
                    f"Unsupported link parent element: {ml_elem}. Only C and N are supported."
                )
            vec = atoms_real[mm_i].position - atoms_real[ml_i].position
            R = np.linalg.norm(vec)
            if R < 1e-6:
                continue
            u = vec / R
            H_pos = atoms_real[ml_i].position + u * dist
            link_idx_in_model_LH = base_model_len + k
            atoms_model_LH[link_idx_in_model_LH].position = H_pos
            added_link_atoms.append((link_idx_in_model_LH, ml_i, mm_i, dist))

        freeze_model: List[int] = []
        if self.freeze_atoms:
            atoms_real.set_constraint(FixAtoms(indices=self.freeze_atoms))
            real_to_model = self._idx_map_real_to_model
            freeze_model = [real_to_model[i] for i in self.freeze_atoms if i in real_to_model]
            if freeze_model:
                atoms_model.set_constraint(FixAtoms(indices=freeze_model))
                atoms_model_LH.set_constraint(FixAtoms(indices=freeze_model))

        return atoms_real, atoms_model, atoms_model_LH, added_link_atoms, freeze_model

    @staticmethod
    def _jacobian_blocks_numpy(r_ml: np.ndarray, r_mm: np.ndarray, dist: float) -> Optional[np.ndarray]:
        vec = r_mm - r_ml
        R = np.linalg.norm(vec)
        if R < 1e-12:
            return None
        u = vec / R
        I = np.eye(3)
        du_dQ = (I - np.outer(u, u)) / R
        dR_dQ = I - dist * du_dQ
        dR_dM = dist * du_dQ
        return np.hstack([dR_dQ, dR_dM]).T

    @staticmethod
    def _jacobian_blocks_torch(
        r_ml: torch.Tensor,
        r_mm: torch.Tensor,
        dist: float,
        *,
        dtype: torch.dtype,
        device: torch.device,
    ) -> Optional[torch.Tensor]:
        vec = r_mm - r_ml
        Rlen = torch.norm(vec)
        if float(Rlen) < 1e-12:
            return None
        u = vec / Rlen
        I3 = torch.eye(3, dtype=dtype, device=device)
        du_dQ = (I3 - torch.outer(u, u)) / Rlen
        dR_dQ = I3 - dist * du_dQ
        dR_dM = dist * du_dQ
        return torch.hstack([dR_dQ, dR_dM])

    def _uma_eval(self, atoms: Atoms, need_grad: bool = True):
        atoms.info.update({"charge": self.model_charge, "spin": self.model_mult})
        data = self._AtomicData.from_ase(
            atoms,
            max_neigh=self._uma_max_neigh,
            radius=self._uma_radius,
            r_edges=False,
        ).to(self.ml_device)
        data.dataset = self.uma_task_name
        batch = self._data_list_collater([data], otf_graph=True).to(self.ml_device)
        pos = batch.pos.detach().clone().to(self.ml_device)
        pos.requires_grad_(need_grad)
        batch.pos = pos
        if need_grad:
            res = self.predictor.predict(batch)
        else:
            with torch.no_grad():
                res = self.predictor.predict(batch)
        E = float(res["energy"].squeeze().detach().item())
        F = res["forces"].detach().cpu().numpy()
        return E, F, batch

    def _uma_hessian_analytical(self, batch, n_atoms: int) -> torch.Tensor:
        p_flags = [p.requires_grad for p in self.predictor.model.parameters()]
        for p in self.predictor.model.parameters():
            p.requires_grad_(False)

        self.predictor.model.train()
        try:
            pos = batch.pos

            def energy_fn(flat_pos: torch.Tensor):
                batch.pos = flat_pos.view(-1, 3)
                return self.predictor.predict(batch)["energy"].squeeze()

            H_flat = torch.autograd.functional.hessian(energy_fn, pos.view(-1), vectorize=False)
            H = H_flat.view(n_atoms, 3, n_atoms, 3).to(self.H_dtype).detach()
        finally:
            self.predictor.model.eval()
            for p, flag in zip(self.predictor.model.parameters(), p_flags):
                p.requires_grad_(flag)
            if self.ml_device.type == "cuda":
                torch.cuda.empty_cache()
        return H

    def _uma_hessian_fd(
        self,
        atoms: Atoms,
        freeze_model: Sequence[int],
        *,
        eps_ang: float = 1.0e-3,
    ) -> torch.Tensor:
        n_atoms = len(atoms)
        dof = n_atoms * 3
        dev = self.ml_device
        dtype = self.H_dtype

        frozen_set = set(int(i) for i in freeze_model)
        active_atoms = [i for i in range(n_atoms) if i not in frozen_set]
        active_dof_idx = [3 * i + j for i in active_atoms for j in range(3)]

        H = torch.zeros((dof, dof), device=dev, dtype=dtype)
        coord0 = atoms.get_positions().copy()
        for k in active_dof_idx:
            a = k // 3
            c = k % 3

            atoms.positions = coord0.copy()
            atoms.positions[a, c] = coord0[a, c] + eps_ang
            _, Fp, _ = self._uma_eval(atoms, need_grad=False)

            atoms.positions = coord0.copy()
            atoms.positions[a, c] = coord0[a, c] - eps_ang
            _, Fm, _ = self._uma_eval(atoms, need_grad=False)

            col = - (torch.from_numpy(Fp.reshape(-1)) - torch.from_numpy(Fm.reshape(-1))) / (2.0 * eps_ang)
            H[:, k] = col.to(dev, dtype=dtype)

        return H.view(n_atoms, 3, n_atoms, 3)

    def _eval_ml_high(self, atoms_model_LH: Atoms, freeze_model: Sequence[int], *, return_hessian: bool) -> _MLHighOut:
        local_timing: Dict[str, float | str] = {}
        E_model_high, F_model_high, batch_h = self._uma_eval(atoms_model_LH, need_grad=True)

        H_high = None
        if return_hessian:
            n_mlLH = len(atoms_model_LH)
            if self._ml_hessian_mode == "analytical":
                t0 = time.perf_counter()
                H_high = self._uma_hessian_analytical(batch_h, n_mlLH)
                local_timing["ml_hessian_mode"] = "Analytical"
                local_timing["ml_hessian_s"] = time.perf_counter() - t0
            else:
                t0 = time.perf_counter()
                H_high = self._uma_hessian_fd(atoms_model_LH, freeze_model, eps_ang=1.0e-3)
                local_timing["ml_hessian_mode"] = "FiniteDifference"
                local_timing["ml_hessian_s"] = time.perf_counter() - t0

        return _MLHighOut(E=E_model_high, F=F_model_high, H=H_high, timing=local_timing)

    def _eval_mm_low(self, atoms_real: Atoms, atoms_model: Atoms, *, return_hessian: bool) -> _MMLowOut:
        local_timing: Dict[str, float | str] = {}

        atoms_real.calc = self.calc_real_low
        atoms_model.calc = self.calc_model_low

        E_real_low = atoms_real.get_potential_energy()
        F_real_low = np.double(atoms_real.get_forces())

        E_model_low = atoms_model.get_potential_energy()
        F_model_low = np.double(atoms_model.get_forces())

        H_real_np = None
        H_model_np = None
        active_atoms_from_fd = None

        if return_hessian and self.mm_fd is True:
            info_real = os.path.join(self.mm_fd_dir, "real.log") if self.mm_fd_dir else None
            info_model = os.path.join(self.mm_fd_dir, "model.log") if self.mm_fd_dir else None

            atoms_real_for_hess = atoms_real.copy()
            atoms_real_for_hess.calc = self.calc_real_low
            if self.hess_freeze_atoms:
                atoms_real_for_hess.set_constraint(FixAtoms(indices=self.hess_freeze_atoms))

            t0 = time.perf_counter()
            H_real_np, active_atoms_from_fd = self.calc_real_low.finite_difference_hessian(
                atoms_real_for_hess,
                delta=self.mm_fd_delta,
                info_path=info_real,
                dtype=self.H_np_dtype,
                return_partial_hessian=True,
            )
            local_timing["mm_fd_real_s"] = time.perf_counter() - t0

            t0 = time.perf_counter()
            H_model_np, _ = self.calc_model_low.finite_difference_hessian(
                atoms_model,
                delta=self.mm_fd_delta,
                info_path=info_model,
                dtype=self.H_np_dtype,
                return_partial_hessian=False,
            )
            local_timing["mm_fd_model_s"] = time.perf_counter() - t0
            local_timing["mm_fd_total_s"] = float(local_timing["mm_fd_real_s"]) + float(local_timing["mm_fd_model_s"])

        return _MMLowOut(
            E_real=E_real_low,
            F_real=F_real_low,
            E_model=E_model_low,
            F_model=F_model_low,
            H_real=H_real_np,
            H_model=H_model_np,
            active_atoms_from_fd=active_atoms_from_fd,
            timing=local_timing,
        )

    def compute(
        self,
        coord_ang: np.ndarray,
        *,
        return_forces: bool = False,
        return_hessian: bool = False,
    ) -> Dict:
        timing: Dict[str, float | str] = {}
        hess_total_start: Optional[float] = time.perf_counter() if return_hessian else None
        hess_vram_base_alloc: Optional[float] = None
        hess_vram_base_reserved: Optional[float] = None
        hess_vram_total: Optional[float] = None
        if return_hessian and self.print_vram and self.ml_device.type == "cuda":
            torch.cuda.synchronize(device=self.ml_device)
            hess_vram_base_alloc = float(torch.cuda.memory_allocated(device=self.ml_device))
            hess_vram_base_reserved = float(torch.cuda.memory_reserved(device=self.ml_device))
            hess_vram_total = float(torch.cuda.get_device_properties(self.ml_device).total_memory)
            torch.cuda.reset_peak_memory_stats(device=self.ml_device)

        atoms_real, atoms_model, atoms_model_LH, added_link_atoms, freeze_model = self._prep_3_layer_atoms(coord_ang)
        atoms_real.set_pbc(False)
        atoms_model.set_pbc(False)
        atoms_model_LH.set_pbc(False)

        use_parallel = (self.ml_device.type == "cuda") and (getattr(self.calc_real_low, "device", None) == "cpu")
        if use_parallel:
            with ThreadPoolExecutor(max_workers=2) as executor:
                fut_ml = executor.submit(self._eval_ml_high, atoms_model_LH, freeze_model, return_hessian=return_hessian)
                fut_mm = executor.submit(self._eval_mm_low, atoms_real, atoms_model, return_hessian=return_hessian)
                ml_out = fut_ml.result()
                mm_out = fut_mm.result()
        else:
            ml_out = self._eval_ml_high(atoms_model_LH, freeze_model, return_hessian=return_hessian)
            mm_out = self._eval_mm_low(atoms_real, atoms_model, return_hessian=return_hessian)

        timing.update(ml_out.timing)
        timing.update(mm_out.timing)

        total_E = mm_out.E_real + ml_out.E - mm_out.E_model
        results: Dict = {"energy": total_E}

        if return_forces or return_hessian:
            F_combined = np.copy(mm_out.F_real)
            for i, ridx in enumerate(self.selection_indices):
                F_combined[ridx] += ml_out.F[i] - mm_out.F_model[i]

            real_to_model = self._idx_map_real_to_model
            for link_idx, ml_idx, mm_idx, dist in added_link_atoms:
                ml_model_idx = real_to_model[ml_idx]
                r_ml = atoms_model_LH[ml_model_idx].position
                r_mm = atoms_real[mm_idx].position
                grad_link = ml_out.F[link_idx]
                J = self._jacobian_blocks_numpy(r_ml, r_mm, dist)
                if J is None:
                    continue
                redistributed = J @ grad_link
                F_combined[ml_idx] += redistributed[:3]
                F_combined[mm_idx] += redistributed[3:]
            results["forces"] = F_combined

        if return_hessian:
            n_real = len(atoms_real)
            n_ml = len(self.selection_indices)
            n_hess_active = self.n_hess_active

            if self.mm_fd is True:
                if mm_out.H_real is None or mm_out.H_model is None:
                    raise RuntimeError("MM Hessians were not computed as expected.")

                if mm_out.active_atoms_from_fd is not None:
                    expected = set(self.hess_active_atoms)
                    got = set(mm_out.active_atoms_from_fd.tolist())
                    if expected != got:
                        raise RuntimeError(
                            f"Hessian active atoms mismatch: expected {len(expected)} atoms, got {len(got)}"
                        )

                H = torch.from_numpy(mm_out.H_real).to(self.ml_device, self.H_dtype)
                H = H.view(n_hess_active, 3, n_hess_active, 3)

                H_model = torch.from_numpy(mm_out.H_model).to(self.ml_device, self.H_dtype)
                H_model = H_model.view(n_ml, 3, n_ml, 3)
            else:
                H = torch.zeros((n_hess_active, 3, n_hess_active, 3), dtype=self.H_dtype, device=self.ml_device)
                H_model = torch.zeros((n_ml, 3, n_ml, 3), dtype=self.H_dtype, device=self.ml_device)

            H_high = ml_out.H
            ml_pairs = [
                (i, self.full_to_hess_active[gi_real])
                for i, gi_real in enumerate(self.selection_indices)
                if gi_real in self.full_to_hess_active
            ]
            if ml_pairs:
                ml_sel_idx = torch.as_tensor([p[0] for p in ml_pairs], dtype=torch.long, device=self.ml_device)
                ml_active_idx = torch.as_tensor([p[1] for p in ml_pairs], dtype=torch.long, device=self.ml_device)
            else:
                ml_sel_idx = torch.empty((0,), dtype=torch.long, device=self.ml_device)
                ml_active_idx = torch.empty((0,), dtype=torch.long, device=self.ml_device)

            if H_high is not None and ml_sel_idx.numel() > 0:
                t_asm = time.perf_counter()
                H_high_mm = H_high.index_select(0, ml_sel_idx).index_select(2, ml_sel_idx)
                H_model_mm = H_model.index_select(0, ml_sel_idx).index_select(2, ml_sel_idx)
                delta_mm = H_high_mm - H_model_mm
                H[ml_active_idx[:, None], :, ml_active_idx[None, :], :] += delta_mm.permute(0, 2, 1, 3)
                timing["hess_asm_mlml_s"] = time.perf_counter() - t_asm
            del H_model

            real_to_model = self._idx_map_real_to_model
            link_data: List[Tuple[int, int, int, int, int, float, torch.Tensor]] = []
            for link_idx, ml_idx, mm_idx, dist in added_link_atoms:
                ml_model_idx = real_to_model[ml_idx]
                r_ml_t = torch.tensor(atoms_model_LH[ml_model_idx].position, dtype=self.H_dtype, device=self.ml_device)
                r_mm_t = torch.tensor(atoms_real[mm_idx].position, dtype=self.H_dtype, device=self.ml_device)
                K = self._jacobian_blocks_torch(r_ml_t, r_mm_t, dist, dtype=self.H_dtype, device=self.ml_device)
                if K is None:
                    continue
                ml_active = self.full_to_hess_active.get(ml_idx)
                mm_active = self.full_to_hess_active.get(mm_idx)
                if ml_active is None or mm_active is None:
                    continue
                link_data.append((link_idx, ml_idx, mm_idx, ml_active, mm_active, dist, K))

            F_high_t = torch.as_tensor(ml_out.F, dtype=self.H_dtype, device=self.ml_device)
            has_link_force = bool((F_high_t.abs() > 1e-12).any().item())
            if link_data and (H_high is not None or has_link_force):
                t_asm = time.perf_counter()
                I3 = torch.eye(3, dtype=self.H_dtype, device=self.ml_device)
                for link_idx, ml_idx, mm_idx, ml_active, mm_active, dist, K in link_data:

                    if H_high is not None:
                        H_l = H_high[link_idx, :, link_idx, :]
                        H_self = K.T @ H_l @ K
                        H[ml_active, :, ml_active, :].add_(H_self[0:3, 0:3])
                        H[ml_active, :, mm_active, :].add_(H_self[0:3, 3:6])
                        H[mm_active, :, ml_active, :].add_(H_self[3:6, 0:3])
                        H[mm_active, :, mm_active, :].add_(H_self[3:6, 3:6])

                    f_L = -F_high_t[link_idx]

                    r_ml_t = torch.as_tensor(atoms_model_LH[real_to_model[ml_idx]].position,
                                             dtype=self.H_dtype, device=self.ml_device)
                    r_mm_t = torch.as_tensor(atoms_real[mm_idx].position,
                                             dtype=self.H_dtype, device=self.ml_device)
                    v = r_mm_t - r_ml_t
                    R_sq = torch.dot(v, v)
                    inv_R = torch.rsqrt(torch.clamp(R_sq, min=1.0e-24))
                    inv_R2 = inv_R * inv_R
                    u = v * inv_R

                    alpha = torch.dot(u, f_L)
                    uuT = torch.outer(u, u)
                    ufT = torch.outer(u, f_L)
                    fTu = torch.outer(f_L, u)
                    B = (alpha * (3.0 * uuT - I3) - (ufT + fTu)) * inv_R2

                    H_corr6 = torch.zeros((6, 6), dtype=self.H_dtype, device=self.ml_device)
                    H_corr6[0:3, 0:3] = B
                    H_corr6[3:6, 3:6] = B
                    H_corr6[0:3, 3:6] = -B
                    H_corr6[3:6, 0:3] = -B
                    H_corr6.mul_(dist)

                    H[ml_active, :, ml_active, :].add_(H_corr6[0:3, 0:3])
                    H[ml_active, :, mm_active, :].add_(H_corr6[0:3, 3:6])
                    H[mm_active, :, ml_active, :].add_(H_corr6[3:6, 0:3])
                    H[mm_active, :, mm_active, :].add_(H_corr6[3:6, 3:6])
                timing["hess_asm_link_self_s"] = time.perf_counter() - t_asm

            if H_high is not None and link_data and ml_sel_idx.numel() > 0:
                t_asm = time.perf_counter()
                for link_idx, _ml_idx, _mm_idx, ml_active, mm_active, _dist, K in link_data:
                    H_coup = H_high[link_idx].index_select(1, ml_sel_idx).permute(1, 0, 2).contiguous()  # (K,3,3)
                    H_row = torch.einsum("ac,bcd->bad", K.T, H_coup)  # (K,6,3)
                    H_col = torch.einsum("bca,cd->bad", H_coup, K)    # (K,3,6)

                    # Mixed scalar/tensor indexing in PyTorch returns (3, K, 3) for
                    # H[scalar, :, tensor, :], so align H_row blocks explicitly.
                    H[ml_active, :, ml_active_idx, :].add_(H_row[:, 0:3, :].permute(1, 0, 2))
                    H[mm_active, :, ml_active_idx, :].add_(H_row[:, 3:6, :].permute(1, 0, 2))
                    H[ml_active_idx, :, ml_active, :].add_(H_col[:, :, 0:3])
                    H[ml_active_idx, :, mm_active, :].add_(H_col[:, :, 3:6])
                timing["hess_asm_link_ml_s"] = time.perf_counter() - t_asm

            if H_high is not None and link_data:
                t_asm = time.perf_counter()
                n_links = len(link_data)
                for a in range(n_links):
                    link_idx_a, _ml_a, _mm_a, ml_a_active, mm_a_active, _dist_a, K_a = link_data[a]

                    for b in range(a + 1, n_links):
                        link_idx_b, _ml_b, _mm_b, ml_b_active, mm_b_active, _dist_b, K_b = link_data[b]

                        H_ab = H_high[link_idx_a, :, link_idx_b, :]
                        HAB = K_a.T @ H_ab @ K_b

                        H[ml_a_active, :, ml_b_active, :].add_(HAB[0:3, 0:3])
                        H[ml_a_active, :, mm_b_active, :].add_(HAB[0:3, 3:6])
                        H[mm_a_active, :, ml_b_active, :].add_(HAB[3:6, 0:3])
                        H[mm_a_active, :, mm_b_active, :].add_(HAB[3:6, 3:6])

                        HBA = HAB.T
                        H[ml_b_active, :, ml_a_active, :].add_(HBA[0:3, 0:3])
                        H[ml_b_active, :, mm_a_active, :].add_(HBA[0:3, 3:6])
                        H[mm_b_active, :, ml_a_active, :].add_(HBA[3:6, 0:3])
                        H[mm_b_active, :, mm_a_active, :].add_(HBA[3:6, 3:6])
                timing["hess_asm_link_link_s"] = time.perf_counter() - t_asm

            if self.symmetrize_hessian:
                t_asm = time.perf_counter()
                H_flat = H.view(3 * n_hess_active, 3 * n_hess_active)
                H_flat = (H_flat + H_flat.t()).mul_(0.5)
                H = H_flat.view(n_hess_active, 3, n_hess_active, 3)
                timing["hess_asm_sym_s"] = time.perf_counter() - t_asm

            if self.return_partial_hessian:
                results["hessian"] = H.detach()
                results["within_partial_hessian"] = self._build_within_partial_hessian()
            else:
                t_asm = time.perf_counter()
                H_full = torch.zeros((n_real, 3, n_real, 3), dtype=self.H_dtype, device=self.ml_device)
                active_idx = torch.as_tensor(self.hess_active_atoms, dtype=torch.long, device=self.ml_device)
                if active_idx.numel() > 0:
                    H_full[active_idx[:, None], :, active_idx[None, :], :] = H.permute(0, 2, 1, 3).contiguous()
                results["hessian"] = H_full.detach()
                timing["hess_asm_full_expand_s"] = time.perf_counter() - t_asm
                del H_full

            if hess_total_start is not None:
                timing["hessian_total_s"] = time.perf_counter() - hess_total_start
                results["timing"] = timing
                if self.print_timing:
                    ml_mode = timing.get("ml_hessian_mode")
                    ml_time = timing.get("ml_hessian_s")
                    if ml_mode is not None and ml_time is not None:
                        click.echo(f"[HessianTiming] ML Hessian ({ml_mode}): {ml_time:.2f} s")
                    if "mm_fd_total_s" in timing:
                        click.echo(
                            f"[HessianTiming] MM Hessian: REAL {timing['mm_fd_real_s']:.2f} s | "
                            f"MODEL {timing['mm_fd_model_s']:.2f} s | "
                            f"total {timing['mm_fd_total_s']:.2f} s"
                        )
                    asm_parts = []
                    for key, label in (
                        ("hess_asm_mlml_s", "ML-ML"),
                        ("hess_asm_link_self_s", "link-self"),
                        ("hess_asm_link_ml_s", "link-ML"),
                        ("hess_asm_link_link_s", "link-link"),
                        ("hess_asm_sym_s", "sym"),
                        ("hess_asm_full_expand_s", "full-expand"),
                    ):
                        if key in timing:
                            asm_parts.append(f"{label} {float(timing[key]):.2f} s")
                    if asm_parts:
                        click.echo(f"[HessianTiming] Assembly: {' | '.join(asm_parts)}")
                    click.echo(f"[HessianTiming] Hessian total: {timing['hessian_total_s']:.2f} s")
                if self.print_vram and self.ml_device.type == "cuda":
                    torch.cuda.synchronize(device=self.ml_device)
                    base_alloc = float(hess_vram_base_alloc or 0.0)
                    base_reserved = float(hess_vram_base_reserved or 0.0)
                    peak_alloc = max(
                        float(torch.cuda.max_memory_allocated(device=self.ml_device)) - base_alloc,
                        0.0,
                    ) / 1e9
                    peak_reserved_abs = float(torch.cuda.max_memory_reserved(device=self.ml_device))
                    peak_reserved = max(
                        peak_reserved_abs - base_reserved,
                        0.0,
                    ) / 1e9
                    total_vram = float(hess_vram_total or torch.cuda.get_device_properties(self.ml_device).total_memory) / 1e9
                    remaining_vram = max((total_vram * 1e9) - peak_reserved_abs, 0.0) / 1e9
                    click.echo(
                        f"[HessianVRAM] total={total_vram:.3f} GB | "
                        f"peak_allocated={peak_alloc:.3f} GB | "
                        f"peak_reserved={peak_reserved:.3f} GB | "
                        f"remaining={remaining_vram:.3f} GB"
                    )

            del H, H_high
            if self.ml_device.type == "cuda":
                torch.cuda.empty_cache()

        return results


# ======================================================================
#                ASE Calculator wrapper for ML/MM (ONIOM)
# ======================================================================

class MLMMASECalculator(Calculator):
    """ASE Calculator wrapping MLMMCore for use with DMF and other ASE-based methods.

    The underlying MLMMCore takes full-system coordinates (Angstrom) and
    returns energy in eV and forces in eV/Angstrom, which matches ASE conventions.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, core: "MLMMCore", **kwargs):
        super().__init__(**kwargs)
        self.core = core

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        coord_ang = atoms.get_positions().astype(float)
        want_forces = "forces" in properties
        res = self.core.compute(coord_ang, return_forces=want_forces, return_hessian=False)
        self.results = {
            "energy": float(res["energy"]),
        }
        if want_forces:
            self.results["forces"] = res["forces"].reshape(-1, 3)


# ======================================================================
#                     PySisyphus Calculator (ML/MM)
# ======================================================================

from pysisyphus.calculators.Calculator import Calculator as PySiCalc


class mlmm(PySiCalc):
    implemented_properties = ["energy", "forces", "hessian"]

    def __init__(
        self,
        input_pdb: str = None,
        real_parm7: str = None,
        model_pdb: str = None,
        *,
        model_charge: int = 0,
        model_mult: int = 1,
        link_mlmm: List[Tuple[str, str]] | None = None,
        uma_model: str = "uma-s-1p1",
        uma_task_name: str = "omol",
        mm_fd: bool = True,
        mm_fd_dir: Optional[str] = None,
        mm_fd_delta: float = 1e-3,
        symmetrize_hessian: bool = True,
        out_hess_torch: bool = True,
        H_double: bool = False,
        ml_hessian_mode: str = "Analytical",
        hessian_calc_mode: Optional[str] = None,
        ml_device: str = "auto",
        ml_cuda_idx: int = 0,
        mm_device: str = "cpu",
        mm_cuda_idx: int = 0,
        mm_threads: int = 16,
        mm_backend: str = "hessian_ff",
        freeze_atoms: List[int] | None = None,
        return_partial_hessian: bool = True,
        print_timing: bool = True,
        print_vram: bool = True,
        hess_cutoff: Optional[float] = None,
        movable_cutoff: Optional[float] = None,
        use_bfactor_layers: bool = False,
        hess_mm_atoms: Optional[List[int]] = None,
        movable_mm_atoms: Optional[List[int]] = None,
        frozen_mm_atoms: Optional[List[int]] = None,
        **kwargs,
    ):
        self._freeze_atoms = [] if freeze_atoms is None else list(freeze_atoms)
        super().__init__(charge=model_charge, mult=model_mult, **kwargs)

        self.core = MLMMCore(
            input_pdb=input_pdb,
            real_parm7=real_parm7,
            model_pdb=model_pdb,
            model_charge=model_charge,
            model_mult=model_mult,
            link_mlmm=link_mlmm,
            uma_model=uma_model,
            uma_task_name=uma_task_name,
            mm_fd=mm_fd,
            mm_fd_dir=mm_fd_dir,
            mm_fd_delta=mm_fd_delta,
            symmetrize_hessian=symmetrize_hessian,
            H_double=H_double,
            ml_device=ml_device,
            ml_cuda_idx=ml_cuda_idx,
            mm_device=mm_device,
            mm_cuda_idx=mm_cuda_idx,
            mm_threads=mm_threads,
            mm_backend=mm_backend,
            freeze_atoms=self._freeze_atoms,
            ml_hessian_mode=ml_hessian_mode,
            hessian_calc_mode=hessian_calc_mode,
            return_partial_hessian=return_partial_hessian,
            print_timing=print_timing,
            print_vram=print_vram,
            hess_cutoff=hess_cutoff,
            movable_cutoff=movable_cutoff,
            use_bfactor_layers=use_bfactor_layers,
            hess_mm_atoms=hess_mm_atoms,
            movable_mm_atoms=movable_mm_atoms,
            frozen_mm_atoms=frozen_mm_atoms,
        )

        self.out_hess_torch = bool(out_hess_torch)
        self.hess_torch_double = bool(H_double)
        self._hess_scale = EV2AU / ANG2BOHR / ANG2BOHR

    @property
    def freeze_atoms(self) -> List[int] | None:
        return self.core.freeze_atoms

    @freeze_atoms.setter
    def freeze_atoms(self, indices: List[int] | None):
        self._freeze_atoms = [] if indices is None else list(indices)
        self.core.freeze_atoms = self._freeze_atoms
        self.core._update_active_dof_mappings()

    def _run_core(self, coords, *, want_forces: bool, want_hessian: bool):
        coord_ang = np.asarray(coords).reshape(-1, 3) * BOHR2ANG
        res = self.core.compute(coord_ang, return_forces=want_forces or want_hessian, return_hessian=want_hessian)
        out = {"energy": res["energy"] * EV2AU}
        if want_forces or want_hessian:
            out["forces"] = (res["forces"] * (EV2AU / ANG2BOHR)).flatten()
        if want_hessian:
            H = res.pop("hessian")
            H = H.view(H.size(0) * 3, H.size(2) * 3)
            H.mul_(self._hess_scale)
            if self.out_hess_torch:
                target_dtype = torch.float64 if self.hess_torch_double else torch.float32
                out["hessian"] = H.to(target_dtype).detach().requires_grad_(False)
            else:
                out["hessian"] = H.detach().cpu().numpy()
            if "within_partial_hessian" in res:
                out["within_partial_hessian"] = res["within_partial_hessian"]
        return out

    def get_energy(self, elem, coords):
        return self._run_core(coords, want_forces=False, want_hessian=False)

    def get_forces(self, elem, coords):
        return self._run_core(coords, want_forces=True, want_hessian=False)

    def get_hessian(self, elem, coords):
        return self._run_core(coords, want_forces=True, want_hessian=True)


# ======================================================================
#                     PySisyphus Calculator (MM-only)
# ======================================================================


class mlmm_mm_only(PySiCalc):
    """PySisyphus calculator that returns MM-only energy and forces (F_real_mm).

    Used for microiteration: relaxes the MM region without ML computation.
    Shares the MLMMCore from an existing ``mlmm`` calculator to avoid
    re-initializing topology and force field objects.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, core: "MLMMCore", *, freeze_atoms: list[int] | None = None, **kwargs):
        super().__init__(charge=core.model_charge, mult=core.model_mult, **kwargs)
        self.core = core
        self._freeze_atoms = list(freeze_atoms) if freeze_atoms else []

    def _run_core(self, coords, *, want_forces: bool):
        coord_ang = np.asarray(coords).reshape(-1, 3) * BOHR2ANG
        atoms_real = self.core._atoms_real_tpl.copy()
        atoms_real.set_positions(coord_ang)
        atoms_real.set_pbc(False)
        atoms_real.calc = self.core.calc_real_low
        E_real = float(atoms_real.get_potential_energy())
        out = {"energy": E_real * EV2AU}
        if want_forces:
            F_real = np.double(atoms_real.get_forces())
            # Zero forces on frozen atoms
            for i in self._freeze_atoms:
                if 0 <= i < F_real.shape[0]:
                    F_real[i, :] = 0.0
            out["forces"] = (F_real * (EV2AU / ANG2BOHR)).flatten()
        return out

    def get_energy(self, elem, coords):
        return self._run_core(coords, want_forces=False)

    def get_forces(self, elem, coords):
        return self._run_core(coords, want_forces=True)

    def get_hessian(self, elem, coords):
        raise NotImplementedError("MM-only calculator does not support Hessian computation.")


# ======================================================================
#                           CLI registration
# ======================================================================

from pysisyphus import run as _run


def run_pysis_mlmm():
    _run.CALC_DICT["mlmm"] = mlmm
    _run.run()


if __name__ == "__main__":
    run_pysis_mlmm()
