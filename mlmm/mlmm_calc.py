"""
ML/MM engine (REAL-MM / MODEL-MM / MODEL-High) with
ML/MM link-atom Jacobian redistribution and optional Hessian support.

* High-level ML part : uma-s predictor (fairchem) **or** AIMNet2
* Low-level MM part  : OpenMM (Amber ff14SB)

Select the ML backend with the argument ``backend`` ("uma" | "aimnet2").
The default is "uma".
"""

import os
import shutil
import tempfile
from typing import List, Tuple, Dict

import numpy as np
import torch
import torch.nn as nn

from ase import Atoms
from ase.io import read
from ase.calculators.calculator import Calculator, all_changes
from ase.data import atomic_numbers
from ase.constraints import FixAtoms

import parmed as pmd
import openmm as mm
from openmm import app, unit, Platform
from openmm.unit import ScaledUnit, Unit, joule

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

from mlmm.hessian_calc import hessian_calc

# UMA / fairchem and AIMNet2 is imported lazily inside the class to avoid hard dependency

# ---------------------------------------------------------------------
# OpenMM → ASE calculator wrapper
# ---------------------------------------------------------------------
ev_base_unit = ScaledUnit(1.602176634e-19, joule, "electron volt", "eV")
eV = Unit({ev_base_unit: 1.0})


class OpenMMCalculator(Calculator):
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

        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"

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

        self.prmtop = app.AmberPrmtopFile(parm7)
        inpcrd = app.AmberInpcrdFile(rst7)

        self.system = self.prmtop.createSystem(
            nonbondedMethod=app.NoCutoff, rigidWater=False
        )
        self.integrator = mm.VerletIntegrator(0 * unit.femtoseconds)
        self.context = mm.Context(self.system, self.integrator, platform, properties)
        self.context.setPositions(inpcrd.positions)

    # ---------------------------------------------------------------------
    def calculate(self, atoms: Atoms = None, properties=None, system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)

        self.context.setPositions(atoms.get_positions() * unit.angstrom)
        state = self.context.getState(getEnergy=True, getForces=True)

        energy = state.getPotentialEnergy().value_in_unit(eV / unit.item)
        forces = state.getForces(asNumpy=True).value_in_unit(eV / unit.angstrom / unit.item)

        self.results = {"energy": energy, "forces": forces}


# ---------------------------------------------------------------------
# MLMMCore – (UMA or AIMNet2)
# ---------------------------------------------------------------------
class MLMMCore:
    """
    Three-layer ML/MM engine with selectable ML backend.

      • REAL   layer : full MM (OpenMM)
      • MODEL-MM     : ML subset, MM evaluation (OpenMM)
      • MODEL-High   : ML + link-H, high-level ML

    Link-atom forces/Hessians are redistributed to ML/MM atoms through
    Jacobian transformation.
    """

    # ---------------------------------------------------------------------
    def __init__(
        self,
        *,
        # === input files =================================================
        real_pdb: str,
        real_parm7: str,
        real_rst7: str,
        model_pdb: str,
        # === model options ==============================================
        model_charge: int | None = None,
        model_mult: int = 1,
        link_mlmm: List[Tuple[str, str]] | None = None,
        backend: str = "uma",  # "uma" or "aimnet2"
        # === hessian / vib =============================================
        vib_run: bool = False,
        vib_dir: str | None = None,
        H_double: bool = False,
        # === device selection ===========================================
        ml_device: str = "auto",
        ml_cuda_idx: int = 0,
        mm_device: str = "cpu",
        mm_cuda_idx: int = 0,
        mm_threads: int = 16,
        freeze_atoms: List[int] | None = None,
    ):
        """
        Args:
            real_pdb (str): Path to the PDB file for the real system.
            real_parm7 (str): Path to the Amber .parm7 file for the real system.
            real_rst7 (str): Path to the Amber .rst7 file for the real system.
            model_pdb (str): Path to the PDB file for the model system.

            model_charge (int): Charge of the model system. If None, charge of model system is calculated automatically with RDKit.
            model_mult (int): Multiplicity of the model system. Default is 1 (singlet). UMA backend only.
            link_mlmm (List[Tuple[str, str]] | None): List of tuples specifying the link atoms between ML and MM regions. e.g.) [("CB  ARG   294", "CA  ARG   294")]. If None, link atoms are determined automatically based on distance and element type.

            vib_run (bool): Whether to run vibrational analysis.

            ml_device (str): Device to use for ML calculations. Options are "auto", "cpu", or "cuda".
            ml_cuda_idx (int): CUDA device index if using GPU.

            mm_device (str): Device to use for calculations. Options are "auto", "cpu", or "cuda".
            mm_cuda_idx (int): CUDA device index if using GPU.
            mm_threads (int): Number of threads to use for CPU calculations.
            freeze_atoms (List[int] | None): 0-based indices of atoms to freeze during MM Hessian calculations.
            H_double (bool): Use double precision for Hessian-related tensors.
        """
        self.backend = backend.lower()
        if self.backend not in ("uma", "aimnet2"):
            raise ValueError("backend must be 'uma' or 'aimnet2'")

        # prepare sandbox dir
        # ---------------------------------------------------------------------
        self._tmpdir_obj = tempfile.TemporaryDirectory()
        self.tmpdir: str = self._tmpdir_obj.name

        for src, dst in [
            (real_pdb, "real.pdb"),
            (real_parm7, "real.parm7"),
            (real_rst7, "real.rst7"),
            (model_pdb, "model.pdb"),
        ]:
            shutil.copy(src, os.path.join(self.tmpdir, dst))

        # standardized filenames
        self.real_pdb = os.path.join(self.tmpdir, "real.pdb")
        self.real_parm7 = os.path.join(self.tmpdir, "real.parm7")
        self.real_rst7 = os.path.join(self.tmpdir, "real.rst7")
        self.model_pdb = os.path.join(self.tmpdir, "model.pdb")
        self.model_parm7 = os.path.join(self.tmpdir, "model.parm7")
        self.model_rst7 = os.path.join(self.tmpdir, "model.rst7")
        
        mol = pmd.load_file(self.real_parm7, self.real_rst7)
        mol.box = None
        mol.save(self.real_parm7, overwrite=True)
        mol.save(self.real_rst7, overwrite=True)

        self.link_mlmm = link_mlmm

        if backend == "aimnet2" and model_mult != 1:
            raise ValueError("Multiplicity except 1 is not supported for AIMNet2 backend.")

        # derive ML region + link pairs
        # ---------------------------------------------------------------------
        self.ml_ID, self.mlmm_links = self._ml_prep()
        self.selection_indices = self._mk_model_parm7()

        self.freeze_atoms = [] if freeze_atoms is None else list(freeze_atoms)

        self.H_double = bool(H_double)
        self.H_dtype = torch.float64 if self.H_double else torch.float32
        self.H_np_dtype = np.float64 if self.H_double else np.float32

        # (optional) vib directory
        # ---------------------------------------------------------------------
        self.vib_run = vib_run
        self.vib_dir = vib_dir or self.tmpdir
        if self.vib_run and not os.path.exists(self.vib_dir):
            os.makedirs(self.vib_dir, exist_ok=True)

        # device selection
        # ---------------------------------------------------------------------
        if ml_device == "auto":
            ml_device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device_str = ml_device
        self.ml_device = torch.device(
            f"cuda:{ml_cuda_idx}" if ml_device == "cuda" else "cpu"
        )

        if mm_device == "auto":
            mm_device = "cuda" if torch.cuda.is_available() else "cpu"
        self.mm_device = mm_device

        # High-level predictor setup
        # ---------------------------------------------------------------------
        # lazy import to avoid mandatory dependency
        if self.backend == "uma":          
            try:
                from fairchem.core import pretrained_mlip
                from fairchem.core.datasets.atomic_data import AtomicData
                from fairchem.core.datasets import data_list_collater
            except ImportError:
                raise ImportError(
                    "Please install fairchem-core to use the 'uma' backend. \n"
                    "pip install fairchem-core\n"
                    "huggingface-cli login\n"
                )

            self._AtomicData = AtomicData
            self._data_list_collater = data_list_collater

            self.predictor = pretrained_mlip.get_predict_unit("uma-s-1", device=self.device_str)

            self.predictor.model.eval()
            for m in self.predictor.model.modules():
                if isinstance(m, nn.Dropout):
                    m.p = 0.0 # set dropout rate to 0.0 for Hessian evaluation
                    m.eval()

        else:  # AIMNet2
            try:
                from aimnet.calculators import AIMNet2Calculator
            except ImportError:
                raise ImportError(
                    "Please install AIMNet2 to use the 'aimnet2' backend. \n"
                    "pip install git+https://github.com/isayevlab/aimnetcentral.git\n"
                )

            model = AIMNet2Calculator("aimnet2")#, device=self.ml_device)
            model.set_lrcoulomb_method("simple")
            self.calc_model_high = model

        # total charge of MODEL system
        self.model_charge = (model_charge if model_charge is not None else self._calc_model_charge())
        self.model_mult = model_mult

        if backend == "aimnet2":
            if self.model_mult != 1:
                raise ValueError("AIMNet2 backend does not support multiplicity other than 1.")

        # Low-level OpenMM calculators
        # ---------------------------------------------------------------------
        self.calc_real_low = OpenMMCalculator(
            parm7=self.real_parm7,
            rst7=self.real_rst7,
            device=self.mm_device,
            cuda_idx=mm_cuda_idx,
            threads=mm_threads,
        )
        self.calc_model_low = OpenMMCalculator(
            parm7=self.model_parm7,
            rst7=self.model_rst7,
            device=self.mm_device,
            cuda_idx=mm_cuda_idx,
            threads=mm_threads,
        )

    # ==================================================================
    #                      INTERNAL HELPERS
    # ==================================================================
    # The following two helpers are backend-specific; we keep both and use
    # the appropriate one at runtime.
    def _ase_to_batch(self, atoms: Atoms):
        """Convert ASE Atoms → UMA AtomicData(Batch). Only used when backend='uma'."""
        atoms.info.update({"charge": self.model_charge, "spin": self.model_mult})
        data = self._AtomicData.from_ase(
            atoms,
            max_neigh=self.predictor.model.module.backbone.max_neighbors,
            radius=self.predictor.model.module.backbone.cutoff,
            r_edges=False,
        ).to(self.ml_device)
        data.dataset = "omol"
        return self._data_list_collater([data], otf_graph=True).to(self.ml_device)

    def _prepare_input(self, elem, coord):
        """Prepare AIMNet2 input. Only used when backend='aimnet2'."""
        numbers = torch.as_tensor(
            [atomic_numbers[symbol] for symbol in elem], 
            dtype=torch.long,
            device=self.ml_device
        )
        coord_t = torch.as_tensor(coord, dtype=torch.float, device=self.ml_device)
        charge = torch.as_tensor([self.model_charge], dtype=torch.float, device=self.ml_device)
        mult = torch.as_tensor([1], dtype=torch.float, device=self.ml_device)
        return dict(coord=coord_t, numbers=numbers, charge=charge, mult=mult)

    # ---------------------------------------------------------------------
    def _ml_prep(self) -> Tuple[List[str], List[Tuple[int, int]]]:
        """Detect ML atom IDs and link pairs (ML idx, MM idx)."""
        # (1) ML ids from model_pdb (residue-based string)
        ml_region = set()
        with open(self.model_pdb) as fh:
            for ln in fh:
                if ln.startswith(("ATOM", "HETATM")):
                    ml_region.add(
                        f"{ln[12:16].strip()} {ln[17:20].strip()} {ln[22:26].strip()}"
                    )

        # (2) parse REAL pdb
        leap_atoms: List[Dict] = []
        with open(self.real_pdb) as fh:
            for ln in fh:
                if not ln.startswith(("ATOM", "HETATM")):
                    continue
                leap_atoms.append(
                    {
                        "idx": int(ln[6:11]),
                        "id": f"{ln[12:16].strip()} {ln[17:20].strip()} {ln[22:26].strip()}",
                        "elem": ln[76:78].strip(),
                        "coord": np.array(
                            [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]
                        ),
                    }
                )

        # (3) ML indices
        ml_ID = [str(a["idx"]) for a in leap_atoms if a["id"] in ml_region]

        # (4) link atom pairs
        if self.link_mlmm:  # user supplied
            processed = [
                (" ".join(q.split()[:3]), " ".join(m.split()[:3])) for q, m in self.link_mlmm
            ]
            ml_indices, mm_indices = [], []
            for a in leap_atoms:
                for qnm, mnm in processed:
                    if a["id"] == qnm:
                        ml_indices.append(a["idx"])
                    elif a["id"] == mnm:
                        mm_indices.append(a["idx"])
            if len(set(ml_indices)) != len(ml_indices) or len(set(mm_indices)) != len(mm_indices):
                raise ValueError("Duplicated ML or MM indices in link specification.")
            mlmm_links = list(zip(ml_indices, mm_indices))
        else:  # automatic neighbour detection
            threshold = 1.7
            ml_set = {a["idx"] for a in leap_atoms if a["id"] in ml_region}
            coords = {a["idx"]: a["coord"] for a in leap_atoms}
            elem = {a["idx"]: a["elem"] for a in leap_atoms}

            ml_indices, mm_indices = [], []
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
                    "Automatic link detection produced duplicate pairs. "
                    "Specify 'link_mlmm' manually."
                )
            mlmm_links = list(zip(ml_indices, mm_indices))
        return ml_ID, mlmm_links

    # ---------------------------------------------------------------------
    def _mk_model_parm7(self):
        """Extract MODEL-ML subset topology."""
        real = pmd.load_file(self.real_parm7, self.real_rst7)
        real.box = None

        ml_atoms = [real.atoms[int(i) - 1] for i in self.ml_ID]
        selection = [a.idx for a in ml_atoms]

        # ML = REAL ?
        if len(selection) == len(real.atoms):
            shutil.copy(self.real_parm7, self.model_parm7)
            shutil.copy(self.real_rst7, self.model_rst7)
            return selection

        model = real[selection]
        # model.recalculate_LJ()
        # model.remake_parm()
        model.box = None
        model.save(self.model_parm7, overwrite=True)
        model.save(self.model_rst7, overwrite=True)
        return selection

    # ---------------------------------------------------------------------
    def _prep_3_layer_atoms(self, real_coord):
        """Build REAL / MODEL(MM) / MODEL(High+linkH) Atoms objects."""
        atoms_real = read(self.real_pdb)
        atoms_real.set_positions(real_coord)

        atoms_model = read(self.model_pdb)
        for i, ridx in enumerate(self.ml_ID):
            atoms_model[i].position = atoms_real[int(ridx) - 1].position

        atoms_model_LH = read(self.model_pdb)
        for i, ridx in enumerate(self.ml_ID):
            atoms_model_LH[i].position = atoms_real[int(ridx) - 1].position

        added_link_atoms = []
        for ml_idx, mm_idx in self.mlmm_links:
            ml_i = ml_idx - 1
            mm_i = mm_idx - 1

            ml_elem = atoms_real[ml_i].symbol
            if ml_elem == "C":
                dist = 1.09
            elif ml_elem == "N":
                dist = 1.01
            else:
                raise ValueError(
                    f"Unsupported link atom element: {ml_elem}. "
                    "Only C and N are supported for parent of link atoms."
                )

            vec = atoms_real[mm_i].position - atoms_real[ml_i].position
            R = np.linalg.norm(vec)
            if R < 1e-6:
                continue  # skip if the atoms are at the same position
            u = vec / np.linalg.norm(vec)
            H_pos = atoms_real[ml_i].position + u * dist
            atoms_model_LH += Atoms("H", positions=[H_pos])
            # (link_idx_in_model_LH , ml_idx_in_REAL , mm_idx_in_REAL , dist)
            added_link_atoms.append((len(atoms_model_LH) - 1, ml_i, mm_i, dist))

        if self.freeze_atoms:
            atoms_real.set_constraint(FixAtoms(indices=self.freeze_atoms))
            idx_map = {idx: pos for pos, idx in enumerate(self.selection_indices)}
            freeze_model = [idx_map[i] for i in self.freeze_atoms if i in idx_map]
            if freeze_model:
                atoms_model.set_constraint(FixAtoms(indices=freeze_model))
                atoms_model_LH.set_constraint(FixAtoms(indices=freeze_model))

        return atoms_real, atoms_model, atoms_model_LH, added_link_atoms

    # ---------------------------------------------------------------------
    def _calc_model_charge(self) -> int:
        """ Calculate Formal charge of MODEL (+link-H) system with RDkit. """
        atoms = read(self.real_pdb)
        real_coord = atoms.get_positions()
        _, _, atoms_model_LH, _ = self._prep_3_layer_atoms(real_coord)
        xyz_lines = [f"{len(atoms_model_LH)}", "MODEL with link-H"]
        for a in atoms_model_LH:
            x, y, z = a.position
            xyz_lines.append(f"{a.symbol:2s} {x:15.8f} {y:15.8f} {z:15.8f}")
        mol = Chem.MolFromXYZBlock("\n".join(xyz_lines) + "\n")
        rdDetermineBonds.DetermineBonds(mol)
        return sum(at.GetFormalCharge() for at in mol.GetAtoms())

    @staticmethod
    def symmetrize_4idx(H: torch.Tensor) -> torch.Tensor:
        n  = H.shape[0]
        H = H.reshape(3 * n, 3 * n)
        H = 0.5 * (H + H.T)
        return H.view(n, 3, n, 3).requires_grad_(False)

    # ==================================================================
    #                        MAIN API
    # ==================================================================
    def compute(
        self,
        coord_ang: np.ndarray,
        *,
        return_forces: bool = False,
        return_hessian: bool = False,
        return_charges: bool = False,
    ) -> Dict:
        """
        Energy / forces / Hessian evaluation.

        Parameters
        ----------
        coord_ang : (N,3) ndarray  Cartesian coordinates in Å (REAL ordering).
        return_*  : bool           Toggles.

        Returns
        -------
        dict
            keys present according to request flags.
        """

        rev_index = {idx: pos for pos, idx in enumerate(self.selection_indices)}

        # 1. build Atoms objects
        # ---------------------------------------------------------------------
        atoms_real, atoms_model, atoms_model_LH, added_link_atoms = self._prep_3_layer_atoms(
            coord_ang
        )

        # 2. high-level evaluation
        # ---------------------------------------------------------------------
        if self.backend == "uma":
            self.predictor.model.eval()

            # prepare UMA batch
            batch_h = self._ase_to_batch(atoms_model_LH)
            n_ml = len(atoms_model_LH)
            pos = batch_h.pos.detach().clone().requires_grad_(True)
            batch_h.pos = pos

            res_pred = self.predictor.predict(batch_h)
            E_tensor = res_pred["energy"].squeeze()
            E_model_high = E_tensor.item()
            F_model_high = res_pred["forces"].detach().cpu().numpy()

            if return_hessian:
                def energy_fn(flat_pos: torch.Tensor):
                    batch_h.pos = flat_pos.view(-1, 3)
                    return self.predictor.predict(batch_h)["energy"].squeeze()
                self.predictor.model.train()
                H_flat = torch.autograd.functional.hessian(energy_fn, pos.view(-1))
                H_high = H_flat.view(n_ml, 3, n_ml, 3).to(self.H_dtype).detach()
                self.predictor.model.eval()
            else:
                H_high = None
            results_h = {}
        else:  # AIMNet2
            elem_ml = [a.symbol for a in atoms_model_LH]
            coord_ml = atoms_model_LH.get_positions()
            data_h = self._prepare_input(elem_ml, coord_ml)
            results_h = self.calc_model_high.eval(
                data_h, forces=True, hessian=return_hessian
            )
            E_model_high = results_h["energy"].item()
            F_model_high = (
                results_h["forces"].to(torch.double).detach().cpu().numpy()
            )
            H_high = (
                results_h["hessian"].to(self.H_dtype)
                if return_hessian and "hessian" in results_h
                else None
            )

        # 3. low-level OpenMM energies
        # ---------------------------------------------------------------------
        atoms_real.set_pbc(False)
        atoms_model.set_pbc(False)
        atoms_model_LH.set_pbc(False)
        atoms_real.calc = self.calc_real_low
        atoms_model.calc = self.calc_model_low

        E_real_low = atoms_real.get_potential_energy()
        F_real_low = np.double(atoms_real.get_forces())

        E_model_low = atoms_model.get_potential_energy()
        F_model_low = np.double(atoms_model.get_forces())

        # 4. combine energies
        # ---------------------------------------------------------------------
        total_E = E_real_low + E_model_high - E_model_low
        results: Dict = {"energy": total_E}

        # 5. combine forces
        # ---------------------------------------------------------------------
        if return_forces or return_hessian:
            F_combined = np.copy(F_real_low)
            for i, ridx in enumerate(self.selection_indices):
                F_combined[ridx] += F_model_high[i] - F_model_low[i]

            for link_idx, ml_idx, mm_idx, dist in added_link_atoms:
                ml_model_idx = rev_index[ml_idx]
                r_ml = atoms_model_LH[ml_model_idx].position
                r_mm = atoms_real[mm_idx].position
                grad_link = F_model_high[link_idx]
                vec = r_mm - r_ml
                R = np.linalg.norm(vec)
                u = vec / R
                I = np.eye(3)
                du_dQ = (I - np.outer(u, u)) / R
                dR_dQ = I - dist * du_dQ
                dR_dM = dist * du_dQ
                J = np.hstack([dR_dQ, dR_dM]).T
                redistributed = J @ grad_link
                F_combined[ml_idx] += redistributed[:3]
                F_combined[mm_idx] += redistributed[3:]
            results["forces"] = F_combined


        # 6. Hessian assembly
        # ---------------------------------------------------------------------
        if return_hessian:
            n_real = len(atoms_real)

            # MM finite-difference Hessians
            # ---------------------------------------------------------------------
            if self.vib_run:
                H_tot = hessian_calc(
                    atoms_real, self.calc_real_low, delta=0.01,
                    info_path=os.path.join(self.vib_dir, "real.log"),
                    dtype=self.H_np_dtype,
                )
                H_tot = torch.from_numpy(H_tot).to(self.ml_device).to(self.H_dtype).reshape(n_real, 3, n_real, 3)

                H_model = hessian_calc(
                    atoms_model, self.calc_model_low, delta=0.01,
                    info_path=os.path.join(self.vib_dir, "model.log"),
                    dtype=self.H_np_dtype,
                )
                H_model = torch.from_numpy(H_model).to(self.ml_device).to(self.H_dtype).reshape(len(atoms_model), 3, len(atoms_model), 3)
            else:
                H_tot = torch.zeros((n_real, 3, n_real, 3),
                                     dtype=self.H_dtype, device=self.ml_device)
                n_ml = len(self.selection_indices)
                H_model = torch.zeros((n_ml, 3, n_ml, 3),
                                      dtype=self.H_dtype, device=self.ml_device)

            n_ml = len(self.selection_indices)

            # ML-ML block (ONIOM diagonal)
            # ---------------------------------------------------------------------
            if H_high is not None:
                for i in range(n_ml):
                    gi = self.selection_indices[i]
                    for j in range(n_ml):
                        gj = self.selection_indices[j]
                        H_tot[gi, :, gj, :].add_(H_high[i, :, j, :]).sub_(H_model[i, :, j, :])
            del H_model

            # ---------------------------------------------------------------------
            #  link metadata & Jacobian
            # ---------------------------------------------------------------------
            link_data: List[Tuple[int, int, int, float, torch.Tensor]] = []

            for link_idx, ml_idx, mm_idx, dist in added_link_atoms:
                ml_model_idx = rev_index[ml_idx]

                r_ml = torch.tensor(atoms_model_LH[ml_model_idx].position,
                                    dtype=self.H_dtype, device=self.ml_device)
                r_mm = torch.tensor(atoms_real[mm_idx].position,
                                    dtype=self.H_dtype, device=self.ml_device)
                vec = r_mm - r_ml
                Rlen = torch.norm(vec)
                u = vec / Rlen
                I3 = torch.eye(3, dtype=self.H_dtype, device=self.ml_device)
                du_dQ = (I3 - torch.outer(u, u)) / Rlen
                dR_dQ = I3 - dist * du_dQ
                dR_dM = dist * du_dQ
                K = torch.hstack([dR_dQ, dR_dM])           # (3×6)

                # (link_idx_REAL , ml_idx_REAL , mm_idx_REAL , dist , K)
                link_data.append((link_idx, ml_idx, mm_idx, dist, K))

            # ---------------------------------------------------------------------
            #  (0) self-diagonal Jᵀ·H·J  and 2-nd term
            # ---------------------------------------------------------------------
            F_high_torch = torch.as_tensor(F_model_high,
                                           dtype=self.H_dtype,
                                           device=self.ml_device)

            for link_idx, ml_idx, mm_idx, dist, K in link_data:
                # Jᵀ H J
                # ---------------------------------------------------------------------
                if H_high is not None:
                    H_l = H_high[link_idx, :, link_idx, :]    # (3×3)
                    H_self = K.T @ H_l @ K                   # (6×6)

                    H_tot[ml_idx, :, ml_idx, :].add_(H_self[0:3, 0:3])
                    H_tot[ml_idx, :, mm_idx, :].add_(H_self[0:3, 3:6])
                    H_tot[mm_idx, :, ml_idx, :].add_(H_self[3:6, 0:3])
                    H_tot[mm_idx, :, mm_idx, :].add_(H_self[3:6, 3:6])

                # 2-nd term  Σ ∂Jᵀ/∂x · f_L
                # ---------------------------------------------------------------------
                f_L = -F_high_torch[link_idx]          # (3,)

                # skip if link force negligible
                if torch.all(torch.abs(f_L) < 1e-12):
                    continue

                r_ml_np = atoms_model_LH[rev_index[ml_idx]].position
                r_mm_np = atoms_real[mm_idx].position
                pos = torch.as_tensor(
                    np.concatenate([r_ml_np, r_mm_np]),
                    dtype=self.H_dtype,
                    device=self.ml_device
                )
                pos.requires_grad_(True)

                def g(p):
                    rq = p[0:3]
                    rm = p[3:6]
                    v = rm - rq
                    uvec = v / torch.norm(v)
                    rL = rq + dist * uvec
                    return torch.dot(f_L, rL)  # scalar

                H_corr6 = torch.autograd.functional.hessian(g, pos).detach()  # (6×6)
                H_corr6 = 0.5 * (H_corr6 + H_corr6.T)  # symmetrize

                H_tot[ml_idx, :, ml_idx, :].add_(H_corr6[0:3, 0:3])
                H_tot[ml_idx, :, mm_idx, :].add_(H_corr6[0:3, 3:6])
                H_tot[mm_idx, :, ml_idx, :].add_(H_corr6[3:6, 0:3])
                H_tot[mm_idx, :, mm_idx, :].add_(H_corr6[3:6, 3:6])

            # ---------------------------------------------------------------------
            #  (a) link–ML off-diagonals  Jᵀ H J
            # ---------------------------------------------------------------------
            if H_high is not None:
                for link_idx, ml_idx, mm_idx, dist, K in link_data:
                    for j_ml, gj in enumerate(self.selection_indices):
                        H_coup = H_high[link_idx, :, j_ml, :]   # 3×3
                        H_row = K.T @ H_coup                    # (6×3)
                        H_col = H_coup.T @ K                    # (3×6)

                        H_tot[ml_idx, :, gj, :].add_(H_row[0:3, :])
                        H_tot[mm_idx, :, gj, :].add_(H_row[3:6, :])
                        H_tot[gj, :, ml_idx, :].add_(H_col[:, 0:3])
                        H_tot[gj, :, mm_idx, :].add_(H_col[:, 3:6])

            # ---------------------------------------------------------------------
            #  (b) link–link couplings  Jᵀ H J
            # ---------------------------------------------------------------------
            if H_high is not None:
                n_links = len(link_data)
                for a in range(n_links):
                    link_idx_a, ml_a, mm_a, dist_a, K_a = link_data[a]
                    for b in range(a + 1, n_links):
                        link_idx_b, ml_b, mm_b, dist_b, K_b = link_data[b]
                        H_ab = H_high[link_idx_a, :, link_idx_b, :]   # 3×3
                        H_tr = K_a.T @ H_ab @ K_b                    # (6×6)

                        H_tot[ml_a, :, ml_b, :].add_(H_tr[0:3, 0:3])
                        H_tot[ml_a, :, mm_b, :].add_(H_tr[0:3, 3:6])
                        H_tot[mm_a, :, ml_b, :].add_(H_tr[3:6, 0:3])
                        H_tot[mm_a, :, mm_b, :].add_(H_tr[3:6, 3:6])

            results["hessian"] = H_tot.detach()
            del H_high, F_high_torch, H_tot; torch.cuda.empty_cache()

        # 7. Charges
        # ---------------------------------------------------------------------
        if return_charges and self.backend == "aimnet2" and "charges" in results_h:
            n_real = len(atoms_real)
            charges = np.zeros(n_real)
            for i, idx in enumerate(self.selection_indices):
                charges[idx] = results_h["charges"][i].item()
            results["charges"] = charges

        return results