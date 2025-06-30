"""PySisyphus calculator wrapping MLMMCore for hybrid ML/MM computations.

Provides energy, force and Hessian evaluations through the PySisyphus
framework.
"""

import numpy as np
from typing import List, Tuple
import torch
import sys

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AU2EV
from pysisyphus import run

from mlmm.mlmm_calc import MLMMCore

EV2AU = 1 / AU2EV   # eV → Hartree

# Hessian dtype
H_dtype = torch.float32

class mlmm(Calculator):
    implemented_properties = ['energy', 'forces', 'hessian', 'charges']

    def __init__(self,
                 real_pdb: str = None,
                 real_parm7: str = None,
                 real_rst7: str = None,
                 model_pdb: str = None,
                 *,
                 model_charge: int = None,
                 model_mult: int = 1,
                 link_mlmm: List[Tuple[str, str]] | None = None,
                 dist_link: float = 1.09,
                 backend: str = "uma",  # "uma" or "aimnet2"

                 vib_run: bool = False,
                 vib_dir: str = None,
                 out_hess_torch: bool = False,
                 H_double: bool = False,

                 ml_device: str = "auto",
                 ml_cuda_idx: int = 0,

                 mm_device: str = "cpu",
                 mm_cuda_idx: int = 0,
                 mm_threads: int = 16,
                 freeze_atoms: List[int] | None = None,
                 **kwargs):
        """
        ML/MM calculator for Pysisyphus
        Args:
            real_pdb (str): Path to the PDB file for the real system.
            real_parm7 (str): Path to the Amber .parm7 file for the real system.
            real_rst7 (str): Path to the Amber .rst7 file for the real system.
            model_pdb (str): Path to the PDB file for the model system.

            model_charge (int): Charge of the model system. If None, charge of model system is calculated automatically with RDKit.
            model_mult (int): Multiplicity of the model system. Default is 1 (singlet).
            link_mlmm (List[Tuple[str, str]] | None): List of tuples specifying the link atoms between ML and MM regions. e.g.) [("CB  ARG   294", "CA  ARG   294")]. If None, link atoms are determined automatically based on distance and element type.
            dist_link (float): Distance for link atoms.
            backend (str): ML backend to use. Options are "uma" or "aimnet2".

            vib_run (bool): Whether to run vibrational analysis.
            out_hess_torch (bool): Whether to output Hessian in torch format. True: torch.Tensor (N,3,N,3) on device, False: numpy.ndarray (N*3,N*3) on cpu.
            H_double (bool): Use double precision for Hessian-related tensors.

            ml_device (str): Device to use for ML calculations. Options are "auto", "cpu", or "cuda".
            ml_cuda_idx (int): CUDA device index if using GPU.

            mm_device (str): Device to use for calculations. Options are "auto", "cpu", or "cuda".
            mm_cuda_idx (int): CUDA device index if using GPU.
            mm_threads (int): Number of threads to use for CPU calculations. {7950X3D/4.20GHz 16 threads faster than RTX 3090 (2x faster for 8->16, Almos same for 16->32)}
            freeze_atoms (List[int] | None): 0-based indices of atoms to freeze during MM Hessian calculations.
        """
        self.model_charge_init = model_charge if model_charge is not None else 0
        self._freeze_atoms = [] if freeze_atoms is None else list(freeze_atoms)
        super().__init__(charge=self.model_charge_init, mult=1, **kwargs)
        self.core = MLMMCore(
                 real_pdb = real_pdb,
                 real_parm7 = real_parm7,
                 real_rst7 = real_rst7,
                 model_pdb = model_pdb,

                 model_charge = model_charge,
                 link_mlmm = link_mlmm,
                 dist_link = dist_link,
                 backend = backend,

                 vib_run = vib_run,
                 vib_dir = vib_dir,
                 
                 ml_device = ml_device,
                 ml_cuda_idx = ml_cuda_idx,

                 mm_device = mm_device,
                 mm_cuda_idx = mm_cuda_idx,
                 mm_threads = mm_threads,
                 freeze_atoms = self._freeze_atoms,
                 H_double = H_double)

        self.out_hess_torch = out_hess_torch
        self.hess_torch_double = H_double
        self.freeze_atoms = freeze_atoms

    # ------------------------------------------------------------------
    @property
    def freeze_atoms(self) -> List[int] | None:
        """Indices of frozen atoms."""
        return self.core.freeze_atoms

    @freeze_atoms.setter
    def freeze_atoms(self, indices: List[int] | None):
        self._freeze_atoms = [] if indices is None else list(indices)
        self.core.freeze_atoms = self._freeze_atoms

    @staticmethod
    def _results_get_energy(results):
        return results['energy'] * EV2AU
    
    @staticmethod
    def _results_get_forces(results):
        return (results['forces'] * (EV2AU / ANG2BOHR)).flatten()
    
    def _results_get_hessian(self, results):
        scale = EV2AU / ANG2BOHR / ANG2BOHR

        H = results.pop("hessian")

        if self.out_hess_torch:
            target_dtype = torch.float64 if self.hess_torch_double else H_dtype
            H = H.view(H.size(0)*3, H.size(2)*3).to(target_dtype)
            H.mul_(scale)
            H = H.detach().requires_grad_(False)
        else:
            H = H.view(H.size(0)*3, H.size(2)*3)
            H.mul_(scale)
            H = H.detach().cpu().numpy()

        del results; torch.cuda.empty_cache()
        return H

    # API for Pysisyphus
    # ---------------------------------------------------------------------
    def get_energy(self, elem, coords):
        coord_ang = np.asarray(coords).reshape(-1, 3) * BOHR2ANG
        res = self.core.compute(coord_ang, return_forces=False, return_hessian=False)
        energy = self._results_get_energy(res)
        del res
        return dict(energy=energy)

    def get_forces(self, elem, coords):
        coord_ang = np.asarray(coords).reshape(-1, 3) * BOHR2ANG
        res = self.core.compute(coord_ang, return_forces=True, return_hessian=False)
        energy = self._results_get_energy(res)
        forces = self._results_get_forces(res)
        del res
        return dict(energy=energy, forces=forces)

    def get_hessian(self, elem, coords):
        coord_ang = np.asarray(coords).reshape(-1, 3) * BOHR2ANG
        res = self.core.compute(coord_ang, return_forces=True, return_hessian=True)
        energy = self._results_get_energy(res)
        forces = self._results_get_forces(res)
        hessian = self._results_get_hessian(res)
        del res; torch.cuda.empty_cache()
        return dict(energy=energy, forces=forces, hessian=hessian)

def run_pysis():
    run.CALC_DICT['mlmm'] = mlmm
    run.run()
