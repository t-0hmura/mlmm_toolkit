"""ASE Calculator adapter exposing ML/MM functionality via MLMMCore.

This class allows ASE to evaluate energies and forces using the hybrid
machine-learning/molecular-mechanics model.

Author
------
Takuto Ohmura
"""

import numpy as np
from typing import List, Tuple
from ase.calculators.calculator import Calculator, all_changes
from mlmm.mlmm_calc import MLMMCore


class mlmm_ase(Calculator):
    implemented_properties = ['energy', 'forces']

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

                 ml_device: str = "auto",
                 ml_cuda_idx: int = 0,

                 mm_device: str = "cpu",
                 mm_cuda_idx: int = 0,
                 mm_threads: int = 16,
                 freeze_atoms: List[int] | None = None,
                 H_double: bool = False):
        """
        ML/MM calculator for ASE

        Parameters
        ----------
        H_double : bool
            Use double precision for Hessian-related tensors.
        """
        super().__init__()
        self.core = MLMMCore(
                 real_pdb = real_pdb,
                 real_parm7 = real_parm7,
                 real_rst7 = real_rst7,
                 model_pdb = model_pdb,

                 model_charge = model_charge,
                 model_mult = model_mult,
                 link_mlmm = link_mlmm,
                 dist_link = dist_link,
                 backend = backend,

                 vib_run = False,

                 ml_device = ml_device,
                 ml_cuda_idx = ml_cuda_idx,

                 mm_device = mm_device,
                 mm_cuda_idx = mm_cuda_idx,
                 mm_threads = mm_threads,
                 freeze_atoms = freeze_atoms,
                 H_double = H_double)

    # ---------------------------------------------------------------------
    def calculate(self, atoms, properties, system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        
        coord_ang = atoms.get_positions()

        res = self.core.compute(coord_ang, return_forces=True)

        self.results['energy'] = res['energy']
        self.results['forces'] = res['forces']
