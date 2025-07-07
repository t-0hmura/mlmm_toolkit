"""Public API exports for the ML/MM calculator package."""

from .mlmm_calc  import MLMMCore
from .mlmm_pysis import mlmm
from .mlmm_ase   import mlmm_ase
from .utils      import add_elem_info, get_freeze_indices
from .def_ml_region import def_ml_region
from .hessian_calc import hessian_calc, calc_freq_from_hessian, write_vib_traj_xyz
from .cartesian_bond_scan import CartesianBondScan
from .partial_hessian_dimer import PartialHessianDimer
from .energy_summary import analyze_single_structure, analyze_three_structures

__version__ = "0.1.0"

__all__ = [
    "__version__", 
    "MLMMCore",
    "mlmm",
    "mlmm_ase",
    "add_elem_info",
    "get_freeze_indices",
    "def_ml_region",
    "hessian_calc",
    "calc_freq_from_hessian",
    "write_vib_traj_xyz",
    "CartesianBondScan",
    "PartialHessianDimer",
    "analyze_single_structure",
    "analyze_three_structures",
]
