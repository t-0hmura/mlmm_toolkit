"""Example: geometry optimization using the ASE interface."""

from ase.io import read, write
from ase.optimize import LBFGS
from mlmm_toolkit.mlmm_calc import MLMMCore, MLMMASECalculator

core = MLMMCore(
    input_pdb="r_complex.pdb",
    real_parm7="p_complex.parm7",
    model_pdb="ml_region_r.pdb",
    model_charge=-1,
    model_mult=1,
)

atoms = read("r_complex_layered.pdb")
atoms.calc = MLMMASECalculator(core)

opt = LBFGS(atoms, logfile="opt_ase.log")
opt.run(fmax=0.02, steps=10000)

write("final_ase.pdb", atoms)
print("Done.")
