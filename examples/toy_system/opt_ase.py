"""Example: geometry optimization using the ASE interface."""

from ase.io import read, write
from ase.optimize import LBFGS
from mlmm_toolkit.mlmm_calc import MLMMASECalculator

atoms = read("structure.pdb")
atoms.calc = MLMMASECalculator(
    input_pdb="complex.pdb",
    real_parm7="complex.parm7",
    model_pdb="ml_region.pdb",
    model_charge=-1,
    model_mult=1,
)

opt = LBFGS(atoms, logfile="opt_ase.log")
opt.run(fmax=0.02, steps=10000)

write("final_ase.pdb", atoms)
print("Done.")
