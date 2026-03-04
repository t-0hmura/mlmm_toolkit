"""Example: geometry optimization using the PySisyphus interface."""

from pysisyphus.io.pdb import geom_from_pdb
from pysisyphus.optimizers.LBFGS import LBFGS
from mlmm_toolkit.mlmm_calc import mlmm

geom = geom_from_pdb("structure.pdb")
geom.set_calculator(mlmm(
    input_pdb="complex.pdb",
    real_parm7="complex.parm7",
    model_pdb="ml_region.pdb",
    model_charge=-1,
    model_mult=1,
))

opt = LBFGS(geom, max_cycles=10000, thresh="gau")
opt.run()

with open("final_pysis.xyz", "w") as fp:
    fp.write(geom.as_xyz() + "\n")
print("Done.")
