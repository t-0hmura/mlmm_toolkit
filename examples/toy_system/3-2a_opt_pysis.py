from pysisyphus.io.pdb import geom_from_pdb
from pysisyphus.optimizers.LBFGS import LBFGS
from mlmm import mlmm            # PySiSyphus calculator

mlmm_kwargs = dict(
    real_pdb     = "complex.pdb",
    real_parm7   = "complex.parm7",
    real_rst7    = "complex.rst7",
    model_pdb    = "ml_region.pdb",
    model_charge = -1,
    model_mult   = 1,
    backend      = "uma",
    uma_model    = "uma-s-1p1",
    ml_device    = "auto",
    ml_cuda_idx  = 0,
    mm_device    = "cpu",
    mm_cuda_idx  = 0,
    mm_threads   = 16,
    mem          = 10000,      # MB – PySiSyphus scratch memory
)

geom = geom_from_pdb("structure.pdb")
geom.set_calculator(mlmm(**mlmm_kwargs))

opt = LBFGS(geom, max_cycles=10000, thresh='gau')
opt.run()

with open("final.xyz", "w") as fp:
    fp.write(geom.as_xyz() + "\n")
