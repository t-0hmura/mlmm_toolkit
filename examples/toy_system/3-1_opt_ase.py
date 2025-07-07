from ase.io import read, write
from ase.optimize import LBFGS
from mlmm import mlmm_ase          # ASE wrapper

mlmm_kwargs = dict(
    real_pdb     = "complex.pdb",
    real_parm7   = "complex.parm7",
    real_rst7    = "complex.rst7",
    model_pdb    = "ml_region.pdb",
    model_charge = 0,              # Charge of ML region including link atoms
    model_mult   = 1,              # Multiplicity of ML region
    backend      = "uma",          # "uma" or "aimnet2"
    uma_model    = "uma-s-1"
    ml_device    = "auto",         # "auto" | "cuda" | "cpu"
    ml_cuda_idx  = 0,
    mm_device    = "cpu",
    mm_cuda_idx  = 0,
    mm_threads   = 16,
)

atoms = read("structure.pdb")
atoms.calc = mlmm_ase(**mlmm_kwargs)

opt = LBFGS(atoms, logfile="opt.log")
opt.run(fmax=0.01, steps=10000)

write("final.pdb", atoms)
