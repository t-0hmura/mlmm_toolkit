from mlmm import MLMMCore

core = MLMMCore(
    real_pdb       = "complex.pdb",    # Full system PDB (protein + substrate + solvent)
    real_parm7     = "complex.parm7",  # Amber topology for the full system
    real_rst7      = "complex.rst7",   # Amber coordinates for the full system
    model_pdb      = "ml_region.pdb",  # ML region only (trimmed PDB)
    model_charge   = 0,               # Formal charge of the ML region including link H atoms
    model_mult     = 1,                # Spin multiplicity of the ML region (used by UMA only)
    link_mlmm      = None,             # default: None, link atom pairs are auto determined.
    dist_link      = 1.09,             # Bond length (Å) between link atom and boundary atom
    backend        = "uma",            # ML backend: "uma" or "aimnet2"
    vib_run        = False,            # Whether to compute numerical Hessian (True = finite difference)
    ml_device      = "auto",           # ML backend device: "auto", "cuda", or "cpu"
    ml_cuda_idx    = 0,                # GPU index for ML backend (if using CUDA)
    mm_device      = "cpu",            # MM backend device: "auto", "cuda", or "cpu"
    mm_cuda_idx    = 0,                # GPU index for MM backend (if using CUDA)
    mm_threads     = 16,               # Number of CPU threads for MM force evaluation
)

from ase.io import read; atoms = read("structure.pdb")
coord_ang = atoms.get_positions()

# Coordinates in Å, shape (N, 3) NumPy array
results = core.compute(coord_ang, return_forces=True, return_hessian=True)

energy   = results["energy"]       # float, eV
forces   = results["forces"]       # ndarray (N, 3), eV Å-1
hessian  = results["hessian"]      # torch.Tensor (3N, 3N), eV Å-2

print(f"Energy: {energy:.6f} eV")
print(f"Max force: {forces.max():.6f} eV/Å")
print(f"Hessian shape: {hessian.shape}")
