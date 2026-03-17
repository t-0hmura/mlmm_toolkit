"""Example: direct MLMMCore API usage (energy, forces, Hessian)."""

from ase.io import read
from mlmm.mlmm_calc import MLMMCore

core = MLMMCore(
    input_pdb="r_complex.pdb",
    real_parm7="p_complex.parm7",
    model_pdb="ml_region_r.pdb",
    model_charge=-1,
    model_mult=1,
)

atoms = read("r_complex_layered.pdb")
coord_ang = atoms.get_positions()

results = core.compute(coord_ang, return_forces=True, return_hessian=True)

energy = results["energy"]
forces = results["forces"]
hessian = results["hessian"]

print(f"Energy:        {energy:.6f} eV")
print(f"Max |force|:   {forces.max():.6f} eV/Ang")
print(f"Hessian shape: {hessian.shape}")
