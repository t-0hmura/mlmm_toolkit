#!/usr/bin/env python3
"""Required native MACE-POLAR embedding lane with S/M/L real models."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import torch
from ase import Atoms
from mace.calculators.foundations_models import mace_polar

from mlmm.backends.mace_embedcharge import MACEPolarEmbedding


def core_gate():
    """A covalent ML/MM cut plus two distant ions exercises link-H assembly."""
    import parmed as pmd
    from mlmm.backends.mlmm_calc import MLMMCore
    structure = pmd.Structure()
    coords = np.array([[0.,0.,0.],[-.36,1.03,0.],[-.36,-.51,.89],[-.36,-.51,-.89],
                       [1.53,0.,0.],[1.89,1.03,0.],[1.89,-.51,.89],[1.89,-.51,-.89],
                       [4.,2.,1.],[-3.,-2.,.5]])
    for i in range(10):
        z,mass,charge = (6,12.011,.1) if i in (0,4) else (1,1.008,-.1/3)
        if i >= 8:
            z,mass,charge = 17,35.45,(.35 if i==8 else -.4)
        atom_type = pmd.AtomType("T"+str(i),None,mass,z)
        atom_type.set_lj_params(.01,1.2)
        atom = pmd.Atom(name="A"+str(i),type=atom_type.name,charge=charge,mass=mass,atomic_number=z)
        atom.atom_type = atom_type
        structure.add_atom(atom,"MLR" if i<4 else "MMR",1 if i<4 else 2)
    bond_type = pmd.BondType(200.,1.09)
    carbon_bond = pmd.BondType(200.,1.53)
    structure.bond_types.extend([bond_type,carbon_bond])
    for a,b in [(0,1),(0,2),(0,3),(4,5),(4,6),(4,7),(0,4)]:
        structure.bonds.append(pmd.Bond(structure.atoms[a],structure.atoms[b],
                                        type=carbon_bond if (a,b)==(0,4) else bond_type))
    structure.coordinates = coords
    pdb,parm,model = Path("native_toy.pdb"),Path("native_toy.parm7"),Path("native_model.pdb")
    structure.save(str(pdb),overwrite=True)
    pmd.amber.AmberParm.from_structure(structure).save(str(parm),overwrite=True)
    structure["@1,2,3,4"].save(str(model),overwrite=True)
    core = MLMMCore(input_pdb=str(pdb),real_parm7=str(parm),model_pdb=str(model),
                     model_charge=0,model_mult=1,backend="mace",mace_model="mace-polar-1-m",
                     mace_dtype="float64",ml_device="cuda",mm_device="cpu",mm_threads=1,
                     mm_fd=False,hessian_calc_mode="FiniteDifference",embedcharge=True,
                     xtb_cmd="xtb-must-not-be-called",print_timing=False,print_vram=False,
                     use_bfactor_layers=False,return_partial_hessian=False)
    try:
        assert core.embedding_method == "mace-polar"
        assert len(core.mlmm_links) == 1, core.mlmm_links
        embedded = core.compute(coords,return_forces=True)
        atoms_real,_,atoms_ml,_,_ = core._prep_3_layer_atoms(coords)
        indices = list(core.embedding_point_charge_indices(atoms_real))
        assert indices == [8,9], indices
        de = core._embed_correction.compute_correction(atoms_ml.get_chemical_symbols(),
             atoms_ml.positions,coords[indices],core._get_mm_charges(indices),0,1)[0]
        classical = core.embedding_coulomb(atoms_real,indices,need_forces=False)[0]
        core.embedcharge = False
        gas = core.compute(coords)["energy"]
        core.embedcharge = True
        np.testing.assert_allclose(embedded["energy"]-gas,de-classical,atol=1e-8)
        # Check ML, link host and MM-charge derivatives on the assembled PES.
        errors = []
        for component in [0,1,12,13,24,25,27,28]:
            plus,minus = coords.copy(),coords.copy()
            plus.flat[component] += 1e-4
            minus.flat[component] -= 1e-4
            fd = -(core.compute(plus)["energy"]-core.compute(minus)["energy"])/2e-4
            errors.append(abs(fd-embedded["forces"].flat[component]))
        if max(errors) > 1e-4:
            raise AssertionError(f"ONIOM/link/MM force mismatch: {max(errors)}")
        result = core.compute(coords,return_forces=True,return_hessian=True)
        hessian = result["hessian"]
        if isinstance(hessian,torch.Tensor):
            hessian = hessian.detach().cpu().numpy()
        if not np.isfinite(hessian).all():
            raise AssertionError("Non-finite assembled embedding Hessian")
        return {"method":core.embedding_method,"link_atoms":len(core.mlmm_links),
                 "force_error_ev_angstrom":max(errors),"hessian_shape":list(hessian.shape)}
    finally:
        core.cleanup()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path("mace_polar_embedding_gate.json"))
    args = parser.parse_args()
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is required for the MACE-POLAR embedding gate")
    symbols = ["O","H","H"]
    ml = np.array([[0.,0.,0.],[.95,.03,0.],[-.25,.91,.1]])
    mm = np.array([[3.,1.,.4],[-2.8,2.1,-.5]])
    charges = np.array([.35,-.4])
    records = []
    for size in ("s","m","l"):
        calc = mace_polar(model="polar-1-"+size,device="cuda",default_dtype="float64")
        adapter = MACEPolarEmbedding(calc)
        atoms = Atoms(symbols,positions=ml)
        atoms.info.update(charge=0,spin=1)
        atoms.calc = calc
        gas = atoms.get_potential_energy()
        x = np.concatenate((ml,mm))
        def evaluate(coords, q=charges, **kw):
            return adapter.compute_correction(symbols,coords[:3],coords[3:],q,0,1,**kw)
        e,f,_ = evaluate(x,need_forces=True)
        fd = np.zeros_like(x)
        step = 1e-4
        for i in range(x.size):
            plus,minus = x.copy(),x.copy()
            plus.flat[i] += step
            minus.flat[i] -= step
            fd.flat[i] = -(evaluate(plus)[0]-evaluate(minus)[0])/(2*step)
        error = float(np.max(np.abs(fd-f)))
        if error > 1e-5:
            raise AssertionError(f"ML/MM energy-force mismatch: {error}")
        np.testing.assert_allclose(f.sum(0),0,atol=1e-8)
        zero_e,zero_f,_ = evaluate(x,q=np.zeros(2),need_forces=True)
        np.testing.assert_allclose(zero_e,0,atol=1e-10)
        np.testing.assert_allclose(zero_f,0,atol=1e-9)
        np.testing.assert_allclose(evaluate(x+100)[0],e,atol=1e-9)
        atoms.calc.reset()
        np.testing.assert_allclose(atoms.get_potential_energy(),gas,atol=1e-10)
        record = {"model":"polar-1-"+size,"delta_energy_ev":e,"force_error_ev_angstrom":error}
        if size == "m":
            _,_,h = evaluate(x,need_hessian=True)
            adapter.hessian_step *= .5
            _,_,h_half = evaluate(x,need_hessian=True)
            error_h = float(np.max(np.abs(h-h_half)))
            if error_h > 1e-3:
                raise AssertionError(f"Mixed Hessian step sensitivity: {error_h}")
            if not np.any(np.abs(h[:9,9:]) > 1e-6):
                raise AssertionError("Missing ML-MM mixed Hessian block")
            partial = evaluate(x,need_hessian=True,hessian_mm_rows=[1])[2]
            dofs = list(range(9))+[12,13,14]
            np.testing.assert_allclose(partial,h_half[np.ix_(dofs,dofs)],atol=1e-9)
            record["hessian_step_error_ev_angstrom2"] = error_h
        records.append(record)
    core = core_gate()
    args.output.write_text(json.dumps({"status":"passed","models":records,"oniom_link_test":core},indent=2)+"\n")
    print("PASS MACE-POLAR S/M/L native ML/MM derivatives, zero-charge limit and mixed Hessian")


if __name__ == "__main__":
    main()
