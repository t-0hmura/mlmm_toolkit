import numpy as np 
import h5py

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import eigval_to_wavenumber


def save_hessian(h5_fn, geom, cart_hessian=None, energy=None, mult=None):
    raise NotImplementedError(
        "Dense Hessian HDF5 export is disabled. Use the workflow's identified "
        "Hessian file or an explicit NumPy export only when needed."
    )


def save_third_deriv(h5_fn, geom, third_deriv_result, H_mw, H_proj):
    raise NotImplementedError(
        "Dense third-derivative HDF5 export is disabled."
    )


def geom_from_hessian(h5_fn, **geom_kwargs):
    with h5py.File(h5_fn, "r") as handle:
        try:
            atoms = [atom.capitalize() for atom in handle.attrs["atoms"]]
        except (KeyError, TypeError):
            atoms = [atom.capitalize() for atom in handle["atoms"][:]]
        coords3d = handle["coords3d"][:]
        energy = handle.attrs["energy"]
        cart_hessian = handle["hessian"][:]

    geom = Geometry(atoms=atoms, coords=coords3d, **geom_kwargs)
    geom.cart_hessian = cart_hessian
    geom.energy = energy
    return geom
