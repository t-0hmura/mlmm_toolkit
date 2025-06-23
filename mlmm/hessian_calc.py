"""Hessian utilities and vibrational analysis helpers for the ML/MM workflow.

These routines construct Hessians, compute frequencies and handle I/O.

Author
------
Takuto Ohmura
"""

import os
import numpy as np
import time
import torch
from typing import Sequence, Union
import ase.units as units
from ase.data import atomic_masses
from ase.constraints import FixAtoms
from ase.atoms import Atoms
from ase.io import write
from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AU2EV, AMU2AU


# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------
def clone_atoms(atoms):
    """Return a shallow copy of *atoms* keeping constraints (no PBC)."""
    new = Atoms(
        symbols=atoms.get_chemical_symbols(),
        positions=atoms.get_positions(),
        pbc=False,
    )

    if atoms.constraints:
        fix_idx = [
            i for c in atoms.constraints if isinstance(c, FixAtoms)
            for i in c.get_indices()
        ]
        if fix_idx:
            new.set_constraint(FixAtoms(indices=fix_idx))

    return new


# ---------------------------------------------------------------------
# Finite-difference Hessian
# ---------------------------------------------------------------------
def hessian_calc(
    atoms,
    calc,
    delta: float = 0.01,
    info_path: str | None = None,
    *,
    dtype: np.dtype = np.float64,
):
    """
    Numerically compute the full Cartesian Hessian (second derivatives).

    Finite differences are taken by displacing each *unfixed* atom by
    ``±delta`` Å along x, y, z.  Central differences give the sub-Hessian,
    which is **always expanded to the full (3N, 3N) matrix** by padding
    zeros for fixed atoms before returning.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure to differentiate.
    calc : ase.Calculator
        Calculator providing forces.
    delta : float, default 0.01 Å
        Displacement size.
    info_path : str | None
        If given, progress info is appended to this file.
    dtype : numpy dtype, default float64
        Data type of the returned Hessian.

    Returns
    -------
    ndarray
        Cartesian Hessian of shape ``(3N, 3N)`` where
        ``N == len(atoms)`` (fixed-atom rows/cols are zero).
    """
    # movable / fixed atom lists
    # ---------------------------------------------------------------------
    fixed = {
        i for c in atoms.constraints if isinstance(c, FixAtoms)
        for i in c.get_indices()
    }
    movable = np.asarray([i for i in range(len(atoms)) if i not in fixed])
    m = len(movable)
    n_dof = 3 * m

    # Short-circuit: all atoms frozen → return zeros
    if m == 0:
        N = len(atoms)
        return np.zeros((3 * N, 3 * N))

    H_sub = np.empty((n_dof, n_dof), dtype=dtype)
    row = 0

    # progress logging
    # ---------------------------------------------------------------------
    if info_path is not None:
        os.makedirs(os.path.dirname(info_path), exist_ok=True)
        log = open(info_path, "w", encoding="utf-8")
        log.write("Hessian calculation by numerical differentiation\n")
        log.write("------------------------------------------------\n")
        t0 = time.time()
        checkpoints = [int(m * k / 10) for k in range(1, 11)]  # every 10 %

    # finite-difference loop
    # ---------------------------------------------------------------------
    for count, a in enumerate(movable, start=1):
        for coord in range(3):
            plus = clone_atoms(atoms)
            plus.calc = calc
            plus.positions[a, coord] += delta
            Fp = plus.calc.get_forces(plus)

            minus = clone_atoms(atoms)
            minus.calc = calc
            minus.positions[a, coord] -= delta
            Fm = minus.calc.get_forces(minus)

            # Central difference:  −∂F/∂x = ∂²E/∂x²
            H_sub[row] = (Fm - Fp)[movable].ravel() / (2 * delta)
            row += 1

        # logging
        # ---------------------------------------------------------------------
        if info_path is not None and count in checkpoints:
            pct = checkpoints.index(count) * 10 + 10
            elapsed = time.time() - t0
            total_est = elapsed * 100.0 / pct
            remaining = total_est - elapsed
            me, se = divmod(int(elapsed), 60)
            mr, sr = divmod(int(remaining), 60)
            log.write(
                f"{pct:3d}% ({count}/{m})   Elapsed {me}m{se}s   ETA {mr}m{sr}s\n"
            )
            log.flush()

    if info_path is not None:
        elapsed = time.time() - t0
        me, se = divmod(int(elapsed), 60)
        log.write(f"Done in {me}m{se}s\n")
        log.close()

    # assemble full (3N, 3N) Hessian
    # ---------------------------------------------------------------------
    N = len(atoms)
    H = np.zeros((3 * N, 3 * N), dtype=dtype)
    for i_local, i_atom in enumerate(movable):
        for j_local, j_atom in enumerate(movable):
            H[
                3 * i_atom : 3 * i_atom + 3,
                3 * j_atom : 3 * j_atom + 3,
            ] = H_sub[
                3 * i_local : 3 * i_local + 3,
                3 * j_local : 3 * j_local + 3,
            ]
    return H


# ---------------------------------------------------------------------
# Vibrational analysis utilities
# ---------------------------------------------------------------------
def _get_masses(Z):
    """Return atomic masses (amu) for a list/array of atomic numbers."""
    return np.array([atomic_masses[z] for z in Z])


def _build_tr_basis(coords_bohr: torch.Tensor, masses_amu: torch.Tensor) -> torch.Tensor:
    """Return (3N,6) mass-weighted basis for Tx,Ty,Tz,Rx,Ry,Rz."""
    device, dtype = coords_bohr.device, coords_bohr.dtype
    N = coords_bohr.shape[0]

    m_au = masses_amu.to(dtype=dtype, device=device) * AMU2AU
    m_sqrt = torch.sqrt(m_au).reshape(-1, 1)

    com = (m_au.reshape(-1, 1) * coords_bohr).sum(0) / m_au.sum()
    x = coords_bohr - com

    eye3 = torch.eye(3, dtype=dtype, device=device)
    cols = []

    for i in range(3):                                   # Tx, Ty, Tz
        cols.append((eye3[i].repeat(N, 1) * m_sqrt).reshape(-1, 1))
    for i in range(3):                                   # Rx, Ry, Rz
        rot = torch.cross(x, eye3[i].expand_as(x), dim=1) * m_sqrt
        cols.append(rot.reshape(-1, 1))

    return torch.cat(cols, dim=1)


def calc_freq_from_hessian(
    H: torch.Tensor,
    elem: Sequence[int],
    tol: float = 1e-6,
    project_tr: bool = True,
    coords_bohr: np.ndarray | None = None,
    masses_amu: np.ndarray | None = None,
    verbose: bool = True,
):
    """
    Compute vibrational frequencies from a (possibly sparse) Cartesian Hessian.

    Parameters
    ----------
    H : torch.Tensor
        Full Cartesian Hessian **in Hartree / Bohr²**.  Rows / cols that are
        *exactly* zero (e.g. fixed DOF padded with zeros) are automatically
        dropped.
    elem : Sequence[int]
        Atomic numbers (length N, unchanged from the original system).
    tol : float
        Threshold for discriminating (near-)zero eigenvalues.
    project_tr : bool
        If True, translational / rotational modes are projected out *for the
        active atoms only (frozen atoms are not included in the point-mass
        projection).
    coords_bohr, masses_amu : ndarray | None
        Only required when `project_tr` is True.
        **coords_bohr must be shape (N,3) for the *original* system** — the
        function takes the subset internally, so you can pass it as-is.

    Returns
    -------
    freqs_cm : ndarray
        Vibrational frequencies (cm⁻¹). Imaginary modes are negative.
    hnu_eV : ndarray
        Zero-point energies hν (eV) for the same modes.
    modes_t : torch.Tensor
        Mass-weighted eigenvectors, shape ``(nmode, 3N)`` kept on the
        same device as ``H``.  Frozen DOF are zero.
    """

    # ---------------------------------------------------------------------
    # 0)   Auto-detect active DOF and reduce if necessary
    # ---------------------------------------------------------------------
    #    Determined by row norm == 0 (columns will also be zero)
    #    e.g. when all three components are zero in freeze_atoms
    # ---------------------------------------------------------------------
    row_nonzero = torch.linalg.norm(H, dim=1) > tol
    active_dof  = torch.nonzero(row_nonzero, as_tuple=False).squeeze()

    n_tot = H.shape[0]           # = 3N
    if active_dof.numel() == 0:
        raise RuntimeError("All DOF are zero – no vibrational analysis possible.")

    reduced = active_dof.numel() != n_tot
    if reduced:
        H = H[active_dof][:, active_dof]         # (3N_act, 3N_act)

        # Map DOF to atom ID (0-based)
        atom_map_full = (active_dof // 3).cpu().numpy()
        active_atoms  = np.unique(atom_map_full)
        elem_act      = np.asarray(elem)[active_atoms]
        if coords_bohr is not None:
            coords_bohr_act = coords_bohr[active_atoms]
        if masses_amu is not None:
            masses_amu_act = masses_amu[active_atoms]
    else:
        elem_act        = elem
        coords_bohr_act = coords_bohr
        masses_amu_act  = masses_amu

    # ---------------------------------------------------------------------
    # 1)   TR projection (optional)
    # ---------------------------------------------------------------------
    if project_tr and coords_bohr_act is not None and masses_amu_act is not None:
        coords = torch.as_tensor(coords_bohr_act, dtype=H.dtype, device=H.device)
        masses = torch.as_tensor(masses_amu_act, dtype=H.dtype, device=H.device)
        B = _build_tr_basis(coords, masses)
        Bt = B.T
        BtB_inv = torch.linalg.inv(Bt @ B)
        with torch.no_grad():
            G = BtB_inv @ (Bt @ H)
            H.sub_(B @ G)
            HB = H @ B
            H.sub_(HB @ BtB_inv @ Bt)

    # ---------------------------------------------------------------------
    # 2)   Mass-weighting and diagonalization
    # ---------------------------------------------------------------------
    m_vec = np.repeat(_get_masses(elem_act), 3)   # amu
    m     = torch.as_tensor(m_vec, dtype=H.dtype, device=H.device)
    inv_sqrt_m = torch.sqrt(1.0 / m)
    H *= inv_sqrt_m[:, None]
    H *= inv_sqrt_m
    H.requires_grad_(False)
    omega2, modes = torch.linalg.eigh(H)
    del H, inv_sqrt_m, m; torch.cuda.empty_cache()

    # ---------------------------------------------------------------------
    # 3)   Remove low frequencies and report eigenvalues
    # ---------------------------------------------------------------------
    if verbose:
        neg_idx = torch.where(omega2 < -tol)[0]
        if neg_idx.numel() == 0:
            print("No negative eigenvalue.")
        elif neg_idx.numel() > 1:
            vals = np.sort(omega2[neg_idx].cpu().numpy())
            print(f"{neg_idx.numel()} negative eigenvalues detected:")
            print(vals)

    sel    = torch.abs(omega2) > tol
    omega2 = omega2[sel]
    modes  = modes[:, sel].T        # (nmode, 3N_act)

    # ---------------------------------------------------------------------
    # 4)   Pad eigenvectors to length 3N
    # ---------------------------------------------------------------------
    if reduced:
        nmode = modes.shape[0]
        full = torch.zeros((nmode, n_tot), dtype=modes.dtype, device=modes.device)
        full[:, active_dof] = modes
        modes = full

    # ---------------------------------------------------------------------
    # 5)   Convert to frequencies and hν
    # ---------------------------------------------------------------------
    #  (Hartree / Bohr² / amu)  →  hν [eV]
    
    #  hν =  ħ ω
    #  ω  = √ω²
    #  ω² = (Hartree / Bohr²) × (J/Hartree) / amu × (amu/kg) / (Bohr²→m²)
    
    #  s_new = ħ·10¹⁰ / √(e·amu) · √(AU2EV) / BOHR2ANG
    #      = 0.6373391776792675  [eV·√(amu·Bohr²/Hartree)]⁻¹
    
    s_new = (units._hbar * 1e10
             / np.sqrt(units._e * units._amu)
             * np.sqrt(AU2EV) / BOHR2ANG)

    hnu = s_new * torch.sqrt(torch.abs(omega2))
    hnu = torch.where(omega2 < 0, -hnu, hnu)   # keep sign

    freqs_cm = (hnu / units.invcm).cpu().numpy()
    hnu_eV   = hnu.cpu().numpy()

    return freqs_cm, hnu_eV, modes


# ---------------------------------------------------------------------
# XYZ animation writer
# ---------------------------------------------------------------------
def write_vib_traj_xyz(
    atoms: Atoms,
    mode: Union[np.ndarray, torch.Tensor, Sequence[float]],
    filename: str,
    amplitude: float = 0.25,
    n_frames: int = 20,
    comment: str = "vibrational mode",
) -> None:
    """Save an XYZ trajectory that animates a single vibrational mode."""
    mode = np.asarray(mode, dtype=float)
    if mode.ndim == 1:
        mode = mode.reshape(-1, 3)
    if mode.shape != (len(atoms), 3):
        raise ValueError(f"mode shape {mode.shape} != ({len(atoms)}, 3)")

    mode /= np.linalg.norm(mode)
    mode *= amplitude

    images = []
    for i in range(n_frames):
        phase = np.sin(2.0 * np.pi * i / n_frames)
        img = atoms.copy()
        img.positions += phase * mode
        img.info["comment"] = f"{comment}  frame={i+1}/{n_frames}"
        images.append(img)

    write(filename, images, format="xyz")
    print(f"Wrote {n_frames} frames to '{filename}'.")