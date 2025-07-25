"""
Partial Hessian Guided Dimer (PHG-DImer) – ML/MM dimer search with partial Hessian & mass-scaled multi-mode flattening
================================================================================================

Workflow
--------
0. **Pre-opt (ML region frozen)** – relaxes starting structure.
1. **Loose-threshold Dimer**  
   * Build *partial* Hessian (ML + MM ≤ cutoff) using dynamic `freeze_atoms`.  
   * Extract most-negative mode → `mode.dat`.  
   * Run Dimer (``thresh_loose``) for `update_interval_hessian` steps.  
   * Recompute partial Hessian & update mode until loose convergence.
2. **Final-threshold Dimer**  
   * Same loop with tighter ``thresh`` → converge to a TS that still
     contains extra imaginary modes.
3. **Full Hessian & flatten loop**  
   * Compute a **full** Hessian (no padding).  
   * If > 1 imaginary mode:  
       • Pick the 2nd most-negative mode, greedily add spatially separated
         modes, displace atoms by  

           Δr = ± ``flatten_amp_ang`` × √(m_C/m_i) × e_i  

         (carbon moves exactly ``flatten_amp_ang`` Å).  
       • After displacement **rerun the Final-threshold Dimer loop** with
         partial Hessian.  
   * Repeat until ≤ 1 imaginary mode or ``flatten_max_iter``.
4. **Export** final geometry and surviving imaginary mode(s).
"""

from __future__ import annotations

# stdlib
# ---------------------------------------------------------------------
import time, tempfile, os
from pathlib import Path
from typing import Optional, List, Tuple, Union

# 3rd-party
# ---------------------------------------------------------------------
import numpy as np
import torch
from ase import Atoms
from ase.io import read, write
from ase.data import atomic_masses

from mlmm import mlmm as MLMM
from .hessian_calc import (
    calc_freq_from_hessian, write_vib_traj_xyz, _build_tr_basis
)
from pysisyphus.io.xyz import geom_from_xyz
from pysisyphus.io.pdb import geom_from_pdb
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.calculators.Dimer import Dimer
from pysisyphus.constants import BOHR2ANG, AMU2AU


class PartialHessianDimer:
# ---------------------------------------------------------------------
    # Constructor
# ---------------------------------------------------------------------
    def __init__(
        # files & geometry
        # ---------------------------------------------------------------------
        self,
        fn: str = "./coord/gs_hei.xyz",
        real_pdb: str = "./parm/complex.pdb",
        freeze_atoms: Optional[list[int]] = None,
        out_dir: str = "./dump/dimer/",
        vib_dir: str = "./dump/vib/",
        # thresholds & cycles
        # ---------------------------------------------------------------------
        thresh_loose: str = "gau_loose",
        thresh: str = "baker",
        update_interval_hessian: int = 50,
        lobpcg: bool = True,
        max_cycles: int = 100_000,
        max_cycles_preopt: Optional[int] = None, # if None, same as `max_cycles`. if 0, no preopt.
        dump: bool = False,
        mem: int = 100_000,
        neg_freq_thresh: float = 10.0,
        partial_mm_cutoff: float = 0.0,
        # flatten loop
        # ---------------------------------------------------------------------
        flatten_amp_ang: float = 0.20,
        flatten_max_iter: int = 20,
        flatten_sep_cutoff: float = 2.0,
        flatten_k: int = 10,
        # MLMM / Dimer
        # ---------------------------------------------------------------------
        mlmm_kwargs: Optional[dict] = None,
        dimer_kwargs: Optional[dict] = None,
    ):
        # basic paths
        # ---------------------------------------------------------------------
        self.fn = fn
        self.real_pdb = real_pdb
        self.out_dir = Path(out_dir)
        self.vib_dir = Path(vib_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.vib_dir.mkdir(parents=True, exist_ok=True)

        # user parameters
        # ---------------------------------------------------------------------
        self.thresh_loose = thresh_loose
        self.thresh = thresh
        self.update_interval_hessian = int(update_interval_hessian)
        self.lobpcg = lobpcg
        self.max_cycles = int(max_cycles)
        self.max_cycles_preopt = (
            self.max_cycles if max_cycles_preopt is None else int(max_cycles_preopt)
        )
        self.dump = dump
        self.mem = int(mem)
        self.neg_freq_thresh = float(neg_freq_thresh)
        self.partial_mm_cutoff = float(partial_mm_cutoff)

        self.flatten_amp_ang = float(flatten_amp_ang)
        self.flatten_max_iter = int(flatten_max_iter)
        self.flatten_sep_cutoff = float(flatten_sep_cutoff)
        self.flatten_k = int(flatten_k)

        # freeze lists
        # ---------------------------------------------------------------------
        self.freeze_atoms_static = [] if freeze_atoms is None else list(freeze_atoms)

        # MLMM kwargs
        # ---------------------------------------------------------------------
        self.mlmm_kwargs = {} if mlmm_kwargs is None else dict(mlmm_kwargs)
        self.dimer_kwargs = {} if dimer_kwargs is None else dict(dimer_kwargs)

        if self.freeze_atoms_static and "freeze_atoms" not in self.mlmm_kwargs:
            self.mlmm_kwargs["freeze_atoms"] = self.freeze_atoms_static.copy()

        self.model_pdb = self.mlmm_kwargs.get("model_pdb")
        if not self.model_pdb:
            raise ValueError("`mlmm_kwargs` must include 'model_pdb'.")

        # dtype / device
        # ---------------------------------------------------------------------
        if self.mlmm_kwargs.get('H_double', False) == True:
            self.H_dtype = torch.float64
        else:
            self.H_dtype = torch.float32

        ml_device = self.mlmm_kwargs.get('ml_device', 'auto')
        ml_cuda_idx = self.mlmm_kwargs.get('ml_cuda_idx', 0)
        if ml_device == "auto":
            ml_device = "cuda" if torch.cuda.is_available() else "cpu"
        self.H_device = torch.device(
            f"cuda:{ml_cuda_idx}" if ml_device == "cuda" else "cpu"
        )

        # workspace
        # ---------------------------------------------------------------------
        self._tmpdir_obj = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self._tmpdir_obj.name)
        self.mode_path = str(self.tmpdir / "mode.dat")

        # geometry
        # ---------------------------------------------------------------------
        if fn.endswith(".xyz"):
            self.geom = geom_from_xyz(fn)
        elif fn.endswith(".pdb"):
            self.geom = geom_from_pdb(fn)
        else:
            raise ValueError("`fn` must be .xyz or .pdb")

        self.geom.freeze_atoms = self.freeze_atoms_static.copy()

        # masses
        # ---------------------------------------------------------------------
        masses_amu = np.array([atomic_masses[z] for z in self.geom.atomic_numbers])
        self.masses_amu_np = masses_amu
        self.masses_au_t = torch.tensor(masses_amu * AMU2AU,
                                        dtype=self.H_dtype, device=self.H_device)
        # for mass-scaled flatten
        self.mass_scale = np.sqrt(12.011 / masses_amu)[:, None]

        # template PDB
        # ---------------------------------------------------------------------
        self._template_real_atoms = read(self.real_pdb)

        # ML region indices (0-based w.r.t real_pdb ordering)
        # -----------------------------------------------------------------
        self.ml_indices = self._extract_ml_indices()

        # clean old traj
        # ---------------------------------------------------------------------
        (self.out_dir / "optimization_all.trj").unlink(missing_ok=True)

    # ================================================================
    # helper – ML region indices
    # ================================================================
    def _extract_ml_indices(self) -> List[int]:
        """Return indices of atoms belonging to the ML region."""
        model_ids: set[str] = set()
        with open(self.model_pdb, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    model_ids.add(
                        f"{line[12:16].strip()} {line[17:20].strip()} {line[22:26].strip()}"
                    )

        ml_idx: list[int] = []
        with open(self.real_pdb, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    idx = int(line[6:11]) - 1
                    atom_id = f"{line[12:16].strip()} {line[17:20].strip()} {line[22:26].strip()}"
                    if atom_id in model_ids:
                        ml_idx.append(idx)
        return ml_idx

    # ================================================================
    # helper – dynamic freeze list
    # ================================================================
    def _compute_dynamic_freeze(self) -> List[int]:
        if self.partial_mm_cutoff <= 0.1:
            return self.freeze_atoms_static.copy()

        atoms = self._template_real_atoms.copy()
        atoms.set_positions(self.geom.coords.reshape(-1, 3) * BOHR2ANG)
        coords = atoms.get_positions()

        coords_t = torch.as_tensor(coords, dtype=self.H_dtype, device=self.H_device)
        ml_coords_t = coords_t[self.ml_indices]

        dmat = torch.cdist(ml_coords_t, coords_t)
        min_dist = torch.min(dmat, dim=0).values

        freeze_mask = min_dist > self.partial_mm_cutoff
        freeze_mask[self.ml_indices] = False

        frozen = set(self.freeze_atoms_static)
        frozen.update(torch.nonzero(freeze_mask).view(-1).cpu().tolist())

        return sorted(frozen)

    # ================================================================
    # helper – TR projection (in-place)
    # ================================================================
    def _project_out_tr(self, H: torch.Tensor, coords_bohr_t: torch.Tensor) -> torch.Tensor:
        if coords_bohr_t.ndim == 1:
            coords_bohr_t = coords_bohr_t.view(-1, 3)
        elif coords_bohr_t.shape[-1] != 3:
            coords_bohr_t = coords_bohr_t.reshape(-1, 3)

        B       = _build_tr_basis(coords_bohr_t, self.masses_au_t)      # (3N,6)
        Bt      = B.T
        BtB_inv = torch.linalg.inv(Bt @ B)                              # (6,6)

        BtH = Bt @ H          # (6,3N) = Bt H₀
        HB  = H  @ B          # (3N,6) = H₀ B

        # --- 1)  H  ←  H₀  -  B (BtB)⁻¹ (Bt H₀) ------------------------------
        H.sub_(B @ (BtB_inv @ BtH))

        # --- 2)  H  ←  H  -  H₀ B (BtB)⁻¹ Bt ---------------------------------
        H.sub_(HB @ BtB_inv @ Bt)

        # --- 3)  H  ←  H  +  B (BtB)⁻¹ (Bt H₀ B) (BtB)⁻¹ Bt ------------------
        H.add_(B @ (BtB_inv @ (BtH @ B)) @ BtB_inv @ Bt)

        del B, Bt, BtB_inv, BtH, HB; torch.cuda.empty_cache()

        return H

    # ================================================================
    # helper – lowest-mode writer
    # ================================================================
    def _write_lowest_mode(self, H: torch.Tensor, freeze_atoms: Optional[List[int]] = None) -> np.ndarray:
        """Diagonalise ``H`` and save the most-negative mode.

        Zero-padded blocks (e.g. frozen atoms) are automatically removed prior to
        translational/rotational projection.  The returned mode is padded back to
        the full length ``3N`` and written to ``mode.dat``.
        """
        # --------------------------------------------------------------
        # 0)   Detect active degrees of freedom and reduce Hessian
        # --------------------------------------------------------------
        tol = 1e-6
        row_nonzero = torch.linalg.norm(H, dim=1) > tol
        active_dof = torch.nonzero(row_nonzero, as_tuple=False).squeeze()

        n_tot = H.shape[0]  # = 3N
        if active_dof.numel() == 0:
            raise RuntimeError("All DOF are zero – no mode available.")

        reduced = active_dof.numel() != n_tot
        if reduced:
            H = H[active_dof][:, active_dof]

            # map to active atoms
            atom_map_full = (active_dof // 3).cpu().numpy()
            active_atoms = np.unique(atom_map_full)
            coords_full = torch.tensor(
                self.geom.coords, dtype=self.H_dtype, device=self.H_device
            ).view(-1, 3)
            coords_t = coords_full[active_atoms]
            masses_t = self.masses_au_t[active_atoms]
        else:
            coords_t = torch.tensor(
                self.geom.coords, dtype=self.H_dtype, device=self.H_device
            ).view(-1, 3)
            masses_t = self.masses_au_t

        # --------------------------------------------------------------
        # 1)   Project out translational/rotational modes (active atoms)
        # --------------------------------------------------------------
        B = _build_tr_basis(coords_t, masses_t)  # (3N_act,6)
        Bt = B.T
        BtB_inv = torch.linalg.inv(Bt @ B)

        BtH = Bt @ H
        HB = H @ B

        H.sub_(B @ (BtB_inv @ BtH))
        H.sub_(HB @ BtB_inv @ Bt)
        H.add_(B @ (BtB_inv @ (BtH @ B)) @ BtB_inv @ Bt)

        del B, Bt, BtB_inv, BtH, HB
        torch.cuda.empty_cache()

        # --------------------------------------------------------------
        # 2)   Diagonalise and extract the lowest eigenvector
        # --------------------------------------------------------------
        if self.lobpcg:
            eigvals, eigvecs = torch.lobpcg(H, k=1, largest=False)
            mode_vec_red = eigvecs[:, 0]
        else:
            eigvals, eigvecs = torch.linalg.eigh(H)
            mode_vec_red = eigvecs[:, torch.argmin(eigvals)]

        # --------------------------------------------------------------
        # 3)   Pad eigenvector back to length 3N
        # --------------------------------------------------------------
        if reduced:
            mode_vec = torch.zeros(n_tot, dtype=mode_vec_red.dtype, device=mode_vec_red.device)
            mode_vec[active_dof] = mode_vec_red
        else:
            mode_vec = mode_vec_red

        # --------------------------------------------------------------
        # 4)   Zero out frozen atoms (vectorised, no Python loop)
        # --------------------------------------------------------------
        if freeze_atoms is None:
            freeze_atoms = self.freeze_atoms_static

        if freeze_atoms:
            freeze_idx = torch.as_tensor(freeze_atoms, dtype=torch.long, device=mode_vec.device)
            flat_idx = (freeze_idx.unsqueeze(1) * 3 + torch.arange(3, device=mode_vec.device)).reshape(-1)
            mode_vec[flat_idx] = 0.0

        # --------------------------------------------------------------
        # 5)   Normalise, reshape, save, return
        # --------------------------------------------------------------
        mode = (mode_vec / torch.linalg.norm(mode_vec)).reshape(-1, 3).cpu().numpy()
        np.savetxt(self.mode_path, mode, fmt="%.10f")

        return mode

    # ================================================================
    # helper – partial Hessian (dynamic freeze)
    # ================================================================
    def _partial_hessian(self) -> Tuple[torch.Tensor, List[int]]:
        freeze = self._compute_dynamic_freeze()
        vib_kw = self.mlmm_kwargs.copy()
        vib_kw["freeze_atoms"] = freeze

        if self.partial_mm_cutoff <= 0.1:
            calc = MLMM(out_hess_torch=True, vib_run=False, **vib_kw)
        else:
            calc = MLMM(out_hess_torch=True, vib_run=True, **vib_kw)

        H = calc.get_hessian(self.geom.atom_types, self.geom.coords)["hessian"].to(dtype=self.H_dtype, device=self.H_device)
        del calc; torch.cuda.empty_cache()
        return H, freeze

    # ================================================================
    # helper – full Hessian (zero out freeze_atoms)
    # ================================================================
    def _full_hessian(self) -> Tuple[np.ndarray, torch.Tensor]:
        """Calculate the full Hessian and perform vibrational analysis after
        zeroing all degrees of freedom corresponding to ``freeze_atoms`` (the
        3×3 blocks and their interactions with other atoms).

        Returns
        -------
        freqs : np.ndarray
            Frequencies in cm⁻¹.
        modes : torch.Tensor
            Mass-weighted eigenvectors with shape ``(nmode, 3N)`` on GPU.
        """
        # 1. obtain the Hessian
        # ---------------------------------------------------------------------
        calc = MLMM(out_hess_torch=True, vib_run=True, **self.mlmm_kwargs)
        H = calc.get_hessian(
            self.geom.atom_types,
            self.geom.coords
        )["hessian"].to(dtype=self.H_dtype, device=self.H_device)

        # 2. zero elements belonging to freeze_atoms
        # ---------------------------------------------------------------------
        freeze_set = set(self.freeze_atoms_static) | set(
            self.mlmm_kwargs.get("freeze_atoms", [])
        )
        if freeze_set:
            # build 1D DoF indices (three consecutive for x, y, z)
            dof_idx: list[int] = []
            for i in sorted(freeze_set):
                base = 3 * i
                dof_idx.extend((base, base + 1, base + 2))
            # zero out rows and columns
            H[dof_idx, :] = 0.0
            H[:, dof_idx] = 0.0

        # 3. vibrational analysis
        # ---------------------------------------------------------------------
        freqs, _, modes = calc_freq_from_hessian(
            H,
            self.geom.atomic_numbers,
            project_tr=True,
            coords_bohr=self.geom.coords.reshape(-1, 3),
            masses_amu=self.masses_amu_np,
        )

        # 4. cleanup
        # ---------------------------------------------------------------------
        del calc, H; torch.cuda.empty_cache()
        return freqs, modes

    # ================================================================
    # helper – run Dimer for *n* steps
    # ================================================================
    def _dimer_segment(self, threshold: str, n_steps: int) -> int:
        calc = MLMM(out_hess_torch=False, **self.mlmm_kwargs)
        dimer = Dimer(
            calculator=calc, N_raw=self.mode_path, mem=self.mem,
            write_orientations=False, seed=0, **self.dimer_kwargs
        )
        self.geom.freeze_atoms = self.freeze_atoms_static.copy()
        self.geom.set_calculator(dimer)
        opt = LBFGS(
            self.geom, max_cycles=n_steps, thresh=threshold, out_dir=self.out_dir,
            dump=self.dump, line_search=True
        )
        opt.run()
        if self.dump:
            with (self.out_dir / "optimization_all.trj").open("a") as f_all, \
                 (self.out_dir / "optimization.trj").open() as f_part:
                f_all.write(f_part.read())
        steps = opt.cur_cycle
        converged = opt.is_converged
        self.geom.set_calculator(None)
        del calc, dimer, opt
        torch.cuda.empty_cache()
        return steps, converged

    # ================================================================
    # helper – Dimer loop (loose or final)
    # ================================================================
    def _dimer_loop(self, threshold: str) -> int:
        total = 0
        while True:
            steps, ok = self._dimer_segment(threshold, self.update_interval_hessian)
            total += steps
            if ok:
                break
            # re-diagonalise partial Hessian & update mode
            H, frz = self._partial_hessian()
            self._write_lowest_mode(H, frz)
            del H; torch.cuda.empty_cache()
        return total

    # ================================================================
    # helper – flatten imaginary modes
    # ================================================================
    @staticmethod
    def _representative_atoms(mode_vec: torch.Tensor, k: int = 10) -> np.ndarray:
        """Return indices of the k atoms with the largest displacements."""
        vec = mode_vec.view(-1, 3)
        norms = torch.linalg.norm(vec, dim=1)
        return torch.argsort(norms)[-k:].cpu().numpy()

    def _flatten_once(
        self, freqs: np.ndarray, modes: torch.Tensor
    ) -> bool:  # returns True if flatten performed
        neg_idx = np.where(freqs < -abs(self.neg_freq_thresh))[0]
        if len(neg_idx) <= 1:
            return False
        # sort ascending (most negative first)
        sorted_neg = neg_idx[np.argsort(freqs[neg_idx])]
        coords_ang = torch.as_tensor(
            self.geom.coords.reshape(-1, 3) * BOHR2ANG,
            dtype=modes.dtype,
            device=modes.device,
        )
        targets = [int(sorted_neg[1])]  # start from 2nd most negative
        reps = [self._representative_atoms(modes[sorted_neg[1]], k=self.flatten_k)]
        for idx in sorted_neg[2:]:
            rep = self._representative_atoms(modes[idx], k=self.flatten_k)
            if all(
                torch.min(
                    torch.linalg.norm(
                        coords_ang[rep][:, None, :] - coords_ang[rp][None, :, :],
                        dim=2,
                    )
                ).item() >= self.flatten_sep_cutoff
                for rp in reps
            ):
                targets.append(int(idx))
                reps.append(rep)

        sp_calc = MLMM(out_hess_torch=False, **self.mlmm_kwargs)
        amp_bohr = self.flatten_amp_ang / BOHR2ANG

        for idx in targets:
            v = modes[idx].view(-1, 3).detach().cpu().numpy()
            v /= np.linalg.norm(v)
            disp = amp_bohr * self.mass_scale * v  # Bohr
            ref = self.geom.coords.reshape(-1, 3)
            plus, minus = ref + disp, ref - disp

            self.geom.coords = plus.flatten()
            self.geom.set_calculator(sp_calc)
            E_plus = self.geom.energy

            self.geom.coords = minus.flatten()
            self.geom.set_calculator(sp_calc)
            E_minus = self.geom.energy

            self.geom.coords = (plus if E_plus < E_minus else minus).flatten()
        self.geom.set_calculator(None)
        del sp_calc
        torch.cuda.empty_cache()
        return True

    # ================================================================
    # run
    # ================================================================
    def run(self):
        wall0 = time.time()

        if self.max_cycles_preopt == 0:
            # 0  no pre-opt
            # ---------------------------------------------------------------------
            print("\n>>> No pre-optimization\n")
            self.geom.freeze_atoms = self.freeze_atoms_static.copy()
            torch.cuda.empty_cache()
        else:
            # 0  pre-opt
            # ---------------------------------------------------------------------
            print("\n>>> Pre-optimization\n")
            calc_pre = MLMM(out_hess_torch=False, **self.mlmm_kwargs)
            ml_idx = list(getattr(calc_pre, "selection_indices",
                                getattr(calc_pre, "core").selection_indices))
            self.ml_indices = ml_idx
            self.geom.freeze_atoms = sorted(set(self.freeze_atoms_static) | set(ml_idx))
            self.geom.set_calculator(calc_pre)
            LBFGS(
                self.geom, max_cycles=self.max_cycles_preopt, thresh=self.thresh,
                out_dir=self.out_dir, dump=self.dump, line_search=True
            ).run()
            self.geom.freeze_atoms = self.freeze_atoms_static.copy()
            torch.cuda.empty_cache()

        # 1  build partial Hessian & loose loop
        # ---------------------------------------------------------------------
        print("\n>>> Loose-threshold Dimer with partial Hessian\n")
        H, frz = self._partial_hessian()
        self._write_lowest_mode(H, frz)
        del H; torch.cuda.empty_cache()
        cycles_loose = self._dimer_loop(self.thresh_loose)

        # 2  final-threshold Dimer loop
        # ---------------------------------------------------------------------
        print("\n>>> Final-threshold Dimer\n")
        H, frz = self._partial_hessian()
        self._write_lowest_mode(H, frz)
        del H; torch.cuda.empty_cache()
        cycles_final = self._dimer_loop(self.thresh)

        total_cycles = cycles_loose + cycles_final

        # 3  flatten loop
        # ---------------------------------------------------------------------
        print("\n>>> Flatten extra imaginary modes\n")
        for it in range(self.flatten_max_iter):
            freqs, modes = self._full_hessian()
            n_imag = np.sum(freqs < -abs(self.neg_freq_thresh))
            print(f"  flatten iter {it:02d}: n_imag = {n_imag}")
            if n_imag <= 1:
                break
            did_flatten = self._flatten_once(freqs, modes)
            del freqs, modes; torch.cuda.empty_cache()
            if not did_flatten:
                break
            # after displacement → run final-threshold Dimer again
            H, frz = self._partial_hessian()
            self._write_lowest_mode(H, frz)
            del H; torch.cuda.empty_cache()
            total_cycles += self._dimer_loop(self.thresh)
        else:
            print("  Warning: flatten_max_iter reached.")

        # 4  export final structure
        # ---------------------------------------------------------------------
        final_xyz = self.out_dir / "final_geometry.xyz"
        final_pdb = self.out_dir / "final_geometry.pdb"

        init = read(self.real_pdb); final_geom = read(final_xyz)
        atoms = init.copy()             
        atoms.set_positions(final_geom.get_positions()) 
        write(final_pdb, atoms)

        print("Saved final geometry:", final_pdb)

        if os.path.exists(self.out_dir / "optimization.trj"):
            os.remove(self.out_dir / "optimization.trj")
        if os.path.exists(self.out_dir / "optimization.h5"):
            os.remove(self.out_dir / "optimization.h5")

        # 5  final vibrational analysis
        # ---------------------------------------------------------------------
        freqs, modes = self._full_hessian()
        neg_idx = np.where(freqs < -abs(self.neg_freq_thresh))[0]
        print(f"Final n_imag = {len(neg_idx)}: ", freqs[neg_idx] if neg_idx.size else "none")

        atoms = Atoms(
            symbols=self.geom.atoms,
            positions=(self.geom.coords * BOHR2ANG).reshape(-1, 3),
            pbc=False,
        )
        for rank, idx in enumerate(neg_idx):
            freq = freqs[idx]
            mode_vec = modes[idx].view(len(self.geom.atomic_numbers), 3).cpu().numpy()
            out_xyz = self.vib_dir / f"mode{rank:02d}_{freq:+.2f}cm-1.xyz"
            write_vib_traj_xyz(
                atoms=atoms, mode=mode_vec, filename=str(out_xyz),
                amplitude=1.0, n_frames=20,
                comment=f"mode {rank}  freq={freq:+.1f} cm^-1"
            )
            init = read(self.real_pdb)
            traj = read(out_xyz, index=":", format="xyz")
            for i, frame in enumerate(traj):
                pdb_frame = init.copy()
                pdb_frame.set_positions(frame.get_positions())
                write(out_xyz.with_suffix(".pdb"), pdb_frame, append=i != 0)

        # 6  summary
        # ---------------------------------------------------------------------
        h, rem = divmod(int(time.time() - wall0), 3600)
        m, s = divmod(rem, 60)
        print(f"Total Dimer cycles      : {total_cycles}")
        print(f"Elapsed time            : {h} h {m} m {s} s")
