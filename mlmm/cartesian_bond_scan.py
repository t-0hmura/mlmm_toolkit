"""
CartesianBondScan – ML/MM geometry scan along a single bond
-----------------------------------------------------------

* High-level ML / Low-level MM handled by MLMM calculator.
* Performs an inward (shortening) and outward (elongating) scan.
* Each scan point is locally optimized with LBFGS.
"""
from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Tuple, Sequence, Union, Optional

import numpy as np
from mlmm import mlmm as MLMM
from pysisyphus.io.xyz import geom_from_xyz
from pysisyphus.io.pdb import geom_from_pdb
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.constants import BOHR2ANG, ANG2BOHR


class CartesianBondScan:
    """Wrapper around the original procedural scan script."""

    # ---------------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------------
    def __init__(
        self,
        # scan settings
        # ---------------------------------------------------------------------
        input_path: str,
        scan_bond: Tuple[int, int],
        scan_step: float,                 # Å
        scan_range: Tuple[float, float], # Å (min,max)
        init_thresh: str,
        thresh: str,
        max_cycles: int,
        out_dir: str,
        freeze_atoms: List[int] = [],
        out_fn: str = "final_geometries.trj",
        # MLMM calculator
        # ---------------------------------------------------------------------
        mlmm_kwargs: Optional[dict] = None,
    ):
        self.input_path = input_path
        self.freeze_atoms = freeze_atoms
        self.scan_bond = scan_bond
        self.scan_step = scan_step
        self.scan_range = scan_range
        self.init_thresh = init_thresh
        self.thresh = thresh
        self.max_cycles = max_cycles
        self.out_dir = out_dir
        self.out_fn = out_fn

        # calculator
        self.mlmm_kwargs = {} if mlmm_kwargs is None else mlmm_kwargs
        self.calc = MLMM(**self.mlmm_kwargs)

        # internal paths
        self.out_path = os.path.join(self.out_dir, self.out_fn)
        self.scan_step_bohr = self.scan_step * ANG2BOHR
        self.scan_range_bohr = [x * ANG2BOHR for x in self.scan_range]

        # ensure output directory
        Path(self.out_dir).mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------------
    # Static helpers
    # ---------------------------------------------------------------------
    @staticmethod
    def _read_geom(path: str):
        if path.endswith(".pdb"):
            return geom_from_pdb(path)
        elif path.endswith(".xyz"):
            return geom_from_xyz(path)
        raise ValueError("Unsupported file format (must be .xyz or .pdb)")

    @staticmethod
    def _bond_length(coords: np.ndarray, idx_pair: Tuple[int, int]) -> float:
        return np.linalg.norm(coords[idx_pair[0]] - coords[idx_pair[1]])

    # ---------------------------------------------------------------------
    # Core routine
    # ---------------------------------------------------------------------
    def run(self):
        start = time.time()

        # 1. initial optimization with frozen scan atoms -----------------
        geom = self._read_geom(self.input_path)
        geom.set_calculator(self.calc)
        geom.freeze_atoms = sorted(list(self.scan_bond) + list(self.freeze_atoms))
        opt = LBFGS(geom, max_cycles=self.max_cycles,
                    thresh=self.init_thresh, out_dir=self.out_dir)
        opt.run()

        # write first structure
        with open(self.out_path, "w") as f:
            f.write(geom.as_xyz() + "\n")

        # reference data
        direction = self._unit_bond_vector(geom.coords.reshape(-1, 3))
        # separate inward / outward geoms so they keep independent histories
        self._scan_branch(
            base_geom_path=self.input_path,
            direction_vector= direction,
            inward=True
        )
        self._scan_branch(
            base_geom_path=self.input_path,
            direction_vector= direction,
            inward=False
        )

        # ---------------------------------------------------------------------
        elapsed = time.time() - start
        h, m, s = map(int, [elapsed // 3600, (elapsed % 3600) // 60, elapsed % 60])
        print(f"Cartesian Scan took: {h}h {m}m {s}s")

    # ---------------------------------------------------------------------
    # Helper to get unit bond vector
    # ---------------------------------------------------------------------
    def _unit_bond_vector(self, coords: np.ndarray) -> np.ndarray:
        v = coords[self.scan_bond[1]] - coords[self.scan_bond[0]]
        v /= np.linalg.norm(v)
        return v

    # ---------------------------------------------------------------------
    # Branch scan (inward / outward)
    # ---------------------------------------------------------------------
    def _scan_branch(
        self,
        base_geom_path: str,
        direction_vector: np.ndarray,
        inward: bool = True,
        max_steps: int = 10_000,
    ):
        geom = self._read_geom(base_geom_path)
        geom.set_calculator(self.calc)

        # select sign for displacement
        sign = +1 if inward else -1

        for _ in range(max_steps):
            coords = geom.coords.reshape(-1, 3)

            # displace scan bond atoms by half-step each
            disp = sign * direction_vector * self.scan_step_bohr / 2.0
            coords[self.scan_bond[0]] += disp
            coords[self.scan_bond[1]] -= disp
            geom.coords = coords.flatten()

            current_len = self._bond_length(coords, self.scan_bond)
            if ((inward and current_len < self.scan_range_bohr[0]) or
                (not inward and current_len > self.scan_range_bohr[1])):
                break

            # optimize
            geom.freeze_atoms = sorted(list(self.scan_bond) + list(self.freeze_atoms))
            opt = LBFGS(geom, max_cycles=self.max_cycles,
                        thresh=self.thresh, out_dir=self.out_dir)
            opt.run()

            # save trajectory (prepend for inward, append for outward)
            mode = "r+" if inward else "a"
            with open(self.out_path, mode) as f:
                if inward:                    # write at beginning
                    existing = f.read()
                    f.seek(0, 0)
                    f.write(geom.as_xyz() + "\n" + existing)
                else:
                    f.write(geom.as_xyz() + "\n")

"""
Example
-------
scanner = CartesianBondScan(
    input_path="./dump/opt1/final_geometry.xyz",
    scan_bond=(100, 110),
    scan_step=0.05,                 # Å
    scan_range=(1.4, 4.0),         # Å (min, max)
    init_thresh="gau",
    thresh="gau",
    max_cycles=10_000,
    out_dir="./dump/cart_scan/",
    out_fn="final_geometries.trj",
    mlmm_kwargs=dict(
        real_pdb="./parm/complex.pdb",
        real_parm7="./parm/complex.parm7",
        real_rst7="./parm/complex.rst7",
        model_pdb="./parm/ml_region.pdb",
        model_charge=1,
        model_mult=1,
        backend="aimnet2",
        ml_device="auto",
        ml_cuda_idx=0,
        mm_device="cpu",
        mm_cuda_idx=0,
        mm_threads=16,
        mem=100_000,
    ),
)
scanner.run()
"""
