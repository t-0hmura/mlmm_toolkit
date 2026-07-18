#!/usr/bin/env python3
"""Strict real-model analytical-Hessian smoke for an ML/MM MLIP backend."""

from __future__ import annotations

import argparse
import json

import numpy as np
import torch
from ase import Atoms


def _finite_difference_hessian(backend, atoms: Atoms, eps_ang: float = 1.0e-3) -> np.ndarray:
    positions = np.asarray(atoms.positions, dtype=np.float64)
    dof = positions.size
    hessian = np.empty((dof, dof), dtype=np.float64)
    for column in range(dof):
        atom, axis = divmod(column, 3)
        plus = atoms.copy()
        minus = atoms.copy()
        plus.positions[atom, axis] += eps_ang
        minus.positions[atom, axis] -= eps_ang
        force_plus = np.asarray(backend.eval(plus, need_grad=False)[1], dtype=np.float64).reshape(-1)
        force_minus = np.asarray(backend.eval(minus, need_grad=False)[1], dtype=np.float64).reshape(-1)
        hessian[:, column] = -(force_plus - force_minus) / (2.0 * eps_ang)
    return 0.5 * (hessian + hessian.T)


def _metrics(analytical: np.ndarray, finite_difference: np.ndarray) -> dict[str, float]:
    difference = analytical - finite_difference
    fd_norm = max(float(np.linalg.norm(finite_difference)), 1.0e-12)
    fd_max = max(float(np.max(np.abs(finite_difference))), 1.0e-12)
    return {
        "relative_frobenius_error": float(np.linalg.norm(difference)) / fd_norm,
        "relative_max_error": float(np.max(np.abs(difference))) / fd_max,
        "analytical_asymmetry": float(np.max(np.abs(analytical - analytical.T))),
        "fd_asymmetry": float(np.max(np.abs(finite_difference - finite_difference.T))),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("backend", choices=("uma", "orb", "mace", "aimnet2"))
    parser.add_argument("--max-relative-frobenius", type=float, default=0.20)
    parser.add_argument("--max-relative-element", type=float, default=0.35)
    args = parser.parse_args()

    if not torch.cuda.is_available():
        raise SystemExit("CUDA is required for the real-model analytical-Hessian lane.")

    if args.backend == "uma":
        from mlmm.backends.mlmm_calc import _UMABackend

        backend = _UMABackend(
            uma_model="uma-s-1p2",
            uma_task_name="omol",
            precision="fp32",
            workers=1,
            workers_per_node=1,
            model_charge=0,
            model_mult=1,
            ml_device=torch.device("cuda"),
        )
    elif args.backend == "orb":
        from mlmm.backends.mlmm_calc import _OrbBackend

        backend = _OrbBackend(
            orb_model="orb_v3_conservative_omol",
            orb_precision="float64",
            model_charge=0,
            model_mult=1,
            ml_device=torch.device("cuda"),
        )
    elif args.backend == "mace":
        from mlmm.backends.mlmm_calc import _MACEBackend

        backend = _MACEBackend(
            mace_model="MACE-OMOL-0",
            mace_dtype="float64",
            model_charge=0,
            model_mult=1,
            ml_device=torch.device("cuda"),
        )
    else:
        from mlmm.backends.mlmm_calc import _AIMNet2Backend

        backend = _AIMNet2Backend(
            aimnet2_model="aimnet2",
            model_charge=0,
            model_mult=1,
            ml_device=torch.device("cuda"),
        )

    atoms = Atoms("H2", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
    _, _, opaque = backend.eval(atoms, need_grad=True)
    analytical = (
        backend.hessian_analytical(opaque, len(atoms), dtype=torch.float64)
        .detach()
        .cpu()
        .numpy()
        .reshape(6, 6)
    )
    finite_difference = _finite_difference_hessian(backend, atoms)

    if not np.isfinite(analytical).all() or not np.isfinite(finite_difference).all():
        raise SystemExit("Hessian contains non-finite values.")
    metrics = _metrics(analytical, finite_difference)
    # The analytical Hessian is formed by GPU autograd double-backward, so its raw
    # asymmetry is bounded by the backend's working precision (fp32 UMA lands near
    # 1e-8), not by an fp64-scale 1e-10. Check the asymmetry relative to the Hessian
    # magnitude; the pipeline symmetrizes before eigendecomposition anyway.
    analytical_scale = max(float(np.max(np.abs(analytical))), 1.0e-12)
    if metrics["analytical_asymmetry"] > 1.0e-5 * analytical_scale:
        raise SystemExit(f"Analytical Hessian is not symmetric: {metrics}")
    if metrics["relative_frobenius_error"] > args.max_relative_frobenius:
        raise SystemExit(f"Analytical/FD Frobenius mismatch: {metrics}")
    if metrics["relative_max_error"] > args.max_relative_element:
        raise SystemExit(f"Analytical/FD element mismatch: {metrics}")

    print(json.dumps({"backend": args.backend, "shape": [6, 6], **metrics}, sort_keys=True))


if __name__ == "__main__":
    main()
