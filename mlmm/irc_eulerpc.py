"""
Mass-weighted Euler Predictor–Corrector IRC search (ML/MM).

Reference
---------
EulerPC itself.
    Hratchian, H. P., & Schlegel, H. B. (2004). Accurate reaction paths using a Hessian based predictor-corrector integrator. Journal of Chemical Physics, 120(21), 9918–9924. https://doi.org/10.1063/1.1724823
    Hratchian, H. P., Frisch, M. J., & Schlegel, H. B. (2010). Steepest descent reaction path integration using a first-order predictor-corrector method. Journal of Chemical Physics, 133(22). https://doi.org/10.1063/1.3514202

Pysisyphus
    Steinmetzer, J., Kupfer, S., & Gräfe, S. (2021). pysisyphus: Exploring potential energy surfaces in ground and excited states. International Journal of Quantum Chemistry, 121(3). https://doi.org/10.1002/qua.26390

This script reimplements the EulerPC algorithm for IRC based on the pysisyphus.irc.EulerPC.

Usage
-----------
# YAML-driven
$ irc_eulerpc input.yaml

YAML example:
  geom:
    fn: ./dump/dimer/final_geometry.xyz
    freeze_atoms: [0, 1]        # optional
  irc:
    max_steps:   25
    step_size:   0.10           # amu½·Bohr
    pc_steps:    500
    out_dir:     ./dump/irc/
    device:      auto           # "auto" | "cpu" | "cuda:0" …
  calc:                         # forwarded to MLMM(**calc)
    real_pdb:   ./parm/complex.pdb
    real_parm7: ./parm/complex.parm7
    ...

Outputs
-------
Trajectory  →  ./dump/irc/final_geometries.trj
(change via YAML key `irc.out_dir`)
"""
from __future__ import annotations

# ------------------ stdlib │ 3rd-party -----------------------------
import argparse
from pathlib import Path
from typing import Final, List, Optional

import time
import numpy as np
import torch
import yaml
from ase.data import atomic_masses
from ase.io import read, write
from mlmm.mlmm_pysis import mlmm as MLMM
from pysisyphus.constants import AMU2AU, ANG2BOHR, BOHR2ANG

# ------------------- constants -------------------------------------
EPS: Final[float] = 1.0e-12           # numerical guard
RICHARDSON_MAX_K: Final[int] = 15    # 2**15 ≈ 65 k Richardson points

# ------------------- utilities -------------------------------------
def _sr1_update(z: torch.Tensor, dx: torch.Tensor) -> torch.Tensor:
    return torch.outer(z, z) / (torch.dot(z, dx) + EPS)


def _psb_update(z: torch.Tensor, dx: torch.Tensor) -> torch.Tensor:
    inv_dx2 = 1.0 / (torch.dot(dx, dx) + EPS)
    return ((torch.outer(dx, z) + torch.outer(z, dx)) * inv_dx2
            - torch.dot(dx, z) * torch.outer(dx, dx) * inv_dx2**2)


def bofill_update(H: torch.Tensor, dx: torch.Tensor, dg: torch.Tensor) -> None:
    """In-place Bofill update (mass-weighted)."""
    z   = dg - H @ dx
    mix = (torch.dot(z, dx) ** 2) / ((torch.dot(z, z) + EPS) *
                                     (torch.dot(dx, dx) + EPS))
    H += mix * _sr1_update(z, dx) + (1.0 - mix) * _psb_update(z, dx)

# ------------------- distance-weighted fit ------------------
class DWI:
    """Two-point distance-weighted interpolator (GPU-friendly)."""

    def __init__(
        self,
        power: int = 4,
        *,
        device: torch.device | None = None,
        dtype: torch.dtype = torch.float64,
    ):
        if power <= 0 or power % 2:
            raise ValueError("`power` must be a positive even integer.")
        self.n      = int(power)
        self.device = torch.device("cpu") if device is None else device
        self.dtype  = dtype
        self.coords:    list[torch.Tensor] = []
        self.energies:  list[float]        = []
        self.gradients: list[torch.Tensor] = []
        self.hessians:  list[torch.Tensor] = []

    # --------------------- storage ---------------------
    def update(
        self,
        coords: torch.Tensor,
        energy: float,
        gradient: torch.Tensor,
        hessian: torch.Tensor,
    ) -> None:
        """Keep only the two most recent reference points."""
        if len(self.coords) == 2:
            self.coords.pop(0)
            self.energies.pop(0)
            self.gradients.pop(0)
            self.hessians.pop(0)

        self.coords.append(coords.to(self.device, dtype=self.dtype))
        self.energies.append(float(energy))
        self.gradients.append(gradient.to(self.device, dtype=self.dtype))
        self.hessians.append(hessian.to(self.device, dtype=self.dtype))

    # --------------------- interpolation ---------------------
    @staticmethod
    def _taylor(
        energy: float, grad: torch.Tensor, hess: torch.Tensor, step: torch.Tensor
    ) -> torch.Tensor:
        """Second-order Taylor polynomial of *E* at a displaced point."""
        return energy + torch.dot(step, grad) + 0.5 * torch.dot(step, hess @ step)

    @staticmethod
    def _taylor_grad(
        grad: torch.Tensor, hess: torch.Tensor, step: torch.Tensor
    ) -> torch.Tensor:
        """Gradient of the Taylor polynomial."""
        return grad + hess @ step

    def interpolate(
        self, coords: torch.Tensor, *, gradient: bool = False
    ) -> float | tuple[float, torch.Tensor]:
        """Energy (and optionally gradient) at *coords* by DWI."""
        if len(self.coords) < 2:
            raise RuntimeError("DWI requires two stored points.")

        c1, c2   = self.coords
        dx1, dx2 = coords - c1, coords - c2
        d1, d2   = torch.linalg.norm(dx1), torch.linalg.norm(dx2)

        # exact match → return cached value
        if d1 < 1.0e-16:
            return (self.energies[0] if not gradient
                    else (self.energies[0], self.gradients[0].clone()))
        if d2 < 1.0e-16:
            return (self.energies[1] if not gradient
                    else (self.energies[1], self.gradients[1].clone()))

        # distance weights
        w1, w2 = d2**self.n, d1**self.n
        denom  = w1 + w2
        w1, w2 = w1 / denom, w2 / denom

        e1, e2 = self.energies
        g1, g2 = self.gradients
        h1, h2 = self.hessians

        t1 = self._taylor(e1, g1, h1, dx1)
        t2 = self._taylor(e2, g2, h2, dx2)
        energy = w1 * t1 + w2 * t2

        if not gradient:
            return energy.item()

        t1_g = self._taylor_grad(g1, h1, dx1)
        t2_g = self._taylor_grad(g2, h2, dx2)

        dw1 = (self.n * d2 ** (self.n - 2)) * dx2 * d1**self.n \
            - (self.n * d1 ** (self.n - 2)) * dx1 * d2**self.n
        dw1 /= denom**2
        grad_dwi = dw1 * t1 + w1 * t1_g - dw1 * t2 + w2 * t2_g
        return energy.item(), grad_dwi

# ------------------- IRC driver -----------------------------------
class EulerPC:
    """
    Euler Predictor–Corrector IRC driver in mass-weighted space.
    Traces forward (products) and backward (reactants) branches.
    """

    # ------------- init -------------
    def __init__(
        self,
        ts_xyz: str | Path,
        *,
        step_length: float = 0.10,       # amu½·Bohr
        max_cycles: int = 25,
        max_pred_steps: int = 500,
        out_dir: str | Path = "./dump/irc/",
        freeze_atoms: Optional[List[int]] = None,
        mlmm_kwargs: Optional[dict] = None,
        H_dtype: torch.dtype = torch.float32,
        device: str | torch.device = "auto",
    ):
        # geometry
        self.atoms  = read(str(ts_xyz))
        self.natom  = len(self.atoms)

        # device / dtype
        self.device = (torch.device("cuda") if device == "auto" and
                       torch.cuda.is_available()
                       else torch.device(device))
        self.dtype  = H_dtype

        # Cartesian → Bohr
        self.coords0 = torch.as_tensor(self.atoms.get_positions() * ANG2BOHR,
                                       dtype=self.dtype, device=self.device).view(-1)

        # masses
        masses_au   = torch.as_tensor(
            np.array([atomic_masses[z] for z in self.atoms.get_atomic_numbers()])
            * AMU2AU, dtype=self.dtype, device=self.device)
        self.m_sqrt = torch.sqrt(masses_au)               # (N,)
        self.m_sqrt_3N = self.m_sqrt.repeat_interleave(3) # (3 N,)

        # freeze mask
        self.freeze_atoms = [] if freeze_atoms is None else list(freeze_atoms)
        self.freeze_mask  = torch.ones_like(self.coords0, dtype=self.dtype)
        for idx in self.freeze_atoms:
            self.freeze_mask[3*idx : 3*idx+3] = 0.0

        # run parameters
        self.step_length    = float(step_length)
        self.max_cycles     = int(max_cycles)
        self.max_pred_steps = int(max_pred_steps)
        self.out_dir        = Path(out_dir).resolve()
        self.out_dir.mkdir(parents=True, exist_ok=True)

        self.mlmm_kwargs = {} if mlmm_kwargs is None else dict(mlmm_kwargs)

    # ------------- mass-weight helpers -------------
    def _mw(self,  x: torch.Tensor) -> torch.Tensor: return x * self.m_sqrt_3N
    def _unmw(self, x: torch.Tensor) -> torch.Tensor: return x / self.m_sqrt_3N
    def _grad_mw(self, g: torch.Tensor) -> torch.Tensor: return g / self.m_sqrt_3N

    def _conv_fact(self, G: torch.Tensor, min_fact: float = 2.0) -> float:
        """Scale predictor sub-step count by gradient norm ratio."""
        return max(min_fact,
                   torch.linalg.norm(G * self.m_sqrt_3N).item()
                   / (torch.linalg.norm(G).item() + EPS))

    # ------------- ML/MM backend -------------
    def _build_calc(self, *, out_hess: bool = False):
        cfg = dict(self.mlmm_kwargs)
        cfg["out_hess_torch"] = out_hess
        cfg.setdefault("vib_run", True)
        if self.freeze_atoms and "freeze_atoms" not in cfg:
            cfg["freeze_atoms"] = self.freeze_atoms
        return MLMM(**cfg)

    def _compute(
        self, coords_bohr: torch.Tensor, *, hessian: bool = False
    ) -> tuple[float, torch.Tensor, Optional[torch.Tensor]]:
        """Energy, gradient, (optional) Hessian at *coords_bohr*."""
        coords_np = coords_bohr.view(self.natom, 3).detach().cpu().numpy()
        calc      = self._build_calc(out_hess=hessian)

        res = (calc.get_hessian if hessian else calc.get_forces)(
            self.atoms.get_chemical_symbols(), coords_np)

        energy = float(res["energy"])
        grad   = -torch.as_tensor(res["forces"], dtype=self.dtype,
                                  device=self.device).view(-1) * self.freeze_mask

        H_MW   = None
        if hessian:
            H_cart = torch.as_tensor(res["hessian"], dtype=self.dtype,
                                     device=self.device)
            H_cart.div_(self.m_sqrt_3N.unsqueeze(0)).div_(self.m_sqrt_3N.unsqueeze(1))
            frozen = torch.where(self.freeze_mask == 0.0)[0]
            if len(frozen):
                H_cart[frozen, :] = 0.0
                H_cart[:, frozen] = 0.0
            H_MW = H_cart

        del calc  # free GPU
        return energy, grad, H_MW

    # ------------- single-branch propagation -------------
    def _propagate(self, sign: int) -> list[np.ndarray]:
        """
        sign = +1 → products, -1 → reactants branch
        """
        # transition state
        E_ts, g_ts, H_ts = self._compute(self.coords0, hessian=True)
        Q_ts, G_ts       = self._mw(self.coords0), self._grad_mw(g_ts)

        # imaginary mode
        eigvals, eigvecs = torch.linalg.eigh(H_ts)
        ev = eigvecs[:, torch.argmin(eigvals)]
        ev /= torch.linalg.norm(ev) + EPS

        # first displaced point
        Q = Q_ts + sign * self.step_length * ev
        E, g, _ = self._compute(self._unmw(Q), hessian=False)
        G       = self._grad_mw(g)

        # DWI initialisation
        dwi = DWI(device=self.device, dtype=self.dtype)
        dwi.update(Q_ts, E_ts, G_ts, H_ts)
        dwi.update(Q,    E,    G,    H_ts)

        # trajectory container
        path = [self.coords0.view(-1, 3).cpu().numpy(),
                self._unmw(Q).view(-1, 3).cpu().numpy()]

        H = H_ts.clone()  # working Hessian

        for cycle in range(self.max_cycles):
            Q0 = Q.clone()

            # predictor (explicit Euler, MW)
            ds        = self.step_length / (self.max_pred_steps / self._conv_fact(G))
            Q_pred, G_pred = Q.clone(), G.clone()

            for _ in range(self.max_pred_steps):
                Q_pred -= ds * G_pred / (torch.linalg.norm(G_pred) + EPS)
                delta_Q = Q_pred - Q0
                G_pred  = G + H @ delta_Q
                if torch.linalg.norm(delta_Q).item() >= self.step_length:
                    break

            # evaluate predictor
            E_pred, g_pred, _ = self._compute(self._unmw(Q_pred), hessian=False)
            G_new = self._grad_mw(g_pred)

            # Bofill update
            bofill_update(H, Q_pred - Q0, G_new - G)

            # corrector (Richardson)
            dwi.update(Q_pred, E_pred, G_new, H)
            Q_corr = self._corrector_step(Q0, self.step_length, dwi)

            E_corr, g_corr, _ = self._compute(self._unmw(Q_corr), hessian=False)
            G_corr = self._grad_mw(g_corr)
            dwi.update(Q_corr, E_corr, G_corr, H)

            # advance
            Q, G = Q_corr, G_corr
            path.append(self._unmw(Q).view(-1, 3).cpu().numpy())
        return path

    # ------------- Richardson corrector -------------
    def _corrector_step(
        self, Q0: torch.Tensor, step_length: float, dwi: DWI
    ) -> torch.Tensor:
        """Return *Q* such that ‖Q − Q₀‖ = *step_length* (MW norm)."""
        richardson: dict[tuple[int, int], torch.Tensor] = {}
        for k in range(RICHARDSON_MAX_K + 1):
            n_points = 20 * 2**k
            ds       = step_length / (n_points - 1)

            Q = Q0.clone()
            while True:
                if abs(step_length - torch.linalg.norm(Q - Q0).item()) < 0.5 * ds:
                    break
                _, G = dwi.interpolate(Q, gradient=True)
                Q   -= ds * G / (torch.linalg.norm(G) + EPS)

            richardson[(k, 0)] = Q
            for j in range(1, k + 1):
                richardson[(k, j)] = ((2**j) * richardson[(k, j-1)]
                                      - richardson[(k-1, j-1)]) / (2**j - 1)

            if k > 0 and torch.linalg.norm(
                    richardson[(k, k)] - richardson[(k, k-1)]).item() <= 1.0e-5:
                break
        return richardson[(k, k)]

    # ------------- public API -------------
    def run(self) -> None:
        """Compute both IRC branches and write the merged trajectory."""
        start_time = time.time()

        print("Starting IRC path search with Euler Predictor–Corrector method...")

        print('Backward path search (TS → Reactants)')
        path_bwd = self._propagate(sign=-1)   # TS → Reactants
        torch.cuda.empty_cache()
        print('Forward path search (TS → Products)')
        path_fwd = self._propagate(sign=+1)   # TS → Products
        torch.cuda.empty_cache()

        merged   = path_bwd[::-1] + path_fwd[1:]  # drop duplicated TS
        traj_file = self.out_dir / "final_geometries.trj"
        for i, frame in enumerate(merged):
            at = self.atoms.copy()
            at.set_positions(frame * BOHR2ANG)
            write(traj_file, at, append=i != 0)

        elapsed = time.time() - start_time
        hours, minutes, seconds = int(elapsed // 3600), int((elapsed % 3600) // 60), elapsed % 60
        print(f"IRC path search completed in {hours}h {minutes}m {seconds:.2f}s.")

# ------------------- CLI helpers ----------------------------------
def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="irc_eulerpc",
        description="Mass-weighted Euler Predictor–Corrector IRC search (ML/MM)",
    )
    p.add_argument(
        "config",
        metavar="input.yaml",
        type=str,
        help="YAML configuration file describing geometry, IRC parameters and ML/MM backend.",
    )
    return p


def _load_yaml(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)

# ------------------- script entry ---------------------------------
def main() -> None:
    args = _build_parser().parse_args()

    cfg = _load_yaml(args.config)

    geom = cfg.get("geom", {})
    irc  = cfg.get("irc",  {})
    mlmm = cfg.get("calc", {})

    EulerPC(
        ts_xyz      = geom.get("fn", "./coord/ts.xyz"),
        step_length = irc.get("step_size", 0.10),
        max_cycles  = irc.get("max_steps", 25),
        max_pred_steps = irc.get("pc_steps", 500),
        out_dir     = irc.get("out_dir", "./dump/irc/"),
        device      = irc.get("device", "auto"),
        freeze_atoms= geom.get("freeze_atoms", []),
        mlmm_kwargs = mlmm,
    ).run()


if __name__ == "__main__":
    main()