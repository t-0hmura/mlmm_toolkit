from __future__ import annotations

import torch
from torch import nn


class BondTerm(nn.Module):
    """AMBER bond term: E = sum k (r - r0)^2"""

    def __init__(self, i: torch.Tensor, j: torch.Tensor, k: torch.Tensor, r0: torch.Tensor):
        super().__init__()
        self.register_buffer("i", i.long())
        self.register_buffer("j", j.long())
        self.register_buffer("k", k)
        self.register_buffer("r0", r0)

    def forward(self, coords: torch.Tensor) -> torch.Tensor:
        if self.i.numel() == 0:
            return coords.new_zeros(())
        rij = coords[self.j] - coords[self.i]
        r = torch.linalg.norm(rij, dim=-1)
        return torch.sum(self.k * (r - self.r0) ** 2)

    def energy_force(self, coords: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        """Return bond energy and analytical force."""
        force = torch.zeros_like(coords)
        if self.i.numel() == 0:
            return coords.new_zeros(()), force

        rij = coords[self.j] - coords[self.i]
        r2 = torch.sum(rij * rij, dim=-1).clamp_min(1.0e-24)
        inv_r = torch.rsqrt(r2)
        r = r2 * inv_r

        dr = r - self.r0
        e = torch.sum(self.k * dr * dr)

        # E = k (r-r0)^2, dE/dr = 2k(r-r0), F_i = dE/dr * (r_ij/r)
        fscale = 2.0 * self.k * dr * inv_r
        fij = fscale.unsqueeze(-1) * rij

        force.index_add_(0, self.i, fij)
        force.index_add_(0, self.j, -fij)
        return e, force
