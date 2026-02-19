from __future__ import annotations

import torch
from torch import nn


class AngleTerm(nn.Module):
    """AMBER angle term: E = sum k (theta - theta0)^2

    theta is computed via atan2(|u x v|, u·v) for numerical stability.
    """

    def __init__(
        self,
        i: torch.Tensor,
        j: torch.Tensor,
        k: torch.Tensor,
        k_theta: torch.Tensor,
        theta0: torch.Tensor,
    ):
        super().__init__()
        self.register_buffer("i", i.long())
        self.register_buffer("j", j.long())
        self.register_buffer("k", k.long())
        self.register_buffer("k_theta", k_theta)
        self.register_buffer("theta0", theta0)

    def forward(self, coords: torch.Tensor) -> torch.Tensor:
        if self.i.numel() == 0:
            return coords.new_zeros(())
        u = coords[self.i] - coords[self.j]
        v = coords[self.k] - coords[self.j]
        dot = torch.sum(u * v, dim=-1)
        cross = torch.linalg.norm(torch.cross(u, v, dim=-1), dim=-1)
        theta = torch.atan2(cross, dot)
        return torch.sum(self.k_theta * (theta - self.theta0) ** 2)

    def energy_force(self, coords: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        """Return angle energy and analytical force."""
        force = torch.zeros_like(coords)
        if self.i.numel() == 0:
            return coords.new_zeros(()), force

        # OpenMM-style geometry for numerical stability.
        d0 = coords[self.j] - coords[self.i]
        d1 = coords[self.j] - coords[self.k]
        p = torch.cross(d0, d1, dim=-1)

        r20 = torch.sum(d0 * d0, dim=-1).clamp_min(1.0e-24)
        r21 = torch.sum(d1 * d1, dim=-1).clamp_min(1.0e-24)
        rp = torch.linalg.norm(p, dim=-1).clamp_min(1.0e-12)
        dot = torch.sum(d0 * d1, dim=-1)
        cos_theta = dot / torch.sqrt(r20 * r21)
        cos_theta = torch.clamp(cos_theta, -1.0, 1.0)
        theta = torch.acos(cos_theta)

        dtheta = theta - self.theta0
        e = torch.sum(self.k_theta * dtheta * dtheta)

        # E = k (theta-theta0)^2, dE/dtheta = 2k(theta-theta0)
        dE_dtheta = 2.0 * self.k_theta * dtheta

        term_i = dE_dtheta / (r20 * rp)
        term_k = -dE_dtheta / (r21 * rp)

        fi = torch.cross(d0, p, dim=-1) * term_i.unsqueeze(-1)
        fk = torch.cross(d1, p, dim=-1) * term_k.unsqueeze(-1)
        fj = -(fi + fk)

        force.index_add_(0, self.i, fi)
        force.index_add_(0, self.j, fj)
        force.index_add_(0, self.k, fk)
        return e, force
