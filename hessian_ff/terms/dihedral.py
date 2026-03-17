from __future__ import annotations

import torch
from torch import nn


def _dihedral_angle_from_vectors(
    v1: torch.Tensor,
    v2: torch.Tensor,
    v3: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """OpenMM-style signed dihedral from 3 bond-like vectors.

    Returns
    -------
    phi, cross12, cross23
    """
    cross12 = torch.cross(v1, v2, dim=-1)
    cross23 = torch.cross(v2, v3, dim=-1)

    n12 = torch.linalg.norm(cross12, dim=-1).clamp_min(1.0e-12)
    n23 = torch.linalg.norm(cross23, dim=-1).clamp_min(1.0e-12)
    cos_phi = torch.sum(cross12 * cross23, dim=-1) / (n12 * n23)
    cos_phi = torch.clamp(cos_phi, -1.0, 1.0)
    phi = torch.acos(cos_phi)

    sign = torch.where(
        torch.sum(v1 * cross23, dim=-1) < 0.0,
        -torch.ones_like(phi),
        torch.ones_like(phi),
    )
    phi = phi * sign
    return phi, cross12, cross23


def _dihedral_angle(p0: torch.Tensor, p1: torch.Tensor, p2: torch.Tensor, p3: torch.Tensor) -> torch.Tensor:
    """Return OpenMM-compatible signed dihedral angle (radians) for batched points."""
    v1 = p0 - p1
    v2 = p2 - p1
    v3 = p2 - p3
    phi, _, _ = _dihedral_angle_from_vectors(v1, v2, v3)
    return phi


def _accumulate_dihedral_forces(
    force: torch.Tensor,
    coords: torch.Tensor,
    i: torch.Tensor,
    j: torch.Tensor,
    k: torch.Tensor,
    l: torch.Tensor,
    dE_dphi: torch.Tensor,
) -> None:
    """Accumulate torsion force for given dE/dphi following OpenMM reference formulas."""
    if i.numel() == 0:
        return

    v1 = coords[i] - coords[j]  # A-B
    v2 = coords[k] - coords[j]  # C-B
    v3 = coords[k] - coords[l]  # C-D

    cross12 = torch.cross(v1, v2, dim=-1)
    cross23 = torch.cross(v2, v3, dim=-1)

    norm_cross12_sq = torch.sum(cross12 * cross12, dim=-1).clamp_min(1.0e-24)
    norm_cross23_sq = torch.sum(cross23 * cross23, dim=-1).clamp_min(1.0e-24)
    norm_v2 = torch.linalg.norm(v2, dim=-1).clamp_min(1.0e-12)
    norm_v2_sq = torch.sum(v2 * v2, dim=-1).clamp_min(1.0e-24)

    f0 = (-dE_dphi * norm_v2) / norm_cross12_sq
    f3 = (dE_dphi * norm_v2) / norm_cross23_sq
    f1 = torch.sum(v1 * v2, dim=-1) / norm_v2_sq
    f2 = torch.sum(v3 * v2, dim=-1) / norm_v2_sq

    ff0 = f0.unsqueeze(-1) * cross12
    ff3 = f3.unsqueeze(-1) * cross23
    s = f1.unsqueeze(-1) * ff0 - f2.unsqueeze(-1) * ff3
    ff1 = ff0 - s
    ff2 = ff3 + s

    force.index_add_(0, i, ff0)
    force.index_add_(0, j, -ff1)
    force.index_add_(0, k, -ff2)
    force.index_add_(0, l, ff3)


class DihedralTerm(nn.Module):
    """AMBER torsion term: E = sum k (1 + cos(n*phi - phase))."""

    def __init__(
        self,
        i: torch.Tensor,
        j: torch.Tensor,
        k: torch.Tensor,
        l: torch.Tensor,
        force: torch.Tensor,
        period: torch.Tensor,
        phase: torch.Tensor,
    ):
        super().__init__()
        self.register_buffer("i", i.long())
        self.register_buffer("j", j.long())
        self.register_buffer("k", k.long())
        self.register_buffer("l", l.long())
        self.register_buffer("force", force)
        self.register_buffer("period", period)
        self.register_buffer("phase", phase)

    def forward(self, coords: torch.Tensor) -> torch.Tensor:
        if self.i.numel() == 0:
            return coords.new_zeros(())
        p0 = coords[self.i]
        p1 = coords[self.j]
        p2 = coords[self.k]
        p3 = coords[self.l]
        phi = _dihedral_angle(p0, p1, p2, p3)

        # Periodicity can be negative in some prmtops; AMBER uses its absolute value.
        n = torch.abs(self.period)
        return torch.sum(self.force * (1.0 + torch.cos(n * phi - self.phase)))

    def energy_force(self, coords: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        """Return torsion energy and analytical force."""
        force = torch.zeros_like(coords)
        if self.i.numel() == 0:
            return coords.new_zeros(()), force

        p0 = coords[self.i]
        p1 = coords[self.j]
        p2 = coords[self.k]
        p3 = coords[self.l]
        phi = _dihedral_angle(p0, p1, p2, p3)

        n = torch.abs(self.period)
        delta = n * phi - self.phase
        energy = torch.sum(self.force * (1.0 + torch.cos(delta)))

        # dE/dphi
        dE_dphi = -self.force * n * torch.sin(delta)
        _accumulate_dihedral_forces(force, coords, self.i, self.j, self.k, self.l, dE_dphi)
        return energy, force
