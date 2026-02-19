from __future__ import annotations

import math
from typing import Tuple

import torch
from torch import nn

from .dihedral import _accumulate_dihedral_forces, _dihedral_angle

_TWO_PI = 2.0 * math.pi

# References:
# - OpenMM CMAP coefficient/derivative setup:
#   https://github.com/openmm/openmm/blob/4768436/openmmapi/src/CMAPTorsionForceImpl.cpp
# - OpenMM periodic spline routines:
#   https://github.com/openmm/openmm/blob/4768436/openmmapi/src/SplineFitter.cpp
# - OpenMM reference CMAP patch evaluation:
#   https://github.com/openmm/openmm/blob/4768436/platforms/reference/src/SimTKReference/ReferenceCMAPTorsionIxn.cpp

# OpenMM CMAP bicubic coefficient matrix (CMAPTorsionForceImpl.cpp, wt[k+16*m]).
_CMAP_WT = (
    torch.tensor(
        [
            1, 0, -3, 2, 0, 0, 0, 0, -3, 0, 9, -6, 2, 0, -6, 4,
            0, 0, 0, 0, 0, 0, 0, 0, 3, 0, -9, 6, -2, 0, 6, -4,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, -6, 0, 0, -6, 4,
            0, 0, 3, -2, 0, 0, 0, 0, 0, 0, -9, 6, 0, 0, 6, -4,
            0, 0, 0, 0, 1, 0, -3, 2, -2, 0, 6, -4, 1, 0, -3, 2,
            0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 3, -2, 1, 0, -3, 2,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 2, 0, 0, 3, -2,
            0, 0, 0, 0, 0, 0, 3, -2, 0, 0, -6, 4, 0, 0, 3, -2,
            0, 1, -2, 1, 0, 0, 0, 0, 0, -3, 6, -3, 0, 2, -4, 2,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 3, -6, 3, 0, -2, 4, -2,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, 2, -2,
            0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 3, -3, 0, 0, -2, 2,
            0, 0, 0, 0, 0, 1, -2, 1, 0, -2, 4, -2, 0, 1, -2, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 2, -1, 0, 1, -2, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, -1, 1,
            0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 2, -2, 0, 0, -1, 1,
        ],
        dtype=torch.float64,
    ).view(16, 16).t()
)


def _solve_tridiagonal(
    a: torch.Tensor, b: torch.Tensor, c: torch.Tensor, rhs: torch.Tensor
) -> torch.Tensor:
    """Solve a tridiagonal linear system (Thomas algorithm)."""
    n = int(a.numel())
    if n == 0:
        return rhs.clone()
    if n == 1:
        return rhs / b

    gamma = torch.zeros(n, dtype=torch.float64)
    sol = torch.zeros(n, dtype=torch.float64)

    beta = b[0]
    sol[0] = rhs[0] / beta
    for i in range(1, n):
        gamma[i] = c[i - 1] / beta
        beta = b[i] - a[i] * gamma[i]
        sol[i] = (rhs[i] - a[i] * sol[i - 1]) / beta
    for i in range(n - 2, -1, -1):
        sol[i] = sol[i] - gamma[i + 1] * sol[i + 1]
    return sol


def _create_periodic_spline(x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    """Match OpenMM SplineFitter::createPeriodicSpline behavior."""
    n = int(x.numel())
    if n < 3:
        raise ValueError("Periodic spline requires at least 3 points")
    if y.numel() != n:
        raise ValueError("x/y size mismatch in periodic spline setup")

    a = torch.zeros(n - 1, dtype=torch.float64)
    b = torch.zeros(n - 1, dtype=torch.float64)
    c = torch.zeros(n - 1, dtype=torch.float64)
    rhs = torch.zeros(n - 1, dtype=torch.float64)

    a[0] = x[n - 1] - x[n - 2]
    b[0] = 2.0 * (x[1] - x[0] + x[n - 1] - x[n - 2])
    c[0] = x[1] - x[0]
    rhs[0] = 6.0 * (
        (y[1] - y[0]) / (x[1] - x[0]) - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2])
    )
    for i in range(1, n - 1):
        a[i] = x[i] - x[i - 1]
        b[i] = 2.0 * (x[i + 1] - x[i - 1])
        c[i] = x[i + 1] - x[i]
        rhs[i] = 6.0 * (
            (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1])
        )

    beta = a[0]
    alpha = c[n - 2]
    gamma = -b[0]

    ntri = n - 1
    b[0] = b[0] - gamma
    b[ntri - 1] = b[ntri - 1] - alpha * beta / gamma
    deriv = _solve_tridiagonal(a, b, c, rhs)

    u = torch.zeros(ntri, dtype=torch.float64)
    u[0] = gamma
    u[ntri - 1] = alpha
    z = _solve_tridiagonal(a, b, c, u)

    scale = (deriv[0] + beta * deriv[ntri - 1] / gamma) / (
        1.0 + z[0] + beta * z[ntri - 1] / gamma
    )
    deriv = deriv - scale * z

    out = torch.zeros(n, dtype=torch.float64)
    out[:ntri] = deriv
    out[ntri] = deriv[0]
    return out


def _eval_spline_derivative(
    x: torch.Tensor, y: torch.Tensor, second_deriv: torch.Tensor, t: torch.Tensor
) -> torch.Tensor:
    """Equivalent to OpenMM SplineFitter::evaluateSplineDerivative."""
    lower = 0
    upper = int(x.numel()) - 1
    while upper - lower > 1:
        middle = (upper + lower) // 2
        if x[middle] > t:
            upper = middle
        else:
            lower = middle
    delta = x[upper] - x[lower]
    a = (x[upper] - t) / delta
    b = (t - x[lower]) / delta
    return (-1.0 / delta) * (y[lower] - y[upper]) + (
        (1.0 - 3.0 * a * a) * second_deriv[lower] + (3.0 * b * b - 1.0) * second_deriv[upper]
    ) * delta / 6.0


def _calc_map_derivatives(energy: torch.Tensor, size: int) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Match OpenMM CMAPTorsionForceImpl::calcMapDerivatives."""
    d1 = torch.zeros(size * size, dtype=torch.float64)
    d2 = torch.zeros(size * size, dtype=torch.float64)
    d12 = torch.zeros(size * size, dtype=torch.float64)

    x = torch.arange(size + 1, dtype=torch.float64) * (_TWO_PI / float(size))

    # d/dphi
    for i in range(size):
        y = torch.zeros(size + 1, dtype=torch.float64)
        for j in range(size):
            y[j] = energy[j + size * i]
        y[size] = energy[size * i]
        sec = _create_periodic_spline(x, y)
        for j in range(size):
            d1[j + size * i] = _eval_spline_derivative(x, y, sec, x[j])

    # d/dpsi
    for i in range(size):
        y = torch.zeros(size + 1, dtype=torch.float64)
        for j in range(size):
            y[j] = energy[i + size * j]
        y[size] = energy[i]
        sec = _create_periodic_spline(x, y)
        for j in range(size):
            d2[i + size * j] = _eval_spline_derivative(x, y, sec, x[j])

    # d2/(dphi dpsi)
    for i in range(size):
        y = torch.zeros(size + 1, dtype=torch.float64)
        for j in range(size):
            y[j] = d2[j + size * i]
        y[size] = d2[size * i]
        sec = _create_periodic_spline(x, y)
        for j in range(size):
            d12[j + size * i] = _eval_spline_derivative(x, y, sec, x[j])

    return d1, d2, d12


def _calc_map_coefficients(map_values: torch.Tensor, size: int) -> torch.Tensor:
    """Create bicubic patch coefficients identical to OpenMM CMAP setup."""
    if map_values.numel() != size * size:
        raise ValueError(f"CMAP map size mismatch: expected {size*size}, got {map_values.numel()}")
    energy = map_values.detach().to(dtype=torch.float64, device="cpu").reshape(-1)
    d1, d2, d12 = _calc_map_derivatives(energy, size=size)
    coeff = torch.zeros((size * size, 16), dtype=torch.float64)
    delta = _TWO_PI / float(size)

    for i in range(size):
        i1 = (i + 1) % size
        for j in range(size):
            j1 = (j + 1) % size
            k = i + size * j
            rhs = torch.tensor(
                [
                    energy[k],
                    energy[i1 + size * j],
                    energy[i1 + size * j1],
                    energy[i + size * j1],
                    d1[k] * delta,
                    d1[i1 + size * j] * delta,
                    d1[i1 + size * j1] * delta,
                    d1[i + size * j1] * delta,
                    d2[k] * delta,
                    d2[i1 + size * j] * delta,
                    d2[i1 + size * j1] * delta,
                    d2[i + size * j1] * delta,
                    d12[k] * delta * delta,
                    d12[i1 + size * j] * delta * delta,
                    d12[i1 + size * j1] * delta * delta,
                    d12[i + size * j1] * delta * delta,
                ],
                dtype=torch.float64,
            )
            coeff[k] = _CMAP_WT @ rhs
    return coeff


class CMapTerm(nn.Module):
    """CMAP torsion correction implemented with Torch bicubic interpolation.

    The coefficient generation follows OpenMM's CMAP setup logic. The forward
    path is composed of standard Torch ops, so both gradients and Hessians are
    available through autograd.
    """

    def __init__(
        self,
        natom: int,
        cmap_type: torch.Tensor,
        cmap_i: torch.Tensor,
        cmap_j: torch.Tensor,
        cmap_k: torch.Tensor,
        cmap_l: torch.Tensor,
        cmap_m: torch.Tensor,
        cmap_resolution: Tuple[int, ...],
        cmap_maps: Tuple[torch.Tensor, ...],
        *,
        precomputed_size: torch.Tensor | None = None,
        precomputed_delta: torch.Tensor | None = None,
        precomputed_offset: torch.Tensor | None = None,
        precomputed_coeff: torch.Tensor | None = None,
    ):
        super().__init__()
        self.natom = int(natom)
        self.register_buffer("cmap_type", cmap_type.long())
        self.register_buffer("cmap_i", cmap_i.long())
        self.register_buffer("cmap_j", cmap_j.long())
        self.register_buffer("cmap_k", cmap_k.long())
        self.register_buffer("cmap_l", cmap_l.long())
        self.register_buffer("cmap_m", cmap_m.long())

        self._enabled = self.cmap_type.numel() > 0
        if not self._enabled:
            self.register_buffer("map_size", torch.zeros((0,), dtype=torch.int64))
            self.register_buffer("map_delta", torch.zeros((0,), dtype=torch.float64))
            self.register_buffer("map_offset", torch.zeros((0,), dtype=torch.int64))
            self.register_buffer("map_coeff", torch.zeros((0, 16), dtype=torch.float64))
            return

        # Use precomputed coefficients if available (avoids redundant computation
        # when AmberSystem already carries cmap_size/delta/offset/coeff).
        if precomputed_coeff is not None:
            self.register_buffer("map_size", precomputed_size)
            self.register_buffer("map_delta", precomputed_delta)
            self.register_buffer("map_offset", precomputed_offset)
            self.register_buffer("map_coeff", precomputed_coeff)
            return

        if len(cmap_resolution) != len(cmap_maps):
            raise ValueError(
                f"CMAP resolution/map count mismatch: {len(cmap_resolution)} vs {len(cmap_maps)}"
            )

        map_size = torch.tensor([int(x) for x in cmap_resolution], dtype=torch.int64)
        coeff_parts = []
        offsets = []
        cursor = 0
        for ngrid, table in zip(cmap_resolution, cmap_maps):
            offsets.append(cursor)
            c = _calc_map_coefficients(table, size=int(ngrid))
            coeff_parts.append(c)
            cursor += int(ngrid) * int(ngrid)

        coeff_cat = torch.cat(coeff_parts, dim=0)
        target_dtype = cmap_maps[0].dtype
        target_device = cmap_maps[0].device

        self.register_buffer("map_size", map_size.to(device=target_device))
        self.register_buffer(
            "map_delta",
            (torch.full_like(map_size, _TWO_PI, dtype=torch.float64) / map_size.to(torch.float64)).to(
                dtype=target_dtype, device=target_device
            ),
        )
        self.register_buffer("map_offset", torch.tensor(offsets, dtype=torch.int64, device=target_device))
        self.register_buffer("map_coeff", coeff_cat.to(dtype=target_dtype, device=target_device))

    def forward(self, coords: torch.Tensor) -> torch.Tensor:
        if not self._enabled:
            return coords.new_zeros(())

        p_i = coords[self.cmap_i]
        p_j = coords[self.cmap_j]
        p_k = coords[self.cmap_k]
        p_l = coords[self.cmap_l]
        p_m = coords[self.cmap_m]

        phi = _dihedral_angle(p_i, p_j, p_k, p_l)
        psi = _dihedral_angle(p_j, p_k, p_l, p_m)

        # OpenMM CMAP interpolation uses [0, 2pi) wrapped angles.
        ang_a = torch.remainder(phi + _TWO_PI, _TWO_PI)
        ang_b = torch.remainder(psi + _TWO_PI, _TWO_PI)

        tmap = self.cmap_type
        size = self.map_size[tmap]
        delta = self.map_delta[tmap].to(dtype=coords.dtype)

        u = ang_a / delta
        v = ang_b / delta
        s = torch.floor(u).to(torch.int64)
        t = torch.floor(v).to(torch.int64)
        s = torch.minimum(s, size - 1)
        t = torch.minimum(t, size - 1)

        da = u - s.to(dtype=coords.dtype)
        db = v - t.to(dtype=coords.dtype)

        patch = s + size * t
        coeff_row = self.map_offset[tmap] + patch
        coeff = self.map_coeff[coeff_row]
        coeff = coeff.to(dtype=coords.dtype).reshape(-1, 4, 4)

        # Horner evaluation in db, then da.
        dbu = db.unsqueeze(-1)
        poly_b = ((coeff[:, :, 3] * dbu + coeff[:, :, 2]) * dbu + coeff[:, :, 1]) * dbu + coeff[:, :, 0]
        e = ((poly_b[:, 3] * da + poly_b[:, 2]) * da + poly_b[:, 1]) * da + poly_b[:, 0]
        return torch.sum(e)

    def energy_force(self, coords: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        """Return CMAP energy and analytical force."""
        force = torch.zeros_like(coords)
        if not self._enabled:
            return coords.new_zeros(()), force

        p_i = coords[self.cmap_i]
        p_j = coords[self.cmap_j]
        p_k = coords[self.cmap_k]
        p_l = coords[self.cmap_l]
        p_m = coords[self.cmap_m]

        phi = _dihedral_angle(p_i, p_j, p_k, p_l)
        psi = _dihedral_angle(p_j, p_k, p_l, p_m)

        ang_a = torch.remainder(phi + _TWO_PI, _TWO_PI)
        ang_b = torch.remainder(psi + _TWO_PI, _TWO_PI)

        tmap = self.cmap_type
        size = self.map_size[tmap]
        delta = self.map_delta[tmap].to(dtype=coords.dtype)

        u = ang_a / delta
        v = ang_b / delta
        s = torch.floor(u).to(torch.int64)
        t = torch.floor(v).to(torch.int64)
        s = torch.minimum(s, size - 1)
        t = torch.minimum(t, size - 1)

        da = u - s.to(dtype=coords.dtype)
        db = v - t.to(dtype=coords.dtype)

        patch = s + size * t
        coeff_row = self.map_offset[tmap] + patch
        coeff = self.map_coeff[coeff_row].to(dtype=coords.dtype).reshape(-1, 4, 4)

        dbu = db.unsqueeze(-1)
        # Row-wise polynomial in db: Pm(db) (m=0..3)
        poly_b = ((coeff[:, :, 3] * dbu + coeff[:, :, 2]) * dbu + coeff[:, :, 1]) * dbu + coeff[:, :, 0]
        # Derivative wrt db: dPm/ddb
        dpoly_b = ((3.0 * coeff[:, :, 3] * dbu + 2.0 * coeff[:, :, 2]) * dbu + coeff[:, :, 1])

        # E = sum_m Pm(db) * da^m
        e = ((poly_b[:, 3] * da + poly_b[:, 2]) * da + poly_b[:, 1]) * da + poly_b[:, 0]
        energy = torch.sum(e)

        # dE/dda and dE/ddb
        dE_dda = (3.0 * poly_b[:, 3] * da + 2.0 * poly_b[:, 2]) * da + poly_b[:, 1]
        dE_ddb = ((dpoly_b[:, 3] * da + dpoly_b[:, 2]) * da + dpoly_b[:, 1]) * da + dpoly_b[:, 0]

        # da = phi/delta - floor(...), db = psi/delta - floor(...)
        # Away from grid boundaries, d(da)/dphi = 1/delta and d(db)/dpsi = 1/delta.
        dE_dphi = dE_dda / delta
        dE_dpsi = dE_ddb / delta

        _accumulate_dihedral_forces(
            force, coords, self.cmap_i, self.cmap_j, self.cmap_k, self.cmap_l, dE_dphi
        )
        _accumulate_dihedral_forces(
            force, coords, self.cmap_j, self.cmap_k, self.cmap_l, self.cmap_m, dE_dpsi
        )
        return energy, force
