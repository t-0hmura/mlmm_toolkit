from __future__ import annotations

from typing import Dict, Sequence, Tuple

import torch

from .constants import TWO_PI
from .system import AmberSystem
from .native.analytical_hessian import (
    bond_hessian_native,
    nonbonded_hessian_native,
    scatter_local_hessian_native,
)
from .native.loader import get_analytical_hessian_extension

Dual = Tuple[torch.Tensor, torch.Tensor, torch.Tensor]
DualVec3 = Tuple[Dual, Dual, Dual]
_NATIVE_ANALYTICAL_AVAILABLE: bool | None = None


def _has_native_analytical_ext() -> bool:
    global _NATIVE_ANALYTICAL_AVAILABLE
    if _NATIVE_ANALYTICAL_AVAILABLE is None:
        _NATIVE_ANALYTICAL_AVAILABLE = (
            get_analytical_hessian_extension(verbose=False, force_rebuild=False) is not None
        )
    return bool(_NATIVE_ANALYTICAL_AVAILABLE)


def _make_active_map(
    natom: int,
    active_atoms: Sequence[int],
    *,
    device: torch.device,
) -> torch.Tensor:
    amap = torch.full((int(natom),), -1, dtype=torch.int64, device=device)
    if len(active_atoms) == 0:
        return amap
    idx = torch.tensor([int(a) for a in active_atoms], dtype=torch.int64, device=device)
    vals = torch.arange(len(active_atoms), dtype=torch.int64, device=device)
    amap.index_put_((idx,), vals)
    return amap


def _slice_by_rank(x: torch.Tensor, rank: int, size: int) -> torch.Tensor:
    if size <= 1 or x.numel() == 0:
        return x
    return x[rank::size]


def _outer(a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
    return a.unsqueeze(-1) * b.unsqueeze(-2)


def _dual_vars(x: torch.Tensor) -> list[Dual]:
    t, m = int(x.shape[0]), int(x.shape[1])
    eye = torch.eye(m, dtype=x.dtype, device=x.device).unsqueeze(0).expand(t, m, m)
    h0 = x.new_zeros((t, m, m))
    return [(x[:, i], eye[:, i, :], h0) for i in range(m)]


def _d_add(a: Dual, b: Dual) -> Dual:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _d_sub(a: Dual, b: Dual) -> Dual:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _d_mul(a: Dual, b: Dual) -> Dual:
    av, ag, ah = a
    bv, bg, bh = b
    v = av * bv
    g = ag * bv.unsqueeze(-1) + bg * av.unsqueeze(-1)
    h = (
        ah * bv.unsqueeze(-1).unsqueeze(-1)
        + bh * av.unsqueeze(-1).unsqueeze(-1)
        + _outer(ag, bg)
        + _outer(bg, ag)
    )
    return v, g, h


def _d_inv(a: Dual) -> Dual:
    av, ag, ah = a
    inv = 1.0 / av
    inv2 = inv * inv
    inv3 = inv2 * inv
    v = inv
    g = -ag * inv2.unsqueeze(-1)
    h = (2.0 * _outer(ag, ag) - ah * av.unsqueeze(-1).unsqueeze(-1)) * inv3.unsqueeze(-1).unsqueeze(-1)
    return v, g, h


def _d_div(a: Dual, b: Dual) -> Dual:
    return _d_mul(a, _d_inv(b))


def _d_clamp_min(a: Dual, lo: float) -> Dual:
    av, ag, ah = a
    v = torch.clamp(av, min=float(lo))
    mask = av < float(lo)
    if bool(torch.any(mask)):
        g = ag.clone()
        h = ah.clone()
        g[mask] = 0.0
        h[mask] = 0.0
        return v, g, h
    return v, ag, ah


def _d_sqrt(a: Dual) -> Dual:
    av, ag, ah = a
    v = torch.sqrt(av)
    inv2v = 0.5 / v
    g = ag * inv2v.unsqueeze(-1)
    h = ah * inv2v.unsqueeze(-1).unsqueeze(-1) - _outer(ag, ag) * (0.25 / (v * v * v)).unsqueeze(-1).unsqueeze(-1)
    return v, g, h


def _d_atan2(y: Dual, x: Dual) -> Dual:
    yv, yg, yh = y
    xv, xg, xh = x
    d = xv * xv + yv * yv
    d = torch.clamp(d, min=1.0e-24)
    v = torch.atan2(yv, xv)

    n = xv.unsqueeze(-1) * yg - yv.unsqueeze(-1) * xg
    g = n / d.unsqueeze(-1)

    n_jac = (
        xv.unsqueeze(-1).unsqueeze(-1) * yh
        - yv.unsqueeze(-1).unsqueeze(-1) * xh
        + _outer(yg, xg)
        - _outer(xg, yg)
    )
    d_grad = 2.0 * xv.unsqueeze(-1) * xg + 2.0 * yv.unsqueeze(-1) * yg
    h = n_jac / d.unsqueeze(-1).unsqueeze(-1) - _outer(n, d_grad) / (d * d).unsqueeze(-1).unsqueeze(-1)
    h = 0.5 * (h + h.transpose(-1, -2))
    return v, g, h


def _dv_sub(a: DualVec3, b: DualVec3) -> DualVec3:
    return (_d_sub(a[0], b[0]), _d_sub(a[1], b[1]), _d_sub(a[2], b[2]))


def _dv_cross(a: DualVec3, b: DualVec3) -> DualVec3:
    c0 = _d_sub(_d_mul(a[1], b[2]), _d_mul(a[2], b[1]))
    c1 = _d_sub(_d_mul(a[2], b[0]), _d_mul(a[0], b[2]))
    c2 = _d_sub(_d_mul(a[0], b[1]), _d_mul(a[1], b[0]))
    return c0, c1, c2


def _dv_dot(a: DualVec3, b: DualVec3) -> Dual:
    return _d_add(_d_add(_d_mul(a[0], b[0]), _d_mul(a[1], b[1])), _d_mul(a[2], b[2]))


def _dv_norm(a: DualVec3, eps: float = 1.0e-24) -> Dual:
    s2 = _dv_dot(a, a)
    return _d_sqrt(_d_clamp_min(s2, float(eps)))


def _dual_angle_from_points(pi: DualVec3, pj: DualVec3, pk: DualVec3) -> Dual:
    u = _dv_sub(pi, pj)
    v = _dv_sub(pk, pj)
    c = _dv_dot(u, v)
    p = _dv_cross(u, v)
    s = _dv_norm(p)
    return _d_atan2(s, c)


def _dual_dihedral_from_points(p0: DualVec3, p1: DualVec3, p2: DualVec3, p3: DualVec3) -> Dual:
    v1 = _dv_sub(p0, p1)  # A-B
    v2 = _dv_sub(p2, p1)  # C-B
    v3 = _dv_sub(p2, p3)  # C-D

    c12 = _dv_cross(v1, v2)
    c23 = _dv_cross(v2, v3)
    x = _dv_dot(c12, c23)
    c = _dv_cross(c12, c23)
    v2n = _dv_norm(v2)
    y = _d_div(_dv_dot(v2, c), v2n)
    return _d_atan2(y, x)


def _local_dof_indices(atom_idx: torch.Tensor, active_map: torch.Tensor) -> torch.Tensor:
    t = int(atom_idx.shape[0])
    na = int(atom_idx.shape[1])
    dof_idx = torch.full((t, 3 * na), -1, dtype=torch.int64, device=active_map.device)
    act = active_map[atom_idx]
    for a in range(na):
        aa = act[:, a]
        mask = aa >= 0
        if not bool(torch.any(mask)):
            continue
        base = 3 * aa[mask]
        dof_idx[mask, 3 * a + 0] = base + 0
        dof_idx[mask, 3 * a + 1] = base + 1
        dof_idx[mask, 3 * a + 2] = base + 2
    return dof_idx


def _scatter_local_hessian(
    h2: torch.Tensor,
    local_h: torch.Tensor,
    dof_idx: torch.Tensor,
) -> None:
    h2.add_(
        scatter_local_hessian_native(
            local_h=local_h,
            dof_idx=dof_idx,
            ndof=int(h2.shape[0]),
        )
    )


def _filter_terms_with_active(active_map: torch.Tensor, atom_cols: Sequence[torch.Tensor]) -> torch.Tensor:
    keep = torch.zeros_like(atom_cols[0], dtype=torch.bool)
    for c in atom_cols:
        keep |= active_map[c] >= 0
    return keep


def _add_bond_hessian_closed_form(
    h2: torch.Tensor,
    system: AmberSystem,
    coords: torch.Tensor,
    active_map: torch.Tensor,
    rank: int,
    size: int,
) -> int:
    ii = _slice_by_rank(system.bond_i, rank, size)
    jj = _slice_by_rank(system.bond_j, rank, size)
    kk = _slice_by_rank(system.bond_k, rank, size)
    r0 = _slice_by_rank(system.bond_r0, rank, size)
    if ii.numel() == 0:
        return 0

    keep = _filter_terms_with_active(active_map, (ii, jj))
    ii, jj, kk, r0 = ii[keep], jj[keep], kk[keep], r0[keep]
    if ii.numel() == 0:
        return 0

    h2.add_(
        bond_hessian_native(
            coords=coords,
            bond_i=ii,
            bond_j=jj,
            bond_k=kk,
            bond_r0=r0,
            active_map=active_map,
            ndof=int(h2.shape[0]),
        )
    )
    return int(ii.numel())


def _add_nonbonded_pairset_closed_form(
    h2: torch.Tensor,
    system: AmberSystem,
    coords: torch.Tensor,
    active_map: torch.Tensor,
    pair_i: torch.Tensor,
    pair_j: torch.Tensor,
    inv_scee: torch.Tensor | None,
    inv_scnb: torch.Tensor | None,
    rank: int,
    size: int,
) -> int:
    ii = _slice_by_rank(pair_i, rank, size)
    jj = _slice_by_rank(pair_j, rank, size)
    scee = _slice_by_rank(inv_scee, rank, size) if inv_scee is not None else None
    scnb = _slice_by_rank(inv_scnb, rank, size) if inv_scnb is not None else None
    if ii.numel() == 0:
        return 0

    keep = _filter_terms_with_active(active_map, (ii, jj))
    ii, jj = ii[keep], jj[keep]
    if scee is not None:
        scee = scee[keep]
    if scnb is not None:
        scnb = scnb[keep]
    n = int(ii.numel())
    if n == 0:
        return 0

    h2.add_(
        nonbonded_hessian_native(
            coords=coords,
            charge=system.charge,
            atom_type=system.atom_type,
            lj_acoef=system.lj_acoef,
            lj_bcoef=system.lj_bcoef,
            hb_acoef=system.hb_acoef,
            hb_bcoef=system.hb_bcoef,
            nb_index=system.nb_index,
            pair_i=ii,
            pair_j=jj,
            inv_scee=scee,
            inv_scnb=scnb,
            active_map=active_map,
            ndof=int(h2.shape[0]),
        )
    )
    return n


def _add_angle_hessian_closed_form(
    h2: torch.Tensor,
    system: AmberSystem,
    coords: torch.Tensor,
    active_map: torch.Tensor,
    rank: int,
    size: int,
    chunk_terms: int = 256,
) -> int:
    ii = _slice_by_rank(system.angle_i, rank, size)
    jj = _slice_by_rank(system.angle_j, rank, size)
    kk = _slice_by_rank(system.angle_k, rank, size)
    ktheta = _slice_by_rank(system.angle_k0, rank, size)
    t0 = _slice_by_rank(system.angle_t0, rank, size)
    if ii.numel() == 0:
        return 0

    keep = _filter_terms_with_active(active_map, (ii, jj, kk))
    ii, jj, kk, ktheta, t0 = ii[keep], jj[keep], kk[keep], ktheta[keep], t0[keep]
    n = int(ii.numel())
    if n == 0:
        return 0

    for s in range(0, n, int(chunk_terms)):
        e = min(s + int(chunk_terms), n)
        ic, jc, kc = ii[s:e], jj[s:e], kk[s:e]
        ktc, t0c = ktheta[s:e], t0[s:e]
        local_xyz = torch.cat((coords[ic], coords[jc], coords[kc]), dim=1)  # [T,9]
        vars_ = _dual_vars(local_xyz)
        pi = (vars_[0], vars_[1], vars_[2])
        pj = (vars_[3], vars_[4], vars_[5])
        pk = (vars_[6], vars_[7], vars_[8])
        th_v, th_g, th_h = _dual_angle_from_points(pi, pj, pk)

        d1 = 2.0 * ktc * (th_v - t0c)
        d2 = 2.0 * ktc
        local_h = d2.unsqueeze(-1).unsqueeze(-1) * _outer(th_g, th_g) + d1.unsqueeze(-1).unsqueeze(-1) * th_h
        local_h = 0.5 * (local_h + local_h.transpose(-1, -2))

        atom_idx = torch.stack((ic, jc, kc), dim=1)
        dof_idx = _local_dof_indices(atom_idx, active_map)
        _scatter_local_hessian(h2, local_h, dof_idx)
    return n


def _add_dihedral_hessian_closed_form(
    h2: torch.Tensor,
    system: AmberSystem,
    coords: torch.Tensor,
    active_map: torch.Tensor,
    rank: int,
    size: int,
    chunk_terms: int = 192,
) -> int:
    ii = _slice_by_rank(system.dihed_i, rank, size)
    jj = _slice_by_rank(system.dihed_j, rank, size)
    kk = _slice_by_rank(system.dihed_k, rank, size)
    ll = _slice_by_rank(system.dihed_l, rank, size)
    force = _slice_by_rank(system.dihed_force, rank, size)
    period = _slice_by_rank(system.dihed_period, rank, size)
    phase = _slice_by_rank(system.dihed_phase, rank, size)
    if ii.numel() == 0:
        return 0

    keep = _filter_terms_with_active(active_map, (ii, jj, kk, ll))
    ii, jj, kk, ll = ii[keep], jj[keep], kk[keep], ll[keep]
    force, period, phase = force[keep], period[keep], phase[keep]
    n = int(ii.numel())
    if n == 0:
        return 0

    for s in range(0, n, int(chunk_terms)):
        e = min(s + int(chunk_terms), n)
        ic, jc, kc, lc = ii[s:e], jj[s:e], kk[s:e], ll[s:e]
        fc = force[s:e]
        pc = torch.abs(period[s:e])
        phc = phase[s:e]
        local_xyz = torch.cat((coords[ic], coords[jc], coords[kc], coords[lc]), dim=1)  # [T,12]
        vars_ = _dual_vars(local_xyz)
        p0 = (vars_[0], vars_[1], vars_[2])
        p1 = (vars_[3], vars_[4], vars_[5])
        p2 = (vars_[6], vars_[7], vars_[8])
        p3 = (vars_[9], vars_[10], vars_[11])
        phi_v, phi_g, phi_h = _dual_dihedral_from_points(p0, p1, p2, p3)

        d = pc * phi_v - phc
        e1 = -fc * pc * torch.sin(d)
        e2 = -fc * pc * pc * torch.cos(d)
        local_h = e2.unsqueeze(-1).unsqueeze(-1) * _outer(phi_g, phi_g) + e1.unsqueeze(-1).unsqueeze(-1) * phi_h
        local_h = 0.5 * (local_h + local_h.transpose(-1, -2))

        atom_idx = torch.stack((ic, jc, kc, lc), dim=1)
        dof_idx = _local_dof_indices(atom_idx, active_map)
        _scatter_local_hessian(h2, local_h, dof_idx)
    return n


def _cmap_patch_derivatives(coeff: torch.Tensor, da: torch.Tensor, db: torch.Tensor) -> tuple[torch.Tensor, ...]:
    # coeff: [T,4,4], da/db: [T]
    db2 = db * db
    da2 = da * da
    da3 = da2 * da

    p = (
        coeff[:, :, 0]
        + coeff[:, :, 1] * db.unsqueeze(-1)
        + coeff[:, :, 2] * db2.unsqueeze(-1)
        + coeff[:, :, 3] * (db2 * db).unsqueeze(-1)
    )
    pd = coeff[:, :, 1] + 2.0 * coeff[:, :, 2] * db.unsqueeze(-1) + 3.0 * coeff[:, :, 3] * db2.unsqueeze(-1)
    pdd = 2.0 * coeff[:, :, 2] + 6.0 * coeff[:, :, 3] * db.unsqueeze(-1)

    e = p[:, 0] + p[:, 1] * da + p[:, 2] * da2 + p[:, 3] * da3
    e_da = p[:, 1] + 2.0 * p[:, 2] * da + 3.0 * p[:, 3] * da2
    e_daa = 2.0 * p[:, 2] + 6.0 * p[:, 3] * da

    e_db = pd[:, 0] + pd[:, 1] * da + pd[:, 2] * da2 + pd[:, 3] * da3
    e_dbb = pdd[:, 0] + pdd[:, 1] * da + pdd[:, 2] * da2 + pdd[:, 3] * da3
    e_dab = pd[:, 1] + 2.0 * pd[:, 2] * da + 3.0 * pd[:, 3] * da2
    return e, e_da, e_db, e_daa, e_dbb, e_dab


def _add_cmap_hessian_closed_form(
    h2: torch.Tensor,
    system: AmberSystem,
    coords: torch.Tensor,
    active_map: torch.Tensor,
    rank: int,
    size: int,
    chunk_terms: int = 96,
) -> int:
    if int(system.cmap_type.numel()) == 0:
        return 0

    tmap = _slice_by_rank(system.cmap_type, rank, size)
    ia = _slice_by_rank(system.cmap_i, rank, size)
    ja = _slice_by_rank(system.cmap_j, rank, size)
    ka = _slice_by_rank(system.cmap_k, rank, size)
    la = _slice_by_rank(system.cmap_l, rank, size)
    ma = _slice_by_rank(system.cmap_m, rank, size)
    if ia.numel() == 0:
        return 0

    keep = _filter_terms_with_active(active_map, (ia, ja, ka, la, ma))
    tmap, ia, ja, ka, la, ma = tmap[keep], ia[keep], ja[keep], ka[keep], la[keep], ma[keep]
    n = int(ia.numel())
    if n == 0:
        return 0

    for s in range(0, n, int(chunk_terms)):
        e = min(s + int(chunk_terms), n)
        tc = tmap[s:e]
        ic, jc, kc, lc, mc = ia[s:e], ja[s:e], ka[s:e], la[s:e], ma[s:e]

        local_xyz = torch.cat((coords[ic], coords[jc], coords[kc], coords[lc], coords[mc]), dim=1)  # [T,15]
        vars_ = _dual_vars(local_xyz)
        pi = (vars_[0], vars_[1], vars_[2])
        pj = (vars_[3], vars_[4], vars_[5])
        pk = (vars_[6], vars_[7], vars_[8])
        pl = (vars_[9], vars_[10], vars_[11])
        pm = (vars_[12], vars_[13], vars_[14])

        phi_v, phi_g, phi_h = _dual_dihedral_from_points(pi, pj, pk, pl)
        psi_v, psi_g, psi_h = _dual_dihedral_from_points(pj, pk, pl, pm)

        size_map = system.cmap_size[tc]
        delta = system.cmap_delta[tc].to(dtype=coords.dtype)
        ang_a = torch.remainder(phi_v + TWO_PI, TWO_PI)
        ang_b = torch.remainder(psi_v + TWO_PI, TWO_PI)
        u = ang_a / delta
        v = ang_b / delta
        su = torch.floor(u).to(torch.int64)
        sv = torch.floor(v).to(torch.int64)
        su = torch.minimum(su, size_map - 1)
        sv = torch.minimum(sv, size_map - 1)
        da = u - su.to(dtype=coords.dtype)
        db = v - sv.to(dtype=coords.dtype)

        patch = su + size_map * sv
        coeff_row = system.cmap_offset[tc] + patch
        coeff = system.cmap_coeff[coeff_row].to(dtype=coords.dtype).reshape(-1, 4, 4)
        _, e_da, e_db, e_daa, e_dbb, e_dab = _cmap_patch_derivatives(coeff, da, db)

        inv_delta = 1.0 / delta
        inv_delta2 = inv_delta * inv_delta
        e_phi = e_da * inv_delta
        e_psi = e_db * inv_delta
        e_phiphi = e_daa * inv_delta2
        e_psipsi = e_dbb * inv_delta2
        e_phipsi = e_dab * inv_delta2

        local_h = (
            e_phiphi.unsqueeze(-1).unsqueeze(-1) * _outer(phi_g, phi_g)
            + e_psipsi.unsqueeze(-1).unsqueeze(-1) * _outer(psi_g, psi_g)
            + e_phipsi.unsqueeze(-1).unsqueeze(-1) * (_outer(phi_g, psi_g) + _outer(psi_g, phi_g))
            + e_phi.unsqueeze(-1).unsqueeze(-1) * phi_h
            + e_psi.unsqueeze(-1).unsqueeze(-1) * psi_h
        )
        local_h = 0.5 * (local_h + local_h.transpose(-1, -2))

        atom_idx = torch.stack((ic, jc, kc, lc, mc), dim=1)
        dof_idx = _local_dof_indices(atom_idx, active_map)
        _scatter_local_hessian(h2, local_h, dof_idx)
    return n


def build_analytical_hessian(
    system: AmberSystem,
    coords: torch.Tensor,
    active_atoms: Sequence[int],
    *,
    mpi_rank: int = 0,
    mpi_size: int = 1,
) -> tuple[torch.Tensor, Dict[str, int]]:
    """Build analytical Hessian without autograd Jacobian/JVP/Hessian.

    Term coverage:
    - Bond: closed-form
    - Angle: closed-form (manual second-order chain rule)
    - Dihedral: closed-form (manual second-order chain rule)
    - CMAP: closed-form (manual second-order chain rule + bicubic patch derivatives)
    - Nonbonded + 1-4 (including NB_INDEX<0 HB 12-10): closed-form
    """
    coords_work = coords.detach()
    system_work = system.to(device=coords_work.device, dtype=coords_work.dtype)

    active_list = [int(a) for a in active_atoms]
    if not _has_native_analytical_ext():
        raise RuntimeError(
            "Native analytical Hessian extension is required. "
            "Automatic torch fallback has been removed."
        )
    active_map = _make_active_map(system_work.natom, active_list, device=coords_work.device)
    ndof = int(3 * len(active_list))
    h2 = torch.zeros((ndof, ndof), dtype=coords_work.dtype, device=coords_work.device)
    nbond = _add_bond_hessian_closed_form(
        h2=h2,
        system=system_work,
        coords=coords_work,
        active_map=active_map,
        rank=int(mpi_rank),
        size=int(mpi_size),
    )
    nangle = _add_angle_hessian_closed_form(
        h2=h2,
        system=system_work,
        coords=coords_work,
        active_map=active_map,
        rank=int(mpi_rank),
        size=int(mpi_size),
    )
    ndihed = _add_dihedral_hessian_closed_form(
        h2=h2,
        system=system_work,
        coords=coords_work,
        active_map=active_map,
        rank=int(mpi_rank),
        size=int(mpi_size),
    )
    ncmap = _add_cmap_hessian_closed_form(
        h2=h2,
        system=system_work,
        coords=coords_work,
        active_map=active_map,
        rank=int(mpi_rank),
        size=int(mpi_size),
    )
    n_nb = _add_nonbonded_pairset_closed_form(
        h2=h2,
        system=system_work,
        coords=coords_work,
        active_map=active_map,
        pair_i=system_work.pair_i,
        pair_j=system_work.pair_j,
        inv_scee=None,
        inv_scnb=None,
        rank=int(mpi_rank),
        size=int(mpi_size),
    )
    n_14 = _add_nonbonded_pairset_closed_form(
        h2=h2,
        system=system_work,
        coords=coords_work,
        active_map=active_map,
        pair_i=system_work.pair14_i,
        pair_j=system_work.pair14_j,
        inv_scee=system_work.pair14_inv_scee,
        inv_scnb=system_work.pair14_inv_scnb,
        rank=int(mpi_rank),
        size=int(mpi_size),
    )

    meta = {
        "bond_pairs_used": int(nbond),
        "angle_terms_used": int(nangle),
        "dihedral_terms_used": int(ndihed),
        "cmap_terms_used": int(ncmap),
        "nonbond_pairs_used": int(n_nb),
        "nonbond14_pairs_used": int(n_14),
        "analytical_force_evals": 0,
    }
    return h2, meta
