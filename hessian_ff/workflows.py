from __future__ import annotations

from dataclasses import dataclass
import time
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Union

import torch

from .analytical_hessian import build_analytical_hessian
from .forcefield import ForceFieldTorch
from .loaders import load_coords, load_system
from .system import AmberSystem
from .terms.nonbonded import NonbondedTerm

PathLike = Union[str, Path]
CoordsLike = Union[PathLike, torch.Tensor, Any]
CoordsBatchLike = Union[torch.Tensor, Sequence[CoordsLike]]


# ---------------------------------------------------------------------------
# Hessian mode synonym map
# ---------------------------------------------------------------------------
_HESSIAN_MODE_MAP = {
    "analytical": "analytical",
    "forcejacobian": "analytical",
    "force_jacobian": "analytical",
    "customautograd": "analytical",
    "custom_autograd": "analytical",
    "autograd": "autograd",
    "finitedifferenceforce": "fd",
    "finite_difference_force": "fd",
    "finitedifference": "fd",
    "finite_difference": "fd",
    "fd": "fd",
}


# ---------------------------------------------------------------------------
# Runtime cache
# ---------------------------------------------------------------------------
@dataclass
class _RuntimeEntry:
    system: AmberSystem
    ff: ForceFieldTorch
    coords_buffer: Optional[torch.Tensor] = None


_RUNTIME_CACHE: Dict[tuple[str, str, str, int, bool], _RuntimeEntry] = {}


def clear_runtime_cache() -> None:
    """Clear cached runtime objects (system/force-field/coord buffer)."""
    _RUNTIME_CACHE.clear()


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def _dtype_from_double(double: bool) -> torch.dtype:
    return torch.float64 if bool(double) else torch.float32


def _precision_info(double: bool) -> Dict[str, Any]:
    dtype = _dtype_from_double(double)
    if dtype == torch.float64:
        dtype_name = "float64"
    elif dtype == torch.float32:
        dtype_name = "float32"
    else:
        dtype_name = str(dtype)
    return {"double": bool(double), "torch_dtype": dtype_name}


def _configure_torch_threads(num_threads: Optional[int]) -> Optional[int]:
    if num_threads is None:
        return None
    n = int(num_threads)
    if n < 1:
        raise ValueError(f"num_threads must be >=1, got {num_threads}")
    import os

    for env_key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ[env_key] = str(n)
    torch.set_num_threads(n)
    return n


def _energy_terms_to_float(energy_terms: Dict[str, torch.Tensor]) -> Dict[str, float]:
    return {k: float(v.detach().cpu()) for k, v in energy_terms.items() if k.startswith("E_")}

# ---------------------------------------------------------------------------
# MPI helpers (unified)
# ---------------------------------------------------------------------------
def _get_mpi_context(
    mpi: bool,
) -> tuple[Optional[Any], Optional[Any], int, int]:
    if not mpi:
        return None, None, 0, 1
    try:
        from mpi4py import MPI  # type: ignore
    except Exception as e:
        raise RuntimeError(
            "mpi=True requires mpi4py and a working MPI runtime (libmpi). "
            f"Original import error: {e}"
        ) from e
    comm = MPI.COMM_WORLD
    return MPI, comm, int(comm.Get_rank()), int(comm.Get_size())


def _mpi_reduce_tensor(local: torch.Tensor, mpi_obj: Any, comm: Any) -> torch.Tensor:
    """Allreduce a tensor (scalar, vector, or matrix) with MPI SUM."""
    if comm is None:
        return local
    arr_local = local.detach().cpu().contiguous().reshape(-1).numpy()
    arr_global = arr_local.copy()
    comm.Allreduce(arr_local, arr_global, op=mpi_obj.SUM)
    return torch.from_numpy(arr_global).reshape(local.shape).to(dtype=local.dtype)


# ---------------------------------------------------------------------------
# Result dict builder (eliminates repeated metadata assembly)
# ---------------------------------------------------------------------------
def _build_result_meta(
    *,
    double: bool,
    num_threads: Optional[int] = None,
    mpi_rank: int = 0,
    mpi_size: int = 1,
    **extra: Any,
) -> Dict[str, Any]:
    result: Dict[str, Any] = dict(extra)
    result["eval_backend"] = "torch"
    if num_threads is not None:
        result["num_threads"] = int(num_threads)
    result["mpi_enabled"] = bool(mpi_size > 1)
    result["mpi_rank"] = int(mpi_rank)
    result["mpi_size"] = int(mpi_size)
    result.update(_precision_info(double))
    return result


# ---------------------------------------------------------------------------
# Runtime loading
# ---------------------------------------------------------------------------
def _load_runtime(
    prmtop: PathLike,
    coords: CoordsLike,
    device: Union[str, torch.device],
    double: bool,
    requires_grad: bool = False,
    nonbonded_cpu_fast: bool = True,
) -> tuple[AmberSystem, torch.Tensor, ForceFieldTorch]:
    dtype = _dtype_from_double(double)
    dev = str(torch.device(device))
    cache_key = (
        str(Path(prmtop).resolve()),
        dev,
        str(dtype),
        bool(nonbonded_cpu_fast),
    )

    cached = _RUNTIME_CACHE.get(cache_key)
    if cached is None:
        system = load_system(prmtop, device=device).to(dtype=dtype)
        ff = ForceFieldTorch(
            system,
            nonbonded_cpu_fast=nonbonded_cpu_fast,
        )
        entry = _RuntimeEntry(system=system, ff=ff, coords_buffer=None)
        _RUNTIME_CACHE[cache_key] = entry
    else:
        entry = cached
        system = entry.system
        ff = entry.ff

    xyz_loaded = load_coords(coords, natom=system.natom, device=device, dtype=dtype)
    if requires_grad:
        xyz = xyz_loaded.clone().detach().requires_grad_(True)
        return system, xyz, ff

    buf = entry.coords_buffer
    if (
        buf is None
        or buf.shape != xyz_loaded.shape
        or buf.dtype != xyz_loaded.dtype
        or buf.device != xyz_loaded.device
    ):
        entry.coords_buffer = xyz_loaded.clone().detach()
    else:
        entry.coords_buffer.copy_(xyz_loaded)
    xyz = entry.coords_buffer
    return system, xyz, ff


# ---------------------------------------------------------------------------
# MPI-distributed energy/force on CPU
# ---------------------------------------------------------------------------
def _dist_energy_force_cpu(
    system: AmberSystem,
    coords: torch.Tensor,
    *,
    force_calc_mode: str,
    mpi_obj: Any,
    mpi_comm: Any,
    mpi_rank: int,
    mpi_size: int,
) -> tuple[Dict[str, torch.Tensor], torch.Tensor]:
    if str(force_calc_mode).strip().lower() != "analytical":
        raise ValueError("mpi=True for E/F currently requires force_calc_mode='Analytical'")
    if coords.device.type != "cpu":
        raise ValueError("Distributed CPU E/F path requires CPU tensor")

    e_bond = coords.new_zeros(())
    e_angle = coords.new_zeros(())
    e_dihed = coords.new_zeros(())
    e_cmap = coords.new_zeros(())
    force_bonded = torch.zeros_like(coords)

    if mpi_rank == 0:
        ff_root = ForceFieldTorch(system)
        e_bond, f_bond = ff_root.bond.energy_force(coords)
        e_angle, f_angle = ff_root.angle.energy_force(coords)
        e_dihed, f_dihed = ff_root.dihedral.energy_force(coords)
        e_cmap, f_cmap = ff_root.cmap.energy_force(coords)
        force_bonded = f_bond + f_angle + f_dihed + f_cmap

    def _shard(x: torch.Tensor) -> torch.Tensor:
        if mpi_size <= 1 or x.numel() == 0:
            return x
        return x[mpi_rank::mpi_size]

    nb_term = NonbondedTerm(
        natom=system.natom,
        charge=system.charge,
        atom_type=system.atom_type,
        lj_acoef=system.lj_acoef,
        lj_bcoef=system.lj_bcoef,
        hb_acoef=system.hb_acoef,
        hb_bcoef=system.hb_bcoef,
        nb_index=system.nb_index,
        pair_i=_shard(system.pair_i),
        pair_j=_shard(system.pair_j),
        pair14_i=_shard(system.pair14_i),
        pair14_j=_shard(system.pair14_j),
        pair14_inv_scee=_shard(system.pair14_inv_scee),
        pair14_inv_scnb=_shard(system.pair14_inv_scnb),
    )
    nb_e, f_nb = nb_term.energy_force(coords)

    local_force = force_bonded + f_nb
    local_terms: Dict[str, torch.Tensor] = {
        "E_bond": e_bond,
        "E_angle": e_angle,
        "E_dihedral": e_dihed,
        "E_cmap": e_cmap,
        "E_coul": nb_e.coulomb,
        "E_lj": nb_e.lj,
        "E_coul14": nb_e.coulomb14,
        "E_lj14": nb_e.lj14,
    }

    if mpi_size > 1:
        force_global = _mpi_reduce_tensor(local_force, mpi_obj=mpi_obj, comm=mpi_comm)
        terms_global = {
            k: _mpi_reduce_tensor(v, mpi_obj=mpi_obj, comm=mpi_comm)
            for k, v in local_terms.items()
        }
    else:
        force_global = local_force
        terms_global = local_terms

    e_nb_total = terms_global["E_coul"] + terms_global["E_lj"] + terms_global["E_coul14"] + terms_global["E_lj14"]
    e_total = terms_global["E_bond"] + terms_global["E_angle"] + terms_global["E_dihedral"] + terms_global["E_cmap"] + e_nb_total
    out = dict(terms_global)
    out["E_total"] = e_total
    out["E_nonbonded_total"] = e_nb_total
    return out, force_global


# ---------------------------------------------------------------------------
# Core energy/force computation (shared by torch_energy and torch_force)
# ---------------------------------------------------------------------------
def _compute_energy_force(
    prmtop: PathLike,
    coords: CoordsLike,
    *,
    device: Union[str, torch.device],
    double: bool,
    force_calc_mode: str,
    with_force: bool,
    num_threads: Optional[int],
    mpi: bool,
 ) -> tuple[Dict[str, Any], Optional[torch.Tensor], Optional[int], int, int]:
    """Shared computation for torch_energy/torch_force.

    Returns (energy_terms_float, force_or_None, used_threads, mpi_rank, mpi_size).
    """
    used_threads = _configure_torch_threads(num_threads)
    mpi_obj, mpi_comm, mpi_rank, mpi_size = _get_mpi_context(mpi)

    dev = torch.device(device)
    if mpi_size > 1 and dev.type != "cpu":
        raise ValueError("mpi=True currently supports CPU execution only")

    force_mode_lower = str(force_calc_mode).strip().lower()
    system, xyz, ff = _load_runtime(
        prmtop=prmtop,
        coords=coords,
        device=device,
        double=double,
        requires_grad=(with_force and force_mode_lower == "autograd"),
        nonbonded_cpu_fast=not (torch.device(device).type == "cpu" and force_mode_lower == "autograd"),
    )

    if mpi_size > 1:
        out, force_tensor = _dist_energy_force_cpu(
            system=system,
            coords=xyz,
            force_calc_mode=force_calc_mode,
            mpi_obj=mpi_obj,
            mpi_comm=mpi_comm,
            mpi_rank=mpi_rank,
            mpi_size=mpi_size,
        )
    elif with_force:
        out, force_tensor = ff.energy_force(coords=xyz, force_calc_mode=force_calc_mode)
    else:
        out = ff(xyz)
        force_tensor = None

    energy_dict = _energy_terms_to_float(out)
    return energy_dict, force_tensor, used_threads, mpi_rank, mpi_size


# ---------------------------------------------------------------------------
# Public API: system_summary
# ---------------------------------------------------------------------------
def system_summary(prmtop: PathLike, device: Union[str, torch.device] = "cpu") -> Dict[str, int]:
    """Return basic term/pair counts for a prmtop."""
    system = load_system(prmtop, device=device)
    return {
        "n_atom": int(system.natom),
        "n_bond": int(system.bond_i.numel()),
        "n_angle": int(system.angle_i.numel()),
        "n_dihedral": int(system.dihed_i.numel()),
        "n_cmap": int(system.cmap_type.numel()),
        "n_pair_general": int(system.pair_i.numel()),
        "n_pair_14": int(system.pair14_i.numel()),
    }


# ---------------------------------------------------------------------------
# Public API: torch_energy
# ---------------------------------------------------------------------------
def torch_energy(
    prmtop: PathLike,
    coords: CoordsLike,
    device: Union[str, torch.device] = "cpu",
    with_grad: bool = False,
    double: bool = True,
    force_calc_mode: str = "Analytical",
    num_threads: Optional[int] = None,
    mpi: bool = False,
) -> Dict[str, Any]:
    """Evaluate Torch MM energies, with optional gradient norm output."""
    energy_dict, force_tensor, used_threads, mpi_rank, mpi_size = _compute_energy_force(
        prmtop=prmtop,
        coords=coords,
        device=device,
        double=double,
        force_calc_mode=force_calc_mode,
        with_force=with_grad,
        num_threads=num_threads,
        mpi=mpi,
    )
    result = _build_result_meta(
        double=double,
        num_threads=used_threads,
        mpi_rank=mpi_rank,
        mpi_size=mpi_size,
        **energy_dict,
    )
    if with_grad and force_tensor is not None:
        grad = (-force_tensor).detach().cpu()
        result["grad_shape"] = [int(grad.shape[0]), int(grad.shape[1])]
        result["grad_norm"] = float(torch.linalg.norm(grad))
        result["force_calc_mode"] = str(force_calc_mode)
    return result


# ---------------------------------------------------------------------------
# Public API: torch_force
# ---------------------------------------------------------------------------
def torch_force(
    prmtop: PathLike,
    coords: CoordsLike,
    device: Union[str, torch.device] = "cpu",
    double: bool = True,
    force_calc_mode: str = "Analytical",
    save_force: Optional[PathLike] = None,
    num_threads: Optional[int] = None,
    mpi: bool = False,
) -> Dict[str, Any]:
    """Compute force (negative gradient) from Torch energy."""
    energy_dict, force_tensor, used_threads, mpi_rank, mpi_size = _compute_energy_force(
        prmtop=prmtop,
        coords=coords,
        device=device,
        double=double,
        force_calc_mode=force_calc_mode,
        with_force=True,
        num_threads=num_threads,
        mpi=mpi,
    )
    force = force_tensor.detach().cpu()
    e_total = energy_dict.get("E_total", energy_dict.get("E_total_kcalmol", 0.0))

    result = _build_result_meta(
        double=double,
        num_threads=used_threads,
        mpi_rank=mpi_rank,
        mpi_size=mpi_size,
        E_total_kcalmol=float(e_total) if not isinstance(e_total, float) else e_total,
        force_shape=[int(force.shape[0]), int(force.shape[1])],
        force_norm=float(torch.linalg.norm(force)),
        force_maxabs=float(torch.max(torch.abs(force))),
        force_calc_mode=str(force_calc_mode),
    )
    if save_force is not None:
        save_path = Path(save_force)
        torch.save(force, save_path)
        result["force_file"] = str(save_path)
    return result


# ---------------------------------------------------------------------------
# Core batch computation (shared by torch_energy_batch / torch_force_batch)
# ---------------------------------------------------------------------------
def _prepare_batch_ff(
    prmtop: PathLike,
    coords_batch: CoordsBatchLike,
    *,
    device: Union[str, torch.device],
    double: bool,
    force_calc_mode: str,
    with_force: bool,
) -> tuple[AmberSystem, torch.Tensor, ForceFieldTorch]:
    dtype = _dtype_from_double(double)
    system = load_system(prmtop, device=device).to(dtype=dtype)
    if isinstance(coords_batch, torch.Tensor):
        xb = coords_batch.to(device=device, dtype=dtype)
    else:
        frames = [load_coords(p, natom=system.natom, device=device, dtype=dtype) for p in coords_batch]
        if not frames:
            raise ValueError("coords_batch is empty")
        xb = torch.stack(frames, dim=0)
    if xb.ndim != 3:
        raise ValueError(f"coords_batch must have shape [B,N,3], got {tuple(xb.shape)}")
    if int(xb.shape[1]) != int(system.natom) or int(xb.shape[2]) != 3:
        raise ValueError(
            f"coords_batch shape mismatch: expected [B,{system.natom},3], got {tuple(xb.shape)}"
        )
    ff = ForceFieldTorch(
        system,
        nonbonded_cpu_fast=not (
            torch.device(device).type == "cpu"
            and with_force
            and str(force_calc_mode).strip().lower() == "autograd"
        ),
    )
    return system, xb, ff


# ---------------------------------------------------------------------------
# Public API: torch_energy_batch
# ---------------------------------------------------------------------------
def torch_energy_batch(
    prmtop: PathLike,
    coords_batch: CoordsBatchLike,
    device: Union[str, torch.device] = "cpu",
    with_grad: bool = False,
    double: bool = True,
    force_calc_mode: str = "Analytical",
    batch_mode: str = "vmap",
    microbatch_size: Optional[int] = None,
    num_threads: Optional[int] = None,
) -> Dict[str, Any]:
    """Evaluate energies for batched coordinates [B,N,3]."""
    used_threads = _configure_torch_threads(num_threads)
    system, xb, ff = _prepare_batch_ff(
        prmtop, coords_batch,
        device=device, double=double,
        force_calc_mode=force_calc_mode, with_force=with_grad,
    )

    t0 = time.perf_counter()
    if with_grad:
        out, force = ff.energy_force_batch(
            xb, force_calc_mode=force_calc_mode,
            batch_mode=batch_mode, microbatch_size=microbatch_size,
        )
        grad = (-force).detach().cpu()
    else:
        out = ff.forward_batch(xb, batch_mode=batch_mode, microbatch_size=microbatch_size)
    elapsed = time.perf_counter() - t0

    result: Dict[str, Any] = {
        "batch_size": int(xb.shape[0]),
        "n_atom": int(system.natom),
        "batch_mode": str(batch_mode),
        "elapsed_s": float(elapsed),
        "energy_terms_kcalmol": {k: out[k].detach().cpu().tolist() for k in out},
    }
    if with_grad:
        result["grad_shape"] = [int(grad.shape[0]), int(grad.shape[1]), int(grad.shape[2])]
        result["grad_norm"] = float(torch.linalg.norm(grad))
        result["force_calc_mode"] = str(force_calc_mode)
    if microbatch_size is not None:
        result["microbatch_size"] = int(microbatch_size)
    if used_threads is not None:
        result["num_threads"] = int(used_threads)
    result.update(_precision_info(double))
    return result


# ---------------------------------------------------------------------------
# Public API: torch_force_batch
# ---------------------------------------------------------------------------
def torch_force_batch(
    prmtop: PathLike,
    coords_batch: CoordsBatchLike,
    device: Union[str, torch.device] = "cpu",
    double: bool = True,
    force_calc_mode: str = "Analytical",
    batch_mode: str = "vmap",
    microbatch_size: Optional[int] = None,
    save_force: Optional[PathLike] = None,
    num_threads: Optional[int] = None,
) -> Dict[str, Any]:
    """Compute batched force tensors for coordinates [B,N,3]."""
    used_threads = _configure_torch_threads(num_threads)
    system, xb, ff = _prepare_batch_ff(
        prmtop, coords_batch,
        device=device, double=double,
        force_calc_mode=force_calc_mode, with_force=True,
    )

    t0 = time.perf_counter()
    out, force_tensor = ff.energy_force_batch(
        xb, force_calc_mode=force_calc_mode,
        batch_mode=batch_mode, microbatch_size=microbatch_size,
    )
    elapsed = time.perf_counter() - t0
    force = force_tensor.detach().cpu()

    result: Dict[str, Any] = {
        "batch_size": int(xb.shape[0]),
        "n_atom": int(system.natom),
        "batch_mode": str(batch_mode),
        "elapsed_s": float(elapsed),
        "E_total_kcalmol": out["E_total"].detach().cpu().tolist(),
        "force_shape": [int(force.shape[0]), int(force.shape[1]), int(force.shape[2])],
        "force_norm": float(torch.linalg.norm(force)),
        "force_maxabs": float(torch.max(torch.abs(force))),
        "force_calc_mode": str(force_calc_mode),
    }
    if microbatch_size is not None:
        result["microbatch_size"] = int(microbatch_size)
    if save_force is not None:
        save_path = Path(save_force)
        torch.save(force, save_path)
        result["force_file"] = str(save_path)
    if used_threads is not None:
        result["num_threads"] = int(used_threads)
    result.update(_precision_info(double))
    return result


# ---------------------------------------------------------------------------
# Public API: torch_hessian
# ---------------------------------------------------------------------------
def _normalize_active_atoms(natom: int, active_atoms: Sequence[int]) -> list[int]:
    out: list[int] = []
    seen: set[int] = set()
    for a in active_atoms:
        ia = int(a)
        if ia < 0 or ia >= natom:
            raise ValueError(f"active atom index out of range: {ia} (natom={natom})")
        if ia in seen:
            continue
        seen.add(ia)
        out.append(ia)
    if not out:
        raise ValueError("active atom list is empty")
    return out


def torch_hessian(
    prmtop: PathLike,
    coords: CoordsLike,
    device: Union[str, torch.device] = "cpu",
    double: bool = True,
    hessian_calc_mode: str = "FiniteDifferenceForce",
    force_calc_mode: str = "Analytical",
    hessian_delta: float = 1.0e-4,
    fd_column_batch: int = 16,
    partial_hessian: bool = True,
    active_atoms: Optional[Sequence[int]] = None,
    save_hessian: Optional[PathLike] = None,
    num_threads: Optional[int] = None,
    mpi: bool = False,
) -> Dict[str, Any]:
    """Compute Hessian by one of three modes: Analytical, Autograd, or FiniteDifferenceForce."""
    used_threads = _configure_torch_threads(num_threads)
    mpi_obj, mpi_comm, mpi_rank, mpi_size = _get_mpi_context(mpi)

    dev = torch.device(device)
    if mpi_size > 1 and dev.type != "cpu":
        raise ValueError("mpi=True currently supports CPU execution only")

    mode = _HESSIAN_MODE_MAP.get(str(hessian_calc_mode).strip().lower())
    if mode is None:
        raise ValueError(f"Unknown hessian_calc_mode: {hessian_calc_mode!r}")

    force_mode = str(force_calc_mode).strip().lower()
    need_diff_nonbonded = dev.type == "cpu" and (mode == "autograd" or force_mode == "autograd")

    system, xyz, ff = _load_runtime(
        prmtop=prmtop, coords=coords, device=device,
        double=double, requires_grad=False,
        nonbonded_cpu_fast=not need_diff_nonbonded,
    )

    natom = int(system.natom)
    if partial_hessian:
        if active_atoms is None:
            raise ValueError("partial_hessian=true requires active_atoms")
        active_list = _normalize_active_atoms(natom, active_atoms)
    else:
        if active_atoms is not None:
            raise ValueError("active_atoms is only valid when partial_hessian=true")
        active_list = list(range(natom))

    active_idx = torch.tensor(active_list, dtype=torch.int64, device=xyz.device)
    x_active_base = xyz.index_select(0, active_idx).clone().detach()
    ndof = int(x_active_base.numel())

    n_force_eval: Optional[int] = None
    analytical_meta: Optional[Dict[str, int]] = None

    if mode == "analytical":
        if force_mode != "analytical":
            raise ValueError("hessian_calc_mode='Analytical' requires force_calc_mode='Analytical'")
        t0 = time.perf_counter()
        h_local, analytical_meta_local = build_analytical_hessian(
            system=system, coords=xyz, active_atoms=active_list,
            mpi_rank=int(mpi_rank), mpi_size=int(mpi_size),
        )
        elapsed_local = time.perf_counter() - t0
        if mpi_size > 1:
            h2 = _mpi_reduce_tensor(h_local, mpi_obj=mpi_obj, comm=mpi_comm)
            elapsed = float(mpi_comm.allreduce(elapsed_local, op=mpi_obj.MAX))
            meta_sum_keys = {
                "bond_pairs_used", "angle_terms_used", "dihedral_terms_used",
                "cmap_terms_used", "nonbond_pairs_used", "nonbond14_pairs_used",
                "analytical_force_evals",
            }
            analytical_meta = {}
            for k, v in analytical_meta_local.items():
                op = mpi_obj.SUM if k in meta_sum_keys else mpi_obj.MAX
                analytical_meta[k] = int(mpi_comm.allreduce(int(v), op=op))
        else:
            h2 = h_local
            elapsed = elapsed_local
            analytical_meta = analytical_meta_local
        h2 = (0.5 * (h2 + h2.T)).detach().cpu()
        n_force_eval = int(analytical_meta.get("analytical_force_evals", 0))

    elif mode == "autograd":
        if mpi_size > 1:
            raise ValueError("hessian_calc_mode='Autograd' does not support mpi=True")

        def _energy_fn(x_sub: torch.Tensor) -> torch.Tensor:
            x_full = xyz.index_copy(0, active_idx, x_sub)
            return ff(x_full)["E_total"]

        x_active = x_active_base.requires_grad_(True)
        t0 = time.perf_counter()
        h4 = torch.autograd.functional.hessian(_energy_fn, x_active, vectorize=False)
        elapsed = time.perf_counter() - t0
        h2 = h4.reshape(ndof, ndof).detach().cpu()

    else:  # mode == "fd"
        delta = float(hessian_delta)
        fd_batch_cols = max(1, int(fd_column_batch))
        t0 = time.perf_counter()
        if mpi_size > 1:
            h2 = torch.zeros((ndof, ndof), dtype=xyz.dtype, device="cpu")
            col_iter = range(mpi_rank, ndof, mpi_size)
        else:
            h2 = torch.empty((ndof, ndof), dtype=xyz.dtype, device="cpu")
            col_iter = range(ndof)
        cols = list(col_iter)
        n_force_eval_local = 0
        for start in range(0, len(cols), fd_batch_cols):
            sub_cols = cols[start : start + fd_batch_cols]
            bsz = len(sub_cols)
            xb = xyz.unsqueeze(0).repeat(2 * bsz, 1, 1).clone()
            for bi, col in enumerate(sub_cols):
                a = col // 3
                c = col % 3
                atom_idx = int(active_list[a])
                xb[2 * bi, atom_idx, c] = xb[2 * bi, atom_idx, c] + delta
                xb[2 * bi + 1, atom_idx, c] = xb[2 * bi + 1, atom_idx, c] - delta
            _, force_b = ff.energy_force_batch(
                xb, force_calc_mode=force_calc_mode, batch_mode="loop",
                microbatch_size=max(1, min(2 * bsz, 64)),
            )
            for bi, col in enumerate(sub_cols):
                fp = force_b[2 * bi].index_select(0, active_idx).reshape(-1)
                fm = force_b[2 * bi + 1].index_select(0, active_idx).reshape(-1)
                h2[:, col] = (-(fp - fm) / (2.0 * delta)).detach().cpu()
            n_force_eval_local += 2 * bsz
        elapsed_local = time.perf_counter() - t0
        if mpi_size > 1:
            h2 = _mpi_reduce_tensor(h2, mpi_obj=mpi_obj, comm=mpi_comm)
            elapsed = float(mpi_comm.allreduce(elapsed_local, op=mpi_obj.MAX))
            n_force_eval = int(mpi_comm.allreduce(n_force_eval_local, op=mpi_obj.SUM))
        else:
            elapsed = elapsed_local
            n_force_eval = n_force_eval_local

    hasym = torch.max(torch.abs(h2 - h2.T))

    out: Dict[str, Any] = {
        "hessian_calc_mode": str(hessian_calc_mode),
        "force_calc_mode": str(force_calc_mode),
        "partial_hessian": bool(partial_hessian),
        "n_atom_total": natom,
        "active_atom_count": len(active_list),
        "active_dof": ndof,
        "hessian_shape": [ndof, ndof],
        "hessian_maxabs": float(torch.max(torch.abs(h2))),
        "hessian_max_asym": float(hasym),
        "hessian_elapsed_s": float(elapsed),
        "mpi_enabled": bool(mpi_size > 1),
        "mpi_rank": int(mpi_rank),
        "mpi_size": int(mpi_size),
    }
    if used_threads is not None:
        out["num_threads"] = int(used_threads)
    if mode == "fd":
        out["hessian_delta_A"] = float(hessian_delta)
        out["fd_column_batch"] = int(fd_column_batch)
        out["force_evals"] = int(n_force_eval) if n_force_eval is not None else None
    elif mode == "analytical":
        out["hessian_delta_A"] = float(hessian_delta)
        if n_force_eval is not None:
            out["force_evals"] = int(n_force_eval)
        if analytical_meta is not None:
            for k, v in analytical_meta.items():
                out[k] = int(v)
    if partial_hessian:
        out["active_atoms"] = active_list

    if save_hessian is not None:
        save_path = Path(save_hessian)
        if mpi_size == 1 or mpi_rank == 0:
            torch.save(
                {
                    "hessian": h2,
                    "active_atoms": torch.tensor(active_list, dtype=torch.int64),
                    "partial_hessian": bool(partial_hessian),
                    "dtype": _precision_info(double)["torch_dtype"],
                },
                save_path,
            )
        out["hessian_file"] = str(save_path)

    out.update(_precision_info(double))
    return out


# ---------------------------------------------------------------------------
# OpenMM compatibility (lazy import)
# ---------------------------------------------------------------------------
def _openmm_positions_from_coords(coords: CoordsLike, natom: int):
    from openmm import app, unit

    if isinstance(coords, (str, Path)):
        coords_path = Path(coords)
        suffix = coords_path.suffix.lower()
        if suffix in {".rst7", ".inpcrd", ".crd", ".restrt"}:
            return app.AmberInpcrdFile(str(coords_path)).getPositions(asNumpy=True)
        if suffix == ".pdb":
            return app.PDBFile(str(coords_path)).getPositions(asNumpy=True)
        if suffix == ".xyz":
            xyz = load_coords(coords_path, natom=natom, device="cpu", dtype=torch.float64)
            return xyz.detach().cpu().numpy() * unit.angstrom
        raise ValueError(f"Unsupported coords file for OpenMM: {coords_path}")
    xyz = load_coords(coords, natom=natom, device="cpu", dtype=torch.float64)
    return xyz.detach().cpu().numpy() * unit.angstrom


def verify_openmm(
    prmtop: PathLike,
    coords: CoordsLike,
    device: Union[str, torch.device] = "cpu",
    openmm_platform: Optional[str] = None,
    double: bool = True,
) -> Dict[str, Any]:
    """Validate Torch energies/gradients against OpenMM."""
    try:
        from openmm import app, openmm, unit
    except Exception as e:
        raise RuntimeError(
            "OpenMM is not installed. Install it (pip install openmm) and retry."
        ) from e

    prmtop_path = Path(prmtop)
    _, xyz, ff = _load_runtime(
        prmtop=prmtop_path, coords=coords, device=device,
        double=double, requires_grad=True,
        nonbonded_cpu_fast=torch.device(device).type != "cpu",
    )
    e_torch = ff(xyz)
    e_total = e_torch["E_total"]
    e_total.backward()
    grad_torch = xyz.grad.detach().cpu().to(torch.float64)

    prmtop_omm = app.AmberPrmtopFile(str(prmtop_path))
    omm_system = prmtop_omm.createSystem(
        nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False,
    )
    for gi, force in enumerate(omm_system.getForces()):
        force.setForceGroup(gi)

    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    preferred = openmm_platform or "CPU"
    if str(preferred).strip().upper() != "CPU":
        raise ValueError(
            "Only OpenMM CPU platform is allowed in this workflow "
            f"(got openmm_platform={openmm_platform!r})."
        )
    names = [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())]
    if "CPU" not in names:
        raise RuntimeError(f"OpenMM CPU platform is not available. Platforms: {names}")
    selected_platform = "CPU"
    platform = openmm.Platform.getPlatformByName(selected_platform)
    context = openmm.Context(omm_system, integrator, platform)

    natom = int(omm_system.getNumParticles())
    context.setPositions(_openmm_positions_from_coords(coords, natom=natom))

    state = context.getState(getEnergy=True, getForces=True)
    e_openmm = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    f_openmm = state.getForces(asNumpy=True).value_in_unit(
        unit.kilocalories_per_mole / unit.angstrom
    )
    grad_openmm = -torch.tensor(f_openmm, dtype=torch.float64)

    omm_terms: Dict[str, float] = {}
    for gi, force in enumerate(omm_system.getForces()):
        st = context.getState(getEnergy=True, groups=1 << gi)
        omm_terms[type(force).__name__] = float(
            st.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        )

    e_t = float(e_total.detach().cpu())
    gdiff = grad_torch - grad_openmm

    out: Dict[str, Any] = {
        "E_total_torch_kcalmol": e_t,
        "E_total_openmm_kcalmol": float(e_openmm),
        "E_total_diff_kcalmol": float(e_t - float(e_openmm)),
        "grad_rms_kcalmolA": float(torch.sqrt(torch.mean(gdiff**2))),
        "grad_maxabs_kcalmolA": float(torch.max(torch.abs(gdiff))),
        "torch_terms_kcalmol": _energy_terms_to_float(e_torch),
        "openmm_force_terms_kcalmol": omm_terms,
        "notes": [
            "OpenMM NonbondedForce energy includes 1-4 exceptions; Torch prints E_coul14/E_lj14 separately.",
            "Agreement depends on matching coordinate units (Angstrom) and using NoCutoff + constraints=None + rigidWater=False.",
        ],
    }
    out["openmm_platform"] = platform.getName()

    out.update(_precision_info(double))
    return out
