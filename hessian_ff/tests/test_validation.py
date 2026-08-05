from __future__ import annotations

from pathlib import Path

import pytest
import torch

from hessian_ff.forcefield import ForceFieldTorch
from hessian_ff.loaders import load_coords, load_system


_DATA_DIR = Path(__file__).resolve().parent / "data" / "small"
_PRMTOP = _DATA_DIR / "complex.parm7"
_COORDS = _DATA_DIR / "complex.rst7"

pytestmark = pytest.mark.skipif(
    (not _PRMTOP.exists()) or (not _COORDS.exists()),
    reason="benchmark/data/small test fixtures are not available in this environment",
)


def _load_ff_and_coords():
    system = load_system(_PRMTOP, device="cpu").to(dtype=torch.float64)
    ff = ForceFieldTorch(system, nonbonded_cpu_fast=True)
    return system, ff


def test_forcefield_explicit_dtype_validation_single() -> None:
    system, ff = _load_ff_and_coords()
    coords = load_coords(_COORDS, natom=system.natom, device="cpu", dtype=torch.float32)
    with pytest.raises(ValueError, match="coords dtype mismatch"):
        ff.energy_force(coords, force_calc_mode="Analytical")


def test_forcefield_explicit_dtype_validation_batch() -> None:
    system, ff = _load_ff_and_coords()
    coords = load_coords(_COORDS, natom=system.natom, device="cpu", dtype=torch.float32)
    coords_batch = coords.unsqueeze(0).repeat(2, 1, 1)
    with pytest.raises(ValueError, match="coords_batch dtype mismatch"):
        ff.energy_force_batch(coords_batch, force_calc_mode="Analytical")


def test_forcefield_to_keeps_system_and_buffers_synchronized() -> None:
    system, ff = _load_ff_and_coords()

    moved = ff.to(dtype=torch.float32)

    assert moved is ff
    assert ff.system.charge.dtype == torch.float32
    assert ff.nonbonded.charge.dtype == torch.float32
    assert ff.bond.k.dtype == torch.float32
    assert ff.system.charge.device == ff.nonbonded.charge.device


def test_default_batch_mode_avoids_vmap_native_kernels(monkeypatch) -> None:
    system, ff = _load_ff_and_coords()
    coords = load_coords(
        _COORDS, natom=system.natom, device="cpu", dtype=torch.float64
    )
    coords_batch = coords.unsqueeze(0).repeat(2, 1, 1)
    monkeypatch.setattr(
        torch,
        "vmap",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("default batch mode must not use vmap")
        ),
    )

    energies = ff.forward_batch(coords_batch)
    force_energies, forces = ff.energy_force_batch(coords_batch)

    assert energies["E_total"].shape == (2,)
    assert force_energies["E_total"].shape == (2,)
    assert forces.shape == coords_batch.shape


def test_unknown_batch_mode_is_rejected_with_a_deterministic_message() -> None:
    system, ff = _load_ff_and_coords()
    coords = load_coords(
        _COORDS, natom=system.natom, device="cpu", dtype=torch.float64
    )
    coords_batch = coords.unsqueeze(0)

    for call in (ff.forward_batch, ff.energy_force_batch):
        with pytest.raises(ValueError) as excinfo:
            call(coords_batch, batch_mode="loops")
        message = str(excinfo.value)
        # The accepted modes are listed in a fixed order, so the public
        # diagnostic is identical between runs.
        assert "('loop', 'vmap', 'vectorized')" in message
        assert "'loops'" in message
