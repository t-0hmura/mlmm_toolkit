"""Shared fixtures for hessian_ff tests."""

from __future__ import annotations

from pathlib import Path

import pytest
import torch

from hessian_ff.forcefield import ForceFieldTorch
from hessian_ff.loaders import load_coords, load_system

DATA_DIR = Path(__file__).resolve().parent / "data" / "small"
PRMTOP = DATA_DIR / "complex.parm7"
COORDS = DATA_DIR / "complex.rst7"
PDB = DATA_DIR / "complex.pdb"


def _skip_no_data():
    if not (PRMTOP.exists() and COORDS.exists()):
        pytest.skip("hessian_ff test data (data/small/) not available")


def _skip_no_openmm():
    try:
        from openmm import openmm as omm  # noqa: F401
    except ImportError:
        pytest.skip("OpenMM is not installed")


@pytest.fixture(scope="session")
def small_system():
    """Load AmberSystem from small test data (float64, CPU)."""
    _skip_no_data()
    return load_system(PRMTOP, device="cpu").to(dtype=torch.float64)


@pytest.fixture(scope="session")
def small_ff(small_system):
    """ForceFieldTorch for the small test system."""
    return ForceFieldTorch(small_system, nonbonded_cpu_fast=True)


@pytest.fixture(scope="session")
def small_ff_autograd(small_system):
    """ForceFieldTorch without nonbonded_cpu_fast (needed for autograd Hessian)."""
    return ForceFieldTorch(small_system, nonbonded_cpu_fast=False)


@pytest.fixture(scope="session")
def small_coords(small_system):
    """Coordinates tensor for the small test system (float64, CPU)."""
    return load_coords(COORDS, natom=small_system.natom, device="cpu", dtype=torch.float64)


@pytest.fixture(scope="session")
def openmm_context_and_pos(small_system):
    """Create an OpenMM context and base positions for the small system."""
    _skip_no_openmm()
    from openmm import app, openmm as omm, unit

    prmtop = app.AmberPrmtopFile(str(PRMTOP))
    system = prmtop.createSystem(
        nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False,
    )
    integ = omm.VerletIntegrator(1.0 * unit.femtoseconds)
    platform = omm.Platform.getPlatformByName("CPU")
    context = omm.Context(system, integ, platform, {"Threads": "4"})

    coords = load_coords(
        COORDS, natom=small_system.natom, device="cpu", dtype=torch.float64,
    )
    pos_a = coords.detach().cpu().numpy()
    context.setPositions(pos_a * unit.angstrom)
    return context, pos_a
