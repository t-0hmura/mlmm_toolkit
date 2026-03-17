"""hessian_ff

Minimal AMBER (parm7/prmtop) force field evaluation in PyTorch.

Design goals:
  - ONIOM/MM style usage (no periodic boundary conditions, no Ewald/PME)
  - Differentiable energies (PyTorch autograd)
  - Match AMBER prmtop conventions: exclusions and 1-4 scaling

Primary entrypoints:
  - :func:`hessian_ff.load_system`
  - :class:`hessian_ff.forcefield.ForceFieldTorch`
  - :func:`hessian_ff.workflows.torch_energy`
  - :func:`hessian_ff.workflows.torch_force`
  - :func:`hessian_ff.workflows.torch_hessian`
"""

from .system import AmberSystem
from .forcefield import ForceFieldTorch
from .loaders import load_system, load_coords
from .native import build_native_extensions
from .workflows import (
    clear_runtime_cache,
    system_summary,
    torch_energy,
    torch_energy_batch,
    torch_force,
    torch_force_batch,
    torch_hessian,
    verify_openmm,
)

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "AmberSystem",
    "ForceFieldTorch",
    "build_native_extensions",
    "load_system",
    "load_coords",
    "clear_runtime_cache",
    "system_summary",
    "torch_energy",
    "torch_energy_batch",
    "torch_force",
    "torch_force_batch",
    "torch_hessian",
    "verify_openmm",
]
