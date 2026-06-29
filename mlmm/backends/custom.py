# mlmm/backends/custom.py

"""Load a user-supplied ASE Calculator for the ML region (``--calc-file``).

Lets a user drive the ML region of the ML/MM ONIOM coupling with an arbitrary
`ASE <https://wiki.fysik.dtu.dk/ase/>`_ Calculator (GFN-xTB via ``tblite`` /
``xtb-python``, DFTB+, ORCA, Psi4, …) instead of a built-in MLIP backend, with
no change to ``mlmm_toolkit``. The boundary is the standard ASE Calculator
interface (energy in eV, forces in eV/Å); :class:`_CustomBackend` (in
``mlmm_calc.py``) wraps the loaded calculator via the existing
``_ASEMLBackend`` adapter.

The user file must expose a factory callable (default name ``get_calculator``)
that returns an :class:`ase.calculators.calculator.Calculator`::

    # my_calc.py
    from ase.calculators.emt import EMT

    def get_calculator(charge=0, spin=1, device="auto", **kwargs):
        return EMT()

Then::

    mlmm sp -i complex.pdb --parm system.parm7 --calc-file my_calc.py
"""

from __future__ import annotations

import importlib.util
import inspect
import sys
from pathlib import Path
from typing import Any, Dict

_MODULE_COUNTER = 0


def _is_ase_calculator(obj: Any) -> bool:
    """Duck-typed check for an ASE-compatible Calculator."""
    try:
        from ase.calculators.calculator import Calculator as _ASECalc
    except ImportError:
        # ASE not importable here; fall back to the duck type (energy + forces).
        _ASECalc = None
    if _ASECalc is not None and isinstance(obj, _ASECalc):
        return True
    return hasattr(obj, "get_potential_energy") and hasattr(obj, "get_forces")


def _load_user_module(calc_file: str):
    """Import the user's Python file as an isolated module by path."""
    path = Path(calc_file).expanduser()
    if not path.is_file():
        raise ValueError(f"--calc-file not found: {path}")
    global _MODULE_COUNTER
    _MODULE_COUNTER += 1
    mod_name = f"mlmm_calc_file_{_MODULE_COUNTER}"
    spec = importlib.util.spec_from_file_location(mod_name, str(path))
    if spec is None or spec.loader is None:
        raise ValueError(f"Could not load --calc-file as a Python module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = module
    try:
        spec.loader.exec_module(module)
    except Exception as exc:  # surface the user's import-time error clearly
        sys.modules.pop(mod_name, None)
        raise ValueError(f"Error while importing --calc-file {path}: {exc}") from exc
    return module, path


def _call_factory(factory, *, charge: int, spin: int, device: str, extra: Dict[str, Any]):
    """Call the user factory, passing only the kwargs its signature accepts."""
    candidate: Dict[str, Any] = {
        "charge": charge,
        "spin": spin,
        "mult": spin,
        "multiplicity": spin,
        "device": device,
    }
    candidate.update(extra or {})
    try:
        sig = inspect.signature(factory)
    except (TypeError, ValueError):
        # Builtins / C callables expose no signature; pass the full candidate set.
        return factory(**candidate)
    params = sig.parameters.values()
    if any(p.kind == p.VAR_KEYWORD for p in params):
        return factory(**candidate)
    names = {p.name for p in params}
    return factory(**{k: v for k, v in candidate.items() if k in names})


def load_ase_calculator(
    calc_file: str,
    calc_factory: str = "get_calculator",
    *,
    charge: int = 0,
    spin: int = 1,
    device: str = "auto",
    **extra: Any,
):
    """Load and return an ASE Calculator from a user Python file.

    Parameters
    ----------
    calc_file
        Path to a Python file exposing ``calc_factory``.
    calc_factory
        Name of a callable returning an ASE Calculator, or of a module-level
        Calculator instance. Defaults to ``"get_calculator"``.
    charge, spin, device
        Forwarded to the factory when its signature accepts them.
    """
    module, path = _load_user_module(calc_file)
    if not hasattr(module, calc_factory):
        exported = [n for n in vars(module) if not n.startswith("_")]
        raise ValueError(
            f"--calc-file {path} has no attribute '{calc_factory}'. "
            f"Define `def {calc_factory}(charge=0, spin=1, **kwargs)` returning "
            f"an ASE Calculator (or use --calc-factory NAME). "
            f"Found top-level names: {', '.join(exported) or '(none)'}."
        )
    attr = getattr(module, calc_factory)

    # A Calculator instance bound directly to the factory name is accepted as-is.
    if _is_ase_calculator(attr):
        return attr
    if not callable(attr):
        raise ValueError(
            f"--calc-file {path}: '{calc_factory}' is neither a callable factory "
            f"nor an ASE Calculator instance (got {type(attr).__name__})."
        )
    try:
        calc = _call_factory(attr, charge=charge, spin=spin, device=device, extra=extra)
    except Exception as exc:
        raise ValueError(
            f"--calc-file {path}: calling {calc_factory}(...) failed: {exc}"
        ) from exc
    if not _is_ase_calculator(calc):
        raise ValueError(
            f"--calc-file {path}: {calc_factory}(...) returned {type(calc).__name__}, "
            f"not an ASE Calculator (needs get_potential_energy / get_forces)."
        )
    return calc
