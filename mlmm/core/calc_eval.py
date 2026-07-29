"""Shared calculator evaluation helpers."""

from __future__ import annotations

import math
from typing import Any, Dict, Optional

import torch

from mlmm.backends.mlmm_calc import mlmm


def calc_energy(
    geom,
    calc_kwargs: Optional[Dict[str, Any]],
    calc: Optional[Any] = None,
) -> float:
    """Compute the ONIOM energy (Hartree) for ``geom``.

    Parameters
    ----------
    geom:
        Pysisyphus ``Geometry`` instance carrying ``atoms`` and Cartesian
        ``cart_coords``.
    calc_kwargs:
        Keyword arguments for constructing a fresh :class:`mlmm.backends.mlmm`
        when ``calc`` is ``None``. ``out_hess_torch`` is forced to ``False``
        so the temporary calculator does not allocate Hessian tensors for
        an energy-only call.
    calc:
        Pre-built calculator to reuse. When provided the helper does not
        construct or release a calculator, so the caller is responsible for
        teardown.
    """
    owns_calc = calc is None
    if owns_calc:
        kw = dict(calc_kwargs or {})
        kw["out_hess_torch"] = False
        calc = mlmm(**kw)
    try:
        result = calc.get_energy(geom.atoms, geom.cart_coords)
        if "energy" not in result:
            raise KeyError("Calculator result is missing 'energy'.")
        energy = float(result["energy"])
        if not math.isfinite(energy):
            raise ValueError("Calculator energy must be finite.")
    finally:
        if owns_calc:
            del calc
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
    return energy
