"""L4 Infra — MLIP backend adapters + ML/MM ONIOM coupling.

Modules:
- ``mlmm_calc`` — ``MLMMCore`` (ML/MM ONIOM coupling) + per-backend adapters
  (``_UMABackend`` / ``_OrbBackend`` / ``_MACEBackend`` / ``_AIMNet2Backend``),
  the ``_create_ml_backend`` factory, ``MLMMASECalculator`` (ASE), and ``mlmm``
  (pysisyphus Calculator).
- ``xtb_embedcharge_correction`` — ``--embedcharge`` correction (xTB point-charge
  embedding for MM→ML environmental effects).

User-facing API (factory pattern, per-backend kwargs, unified
``--precision fp32|fp64`` option, add-a-backend recipe) — see
``docs/backends.md``.
"""
from __future__ import annotations

import warnings
from typing import Any, Dict

# Backend-specific value domain for the unified ``--precision`` CLI flag.
# Routes the user-facing `fp32` / `fp64` choice to each backend's native
# kwarg via `_create_ml_backend`. ORB binding (`orb_precision`) follows
# orb_models matmul-precision strings; MACE uses numpy-dtype strings;
# AIMNet2 has no precision knob: fp32 is a no-op, fp64 is rejected.
_PRECISION_DISPATCH: Dict[str, Dict[str, tuple]] = {
    "fp32": {
        "uma":  ("uma_precision", "fp32"),
        "orb":  ("orb_precision", "float32-high"),
        "mace": ("mace_dtype", "float32"),
    },
    "fp64": {
        "uma":  ("uma_precision", "fp64"),
        "orb":  ("orb_precision", "float64"),
        "mace": ("mace_dtype", "float64"),
    },
}


def apply_precision_to_calc_cfg(calc_cfg: Dict[str, Any], precision: str) -> None:
    """Route the unified ``--precision`` CLI value into backend-specific kwargs.

    Mutates ``calc_cfg`` in place; for aimnet2 fp32 is a no-op and fp64
    is rejected (its model inputs are cast to float32 upstream).
    Raises ``ValueError`` for invalid values.
    """
    val = str(precision or "").lower()
    if val not in _PRECISION_DISPATCH:
        raise ValueError(
            f"--precision must be 'fp32' or 'fp64', got {precision!r}"
        )
    backend = str(calc_cfg.get("backend") or "uma").strip().lower()
    mapping = _PRECISION_DISPATCH[val]
    if backend not in mapping:
        # AIMNet2 (or any future backend with no precision knob) cannot honour
        # `--precision fp64`: upstream `aimnet` casts model inputs to float32
        # (`aimnet/calculators/calculator.py keys_in: torch.float`), so silently
        # accepting fp64 would lie about the actual numeric precision of the
        # run. Reject the combination loudly; fp32 stays a no-op so users can
        # swap `--backend` without changing scripts.
        if val == "fp64":
            raise ValueError(
                f"--precision fp64 is not supported by backend {backend!r}: "
                f"its model inputs are cast to float32 upstream, so the run "
                f"would not actually be fp64."
            )
        return
    kw_name, kw_val = mapping[backend]
    calc_cfg[kw_name] = kw_val
    if val == "fp64":
        # fp64 model precision implies a fp64 Hessian. Leaving ``H_double``
        # off while the forward pass is fp64 is internally inconsistent: the
        # optimizer / eigen linear algebra would discard exactly the precision
        # the user paid for. ``H_double`` is not a CLI flag (always True by
        # default), so the only way to reach the mismatch is a hand-edited
        # config. Warn that the deliberate inconsistency is being overridden,
        # then force fp64 on so the run stays self-consistent.
        if calc_cfg.get("H_double", True) is False:
            warnings.warn(
                "--precision fp64 forces a fp64 Hessian: overriding the "
                "H_double=False set in the config. A fp32 Hessian under fp64 "
                "model precision is internally inconsistent (the optimizer and "
                "eigen linear algebra would discard the fp64 precision). Use "
                "--precision fp32 if a fp32 Hessian is intended.",
                stacklevel=2,
            )
        calc_cfg["H_double"] = True
