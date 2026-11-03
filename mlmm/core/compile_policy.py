"""Whether UMA and ORB run their model through torch.compile (``--compile``).

The models run uncompiled unless ``--compile`` or ``calc.compile_model: true``
asks for it: compiling makes the first call slower and later calls only
slightly faster. The analytical Hessian is never compiled, because it needs a
double backward that the compiled model does not support. A compiled run checks
for the compiler before the model loads and stops with the install guidance of
``core.prerequisites`` when there is none.
"""

from __future__ import annotations

import os
import shutil
from typing import Any, Mapping, MutableMapping, Optional

import click

from mlmm.core.prerequisites import (
    KIND_C_COMPILER,
    KIND_CXX_COMPILER,
    MissingPrerequisiteError,
)

COMPILE_BACKENDS = ("uma", "orb")

_BACKEND_ERROR = "--compile applies to the uma and orb backends"
_ANALYTICAL_NOTE = "[compile] --compile is ignored for the analytical Hessian."
_analytical_noted = False
_NITRO_NOTE = (
    "[compile] --compile is ignored with --precision tf32 on UMA "
    "(UMA Nitro does not use torch.compile)."
)
_nitro_noted = False


def _on_cuda(device: Any) -> bool:
    text = "auto" if device is None else str(device).strip().lower()
    if text == "auto":
        import torch

        return bool(torch.cuda.is_available())
    return text.startswith("cuda")


def find_compiler(device: Any) -> Optional[str]:
    """Path of the compiler torch.compile needs on ``device``, or ``None``.

    On CUDA, Triton builds the kernels with a C compiler ($CC, then cc, gcc,
    clang); on a CPU, Inductor builds them with a C++ compiler ($CXX, then c++,
    g++, clang++).
    """
    if _on_cuda(device):
        variable, names = "CC", ("cc", "gcc", "clang")
    else:
        variable, names = "CXX", ("c++", "g++", "clang++")
    chosen = os.environ.get(variable, "").strip()
    for name in ((chosen,) if chosen else ()) + names:
        path = shutil.which(name)
        if path:
            return path
    return None


def resolve_compile(
    requested: Any,
    *,
    backend: str,
    device: Any,
    analytical_hessian: bool,
    uma_tf32_single_worker: bool = False,
) -> bool:
    """Whether one calculator compiles its model, given the resolved request."""
    global _analytical_noted, _nitro_noted
    if not requested:
        return False
    if str(backend).strip().lower() not in COMPILE_BACKENDS:
        raise click.UsageError(_BACKEND_ERROR)
    if analytical_hessian:
        if not _analytical_noted:
            _analytical_noted = True
            click.echo(_ANALYTICAL_NOTE, err=True)
        return False
    if uma_tf32_single_worker:
        if not _nitro_noted:
            _nitro_noted = True
            click.echo(_NITRO_NOTE, err=True)
        return False
    cuda = _on_cuda(device)
    if find_compiler("cuda" if cuda else "cpu") is None:
        raise MissingPrerequisiteError(KIND_C_COMPILER if cuda else KIND_CXX_COMPILER)
    return True


def apply_compile_option(calc_cfg: MutableMapping[str, Any], cli_value: Optional[bool]) -> None:
    """Resolve ``calc.compile_model`` from ``--compile/--no-compile`` and the config.

    Call once the backend is known. A CLI value and a config value that are
    both given and differ are a usage error; a single given value wins; neither
    leaves compilation off. ``--compile`` with a backend other than uma or orb
    is a usage error, while ``--no-compile`` is accepted there.
    """
    config_value = calc_cfg.get("compile_model")
    if config_value is not None and not isinstance(config_value, bool):
        raise click.BadParameter(
            f"calc.compile_model must be true or false, got {config_value!r}."
        )
    if cli_value is not None and config_value is not None and bool(cli_value) != config_value:
        flag = "--compile" if cli_value else "--no-compile"
        raise click.UsageError(
            f"{flag} and calc.compile_model: {str(config_value).lower()} conflict."
        )
    value = bool(cli_value) if cli_value is not None else bool(config_value)
    if value and str(calc_cfg.get("backend") or "uma").strip().lower() not in COMPILE_BACKENDS:
        raise click.UsageError(_BACKEND_ERROR)
    calc_cfg["compile_model"] = value


def runs_compiled(calc_cfg: Mapping[str, Any]) -> bool:
    """Whether the calculator built from a resolved ``calc_cfg`` compiles its model."""
    backend = str(calc_cfg.get("backend") or "uma").strip().lower()
    return (
        calc_cfg.get("compile_model") is True
        and backend in COMPILE_BACKENDS
        and str(calc_cfg.get("hessian_calc_mode") or "").strip().lower() != "analytical"
        and not (
            backend == "uma"
            and str(calc_cfg.get("uma_precision") or "").strip().lower() == "tf32"
            and int(calc_cfg.get("workers") or 1) <= 1
        )
    )


__all__ = [
    "COMPILE_BACKENDS",
    "apply_compile_option",
    "find_compiler",
    "resolve_compile",
    "runs_compiled",
]
