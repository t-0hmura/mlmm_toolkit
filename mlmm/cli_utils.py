"""CLI utilities for standardized exception handling and shared boilerplate."""

from __future__ import annotations

import argparse
import gc
import os
import shutil
import sys
import textwrap
import time
import traceback
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Tuple, Type

import click

from .utils import deep_update, load_yaml_dict


# ---------------------------------------------------------------------------
# Click parameter source helper
# ---------------------------------------------------------------------------

def make_is_param_explicit(ctx: "click.Context"):
    """Return a helper that checks whether a Click parameter was set explicitly."""
    from click.core import ParameterSource
    def _is_param_explicit(name: str) -> bool:
        try:
            source = ctx.get_parameter_source(name)
            return source not in (None, ParameterSource.DEFAULT)
        except Exception:
            return False
    return _is_param_explicit


_TRUE_VALUES = {"true", "1", "yes", "y", "t"}
_FALSE_VALUES = {"false", "0", "no", "n", "f"}


def parse_bool(value: Any) -> bool:
    """Parse common boolean strings into bool; raise ValueError on invalid input."""
    if value is None:
        raise ValueError("Invalid boolean value: None. Use True/False.")
    text = str(value).strip().lower()
    if text in _TRUE_VALUES:
        return True
    if text in _FALSE_VALUES:
        return False
    raise ValueError(f"Invalid boolean value: {value!r}. Use True/False.")


def argparse_bool(value: str) -> bool:
    """argparse-compatible boolean parser using parse_bool()."""
    try:
        return parse_bool(value)
    except ValueError as e:
        raise argparse.ArgumentTypeError(str(e))


# ---------------------------------------------------------------------------
# YAML source resolution
# ---------------------------------------------------------------------------

def resolve_yaml_sources(
    config_yaml: Optional[Path],
    override_yaml: Optional[Path],
    args_yaml_legacy: Optional[Path],
) -> Tuple[Optional[Path], Optional[Path], bool]:
    """Resolve which YAML files to use, raising on conflicting options."""
    if override_yaml is not None and args_yaml_legacy is not None:
        raise click.BadParameter(
            "Use a single YAML source option."
        )
    if args_yaml_legacy is not None:
        return config_yaml, args_yaml_legacy, True
    return config_yaml, override_yaml, False


def load_merged_yaml_cfg(
    config_yaml: Optional[Path],
    override_yaml: Optional[Path],
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """Load and merge YAML config and override files.

    Returns ``(merged, config_dict, override_dict)`` so that callers can
    use the individual layers for staged ``apply_yaml_overrides`` and the
    merged dict for ``show_config`` display without re-reading the files.
    """
    config_dict = load_yaml_dict(config_yaml)
    override_dict = load_yaml_dict(override_yaml)
    merged: Dict[str, Any] = {}
    deep_update(merged, config_dict)
    deep_update(merged, override_dict)
    return merged, config_dict, override_dict


# ---------------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------------

def link_or_copy_file(src: Path, dst: Path) -> bool:
    """Create a symlink when possible; fall back to copy."""
    try:
        if dst.exists() or dst.is_symlink():
            if dst.is_dir():
                return False
            dst.unlink()
        rel = os.path.relpath(src, start=dst.parent)
        dst.symlink_to(rel)
        return True
    except OSError:
        try:
            shutil.copy2(src, dst)
            return True
        except OSError:
            return False


def _write_error_json(
    out_dir: Path, command: str, exc: Exception, label: str,
    time_start: Optional[float] = None,
) -> None:
    """Write a result.json with error status when a subcommand fails."""
    try:
        from .utils import write_result_json
        elapsed = time.perf_counter() - time_start if time_start else None
        write_result_json(
            out_dir,
            {"status": "error", "error": str(exc), "error_type": type(exc).__name__},
            command=command,
            elapsed_seconds=elapsed,
        )
    except Exception:
        pass  # Best-effort; don't mask the original error


# ---------------------------------------------------------------------------
# CLI exception wrapper
# ---------------------------------------------------------------------------

def run_cli(
    fn: Callable[[], None],
    *,
    label: str,
    zero_step_exc: Optional[Type[BaseException]] = None,
    zero_step_msg: Optional[str] = None,
    opt_exc: Optional[Type[BaseException]] = None,
    opt_msg: Optional[str] = None,
    out_dir: Optional[Path] = None,
    command: Optional[str] = None,
    time_start: Optional[float] = None,
) -> None:
    """Standard CLI exception handling with consistent messaging."""
    try:
        fn()
    except KeyboardInterrupt:
        click.echo("Interrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        if out_dir and command:
            _write_error_json(out_dir, command, e, label, time_start)
        if zero_step_exc is not None and isinstance(e, zero_step_exc):
            click.echo(
                zero_step_msg
                or "ERROR: Proposed step length dropped below the minimum allowed (ZeroStepLength).",
                err=True,
            )
            sys.exit(2)
        if opt_exc is not None and isinstance(e, opt_exc):
            msg = opt_msg or "ERROR: Optimization failed - {e}"
            click.echo(msg.format(e=e), err=True)
            sys.exit(3)
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo(
            f"Unhandled error during {label}:\n" + textwrap.indent(tb, "  "),
            err=True,
        )
        sys.exit(1)
    finally:
        # Release GPU memory (model + Hessian) after CLI command finishes
        # so that subsequent pipeline stages (e.g. tsopt → irc) don't OOM.
        # gc.collect() breaks cyclic refs inside torch.nn.Module.
        gc.collect()
        try:
            import torch
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        except ImportError:
            pass
