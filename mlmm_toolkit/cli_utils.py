"""CLI utilities for standardized exception handling."""

from __future__ import annotations

import argparse
import sys
import textwrap
import traceback
from typing import Any, Callable, Optional, Type

import click


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


def run_cli(
    fn: Callable[[], None],
    *,
    label: str,
    zero_step_exc: Optional[Type[BaseException]] = None,
    zero_step_msg: Optional[str] = None,
    opt_exc: Optional[Type[BaseException]] = None,
    opt_msg: Optional[str] = None,
) -> None:
    """Standard CLI exception handling with consistent messaging."""
    try:
        fn()
    except KeyboardInterrupt:
        click.echo("Interrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
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
