"""Dependency-light adapter for MLMM's private console presentation tags."""

from __future__ import annotations

from typing import Any

import click


_TAG_AWARE_MARKER = "__mlmm_private_echo_tags__"


def emit(
    message: Any = "",
    *,
    narrative: bool = True,
    detail: bool = False,
    **kwargs: Any,
) -> None:
    """Emit through the CLI gate when installed, or native Click otherwise.

    ``narrative``, ``detail``, ``force``, and ``raw_path`` are MLMM-private
    presentation metadata.  Programmatic/library use before CLI bootstrap
    consumes those tags locally so they never reach native ``click.echo``.
    """

    if getattr(click.echo, _TAG_AWARE_MARKER, False):
        click.echo(message, narrative=narrative, detail=detail, **kwargs)
        return
    kwargs.pop("force", None)
    kwargs.pop("raw_path", None)
    click.echo(message, **kwargs)
