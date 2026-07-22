"""Dependency-light adapter for MLMM's private console presentation tags."""

from __future__ import annotations

from typing import Any

import click


_TAG_AWARE_MARKER = "__mlmm_private_echo_tags__"


def emit(
    message: str = "",
    *,
    narrative: bool | None = None,
    detail: bool = False,
    force: bool = False,
    raw_path: bool = False,
    **kwargs: Any,
) -> None:
    """Echo one line while safely consuming private tags before bootstrap."""

    if narrative is None:
        narrative = not (detail or force or raw_path)

    if not getattr(click.echo, _TAG_AWARE_MARKER, False):
        click.echo(message, **kwargs)
        return
    click.echo(
        message,
        narrative=narrative,
        detail=detail,
        force=force,
        raw_path=raw_path,
        **kwargs,
    )
