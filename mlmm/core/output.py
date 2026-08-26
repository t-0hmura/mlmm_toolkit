"""Dependency-light adapter for MLMM's private console presentation tags."""

from __future__ import annotations

import re
from typing import Any

import click


_TAG_AWARE_MARKER = "__mlmm_private_echo_tags__"


def mlip_model_label(backend: Any, model: Any, task_name: Any = None) -> str:
    """Return a concise publication-facing MLIP model label."""
    backend_token = str(backend or "").strip().lower()
    model_token = str(model or "").strip()
    if not model_token:
        return "-"
    if backend_token == "uma":
        match = re.fullmatch(r"uma-([a-z]+)-([0-9]+(?:p[0-9]+)*)", model_token.lower())
        if match:
            model_label = "UMA-%s-%s" % (
                match.group(1).upper(), match.group(2).replace("p", ".")
            )
        else:
            model_label = model_token
        task_token = str(task_name or "omol").strip().lower()
        task_label = {
            "omol": "OMol", "omat": "OMat", "oc20": "OC20",
            "odac": "ODAC", "odac23": "ODAC23", "omc": "OMC",
            "omc25": "OMC25",
        }.get(task_token, str(task_name or "").strip())
        return f"{model_label} ({task_label})" if task_label else model_label
    if backend_token == "orb":
        label = re.sub(r"(?i)^orb[_-]", "ORB-", model_token).replace("_", "-")
        return label.replace("-omol", "-OMol").replace("-omat", "-OMat")
    if backend_token == "mace":
        off_match = re.fullmatch(r"(?i)off:(small|medium|large)", model_token)
        if off_match:
            return f"MACE-OFF23-{off_match.group(1).lower()}"
        label = re.sub(r"(?i)^mace-", "MACE-", model_token)
        for family in ("omol", "off23", "omat"):
            prefix = f"MACE-{family}"
            if label.lower().startswith(prefix.lower()):
                label = prefix.upper() + label[len(prefix):]
                break
        return label
    return model_token


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
