"""Release policy for the retired experimental electronic-embedding path."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any


EMBEDCHARGE_UNAVAILABLE = (
    "Electronic embedding is unavailable in v0.3.3. The former experimental "
    "implementation double-counted ML/MM electrostatics, and its xTB correction "
    "used an inconsistent uncapped boundary model. Use mechanical embedding "
    "(`--no-embedcharge`) and rerun prior embedcharge results."
)


class EmbedChargeUnavailableError(ValueError):
    """Raised before a retired electronic-embedding calculation can start."""


def validate_retired_embedcharge(
    *,
    embedcharge: Any,
    cutoff_requested: bool = False,
) -> None:
    """Fail closed for every request that would activate or configure embedding."""
    if bool(embedcharge) or bool(cutoff_requested):
        raise EmbedChargeUnavailableError(EMBEDCHARGE_UNAVAILABLE)


def reject_retired_embedcharge_cli(
    calc_cfg: Mapping[str, Any],
    *,
    cutoff_requested: bool = False,
) -> None:
    """Apply the policy to an effective CLI calculator mapping."""
    try:
        validate_retired_embedcharge(
            embedcharge=calc_cfg.get("embedcharge", False),
            cutoff_requested=cutoff_requested,
        )
    except EmbedChargeUnavailableError as exc:
        import click

        raise click.UsageError(str(exc)) from exc
