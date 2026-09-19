"""CLI utilities for standardized exception handling and shared boilerplate."""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
import logging
import sys
import textwrap
import time
import traceback
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import click

from mlmm.core.utils import deep_update, load_yaml_dict

logger = logging.getLogger(__name__)


def canonicalize_calculator_section(config: Dict[str, Any]) -> Dict[str, Any]:
    """Merge the accepted ``mlmm`` section into ``calc`` safely."""

    normalized = deepcopy(config)
    canonical = normalized.get("calc")
    legacy = normalized.get("mlmm")
    if canonical is not None and not isinstance(canonical, Mapping):
        raise click.BadParameter("YAML section 'calc' must be a mapping.")
    if legacy is not None and not isinstance(legacy, Mapping):
        raise click.BadParameter("YAML section 'mlmm' must be a mapping.")
    if isinstance(legacy, Mapping):
        merged = deepcopy(dict(legacy))
        for key, value in dict(canonical or {}).items():
            if key in merged and merged[key] != value:
                raise click.BadParameter(
                    f"Conflicting YAML values for calc.{key} and mlmm.{key}."
                )
            merged[key] = deepcopy(value)
        normalized["calc"] = merged
        normalized.pop("mlmm", None)
    return normalized


def resolve_model_indices_setting(
    ctx: click.Context,
    yaml_cfg: Mapping[str, Any],
    model_indices: Optional[str],
    legacy_one_based: bool,
) -> tuple[Optional[str], bool]:
    """Resolve CLI/YAML ML-region indices with a fixed 1-based CLI contract."""

    explicit = make_is_param_explicit(ctx)
    if explicit("model_indices_str"):
        return model_indices, (
            bool(legacy_one_based)
            if explicit("model_indices_one_based")
            else True
        )
    normalized = canonicalize_calculator_section(dict(yaml_cfg))
    calc = normalized.get("calc")
    if not isinstance(calc, Mapping) or calc.get("model_indices") is None:
        return model_indices, True
    raw = calc.get("model_indices")
    if isinstance(raw, str):
        rendered = raw.strip()
    elif isinstance(raw, (list, tuple)):
        if not raw or any(isinstance(value, bool) for value in raw):
            raise click.BadParameter("calc.model_indices must contain atom indices.")
        try:
            rendered = ",".join(str(int(value)) for value in raw)
        except (TypeError, ValueError) as exc:
            raise click.BadParameter(
                "calc.model_indices must be a string or a sequence of integers."
            ) from exc
    else:
        raise click.BadParameter(
            "calc.model_indices must be a string or a sequence of integers."
        )
    if not rendered:
        raise click.BadParameter("calc.model_indices must not be empty.")
    base = calc.get("model_indices_base", 1)
    if isinstance(base, bool) or base not in (0, 1):
        raise click.BadParameter("calc.model_indices_base must be 0 or 1.")
    return rendered, bool(base == 1)


def make_is_param_explicit(ctx: "click.Context"):
    """Return a helper that checks whether a Click parameter was set explicitly."""
    from click.core import ParameterSource
    def _is_param_explicit(name: str) -> bool:
        try:
            source = ctx.get_parameter_source(name)
            return source not in (None, ParameterSource.DEFAULT)
        except Exception as exc:
            # Unknown param name: treat as not explicit (fall through to
            # YAML/defaults) rather than letting a Click default outrank a
            # --config entry. Log so a genuine typo is diagnosable.
            logger.debug("make_is_param_explicit: failed to query source for %r: %s", name, exc)
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


# Every top-level section any subcommand reads via apply_yaml_overrides. Deliberately the UNION
# over all commands: a per-command list would fire on a legitimate shared config file and become
# noise. `mlmm` is an accepted alias of `calc`.
_KNOWN_YAML_SECTIONS = frozenset({
    "geom", "calc", "mlmm", "opt", "lbfgs", "rfo", "rsirfo", "stopt", "gs", "dmf",
    "sp", "freq", "thermo", "microiter", "hessian_dimer", "irc", "dft", "bond",
    "bias", "search",
})


def _warn_unknown_yaml_sections(merged: Mapping[str, Any]) -> None:
    """Tell the user about YAML sections that will be silently ignored.

    A misspelled section name (``clac:`` for ``calc:``) is dropped without a word, so the run
    completes on stock defaults while the user believes the file took effect.
    """
    bad = sorted(k for k in merged if k not in _KNOWN_YAML_SECTIONS)
    if bad:
        click.echo(
            f"[config] WARNING: YAML section(s) {', '.join(bad)} are not recognized and were "
            f"ignored. Known sections: {', '.join(sorted(_KNOWN_YAML_SECTIONS))}.",
            err=True,
        )


def load_merged_yaml_cfg(
    config_yaml: Optional[Path],
    override_yaml: Optional[Path],
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """Load and merge YAML config and override files.

    Returns ``(merged, config_dict, override_dict)`` so that callers can
    use the individual layers for staged ``apply_yaml_overrides`` and the
    merged dict for ``show_config`` display without re-reading the files.
    """
    try:
        config_dict = canonicalize_calculator_section(load_yaml_dict(config_yaml))
        override_dict = canonicalize_calculator_section(load_yaml_dict(override_yaml))
    except ValueError as exc:
        raise click.BadParameter(str(exc), param_hint="--config/--override") from exc
    merged: Dict[str, Any] = {}
    deep_update(merged, config_dict)
    deep_update(merged, override_dict)
    merged = canonicalize_calculator_section(merged)
    _warn_unknown_yaml_sections(merged)
    return merged, config_dict, override_dict


def _write_error_json(
    out_dir: Path, command: str, exc: Exception, label: str,
    time_start: Optional[float] = None,
) -> None:
    """Write result.json + summary.json with a structured error envelope.

    Envelope fields beyond the existing ``status`` / ``error`` /
    ``error_type`` (kept for backward compatibility) include
    ``error_class_chain`` (MRO names so consumers can match
    ``OptimizationError`` → ``RuntimeError`` → ``Exception``),
    ``error_module`` (defining module, useful when the same class name
    is reused), and ``error_label`` (the high-level CLI stage label).
    """
    try:
        from mlmm.core.utils import write_result_json
        elapsed = time.perf_counter() - time_start if time_start else None
        klass = type(exc)
        write_result_json(
            out_dir,
            {
                "status": "error",
                "error": str(exc),
                "error_type": klass.__name__,
                "error_class_chain": [
                    c.__name__
                    for c in klass.__mro__
                    if c is not object
                ],
                "error_module": klass.__module__,
                "error_label": label,
            },
            command=command,
            elapsed_seconds=elapsed,
        )
    except OSError:
        pass  # Best-effort; don't mask the original error



def render_cli_exception(
    e,
    *,
    label: str,
    out_dir=None,
    command=None,
    time_start=None,
):
    """Shared terminal renderer for a subcommand's top-level exception.

    User-input errors (Click UsageError/BadParameter/ClickException, a
    malformed --config/override YAML) print a clean one-line ``Error:`` and exit with the conventional
    code (2). Everything else falls through to the full traceback (exit 1)
    so genuine internal bugs stay visible — no exception masking. Always
    calls sys.exit (never returns).
    """
    if out_dir is not None and command:
        # _write_error_json is best-effort and swallows its own errors
        # internally (see its definition), so it cannot mask the real error.
        _write_error_json(Path(out_dir).resolve(), command, e, label, time_start)
    if isinstance(e, click.ClickException):
        e.show()
        click.echo(
            f"Try 'mlmm {command} -h' for help." if command else "Try 'mlmm -h' for help.",
            err=True,
        )
        sys.exit(e.exit_code)
    try:
        import yaml as _yaml
        _is_yaml_err = isinstance(e, _yaml.YAMLError)
    except Exception:
        _is_yaml_err = False
    if _is_yaml_err:
        click.echo(
            f"Error: invalid YAML in a configuration/override file: {e}",
            err=True,
        )
        click.echo(
            f"Try 'mlmm {command} -h' for help." if command else "Try 'mlmm -h' for help.",
            err=True,
        )
        sys.exit(2)
    tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
    click.echo(
        f"Unhandled error during {label}:\n" + textwrap.indent(tb, "  "),
        err=True,
    )
    sys.exit(1)
