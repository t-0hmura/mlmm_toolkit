"""Public result fields and command completion codes."""

from functools import wraps

import click

EXECUTION_STATUS_VALUES = ("completed", "failed")
SCIENTIFIC_STATUS_VALUES = ("success", "partial", "failed")

_COMPLETION_KEY = "mlmm.completion"
_EXECUTION_FAILURE_KEY = "mlmm.execution_failure"


def result_fields(data, *, command=None):
    """Keep execution and scientific outcomes distinct from engine diagnostics."""
    result = dict(data)
    terminal = result.pop("status", None)
    result.pop("status_reasons", None)
    if command in {"opt", "tsopt", "path-opt"} and terminal in {
        "converged", "not_converged", "stalled", "unknown", "completed",
    }:
        result.setdefault("optimization_status", terminal)
    execution = result.get("execution_status")
    if execution is None:
        execution = "failed" if terminal == "error" or result.get("error_type") else "completed"
    if any(point.get("executed") is False for point in result.get("point_outcomes", [])):
        execution = "failed"
    scientific = result.get("scientific_status")
    if scientific is None:
        if terminal in {"success", "partial", "failed"}:
            scientific = terminal
        elif terminal in {"ok", "completed", "converged"}:
            scientific = "success"
        else:
            scientific = "failed"
    if result.get("hessian_status") == "failed":
        execution = "failed"
        if scientific == "success":
            scientific = "partial"
    if execution not in EXECUTION_STATUS_VALUES:
        raise ValueError(f"Invalid execution_status: {execution!r}")
    if scientific not in SCIENTIFIC_STATUS_VALUES:
        raise ValueError(f"Invalid scientific_status: {scientific!r}")
    result["execution_status"] = execution
    result["scientific_status"] = scientific
    return result


def record_completion(data, *, command=None):
    """Record a verdict even when JSON output was not requested."""
    result = result_fields(data, command=command)
    context = click.get_current_context(silent=True)
    if context is not None:
        if context.meta.get(_EXECUTION_FAILURE_KEY):
            result["execution_status"] = "failed"
        context.meta[_COMPLETION_KEY] = result
    return result


def record_child_failure(verdict=None):
    """Carry a caught child exception into the parent execution outcome."""
    context = click.get_current_context(silent=True)
    if context is not None and (verdict is None or verdict.get("execution_status") == "failed"):
        context.meta[_EXECUTION_FAILURE_KEY] = True


def completion_code(data):
    return int(
        data["execution_status"] == "failed"
        or data["scientific_status"] == "failed"
    )


def completion_guard(callback):
    """Apply the result verdict after the command has saved outputs and cleaned up."""
    @wraps(callback)
    def run(*args, **kwargs):
        context = click.get_current_context(silent=True)
        if context is None:
            return callback(*args, **kwargs)
        previous = context.meta.pop(_COMPLETION_KEY, None)
        previous_failure = context.meta.pop(_EXECUTION_FAILURE_KEY, None)
        try:
            try:
                value = callback(*args, **kwargs)
            except SystemExit as exc:
                verdict = context.meta.get(_COMPLETION_KEY)
                if verdict is not None:
                    exc.completion_result = verdict
                raise
            except KeyboardInterrupt as exc:
                raise SystemExit(130) from exc
            verdict = context.meta.get(_COMPLETION_KEY)
            if verdict is not None and completion_code(verdict):
                stopped = SystemExit(1)
                stopped.completion_result = verdict
                raise stopped
            return value
        finally:
            context.meta.pop(_COMPLETION_KEY, None)
            context.meta.pop(_EXECUTION_FAILURE_KEY, None)
            if previous_failure is not None:
                context.meta[_EXECUTION_FAILURE_KEY] = previous_failure
            if previous is not None:
                context.meta[_COMPLETION_KEY] = previous
    return run
