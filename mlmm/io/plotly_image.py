"""Bounded, process-isolated Plotly static-image export."""

from __future__ import annotations

import multiprocessing as mp
import os
import signal
import uuid
from pathlib import Path
from typing import Any, Mapping


PLOTLY_IMAGE_TIMEOUT_SECONDS = 600.0


class PlotlyImageTimeoutError(TimeoutError):
    """Raised when a Plotly/Kaleido image export exceeds its wall-time bound."""


def _plotly_image_worker(
    figure_json: str,
    temporary_path: str,
    output_path: str,
    image_format: str,
    kwargs: Mapping[str, Any],
    result_connection: Any,
) -> None:
    """Render in a new process group so Chrome descendants can be terminated."""

    if os.name == "posix":
        try:
            os.setsid()
        except OSError:
            pass

    try:
        from plotly import io as pio

        figure = pio.from_json(figure_json)
        figure.write_image(
            temporary_path,
            format=image_format,
            **dict(kwargs),
        )
        os.replace(temporary_path, output_path)
        result_connection.send(("ok", ""))
    except BaseException as exc:
        result_connection.send(
            ("error", f"{type(exc).__name__}: {exc}")
        )
    finally:
        result_connection.close()


def _terminate_process_tree(process: Any) -> None:
    """Terminate the renderer and its Chrome process group, then reap it."""

    if not process.is_alive():
        process.join()
        return

    terminated_group = False
    if os.name == "posix" and process.pid is not None:
        try:
            os.killpg(process.pid, signal.SIGTERM)
            terminated_group = True
        except (ProcessLookupError, PermissionError):
            pass
    if not terminated_group:
        process.terminate()
    process.join(5.0)

    if process.is_alive():
        killed_group = False
        if os.name == "posix" and process.pid is not None:
            try:
                os.killpg(process.pid, signal.SIGKILL)
                killed_group = True
            except (ProcessLookupError, PermissionError):
                pass
        if not killed_group:
            process.kill()
        process.join()


def write_plotly_image(
    figure: Any,
    output_path: str | Path,
    *,
    timeout_seconds: float = PLOTLY_IMAGE_TIMEOUT_SECONDS,
    **kwargs: Any,
) -> Path:
    """Write one Plotly image atomically, failing after ten minutes by default."""

    timeout = float(timeout_seconds)
    if timeout <= 0:
        raise ValueError("timeout_seconds must be positive")

    destination = Path(output_path).expanduser().resolve()
    image_format = str(kwargs.pop("format", destination.suffix.lstrip("."))).lower()
    if not image_format:
        raise ValueError("Plotly image output requires a filename extension or format")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.unlink(missing_ok=True)
    temporary = destination.with_name(
        f".{destination.name}.{uuid.uuid4().hex}.tmp"
    )

    context = mp.get_context("spawn")
    result_receiver, result_sender = context.Pipe(duplex=False)
    process = context.Process(
        target=_plotly_image_worker,
        args=(
            figure.to_json(),
            str(temporary),
            str(destination),
            image_format,
            dict(kwargs),
            result_sender,
        ),
        daemon=True,
    )
    process.start()
    result_sender.close()
    process.join(timeout)

    if process.is_alive():
        _terminate_process_tree(process)
        result_receiver.close()
        temporary.unlink(missing_ok=True)
        destination.unlink(missing_ok=True)
        raise PlotlyImageTimeoutError(
            f"Plotly image export timed out after {timeout:g} s: {destination}"
        )

    result = result_receiver.recv() if result_receiver.poll() else None
    result_receiver.close()
    if result is None or result[0] != "ok" or process.exitcode != 0:
        temporary.unlink(missing_ok=True)
        destination.unlink(missing_ok=True)
        detail = result[1] if result is not None else f"worker exit code {process.exitcode}"
        raise RuntimeError(f"Plotly image export failed: {detail}")

    return destination
