from __future__ import annotations

from pathlib import Path
import signal

import pytest

from mlmm.io import plotly_image


class _Figure:
    def to_json(self) -> str:
        return "{}"


class _Connection:
    def __init__(self) -> None:
        self.closed = False

    def close(self) -> None:
        self.closed = True

    def poll(self) -> bool:
        return False


class _HungProcess:
    pid = 12345
    exitcode = None

    def __init__(self, **kwargs) -> None:
        self.kwargs = kwargs
        self.alive = True
        self.join_timeouts: list[float] = []

    def start(self) -> None:
        pass

    def join(self, timeout=None) -> None:
        self.join_timeouts.append(timeout)

    def is_alive(self) -> bool:
        return self.alive


class _Context:
    def __init__(self) -> None:
        self.receiver = _Connection()
        self.sender = _Connection()
        self.process: _HungProcess | None = None

    def Pipe(self, duplex=False):
        assert duplex is False
        return self.receiver, self.sender

    def Process(self, **kwargs):
        self.process = _HungProcess(**kwargs)
        return self.process


def test_plotly_image_timeout_is_ten_minutes() -> None:
    assert plotly_image.PLOTLY_IMAGE_TIMEOUT_SECONDS == 600.0


def test_plotly_image_timeout_terminates_worker_and_removes_stale_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    context = _Context()
    monkeypatch.setattr(plotly_image.mp, "get_context", lambda method: context)

    terminated: list[_HungProcess] = []

    def terminate(process: _HungProcess) -> None:
        terminated.append(process)
        process.alive = False

    monkeypatch.setattr(plotly_image, "_terminate_process_tree", terminate)
    output = tmp_path / "plot.png"
    output.write_bytes(b"stale")

    with pytest.raises(plotly_image.PlotlyImageTimeoutError, match="0.01 s"):
        plotly_image.write_plotly_image(
            _Figure(), output, timeout_seconds=0.01, scale=2
        )

    assert context.process is not None
    assert context.process.join_timeouts == [0.01]
    assert terminated == [context.process]
    assert context.sender.closed is True
    assert context.receiver.closed is True
    assert not output.exists()


def test_plotly_image_rejects_nonpositive_timeout(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="must be positive"):
        plotly_image.write_plotly_image(
            _Figure(), tmp_path / "plot.png", timeout_seconds=0
        )


def test_terminate_process_tree_signals_renderer_group(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    process = _HungProcess()
    signals: list[tuple[int, signal.Signals]] = []

    def killpg(pid: int, signum: signal.Signals) -> None:
        signals.append((pid, signum))
        process.alive = False

    monkeypatch.setattr(plotly_image.os, "killpg", killpg)
    plotly_image._terminate_process_tree(process)

    assert signals == [(process.pid, signal.SIGTERM)]
    assert process.join_timeouts == [5.0]
