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


@pytest.mark.parametrize("configured", [False, True])
@pytest.mark.parametrize("render_error", [False, True])
def test_worker_avoids_shared_browser_configuration(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, configured: bool, render_error: bool,
) -> None:
    import os
    import plotly.io as pio

    shared_config = str(tmp_path / "shared-config")
    if configured:
        monkeypatch.setenv("XDG_CONFIG_HOME", shared_config)
    else:
        monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    caller_home = os.environ.get("HOME")
    seen_config: list[Path] = []

    class Figure:
        def write_image(self, path, **kwargs):
            config = os.environ.get("XDG_CONFIG_HOME")
            if config is None or config == shared_config:
                raise TimeoutError("browser configuration lock is held")
            seen_config.append(Path(config))
            assert Path(config).is_dir()
            if render_error:
                raise ValueError("render failed")
            Path(path).write_bytes(b"rendered")

    class Connection:
        def __init__(self):
            self.responses = []
            self.closed = False

        def send(self, value):
            self.responses.append(value)

        def close(self):
            self.closed = True

    monkeypatch.setattr(pio, "from_json", lambda value: Figure())
    monkeypatch.setattr(plotly_image.os, "setsid", lambda: None)
    connection = Connection()
    output = tmp_path / "figure.png"
    plotly_image._plotly_image_worker(
        "{}", str(tmp_path / "temporary.png"), str(output), "png", {}, connection,
    )
    expected = ("error", "ValueError: render failed") if render_error else ("ok", "")
    assert connection.responses == [expected]
    assert connection.closed
    assert bool(seen_config)
    assert all(not path.exists() for path in seen_config)
    if os.name == "posix":
        assert seen_config[0].parent == Path("/tmp")
    assert os.environ.get("XDG_CONFIG_HOME") == (shared_config if configured else None)
    assert os.environ.get("HOME") == caller_home
    if not render_error:
        assert output.read_bytes() == b"rendered"


def _navigate_commands() -> list[str]:
    import asyncio
    from unittest.mock import AsyncMock, Mock

    import choreographer
    from kaleido._kaleido_tab import _KaleidoTab

    async def navigate() -> list[str]:
        loop = asyncio.get_running_loop()
        javascript_ready = loop.create_future()
        page_ready = loop.create_future()
        tab = Mock(spec=choreographer.Tab)
        tab.subscribe_once.side_effect = {
            "Runtime.executionContextCreated": javascript_ready,
            "Page.loadEventFired": page_ready,
        }.__getitem__
        commands: list[str] = []

        async def send_command(command, params=None):
            commands.append(command)
            # Chrome does not replay a load event missed while Page was disabled.
            if command == "Page.navigate" and "Page.enable" in commands:
                page_ready.set_result({})
            elif command == "Runtime.enable":
                javascript_ready.set_result({"params": {"context": {"id": 2}}})
            return {"result": {}}

        tab.send_command = AsyncMock(side_effect=send_command)
        await asyncio.wait_for(_KaleidoTab(tab).navigate("file:///index.html"), timeout=1)
        return commands

    return asyncio.run(navigate())


def test_page_events_are_enabled_before_navigation(monkeypatch: pytest.MonkeyPatch) -> None:
    tab_module = pytest.importorskip("kaleido._kaleido_tab")
    monkeypatch.setattr(
        tab_module._KaleidoTab, "_apply_headers", tab_module._KaleidoTab._apply_headers
    )
    plotly_image._enable_page_events_before_navigation()

    commands = _navigate_commands()

    assert commands.index("Page.enable") < commands.index("Page.navigate")


def test_page_event_order_is_left_to_fixed_kaleido(monkeypatch: pytest.MonkeyPatch) -> None:
    import importlib.metadata

    tab_module = pytest.importorskip("kaleido._kaleido_tab")
    apply_headers = tab_module._KaleidoTab._apply_headers
    monkeypatch.setattr(tab_module._KaleidoTab, "_apply_headers", apply_headers)
    monkeypatch.setattr(importlib.metadata, "version", lambda name: "1.5.1")

    plotly_image._enable_page_events_before_navigation()

    assert tab_module._KaleidoTab._apply_headers is apply_headers
