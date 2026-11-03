"""The dependency log filters drop only their two messages, below -v 3."""

from __future__ import annotations

import logging
import sys

import pytest

from mlmm.core import log_filters, utils
from mlmm.io import plotly_image

_FAIRCHEM_MSG = (
    "Running (unknown commit) with inference settings: "
    "InferenceSettings(tf32=False, activation_checkpointing=False, merge_mole=True)"
)
_KILL_MSG = "Resorting to unclean kill browser."


def _record(name: str, msg: str, level: int = logging.WARNING) -> logging.LogRecord:
    return logging.LogRecord(name, level, __file__, 1, msg, None, None)


def _remove_filters(root: logging.Logger, chrome: logging.Logger) -> None:
    root.filters[:] = [f for f in root.filters if f is not log_filters._FAIRCHEM_FILTER]
    chrome.filters[:] = [
        f for f in chrome.filters if f is not log_filters._CHOREOGRAPHER_FILTER
    ]


@pytest.fixture
def installed():
    import logistro

    root = logging.getLogger()
    chrome = logistro.getLogger(log_filters._CHOREOGRAPHER_LOGGER)
    had_root = log_filters._FAIRCHEM_FILTER in root.filters
    had_chrome = log_filters._CHOREOGRAPHER_FILTER in chrome.filters
    log_filters.install_fairchem_log_filter()
    log_filters.install_choreographer_log_filter()
    log_filters.install_fairchem_log_filter()
    log_filters.install_choreographer_log_filter()
    yield root, chrome
    if not had_root:
        root.removeFilter(log_filters._FAIRCHEM_FILTER)
    if not had_chrome:
        chrome.removeFilter(log_filters._CHOREOGRAPHER_FILTER)


@pytest.mark.parametrize("level", [0, 1, 2])
def test_filters_drop_exactly_the_two_messages_by_default(monkeypatch, level):
    monkeypatch.setattr(utils, "_VERBOSE_LEVEL", level)
    fairchem = log_filters.FairchemInferenceSettingsFilter()
    chrome = log_filters.ChoreographerUncleanKillFilter()

    assert not fairchem.filter(_record("root", _FAIRCHEM_MSG))
    assert not fairchem.filter(
        _record("root", "Running (3f2a9c1) with inference settings: InferenceSettings()")
    )
    assert not chrome.filter(_record(log_filters._CHOREOGRAPHER_LOGGER, _KILL_MSG))
    # Anything else passes.
    assert fairchem.filter(_record("root", "Running the model"))
    assert fairchem.filter(_record("fairchem.core", _FAIRCHEM_MSG))
    assert fairchem.filter(_record("root", "x " + _FAIRCHEM_MSG))
    assert chrome.filter(_record(log_filters._CHOREOGRAPHER_LOGGER, "Browser closed."))
    assert chrome.filter(_record(log_filters._CHOREOGRAPHER_LOGGER, _KILL_MSG + " now"))


def test_filters_keep_both_messages_at_level_three(monkeypatch):
    monkeypatch.setattr(utils, "_VERBOSE_LEVEL", 3)
    assert log_filters.FairchemInferenceSettingsFilter().filter(
        _record("root", _FAIRCHEM_MSG)
    )
    assert log_filters.ChoreographerUncleanKillFilter().filter(
        _record(log_filters._CHOREOGRAPHER_LOGGER, _KILL_MSG)
    )


def test_level_reads_zero_without_mlmm_core_utils(monkeypatch):
    monkeypatch.setattr(utils, "_VERBOSE_LEVEL", 3)
    assert log_filters.cli_verbose_level() == 3
    monkeypatch.delitem(sys.modules, "mlmm.core.utils")
    assert log_filters.cli_verbose_level() == 0
    assert not log_filters.ChoreographerUncleanKillFilter().filter(
        _record(log_filters._CHOREOGRAPHER_LOGGER, _KILL_MSG)
    )


def test_install_is_idempotent_and_keeps_logistro_logger_class(installed):
    root, chrome = installed
    assert root.filters.count(log_filters._FAIRCHEM_FILTER) == 1
    assert chrome.filters.count(log_filters._CHOREOGRAPHER_FILTER) == 1
    # choreographer calls logistro's extra methods on this logger.
    assert hasattr(chrome, "debug2")


@pytest.mark.parametrize("level", [2, 3])
def test_installed_filters_on_the_emitting_loggers(installed, monkeypatch, caplog, level):
    _, chrome = installed
    monkeypatch.setattr(utils, "_VERBOSE_LEVEL", level)

    caplog.set_level(logging.WARNING)
    # FAIR-Chem calls logging.warning() on the root logger.
    logging.warning(_FAIRCHEM_MSG)
    logging.warning("unrelated root warning")
    chrome.warning(_KILL_MSG)
    chrome.warning("unrelated choreographer warning")

    messages = [r.getMessage() for r in caplog.records]
    assert "unrelated root warning" in messages
    assert "unrelated choreographer warning" in messages
    assert (_FAIRCHEM_MSG in messages) is (level >= 3)
    assert (_KILL_MSG in messages) is (level >= 3)


def test_cli_start_installs_the_filters(installed):
    from click.testing import CliRunner

    from mlmm.cli import cli

    root, chrome = installed
    _remove_filters(root, chrome)

    result = CliRunner().invoke(cli, ["opt", "--help"])

    assert result.exit_code == 0, result.output
    assert log_filters._FAIRCHEM_FILTER in root.filters
    assert log_filters._CHOREOGRAPHER_FILTER not in chrome.filters


@pytest.mark.parametrize("hide", [False, True])
def test_png_worker_installs_the_filters_when_told(installed, monkeypatch, tmp_path, hide):
    import plotly.io as pio

    root, chrome = installed
    _remove_filters(root, chrome)

    class Figure:
        def write_image(self, path, **kwargs):
            with open(path, "wb") as handle:
                handle.write(b"rendered")

    class Connection:
        def send(self, value):
            self.value = value

        def close(self):
            pass

    monkeypatch.setattr(pio, "from_json", lambda value: Figure())
    monkeypatch.setattr(plotly_image.os, "setsid", lambda: None)
    connection = Connection()
    plotly_image._plotly_image_worker(
        "{}", str(tmp_path / "temporary.png"), str(tmp_path / "figure.png"), "png", {},
        connection, hide,
    )

    assert connection.value == ("ok", "")
    assert log_filters._FAIRCHEM_FILTER not in root.filters
    assert (log_filters._CHOREOGRAPHER_FILTER in chrome.filters) is hide


@pytest.mark.parametrize("level", [2, 3])
def test_png_export_tells_the_worker_the_parent_level(monkeypatch, tmp_path, level):
    monkeypatch.setattr(utils, "_VERBOSE_LEVEL", level)
    started: list[tuple] = []

    class Connection:
        def close(self):
            pass

        def poll(self):
            return True

        def recv(self):
            return ("ok", "")

    class FinishedProcess:
        exitcode = 0

        def __init__(self, target, args, daemon):
            started.append(args)

        def start(self):
            pass

        def join(self, timeout=None):
            pass

        def is_alive(self):
            return False

    class Context:
        def Pipe(self, duplex=False):
            return Connection(), Connection()

        def Process(self, **kwargs):
            return FinishedProcess(**kwargs)

    class Figure:
        def to_json(self):
            return "{}"

    monkeypatch.setattr(plotly_image.mp, "get_context", lambda method: Context())
    plotly_image.write_plotly_image(Figure(), tmp_path / "plot.png")

    assert started[0][-1] is (level < 3)
