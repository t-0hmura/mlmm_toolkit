"""GSM progress tables retain their segment section headings."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest
from click.testing import CliRunner

from mlmm.cli import cli as root_cli


class _OptimizerStopped(Exception):
    pass


def _record_gsm_order(monkeypatch, module) -> list:
    """Collect the section headings and the start of the string optimizer, in order."""
    events: list = []

    class FakeString:
        def __init__(self, **kwargs):
            pass

    class FakeOptimizer:
        def __init__(self, geometry, **kwargs):
            pass

        def run(self):
            events.append("optimizer.run")
            raise _OptimizerStopped

    monkeypatch.setattr(module, "GrowingString", FakeString)
    monkeypatch.setattr(module, "StringOptimizer", FakeOptimizer)
    monkeypatch.setattr(module, "emit", lambda text="", **tags: events.append(text))
    return events


def test_recursive_gsm_emits_its_segment_section_before_optimization(tmp_path, monkeypatch) -> None:
    from mlmm.workflows import path_search

    events = _record_gsm_order(monkeypatch, path_search)
    endpoint = SimpleNamespace(set_calculator=lambda calculator: None)

    with pytest.raises(_OptimizerStopped):
        path_search._run_gsm_between(
            endpoint, endpoint, object(), {}, {}, tmp_path, "seg_00", None,
        )

    assert events == ["\n====== [seg_00] GSM ======\n", "optimizer.run"]


def test_standalone_gsm_emits_its_section_before_optimization(tmp_path, monkeypatch) -> None:
    from mlmm.workflows import path_opt

    toy = Path(__file__).resolve().parents[1] / "examples" / "toy_system"
    if not (toy / "p_toy.parm7").exists():
        pytest.skip("toy_system example inputs not present")
    events = _record_gsm_order(monkeypatch, path_opt)
    monkeypatch.setattr(path_opt, "mlmm", lambda **kwargs: SimpleNamespace(core=None))

    CliRunner().invoke(root_cli, [
        "path-opt", "-i", str(toy / "r_toy.pdb"), str(toy / "p_toy.pdb"),
        "--parm", str(toy / "p_toy.parm7"), "-q", "-1", "-m", "1", "--no-preopt",
        "--out-dir", str(tmp_path / "out"),
    ])

    heading = "\n====== Growing String optimization ======\n"
    assert heading in events and "optimizer.run" in events
    assert events.index(heading) < events.index("optimizer.run")
