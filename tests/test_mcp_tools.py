"""Argv-contract tests for MLMM MCP tool registration."""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from mlmm.mcp import _tools


class _FakeMCP:
    def __init__(self) -> None:
        self.tools = {}

    def tool(self):
        def register(function):
            self.tools[function.__name__] = function
            return function

        return register


class _FakeResult:
    def __init__(self, payload):
        self.payload = payload

    def to_dict(self):
        return dict(self.payload)


@pytest.fixture
def registry(monkeypatch: pytest.MonkeyPatch):
    calls = []

    def fake_run(argv, **kwargs):
        calls.append((list(argv), dict(kwargs)))
        return _FakeResult({"status": "ok", "argv": list(argv)})

    monkeypatch.setattr(_tools, "run_subcmd", fake_run)
    mcp = _FakeMCP()
    _tools.register_all(mcp)
    return mcp.tools, calls


@pytest.mark.parametrize(
    "extra_args",
    [
        ["-o", "other"],
        ["-oother"],
        ["--out-dir", "other"],
        ["--out-dir=other"],
        ["--out-json"],
        ["--no-out-json"],
    ],
)
def test_summary_tool_rejects_every_managed_output_spelling_before_spawn(
    registry, tmp_path: Path, extra_args: list[str]
) -> None:
    tools, calls = registry
    with pytest.raises(ValueError, match="MCP-managed output"):
        tools["optimize_geometry"](
            "input.pdb",
            "input.parm7",
            0,
            1,
            out_dir=str(tmp_path / "out"),
            extra_args=extra_args,
        )
    assert calls == []


@pytest.mark.parametrize(
    "extra_args",
    [
        ["-o", "other.pdb"],
        ["-oother.pdb"],
        ["--output", "other.pdb"],
        ["--output=other.pdb"],
        ["--out", "other.pdb"],
        ["--out=other.pdb"],
        ["--output-prefix", "other"],
        ["--output-file=other"],
    ],
)
def test_utility_tool_rejects_typed_output_overrides_before_spawn(
    registry, extra_args: list[str]
) -> None:
    tools, calls = registry
    with pytest.raises(ValueError, match="MCP-managed output"):
        tools["define_layer"](
            "input.pdb",
            "output.pdb",
            extra_args=extra_args,
        )
    assert calls == []


def test_tool_argv_preserves_boolean_toggle_syntax(registry, tmp_path: Path) -> None:
    tools, calls = registry
    tools["optimize_geometry"](
        "input.pdb",
        "input.parm7",
        -1,
        2,
        microiter=False,
        embedcharge=False,
        out_dir=str(tmp_path / "out"),
        extra_args=["--thresh", "gau"],
    )
    assert len(calls) == 1
    argv, kwargs = calls[0]
    assert argv[:2] == ["mlmm", "opt"]
    assert "--no-microiter" in argv
    assert "--no-embedcharge" in argv
    assert ["--embedcharge", "False"] != argv[-2:]
    assert argv[-2:] == ["--thresh", "gau"]
    assert kwargs["out_dir"] == tmp_path / "out"


def test_summary_only_tools_do_not_expose_leaf_pair_override(registry, tmp_path: Path) -> None:
    tools, calls = registry
    tools["run_full_pipeline"](
        "reactant.pdb",
        do_tsopt=True,
        do_dft=False,
        do_thermo=True,
        out_dir=str(tmp_path / "all"),
    )
    argv, kwargs = calls[0]
    assert argv[:2] == ["mlmm", "all"]
    assert "--tsopt" in argv
    assert "--no-dft" in argv
    assert "--thermo" in argv
    assert "true" not in argv and "false" not in argv
    assert "expected_primary_filename" not in kwargs


def test_search_paths_always_passes_two_ordered_endpoints(
    registry, tmp_path: Path,
) -> None:
    tools, calls = registry
    signature = inspect.signature(tools["search_paths"])
    assert signature.parameters["product_pdb"].default is inspect.Parameter.empty

    tools["search_paths"](
        "R.pdb",
        "full.parm7",
        -1,
        1,
        product_pdb="P.pdb",
        intermediate_pdbs=["IM1.pdb", "IM2.pdb"],
        out_dir=str(tmp_path / "path-search"),
    )

    argv, _ = calls[-1]
    input_at = argv.index("-i")
    assert argv[input_at + 1 : input_at + 5] == [
        "R.pdb", "IM1.pdb", "IM2.pdb", "P.pdb",
    ]
