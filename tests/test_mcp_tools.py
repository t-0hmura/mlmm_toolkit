"""Argv-contract tests for MLMM MCP tool registration."""

from __future__ import annotations

import asyncio
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
        asyncio.run(
            tools["optimize_geometry"](
                "input.pdb",
                "input.parm7",
                0,
                1,
                out_dir=str(tmp_path / "out"),
                extra_args=extra_args,
            )
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
        ["--out-prefix", "other"],
        ["--out-prefix=other"],
        ["--output-file=other"],
    ],
)
def test_utility_tool_rejects_typed_output_overrides_before_spawn(
    registry, extra_args: list[str]
) -> None:
    tools, calls = registry
    with pytest.raises(ValueError, match="MCP-managed output"):
        asyncio.run(
            tools["define_layer"](
                "input.pdb",
                "output.pdb",
                extra_args=extra_args,
            )
        )
    assert calls == []


def test_tool_argv_preserves_boolean_toggle_syntax(registry, tmp_path: Path) -> None:
    tools, calls = registry
    asyncio.run(
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
    asyncio.run(
        tools["run_full_pipeline"](
            "reactant.pdb",
            do_tsopt=True,
            do_dft=False,
            do_thermo=True,
            out_dir=str(tmp_path / "all"),
        )
    )
    argv, kwargs = calls[0]
    assert argv[:2] == ["mlmm", "all"]
    assert "--tsopt" in argv
    assert "--no-dft" in argv
    assert "--thermo" in argv
    assert "true" not in argv and "false" not in argv
    assert "expected_primary_filename" not in kwargs


def test_single_point_tool_forwards_print_every(registry, tmp_path: Path) -> None:
    tools, calls = registry
    signature = inspect.signature(tools["run_single_point_oniom"])
    assert "print_every" in signature.parameters
    assert "detect_layer" not in signature.parameters

    asyncio.run(
        tools["run_single_point_oniom"](
            "input.pdb",
            "input.parm7",
            charge=0,
            multiplicity=1,
            ref_pdb="topology.pdb",
            print_every=3,
            out_dir=str(tmp_path / "sp"),
        )
    )

    argv, _kwargs = calls[-1]
    option_start = argv.index("--print-every")
    assert argv[option_start : option_start + 2] == ["--print-every", "3"]
    ref_start = argv.index("--ref-pdb")
    assert argv[ref_start : ref_start + 2] == ["--ref-pdb", "topology.pdb"]


def test_bond_change_tool_owns_json_stdout_contract(registry) -> None:
    tools, calls = registry
    asyncio.run(tools["detect_bond_changes"]("R.pdb", "P.pdb"))
    argv, kwargs = calls[-1]
    assert argv == ["mlmm", "bond-summary", "-i", "R.pdb", "P.pdb", "--json"]
    assert kwargs["out_dir"] is None
    assert kwargs["parse_stdout_json"] is True

    with pytest.raises(ValueError, match="MCP-managed output"):
        asyncio.run(
            tools["detect_bond_changes"](
                "R.pdb", "P.pdb", extra_args=["--no-json"]
            )
        )


def test_search_paths_always_passes_two_ordered_endpoints(
    registry, tmp_path: Path,
) -> None:
    tools, calls = registry
    signature = inspect.signature(tools["search_paths"])
    assert signature.parameters["product_pdb"].default is inspect.Parameter.empty

    asyncio.run(
        tools["search_paths"](
            "R.pdb",
            "full.parm7",
            -1,
            1,
            product_pdb="P.pdb",
            intermediate_pdbs=["IM1.pdb", "IM2.pdb"],
            out_dir=str(tmp_path / "path-search"),
        )
    )

    argv, _ = calls[-1]
    input_at = argv.index("-i")
    assert argv[input_at + 1 : input_at + 5] == [
        "R.pdb", "IM1.pdb", "IM2.pdb", "P.pdb",
    ]


@pytest.mark.parametrize(
    ("tool_name", "args", "kwargs"),
    [
        ("compute_frequencies", ("in.pdb", "in.parm7", 0, 1), {}),
        ("scan_1d", ("in.pdb", "in.parm7", 0, 1, "1,2,1.5"), {}),
        ("scan_2d", ("in.pdb", "in.parm7", 0, 1, "1,2,1.5;2,3,2.0"), {}),
        ("scan_3d", ("in.pdb", "in.parm7", 0, 1, "1,2,1.5;2,3,2.0;3,4,2.5"), {}),
        ("optimize_path", ("R.pdb", "P.pdb", "in.parm7", 0, 1), {}),
        (
            "search_paths",
            ("R.pdb", "in.parm7", 0, 1),
            {"product_pdb": "P.pdb"},
        ),
        ("run_full_pipeline", ("R.pdb",), {}),
        ("run_single_point_dft", ("in.pdb", "in.parm7", 0, 1), {}),
    ],
)
def test_stage_tools_forward_typed_mm_controls(
    registry,
    tmp_path: Path,
    tool_name: str,
    args: tuple,
    kwargs: dict,
) -> None:
    tools, calls = registry
    asyncio.run(
        tools[tool_name](
            *args,
            **kwargs,
            link_atom_method="scaled",
            mm_backend="openmm",
            use_cmap=False,
            out_dir=str(tmp_path / tool_name),
        )
    )
    argv, _ = calls[-1]
    assert ["--link-atom-method", "scaled"] == argv[
        argv.index("--link-atom-method") : argv.index("--link-atom-method") + 2
    ]
    assert ["--mm-backend", "openmm"] == argv[
        argv.index("--mm-backend") : argv.index("--mm-backend") + 2
    ]
    assert "--no-cmap" in argv
