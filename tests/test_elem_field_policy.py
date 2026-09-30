"""Element-column policy: the PDB reader and ``add-elem-info`` keep valid element columns."""

from __future__ import annotations

import asyncio
from pathlib import Path

from click.testing import CliRunner
import pytest

from mlmm.cli import cli as root_cli
from mlmm.domain.add_elem_info import assign_elements
from mlmm.mcp import _tools


def _atom(serial: int, name: str, resname: str, element: str, record: str = "HETATM") -> str:
    return (
        f"{record:<6}{serial:5d} {name:4s} {resname:>3s} A{serial:4d}    "
        f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}"
        f"{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
    )


def _write(path: Path, records: list[tuple[str, str, str]]) -> Path:
    path.write_text(
        "".join(_atom(i, name, res, elem) for i, (name, res, elem) in enumerate(records, start=1))
        + "END\n",
        encoding="utf-8",
    )
    return path


def _elements(path: Path) -> list[str]:
    return [
        line[76:78]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]


def test_reader_uses_a_real_element_column_only_when_the_name_yields_x(tmp_path: Path) -> None:
    from pysisyphus.io.pdb import parse_pdb

    records = [
        (" X1 ", "LIG", "CL"),
        (" X2 ", "LIG", "Fe"),
        (" XE ", "LIG", "XE"),
        (" X3 ", "LIG", ""),
        (" X4 ", "LIG", "X"),
        (" X5 ", "LIG", "EP"),
        # Names that determine an element still win over a corrupted column.
        (" CA ", "ALA", "CA"),
        (" NH1", "ARG", "NH"),
        ("ZN  ", "ZN", "N"),
    ]
    atoms, *_ = parse_pdb(str(_write(tmp_path / "x.pdb", records)))
    assert atoms == ["Cl", "Fe", "Xe", "X", "X", "X", "C", "N", "Zn"]


def test_default_keeps_valid_columns_and_fills_empty_or_invalid_ones(tmp_path: Path) -> None:
    source = _write(
        tmp_path / "mixed.pdb",
        [
            (" CA ", "ALA", ""),
            (" X1 ", "LIG", "CL"),
            (" CA ", "ALA", "CA"),
            (" N  ", "ALA", "XX"),
            (" EPW", "HOH", "EP"),
        ],
    )
    before = source.read_text(encoding="utf-8").splitlines()
    target = tmp_path / "fixed.pdb"

    result = CliRunner().invoke(
        root_cli, ["add-elem-info", "-i", str(source), "-o", str(target)], catch_exceptions=False
    )

    assert result.exit_code == 0, result.output
    after = target.read_text(encoding="utf-8").splitlines()
    assert _elements(target) == [" C", "CL", "CA", " N", "EP"]
    assert [after[i] for i in (1, 2, 4)] == [before[i] for i in (1, 2, 4)]
    assert "assigned/updated            : 2" in result.output
    assert "kept existing               : 3" in result.output


@pytest.mark.parametrize(
    "flag",
    [["--overwrite-elem"], ["--overwrite-elem", "true"], ["--overwrite-elem", "True"]],
)
def test_overwrite_elem_reinfers_valid_columns(tmp_path: Path, flag: list[str]) -> None:
    source = _write(tmp_path / "corrupt.pdb", [(" CA ", "ALA", "CA"), (" CD ", "ARG", "CD")])
    target = tmp_path / "fixed.pdb"

    result = CliRunner().invoke(
        root_cli, ["add-elem-info", "-i", str(source), "-o", str(target), *flag]
    )

    assert result.exit_code == 0, result.output
    assert _elements(target) == [" C", " C"]
    assert "kept existing               : 0" in result.output


@pytest.mark.parametrize("flag", [[], ["--no-overwrite-elem"], ["--overwrite-elem", "false"]])
def test_overwrite_elem_off_keeps_valid_columns(tmp_path: Path, flag: list[str]) -> None:
    source = _write(tmp_path / "corrupt.pdb", [(" CA ", "ALA", "CA"), (" CD ", "ARG", "CD")])
    target = tmp_path / "fixed.pdb"

    result = CliRunner().invoke(
        root_cli, ["add-elem-info", "-i", str(source), "-o", str(target), *flag]
    )

    assert result.exit_code == 0, result.output
    assert target.read_bytes() == source.read_bytes()
    assert "kept existing               : 2" in result.output


@pytest.mark.parametrize("flag", ["--overwrite", "--inplace"])
def test_overwrite_still_means_replacing_the_input_file(tmp_path: Path, flag: str) -> None:
    source = _write(tmp_path / "enzyme.pdb", [(" CA ", "ALA", "CA"), (" CB ", "ALA", "")])

    result = CliRunner().invoke(root_cli, ["add-elem-info", "-i", str(source), flag])

    assert result.exit_code == 0, result.output
    assert "Wrote: " + str(source) in result.output
    assert _elements(source) == ["CA", " C"]


def test_python_api_takes_overwrite_elem(tmp_path: Path) -> None:
    source = _write(tmp_path / "corrupt.pdb", [(" CA ", "ALA", "CA")])
    kept = tmp_path / "kept.pdb"
    redone = tmp_path / "redone.pdb"

    assign_elements(str(source), str(kept))
    assign_elements(str(source), str(redone), overwrite_elem=True)

    assert _elements(kept) == ["CA"]
    assert _elements(redone) == [" C"]


def test_all_preflight_runs_only_for_blank_columns_and_keeps_valid_ones(tmp_path: Path) -> None:
    from mlmm.workflows.all import _assign_elem_info, _pdb_needs_elem_fix

    filled = _write(tmp_path / "filled.pdb", [(" CA ", "ALA", "CA"), (" CB ", "ALA", "C")])
    blank = _write(tmp_path / "blank.pdb", [(" CA ", "ALA", "CA"), (" CB ", "ALA", "")])
    assert not _pdb_needs_elem_fix(filled)
    assert _pdb_needs_elem_fix(blank)

    fixed = tmp_path / "fixed.pdb"
    _assign_elem_info(str(blank), str(fixed), overwrite=False)
    assert _elements(fixed) == ["CA", " C"]


def test_advanced_help_lists_the_overwrite_elem_toggle() -> None:
    result = CliRunner().invoke(root_cli, ["add-elem-info", "--help-advanced"])

    assert result.exit_code == 0, result.output
    assert "--overwrite-elem / --no-overwrite-elem" in result.output


class _FakeMCP:
    def __init__(self) -> None:
        self.tools: dict[str, object] = {}

    def tool(self):
        def decorator(func):
            self.tools[func.__name__] = func
            return func

        return decorator


class _FakeResult:
    def to_dict(self) -> dict:
        return {"status": "ok"}


@pytest.mark.parametrize("value, expected", [(False, False), (True, True)])
def test_mcp_add_element_info_passes_overwrite_elem(monkeypatch, value, expected) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(
        _tools, "run_subcmd", lambda argv, **_kwargs: calls.append(list(argv)) or _FakeResult()
    )
    mcp = _FakeMCP()
    _tools.register_all(mcp)

    asyncio.run(mcp.tools["add_element_info"]("raw.pdb", "fixed.pdb", overwrite_elem=value))

    assert ("--overwrite-elem" in calls[-1]) is expected
    assert "--overwrite" not in calls[-1]
