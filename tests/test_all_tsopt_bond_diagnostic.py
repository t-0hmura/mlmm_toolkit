"""Actual parent CLI dispatch with current-run MEP artifacts and no backend.

The child/workspace fixtures follow test_all_mep_pdb.py. A sentinel at TSOPT
entry tests selection, without claiming an ML/MM TS calculation succeeded.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest
from ase import Atoms
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.core import utils
from mlmm.core.result_commit import MLMM_RUN_ID_ENV, with_current_run_id
from mlmm.workflows import all as all_workflow
from mlmm.workflows import dft


class ReachedTSDispatch(BaseException):
    """Escape before backend setup, including outer ``except Exception``."""


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.next")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)
    return path


def _coords(index):
    return [(index * .125, .25, -.5), (1.25 + index * .25, .5, .75)]


def _frame(index):
    atoms = "".join(
        f"{element} {x:.3f} {y:.3f} {z:.3f}\n"
        for element, (x, y, z) in zip(("C", "O"), _coords(index))
    )
    # The HEI is strictly internal; endpoint-HEI handling is not under test.
    return f"2\nE={(0., .1, .02)[index]:.3f} unit=hartree frame{index}\n{atoms}"


def _pdb(index):
    return "".join(
        f"HETATM{serial:5d} {name:<4s} MOL A   1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}          {element:>2s}\n"
        for serial, (name, element, (x, y, z)) in enumerate(
            zip(("C1", "O1"), ("C", "O"), _coords(index)), 1
        )
    ) + "END\n"


@pytest.mark.parametrize(
    "bond_changed,do_tsopt",
    [
        pytest.param(False, True, id="no-bond-change-still-dispatches"),
        pytest.param(True, True, id="bond-change-ordinary-control"),
        pytest.param(False, False, id="explicit-tsopt-opt-out"),
    ],
)
def test_all_tsopt_dispatch_ignores_bond_diagnostic(tmp_path, monkeypatch, bond_changed, do_tsopt):
    monkeypatch.delenv(MLMM_RUN_ID_ENV, raising=False)
    monkeypatch.setattr(utils, "_CONVERT_FILES_ENABLED", True)
    out = tmp_path / "out"
    inputs = [_write(tmp_path / f"input{i}.pdb", _pdb(i)) for i in (0, 2)]
    input_bytes = [path.read_bytes() for path in inputs]
    # Only initial ML/MM workspace construction is replaced. This fixture is
    # not a force field and must never reach a calculator or topology reader.
    parm = _write(tmp_path / "unused.parm7", "dispatch fixture; not a force field\n")
    model = Atoms("CO", positions=_coords(0))
    workspace = SimpleNamespace(
        atoms_model=model, atoms_model_lh=model.copy(), model_pdb=inputs[0],
        link_pairs=[], cleanup=lambda: None,
    )
    monkeypatch.setattr(dft, "_prepare_ml_region_workspace", lambda **kwargs: workspace)
    child_calls = []
    bond_calls = []
    ts_calls = []
    calculator_calls = []

    def child(name, _cli, args, **kwargs):
        assert name == "path_opt"
        child_calls.append(name)
        child_out = Path(args[args.index("--out-dir") + 1])
        _write(child_out / "final_geometries_trj.xyz", "".join(_frame(i) for i in range(3)))
        _write(child_out / "hei.xyz", _frame(1))
        _write(child_out / "hei.pdb", _pdb(1))
        payload = {
            "status": "converged", "converged": True,
            "preopt_requested": False, "preopt_converged": None,
            "path_optimizers": [],
            "stage_outcomes": [{"stage": "path-opt", "item_id": "gsm_mep",
                                "executed": True, "converged": True,
                                "usable": True, "required": True}],
        }
        _write(child_out / "result.json", json.dumps(with_current_run_id(payload)))

    def bond_diagnostic(*args, **kwargs):
        bond_calls.append(True)
        return bond_changed, "Bond formed" if bond_changed else ""

    def diagram(prefix, *, labels, energies_eh, **kwargs):
        return {"name": prefix.stem, "labels": labels, "energies_au": list(energies_eh)}

    def no_calculator(*args, **kwargs):
        calculator_calls.append(True)
        raise AssertionError("The TS entry sentinel must precede calculator construction")

    def stop_at_ts(hei, *args, **kwargs):
        hei = Path(hei)
        assert hei.is_file()
        assert hei.stem == "hei_seg_01"
        ts_calls.append(hei)
        raise ReachedTSDispatch

    monkeypatch.setattr(all_workflow, "_run_cli_main", child)
    monkeypatch.setattr(all_workflow, "_mlmm_calc", no_calculator)
    monkeypatch.setattr(all_workflow, "run_trj2fig", lambda *args, **kwargs: None)
    monkeypatch.setattr(all_workflow, "close_matplotlib_figures", lambda: None)
    monkeypatch.setattr(all_workflow, "_write_segment_energy_diagram", diagram)
    monkeypatch.setattr(all_workflow._path_search, "_has_bond_change", bond_diagnostic)
    monkeypatch.setattr(all_workflow, "_run_tsopt_on_hei", stop_at_ts)
    args = ["all", *[arg for path in inputs for arg in ("-i", str(path))]]
    args += ["-q", "0", "-m", "1", "--parm", str(parm), "--model-pdb", str(inputs[0]),
             "--out-dir", str(out), "--no-preopt", "--convert-files",
             "--no-tsopt-from-mep-tan", "--tsopt" if do_tsopt else "--no-tsopt"]
    if do_tsopt:
        with pytest.raises(ReachedTSDispatch):
            CliRunner().invoke(root_cli, args, catch_exceptions=False)
    else:
        result = CliRunner().invoke(root_cli, args)
        assert result.exit_code == 0, result.output + repr(result.exception)

    assert child_calls == ["path_opt"]
    assert bond_calls == [True]
    assert len(ts_calls) == int(do_tsopt)
    assert calculator_calls == []
    assert [path.read_bytes() for path in inputs] == input_bytes
    summary = json.loads((out / "summary.json").read_text())
    segment, = summary["segments"]
    assert summary["n_segments"] == 1 and summary["n_images"] == 3
    assert segment["kind"] == "seg"
    assert segment["converged"] is True
    expected = "Bond formed" if bond_changed else "(no covalent changes detected)"
    assert expected in str(segment["bond_changes"])
    assert (out / "_work/path_opt/hei_seg_01.pdb").is_file()
    assert (out / "_work/path_opt/mep_seg_01_trj.xyz").is_file()
