"""``all --print-every`` reaches optimizing children as a CLI flag, not via YAML.

path-opt and path-search accept a hidden ``--print-every`` that only sets the
single-structure optimizers (LBFGS/RFO) and wins over YAML without a conflict.
"""

from __future__ import annotations

import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml
from ase import Atoms
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.core import utils
from mlmm.core.result_commit import MLMM_RUN_ID_ENV
from mlmm.workflows import all as all_workflow
from mlmm.workflows import dft
from mlmm.workflows import path_opt as path_opt_workflow

SMOKE = Path(__file__).resolve().parent / "smoke"


class ReachedChild(BaseException):
    """Escape the parent pipeline at the first child dispatch."""


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.next")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)
    return path


def _coords(index: int):
    return [(index * .125, .25, -.5), (1.25 + index * .25, .5, .75)]


def _pdb(index: int) -> str:
    return "".join(
        f"HETATM{serial:5d} {name:<4s} MOL A   1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}          {element:>2s}\n"
        for serial, (name, element, (x, y, z)) in enumerate(
            zip(("C1", "O1"), ("C", "O"), _coords(index)), 1
        )
    ) + "END\n"


def _all_args(tmp_path: Path, monkeypatch, n_inputs: int, extra, print_every):
    monkeypatch.delenv(MLMM_RUN_ID_ENV, raising=False)
    monkeypatch.setattr(utils, "_CONVERT_FILES_ENABLED", True)
    inputs = [_write(tmp_path / f"input{i}.pdb", _pdb(i)) for i in range(n_inputs)]
    parm = _write(tmp_path / "unused.parm7", "dispatch fixture; not a force field\n")
    config = _write(tmp_path / "config.yaml", "lbfgs:\n  print_every: 50\n")
    model = Atoms("CO", positions=_coords(0))
    workspace = SimpleNamespace(
        atoms_model=model, atoms_model_lh=model.copy(), model_pdb=inputs[0],
        link_pairs=[], cleanup=lambda: None,
    )
    monkeypatch.setattr(dft, "_prepare_ml_region_workspace", lambda **kwargs: workspace)

    def no_calculator(*args, **kwargs):
        raise AssertionError("No calculator may be built before the child dispatch")

    monkeypatch.setattr(all_workflow, "_mlmm_calc", no_calculator)
    args = ["all", *[arg for path in inputs for arg in ("-i", str(path))]]
    args += ["-q", "0", "-m", "1", "--parm", str(parm), "--model-pdb", str(inputs[0]),
             "--out-dir", str(tmp_path / "out"), "--no-preopt",
             "--config", str(config), *extra]
    if print_every is not None:
        args += ["--print-every", str(print_every)]
    return args


def _assert_flag(child_args, print_every):
    if print_every is None:
        assert "--print-every" not in child_args
    else:
        assert child_args.count("--print-every") == 1
        assert child_args[child_args.index("--print-every") + 1] == str(print_every)


_CASES = {
    "path_opt": (2, []),
    "path_search": (2, ["--refine-path"]),
    "scan": (1, ["--scan-lists", "[(1,2,1.50)]"]),
    "tsopt": (1, ["--tsopt"]),
}


@pytest.mark.parametrize("print_every", [7, None], ids=["explicit", "omitted"])
@pytest.mark.parametrize("child_name", list(_CASES))
def test_all_forwards_print_every_as_child_cli_flag(
    tmp_path, monkeypatch, child_name, print_every
):
    n_inputs, extra = _CASES[child_name]
    args = _all_args(tmp_path, monkeypatch, n_inputs, extra, print_every)
    captured = {}

    def child(name, _cli, child_args, **kwargs):
        captured["name"] = name
        captured["args"] = list(child_args)
        if "--config" in child_args:
            captured["config"] = yaml.safe_load(
                Path(child_args[child_args.index("--config") + 1]).read_text(encoding="utf-8")
            )
        raise ReachedChild

    monkeypatch.setattr(all_workflow, "_run_cli_main", child)
    with pytest.raises(ReachedChild):
        CliRunner().invoke(root_cli, args, catch_exceptions=False)

    assert captured["name"] == child_name
    _assert_flag(captured["args"], print_every)
    forwarded = captured["config"]
    assert "print_every" not in (forwarded.get("opt") or {})
    assert forwarded["lbfgs"]["print_every"] == 50


@pytest.mark.parametrize("print_every", [7, None], ids=["explicit", "omitted"])
def test_all_tsonly_endpoint_opt_receives_print_every(tmp_path, monkeypatch, print_every):
    args = _all_args(tmp_path, monkeypatch, 1, ["--tsopt"], print_every)

    def geom(energy):
        return SimpleNamespace(energy=energy, cart_coords=np.asarray(_coords(0)).reshape(-1))

    def fake_tsopt(hei_pdb, *a, **k):
        ts_pdb = Path(hei_pdb)
        g_ts = geom(-1.0)
        g_ts._tsopt_result = {"energy_hartree": -1.0, "n_imaginary_modes": 1}
        g_ts._tsopt_continuation = {"continue_irc": True, "reaction_mode_index": 0}
        g_ts._tsopt_result_path = None
        return ts_pdb, None, g_ts

    def fake_irc(**kwargs):
        return {"left_min_geom": geom(-1.2), "right_min_geom": geom(-1.1), "ts_geom": geom(-1.0)}

    def fake_save(g, ref_pdb, out_dir, name):
        return Path(out_dir) / f"{name}.xyz", Path(out_dir) / f"{name}.pdb"

    captured = {}

    def fake_opt(*a, **kwargs):
        captured.update(kwargs)
        raise ReachedChild

    monkeypatch.setattr(all_workflow, "_run_tsopt_on_hei", fake_tsopt)
    monkeypatch.setattr(all_workflow, "_irc_and_match", fake_irc)
    monkeypatch.setattr(all_workflow, "_save_single_geom_for_tools", fake_save)
    monkeypatch.setattr(all_workflow, "_run_opt_for_state", fake_opt)
    with pytest.raises(ReachedChild):
        CliRunner().invoke(root_cli, args, catch_exceptions=False)

    assert captured["print_every"] == print_every


@pytest.mark.parametrize("print_every", [7, None], ids=["explicit", "omitted"])
def test_endpoint_opt_child_argv_carries_print_every(tmp_path, monkeypatch, print_every):
    pdb = _write(tmp_path / "state.pdb", _pdb(0))
    captured = {}

    def child(name, _cli, child_args, **kwargs):
        captured["name"] = name
        captured["args"] = list(child_args)
        raise ReachedChild

    monkeypatch.setattr(all_workflow, "_run_cli_main", child)
    with pytest.raises(ReachedChild):
        all_workflow._run_opt_for_state(
            pdb, 0, 1, tmp_path / "unused.parm7", pdb, False, tmp_path / "opt",
            None, "grad", resolved_calc_template=None, print_every=print_every,
        )

    assert captured["name"] == "opt"
    _assert_flag(captured["args"], print_every)


def _path_invoke(tmp_path: Path, command: str, yaml_text: str, *extra: str):
    config = _write(tmp_path / "config.yaml", yaml_text)
    return CliRunner().invoke(
        root_cli,
        [
            command,
            "-i", str(SMOKE / "r_complex_layered.pdb"), str(SMOKE / "p_complex_layered.pdb"),
            "--parm", str(SMOKE / "p_complex.parm7"),
            "-q", "-1", "-m", "1", "--config", str(config), "--show-config", "-v", "3",
            "--out-dir", str(tmp_path / command), *extra,
        ],
    )


def _block(output: str, title: str) -> dict:
    head = f"\n{title}\n{'-' * len(title)}\n"
    start = output.index(head) + len(head)
    end = output.find("\n\n", start)
    return yaml.safe_load(output[start:] if end < 0 else output[start:end])


def _resolved_blocks(tmp_path, monkeypatch, command, yaml_text, *extra):
    """Resolved ``stopt`` and ``opt.<kind>`` blocks as the command prints them."""
    if command == "path-search":
        result = _path_invoke(tmp_path, command, yaml_text, "--dry-run", *extra)
        assert result.exit_code == 0, result.output
        return lambda title: _block(result.output, title)

    blocks = {}
    real_pretty_block = path_opt_workflow.pretty_block

    def recording_pretty_block(title, payload, *args, **kwargs):
        blocks[title] = dict(payload) if isinstance(payload, dict) else payload
        return real_pretty_block(title, payload, *args, **kwargs)

    def stop_at_calculator(**kwargs):
        raise ReachedChild

    monkeypatch.setattr(path_opt_workflow, "pretty_block", recording_pretty_block)
    monkeypatch.setattr(path_opt_workflow, "mlmm", stop_at_calculator)
    with pytest.raises(ReachedChild):
        CliRunner().invoke(
            root_cli,
            [
                command,
                "-i", str(SMOKE / "r_complex_layered.pdb"), str(SMOKE / "p_complex_layered.pdb"),
                "--parm", str(SMOKE / "p_complex.parm7"),
                "-q", "-1", "-m", "1",
                "--config", str(_write(tmp_path / "config.yaml", yaml_text)),
                "--out-dir", str(tmp_path / command), *extra,
            ],
            catch_exceptions=False,
        )
    return blocks.__getitem__


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
@pytest.mark.parametrize(("opt_mode", "kind"), [("grad", "lbfgs"), ("hess", "rfo")])
@pytest.mark.parametrize(
    "yaml_text",
    ["{}\n", "opt:\n  print_every: 50\n", "lbfgs:\n  print_every: 50\nrfo:\n  print_every: 50\n"],
    ids=["no-yaml", "opt-yaml", "optimizer-yaml"],
)
def test_path_hidden_print_every_sets_single_structure_optimizers(
    tmp_path, monkeypatch, command, opt_mode, kind, yaml_text
):
    block = _resolved_blocks(
        tmp_path, monkeypatch, command, yaml_text,
        "--opt-mode", opt_mode, "--print-every", "7",
    )

    assert block(f"opt.{kind}")["print_every"] == 7
    assert (block("stopt") or {}).get("print_every", 10) == 10


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
def test_path_print_every_without_cli_keeps_yaml_value(tmp_path, monkeypatch, command):
    block = _resolved_blocks(tmp_path, monkeypatch, command, "lbfgs:\n  print_every: 50\n")

    assert block("opt.lbfgs")["print_every"] == 50


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
def test_path_print_every_is_hidden_from_help(command):
    for flag in ("--help", "--help-advanced"):
        result = CliRunner().invoke(root_cli, [command, flag])
        assert result.exit_code == 0, result.output
        assert "--print-every" not in result.output
