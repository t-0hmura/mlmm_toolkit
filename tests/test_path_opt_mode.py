"""Hidden ``--opt-mode`` of path-opt/path-search: grad (LBFGS) or hess (RFO)."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import click
import pytest
from ase import Atoms
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.workflows import path_opt as path_opt_workflow
from mlmm.workflows import path_search as path_search_workflow

SMOKE = Path(__file__).resolve().parent / "smoke"
COMMANDS = ("path-opt", "path-search")
MODULES = {"path-opt": path_opt_workflow, "path-search": path_search_workflow}


def _has_option_row(output: str, option: str) -> bool:
    for line in output.splitlines():
        stripped = line.lstrip()
        if len(line) - len(stripped) > 2 or not stripped.startswith(option):
            continue
        tail = stripped[len(option):]
        if not tail or tail[0].isspace() or tail[0] in {",", "/"}:
            return True
    return False


def _dry_run(tmp_path: Path, monkeypatch, command: str, *extra: str, config: str | None = None):
    blocks = {}

    def capture(title, content, **_kwargs):
        blocks[title] = dict(content)
        return ""

    monkeypatch.setattr(MODULES[command], "pretty_block", capture)
    args = [
        command,
        "-i", str(SMOKE / "r_complex_layered.pdb"), str(SMOKE / "p_complex_layered.pdb"),
        "--parm", str(SMOKE / "p_complex.parm7"),
        "-q", "-1", "-m", "1",
        "--dry-run",
        "--out-dir", str(tmp_path / command),
        *extra,
    ]
    if config is not None:
        config_path = tmp_path / "config.yaml"
        config_path.write_text(config, encoding="utf-8")
        args += ["--config", str(config_path)]
    return CliRunner().invoke(root_cli, args), blocks


@pytest.mark.parametrize("command", COMMANDS)
def test_opt_mode_is_advanced_only_and_defaults_to_grad(command: str) -> None:
    runner = CliRunner()
    primary = runner.invoke(root_cli, [command, "--help"])
    advanced = runner.invoke(root_cli, [command, "--help-advanced"])
    assert primary.exit_code == 0, primary.output
    assert advanced.exit_code == 0, advanced.output
    assert not _has_option_row(primary.output, "--opt-mode")
    assert _has_option_row(advanced.output, "--opt-mode")

    subcommand = root_cli.get_command(click.Context(root_cli), command)
    option = next(p for p in subcommand.params if "--opt-mode" in p.opts)
    assert option.default == "grad"
    assert list(option.type.choices) == ["grad", "hess"]


@pytest.mark.parametrize("command", COMMANDS)
@pytest.mark.parametrize(("mode", "kind"), [("grad", "lbfgs"), ("hess", "rfo")])
def test_dry_run_records_the_selected_mode(tmp_path, monkeypatch, command, mode, kind) -> None:
    result, blocks = _dry_run(tmp_path, monkeypatch, command, "--opt-mode", mode)

    assert result.exit_code == 0, result.output
    assert blocks["dry_run_plan"]["opt_mode"] == mode
    if command == "path-search":
        assert f"opt.{kind}" in blocks
        assert blocks[f"opt.{kind}"]["out_dir_per_tag"].endswith(f"<tag>_{kind}_opt")


def test_path_search_show_config_prints_the_rfo_block(tmp_path: Path) -> None:
    result = CliRunner().invoke(
        root_cli,
        [
            "path-search",
            "-i", str(SMOKE / "r_complex_layered.pdb"), str(SMOKE / "p_complex_layered.pdb"),
            "--parm", str(SMOKE / "p_complex.parm7"),
            "-q", "-1", "-m", "1",
            "--opt-mode", "hess", "--dry-run", "--show-config",
            "--out-dir", str(tmp_path / "path-search"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "\nopt.rfo\n-------\n" in result.output
    assert "\nopt.lbfgs\n" not in result.output


@pytest.mark.parametrize("command", COMMANDS)
@pytest.mark.parametrize(
    "config",
    ["rfo:\n  thresh: gau_tight\n", "opt:\n  rfo:\n    thresh: gau_tight\n"],
)
def test_yaml_rfo_sections_reach_only_the_rfo_run(tmp_path, monkeypatch, command, config) -> None:
    captured = []
    real_apply = MODULES[command].apply_single_opt_yaml_layer

    def record(layer_cfg, **kwargs):
        real_apply(layer_cfg, **kwargs)
        captured.append((dict(kwargs["lbfgs_cfg"]), dict(kwargs["rfo_cfg"])))

    monkeypatch.setattr(MODULES[command], "apply_single_opt_yaml_layer", record)
    result, _ = _dry_run(tmp_path, monkeypatch, command, "--opt-mode", "hess", config=config)

    assert result.exit_code == 0, result.output
    lbfgs_cfg, rfo_cfg = captured[-1]
    assert rfo_cfg["thresh"] == "gau_tight"
    assert lbfgs_cfg["thresh"] != "gau_tight"


@pytest.mark.parametrize("command", COMMANDS)
def test_yaml_conflict_is_checked_for_the_optimizer_that_runs(tmp_path, monkeypatch, command) -> None:
    config = "opt:\n  thresh: gau\nrfo:\n  thresh: gau_tight\n"

    rejected, _ = _dry_run(tmp_path, monkeypatch, command, "--opt-mode", "hess", config=config)
    assert rejected.exit_code != 0
    assert "opt.thresh and rfo.thresh conflict" in rejected.output

    accepted, _ = _dry_run(tmp_path, monkeypatch, command, "--opt-mode", "grad", config=config)
    assert accepted.exit_code == 0, accepted.output


def test_rfo_starts_from_the_calculator_hessian(monkeypatch) -> None:
    seeded, built = [], []
    monkeypatch.setattr(
        path_opt_workflow, "seed_scan_rfo_hessian",
        lambda geom, calc_cfg: seeded.append((geom, dict(calc_cfg))),
    )
    monkeypatch.setattr(
        path_opt_workflow, "RFOptimizer",
        lambda geom, **kwargs: built.append(("rfo", kwargs)) or "rfo",
    )
    monkeypatch.setattr(
        path_opt_workflow, "LBFGS",
        lambda geom, **kwargs: built.append(("lbfgs", kwargs)) or "lbfgs",
    )
    geom = object()
    make = path_opt_workflow._make_single_optimizer

    assert make(geom, "rfo", {"hessian_init": "calc", "max_cycles": 3}, {"model_charge": -1}) == "rfo"
    assert seeded == [(geom, {"model_charge": -1})]
    assert built[-1] == ("rfo", {"hessian_init": "calc", "max_cycles": 3})

    assert make(geom, "rfo", {"hessian_init": "unit"}, None) == "rfo"
    assert make(geom, "lbfgs", {"max_cycles": 3}, None) == "lbfgs"
    assert len(seeded) == 1
    with pytest.raises(ValueError, match="calculator settings"):
        make(geom, "rfo", {"hessian_init": "calc"}, None)
    with pytest.raises(ValueError, match="Unknown"):
        make(geom, "dimer", {}, None)


def test_path_search_single_optimization_uses_the_selected_optimizer(tmp_path, monkeypatch) -> None:
    built = []

    class FakeOptimizer:
        is_converged = True
        final_fn = tmp_path / "missing.xyz"

        def run(self) -> None:
            return None

    def make(geom, kind, args, calc_cfg):
        built.append((kind, args["out_dir"], calc_cfg))
        return FakeOptimizer()

    monkeypatch.setattr(path_search_workflow, "_make_single_optimizer", make)
    geom = SimpleNamespace(set_calculator=lambda _calc: None, coord_type="cart")

    result, converged = path_search_workflow._optimize_single(
        geom, None, "rfo", {"thresh": "gau"}, tmp_path,
        tag="seg_000_left", ref_pdb_path=None, calc_cfg={"model_charge": -1},
    )

    assert result is geom and converged is True
    assert built == [("rfo", str(tmp_path / "seg_000_left_rfo_opt"), {"model_charge": -1})]


def test_recursive_search_passes_the_optimizer_to_hei_neighbours(tmp_path, monkeypatch) -> None:
    def mep():
        energies = [0.0, 1.0, 0.0]
        return SimpleNamespace(
            images=[SimpleNamespace(energy=e) for e in energies],
            energies=energies, hei_idx=1, is_converged=True,
        )

    calls = []

    def optimize(geometry, _shared_calc, kind, _cfg, _out_dir, **kwargs):
        calls.append((kind, kwargs["calc_cfg"]))
        return geometry, True

    changes = iter([True, True, False, False])
    monkeypatch.setattr(path_search_workflow, "_run_mep_between", lambda *_a, **_k: mep())
    monkeypatch.setattr(path_search_workflow, "_refine_between", lambda *_a, **_k: mep())
    monkeypatch.setattr(path_search_workflow, "_optimize_single", optimize)
    monkeypatch.setattr(
        path_search_workflow, "_has_bond_change",
        lambda *_a, **_k: (next(changes), "Bond formed"),
    )
    monkeypatch.setattr(path_search_workflow, "_stitch_paths", lambda parts, **_k: parts[0])
    primary = mep()

    result = path_search_workflow._build_multistep_path(
        primary.images[0], primary.images[-1], None,
        geom_cfg={}, gs_cfg={}, stopt_cfg={}, single_opt_kind="rfo", single_opt_cfg={},
        bond_cfg={}, search_cfg={"max_depth": 1, "stitch_rmsd_thresh": 1e-4,
                                 "bridge_rmsd_thresh": 1e-4},
        refine_mode_kind="peak", mep_mode_kind="gsm", out_dir=tmp_path,
        ref_pdb_path=None, depth=0, seg_counter=[0], branch_tag="pair_00",
        calc_cfg={"model_charge": -1}, dmf_cfg={},
    )

    assert result.single_opt_executed is True
    assert calls == [("rfo", {"model_charge": -1})] * 2


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _pdb(x0: float) -> str:
    return "".join(
        f"HETATM{serial:5d} {name:<4s} LIG A   1    "
        f"{x:8.3f}{8.0:8.3f}{9.0:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
        for serial, name, element, x in ((1, "C1", "C", x0), (2, "O1", "O", x0 + 1.2))
    ) + "END\n"


@pytest.mark.parametrize(("extra", "mode"), [([], "grad"), (["--opt-mode", "hess"], "hess")])
def test_all_forwards_opt_mode_to_path_children(tmp_path, monkeypatch, extra, mode) -> None:
    from mlmm.core.result_commit import MLMM_RUN_ID_ENV, with_current_run_id
    from mlmm.workflows import all as all_workflow
    from mlmm.workflows import dft

    monkeypatch.delenv(MLMM_RUN_ID_ENV, raising=False)
    inputs = [_write(tmp_path / f"input{i}.pdb", _pdb(8.0 + 0.1 * i)) for i in range(2)]
    parm = _write(tmp_path / "unused.parm7", "fixture; never parsed\n")
    model = Atoms("CO", positions=[(8.0, 8.0, 9.0), (9.2, 8.0, 9.0)])
    workspace = SimpleNamespace(
        atoms_model=model, atoms_model_lh=model.copy(), model_pdb=inputs[0],
        link_pairs=[], cleanup=lambda: None,
    )
    monkeypatch.setattr(dft, "_prepare_ml_region_workspace", lambda **_k: workspace)
    forwarded = []

    def fake_child(name, _cli, args, **_kwargs):
        assert name == "path_opt", name
        forwarded.append(args[args.index("--opt-mode") + 1])
        child_out = Path(args[args.index("--out-dir") + 1])
        frame = "2\nE=0.0 unit=hartree\nC 8.0 8.0 9.0\nO 9.2 8.0 9.0\n"
        _write(child_out / "final_geometries_trj.xyz", frame + frame)
        payload = {"stage_outcomes": [{"item_id": "gsm_mep", "converged": True}]}
        _write(child_out / "result.json", json.dumps(with_current_run_id(payload)))

    monkeypatch.setattr(all_workflow, "_run_cli_main", fake_child)
    monkeypatch.setattr(all_workflow, "run_trj2fig", lambda *_a, **_k: None)
    monkeypatch.setattr(all_workflow, "close_matplotlib_figures", lambda: None)
    monkeypatch.setattr(all_workflow, "_write_segment_energy_diagram", lambda *_a, **_k: None)
    args = [arg for path in inputs for arg in ("-i", str(path))]
    args += ["-q", "0", "-m", "1", "--parm", str(parm), "--model-pdb", str(inputs[0]),
             "--no-preopt", "--no-convert-files", "--out-dir", str(tmp_path / "out"), *extra]

    result = CliRunner().invoke(all_workflow.cli, args)

    assert forwarded == [mode], f"{result.output}\n{result.exception!r}"
