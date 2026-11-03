"""``--compile/--no-compile``: the UMA and ORB ML region runs uncompiled unless it is asked for.

The flag and ``calc.compile_model`` resolve as one key (both given and
different is an error), other backends reject ``--compile``, the analytical
Hessian is never compiled, and a compiled run without a compiler stops with
the install guidance before the model loads. No model is loaded here.
"""

from __future__ import annotations

import sys
from pathlib import Path
from types import ModuleType

import click
import pytest
import torch
import yaml
from click.testing import CliRunner

import mlmm.backends.mlmm_calc as mlmm_calc
from mlmm.cli import cli as root_cli
from mlmm.core import compile_policy
from mlmm.core import prerequisites as pre
from mlmm.core.compile_policy import (
    apply_compile_option,
    find_compiler,
    resolve_compile,
    runs_compiled,
)
from mlmm.core.utils import calculator_provenance, calculator_run_label

HELP = (
    "torch.compile the UMA or ORB model: slower first call, slightly faster "
    "afterwards; worth it only for long runs."
)

_TOY = Path(__file__).resolve().parents[1] / "examples" / "toy_system"
_R, _PARM = str(_TOY / "r_toy.pdb"), str(_TOY / "p_toy.parm7")
_SYSTEM = ["--parm", _PARM, "-q", "-1", "-m", "1"]

needs_toy_system = pytest.mark.skipif(
    not (_TOY / "p_toy.parm7").exists(), reason="toy_system example inputs not present"
)


@pytest.fixture(autouse=True)
def _fresh_state(monkeypatch):
    """Reset the once-per-process note and run at the default verbosity with the console gate off."""
    from mlmm.core import utils

    monkeypatch.setattr(compile_policy, "_analytical_noted", False)
    for module, name, value in (
        (utils, "_VERBOSE_LEVEL", 2),
        (utils, "_GATE_ACTIVE", False),
        (utils, "_PIPELINE_MODE", False),
        (utils, "_CHILD_MODE", False),
    ):
        monkeypatch.setattr(module, name, value)


def _which_from(available, probed=None):
    def which(name):
        if probed is not None:
            probed.append(name)
        return f"/usr/bin/{name}" if name in available else None

    return which


# ---------------------------------------------------------------------------
# find_compiler
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("device", "variable", "order"),
    [("cpu", "CXX", ["c++", "g++", "clang++"]), ("cuda", "CC", ["cc", "gcc", "clang"])],
)
def test_find_compiler_probes_the_variable_then_the_names(monkeypatch, device, variable, order):
    monkeypatch.delenv("CC", raising=False)
    monkeypatch.delenv("CXX", raising=False)
    probed = []
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from(set(), probed))
    assert find_compiler(device) is None
    assert probed == order

    monkeypatch.setenv(variable, "my-compiler")
    probed.clear()
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from({"my-compiler", order[-1]}, probed))
    assert find_compiler(device) == "/usr/bin/my-compiler"
    assert probed == ["my-compiler"]


@pytest.mark.parametrize(
    ("device", "found"),
    [("cpu", "clang++"), (torch.device("cpu"), "clang++"), (torch.device("cuda", 0), "clang")],
)
def test_find_compiler_takes_the_first_name_on_path(monkeypatch, device, found):
    monkeypatch.delenv("CC", raising=False)
    monkeypatch.delenv("CXX", raising=False)
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from({found}))
    assert find_compiler(device) == f"/usr/bin/{found}"


# ---------------------------------------------------------------------------
# resolve_compile
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("requested", [None, False])
def test_unrequested_compile_is_off_for_any_backend(monkeypatch, requested):
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from(set()))
    for backend in ("uma", "orb", "mace", "dft"):
        assert resolve_compile(requested, backend=backend, device="cpu", analytical_hessian=False) is False


@pytest.mark.parametrize("backend", ["mace", "aimnet2", "dft", "custom"])
def test_compile_for_another_backend_is_a_usage_error(backend):
    with pytest.raises(click.UsageError, match="^--compile applies to the uma and orb backends$"):
        resolve_compile(True, backend=backend, device="cpu", analytical_hessian=False)


@pytest.mark.parametrize("backend", ["uma", "orb"])
def test_analytical_hessian_is_not_compiled_and_says_so_once(monkeypatch, capsys, backend):
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from(set()))
    for _ in range(2):
        assert resolve_compile(True, backend=backend, device="cpu", analytical_hessian=True) is False
    err = capsys.readouterr().err
    assert err.count("[compile] --compile is ignored for the analytical Hessian.") == 1


@pytest.mark.parametrize(("device", "name"), [("cpu", "g++"), (torch.device("cuda", 0), "gcc")])
def test_compile_with_a_compiler_is_on(monkeypatch, device, name):
    monkeypatch.delenv("CC", raising=False)
    monkeypatch.delenv("CXX", raising=False)
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from({name}))
    assert resolve_compile(True, backend="orb", device=device, analytical_hessian=False) is True


@pytest.mark.parametrize(
    ("device", "kind"),
    [(torch.device("cpu"), pre.KIND_CXX_COMPILER), (torch.device("cuda", 0), pre.KIND_C_COMPILER)],
)
def test_compile_without_a_compiler_stops_with_the_guidance(monkeypatch, capsys, device, kind):
    monkeypatch.delenv("CC", raising=False)
    monkeypatch.delenv("CXX", raising=False)
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from(set()))
    with pytest.raises(pre.MissingPrerequisiteError) as missing:
        resolve_compile(True, backend="uma", device=device, analytical_hessian=False)
    assert missing.value.kind == kind
    with pytest.raises(SystemExit) as stop:
        pre.exit_if_missing_prerequisite(missing.value)
    assert stop.value.code == 1
    assert capsys.readouterr().err.strip() == pre.guidance(kind)


# ---------------------------------------------------------------------------
# apply_compile_option: CLI and config
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("cli", "config", "expected"),
    [
        (None, None, False),
        (True, None, True),
        (False, None, False),
        (None, True, True),
        (None, False, False),
        (True, True, True),
        (False, False, False),
    ],
)
def test_one_explicit_value_wins_and_unset_is_off(cli, config, expected):
    calc_cfg = {"backend": "uma"}
    if config is not None:
        calc_cfg["compile_model"] = config
    apply_compile_option(calc_cfg, cli)
    assert calc_cfg["compile_model"] is expected


@pytest.mark.parametrize(("cli", "config"), [(True, False), (False, True)])
def test_cli_and_config_that_differ_are_a_usage_error(cli, config):
    flag = "--compile" if cli else "--no-compile"
    with pytest.raises(click.UsageError, match=f"^{flag} and calc.compile_model: {str(config).lower()} conflict"):
        apply_compile_option({"backend": "orb", "compile_model": config}, cli)


@pytest.mark.parametrize("backend", ["mace", "aimnet2", "dft", "custom"])
def test_other_backends_reject_compile_and_accept_no_compile(backend):
    with pytest.raises(click.UsageError, match="^--compile applies to the uma and orb backends$"):
        apply_compile_option({"backend": backend}, True)
    with pytest.raises(click.UsageError, match="^--compile applies to the uma and orb backends$"):
        apply_compile_option({"backend": backend, "compile_model": True}, None)
    calc_cfg = {"backend": backend}
    apply_compile_option(calc_cfg, False)
    assert calc_cfg["compile_model"] is False


def test_config_value_must_be_a_bool():
    with pytest.raises(click.BadParameter, match="calc.compile_model must be true or false"):
        apply_compile_option({"backend": "uma", "compile_model": "yes"}, None)


# ---------------------------------------------------------------------------
# The flag on every command that has --precision
# ---------------------------------------------------------------------------


def _commands_with_precision():
    import importlib

    from mlmm.cli.app import _LAZY_SUBCOMMANDS

    commands = {}
    for name, (module_name, attr, _help) in _LAZY_SUBCOMMANDS.items():
        command = getattr(importlib.import_module(module_name), attr)
        names = {param.name for param in command.params}
        if "precision" in names:
            commands[name] = command
    return commands


def test_every_command_with_precision_has_compile():
    commands = _commands_with_precision()
    assert {"opt", "tsopt", "freq", "irc", "sp", "scan", "scan2d", "scan3d",
            "path-opt", "path-search", "all", "trj2fig"} <= set(commands)
    for name, command in commands.items():
        option = next((p for p in command.params if p.name == "compile_model"), None)
        assert option is not None, name
        assert option.opts == ["--compile"] and option.secondary_opts == ["--no-compile"], name
        assert option.default is None and option.show_default == "off", name
        assert option.help == HELP, name


# ---------------------------------------------------------------------------
# What ran: the backend line and mlip_compiled
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("calc_cfg", "label", "compiled"),
    [
        ({"backend": "uma", "uma_precision": "fp64", "compile_model": True},
         "UMA (UMA-S-1.2.1 (OMol), fp64, compiled)", True),
        ({"backend": "uma", "uma_precision": "fp64", "compile_model": False},
         "UMA (UMA-S-1.2.1 (OMol), fp64)", False),
        ({"backend": "uma", "uma_precision": "fp64"}, "UMA (UMA-S-1.2.1 (OMol), fp64)", False),
        ({"backend": "uma", "uma_precision": "fp32", "compile_model": True,
          "hessian_calc_mode": "Analytical"}, "UMA (UMA-S-1.2.1 (OMol), fp32)", False),
        ({"backend": "orb", "orb_precision": "float64", "compile_model": True},
         "ORB (OrbMol-v2, fp64, compiled)", True),
        ({"backend": "mace", "mace_dtype": "float64", "compile_model": False},
         "MACE (MACE-OMOL-0, fp64)", False),
    ],
)
def test_backend_line_and_provenance_record_compilation(calc_cfg, label, compiled):
    assert calculator_run_label(calc_cfg) == label
    assert calculator_provenance(calc_cfg)["mlip_compiled"] is compiled
    assert runs_compiled(calc_cfg) is compiled


def test_custom_calculator_is_recorded_uncompiled(tmp_path):
    from mlmm.backends import apply_calc_file_to_calc_cfg

    calc_file = tmp_path / "toy.py"
    calc_file.write_text("def get_calculator(**kwargs):\n    return None\n")
    calc_cfg = {"backend": "uma"}
    apply_calc_file_to_calc_cfg(calc_cfg, str(calc_file), None)
    apply_compile_option(calc_cfg, None)
    assert calculator_provenance(calc_cfg)["mlip_compiled"] is False


@needs_toy_system
def test_opt_dry_run_prints_the_compiled_backend_line(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    base = ["opt", "-i", _R, *_SYSTEM, "--out-dir", str(tmp_path / "out"),
            "--precision", "fp64", "--dry-run"]

    compiled = CliRunner().invoke(root_cli, [*base, "--compile"])
    assert compiled.exit_code == 0, compiled.output
    assert "[backend] UMA (UMA-S-1.2.1 (OMol), fp64, compiled)" in compiled.output

    default = CliRunner().invoke(root_cli, base)
    assert default.exit_code == 0, default.output
    assert "[backend] UMA (UMA-S-1.2.1 (OMol), fp64)" in default.output
    assert "compiled" not in default.output


@needs_toy_system
@pytest.mark.parametrize(
    ("args", "config", "message"),
    [
        (["--compile", "-b", "mace"], None, "--compile applies to the uma and orb backends"),
        (["--compile"], "calc:\n  compile_model: false\n", "--compile and calc.compile_model: false conflict."),
        (["--no-compile"], "calc:\n  compile_model: true\n", "--no-compile and calc.compile_model: true conflict."),
    ],
)
def test_opt_rejects_compile_misuse(tmp_path, monkeypatch, args, config, message):
    monkeypatch.chdir(tmp_path)
    cmd = ["opt", "-i", _R, *_SYSTEM, "--out-dir", str(tmp_path / "out"), "--dry-run", *args]
    if config is not None:
        path = tmp_path / "config.yaml"
        path.write_text(config, encoding="utf-8")
        cmd += ["--config", str(path)]
    result = CliRunner().invoke(root_cli, cmd)
    assert result.exit_code == 2, result.output
    assert message in result.output


@needs_toy_system
def test_no_compile_is_accepted_by_other_backends(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    result = CliRunner().invoke(root_cli, [
        "opt", "-i", _R, *_SYSTEM, "--out-dir", str(tmp_path / "out"),
        "-b", "mace", "--no-compile", "--dry-run",
    ])
    assert result.exit_code == 0, result.output


@needs_toy_system
def test_sp_with_compile_and_no_compiler_stops_before_the_model_loads(monkeypatch, tmp_path):
    def no_load(*args, **kwargs):
        pytest.fail("the model must not load when the compiler is missing")

    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("CXX", raising=False)
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from(set()))
    monkeypatch.setattr(compile_policy, "_on_cuda", lambda device: False)
    monkeypatch.setattr(mlmm_calc, "_UMABackend", no_load)
    result = CliRunner().invoke(
        root_cli,
        ["sp", "-i", _R, *_SYSTEM, "--compile", "--out-dir", str(tmp_path / "out")],
    )
    assert result.exit_code == 1, result.output
    assert "Missing system prerequisite: C++ compiler (g++)." in result.output
    assert "Traceback" not in result.output and "Unhandled" not in result.output


# ---------------------------------------------------------------------------
# Backends
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("backend", "compile_model", "analytical", "expected"),
    [
        ("uma", None, False, False),
        ("uma", False, False, False),
        ("uma", True, False, True),
        ("uma", True, True, False),
        ("orb", None, False, False),
        ("orb", True, False, True),
        ("orb", True, True, False),
    ],
)
def test_factory_passes_the_resolved_choice_to_the_backend(
    monkeypatch, backend, compile_model, analytical, expected
):
    captured = {}

    def fake_backend(**kwargs):
        captured.update(kwargs)
        return object()

    monkeypatch.delenv("CXX", raising=False)
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from({"g++"}))
    monkeypatch.setattr(mlmm_calc, "_UMABackend" if backend == "uma" else "_OrbBackend", fake_backend)
    mlmm_calc._create_ml_backend(
        backend,
        ml_device=torch.device("cpu"),
        analytical_hessian=analytical,
        compile_model=compile_model,
    )
    assert captured["compile_model"] is expected


def test_factory_checks_the_compiler_before_the_model_loads(monkeypatch):
    def no_load(**kwargs):
        pytest.fail("the model must not load when the compiler is missing")

    monkeypatch.delenv("CXX", raising=False)
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from(set()))
    monkeypatch.setattr(mlmm_calc, "_UMABackend", no_load)
    with pytest.raises(pre.MissingPrerequisiteError) as missing:
        mlmm_calc._create_ml_backend("uma", ml_device=torch.device("cpu"), compile_model=True)
    assert missing.value.kind == pre.KIND_CXX_COMPILER


def test_factory_rejects_compile_for_another_backend():
    with pytest.raises(click.UsageError, match="^--compile applies to the uma and orb backends$"):
        mlmm_calc._create_ml_backend("mace", ml_device=torch.device("cpu"), compile_model=True)


@pytest.mark.parametrize("precision", ["fp32", "fp64"])
@pytest.mark.parametrize("compile_model", [False, True])
def test_uma_backend_builds_the_fast_path(monkeypatch, precision, compile_model):
    captured = {}

    class FakeSettings:
        def __init__(self, **kwargs):
            self.kwargs = kwargs
            self.tf32 = kwargs.get("tf32", False)

    class FakePretrained:
        @staticmethod
        def get_predict_unit(_model, **kwargs):
            captured.update(kwargs)
            return object()

    monkeypatch.setattr(mlmm_calc, "HAS_FAIRCHEM", True)
    monkeypatch.setattr(mlmm_calc, "_UMAInferenceSettings", FakeSettings)
    monkeypatch.setattr(mlmm_calc, "pretrained_mlip", FakePretrained(), raising=False)
    monkeypatch.setattr(mlmm_calc, "AtomicData", object, raising=False)
    monkeypatch.setattr(mlmm_calc, "data_list_collater", object(), raising=False)
    mlmm_calc._UMABackend(
        ml_device=torch.device("cpu"), precision=precision, compile_model=compile_model
    )
    expected = {
        "tf32": False,
        "activation_checkpointing": False,
        "merge_mole": True,
        "compile": compile_model,
        "external_graph_gen": False,
        "internal_graph_gen_version": 2,
    }
    if precision == "fp64":
        expected["base_precision_dtype"] = "float64"
    assert captured["inference_settings"].kwargs == expected


class ConservativeForcefieldRegressor:
    pass


def _install_fake_orb(monkeypatch):
    """Just enough of orb_models for _OrbBackend.__init__ to load a model."""
    calls = []

    def loader(**kwargs):
        calls.append(kwargs)
        return ConservativeForcefieldRegressor()

    pretrained = ModuleType("orb_models.forcefield.pretrained")
    pretrained.ORB_PRETRAINED_MODELS = {"orbmol-v2": loader}
    calculator = ModuleType("orb_models.forcefield.inference.calculator")

    class ORBCalculator:
        def __init__(self, model, *, device=None):
            self.model, self.device = model, device

    calculator.ORBCalculator = ORBCalculator
    inference = ModuleType("orb_models.forcefield.inference")
    inference.calculator = calculator
    forcefield = ModuleType("orb_models.forcefield")
    forcefield.pretrained = pretrained
    forcefield.inference = inference
    package = ModuleType("orb_models")
    package.forcefield = forcefield
    for name, module in (
        ("orb_models", package),
        ("orb_models.forcefield", forcefield),
        ("orb_models.forcefield.pretrained", pretrained),
        ("orb_models.forcefield.inference", inference),
        ("orb_models.forcefield.inference.calculator", calculator),
    ):
        monkeypatch.setitem(sys.modules, name, module)
    monkeypatch.setattr(mlmm_calc, "HAS_ORB", True)
    return calls


@pytest.mark.parametrize(
    ("compile_model", "analytical", "expected"),
    [(None, False, False), (False, False, False), (True, False, True), (True, True, False)],
)
def test_orb_receives_an_explicit_bool(monkeypatch, compile_model, analytical, expected):
    calls = _install_fake_orb(monkeypatch)
    monkeypatch.delenv("CXX", raising=False)
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from({"g++"}))
    mlmm_calc._create_ml_backend(
        "orb",
        ml_device=torch.device("cpu"),
        analytical_hessian=analytical,
        compile_model=compile_model,
    )
    assert calls[0]["compile"] is expected


# ---------------------------------------------------------------------------
# all
# ---------------------------------------------------------------------------


def test_all_child_config_carries_only_an_explicit_value(tmp_path):
    from mlmm.workflows.all import _inject_coord_type_into_args_yaml

    assert _inject_coord_type_into_args_yaml(None, None, compile_model=None) is None
    for value in (True, False):
        path = _inject_coord_type_into_args_yaml(None, None, compile_model=value)
        payload = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        assert payload["calc"]["compile_model"] is value
    source = tmp_path / "config.yaml"
    source.write_text("calc:\n  backend: orb\n", encoding="utf-8")
    path = _inject_coord_type_into_args_yaml(source, None, precision="fp32")
    assert "compile_model" not in yaml.safe_load(Path(path).read_text(encoding="utf-8"))["calc"]


@pytest.mark.parametrize(
    ("config", "cli", "message"),
    [
        ("calc:\n  compile_model: true\n", False, "--no-compile and calc.compile_model: true conflict."),
        ("calc:\n  backend: mace\n", True, "--compile applies to the uma and orb backends"),
    ],
)
def test_all_child_config_checks_the_given_yaml(tmp_path, config, cli, message):
    from mlmm.workflows.all import _inject_coord_type_into_args_yaml

    source = tmp_path / "config.yaml"
    source.write_text(config, encoding="utf-8")
    with pytest.raises(click.UsageError) as error:
        _inject_coord_type_into_args_yaml(source, None, compile_model=cli)
    assert error.value.message == message


@needs_toy_system
@pytest.mark.parametrize("flag", ["--compile", "--no-compile"])
def test_all_passes_an_explicit_value_to_the_child_config(tmp_path, monkeypatch, flag):
    from mlmm.workflows import all as all_workflow

    class Reached(BaseException):
        pass

    seen = []

    def recording_inject(*args, **kwargs):
        seen.append(kwargs.get("compile_model"))
        raise Reached

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(all_workflow, "_inject_coord_type_into_args_yaml", recording_inject)
    with pytest.raises(Reached):
        CliRunner().invoke(
            root_cli,
            ["all", "-i", _R, "-i", str(_TOY / "p_toy.pdb"), *_SYSTEM,
             "--out-dir", str(tmp_path / "out"), "--dry-run", flag],
            catch_exceptions=False,
        )
    assert seen == [flag == "--compile"]


@needs_toy_system
def test_all_rejects_compile_with_another_backend_in_dry_run(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    result = CliRunner().invoke(
        root_cli,
        ["all", "-i", _R, "-i", str(_TOY / "p_toy.pdb"), *_SYSTEM, "-b", "mace",
         "--compile", "--out-dir", str(tmp_path / "out"), "--dry-run"],
    )
    assert result.exit_code == 2, result.output
    assert "--compile applies to the uma and orb backends" in result.output


def test_all_summary_records_mlip_compiled(tmp_path):
    from mlmm.workflows import all as all_workflow

    for compiled in (True, False):
        summary = {"segments": []}
        all_workflow._enrich_summary(
            summary, pipeline_mode="path-opt", out_dir=tmp_path,
            mlip_backend="uma", mlip_model="uma-s-1p2p1", mlip_precision="fp32",
            charge=0, spin=1, config={"tsopt": False, "thermo": False, "dft": False},
            calculator_config={"backend": "uma", "uma_precision": "fp32", "compile_model": compiled},
        )
        assert summary["mlip_compiled"] is compiled


# ---------------------------------------------------------------------------
# trj2fig
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(("compile_model", "expected"), [(None, False), (True, True)])
def test_trj2fig_records_mlip_compiled(compile_model, expected):
    from mlmm.io.trj2fig import _resolved_backend_config

    calc_cfg, provenance = _resolved_backend_config(
        "uma", precision="fp64", compile_model=compile_model
    )
    assert calc_cfg["compile_model"] is expected
    assert provenance["mlip_compiled"] is expected
