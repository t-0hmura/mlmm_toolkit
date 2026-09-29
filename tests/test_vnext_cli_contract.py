from pathlib import Path

import click
from click.core import ParameterSource
from click.testing import CliRunner
import pytest

from mlmm.cli.app import cli as root_cli
from mlmm.cli.decorators import (
    canonicalize_calculator_section,
    resolve_model_indices_setting,
)
from mlmm.core.utils import apply_yaml_overrides


def _option(command: str, name: str) -> click.Option:
    cmd = root_cli.get_command(click.Context(root_cli), command)
    assert cmd is not None
    return next(param for param in cmd.params if param.name == name)


class _Context:
    def __init__(self, **sources):
        self.sources = sources

    def get_parameter_source(self, name):
        return self.sources.get(name, ParameterSource.DEFAULT)


def test_yaml_model_indices_base_and_fixed_cli_base() -> None:
    spec, one_based = resolve_model_indices_setting(
        _Context(), {"calc": {"model_indices": [0, 2], "model_indices_base": 0}},
        None, True,
    )
    assert (spec, one_based) == ("0,2", False)

    spec, one_based = resolve_model_indices_setting(
        _Context(model_indices_str=ParameterSource.COMMANDLINE),
        {"calc": {"model_indices": [0, 2], "model_indices_base": 0}},
        "1,3", False,
    )
    assert (spec, one_based) == ("1,3", True)

    spec, one_based = resolve_model_indices_setting(
        _Context(
            model_indices_str=ParameterSource.COMMANDLINE,
            model_indices_one_based=ParameterSource.COMMANDLINE,
        ), {}, "0,2", False,
    )
    assert (spec, one_based) == ("0,2", False)


def test_calc_and_legacy_section_conflicts_are_rejected() -> None:
    merged = canonicalize_calculator_section(
        {"calc": {"backend": "uma"}, "mlmm": {"workers": 2}}
    )
    assert merged["calc"] == {"backend": "uma", "workers": 2}
    with pytest.raises(click.BadParameter, match="Conflicting YAML values"):
        canonicalize_calculator_section(
            {"calc": {"backend": "uma"}, "mlmm": {"backend": "orb"}}
        )


def test_calc_and_mlmm_sections_merge_in_runtime_override_path() -> None:
    target = {}
    apply_yaml_overrides(
        {"calc": {"backend": "uma"}, "mlmm": {"workers": 2}},
        [(target, (("calc",), ("mlmm",)))],
    )
    assert target == {"backend": "uma", "workers": 2}
    with pytest.raises(click.BadParameter, match="Conflicting YAML values"):
        apply_yaml_overrides(
            {"calc": {"backend": "uma"}, "mlmm": {"backend": "orb"}},
            [(target, (("calc",), ("mlmm",)))],
        )


def test_canonical_names_keep_compatibility_aliases() -> None:
    assert _option("all", "workers").opts == ["--uma-workers", "--workers"]
    assert _option("all", "parm7_override").opts == ["--parm7", "--parm"]
    assert _option("all", "dft_func_basis").opts == ["--func-basis", "--dft-func-basis"]
    assert _option("tsopt", "hess_cutoff").opts == ["--hessian-cutoff", "--radius-hessian", "--hess-cutoff"]
    assert _option("tsopt", "model_indices_one_based").hidden is True


def test_auto_mm_toggles_use_no_form_and_keep_hidden_old_names(
    tmp_path: Path, monkeypatch
) -> None:
    import mlmm.workflows.all as all_workflow

    assert _option("all", "mm_add_ter").secondary_opts == ["--no-auto-mm-add-ter"]
    assert _option("all", "mm_auto_disulfide").secondary_opts == ["--no-auto-mm-disulfide"]
    for name, legacy in (
        ("mm_no_add_ter", "--auto-mm-no-add-ter"),
        ("mm_no_disulfide", "--auto-mm-no-disulfide"),
    ):
        assert _option("all", name).opts == [legacy]
        assert _option("all", name).hidden is True

    smoke = Path(__file__).resolve().parent / "smoke"
    inputs = [smoke / "r_complex_layered.pdb", smoke / "p_complex_layered.pdb"]
    if not all(path.is_file() for path in inputs):
        pytest.skip("smoke inputs are not present")
    captured = {}

    class _StopAtMmParm(Exception):
        pass

    def fake_build_mm_parm7(**kwargs):
        captured.update(kwargs)
        raise _StopAtMmParm

    monkeypatch.setattr(all_workflow, "_missing_ambertools_commands", lambda _paths: [])
    monkeypatch.setattr(all_workflow, "_build_mm_parm7", fake_build_mm_parm7)
    base = [
        "all", "-i", str(inputs[0]), str(inputs[1]),
        "-q", "-1", "-m", "1", "--detect-layer", "--out-dir", str(tmp_path / "out"),
    ]
    CliRunner().invoke(root_cli, [*base, "--auto-mm-no-add-ter", "--auto-mm-no-disulfide"])
    assert captured["add_ter"] is False
    assert captured["auto_disulfide"] is False

    result = CliRunner().invoke(root_cli, [*base, "--auto-mm-add-ter", "--auto-mm-no-add-ter"])
    assert result.exit_code == 2
    assert "Conflicting values were supplied through aliases" in result.output


def test_conflicting_alias_values_are_rejected_before_execution() -> None:
    result = CliRunner().invoke(
        root_cli,
        ["dft", "--scf-tol", "1e-8", "--conv-tol", "1e-7"],
    )
    assert result.exit_code == 2
    assert "Conflicting values were supplied through aliases" in result.output
