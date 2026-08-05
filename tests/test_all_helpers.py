"""Unit tests for mlmm.workflows._all_helpers and the `mlmm all` internals in
mlmm.workflows.all."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest


def test_element_fix_paths_do_not_collide_for_same_basename(tmp_path: Path) -> None:
    from mlmm.workflows.all import _element_fix_path

    first = _element_fix_path(tmp_path, Path("reactant/input.pdb"), 1)
    second = _element_fix_path(tmp_path, Path("product/input.pdb"), 2)

    assert first != second
    assert first.name == "001_input.pdb"
    assert second.name == "002_input.pdb"
import yaml

from mlmm.workflows._all_helpers import (
    build_dft_overrides,
    build_energy_level_dict,
    build_freq_overrides,
    build_path_child_argv,
    build_pipeline_summary_payload,
    build_scan_child_argv,
    build_tsopt_overrides,
    build_thermo_symmetry_provenance,
    copy_path_outputs_to_root,
    promote_diag_for_root,
    resolve_dft_func_basis_forwarding,
    resolve_post_thresh_forwarding,
    has_complete_segment_energy_series,
)


def _path_child_kwargs(*, include_opt_mode: bool) -> dict:
    return {
        "include_opt_mode": include_opt_mode,
        "mep_mode": "dmf",
        "dmf_backend": "cpu",
        "max_nodes": 31,
        "max_cycles": 47,
        "climb": False,
        "opt_mode": "hess",
        "dump": False,
        "pre_opt": False,
        "convert_files": False,
        "thresh": "gau_tight",
        "thresh_gsm": "gau",
        "thresh_dmf": "loose",
    }


_PATH_COMMON_CASES = (
    ("dmf_backend", ["--dmf-backend", "cpu"]),
    ("max_nodes", ["--max-nodes", "31"]),
    ("max_cycles", ["--max-cycles", "47"]),
    ("climb", ["--no-climb"]),
    ("dump", ["--no-dump"]),
    ("pre_opt", ["--no-preopt"]),
    ("convert_files", ["--no-convert-files"]),
    ("thresh", ["--thresh", "gau_tight"]),
    ("thresh_gsm", ["--thresh-gsm", "gau"]),
    ("thresh_dmf", ["--thresh-dmf", "loose"]),
)


@pytest.mark.parametrize(
    ("parameter", "expected"),
    _PATH_COMMON_CASES + (("opt_mode", ["--opt-mode", "hess"]),),
)
def test_path_search_child_forwards_each_explicit_field_once(
    parameter: str,
    expected: list[str],
) -> None:
    argv = build_path_child_argv(
        {parameter},
        **_path_child_kwargs(include_opt_mode=True),
    )
    assert argv == ["--mep-mode", "dmf", *expected]
    assert argv.count(expected[0]) == 1


@pytest.mark.parametrize(("parameter", "expected"), _PATH_COMMON_CASES)
def test_path_opt_child_forwards_each_explicit_field_once(
    parameter: str,
    expected: list[str],
) -> None:
    argv = build_path_child_argv(
        {parameter},
        **_path_child_kwargs(include_opt_mode=False),
    )
    assert argv == ["--mep-mode", "dmf", *expected]
    assert argv.count(expected[0]) == 1


def test_path_opt_explicitly_omits_unsupported_opt_mode() -> None:
    assert build_path_child_argv(
        {"opt_mode"},
        **_path_child_kwargs(include_opt_mode=False),
    ) == ["--mep-mode", "dmf"]


@pytest.mark.parametrize("include_opt_mode", [True, False])
def test_path_defaults_leave_pipeline_owned_and_yaml_tokens_unchanged(
    include_opt_mode: bool,
) -> None:
    pipeline_owned = [
        "-i",
        "state.pdb",
        "-q",
        "-1",
        "--parm",
        "real.parm7",
        "--out-dir",
        "result",
        "--config",
        "effective.yaml",
    ]
    argv = pipeline_owned + build_path_child_argv(
        set(),
        **_path_child_kwargs(include_opt_mode=include_opt_mode),
    )
    assert argv == [*pipeline_owned, "--mep-mode", "dmf"]


@pytest.mark.parametrize(
    ("parameter", "expected"),
    [
        ("convert_files", ["--no-convert-files"]),
        ("thresh", ["--thresh", "gau_tight"]),
    ],
)
def test_scan_child_forwards_each_explicit_field_once(
    parameter: str,
    expected: list[str],
) -> None:
    argv = build_scan_child_argv(
        {parameter},
        convert_files=False,
        thresh="gau_tight",
    )
    assert argv == expected
    assert argv.count(expected[0]) == 1


def test_scan_defaults_leave_pipeline_owned_and_yaml_tokens_unchanged() -> None:
    pipeline_owned = [
        "-i",
        "state.pdb",
        "-q",
        "-1",
        "--parm",
        "real.parm7",
        "--out-dir",
        "result",
        "--config",
        "effective.yaml",
    ]
    argv = pipeline_owned + build_scan_child_argv(
        set(),
        convert_files=False,
        thresh="gau_tight",
    )
    assert argv == pipeline_owned


def _tsopt_override_kwargs() -> dict:
    return {
        "tsopt_max_cycles": None,
        "dump": False,
        "dump_override_requested": False,
        "tsopt_out_dir": None,
        "hessian_calc_mode": None,
        "opt_mode_post_norm": "hess",
        "opt_mode_post_set": False,
        "opt_mode_set": False,
        "tsopt_opt_mode_default": "hess",
        "convert_files": True,
        "convert_files_explicit": False,
        "thresh_post_forward": None,
        "flatten_explicit": False,
        "flatten": False,
        "skip_final_freq": False,
        "skip_final_freq_explicit": False,
    }


def _freq_override_kwargs() -> dict:
    return {
        "freq_max_write": None,
        "freq_amplitude_ang": None,
        "freq_n_frames": None,
        "freq_sort": None,
        "freq_temperature": None,
        "freq_pressure": None,
        "dump_override_requested": False,
        "dump": False,
        "require_thermo_artifact": False,
        "hessian_calc_mode": None,
        "convert_files": True,
        "convert_files_explicit": False,
    }


def _dft_override_kwargs() -> dict:
    return {
        "dft_max_cycle": None,
        "dft_conv_tol": None,
        "dft_grid_level": None,
        "dft_engine": None,
        "dft_func_basis_forward": None,
        "convert_files": True,
        "convert_files_explicit": False,
    }


def test_all_post_stage_overrides_do_not_reemit_parent_defaults() -> None:
    assert build_tsopt_overrides(**_tsopt_override_kwargs()) == {}
    assert build_freq_overrides(**_freq_override_kwargs()) == {}
    assert build_dft_overrides(**_dft_override_kwargs()) == {}


@pytest.mark.parametrize(
    ("updates", "expected"),
    [
        ({"tsopt_max_cycles": 41}, {"max_cycles": 41}),
        ({"dump_override_requested": True}, {"dump": False}),
        ({"tsopt_out_dir": Path("ts")}, {"out_dir": Path("ts")}),
        ({"hessian_calc_mode": "Analytical"}, {"hessian_calc_mode": "Analytical"}),
        (
            {"opt_mode_post_norm": "grad", "opt_mode_post_set": True},
            {"opt_mode": "grad"},
        ),
        (
            {"opt_mode_set": True, "tsopt_opt_mode_default": "grad"},
            {"opt_mode": "grad"},
        ),
        (
            {"convert_files": False, "convert_files_explicit": True},
            {"convert_files": False},
        ),
        ({"thresh_post_forward": "gau_loose"}, {"thresh": "gau_loose"}),
        ({"flatten_explicit": True}, {"flatten": False}),
        ({"skip_final_freq_explicit": True}, {"skip_final_freq": False}),
    ],
)
def test_tsopt_override_builder_covers_each_forwarded_field(
    updates: dict,
    expected: dict,
) -> None:
    kwargs = _tsopt_override_kwargs()
    kwargs.update(updates)
    assert build_tsopt_overrides(**kwargs) == expected


@pytest.mark.parametrize(
    ("updates", "expected"),
    [
        ({"freq_max_write": 7}, {"max_write": 7}),
        ({"freq_amplitude_ang": 0.25}, {"amplitude_ang": 0.25}),
        ({"freq_n_frames": 11}, {"n_frames": 11}),
        ({"freq_sort": "ABS"}, {"sort": "abs"}),
        ({"freq_temperature": 310.0}, {"temperature": 310.0}),
        ({"freq_pressure": 0.9}, {"pressure": 0.9}),
        ({"dump_override_requested": True}, {"dump": False}),
        ({"require_thermo_artifact": True}, {"dump": True}),
        (
            {
                "dump_override_requested": True,
                "dump": False,
                "require_thermo_artifact": True,
            },
            {"dump": True},
        ),
        ({"hessian_calc_mode": "Analytical"}, {"hessian_calc_mode": "Analytical"}),
        (
            {"convert_files": False, "convert_files_explicit": True},
            {"convert_files": False},
        ),
    ],
)
def test_freq_override_builder_covers_each_forwarded_field(
    updates: dict,
    expected: dict,
) -> None:
    kwargs = _freq_override_kwargs()
    kwargs.update(updates)
    assert build_freq_overrides(**kwargs) == expected


def test_thermo_symmetry_provenance_copies_every_complete_valid_state() -> None:
    payload = build_thermo_symmetry_provenance(
        {
            "R": {
                "point_group": "C2",
                "point_group_source": "auto",
                "symmetry_number": 2,
                "symmetry_number_source": "auto",
            },
            "TS": {"symmetry_number": 3, "symmetry_number_source": "config"},
            "P": {"symmetry_number": 1},
            "other": {"symmetry_number": 9, "symmetry_number_source": "test"},
        }
    )

    assert payload == {
        "R": {
            "point_group": "C2",
            "point_group_source": "auto",
            "symmetry_number": 2,
            "symmetry_number_source": "auto",
        },
        "TS": {
            "symmetry_number": 3,
            "symmetry_number_source": "config",
        },
        "other": {
            "symmetry_number": 9,
            "symmetry_number_source": "test",
        },
    }


@pytest.mark.parametrize("invalid", [True, 0, -1, 2.0, "2", None])
def test_thermo_symmetry_provenance_rejects_invalid_values(invalid) -> None:
    assert build_thermo_symmetry_provenance(
        {
            "R": {
                "symmetry_number": invalid,
                "symmetry_number_source": "auto",
            }
        }
    ) == {}


def test_thermo_values_are_independent_of_endpoint_mode_counts() -> None:
    from mlmm.workflows.all import _thermo_correction_ha, _thermo_gibbs_ha

    payload = {
        "num_imag_freq": 3,
        "sum_EE_and_thermal_free_energy_ha": -11.0,
        "thermal_correction_free_energy_ha": 0.2,
    }
    assert _thermo_gibbs_ha(payload) == -11.0
    assert _thermo_correction_ha(payload) == 0.2


def test_all_segment_energy_series_requires_every_reactive_segment() -> None:
    one = [(-10.0, -9.0, -11.0)]
    two = [*one, (-11.0, -10.0, -12.0)]
    assert not has_complete_segment_energy_series(one, expected_segments=2)
    assert has_complete_segment_energy_series(two, expected_segments=2)
    assert not has_complete_segment_energy_series(
        [(-10.0, -9.0)], expected_segments=1
    )


@pytest.mark.parametrize(
    ("updates", "expected"),
    [
        ({"dft_max_cycle": 55}, {"max_cycle": 55}),
        ({"dft_conv_tol": 1e-8}, {"conv_tol": 1e-8}),
        ({"dft_grid_level": 5}, {"grid_level": 5}),
        ({"dft_engine": "cpu"}, {"engine": "cpu"}),
        (
            {"dft_func_basis_forward": "pbe0/def2-svp"},
            {"func_basis": "pbe0/def2-svp"},
        ),
        (
            {"convert_files": False, "convert_files_explicit": True},
            {"convert_files": False},
        ),
    ],
)
def test_dft_override_builder_covers_each_forwarded_field(
    updates: dict,
    expected: dict,
) -> None:
    kwargs = _dft_override_kwargs()
    kwargs.update(updates)
    assert build_dft_overrides(**kwargs) == expected


def test_all_dft_child_relays_success_stderr_without_forwarding_mlip_backend(
    tmp_path: Path, monkeypatch
) -> None:
    import subprocess
    from types import SimpleNamespace
    from mlmm.workflows import all as all_workflow

    commands = []
    emitted = []

    def _run(cmd, **kwargs):
        commands.append(cmd)
        return SimpleNamespace(returncode=0, stdout="", stderr="fallback warning")

    monkeypatch.setattr(subprocess, "run", _run)
    monkeypatch.setattr(
        all_workflow,
        "_echo",
        lambda message, err=False: emitted.append((message, err)),
    )

    all_workflow._run_dft_for_state(
        tmp_path / "state.pdb",
        0,
        1,
        tmp_path / "system.parm7",
        tmp_path / "model.pdb",
        False,
        tmp_path / "dft",
        None,
        backend="uma",
    )

    assert "--backend" not in commands[0]
    assert emitted.count(("fallback warning", True)) == 1
    assert not any("exited with code" in message for message, _ in emitted)


def test_yaml_post_values_are_not_reemitted_as_default_cli_tokens() -> None:
    yaml_cfg = {
        "opt": {"thresh": "gau_tight"},
        "dft": {"func_basis": "pbe0/def2-svp"},
    }
    thresh_forward = resolve_post_thresh_forwarding(
        set(),
        thresh_post="baker",
        yaml_cfg=yaml_cfg,
    )
    func_basis_forward, effective_method = resolve_dft_func_basis_forwarding(
        set(),
        dft_func_basis=None,
        yaml_cfg=yaml_cfg,
    )

    tsopt_kwargs = _tsopt_override_kwargs()
    tsopt_kwargs["thresh_post_forward"] = thresh_forward
    dft_kwargs = _dft_override_kwargs()
    dft_kwargs["dft_func_basis_forward"] = func_basis_forward
    assert "thresh" not in build_tsopt_overrides(**tsopt_kwargs)
    assert "func_basis" not in build_dft_overrides(**dft_kwargs)
    assert effective_method == "pbe0/def2-svp"


def test_explicit_post_values_override_yaml_once() -> None:
    yaml_cfg = {
        "opt": {"thresh": "gau_tight"},
        "dft": {"func_basis": "pbe0/def2-svp"},
    }
    thresh_forward = resolve_post_thresh_forwarding(
        {"thresh_post"},
        thresh_post="baker",
        yaml_cfg=yaml_cfg,
    )
    func_basis_forward, effective_method = resolve_dft_func_basis_forwarding(
        {"dft_func_basis"},
        dft_func_basis="wb97x/def2-tzvp",
        yaml_cfg=yaml_cfg,
    )

    tsopt_kwargs = _tsopt_override_kwargs()
    tsopt_kwargs["thresh_post_forward"] = thresh_forward
    dft_kwargs = _dft_override_kwargs()
    dft_kwargs["dft_func_basis_forward"] = func_basis_forward
    assert build_tsopt_overrides(**tsopt_kwargs) == {"thresh": "baker"}
    assert build_dft_overrides(**dft_kwargs) == {
        "func_basis": "wb97x/def2-tzvp"
    }
    assert effective_method == "wb97x/def2-tzvp"


def test_post_threshold_default_is_forwarded_only_without_yaml_value() -> None:
    assert resolve_post_thresh_forwarding(
        set(),
        thresh_post="baker",
        yaml_cfg={},
    ) == "baker"


def test_all_effective_yaml_canonicalizes_alias_only_calculator(tmp_path: Path) -> None:
    from mlmm.workflows.all import _build_effective_args_yaml

    source = tmp_path / "alias.yaml"
    source.write_text(
        yaml.safe_dump(
            {
                "mlmm": {
                    "backend": "orb",
                    "orb_model": "alias-model",
                    "embedcharge": True,
                }
            }
        ),
        encoding="utf-8",
    )

    effective, payload = _build_effective_args_yaml(
        source, None, tmp_prefix="test_mlmm_alias_"
    )
    assert effective is not None
    assert effective != source
    assert payload["calc"] == payload["mlmm"]
    assert payload["calc"]["backend"] == "orb"
    assert yaml.safe_load(effective.read_text(encoding="utf-8"))["calc"] == payload["calc"]


def test_all_effective_yaml_keeps_canonical_whole_section_precedence(
    tmp_path: Path,
) -> None:
    from mlmm.workflows.all import _build_effective_args_yaml

    source = tmp_path / "both.yaml"
    source.write_text(
        yaml.safe_dump(
            {
                "calc": {"backend": "mace", "mace_model": "canonical-model"},
                "mlmm": {"backend": "orb", "orb_model": "legacy-model"},
            }
        ),
        encoding="utf-8",
    )

    effective, payload = _build_effective_args_yaml(
        source, None, tmp_prefix="test_mlmm_both_"
    )
    assert effective == source
    assert payload["calc"] == {
        "backend": "mace",
        "mace_model": "canonical-model",
    }
    assert payload["mlmm"]["backend"] == "orb"


def test_all_injection_preserves_alias_only_calculator_values(tmp_path: Path) -> None:
    from mlmm.workflows.all import _inject_coord_type_into_args_yaml

    source = tmp_path / "alias.yaml"
    source.write_text(
        yaml.safe_dump(
            {
                "mlmm": {
                    "backend": "orb",
                    "orb_model": "alias-model",
                    "embedcharge": True,
                }
            }
        ),
        encoding="utf-8",
    )

    effective = _inject_coord_type_into_args_yaml(
        source,
        None,
        precision="fp64",
        workers=3,
        workers_per_node=2,
        backend_model="explicit-model",
    )
    assert effective is not None
    payload = yaml.safe_load(effective.read_text(encoding="utf-8"))
    assert payload["calc"]["backend"] == "orb"
    assert payload["calc"]["orb_model"] == "explicit-model"
    assert payload["calc"]["embedcharge"] is True
    assert payload["calc"]["orb_precision"] == "float64"
    assert payload["calc"]["workers"] == 3
    assert payload["calc"]["workers_per_node"] == 2


@pytest.mark.parametrize(
    ("backend", "precision_key", "precision_value", "model_key"),
    [
        ("orb", "orb_precision", "float32-high", "orb_model"),
        ("mace", "mace_dtype", "float32", "mace_model"),
    ],
)
def test_all_injection_routes_cli_values_to_explicit_backend(
    backend: str,
    precision_key: str,
    precision_value: str,
    model_key: str,
) -> None:
    """CLI provenance and the child calculator config describe the same run."""
    from mlmm.workflows.all import (
        _inject_coord_type_into_args_yaml,
        _resolve_calculator_template,
        _resolve_mlip_provenance,
    )

    effective = _inject_coord_type_into_args_yaml(
        None,
        None,
        backend=backend,
        precision="fp32",
        backend_model=f"{backend}-custom",
    )
    assert effective is not None
    payload = yaml.safe_load(effective.read_text(encoding="utf-8"))
    calc = payload["calc"]
    assert calc["backend"] == backend
    assert calc[precision_key] == precision_value
    assert calc[model_key] == f"{backend}-custom"
    assert "uma_precision" not in calc
    assert "uma_model" not in calc

    template = _resolve_calculator_template(
        effective,
        backend=None,
        embedcharge=False,
        embedcharge_explicit=False,
        embedcharge_cutoff=None,
        link_atom_method=None,
        mm_backend=None,
        use_cmap=None,
    ).materialize()
    provenance = _resolve_mlip_provenance(
        backend=None,
        backend_model=None,
        calc_file=None,
        calc_factory=None,
        precision=None,
        merged_yaml_cfg=payload,
    )
    assert template["backend"] == provenance[0] == backend
    assert template[model_key] == provenance[1] == f"{backend}-custom"
    assert provenance[2] == "fp32"


def test_all_injection_enforces_aimnet_precision_contract() -> None:
    from mlmm.workflows.all import _inject_coord_type_into_args_yaml

    with pytest.raises(ValueError, match="fp64 is not supported by backend 'aimnet2'"):
        _inject_coord_type_into_args_yaml(
            None, None, backend="aimnet2", precision="fp64",
        )

    effective = _inject_coord_type_into_args_yaml(
        None,
        None,
        backend="aimnet2",
        precision="fp32",
        backend_model="aimnet-custom",
    )
    assert effective is not None
    calc = yaml.safe_load(effective.read_text(encoding="utf-8"))["calc"]
    assert calc == {
        "backend": "aimnet2",
        "aimnet2_model": "aimnet-custom",
    }


def test_all_injection_translates_yaml_only_generic_backend_aliases(
    tmp_path: Path,
) -> None:
    from mlmm.workflows.all import (
        _inject_coord_type_into_args_yaml,
        _resolve_calculator_template,
        _resolve_mlip_provenance,
    )

    source = tmp_path / "generic-orb.yaml"
    source.write_text(yaml.safe_dump({
        "calc": {
            "backend": "orb",
            "precision": "fp32",
            "backend_model": "orb-custom",
        }
    }))
    effective = _inject_coord_type_into_args_yaml(source, None)
    assert effective is not None and effective != source
    payload = yaml.safe_load(effective.read_text(encoding="utf-8"))
    assert payload["calc"] == {
        "backend": "orb",
        "orb_precision": "float32-high",
        "orb_model": "orb-custom",
    }
    template = _resolve_calculator_template(
        effective,
        backend=None,
        embedcharge=False,
        embedcharge_explicit=False,
        embedcharge_cutoff=None,
        link_atom_method=None,
        mm_backend=None,
        use_cmap=None,
    ).materialize()
    provenance = _resolve_mlip_provenance(
        backend=None,
        backend_model=None,
        calc_file=None,
        calc_factory=None,
        precision=None,
        merged_yaml_cfg=payload,
    )
    assert template["orb_model"] == provenance[1] == "orb-custom"
    assert template["orb_precision"] == "float32-high"
    assert provenance == ("orb", "orb-custom", "fp32")


def test_all_injection_translates_yaml_only_custom_calculator(
    tmp_path: Path,
) -> None:
    from mlmm.workflows.all import (
        _inject_coord_type_into_args_yaml,
        _resolve_calculator_template,
        _resolve_mlip_provenance,
    )

    source = tmp_path / "custom.yaml"
    source.write_text(yaml.safe_dump({
        "calc": {
            "backend": "orb",
            "calc_file": "custom_calc.py",
            "calc_factory": "build_calc",
        }
    }))
    effective = _inject_coord_type_into_args_yaml(source, None)
    assert effective is not None and effective != source
    payload = yaml.safe_load(effective.read_text(encoding="utf-8"))
    assert payload["calc"]["backend"] == "custom"
    assert payload["calc"]["calc_file"] == "custom_calc.py"
    assert payload["calc"]["calc_factory"] == "build_calc"
    template = _resolve_calculator_template(
        effective,
        backend=None,
        embedcharge=False,
        embedcharge_explicit=False,
        embedcharge_cutoff=None,
        link_atom_method=None,
        mm_backend=None,
        use_cmap=None,
    ).materialize()
    provenance = _resolve_mlip_provenance(
        backend=None,
        backend_model=None,
        calc_file=None,
        calc_factory=None,
        precision=None,
        merged_yaml_cfg=payload,
    )
    assert template["backend"] == provenance[0] == "custom"
    assert template["calc_file"] == "custom_calc.py"
    assert provenance[1] == "custom_calc.py:build_calc"


def test_all_calc_file_injection_preserves_alias_only_calculator_values(
    tmp_path: Path,
) -> None:
    from mlmm.workflows.all import _inject_coord_type_into_args_yaml

    source = tmp_path / "alias.yaml"
    source.write_text(
        yaml.safe_dump(
            {
                "mlmm": {
                    "backend": "orb",
                    "embedcharge": True,
                    "embedcharge_cutoff": 8.5,
                }
            }
        ),
        encoding="utf-8",
    )

    effective = _inject_coord_type_into_args_yaml(
        source,
        None,
        calc_file="custom_calc.py",
        calc_factory="make_calculator",
    )
    assert effective is not None
    payload = yaml.safe_load(effective.read_text(encoding="utf-8"))
    assert payload["calc"]["backend"] == "custom"
    assert payload["calc"]["calc_file"] == "custom_calc.py"
    assert payload["calc"]["calc_factory"] == "make_calculator"
    assert payload["calc"]["embedcharge"] is True
    assert payload["calc"]["embedcharge_cutoff"] == 8.5


def test_build_energy_level_dict_zero_referenced_kcal() -> None:
    # AU2KCALPERMOL is approximately 627.5095
    d = build_energy_level_dict(
        labels=["R", "TS", "P"],
        energies_au=[-100.0, -99.5, -100.2],
        ref_energy=-100.0,
        au_to_kcal=627.5095,
        diagram_path="/tmp/diag.png",
        structures={"R": "r.pdb", "TS": "ts.pdb", "P": "p.pdb"},
    )
    assert d["labels"] == ["R", "TS", "P"]
    assert d["energies_au"] == [-100.0, -99.5, -100.2]
    # 0, 0.5*627.5, -0.2*627.5
    assert d["energies_kcal"][0] == 0.0
    assert abs(d["energies_kcal"][1] - 313.75475) < 1e-3
    assert abs(d["energies_kcal"][2] - (-125.5019)) < 1e-3
    assert d["barrier_kcal"] == d["energies_kcal"][1]
    assert d["delta_kcal"] == d["energies_kcal"][-1]
    assert d["diagram"] == "/tmp/diag.png"
    assert d["structures"] == {"R": "r.pdb", "TS": "ts.pdb", "P": "p.pdb"}


def test_build_energy_level_dict_keeps_unassigned_endpoints_directional() -> None:
    d = build_energy_level_dict(
        labels=["E1", "TS", "E2"],
        energies_au=[-100.0, -99.5, -100.2],
        ref_energy=-100.0,
        au_to_kcal=627.5095,
        diagram_path="/tmp/diag.png",
        structures={"E1": "e1.pdb", "TS": "ts.pdb", "E2": "e2.pdb"},
    )
    assert "barrier_kcal" not in d
    assert "delta_kcal" not in d
    assert d["barrier_from_endpoint_1_kcal"] == d["energies_kcal"][1]
    assert d["barrier_from_endpoint_2_kcal"] == (
        d["energies_kcal"][1] - d["energies_kcal"][2]
    )


def test_build_energy_level_dict_does_not_mutate_inputs() -> None:
    labels = ["R", "TS", "P"]
    energies = [-1.0, -0.5, -1.1]
    structs = {"R": "r.pdb"}
    d = build_energy_level_dict(
        labels=labels,
        energies_au=energies,
        ref_energy=-1.0,
        au_to_kcal=627.5,
        diagram_path="/tmp/x.png",
        structures=structs,
    )
    # caller's containers must not be aliased
    d["labels"].append("EXTRA")
    d["structures"]["NEW"] = "n.pdb"
    assert labels == ["R", "TS", "P"]
    assert "NEW" not in structs


def test_promote_diag_for_root_returns_none_on_empty() -> None:
    assert promote_diag_for_root(None, "stem", Path(".")) is None
    assert promote_diag_for_root({}, "stem", Path(".")) is None


def test_promote_diag_for_root_rewrites_name_and_image() -> None:
    original = {"name": "MEP", "image": "old.png", "x": [0, 1]}
    promoted = promote_diag_for_root(original, "energy_diagram_MLIP", Path("/tmp/out"))
    assert promoted is not None
    assert promoted["name"] == "energy_diagram_MLIP_all"
    assert promoted["image"] == "/tmp/out/energy_diagram_MLIP_all.png"
    # caller's dict must not be mutated
    assert original["name"] == "MEP"
    assert original["image"] == "old.png"
    # other keys pass through
    assert promoted["x"] == [0, 1]


def test_copy_path_outputs_skips_missing_files(tmp_path: Path) -> None:
    src = tmp_path / "src"
    dst = tmp_path / "dst"
    src.mkdir()
    dst.mkdir()
    # No files present — best-effort no-op
    warnings: list[str] = []
    copy_path_outputs_to_root(src, dst, warn_fn=warnings.append)
    assert warnings == []
    assert list(dst.iterdir()) == []


def test_copy_path_outputs_copies_known_artefacts(tmp_path: Path) -> None:
    src = tmp_path / "src"
    dst = tmp_path / "dst"
    src.mkdir()
    dst.mkdir()
    (src / "mep_plot.png").write_text("dummy-png")
    (src / "mep.pdb").write_text("dummy-pdb")
    (src / "mep.cif").write_text("dummy-cif")
    (src / "summary.json").write_text("{}")
    (src / "mep_trj.xyz").write_text("xyz")
    (src / "unrelated.tmp").write_text("ignore me")

    copy_path_outputs_to_root(src, dst)

    assert (dst / "mep_plot.png").read_text() == "dummy-png"
    assert (dst / "mep.pdb").read_text() == "dummy-pdb"
    assert (dst / "mep.cif").read_text() == "dummy-cif"
    assert (dst / "summary.json").read_text() == "{}"
    assert (dst / "mep_trj.xyz").read_text() == "xyz"
    assert not (dst / "unrelated.tmp").exists()


def test_build_pipeline_summary_payload_shape() -> None:
    with tempfile.TemporaryDirectory() as d:
        out_dir = Path(d) / "out"
        path_dir = Path(d) / "path"
        out_dir.mkdir()
        path_dir.mkdir()
        summary = {
            "n_images": 5,
            "n_segments": 1,
            "segments": [{"kind": "seg", "bond_changes": "A->B"}],
            "energy_diagrams": [{"name": "MEP", "x": [0, 1]}],
            "status": "partial",
            "status_reasons": ["legacy incomplete"],
            "execution_status": "completed",
            "scientific_status": "failed",
            "scientific_status_reasons": ["endpoint optimization failed"],
        }
        payload = build_pipeline_summary_payload(
            out_dir=out_dir,
            path_dir=path_dir,
            summary=summary,
            refine_path=True,
            thresh="gau_loose",
            thresh_post="gau",
            flatten=False,
            do_tsopt=True,
            do_thermo=False,
            do_dft=False,
            opt_mode_norm="grad",
            opt_mode_post="HESS",
            path_opt_mode="grad",
            post_opt_mode="HESS",
            ts_opt_mode="HESS",
            endpoint_opt_mode="GRAD",
            mep_mode="dmf",
            dmf_backend="cpu",
            dmf_correlated=True,
            command_str="mlmm all -i foo.pdb",
            q_int=-1,
            spin=1,
            post_segment_logs=[{"seg": 1, "status": "ok"}],
        )
    assert payload["pipeline_mode"] == "path-search"
    assert payload["refine_path"] is True
    assert payload["opt_mode"] == "grad"
    assert payload["opt_mode_post"] == "hess"
    assert payload["path_opt_mode"] == "grad"
    assert payload["post_opt_mode"] == "hess"
    assert payload["ts_opt_mode"] == "hess"
    assert payload["endpoint_opt_mode"] == "grad"
    assert payload["mep_mode"] == "dmf"
    assert payload["dmf_backend"] == "cpu"
    assert payload["dmf_correlated"] is True
    assert payload["charge"] == -1
    assert payload["spin"] == 1
    assert payload["mlip_backend"] == "uma"
    assert payload["mlip_model"] is None
    assert payload["mlip_precision"] is None
    assert payload["status"] == "partial"
    assert payload["status_reasons"] == ["legacy incomplete"]
    assert payload["execution_status"] == "completed"
    assert payload["scientific_status"] == "failed"
    assert payload["scientific_status_reasons"] == [
        "endpoint optimization failed"
    ]
    assert payload["mep"]["n_images"] == 5
    assert payload["mep"]["diagram"]["name"] == "MEP"
    assert payload["post_segments"] == [{"seg": 1, "status": "ok"}]
    assert payload["energy_diagrams"] == summary["energy_diagrams"]
def test_irc_endpoint_topology_tie_uses_rmsd_and_records_provenance(
    monkeypatch,
):
    from types import SimpleNamespace

    import numpy as np

    from mlmm.workflows import all as workflow

    left = SimpleNamespace(coords=np.array([10.0, 0.0, 0.0]))
    right = SimpleNamespace(coords=np.array([0.0, 0.0, 0.0]))
    mep_left = SimpleNamespace(coords=np.array([0.0, 0.0, 0.0]))
    mep_right = SimpleNamespace(coords=np.array([10.0, 0.0, 0.0]))
    monkeypatch.setattr(
        workflow._path_search,
        "_has_bond_change",
        lambda *_args, **_kwargs: (False, ""),
    )

    (
        oriented_left,
        oriented_right,
        _left_tag,
        _right_tag,
        reversed_irc,
        assignment,
    ) = workflow._orient_irc_endpoint_geometries(
        left,
        right,
        mep_left,
        mep_right,
    )

    assert (oriented_left, oriented_right) == (right, left)
    assert reversed_irc is True
    assert assignment["method"] == "rmsd_topology_tie"
    assert assignment["rmsd_swapped"] < assignment["rmsd_direct"]
    assert assignment["connectivity_validated"] is True
