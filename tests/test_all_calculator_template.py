"""Resolved-calculator identity contracts for the ``all`` orchestrator."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from mlmm.workflows.all import (
    _inject_coord_type_into_args_yaml,
    _resolve_calculator_template,
    _resolve_mlip_provenance,
    _stage_calc_kwargs,
)
from mlmm.workflows._all_helpers import append_backend_forwarding_args


def _resolve(path: Path, **overrides):
    options = {
        "backend": None,
        "embedcharge": False,
        "embedcharge_explicit": False,
        "embedcharge_cutoff": None,
        "link_atom_method": None,
        "mm_backend": None,
        "use_cmap": None,
    }
    options.update(overrides)
    return _resolve_calculator_template(path, **options)


@pytest.mark.parametrize(
    ("backend", "model_key", "model_value", "precision_key", "precision_value"),
    [
        ("orb", "orb_model", "orb-sentinel", "orb_precision", "float32-highest"),
        ("mace", "mace_model", "mace-sentinel", "mace_dtype", "float32"),
        ("uma", "uma_model", "uma-sentinel", "uma_precision", "fp64"),
    ],
)
def test_stage_template_preserves_backend_specific_and_mm_fields(
    tmp_path: Path,
    backend: str,
    model_key: str,
    model_value: str,
    precision_key: str,
    precision_value: str,
) -> None:
    config = tmp_path / f"{backend}.yaml"
    config.write_text(
        yaml.safe_dump(
            {
                "calc": {
                    "backend": backend,
                    model_key: model_value,
                    precision_key: precision_value,
                    "workers": 7,
                    "workers_per_node": 3,
                    "embedcharge": False,
                    "mm_backend": "openmm",
                    "mm_threads": 5,
                    "mm_fd_delta": 0.004,
                }
            }
        ),
        encoding="utf-8",
    )

    template = _resolve(config)
    first = _stage_calc_kwargs(
        template,
        input_pdb="first.pdb",
        real_parm7="real.parm7",
        model_pdb="model.pdb",
        charge=-1,
        spin=2,
        use_bfactor_layers=False,
    )
    second = _stage_calc_kwargs(
        template,
        input_pdb="second.pdb",
        real_parm7="other.parm7",
        model_pdb="other-model.pdb",
        charge=1,
        spin=1,
        use_bfactor_layers=True,
    )

    for derived in (first, second):
        assert derived["backend"] == backend
        assert derived[model_key] == model_value
        assert derived[precision_key] == precision_value
        assert derived["workers"] == 7
        assert derived["workers_per_node"] == 3
        assert derived["embedcharge"] is False
        assert derived["mm_backend"] == "openmm"
        assert derived["mm_threads"] == 5
        assert derived["mm_fd_delta"] == pytest.approx(0.004)

    assert first["input_pdb"] == "first.pdb"
    assert second["input_pdb"] == "second.pdb"
    first["freeze_atoms"].append(99)
    assert 99 not in second["freeze_atoms"]
    assert 99 not in template.materialize()["freeze_atoms"]


def test_custom_calculator_and_explicit_overlays_survive_stage_derivation(
    tmp_path: Path,
) -> None:
    config = tmp_path / "custom.yaml"
    config.write_text(
        yaml.safe_dump(
            {
                "mlmm": {
                    "backend": "custom",
                    "calc_file": "custom_calc.py",
                    "calc_factory": "build_calc",
                    "embedcharge": True,
                    "link_atom_method": "scaled",
                    "use_cmap": False,
                }
            }
        ),
        encoding="utf-8",
    )

    template = _resolve(
        config,
        embedcharge=False,
        embedcharge_explicit=True,
        link_atom_method="fixed",
        mm_backend="openmm",
        use_cmap=True,
    )
    derived = _stage_calc_kwargs(
        template,
        input_pdb="state.pdb",
        real_parm7="real.parm7",
        model_pdb="model.pdb",
        charge=0,
        spin=1,
        use_bfactor_layers=True,
    )

    assert derived["backend"] == "custom"
    assert derived["calc_file"] == "custom_calc.py"
    assert derived["calc_factory"] == "build_calc"
    assert derived["embedcharge"] is False
    assert derived["link_atom_method"] == "fixed"
    assert derived["mm_backend"] == "openmm"
    assert derived["use_cmap"] is True


def test_yaml_embedding_and_cli_cutoff_forward_as_one_effective_state(
    tmp_path: Path,
) -> None:
    config = tmp_path / "embedding.yaml"
    config.write_text(
        yaml.safe_dump(
            {"calc": {"embedcharge": True, "embedcharge_cutoff": 12.0}}
        ),
        encoding="utf-8",
    )
    effective = _resolve(config, embedcharge_cutoff=8.5).materialize()
    argv: list[str] = []

    append_backend_forwarding_args(
        argv,
        backend=None,
        embedcharge=bool(effective["embedcharge"]),
        embedcharge_cutoff=float(effective["embedcharge_cutoff"]),
        embedcharge_explicit=False,
        link_atom_method=None,
        mm_backend=None,
        use_cmap=None,
    )

    assert argv == ["--embedcharge", "--embedcharge-cutoff", "8.5"]


@pytest.mark.parametrize("source_kind", ["config", "cli"])
def test_custom_calculator_outranks_explicit_backend_for_every_stage(
    tmp_path: Path,
    source_kind: str,
) -> None:
    """Child and in-process stages must resolve the same custom calculator."""
    calc_file = tmp_path / "custom_calc.py"
    calc_file.write_text("def build_calc():\n    return None\n", encoding="utf-8")
    source = None
    injected_calc_file = str(calc_file) if source_kind == "cli" else None
    if source_kind == "config":
        source = tmp_path / "config.yaml"
        source.write_text(
            yaml.safe_dump({
                "calc": {
                    "backend": "uma",
                    "calc_file": str(calc_file),
                    "calc_factory": "build_calc",
                }
            }),
            encoding="utf-8",
        )

    effective = _inject_coord_type_into_args_yaml(
        source,
        None,
        backend="orb",
        calc_file=injected_calc_file,
        calc_factory="build_calc",
    )
    assert effective is not None
    payload = yaml.safe_load(effective.read_text(encoding="utf-8"))
    template = _resolve(
        effective,
        backend="orb",
    ).materialize()
    provenance = _resolve_mlip_provenance(
        backend="orb",
        backend_model=None,
        calc_file=injected_calc_file,
        calc_factory="build_calc",
        merged_yaml_cfg=payload,
    )

    assert payload["calc"]["backend"] == "custom"
    assert template["backend"] == provenance[0] == "custom"
    assert template["calc_file"] == str(calc_file)
    assert provenance[1] == f"{calc_file.name}:build_calc"


def test_template_drops_config_only_aliases_and_is_read_only(tmp_path: Path) -> None:
    config = tmp_path / "aliases.yaml"
    config.write_text(
        yaml.safe_dump(
            {
                "calc": {
                    "backend": "orb",
                    "charge": 8,
                    "spin": 9,
                    "precision": "fp32",
                    "backend_model": "legacy-model",
                }
            }
        ),
        encoding="utf-8",
    )

    template = _resolve(config)
    with pytest.raises(TypeError):
        template.values["backend"] = "mace"  # type: ignore[index]

    derived = _stage_calc_kwargs(
        template,
        input_pdb="state.pdb",
        real_parm7="real.parm7",
        model_pdb="model.pdb",
        charge=-2,
        spin=3,
        use_bfactor_layers=False,
    )
    assert "charge" not in derived
    assert "spin" not in derived
    assert "precision" not in derived
    assert "backend_model" not in derived
    assert derived["model_charge"] == -2
    assert derived["model_mult"] == 3
