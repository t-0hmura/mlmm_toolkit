"""Strict SP ML-region selection and provenance contracts."""

from __future__ import annotations

from pathlib import Path

import click
from click.testing import CliRunner
import pytest

from mlmm.core.defaults import MLMM_CALC_KW
from mlmm.core.utils import read_bfactors_from_pdb
from mlmm.workflows import sp


def _write_pdb(path: Path, bfactors) -> Path:
    lines = []
    for serial, bfactor in enumerate(bfactors, start=1):
        lines.append(
            f"ATOM  {serial:5d}  H   MOL A{serial:4d}    "
            f"{float(serial - 1):8.3f}{0.0:8.3f}{0.0:8.3f}"
            f"{1.0:6.2f}{float(bfactor):6.2f}          H \n"
        )
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")
    assert read_bfactors_from_pdb(path) == [float(value) for value in bfactors]
    return path


def _calc_cfg(**updates):
    config = dict(MLMM_CALC_KW)
    config.update(updates)
    return config


def test_mixed_out_of_range_indices_fail_without_full_system_fallback(
    tmp_path: Path,
) -> None:
    source = _write_pdb(tmp_path / "source.pdb", [50.0, 50.0, 50.0])
    config = _calc_cfg(model_pdb=None, use_bfactor_layers=False)

    with pytest.raises(click.BadParameter, match="outside the input atom bounds"):
        sp._resolve_sp_ml_region(
            source_path=source,
            out_dir_path=tmp_path / "out",
            calc_cfg=config,
            model_indices_str="1,4",
            model_indices_one_based=True,
        )

    assert config["model_pdb"] is None
    assert not (tmp_path / "out").exists()


def test_absent_selection_on_unlayered_input_fails_strictly(tmp_path: Path) -> None:
    source = _write_pdb(tmp_path / "source.pdb", [50.0, 50.0])
    with pytest.raises(click.ClickException, match="Invalid or missing layer B-factors"):
        sp._resolve_sp_ml_region(
            source_path=source,
            out_dir_path=tmp_path / "out",
            calc_cfg=_calc_cfg(model_pdb=None, use_bfactor_layers=True),
            model_indices_str=None,
            model_indices_one_based=True,
        )


def test_explicit_all_atom_indices_are_allowed_and_recorded(tmp_path: Path) -> None:
    source = _write_pdb(tmp_path / "source.pdb", [50.0, 50.0, 50.0])
    config = _calc_cfg(model_pdb=None, use_bfactor_layers=True)
    provenance = sp._resolve_sp_ml_region(
        source_path=source,
        out_dir_path=tmp_path / "out",
        calc_cfg=config,
        model_indices_str="1-3",
        model_indices_one_based=True,
    )

    assert provenance["ml_region_source"] == "model_indices"
    assert provenance["ml_region_atom_count"] == 3
    assert provenance["full_system_ml"] is True
    assert provenance["ml_region_indices"] == [0, 1, 2]
    assert Path(provenance["ml_region_model_pdb"]).exists()
    assert config["use_bfactor_layers"] is True


def test_explicit_model_pdb_has_precedence_and_path_provenance(tmp_path: Path) -> None:
    source = _write_pdb(tmp_path / "source.pdb", [50.0, 50.0])
    model = _write_pdb(tmp_path / "model.pdb", [0.0, 0.0])
    config = _calc_cfg(model_pdb=str(model), use_bfactor_layers=True)
    provenance = sp._resolve_sp_ml_region(
        source_path=source,
        out_dir_path=tmp_path / "out",
        calc_cfg=config,
        # The established precedence ignores indices when model_pdb is explicit.
        model_indices_str="999",
        model_indices_one_based=True,
    )

    assert provenance == {
        "ml_region_source": "model_pdb",
        "ml_region_atom_count": 2,
        "full_system_ml": True,
        "ml_region_model_pdb": str(model),
        "ml_region_indices": None,
    }
    assert config["model_pdb"] == str(model)
    assert config["use_bfactor_layers"] is True


def test_valid_bfactor_subset_records_layer_provenance(tmp_path: Path) -> None:
    source = _write_pdb(tmp_path / "source.pdb", [0.0, 10.0, 20.0])
    provenance = sp._resolve_sp_ml_region(
        source_path=source,
        out_dir_path=tmp_path / "out",
        calc_cfg=_calc_cfg(model_pdb=None, use_bfactor_layers=True),
        model_indices_str=None,
        model_indices_one_based=True,
    )
    assert provenance["ml_region_source"] == "bfactor"
    assert provenance["ml_region_atom_count"] == 1
    assert provenance["full_system_ml"] is False
    assert provenance["ml_region_indices"] is None


def test_cli_rejects_invalid_indices_before_calculator_construction(
    tmp_path: Path, monkeypatch
) -> None:
    source = _write_pdb(tmp_path / "source.pdb", [50.0, 50.0])
    parm = tmp_path / "real.parm7"
    parm.write_text("placeholder\n", encoding="utf-8")
    constructed = []

    def _forbidden_constructor(**kwargs):
        constructed.append(kwargs)
        pytest.fail("calculator construction must follow ML-region validation")

    monkeypatch.setattr(sp, "mlmm", _forbidden_constructor)
    result = CliRunner().invoke(
        sp.cli,
        [
            "-i",
            str(source),
            "--parm",
            str(parm),
            "-q",
            "0",
            "--no-detect-layer",
            "--model-indices",
            "1,3",
            "--out-dir",
            str(tmp_path / "result"),
        ],
    )

    assert result.exit_code != 0
    assert "outside the input atom bounds" in result.output
    assert constructed == []
