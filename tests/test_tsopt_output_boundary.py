"""TSOPT output ownership and mode-provenance regressions."""

from pathlib import Path

import click
import pytest

from mlmm.workflows.tsopt import (
    _heavy_mode_label,
    _post_analysis_hessian_config,
    _prepare_tsopt_output_dir,
)


@pytest.mark.parametrize(
    ("mode", "label"),
    [
        ("rsirfo", "RS-I-RFO"),
        ("rsprfo", "RS-P-RFO"),
        ("trim", "TRIM"),
    ],
)
def test_heavy_mode_labels_are_specific(mode: str, label: str) -> None:
    assert _heavy_mode_label(mode) == label


def test_full_hessian_flatten_changes_post_analysis_storage_policy() -> None:
    cfg = {"return_partial_hessian": True, "hess_cutoff": 8.0}

    resolved = _post_analysis_hessian_config(cfg, partial=False)

    assert resolved["return_partial_hessian"] is False
    assert resolved["hess_cutoff"] == 8.0
    assert cfg["return_partial_hessian"] is True


def test_prepare_tsopt_output_removes_only_command_owned_artifacts(
    tmp_path: Path,
) -> None:
    vib = tmp_path / "vib"
    vib.mkdir()
    stale = [
        tmp_path / "final_geometry.xyz",
        tmp_path / "result.json",
        vib / "imaginary_mode_1.xyz",
        vib / "imaginary_mode_1.pdb",
    ]
    for path in stale:
        path.write_text("stale", encoding="utf-8")
    retained = tmp_path / "notes.txt"
    retained.write_text("keep", encoding="utf-8")

    _prepare_tsopt_output_dir(tmp_path)

    assert not any(path.exists() for path in stale)
    assert retained.read_text(encoding="utf-8") == "keep"


def test_prepare_tsopt_output_rejects_input_collision(tmp_path: Path) -> None:
    source = tmp_path / "final_geometry.xyz"
    source.write_text("input", encoding="utf-8")

    with pytest.raises(click.UsageError, match="collides"):
        _prepare_tsopt_output_dir(tmp_path, protected_inputs=(source,))

    assert source.read_text(encoding="utf-8") == "input"
