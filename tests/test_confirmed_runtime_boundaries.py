"""Focused regressions for confirmed runtime boundary defects."""

from __future__ import annotations

from pathlib import Path

import click
import pytest


def test_opt_collision_helper_detects_hardlink(tmp_path: Path) -> None:
    from mlmm.workflows.opt import _reject_opt_output_collisions

    source = tmp_path / "input.pdb"
    source.write_text("input", encoding="utf-8")
    destination = tmp_path / "final_geometry.pdb"
    destination.hardlink_to(source)

    with pytest.raises(click.UsageError, match="aliases consumed input"):
        _reject_opt_output_collisions(tmp_path, (source,))


def test_opt_optional_invalidation_preserves_unrelated_files(tmp_path: Path) -> None:
    from mlmm.workflows.opt import _invalidate_opt_optional_outputs

    stale = tmp_path / "result.json"
    unrelated = tmp_path / "notes.txt"
    stale.write_text("old", encoding="utf-8")
    unrelated.write_text("keep", encoding="utf-8")

    _invalidate_opt_optional_outputs(tmp_path)

    assert not stale.exists()
    assert unrelated.read_text(encoding="utf-8") == "keep"


def test_opt_resolves_model_pdb_before_collision_validation() -> None:
    from mlmm.workflows import opt

    source = Path(opt.__file__).read_text(encoding="utf-8")
    cli_source = source[source.index("def cli("):]
    assignment = cli_source.index('model_pdb_cfg = calc_cfg.get("model_pdb")')
    collision_check = cli_source.index("        _reject_opt_output_collisions(")

    assert assignment < collision_check


def test_path_search_collision_helper_detects_symlink(tmp_path: Path) -> None:
    from mlmm.workflows.path_search import _reject_path_search_output_collisions

    source = tmp_path / "input.xyz"
    source.write_text("input", encoding="utf-8")
    (tmp_path / "mep_trj.xyz").symlink_to(source)

    with pytest.raises(click.UsageError, match="aliases consumed input"):
        _reject_path_search_output_collisions(tmp_path, (source,))


def test_global_segment_labels_preserve_noncontiguous_ids() -> None:
    from mlmm.workflows.all import _build_global_segment_labels

    assert _build_global_segment_labels([2, 4]) == [
        "R",
        "TS2",
        "IM2_1",
        "IM2_2",
        "TS4",
        "P",
    ]
