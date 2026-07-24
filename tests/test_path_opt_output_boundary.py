"""Path-optimization output-boundary regressions."""

from __future__ import annotations

from pathlib import Path


def test_prepare_path_output_dir_invalidates_prior_envelopes(
    tmp_path: Path,
) -> None:
    from mlmm.workflows.path_opt import _prepare_path_output_dir

    out_dir = tmp_path / "path"
    out_dir.mkdir()
    for name in ("result.json", "summary.json"):
        (out_dir / name).write_text('{"status": "error"}\n', encoding="utf-8")
    unrelated = out_dir / "notes.txt"
    unrelated.write_text("keep\n", encoding="utf-8")

    resolved = _prepare_path_output_dir(out_dir)

    assert resolved == out_dir.resolve()
    assert not (out_dir / "result.json").exists()
    assert not (out_dir / "summary.json").exists()
    assert unrelated.read_text(encoding="utf-8") == "keep\n"
