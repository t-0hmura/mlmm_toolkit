"""``all --resume-segment`` continues with the parm7 the completed MEP saved.

A rebuilt parm7 need not reproduce its AM1-BCC charges bit for bit, so the
resume no longer rebuilds and compares it. The resume still refuses a changed
input, a different ``--parm`` and a saved parm7 that no longer matches the
recorded hash.
"""

from __future__ import annotations

from pathlib import Path

from mlmm.workflows._segment_resume import parm7_content_sha256
from mlmm.workflows import all as all_workflow

from _all_stub_harness import option, run_all, write

RESUME = ("--tsopt", "--resume-segment", "1")


def _stub_mm_parm(monkeypatch, builds):
    """Every build writes other charges, as sqm under another load could."""

    def build(*, pdb, out_dir, name=None, **_kwargs):
        builds.append(Path(pdb))
        out_dir.mkdir(parents=True, exist_ok=True)
        stem = name or Path(pdb).stem
        parm7, rst7 = out_dir / f"{stem}.parm7", out_dir / f"{stem}.rst7"
        parm7.write_text(
            f"%VERSION  build {len(builds)}\n%FLAG CHARGE\n{0.1 * len(builds):.4f}\n",
            encoding="utf-8",
        )
        rst7.write_text("fixture rst7\n", encoding="utf-8")
        return parm7.resolve(), rst7.resolve()

    monkeypatch.setattr(all_workflow, "_build_mm_parm7", build)
    monkeypatch.setattr(all_workflow, "_missing_ambertools_commands", lambda *_a, **_k: [])


def _completed_mep(tmp_path, monkeypatch, builds):
    first = run_all(tmp_path, monkeypatch, extra=("--tsopt",), parm=False)
    assert first.result.exit_code == 0, first.result.output[-3000:]
    saved = (first.out / "mm_parm" / "input0.parm7").resolve()
    recorded = first.summary["resume_identity"]["path_request"]["parm7_sha256"]
    assert len(builds) == 1 and parm7_content_sha256(saved) == recorded
    return first, saved


def test_resume_uses_the_saved_parm7_not_a_rebuilt_one(tmp_path, monkeypatch) -> None:
    builds = []
    _stub_mm_parm(monkeypatch, builds)
    _first, saved = _completed_mep(tmp_path, monkeypatch, builds)

    resumed = run_all(tmp_path, monkeypatch, extra=RESUME, parm=False)

    assert resumed.result.exit_code == 0, resumed.result.output[-3000:]
    assert "Reusing the verified MEP" in resumed.result.output
    assert len(builds) == 1  # a rebuild would have given other charges
    tsopt_args = resumed.child_args("tsopt")
    assert tsopt_args and {option(args, "--parm") for args in tsopt_args} == {str(saved)}


def test_resume_still_refuses_a_changed_input(tmp_path, monkeypatch) -> None:
    builds = []
    _stub_mm_parm(monkeypatch, builds)
    _completed_mep(tmp_path, monkeypatch, builds)

    resumed = run_all(tmp_path, monkeypatch, frames=(0.0, 2.1), extra=RESUME, parm=False)

    assert resumed.result.exit_code != 0
    assert "cannot reuse the saved MEP because inputs changed" in resumed.result.output


def test_resume_still_refuses_a_different_parm(tmp_path, monkeypatch) -> None:
    builds = []
    _stub_mm_parm(monkeypatch, builds)
    _completed_mep(tmp_path, monkeypatch, builds)
    other = write(tmp_path / "other.parm7", "%VERSION  other\n%FLAG CHARGE\n0.5000\n")

    resumed = run_all(
        tmp_path, monkeypatch, extra=(*RESUME, "--parm", str(other)), parm=False,
    )

    assert resumed.result.exit_code != 0
    assert "cannot reuse the saved MEP because path_request changed" in resumed.result.output


def test_resume_refuses_a_saved_parm7_that_no_longer_matches(tmp_path, monkeypatch) -> None:
    builds = []
    _stub_mm_parm(monkeypatch, builds)
    _first, saved = _completed_mep(tmp_path, monkeypatch, builds)
    saved.write_text("%VERSION  edited\n%FLAG CHARGE\n0.9000\n", encoding="utf-8")

    resumed = run_all(tmp_path, monkeypatch, extra=RESUME, parm=False)

    assert resumed.result.exit_code != 0
    assert "Saved parm7 does not match the completed MEP" in resumed.result.output
    assert len(builds) == 1
