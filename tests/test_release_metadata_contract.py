"""Contract for release metadata, landing substitution, and licensing.

Positive: the real repository's landing headers render ``v{{ release }}`` and the
CFF/pyproject SPDX license identities agree on ``GPL-3.0-or-later``. Negative: a
hardcoded landing literal or a mismatched CFF license fails.

The ``cffconvert --validate`` leg runs when its console entry point is
available, matching the release workflow invocation.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = REPO_ROOT / ".github" / "scripts"
for _p in (str(REPO_ROOT), str(SCRIPTS)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import check_release_versions as crv  # noqa: E402


def test_release_checker_passes_on_real_repo(monkeypatch) -> None:
    monkeypatch.setattr(sys, "argv", ["check_release_versions.py"])
    assert crv.main() == 0


def test_real_landing_pages_use_substitution() -> None:
    assert crv._check_landing_pages() == []


def test_real_license_identities_agree_on_gpl3_or_later() -> None:
    assert crv._check_license() == []
    assert crv._cff_license() == "GPL-3.0-or-later"
    assert crv._pyproject_license() == "GPL-3.0-or-later"


def test_primary_source_title_is_identical_across_citation_surfaces() -> None:
    surfaces = (
        "CITATION.cff",
        "README.md",
        "docs/index.md",
        "docs/ja/index.md",
        "mlmm/io/summary.py",
        "examples/mlmm_colab.ipynb",
    )
    for relative in surfaces:
        text = (REPO_ROOT / relative).read_text(encoding="utf-8")
        assert "Toward Accelerated" in text, relative
        assert "Towards Accelerated" not in text, relative


def test_hardcoded_landing_literal_fails(tmp_path, monkeypatch) -> None:
    good = tmp_path / "index.md"
    good.write_text("# Docs\n\n*Version: v{{ release }}* — text\n", encoding="utf-8")
    monkeypatch.setattr(crv, "LANDING_PAGES", (good,))
    assert crv._check_landing_pages() == []

    bad = tmp_path / "bad.md"
    bad.write_text("# Docs\n\n*Version: v0.3.3* — text\n", encoding="utf-8")
    monkeypatch.setattr(crv, "LANDING_PAGES", (bad,))
    errors = crv._check_landing_pages()
    assert errors
    assert any("hardcode a version literal" in e for e in errors)


def test_mismatched_cff_license_fails(monkeypatch) -> None:
    monkeypatch.setattr(crv, "_cff_license", lambda: "GPL-3.0")
    errors = crv._check_license()
    assert errors
    assert any("!=" in e for e in errors)


def test_cffconvert_accepts_license_when_available() -> None:
    # Exercise the same console entry point used by the release workflow.
    executable = shutil.which("cffconvert")
    if executable is None:
        assert (
            crv._cff_license()
            == crv._pyproject_license()
            == "GPL-3.0-or-later"
        )
        return
    completed = subprocess.run(
        [executable, "--validate"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
