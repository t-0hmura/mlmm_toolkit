"""Contract for release metadata: landing substitution + license (M67, P03).

Positive: the real repository's landing headers render ``v{{ release }}`` and the
CFF/pyproject SPDX license identities agree on ``GPL-3.0-only``. Negative: a
hardcoded landing literal or a mismatched CFF license fails.

NOTE (P03): ``cffconvert`` is not installed in the current environment, so the
``cffconvert --validate`` leg is exercised only when available. Freeze-time
confirmation that ``cffconvert --validate`` accepts ``GPL-3.0-only`` is still
required before release.
"""

from __future__ import annotations

import importlib.util
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


def test_real_license_identities_agree_on_gpl3_only() -> None:
    assert crv._check_license() == []
    assert crv._cff_license() == "GPL-3.0-only"
    assert crv._pyproject_license() == "GPL-3.0-only"


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
    # Exercised only if cffconvert is installed; otherwise the exact-equality
    # checks above stand in until freeze-time confirmation.
    if importlib.util.find_spec("cffconvert") is None:
        assert crv._cff_license() == crv._pyproject_license() == "GPL-3.0-only"
        return
    completed = subprocess.run(
        [sys.executable, "-m", "cffconvert", "--validate", "-i", str(REPO_ROOT / "CITATION.cff")],
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
