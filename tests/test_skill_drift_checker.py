from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / ".github" / "scripts" / "check_skill_drift.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("mlmm_check_skill_drift", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_skill_drift_accepts_live_flag_and_rejects_fabricated_flag(
    tmp_path, monkeypatch
) -> None:
    checker = _load_checker()
    skill_dir = tmp_path / "skills" / "mlmm-cli"
    skill_dir.mkdir(parents=True)
    skill = skill_dir / "example.md"
    monkeypatch.setattr(checker, "REPO_ROOT", tmp_path)

    flags = checker._collect_flag_union()
    skill.write_text("Use `--max-cycles`.\n", encoding="utf-8")
    assert checker._scan_file(skill, flags) == []

    skill.write_text("Use `--definitely-not-a-real-mlmm-flag`.\n", encoding="utf-8")
    warnings = checker._scan_file(skill, flags)
    assert len(warnings) == 1
    assert "unknown flag --definitely-not-a-real-mlmm-flag" in warnings[0]


def test_skill_drift_main_fails_on_detected_drift(tmp_path, monkeypatch) -> None:
    checker = _load_checker()
    skill_dir = tmp_path / "skills" / "mlmm-cli"
    skill_dir.mkdir(parents=True)
    (skill_dir / "example.md").write_text(
        "Use `--definitely-not-a-real-mlmm-flag`.\n", encoding="utf-8"
    )
    monkeypatch.setattr(checker, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(checker, "SKILLS_DIR", tmp_path / "skills")

    assert checker.main() == 1
