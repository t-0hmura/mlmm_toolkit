"""Contract for the expanded public-Markdown link checker (M66).

Positive: the real repository passes and the explicit roots include README, a
skills page, and docs pages. Negative: an independent broken link in each of a
README, a skill, and an example fixture is reported.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = REPO_ROOT / ".github" / "scripts"
for _p in (str(REPO_ROOT), str(SCRIPTS)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import check_markdown_links as cml  # noqa: E402


def test_link_checker_passes_on_real_repo() -> None:
    assert cml.main() == 0


def test_public_roots_include_readme_docs_and_skills() -> None:
    pages = cml.public_markdown_paths()
    rels = {str(p.relative_to(REPO_ROOT)) for p in pages}
    assert "README.md" in rels
    assert any(r.startswith("docs/") for r in rels)
    assert any(r.startswith("skills/") for r in rels)


def test_broken_links_reported_independently_per_public_page(tmp_path: Path) -> None:
    fixtures = {
        "README.md": "[missing](does_not_exist_readme.md)\n",
        "skill.md": "[missing](does_not_exist_skill.md)\n",
        "example.md": "[missing](does_not_exist_example.md)\n",
    }
    errors: list[str] = []
    for name, text in fixtures.items():
        page = tmp_path / name
        page.write_text(text, encoding="utf-8")
        cml._check_path(page, errors)
    assert len(errors) == 3
    assert all("broken local link" in e for e in errors)


def test_toctree_only_checked_on_docs_pages(tmp_path: Path) -> None:
    # A non-docs page with a toctree-looking block must not raise a toctree error.
    page = tmp_path / "skill.md"
    page.write_text(
        "```{toctree}\nmissing_target\n```\n", encoding="utf-8"
    )
    errors: list[str] = []
    cml._check_path(page, errors)
    assert errors == []
