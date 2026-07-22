from __future__ import annotations

from pathlib import Path


RELEASE_WORKFLOW = Path(__file__).parents[1] / ".github" / "workflows" / "release.yml"


def test_release_workflow_publishes_and_attaches_both_distributions() -> None:
    text = RELEASE_WORKFLOW.read_text(encoding="utf-8")

    assert "release:\n    types:\n      - published" in text
    assert "contents: write" in text
    assert "id-token: write" in text
    assert text.count("pypa/gh-action-pypi-publish@release/v1") == 1
    assert text.count('gh release upload "${{ github.event.release.tag_name }}"') == 1
    assert "dist/*.whl dist/*.tar.gz --clobber" in text
    assert text.index("Attach distribution to GitHub Release") < text.index(
        "Publish to PyPI (Trusted Publishing)"
    )
    assert text.rsplit("      - name: ", 1)[1].startswith(
        "Publish to PyPI (Trusted Publishing)"
    )
    assert "attestations: false" not in text
    assert "skip-existing: true" not in text
