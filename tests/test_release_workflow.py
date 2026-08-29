from __future__ import annotations

from pathlib import Path


RELEASE_WORKFLOW = Path(__file__).parents[1] / ".github" / "workflows" / "release.yml"
MANIFEST = Path(__file__).parents[1] / "MANIFEST.in"


def test_release_workflow_publishes_and_attaches_both_distributions() -> None:
    text = RELEASE_WORKFLOW.read_text(encoding="utf-8")

    assert "release:\n    types:\n      - published" in text
    assert "workflow_dispatch:" in text
    assert "Existing published release tag to build and publish" in text
    assert "contents: write" in text
    assert "id-token: write" in text
    assert text.count("pypa/gh-action-pypi-publish@release/v1") == 1
    assert text.count('gh release upload "$RELEASE_TAG"') == 1
    assert "dist/*.whl dist/*.tar.gz --clobber" in text
    assert 'ref: ${{ env.RELEASE_TAG }}' in text
    assert 'git merge-base --is-ancestor "$tag_commit" "$main_commit"' in text
    assert '--no-deps' not in text
    assert 'python -m pip install --force-reinstall "${wheels[0]}"' in text
    assert text.index("Attach distribution to GitHub Release") < text.index(
        "Publish to PyPI (Trusted Publishing)"
    )
    assert text.rsplit("      - name: ", 1)[1].startswith(
        "Publish to PyPI (Trusted Publishing)"
    )
    assert "attestations: false" not in text
    assert "skip-existing: true" not in text


def test_release_workflow_rejects_a_title_that_differs_from_the_tag() -> None:
    text = RELEASE_WORKFLOW.read_text(encoding="utf-8")

    assert (
        "RELEASE_TAG: ${{ github.event_name == 'release' "
        "&& github.event.release.tag_name || inputs.tag }}"
    ) in text
    assert "RELEASE_NAME: ${{ github.event.release.name }}" in text
    assert '[[ "$RELEASE_NAME" != "$RELEASE_TAG" ]]' in text
    assert 'gh release view "$RELEASE_TAG"' in text
    assert '--json isDraft --jq .isDraft' in text
    assert '--json isPrerelease --jq .isPrerelease' in text
    assert text.index("Verify release title matches tag") < text.index(
        "- name: Checkout"
    )


def test_sdist_prunes_local_smoke_products() -> None:
    text = MANIFEST.read_text(encoding="utf-8")

    assert "recursive-include tests *" in text
    assert "\nprune tests/smoke\n" in text
    assert "recursive-exclude tests/smoke" not in text
