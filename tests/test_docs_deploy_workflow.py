"""Guard the automatic documentation publication dependency."""

from __future__ import annotations

from pathlib import Path


WORKFLOW = Path(__file__).parents[1] / ".github" / "workflows" / "docs.yml"


def test_pages_deploy_waits_for_docs_quality_on_the_same_main_sha() -> None:
    text = WORKFLOW.read_text(encoding="utf-8")

    assert 'workflows: ["Docs Quality"]' in text
    assert "types: [completed]" in text
    assert "github.event.workflow_run.conclusion == 'success'" in text
    assert "github.event.workflow_run.event == 'push'" in text
    assert "github.event.workflow_run.head_branch == 'main'" in text
    assert "github.event.workflow_run.head_sha == github.sha" in text
    assert (
        "github.event.workflow_run.head_sha || github.sha"
        in text
    )
    assert "needs: build" in text
    assert "\n  push:\n" not in text


def test_manual_dispatch_remains_an_explicit_operator_override() -> None:
    text = WORKFLOW.read_text(encoding="utf-8")

    assert "workflow_dispatch:" in text
    assert "github.event_name == 'workflow_dispatch'" in text
