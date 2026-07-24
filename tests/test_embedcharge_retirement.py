"""Fail-closed contract for the retired experimental embedding path."""

from __future__ import annotations

import re
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.core.embedcharge_policy import (
    EMBEDCHARGE_UNAVAILABLE,
    EmbedChargeUnavailableError,
    reject_retired_embedcharge_cli,
    validate_retired_embedcharge,
)


_WORKFLOW_MODULES = (
    "all",
    "dft",
    "freq",
    "irc",
    "opt",
    "path_opt",
    "path_search",
    "scan",
    "scan2d",
    "scan3d",
    "sp",
    "tsopt",
)
_FIXTURE_DIR = (
    Path(__file__).resolve().parents[1] / "hessian_ff" / "tests" / "data" / "small"
)


@pytest.mark.parametrize(
    ("embedcharge", "cutoff_requested"),
    [(True, False), (False, True), (True, True)],
)
def test_policy_rejects_activation_and_cutoff(embedcharge, cutoff_requested) -> None:
    with pytest.raises(EmbedChargeUnavailableError, match="double-counted"):
        validate_retired_embedcharge(
            embedcharge=embedcharge,
            cutoff_requested=cutoff_requested,
        )


def test_policy_allows_mechanical_embedding() -> None:
    validate_retired_embedcharge(embedcharge=False, cutoff_requested=False)
    reject_retired_embedcharge_cli({"embedcharge": False})


@pytest.mark.parametrize("module_name", _WORKFLOW_MODULES)
def test_every_compute_workflow_applies_effective_config_gate(module_name) -> None:
    source_path = (
        Path(__file__).resolve().parents[1]
        / "mlmm"
        / "workflows"
        / f"{module_name}.py"
    )
    source = source_path.read_text(encoding="utf-8")
    assert "reject_retired_embedcharge_cli(" in source
    assert "Unavailable in v0.3.3" in source


def _all_args(*extra: str) -> list[str]:
    structure = _FIXTURE_DIR / "complex.pdb"
    topology = _FIXTURE_DIR / "complex.parm7"
    return [
        "all",
        "-i",
        str(structure),
        str(structure),
        "--parm",
        str(topology),
        "-q",
        "0",
        "-m",
        "1",
        "--dry-run",
        *extra,
    ]


@pytest.mark.parametrize("extra", [("--embedcharge",), ("--embedcharge-cutoff", "8.0")])
def test_all_cli_rejects_retired_options_before_execution(extra) -> None:
    result = CliRunner().invoke(root_cli, _all_args(*extra))
    assert result.exit_code != 0
    assert EMBEDCHARGE_UNAVAILABLE in result.output


def test_yaml_activation_is_rejected_but_explicit_no_overrides_it(tmp_path) -> None:
    config = tmp_path / "embed.yaml"
    config.write_text(
        yaml.safe_dump({"calc": {"embedcharge": True}}),
        encoding="utf-8",
    )
    rejected = CliRunner().invoke(root_cli, _all_args("--config", str(config)))
    assert rejected.exit_code != 0
    assert EMBEDCHARGE_UNAVAILABLE in rejected.output

    from mlmm.workflows.all import _resolve_calculator_template

    mechanical = _resolve_calculator_template(
        config,
        backend=None,
        embedcharge=False,
        embedcharge_explicit=True,
        embedcharge_cutoff=None,
        link_atom_method=None,
        mm_backend=None,
        use_cmap=None,
    )
    reject_retired_embedcharge_cli(mechanical.materialize())


def test_mlmm_core_rejects_before_temporary_workspace(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from mlmm.backends import mlmm_calc

    touched = False

    def forbidden_workspace(*args, **kwargs):
        nonlocal touched
        touched = True
        raise AssertionError("workspace allocation must not run")

    monkeypatch.setattr(mlmm_calc.tempfile, "TemporaryDirectory", forbidden_workspace)
    with pytest.raises(EmbedChargeUnavailableError, match="double-counted"):
        mlmm_calc.MLMMCore(
            input_pdb="missing.pdb",
            real_parm7="missing.parm7",
            model_pdb="missing-model.pdb",
            embedcharge=True,
        )
    assert touched is False


def test_current_docs_and_skills_never_instruct_activation() -> None:
    repo = Path(__file__).resolve().parents[1]
    texts: list[tuple[Path, str]] = []
    for root in (repo / "docs", repo / "skills"):
        for path in root.rglob("*.md"):
            if path.name.lower().startswith("changelog"):
                continue
            texts.append((path, path.read_text(encoding="utf-8")))

    forbidden_phrases = (
        "Add `--embedcharge` to enable",
        "When `--embedcharge` is enabled",
        "With `--embedcharge`, xTB",
        "`--embedcharge` を追加すると",
        "`--embedcharge` で xTB 点電荷埋め込み補正を有効",
        "`embedcharge: true` enables",
    )
    for path, text in texts:
        for phrase in forbidden_phrases:
            assert phrase not in text, f"{path}: active instruction {phrase!r}"
        assert not re.search(
            r"(?m)^\s*mlmm\s+\S+.*(?:^|\s)--embedcharge(?:\s|$)",
            text,
        ), f"{path}: runnable embedcharge command remains"
        assert not re.search(
            r"(?m)^\s*embedcharge:\s*true(?:\s|#|$)",
            text,
        ), f"{path}: activating YAML example remains"
