"""Contract for the curated starter-snapshot reference page (M71, P02).

Positive: the generated page is a curated, non-exhaustive snapshot that links
the full YAML reference, the JA navigation labels it a starter snapshot, and the
runtime-owner parity validator accepts every current scalar. Negative: mutating
a scalar, adding an ownerless scalar, or declaring an unused owner each fails and
names the exact dotted path.
"""

from __future__ import annotations

import copy
import sys
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = REPO_ROOT / ".github" / "scripts"
for _p in (str(REPO_ROOT), str(SCRIPTS)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import generate_reference as gr  # noqa: E402

GENERATED = REPO_ROOT / "docs" / "reference" / "yaml.md"
JA_INDEX = REPO_ROOT / "docs" / "ja" / "index.md"


def _owner_values():
    return gr._resolve_owner_values(gr.root_cli)


def _template():
    return yaml.safe_load(gr._ALL_TEMPLATE) or {}


# --------------------------------------------------------------------------- #
# Positive controls
# --------------------------------------------------------------------------- #
def test_generated_page_is_curated_non_exhaustive_and_links_full_reference() -> None:
    text = GENERATED.read_text(encoding="utf-8")
    assert text.startswith("# Curated `mlmm all` Starter Snapshot")
    assert "non-exhaustive" in text
    assert "(../yaml-reference.md)" in text
    assert "## Included Sections" in text
    assert "YAML Schema" not in text
    assert "Top-level Keys" not in text


def test_ja_navigation_labels_starter_snapshot_not_schema() -> None:
    line = next(
        ln
        for ln in JA_INDEX.read_text(encoding="utf-8").splitlines()
        if "../reference/yaml.md" in ln
    )
    assert "スキーマ" not in line
    assert "スターター" in line


def test_current_starter_snapshot_passes_parity() -> None:
    gr._validate_starter_snapshot(_template(), _owner_values())


def test_every_scalar_has_an_owner_and_equals_it() -> None:
    owners = _owner_values()
    template_items = dict(gr._iter_scalar_items(_template()))
    assert set(template_items) == set(owners)
    for path, value in template_items.items():
        _label, owner_value = owners[path]
        assert gr._scalar_equal(value, owner_value), path


# --------------------------------------------------------------------------- #
# Negative controls
# --------------------------------------------------------------------------- #
def test_mutated_scalar_value_fails_with_exact_path() -> None:
    owners = _owner_values()
    data = copy.deepcopy(_template())
    data["calc"]["backend"] = "not-a-backend"
    with pytest.raises(RuntimeError) as exc:
        gr._validate_starter_snapshot(data, owners)
    assert "calc.backend" in str(exc.value)


def test_mutated_thermo_temperature_fails_with_exact_path() -> None:
    owners = _owner_values()
    data = copy.deepcopy(_template())
    data["thermo"]["temperature"] = 300.0
    with pytest.raises(RuntimeError) as exc:
        gr._validate_starter_snapshot(data, owners)
    assert "thermo.temperature" in str(exc.value)


def test_ownerless_scalar_fails() -> None:
    owners = _owner_values()
    data = copy.deepcopy(_template())
    data["freq"]["brand_new_knob"] = 7
    with pytest.raises(RuntimeError) as exc:
        gr._validate_starter_snapshot(data, owners)
    assert "freq.brand_new_knob" in str(exc.value)
    assert "ownerless" in str(exc.value)


def test_unused_owner_declaration_fails() -> None:
    owners = dict(_owner_values())
    owners["scan.nonexistent"] = ("`SOME_KW[\"x\"]`", 1)
    with pytest.raises(RuntimeError) as exc:
        gr._validate_starter_snapshot(_template(), owners)
    assert "scan.nonexistent" in str(exc.value)
    assert "stale owner declaration" in str(exc.value)
