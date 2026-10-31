"""Contract for the curated starter-snapshot reference page.

Positive: the runtime-owner parity validator accepts every current scalar (the
page's title, its link to the full YAML reference and the JA navigation label are
checked by ``.github/scripts/check_docs_contract.py``). Negative: mutating a
scalar, adding an ownerless scalar, or declaring an unused owner each fails and
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


def _owner_values():
    return gr._resolve_owner_values(gr.root_cli)


def _template():
    return yaml.safe_load(gr._ALL_TEMPLATE) or {}


# --------------------------------------------------------------------------- #
# Positive controls
# --------------------------------------------------------------------------- #
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
