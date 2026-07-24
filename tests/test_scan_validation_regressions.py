"""Fail-closed validation for scan axes and scientific eligibility."""

import click
import pytest

from mlmm.core.utils import (
    parse_scan_list_quads,
    parse_scan_list_triples,
    parse_scan_spec_stages,
)


@pytest.mark.parametrize(
    "raw",
    [
        "[]",
        "[(1, 1, 1.5)]",
        "[(1, 2, 1.5), (2, 1, 1.6)]",
    ],
)
def test_scan_stage_rejects_empty_self_and_duplicate_pairs(raw) -> None:
    with pytest.raises(click.BadParameter):
        parse_scan_list_triples(
            raw,
            one_based=False,
            atom_meta=None,
            option_name="--scan-lists",
        )


@pytest.mark.parametrize(
    "raw",
    [
        "[(1, 1, 1.0, 2.0), (2, 3, 1.0, 2.0)]",
        "[(1, 2, 1.0, 2.0), (2, 1, 1.0, 2.0)]",
    ],
)
def test_grid_scan_rejects_self_and_duplicate_axes(raw) -> None:
    with pytest.raises(click.BadParameter):
        parse_scan_list_quads(
            raw,
            expected_len=2,
            one_based=False,
            atom_meta=None,
            option_name="--scan-lists",
        )


def test_scan_spec_expands_bidirectional_stage_with_reset_markers(tmp_path) -> None:
    spec = tmp_path / "scan.yaml"
    spec.write_text(
        "one_based: true\n"
        "stages:\n"
        "  - [[1, 2, 1.2, 1.8]]\n",
        encoding="utf-8",
    )

    stages, one_based, snapshots, resets = parse_scan_spec_stages(
        spec,
        one_based_default=False,
        atom_meta=None,
        return_bidirectional_markers=True,
    )

    assert one_based is True
    assert stages == [[(0, 1, 1.2)], [(0, 1, 1.8)]]
    assert snapshots == frozenset({0})
    assert resets == frozenset({1})
