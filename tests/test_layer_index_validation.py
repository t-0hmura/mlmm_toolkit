"""Validation of YAML layer atom-index normalization."""

from __future__ import annotations

import click
import pytest

from mlmm.workflows._opt_freq_common import _convert_yaml_layer_atoms_1to0


def test_layer_indices_are_deduplicated_and_converted() -> None:
    config = {"hess_mm_atoms": [3, 1, 3]}
    _convert_yaml_layer_atoms_1to0(config)
    assert config["hess_mm_atoms"] == [0, 2]


@pytest.mark.parametrize(
    "value",
    [
        "1,2,3",
        [1, "2"],
        [1.9, 3],
        {"one": 1},
        [0],
        [-1],
        [True],
    ],
)
def test_invalid_layer_indices_fail_closed(value) -> None:
    config = {"movable_mm_atoms": value}
    with pytest.raises(click.BadParameter, match=r"calc\.movable_mm_atoms"):
        _convert_yaml_layer_atoms_1to0(config)
