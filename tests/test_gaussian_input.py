"""Ordinary Gaussian parsing and bond-summary integration regressions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

from mlmm.cli import cli as root_cli


def _write_ordinary_gjf(path: Path, coordinates: list[tuple[str, float, float, float]]) -> None:
    rows = "\n".join(
        f"{element:<2s} {x:10.4f} {y:10.4f} {z:10.4f}"
        for element, x, y, z in coordinates
    )
    path.write_text(
        "%mem=1GB\n"
        "#p hf/sto-3g\n\n"
        "ordinary molecule\n\n"
        f"-1 2\n{rows}\n\n",
        encoding="utf-8",
    )


def test_parse_ordinary_gaussian_preserves_order_coordinates_and_state(
    tmp_path: Path,
) -> None:
    from mlmm.io.gaussian_input import parse_gaussian_input

    source = tmp_path / "ordinary.gjf"
    _write_ordinary_gjf(
        source,
        [("C", 0.0, 0.0, 0.0), ("Cl", 1.7, 0.0, 0.0)],
    )

    parsed = parse_gaussian_input(source)

    assert parsed.elements == ("C", "Cl")
    assert np.allclose(parsed.coordinates, [[0.0, 0.0, 0.0], [1.7, 0.0, 0.0]])
    assert (parsed.charge, parsed.multiplicity) == (-1, 2)


def test_bond_summary_cli_reads_ordinary_gaussian_inputs(tmp_path: Path) -> None:
    reactant = tmp_path / "reactant.gjf"
    product = tmp_path / "product.gjf"
    spectators = [
        ("H", 20.0, 0.0, 0.0),
        ("H", 25.0, 0.0, 0.0),
        ("H", 30.0, 0.0, 0.0),
    ]
    _write_ordinary_gjf(
        reactant,
        [("C", 0.0, 0.0, 0.0), *spectators, ("F", 5.0, 0.0, 0.0), ("Cl", 1.7, 0.0, 0.0)],
    )
    _write_ordinary_gjf(
        product,
        [("C", 0.0, 0.0, 0.0), *spectators, ("F", 1.3, 0.0, 0.0), ("Cl", 5.0, 0.0, 0.0)],
    )

    result = CliRunner().invoke(
        root_cli,
        ["bond-summary", "-i", str(reactant), "-i", str(product), "--device", "cpu"],
    )

    assert result.exit_code == 0, result.output
    assert "C1-F5" in result.output
    assert "C1-Cl6" in result.output


def test_specialized_oniom_parser_remains_separate(tmp_path: Path) -> None:
    from mlmm.io.gaussian_input import parse_gaussian_input
    from mlmm.workflows.oniom_import import _parse_gaussian_oniom

    source = tmp_path / "oniom.gjf"
    source.write_text(
        "#p oniom(hf/sto-3g:amber)\n\n"
        "ONIOM molecule\n\n"
        "0 1 0 1 0 1\n"
        "C-CT-0.0 0 0.0 0.0 0.0 H\n"
        "H-HC-0.0 0 1.0 0.0 0.0 L\n\n",
        encoding="utf-8",
    )

    coords, elements, qm_indices, movable_indices, charge, multiplicity = (
        _parse_gaussian_oniom(source)
    )

    assert elements == ["C", "H"]
    assert np.allclose(coords, [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    assert qm_indices == {0}
    assert movable_indices == {0, 1}
    assert (charge, multiplicity) == (0, 1)
    with pytest.raises(ValueError, match="ordinary two-integer"):
        parse_gaussian_input(source)
