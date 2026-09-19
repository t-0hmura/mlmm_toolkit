from __future__ import annotations

import numpy as np
from click.testing import CliRunner

from pysisyphus.constants import ANG2BOHR

from mlmm.core.utils import parse_scan_list_quads, parse_scan_list_triples
from mlmm.domain.scan_coordinates import coordinate_delta, coordinate_value
from mlmm.workflows.all import _parse_scan_lists_literals
from mlmm.workflows.restraints import harmonic_internal_energy_forces_hessian
from mlmm.workflows.scan import cli as scan_cli
from mlmm.workflows.scan2d import cli as scan2d_cli
from mlmm.workflows.scan3d import cli as scan3d_cli
from mlmm.cli import cli as root_cli


def _coords() -> np.ndarray:
    return np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
         [1.0, 1.0, 0.0], [2.0, 1.0, 1.0]], dtype=float,
    ) * ANG2BOHR


def test_target_and_range_grammars_cover_all_coordinate_kinds() -> None:
    targets, _ = parse_scan_list_triples(
        "[(1,2,1.2),(1,2,3,100),(1,2,3,4,-60)]",
        one_based=True, atom_meta=None, option_name="--scan-lists",
    )
    assert targets == [(0, 1, 1.2), (0, 1, 2, 100.0), (0, 1, 2, 3, -60.0)]
    ranges, _ = parse_scan_list_quads(
        "[(1,2,1.0,1.5),(1,2,3,80,120),(1,2,3,4,-180,180)]",
        expected_len=3, one_based=True, atom_meta=None,
        option_name="--scan-lists",
    )
    assert [len(entry) for entry in ranges] == [4, 5, 6]
    assert _parse_scan_lists_literals(
        ("[(1,2,3,100),(1,2,3,4,-60)]",), one_based=True,
    )[0] == [(1, 2, 3, 100.0), (1, 2, 3, 4, -60.0)]


def test_internal_coordinate_bias_force_and_hessian() -> None:
    coords = _coords()
    np.testing.assert_allclose(coordinate_value(coords, (0, 1, 1.0)), 1.0)
    np.testing.assert_allclose(coordinate_value(coords, (0, 1, 2, 90.0)), 90.0)
    np.testing.assert_allclose(coordinate_value(coords, (0, 1, 2, 3, 135.0)), 135.0)
    restraints = [(0, 1, 2, 100.0), (0, 1, 2, 3, 45.0)]
    energy, forces, hessian = harmonic_internal_energy_forces_hessian(coords, 5.0, restraints)
    assert energy > 0.0 and hessian is not None
    np.testing.assert_allclose(hessian, hessian.T, atol=1e-10)
    numeric_gradient = np.zeros(coords.size)
    flat = coords.reshape(-1)
    step = 1.0e-6
    for index in range(flat.size):
        plus = flat.copy(); plus[index] += step
        minus = flat.copy(); minus[index] -= step
        e_plus = harmonic_internal_energy_forces_hessian(plus.reshape(-1, 3), 5.0, restraints, need_hessian=False)[0]
        e_minus = harmonic_internal_energy_forces_hessian(minus.reshape(-1, 3), 5.0, restraints, need_hessian=False)[0]
        numeric_gradient[index] = (e_plus - e_minus) / (2.0 * step)
    np.testing.assert_allclose(forces, -numeric_gradient, atol=2e-8)


def test_dihedral_delta_crosses_periodic_boundary_in_short_direction() -> None:
    assert coordinate_delta("dihedral", -170.0, 170.0) == 20.0
    assert coordinate_delta("dihedral", 170.0, -170.0) == -20.0


def test_scan_help_publishes_internal_coordinate_ranges() -> None:
    runner = CliRunner()
    for command_name in ("scan", "scan2d", "scan3d"):
        result = runner.invoke(root_cli, [command_name, "--help-advanced"])
        assert result.exit_code == 0
        assert "--max-angle-step-size" in result.output
        assert "--max-dihedral-step-size" in result.output
