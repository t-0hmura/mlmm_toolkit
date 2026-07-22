from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / ".github" / "scripts" / "check_skill_commands.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("mlmm_check_skill_commands", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_command_checker_requires_topology_for_compute_examples() -> None:
    checker = _load_checker()
    contracts = checker._collect_subcommand_contracts()

    issues = checker._check_command("mlmm sp -i system.pdb", contracts)
    assert any("missing topology option" in issue for issue in issues)
    assert checker._check_command(
        "mlmm sp -i system.pdb --parm system.parm7 -q 0", contracts
    ) == []


def test_command_checker_requires_charge_for_compute_examples() -> None:
    checker = _load_checker()
    contracts = checker._collect_subcommand_contracts()

    issues = checker._check_command(
        "mlmm opt -i system.pdb --parm system.parm7", contracts
    )
    assert any("missing charge option" in issue for issue in issues)
    assert checker._check_command(
        "mlmm opt -i system.pdb --parm system.parm7 -q 0", contracts
    ) == []
    assert checker._check_command(
        "mlmm opt -i system.pdb --parm system.parm7 --config run.yaml", contracts
    ) == []


def test_command_checker_requires_matching_references_for_xyz_inputs() -> None:
    checker = _load_checker()
    contracts = checker._collect_subcommand_contracts()

    issues = checker._check_command(
        "mlmm path-opt -i reactant.xyz product.xyz --parm system.parm7 "
        "--ref-pdb reactant.pdb -q 0",
        contracts,
    )
    assert "XYZ input requires 2 corresponding --ref-pdb value(s)" in issues
    assert checker._check_command(
        "mlmm path-opt -i reactant.xyz product.xyz --parm system.parm7 "
        "--ref-pdb reactant.pdb --ref-pdb product.pdb -q 0",
        contracts,
    ) == []


def test_command_checker_accepts_scan3d_plot_only_mode() -> None:
    checker = _load_checker()
    contracts = checker._collect_subcommand_contracts()

    assert checker._check_command(
        "mlmm scan3d --csv surface.csv -o result_plot", contracts
    ) == []
