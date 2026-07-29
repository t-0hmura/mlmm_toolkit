from __future__ import annotations

import shlex
import sys
from pathlib import Path

from mlmm.core.utils import load_pdb_atom_metadata
from mlmm.workflows.all import _parse_scan_lists_literals


BEZA = Path(__file__).parents[1] / "examples" / "beza"
SCRIPTS = Path(__file__).parents[1] / ".github" / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import docs_command_contract as dc  # noqa: E402


STATES = ("1.R.pdb", "2.IM.pdb", "3.P.pdb")
IDENTITY_FIELDS = (
    "is_hetatm",
    "chain",
    "resname",
    "resseq",
    "icode",
    "name",
    "altloc",
    "element",
)


def _identity(path: Path) -> list[tuple[object, ...]]:
    return [
        tuple(atom.get(field) for field in IDENTITY_FIELDS)
        for atom in load_pdb_atom_metadata(path)
    ]


def _coordinates(path: Path) -> list[str]:
    return [
        line[30:54]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]


def test_beza_states_share_one_ordered_topology() -> None:
    identities = [_identity(BEZA / name) for name in STATES]
    assert len(identities[0]) == 9215
    assert identities[1:] == [identities[0], identities[0]]

    residues = {identity[2] for identity in identities[0]}
    assert {"SAM", "GPP", "MG", "GLU"} <= residues


def test_beza_states_have_distinct_coordinates() -> None:
    coordinates = [_coordinates(BEZA / name) for name in STATES]
    assert coordinates[0] != coordinates[1]
    assert coordinates[0] != coordinates[2]
    assert coordinates[1] != coordinates[2]


def test_beza_runner_preserves_endpoint_and_scan_contract() -> None:
    commands = dc.extract_shell_commands([BEZA / "run.sh"])
    assert len(commands) == 2
    endpoint, scan = (shlex.split(command.text) for command in commands)

    assert endpoint[endpoint.index("-i") + 1 : endpoint.index("-c")] == [
        "$script_dir/1.R.pdb",
        "$script_dir/3.P.pdb",
    ]
    assert scan[scan.index("-i") + 1 : scan.index("-c")] == [
        "$script_dir/1.R.pdb",
    ]
    for tokens, out_dir in ((endpoint, "result_mep"), (scan, "result_scan")):
        assert tokens[tokens.index("-c") + 1] == "SAM,GPP,MG"
        assert tokens[tokens.index("-l") + 1] == "SAM:1,GPP:-3"
        assert "--tsopt" in tokens
        assert "--thermo" in tokens
        assert tokens[tokens.index("--out-dir") + 1] == out_dir

    scan_index = scan.index("--scan-lists")
    scan_stages = scan[scan_index + 1 : scan_index + 3]
    metadata = load_pdb_atom_metadata(BEZA / "1.R.pdb")
    stages = _parse_scan_lists_literals(
        scan_stages,
        atom_meta=metadata,
        one_based=True,
    )
    assert [len(stage) for stage in stages] == [2, 2]
