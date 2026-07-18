"""Regression tests for coherent PDB alternate-location selection."""

from pathlib import Path

import pytest

from mlmm.io.pdb_fix import process_block


def _atoms(lines: list[str]) -> list[str]:
    return [line for line in process_block(lines) if line.startswith("ATOM")]


def test_altloc_selection_does_not_create_hybrid_residue() -> None:
    lines = [
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N\n",
        "ATOM      2  CA AALA A   1       1.000   0.000   0.000  0.60 10.00           C\n",
        "ATOM      3  CG AALA A   1       2.500   1.000   0.000  0.60 10.00           C\n",
        "ATOM      4  CA BALA A   1       1.100   0.000   0.000  0.40 10.00           C\n",
        "ATOM      5  CD BALA A   1       2.600  -1.000   0.000  0.40 10.00           C\n",
        "END\n",
    ]

    atoms = _atoms(lines)

    assert [line[12:16].strip() for line in atoms] == ["N", "CA", "CG"]
    assert all(line[16] == " " for line in atoms)


def test_altloc_selection_uses_residue_mean_occupancy() -> None:
    lines = [
        "ATOM      1  CA AALA A   1       1.000   0.000   0.000  0.70 10.00           C\n",
        "ATOM      2  CB AALA A   1       1.500   1.000   0.000  0.50 10.00           C\n",
        "ATOM      3  CA BALA A   1       1.100   0.000   0.000  0.30 10.00           C\n",
        "ATOM      4  CB BALA A   1       1.600  -1.000   0.000  0.60 10.00           C\n",
    ]

    assert [int(line[6:11]) for line in _atoms(lines)] == [1, 2]


def test_altloc_selection_prefers_parsed_mean_over_missing_only_label() -> None:
    missing_a = (
        "ATOM      1  CA AALA A   1       1.000   0.000   0.000  0.70 10.00           C\n"
    )
    missing_a = missing_a[:54] + "      " + missing_a[60:]
    parsed_b = (
        "ATOM      2  CA BALA A   1       1.100   0.000   0.000  0.10 10.00           C\n"
    )

    assert [int(line[6:11]) for line in _atoms([missing_a, parsed_b])] == [2]


@pytest.mark.parametrize(
    ("occupancy_a", "occupancy_b", "expected_serial"),
    [
        (None, 0.10, 2),
        (0.00, None, 1),
        (0.25, 0.25, 1),
        (None, None, 1),
    ],
)
def test_typed_and_text_altloc_paths_share_missing_occupancy_policy(
    tmp_path: Path,
    occupancy_a: float | None,
    occupancy_b: float | None,
    expected_serial: int,
) -> None:
    from mlmm.io.structure_formats import read_pdb_atom_sites

    def line(serial: int, label: str, x: float, occupancy: float | None) -> str:
        occupancy_field = "      " if occupancy is None else f"{occupancy:6.2f}"
        return (
            f"ATOM  {serial:5d}  CA {label}ALA A   1    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{occupancy_field}{10.0:6.2f}"
            "           C\n"
        )

    lines = [line(1, "A", 1.0, occupancy_a), line(2, "B", 2.0, occupancy_b)]
    source = tmp_path / "altloc.pdb"
    source.write_text("".join(lines) + "END\n", encoding="utf-8")

    typed, _nonstandard = read_pdb_atom_sites(source)
    text = _atoms(lines)

    assert len(typed) == len(text) == 1
    assert typed[0].x == pytest.approx(float(expected_serial))
    assert int(text[0][6:11]) == expected_serial
    assert text[0][16] == " "


@pytest.mark.parametrize("nonfinite", ["nan", "inf", "-inf"])
def test_nonfinite_occupancy_is_missing_in_both_altloc_paths(
    tmp_path: Path,
    nonfinite: str,
) -> None:
    from mlmm.io.structure_formats import read_pdb_atom_sites

    def line(serial: int, label: str, x: float, occupancy_field: str) -> str:
        return (
            f"ATOM  {serial:5d}  CA {label}ALA A   1    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{occupancy_field:>6s}{10.0:6.2f}"
            "           C\n"
        )

    lines = [line(1, "A", 1.0, nonfinite), line(2, "B", 2.0, "0.00")]
    source = tmp_path / "nonfinite.pdb"
    source.write_text("".join(lines) + "END\n", encoding="utf-8")

    typed, _ = read_pdb_atom_sites(source)
    text = _atoms(lines)

    assert typed[0].x == pytest.approx(2.0)
    assert int(text[0][6:11]) == 2


def test_labelled_atom_supersedes_high_occupancy_blank_in_both_paths(
    tmp_path: Path,
) -> None:
    from mlmm.io.structure_formats import read_pdb_atom_sites

    def line(serial: int, label: str, x: float, occupancy: float) -> str:
        return (
            f"ATOM  {serial:5d}  CA {label}ALA A   1    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{occupancy:6.2f}{10.0:6.2f}"
            "           C\n"
        )

    lines = [
        line(1, " ", 1.0, 1.00),
        line(2, "A", 2.0, 0.10),
        line(3, "B", 3.0, 0.05),
    ]
    source = tmp_path / "shared-and-altloc.pdb"
    source.write_text("".join(lines) + "END\n", encoding="utf-8")

    typed, _ = read_pdb_atom_sites(source)
    text = _atoms(lines)

    assert typed[0].x == pytest.approx(2.0)
    assert int(text[0][6:11]) == 2
