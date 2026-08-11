"""Contracts for strict XYZ frame and energy parsing."""

from __future__ import annotations

from pathlib import Path

import pytest

from mlmm.io.xyz_trajectory import (
    parse_xyz_energy_comment,
    read_xyz_trajectory,
)


def test_energy_comment_prefers_explicit_key_and_converts_units() -> None:
    assert parse_xyz_energy_comment("frame 8 E=-1.25 Ha") == (
        -1.25,
        "explicit-Ha",
    )
    energy, source = parse_xyz_energy_comment("Energy=27.211386245988 eV")
    assert energy == pytest.approx(1.0)
    assert source == "explicit-eV"
    energy, source = parse_xyz_energy_comment("Energy=27.211386245988 [eV]")
    assert energy == pytest.approx(1.0)
    assert source == "explicit-eV"
    energy, source = parse_xyz_energy_comment("E=627.5094740631 (kcal/mol)")
    assert energy == pytest.approx(1.0)
    assert source == "explicit-kcal/mol"
    energy, source = parse_xyz_energy_comment(
        "frame=3 E=27.211386245988; units=eV"
    )
    assert energy == pytest.approx(1.0)
    assert source == "explicit-eV"
    energy, source = parse_xyz_energy_comment(
        "frame=3 E=27.211386245988, units=eV"
    )
    assert energy == pytest.approx(1.0)
    assert source == "explicit-eV"
    energy, source = parse_xyz_energy_comment("E=1.0 Eh; unit=Ha")
    assert energy == pytest.approx(1.0)
    assert source == "explicit-Ha"
    energy, source = parse_xyz_energy_comment(
        'Lattice="10 0 0 0 10 0 0 0 10" Properties=species:S:1:pos:R:3 energy=27.211386245988'
    )
    assert energy == pytest.approx(1.0)
    assert source == "explicit-extxyz-eV"


def test_energy_comment_never_uses_integer_counter_as_energy() -> None:
    assert parse_xyz_energy_comment("-100.0") == (-100.0, "bare-assumed-Ha")
    assert parse_xyz_energy_comment("0 -100.0")[0] is None
    assert parse_xyz_energy_comment("frame 0")[0] is None
    assert parse_xyz_energy_comment("-100.0 -99.0")[0] is None
    assert parse_xyz_energy_comment("mode 1 -123.45 cm-1 frame=4")[0] is None


def test_energy_comment_accepts_bundled_pysisyphus_trajectory_format() -> None:
    assert parse_xyz_energy_comment("       -599.97249452 , ") == (
        -599.97249452,
        "pysisyphus-legacy-Ha",
    )
    assert parse_xyz_energy_comment(" -1.25 , cycle 3") == (
        -1.25,
        "pysisyphus-legacy-Ha",
    )
    assert parse_xyz_energy_comment("frame 3, -1.25")[0] is None


def test_conflicting_explicit_energies_are_structure_only() -> None:
    assert parse_xyz_energy_comment("E=-1.0 Ha Energy=-2.0 Ha") == (
        None,
        "conflicting-explicit",
    )
    assert parse_xyz_energy_comment("E=1e999 Ha") == (
        None,
        "nonfinite-explicit",
    )


@pytest.mark.parametrize(
    "comment, provenance",
    [
        ("E=-1.0 kJ/mol", "unsupported-unit"),
        ("Energy=-1.0 J/mol", "unsupported-unit"),
        ("E=-1.0foo", "malformed-explicit"),
        ("Energy=27.2 bananas", "unsupported-unit"),
        ("E=-1.0 [kJ/mol]", "unsupported-unit"),
        ("E=-1.0 [eV", "malformed-explicit"),
        ("E=-1.0 (eV]", "malformed-explicit"),
        ("E=-1.0; unit=kJ/mol", "unsupported-unit"),
        ("E=-1.0; unit=", "malformed-explicit"),
        ("E=-1.0 Ha; units=eV", "conflicting-explicit"),
        ("E=-1.0e+", "malformed-explicit"),
    ],
)
def test_energy_comment_rejects_unsupported_or_malformed_suffixes(
    comment: str,
    provenance: str,
) -> None:
    assert parse_xyz_energy_comment(comment) == (None, provenance)


def test_xyz_reader_validates_blocks_and_can_keep_missing_energies(
    tmp_path: Path,
) -> None:
    trajectory = tmp_path / "trajectory.xyz"
    trajectory.write_text(
        "1\nframe 0\nH 0.0 0.0 0.0\n"
        "1\nEnergy=-99.5\nH 0.1 0.0 0.0\n",
        encoding="utf-8",
    )
    parsed = read_xyz_trajectory(trajectory)
    assert len(parsed["frames"]) == 2
    assert parsed["energies_ha"] == [None, -99.5]
    assert parsed["energy_provenance"] == [
        "ambiguous",
        "explicit-assumed-Ha",
    ]
    with pytest.raises(RuntimeError, match="frame 1"):
        read_xyz_trajectory(trajectory, require_energies=True)


@pytest.mark.parametrize(
    "contents, message",
    [
        ("-1\ncomment\n", "invalid atom count"),
        ("2\ncomment\nH 0 0 0\n", "Incomplete XYZ frame"),
        ("1\ncomment\nH x 0 0\n", "Non-numeric coordinate"),
        ("1\ncomment\nH nan 0 0\n", "Non-finite coordinate"),
        ("1\ncomment\nH 0 inf 0\n", "Non-finite coordinate"),
        ("1\ncomment\nH 0 0 -inf\n", "Non-finite coordinate"),
        ("not-an-xyz\n", "Malformed XYZ header"),
        (
            "1\n-1.0\nH 0 0 0\n1\n-0.5\nC 0 0 0\n",
            "atom identity/order changes",
        ),
    ],
)
def test_xyz_reader_rejects_malformed_or_incomplete_frames(
    tmp_path: Path,
    contents: str,
    message: str,
) -> None:
    trajectory = tmp_path / "bad.xyz"
    trajectory.write_text(contents, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        read_xyz_trajectory(trajectory)
