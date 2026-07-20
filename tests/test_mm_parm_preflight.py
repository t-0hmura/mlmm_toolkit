"""Preflight checks for AmberTools dependency detection in mm_parm."""

from __future__ import annotations

from mlmm.workflows import mm_parm


def _atom(serial: int, name: str, resname: str, resseq: int, x: float, y: float, z: float) -> str:
    elem = name.strip()[0]
    return (
        f"ATOM  {serial:5d} {name:^4s} {resname:>3s} A{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}\n"
    )


def test_missing_ambertools_commands_with_partial_paths() -> None:
    paths = {
        "tleap": "/usr/bin/tleap",
        "antechamber": None,
        "parmchk2": None,
    }
    missing = mm_parm.missing_ambertools_commands(paths)
    assert missing == ["antechamber", "parmchk2"]


def test_ambertools_available_uses_missing_command_detection(monkeypatch) -> None:
    def fake_which(cmd: str):
        if cmd == "antechamber":
            return None
        return f"/opt/amber/bin/{cmd}"

    monkeypatch.setattr(mm_parm, "which", fake_which)

    paths = mm_parm.ambertools_command_paths()
    assert paths["tleap"] == "/opt/amber/bin/tleap"
    assert paths["antechamber"] is None
    assert paths["parmchk2"] == "/opt/amber/bin/parmchk2"
    assert mm_parm.ambertools_available() is False


def test_add_ter_keeps_connected_peptide_block(tmp_path) -> None:
    src = tmp_path / "connected.pdb"
    dst = tmp_path / "connected_with_ter.pdb"
    src.write_text(
        "".join(
            [
                _atom(1, "N", "LEU", 1, 0.0, 0.0, 0.0),
                _atom(2, "CA", "LEU", 1, 1.0, 0.0, 0.0),
                _atom(3, "C", "LEU", 1, 2.0, 0.0, 0.0),
                _atom(4, "N", "MET", 2, 3.3, 0.0, 0.0),
                _atom(5, "CA", "MET", 2, 4.3, 0.0, 0.0),
            ]
        )
    )

    mm_parm.insert_ter_around_special_residues(src, dst, set())

    assert "TER\n" not in dst.read_text()


def test_add_ter_splits_distant_peptide_block(tmp_path) -> None:
    src = tmp_path / "broken.pdb"
    dst = tmp_path / "broken_with_ter.pdb"
    src.write_text(
        "".join(
            [
                _atom(1, "N", "LEU", 1, 0.0, 0.0, 0.0),
                _atom(2, "CA", "LEU", 1, 1.0, 0.0, 0.0),
                _atom(3, "C", "LEU", 1, 2.0, 0.0, 0.0),
                _atom(4, "N", "MET", 2, 40.0, 0.0, 0.0),
                _atom(5, "CA", "MET", 2, 41.0, 0.0, 0.0),
            ]
        )
    )

    mm_parm.insert_ter_around_special_residues(src, dst, set())

    out = dst.read_text().splitlines()
    assert out[2].startswith("ATOM")
    assert out[3] == "TER"
    assert out[4].startswith("ATOM")


def test_disulfide_cys_pair_is_renamed_to_cyx(tmp_path) -> None:
    """A CYS taking part in a detected S-S must become CYX before loadpdb.

    tleap's `bond` adds the S-S connection but does not strip the CYS template's
    HG, so a CYS-named disulfide cysteine would end up with a hypervalent SG.
    """
    src = tmp_path / "in.pdb"
    dst = tmp_path / "out.pdb"
    src.write_text(
        "".join(
            [
                _atom(1, "CB", "CYS", 10, 1.5, 0.0, 0.0),
                _atom(2, "SG", "CYS", 10, 2.5, 0.0, 0.0),
                _atom(3, "HG", "CYS", 10, 3.0, 0.9, 0.0),
                _atom(4, "CB", "CYS", 20, 5.5, 0.0, 0.0),
                _atom(5, "SG", "CYS", 20, 4.55, 0.0, 0.0),
            ]
        )
    )

    pairs = mm_parm.detect_disulfides_from_pdb(src, cutoff=mm_parm.DISULFIDE_CUTOFF)
    assert pairs == [(("A", 10), ("A", 20))]

    renamed = mm_parm.rename_disulfide_cys_to_cyx(src, dst, pairs)

    assert renamed == 2  # residues, not atom lines
    out = dst.read_text()
    assert "CYS" not in out
    assert out.count("CYX") == 5


def test_isolated_cys_and_disulfide_cym_are_left_alone(tmp_path) -> None:
    """Only bonded CYS is renamed; an isolated CYS and a CYM keep their names.

    CYM is the thiolate (formal -1); renaming it would silently change the net
    charge, so it is reported rather than converted.
    """
    src = tmp_path / "in.pdb"
    dst = tmp_path / "out.pdb"
    src.write_text(
        "".join(
            [
                _atom(1, "SG", "CYS", 30, 20.0, 0.0, 0.0),
                _atom(2, "SG", "CYM", 40, 40.0, 0.0, 0.0),
                _atom(3, "SG", "CYM", 41, 42.05, 0.0, 0.0),
            ]
        )
    )

    pairs = mm_parm.detect_disulfides_from_pdb(src, cutoff=mm_parm.DISULFIDE_CUTOFF)
    assert pairs == [(("A", 40), ("A", 41))]

    renamed = mm_parm.rename_disulfide_cys_to_cyx(src, dst, pairs)

    assert renamed == 0
    out = dst.read_text()
    assert "CYX" not in out
    assert out.count("CYM") == 2
    assert out.count("CYS") == 1


def test_rename_without_disulfides_is_a_noop(tmp_path) -> None:
    src = tmp_path / "in.pdb"
    dst = tmp_path / "out.pdb"
    body = "".join(
        [
            _atom(1, "SG", "CYS", 10, 0.0, 0.0, 0.0),
            _atom(2, "SG", "CYS", 20, 20.0, 0.0, 0.0),
        ]
    )
    src.write_text(body)

    renamed = mm_parm.rename_disulfide_cys_to_cyx(src, dst, [])

    assert renamed == 0
    assert dst.read_text() == body


def test_auto_disulfide_off_bonds_only_explicit_cyx(tmp_path) -> None:
    """--no-auto-disulfide trusts the input naming: only CYX pairs are bonded.

    A CYS that merely sits within the cutoff of another one is then left alone,
    both as a bond candidate and as a rename target.
    """
    src = tmp_path / "in.pdb"
    src.write_text(
        "".join(
            [
                _atom(1, "SG", "CYS", 10, 0.0, 0.0, 0.0),
                _atom(2, "SG", "CYS", 20, 2.05, 0.0, 0.0),
                _atom(3, "SG", "CYX", 30, 20.0, 0.0, 0.0),
                _atom(4, "SG", "CYX", 40, 22.05, 0.0, 0.0),
            ]
        )
    )

    auto = mm_parm.detect_disulfides_from_pdb(src, cutoff=mm_parm.DISULFIDE_CUTOFF)
    explicit = mm_parm.detect_disulfides_from_pdb(
        src, cutoff=mm_parm.DISULFIDE_CUTOFF, cyx_only=True
    )

    assert auto == [(("A", 10), ("A", 20)), (("A", 30), ("A", 40))]
    assert explicit == [(("A", 30), ("A", 40))]

    dst = tmp_path / "out.pdb"
    renamed = mm_parm.rename_disulfide_cys_to_cyx(src, dst, explicit)

    assert renamed == 0
    assert dst.read_text().count("CYS") == 2
