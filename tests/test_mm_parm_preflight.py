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
