"""Preflight checks for AmberTools dependency detection in mm_parm."""

from __future__ import annotations

from mlmm_toolkit import mm_parm


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
