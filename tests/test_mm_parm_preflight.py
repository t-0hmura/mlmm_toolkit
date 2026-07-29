"""mm_parm preflight: AmberTools dependency detection plus the surrounding
input/topology checks."""

from __future__ import annotations

import stat

import pytest

from mlmm.workflows import mm_parm


def _atom(
    serial: int,
    name: str,
    resname: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    icode: str = "",
) -> str:
    elem = name.strip()[0]
    return (
        f"ATOM  {serial:5d} {name:^4s} {resname:>3s} A{resseq:4d}{icode:1s}   "
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


def test_run_executes_the_resolved_environment_path(monkeypatch) -> None:
    seen = []

    class FakePopen:
        def __init__(self, argv, **_kwargs):
            seen.append(argv)
            self.stdout = []

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        @staticmethod
        def wait():
            return 0

    monkeypatch.setattr(
        mm_parm, "which", lambda cmd: f"/env/bin/{cmd}",
    )
    monkeypatch.setattr(mm_parm.subprocess, "Popen", FakePopen)

    assert mm_parm.run(["tleap", "-f", "input.in"]) == 0
    assert seen == [["/env/bin/tleap", "-f", "input.in"]]


def test_second_tleap_pass_cannot_reuse_first_pass_outputs(
    tmp_path, monkeypatch
) -> None:
    source = tmp_path / "input.pdb"
    source.write_text(
        _atom(1, "C1", "LIG", 1, 0.0, 0.0, 0.0),
        encoding="utf-8",
    )
    calls = []

    def fake_run(cmd, cwd=None, logfile=None):
        calls.append(list(cmd))
        work = tmp_path if cwd is None else cwd
        if len(calls) == 1:
            for name in mm_parm._TLEAP_COMPLEX_OUTPUTS:
                (work / name).write_text("pass-1", encoding="utf-8")
            if logfile is not None:
                logfile.write_text("unknown LIG", encoding="utf-8")
            return 0
        assert all(not (work / name).exists() for name in mm_parm._TLEAP_COMPLEX_OUTPUTS)
        if logfile is not None:
            logfile.write_text("pass 2 failed", encoding="utf-8")
        return 1

    monkeypatch.setattr(mm_parm, "run", fake_run)
    monkeypatch.setattr(
        mm_parm, "parse_tleap_unknown_residues", lambda _path: {"LIG"}
    )
    monkeypatch.setattr(mm_parm, "extract_first_residue_pdb", lambda *_args: True)
    monkeypatch.setattr(
        mm_parm,
        "antechamber_parametrize",
        lambda *_args: (tmp_path / "LIG.mol2", tmp_path / "LIG.frcmod"),
    )

    with pytest.raises(RuntimeError, match="pass 2 exited with code 1"):
        mm_parm.ambertools_route(
            source,
            str(tmp_path / "out"),
            {},
            {},
            False,
            tmp_path,
            "ff19SB",
            False,
        )

    assert len(calls) == 2
    assert not (tmp_path / "out.parm7").exists()
    assert not (tmp_path / "out.rst7").exists()


def test_leap_pdb_export_populates_elements_without_reordering(tmp_path) -> None:
    source = tmp_path / "leap.pdb"
    destination = tmp_path / "system.pdb"
    source.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00\n"
        "HETATM    2 MG    MG A   2       2.000   0.000   0.000  1.00  0.00\n"
        "TER\nEND\n",
        encoding="utf-8",
    )

    assigned, unresolved = mm_parm.copy_pdb_with_element_fields(source, destination)

    lines = destination.read_text(encoding="utf-8").splitlines()
    assert (assigned, unresolved) == (2, 0)
    assert [line[:6].strip() for line in lines] == ["ATOM", "HETATM", "TER", "END"]
    assert lines[0][76:78].strip() == "C"
    assert lines[1][76:78].strip() == "Mg"


def test_leap_pdb_export_preserves_line_endings_and_mode(tmp_path) -> None:
    source = tmp_path / "leap_mixed_endings.pdb"
    destination = tmp_path / "system.pdb"
    source.write_bytes(
        b"ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00\r\n"
        b"HETATM    2 MG    MG A   2       2.000   0.000   0.000  1.00  0.00\r"
        b"TER\nEND"
    )
    source.chmod(0o640)

    assigned, unresolved = mm_parm.copy_pdb_with_element_fields(source, destination)

    payload = destination.read_bytes()
    assert (assigned, unresolved) == (2, 0)
    assert payload.count(b"\r\n") == 1
    assert payload.count(b"\r") == 2
    assert payload.count(b"\n") == 2
    assert payload.endswith(b"END")
    assert stat.S_IMODE(destination.stat().st_mode) == 0o640


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
    assert pairs == [(("A", 10, ""), ("A", 20, ""))]

    renamed = mm_parm.rename_disulfide_cys_to_cyx(src, dst, pairs)

    assert renamed == 2  # residues, not atom lines
    out = dst.read_text()
    assert "CYS" not in out
    assert out.count("CYX") == 5


def test_auto_disulfide_excludes_cym_thiolate_pairs(tmp_path) -> None:
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
    assert pairs == []

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

    assert auto == [
        (("A", 10, ""), ("A", 20, "")),
        (("A", 30, ""), ("A", 40, "")),
    ]
    assert explicit == [(("A", 30, ""), ("A", 40, ""))]

    dst = tmp_path / "out.pdb"
    renamed = mm_parm.rename_disulfide_cys_to_cyx(src, dst, explicit)

    assert renamed == 0
    assert dst.read_text().count("CYS") == 2


def test_disulfide_identity_retains_insertion_code(tmp_path) -> None:
    source = tmp_path / "insertions.pdb"
    source.write_text(
        "".join(
            [
                _atom(1, "SG", "CYS", 10, 20.0, 0.0, 0.0),
                _atom(2, "SG", "CYS", 10, 0.0, 0.0, 0.0, "A"),
                _atom(3, "SG", "CYS", 20, 2.05, 0.0, 0.0),
            ]
        ),
        encoding="utf-8",
    )

    pairs = mm_parm.detect_disulfides_from_pdb(source)
    assert pairs == [(("A", 10, "A"), ("A", 20, ""))]

    mapping = mm_parm.build_leap_residue_index(source)
    assert mapping[("A", "  10", "")] == 1
    assert mapping[("A", "  10", "A")] == 2
    assert mapping[("A", "  20", "")] == 3


def test_failed_hydrogenated_fallback_preserves_existing_canonical_pdb(
    tmp_path, monkeypatch
) -> None:
    source = tmp_path / "input.pdb"
    source.write_text(_atom(1, "C1", "LIG", 1, 0.0, 0.0, 0.0))
    prefix = tmp_path / "system"
    canonical = prefix.with_suffix(".pdb")
    canonical.write_bytes(b"existing canonical generation\n")

    monkeypatch.setattr(
        mm_parm,
        "ambertools_command_paths",
        lambda: {
            name: f"/bin/{name}"
            for name in mm_parm._AMBERTOOLS_REQUIRED_COMMANDS
        },
    )

    def add_h(_source, destination, _ph):
        destination.write_bytes(b"partial hydrogenated fallback\n")

    monkeypatch.setattr(mm_parm, "add_hydrogens_with_pdbfixer", add_h)
    monkeypatch.setattr(
        mm_parm,
        "ambertools_route",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            RuntimeError("injected")
        ),
    )
    args = mm_parm.Args(
        pdb=source,
        out_prefix=str(prefix),
        ligand_charge={},
        ligand_mult={},
        keep_temp=False,
        add_ter=False,
        auto_disulfide=True,
        add_h=True,
        ph=7.0,
        ff_set="ff19SB",
        out_prefix_given=True,
    )

    with pytest.raises(RuntimeError, match="injected"):
        mm_parm.run_pipeline(args)

    assert canonical.read_bytes() == b"existing canonical generation\n"
