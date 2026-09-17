"""Plain XYZ parsing preserves ASE semantics without reading future frames."""
import gzip
from io import StringIO

import pytest
from ase.io import read as ase_read

from mlmm.core.utils import _iter_plain_xyz_atoms, convert_xyz_to_pdb
from mlmm.io.structure_formats import render_pdb_coordinate_frames

XYZ = "2\n arbitrary unclosed quote ' Properties=ignored\nc 0.123456 0 0 extra\nZN 1.345678 0 0\n"


@pytest.mark.parametrize("compressed", [False, True])
def test_framewise_overlay_matches_plain_xyz_semantics(tmp_path, compressed):
    ref = tmp_path / "ref.pdb"
    ref.write_text(
        "HEADER    retained metadata\n"
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  0.75 18.00           C  \n"
        "TER       2      ALA A   1\n"
        "HETATM    3 ZN    ZN B   7       1.200   0.000   0.000  0.50 21.00          ZN2+\n"
        "CONECT    1    3\nEND\n"
    )
    xyz = tmp_path / ("frames.xyz.gz" if compressed else "frames.xyz")
    content = (XYZ + XYZ.replace("0.123456", "0.654321")).encode()
    xyz.write_bytes(gzip.compress(content) if compressed else content)
    old_frames = ase_read(StringIO(content.decode()), index=":", format="xyz")
    expected = render_pdb_coordinate_frames(
        ref, [f.get_chemical_symbols() for f in old_frames],
        [f.get_positions() for f in old_frames],
    )
    out = tmp_path / "out.pdb"
    convert_xyz_to_pdb(xyz, ref, out)
    assert out.read_text() == expected


def test_reader_defers_late_invalid_frame(tmp_path):
    xyz = tmp_path / "frames.xyz"
    xyz.write_text(XYZ + "2\ntruncated\nC 0 0 0\n")
    frames = _iter_plain_xyz_atoms(xyz)
    assert next(frames).get_chemical_symbols() == ["C", "Zn"]
    with pytest.raises(ValueError, match="Incomplete XYZ frame 2"):
        next(frames)
