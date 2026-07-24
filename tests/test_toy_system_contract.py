from __future__ import annotations

from pathlib import Path


TOY = Path(__file__).parents[1] / "examples" / "toy_system"


def _records(path: Path, *, bfactor: float | None = None):
    records = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if bfactor is not None and float(line[60:66]) != bfactor:
            continue
        records.append(
            (
                line[21:22],
                line[22:26],
                line[26:27],
                line[17:20],
                line[12:16],
            )
        )
    return records


def _coordinates(path: Path):
    return {
        (
            line[21:22],
            line[22:26],
            line[26:27],
            line[17:20],
            line[12:16],
        ): line[30:54]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    }


def test_toy_endpoints_reuse_one_ordered_ml_selection() -> None:
    model_r = _records(TOY / "ml_region_r.pdb")
    model_p = _records(TOY / "ml_region_p.pdb")
    layered_r = _records(TOY / "r_complex_layered.pdb", bfactor=0.0)
    layered_p = _records(TOY / "p_complex_layered.pdb", bfactor=0.0)

    assert model_r
    assert model_p == model_r
    assert layered_r == model_r
    assert layered_p == model_r
    full_r_coords = _coordinates(TOY / "r_complex.pdb")
    full_p_coords = _coordinates(TOY / "p_complex.pdb")
    model_r_coords = _coordinates(TOY / "ml_region_r.pdb")
    model_p_coords = _coordinates(TOY / "ml_region_p.pdb")
    assert all(
        model_r_coords[record] == full_r_coords[record]
        for record in model_r
    )
    assert all(
        model_p_coords[record] == full_p_coords[record]
        for record in model_p
    )


def test_toy_full_endpoint_atom_order_is_identical() -> None:
    assert _records(TOY / "r_complex.pdb") == _records(TOY / "p_complex.pdb")
    assert _records(TOY / "r_complex_layered.pdb") == _records(
        TOY / "p_complex_layered.pdb",
    )
