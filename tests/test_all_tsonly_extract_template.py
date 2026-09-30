"""TS-only ``mlmm all`` with pocket extraction (-c).

The TS is optimized on the full layered system, so every TS/IRC geometry that
is converted to PDB must use a full-system template, not the extracted pocket.
Extraction, layer assignment, TSOPT and IRC are replaced by lightweight fakes.
"""

from __future__ import annotations

import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from click.testing import CliRunner

from mlmm.cli import cli as root_cli
from mlmm.core import utils
from mlmm.core.result_commit import MLMM_RUN_ID_ENV
from mlmm.workflows import all as all_workflow
from mlmm.workflows import dft


class ReachedTemplate(BaseException):
    """Escape the pipeline (including ``except Exception``) once a template is used."""


_FULL = (
    ("HETATM", "C1", "MOL", 1, "C", (0.000, 0.000, 0.000)),
    ("HETATM", "O1", "MOL", 1, "O", (1.200, 0.000, 0.000)),
    ("HETATM", "O", "WAT", 2, "O", (4.000, 0.000, 0.000)),
    ("HETATM", "H1", "WAT", 2, "H", (4.900, 0.300, 0.000)),
)


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.next")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)
    return path


def _pdb(rows) -> str:
    return "".join(
        f"{rec:<6s}{serial:5d} {name:<4s} {resn:>3s} A{resi:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}          {element:>2s}\n"
        for serial, (rec, name, resn, resi, element, (x, y, z)) in enumerate(rows, 1)
    ) + "END\n"


def _n_atoms(pdb: Path) -> int:
    return sum(
        1 for line in Path(pdb).read_text().splitlines()
        if line.startswith(("ATOM", "HETATM"))
    )


def _fake_geom(energy: float) -> SimpleNamespace:
    coords = np.array([xyz for *_, xyz in _FULL], dtype=float)
    return SimpleNamespace(
        energy=energy,
        cart_coords=coords.reshape(-1),
        coords3d=coords,
        atoms=tuple(row[4] for row in _FULL),
    )


@pytest.mark.parametrize("continue_irc", [False, True], ids=["stop-before-irc", "irc"])
def test_tsonly_with_extraction_uses_full_system_template(tmp_path, monkeypatch, continue_irc):
    monkeypatch.delenv(MLMM_RUN_ID_ENV, raising=False)
    monkeypatch.setattr(utils, "_CONVERT_FILES_ENABLED", True)
    out = tmp_path / "out"
    full = _write(tmp_path / "complex.pdb", _pdb(_FULL))
    ml_only = _write(tmp_path / "ml_only.pdb", _pdb(_FULL[:2]))
    parm = _write(tmp_path / "unused.parm7", "template fixture; not a force field\n")

    model = Atoms("CO", positions=[row[5] for row in _FULL[:2]])
    workspace = SimpleNamespace(
        atoms_model=model, atoms_model_lh=model.copy(), model_pdb=ml_only,
        link_pairs=[], cleanup=lambda: None,
    )
    monkeypatch.setattr(dft, "_prepare_ml_region_workspace", lambda **kwargs: workspace)

    pockets = []

    def fake_extract(*, complex_pdb, output, **kwargs):
        for dest in output:
            _write(Path(dest), _pdb(_FULL[:2]))
            pockets.append(Path(dest).resolve())
        return {"charge_summary": {"total_charge": 0.0}}

    def fake_define_layers(*, input_pdb, output_pdb, model_pdb, **kwargs):
        _write(Path(output_pdb), Path(input_pdb).read_text())
        return {"ml_indices": [0, 1], "movable_mm_indices": [2, 3], "frozen_indices": []}

    ts_inputs = []

    def fake_tsopt(hei_pdb, *args, **kwargs):
        ts_inputs.append(Path(hei_pdb))
        ts_dir = Path(args[6])
        ts_pdb = _write(ts_dir / "final_geometry.pdb", Path(hei_pdb).read_text())
        ts_xyz = _write(ts_dir / "final_geometry.xyz", "4\n\n" + "".join(
            f"{el} {x} {y} {z}\n" for *_, el, (x, y, z) in _FULL
        ))
        g_ts = _fake_geom(-1.0)
        g_ts._tsopt_result = {"energy_hartree": -1.0, "n_imaginary_modes": 1}
        g_ts._tsopt_continuation = {
            "continue_irc": continue_irc,
            "reason": "fixture",
            "reaction_mode_index": 0,
            "n_imaginary_modes": 1,
        }
        g_ts._tsopt_result_path = None
        return ts_pdb, ts_xyz, g_ts

    irc_templates = []

    def fake_irc(**kwargs):
        irc_templates.append(Path(kwargs["seg_pocket_pdb"]))
        return {
            "left_min_geom": _fake_geom(-1.2),
            "right_min_geom": _fake_geom(-1.1),
            "ts_geom": _fake_geom(-1.0),
        }

    saved = []

    def fake_save(g, ref_pdb, out_dir, name):
        saved.append((name, Path(ref_pdb)))
        raise ReachedTemplate

    def no_calculator(*args, **kwargs):
        raise AssertionError("No calculator may be built in this fixture")

    monkeypatch.setattr(all_workflow, "extract_api", fake_extract)
    monkeypatch.setattr(all_workflow, "_define_layers", fake_define_layers)
    monkeypatch.setattr(all_workflow, "_run_tsopt_on_hei", fake_tsopt)
    monkeypatch.setattr(all_workflow, "_irc_and_match", fake_irc)
    monkeypatch.setattr(all_workflow, "_save_single_geom_for_tools", fake_save)
    monkeypatch.setattr(all_workflow, "_mlmm_calc", no_calculator)

    args = [
        "all", "-i", str(full), "-c", "MOL", "-q", "0", "-m", "1",
        "--parm", str(parm), "--model-pdb", str(ml_only),
        "--out-dir", str(out), "--tsopt", "--convert-files",
    ]
    with pytest.raises(ReachedTemplate):
        CliRunner().invoke(root_cli, args, catch_exceptions=False)

    assert len(pockets) == 1 and _n_atoms(pockets[0]) == 2
    ts_input, = ts_inputs
    assert _n_atoms(ts_input) == len(_FULL)

    name, template = saved[-1]
    assert name == ("reactant_irc" if continue_irc else "ts")
    assert template.resolve() == ts_input.resolve()
    assert _n_atoms(template) == len(_FULL)
    if continue_irc:
        irc_template, = irc_templates
        assert irc_template.resolve() == ts_input.resolve()
    else:
        assert irc_templates == []
