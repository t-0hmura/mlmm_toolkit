"""Flags and settings that do nothing must say so.

Requested settings that cannot be honored must emit a diagnostic instead of
appearing to have been applied. Output-only guards do not change behavior,
thresholds, or defaults.
"""

from pathlib import Path

import numpy as np
import pytest

import mlmm.workflows.freq as freq_mod
import mlmm.workflows.mm_parm as mm_parm_mod
import mlmm.workflows.tsopt as tsopt_mod


def test_unresolvable_layer_sets_are_reported(monkeypatch, capsys) -> None:
    # Empty layer sets collapse every --active-dof-mode to ALL atoms, and the
    # only previous signal was a debug-level log line. The warning must not
    # change what callers receive.
    def failing_calculator(**_kwargs):
        raise RuntimeError("no topology")

    monkeypatch.setattr(freq_mod, "mlmm", failing_calculator)

    layer_sets = freq_mod._collect_layer_atom_sets({})

    assert layer_sets == {"ml": set(), "hess_mm": set(), "movable_mm": set(), "frozen_mm": set()}
    err = capsys.readouterr().err
    assert "could not resolve ML/MM layer sets (no topology)" in err
    assert "falls back to ALL atoms regardless of --active-dof-mode" in err


@pytest.mark.parametrize("module", [freq_mod, tsopt_mod], ids=["freq", "tsopt"])
def test_mode_pdb_fallback_is_reported(module, tmp_path, monkeypatch, capsys) -> None:
    # The file still appears, but without the reference topology its atom names
    # and residues are not the input's.
    from pysisyphus.Geometry import Geometry

    def failing_converter(*_args, **_kwargs):
        raise RuntimeError("template mismatch")

    monkeypatch.setattr(module, "is_convert_file_enabled", lambda: True)
    monkeypatch.setattr(module, "convert_xyz_to_pdb", failing_converter)
    geometry = Geometry(("H", "H"), np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.4]), coord_type="cart")
    reference = tmp_path / "input.pdb"
    reference.write_text("END\n", encoding="utf-8")
    out_pdb = tmp_path / "mode.pdb"

    module._write_mode_trj_and_pdb(
        geometry, np.array([1.0, 0.0, 0.0, -1.0, 0.0, 0.0]), tmp_path / "mode_trj.xyz",
        out_pdb, n_frames=2, ref_pdb=reference,
    )

    assert out_pdb.is_file()
    assert (
        "mode PDB fell back to plain ASE output without the reference topology: "
        "template mismatch"
    ) in capsys.readouterr().err


def test_unused_ligand_charge_entries_are_reported(tmp_path, monkeypatch, capsys) -> None:
    # A --ligand-charge entry for a residue tleap already knows is never read;
    # a mistyped residue name is exactly how the wrong net charge gets used.
    source = tmp_path / "input.pdb"
    source.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )

    def tleap_knows_everything(cmd, cwd=None, logfile=None):
        for name in mm_parm_mod._TLEAP_REQUIRED_OUTPUTS:
            (Path(cwd) / name).write_text("fixture\n", encoding="utf-8")
        return 0

    monkeypatch.setattr(mm_parm_mod, "run", tleap_knows_everything)
    monkeypatch.setattr(mm_parm_mod, "parse_tleap_unknown_residues", lambda _log: set())

    mm_parm_mod.ambertools_route(
        source, str(tmp_path / "out"), {"ZZZ": 1}, {}, False, tmp_path, "ff19SB", False,
    )

    assert (
        "--ligand-charge/--ligand-mult entry ZZZ matched no residue needing parameters"
    ) in capsys.readouterr().err


@pytest.mark.parametrize("dry_run", [True, False], ids=["dry-run", "run"])
@pytest.mark.parametrize(("calc_yaml", "expected"), [("", True), ("calc:\n  use_cmap: false\n", False)],
                         ids=["default", "off"])
def test_dft_passes_the_cmap_switch_to_the_model_workspace(
    tmp_path, monkeypatch, dry_run, calc_yaml, expected
) -> None:
    """The CMAP policy must reach the ordinary ML/MM core used by dft.

    `dft.py` used to build its own model workspace and never read the setting,
    so with a backbone-containing ML region the ML/MM and QM/MM paths computed
    E_model_low on different models -- and `result.yaml` recorded
    `use_cmap: False` either way. Both builds of the workspace (the dry-run
    validation and the run) have to receive it.
    """
    from click.testing import CliRunner

    from mlmm.cli import cli as root_cli
    import mlmm.workflows.dft as dft_mod

    toy = Path(__file__).resolve().parents[1] / "examples" / "toy_system"
    if not (toy / "p_toy.parm7").exists():
        pytest.skip("toy_system example inputs not present")

    class WorkspaceRequested(Exception):
        pass

    requested = {}

    def prepare(**kwargs):
        requested.update(kwargs)
        raise WorkspaceRequested

    monkeypatch.setattr(dft_mod, "_prepare_ml_region_workspace", prepare)
    config = tmp_path / "config.yaml"
    config.write_text(calc_yaml, encoding="utf-8")

    CliRunner().invoke(root_cli, [
        "dft", "-i", str(toy / "p_toy.pdb"), "--parm", str(toy / "p_toy.parm7"),
        "-q", "-1", "-m", "1", "--config", str(config), "-o", str(tmp_path / "out"),
        *(["--dry-run"] if dry_run else []),
    ])

    assert requested["use_cmap"] is expected


def test_dft_cmap_default_matches_the_product_default() -> None:
    from mlmm.core.defaults import MLMM_CALC_KW

    assert MLMM_CALC_KW["use_cmap"] is True
