"""Children of ``all`` that consume the tsopt TS Hessian must see tsopt's topology PDB."""

import json
from pathlib import Path

import pytest


class _StopHere(RuntimeError):
    pass


class _EnergyOnly:
    def get_energy(self, _atoms, _coords):
        return {"energy": 0.0}

    def close(self):
        pass


def test_irc_and_ts_freq_receive_tsopt_topology_pdb(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    # The in-process TS Hessian cache matches evaluator.potential.input_pdb by
    # path, so a different --ref-pdb makes IRC and freq recompute the Hessian.
    from mlmm.workflows import all as all_workflow

    layered = tmp_path / "layered.pdb"
    layered.write_text("END\n", encoding="utf-8")
    hei_pdb = tmp_path / "path_opt" / "hei_seg_01.pdb"
    hei_pdb.parent.mkdir()
    hei_pdb.write_text("END\n", encoding="utf-8")
    hei_pdb.with_suffix(".xyz").write_text("1\nHEI\nHe 0 0 0\n", encoding="utf-8")
    argv: dict[str, list[str]] = {}

    class _Prepared:
        def __init__(self, path):
            self.geom_path = Path(path)
            self.source_path = Path(path)

        def cleanup(self):
            pass

    def _override(prepared, ref_pdb):
        prepared.source_path = Path(ref_pdb)

    def _fake_cli(name, _command, args, **_kwargs):
        argv[name] = list(args)
        if name != "tsopt":
            raise _StopHere
        ts_dir = Path(args[args.index("--out-dir") + 1])
        (ts_dir / "result.json").write_text(
            json.dumps({"optimization_status": "converged", "hessian_status": "completed",
                        "n_imaginary_modes": 1}),
            encoding="utf-8",
        )
        (ts_dir / "final_geometry.xyz").write_text("1\nTS\nHe 0 0 0\n", encoding="utf-8")
        (ts_dir / "final_geometry.pdb").write_text("END\n", encoding="utf-8")

    monkeypatch.setattr(all_workflow, "prepare_input_structure", _Prepared)
    monkeypatch.setattr(all_workflow, "apply_ref_pdb_override", _override)
    monkeypatch.setattr(all_workflow, "_mlmm_calc", lambda **_kwargs: _EnergyOnly())
    monkeypatch.setattr(all_workflow, "_run_cli_main", _fake_cli)
    template = all_workflow._ResolvedCalculatorTemplate.from_mapping({})
    seg_dir = tmp_path / "segments" / "seg_01"

    ts_pdb, ts_xyz, g_ts = all_workflow._run_tsopt_on_hei(
        hei_pdb, 0, 1, tmp_path / "system.parm7", tmp_path / "model.pdb", True,
        None, seg_dir, "hess", resolved_calc_template=template, ref_pdb=layered,
    )
    tsopt_ref = argv["tsopt"][argv["tsopt"].index("--ref-pdb") + 1]
    assert tsopt_ref == str(layered)

    with pytest.raises(_StopHere):
        all_workflow._irc_and_match(
            seg_idx=1, seg_dir=seg_dir, ref_pdb_for_seg=ts_pdb, seg_pocket_pdb=hei_pdb,
            g_ts=g_ts, q_int=0, spin=1, resolved_calc_template=template,
            ts_xyz_path=ts_xyz, real_parm7=tmp_path / "system.parm7",
            model_pdb=tmp_path / "model.pdb", detect_layer=True,
        )
    assert argv["irc"][argv["irc"].index("-i") + 1] == str(ts_xyz)
    assert argv["irc"][argv["irc"].index("--ref-pdb") + 1] == tsopt_ref

    # The TS freq call site in ``all`` forwards the same topology.
    x_ts = seg_dir / "structures" / "ts.xyz"
    x_ts.parent.mkdir()
    x_ts.write_text("1\nTS\nHe 0 0 0\n", encoding="utf-8")
    with pytest.raises(_StopHere):
        all_workflow._run_freq_for_state(
            seg_dir / "structures" / "ts.pdb", 0, 1, tmp_path / "system.parm7",
            tmp_path / "model.pdb", True, seg_dir / "freq" / "TS", None,
            xyz_path=x_ts, topology_pdb=getattr(g_ts, "_tsopt_topology_pdb", None),
        )
    assert argv["freq"][argv["freq"].index("--ref-pdb") + 1] == tsopt_ref

    endpoint_xyz = tmp_path / "reactant.xyz"
    endpoint_xyz.write_text("1\nendpoint\nHe 0 0 0\n", encoding="utf-8")
    with pytest.raises(_StopHere):
        all_workflow._run_opt_for_state(
            tmp_path / "reactant.pdb", 0, 1, tmp_path / "system.parm7",
            tmp_path / "model.pdb", True, tmp_path / "endpoint_opt", None, "hess",
            resolved_calc_template=template, xyz_path=endpoint_xyz,
            topology_pdb=getattr(g_ts, "_tsopt_topology_pdb", None),
        )
    assert argv["opt"][argv["opt"].index("--ref-pdb") + 1] == tsopt_ref
