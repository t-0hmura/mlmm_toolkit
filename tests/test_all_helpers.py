"""Unit tests for mlmm.workflows._all_helpers."""

from __future__ import annotations

import tempfile
from pathlib import Path

from mlmm.workflows._all_helpers import (
    AllContext,
    build_energy_level_dict,
    build_pipeline_summary_payload,
    copy_path_outputs_to_root,
    promote_diag_for_root,
)


def test_build_energy_level_dict_zero_referenced_kcal() -> None:
    # AU2KCALPERMOL is approximately 627.5095
    d = build_energy_level_dict(
        labels=["R", "TS", "P"],
        energies_au=[-100.0, -99.5, -100.2],
        ref_energy=-100.0,
        au_to_kcal=627.5095,
        diagram_path="/tmp/diag.png",
        structures={"R": "r.pdb", "TS": "ts.pdb", "P": "p.pdb"},
    )
    assert d["labels"] == ["R", "TS", "P"]
    assert d["energies_au"] == [-100.0, -99.5, -100.2]
    # 0, 0.5*627.5, -0.2*627.5
    assert d["energies_kcal"][0] == 0.0
    assert abs(d["energies_kcal"][1] - 313.75475) < 1e-3
    assert abs(d["energies_kcal"][2] - (-125.5019)) < 1e-3
    assert d["barrier_kcal"] == d["energies_kcal"][1]
    assert d["delta_kcal"] == d["energies_kcal"][-1]
    assert d["diagram"] == "/tmp/diag.png"
    assert d["structures"] == {"R": "r.pdb", "TS": "ts.pdb", "P": "p.pdb"}


def test_build_energy_level_dict_does_not_mutate_inputs() -> None:
    labels = ["R", "TS", "P"]
    energies = [-1.0, -0.5, -1.1]
    structs = {"R": "r.pdb"}
    d = build_energy_level_dict(
        labels=labels,
        energies_au=energies,
        ref_energy=-1.0,
        au_to_kcal=627.5,
        diagram_path="/tmp/x.png",
        structures=structs,
    )
    # caller's containers must not be aliased
    d["labels"].append("EXTRA")
    d["structures"]["NEW"] = "n.pdb"
    assert labels == ["R", "TS", "P"]
    assert "NEW" not in structs


def test_promote_diag_for_root_returns_none_on_empty() -> None:
    assert promote_diag_for_root(None, "stem", Path(".")) is None
    assert promote_diag_for_root({}, "stem", Path(".")) is None


def test_promote_diag_for_root_rewrites_name_and_image() -> None:
    original = {"name": "MEP", "image": "old.png", "x": [0, 1]}
    promoted = promote_diag_for_root(original, "energy_diagram_UMA", Path("/tmp/out"))
    assert promoted is not None
    assert promoted["name"] == "energy_diagram_UMA_all"
    assert promoted["image"] == "/tmp/out/energy_diagram_UMA_all.png"
    # caller's dict must not be mutated
    assert original["name"] == "MEP"
    assert original["image"] == "old.png"
    # other keys pass through
    assert promoted["x"] == [0, 1]


def test_copy_path_outputs_skips_missing_files(tmp_path: Path) -> None:
    src = tmp_path / "src"
    dst = tmp_path / "dst"
    src.mkdir()
    dst.mkdir()
    # No files present — best-effort no-op
    warnings: list[str] = []
    copy_path_outputs_to_root(src, dst, warn_fn=warnings.append)
    assert warnings == []
    assert list(dst.iterdir()) == []


def test_copy_path_outputs_copies_known_artefacts(tmp_path: Path) -> None:
    src = tmp_path / "src"
    dst = tmp_path / "dst"
    src.mkdir()
    dst.mkdir()
    (src / "mep_plot.png").write_text("dummy-png")
    (src / "mep.pdb").write_text("dummy-pdb")
    (src / "summary.json").write_text("{}")
    (src / "mep_trj.xyz").write_text("xyz")
    (src / "unrelated.tmp").write_text("ignore me")

    copy_path_outputs_to_root(src, dst)

    assert (dst / "mep_plot.png").read_text() == "dummy-png"
    assert (dst / "mep.pdb").read_text() == "dummy-pdb"
    assert (dst / "summary.json").read_text() == "{}"
    assert (dst / "mep_trj.xyz").read_text() == "xyz"
    assert not (dst / "unrelated.tmp").exists()


def test_build_pipeline_summary_payload_shape() -> None:
    with tempfile.TemporaryDirectory() as d:
        out_dir = Path(d) / "out"
        path_dir = Path(d) / "path"
        out_dir.mkdir()
        path_dir.mkdir()
        summary = {
            "n_images": 5,
            "n_segments": 1,
            "segments": [{"kind": "seg", "bond_changes": "A->B"}],
            "energy_diagrams": [{"name": "MEP", "x": [0, 1]}],
        }
        payload = build_pipeline_summary_payload(
            out_dir=out_dir,
            path_dir=path_dir,
            summary=summary,
            refine_path=True,
            thresh="gau_loose",
            thresh_post="gau",
            flatten=False,
            do_tsopt=True,
            do_thermo=False,
            do_dft=False,
            opt_mode_norm="grad",
            opt_mode_post="HESS",
            command_str="mlmm all -i foo.pdb",
            q_int=-1,
            spin=1,
            post_segment_logs=[{"seg": 1, "status": "ok"}],
        )
    assert payload["pipeline_mode"] == "path-search"
    assert payload["refine_path"] is True
    assert payload["opt_mode"] == "grad"
    assert payload["opt_mode_post"] == "hess"
    assert payload["charge"] == -1
    assert payload["spin"] == 1
    assert payload["mep"]["n_images"] == 5
    assert payload["mep"]["diagram"]["name"] == "MEP"
    assert payload["post_segments"] == [{"seg": 1, "status": "ok"}]
    assert payload["energy_diagrams"] == summary["energy_diagrams"]


def test_all_context_fields_match_cli_signature() -> None:
    """Regression guard: AllContext field names must mirror cli() params.

    When someone adds / removes a parameter on `mlmm.workflows.all.cli`
    they must also update `AllContext`; this test enforces the contract
    so the dataclass cannot silently drift out of sync. Excludes
    `ctx: click.Context` since it isn't carried in the bundle.
    """
    import dataclasses
    import inspect
    from mlmm.workflows.all import cli
    # cli is a Click command; the underlying callable is .callback
    callback = cli.callback
    assert callback is not None, "cli() should expose its callback"
    sig = inspect.signature(callback)
    cli_param_names = {n for n in sig.parameters if n != "ctx"}
    ctx_field_names = {f.name for f in dataclasses.fields(AllContext)}
    missing = cli_param_names - ctx_field_names
    extra = ctx_field_names - cli_param_names
    assert not missing, (
        f"AllContext is missing fields present in cli(): {sorted(missing)}"
    )
    assert not extra, (
        f"AllContext has fields NOT in cli(): {sorted(extra)}"
    )


def test_all_context_frozen_and_field_count() -> None:
    # AllContext mirrors the cli() signature; we don't need every field
    # populated to verify the contract here, but the dataclass should
    # be frozen and provide the canonical bundle for future decomp.
    ctx = AllContext(
        input_paths=(),
        center_spec=None,
        out_dir=Path("."),
        radius=5.0,
        radius_het2het=0.0,
        include_h2o=True,
        exclude_backbone=False,
        add_linkh=False,
        selected_resn="",
        modified_residue="",
        ligand_charge=None,
        charge_override=None,
        parm7_override=None,
        model_pdb_override=None,
        mm_ff_set="ff19SB",
        mm_add_ter=True,
        mm_keep_temp=False,
        mm_ligand_mult=None,
        spin=1,
        max_nodes=5,
        max_cycles=100,
        climb=True,
        opt_mode="grad",
        opt_mode_post=None,
        dump=False,
        refine_path=True,
        thresh=None,
        thresh_post="gau",
        config_yaml=None,
        show_config=False,
        dry_run=False,
        pre_opt=True,
        hessian_calc_mode=None,
        detect_layer=True,
        do_tsopt=False,
        do_thermo=False,
        do_dft=False,
        scan_lists_raw=(),
        scan_out_dir=None,
        scan_one_based=None,
        scan_max_step_size=None,
        scan_bias_k=None,
        scan_relax_max_cycles=None,
        scan_preopt_override=None,
        scan_endopt_override=None,
        convert_files=True,
        ref_pdb_cli=None,
        backend=None,
        embedcharge=False,
        embedcharge_cutoff=None,
        link_atom_method=None,
        mm_backend=None,
        use_cmap=None,
        tsopt_max_cycles=None,
        flatten=False,
        skip_final_freq=False,
        tsopt_out_dir=None,
        freq_out_dir=None,
        freq_max_write=None,
        freq_amplitude_ang=None,
        freq_n_frames=None,
        freq_sort=None,
        freq_temperature=None,
        freq_pressure=None,
        dft_out_dir=None,
        dft_func_basis=None,
        dft_max_cycle=None,
        dft_conv_tol=None,
        dft_grid_level=None,
        dft_engine=None,
        cli_coord_type=None,
        precision=None,
        backend_model=None,
    )
    import dataclasses
    assert dataclasses.is_dataclass(ctx)
    assert getattr(ctx.__dataclass_params__, "frozen", False) is True  # frozen
    field_names = {f.name for f in dataclasses.fields(ctx)}
    assert "input_paths" in field_names
    assert "do_dft" in field_names
    assert "precision" in field_names
    assert len(field_names) == 75  # current cli() param count
