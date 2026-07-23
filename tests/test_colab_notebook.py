"""Static contracts for the release-matched Colab GUI notebook."""

from __future__ import annotations

import ast
import csv
import datetime
import glob
import html
import json
import os
import shlex
import types
from pathlib import Path

import pytest


NOTEBOOK = Path(__file__).parents[1] / "examples" / "mlmm_colab.ipynb"


def _execute_app(monkeypatch, tmp_path: Path) -> tuple[dict, list]:
    """Execute the complete app cell with real widgets and captured HTML."""
    rendered: list = []
    import IPython.display as ipd
    monkeypatch.setattr(ipd, "display", lambda *args, **kwargs: rendered.extend(args))
    monkeypatch.setattr(ipd, "clear_output", lambda *args, **kwargs: None)
    monkeypatch.setattr(ipd, "HTML", lambda value: value)
    monkeypatch.setattr(ipd, "Image", lambda *args, **kwargs: (args, kwargs))
    monkeypatch.chdir(tmp_path)
    namespace = {"TOOL": "mlmm", "BACKEND": "mace", "REPO_DIR": "unused"}
    source = _notebook()["cells"][2]["source"]
    exec(compile(source, str(NOTEBOOK), "exec"), namespace)
    return namespace, rendered


def _notebook() -> dict:
    return json.loads(NOTEBOOK.read_text(encoding="utf-8"))


def _output_contract() -> dict:
    """Execute only the pure output-tracking helpers from the app cell."""
    source = _notebook()["cells"][2]["source"]
    wanted = {
        "_normalized_scope_argv", "_effective_out_dir",
        "_click_subcommand_params", "_click_output_default",
        "_trj2fig_output_targets", "_flag_enabled", "_force_dry_run",
        "_grouped_option_values", "_parsed_path",
        "_input_needs_cif_companion", "_exact_output_scope",
        "_output_scope", "_effective_result_root", "_matches_output_scope",
        "_snapshot_files", "_snapshot_output_scope",
        "_output_scope_collision", "_structured_current_paths",
    }
    tree = ast.parse(source)
    module = ast.Module(
        body=[node for node in tree.body
              if isinstance(node, ast.FunctionDef) and node.name in wanted],
        type_ignores=[],
    )
    namespace = {
        "glob": glob, "json": json, "os": os, "shlex": shlex,
        "Path": Path, "CLI": "mlmm",
        "S": {"out_dir": "result", "subcmd": "all"},
    }
    exec(compile(module, str(NOTEBOOK), "exec"), namespace)
    assert wanted <= namespace.keys()
    return namespace


def _viewer_contract() -> dict:
    """Execute pure selection/highlight and trajectory-label helpers."""
    source = _notebook()["cells"][2]["source"]
    wanted = {
        "_pick_text",
        "_resolve_click_meta",
        "_residue_id_selector",
        "_rich_residue_selector",
        "_center_cli_selectors",
        "_input_count_error",
        "_artifact_kind",
        "_csv_preview_html",
        "_text_preview_html",
        "_atom_signatures",
        "_resolve_atom_query",
        "_trajectory_semantics",
        "_stationary",
    }
    tree = ast.parse(source)
    module = ast.Module(
        body=[node for node in tree.body
              if isinstance(node, ast.FunctionDef) and node.name in wanted],
        type_ignores=[],
    )
    namespace = {
        "os": os, "S": {}, "_TRAJ": {}, "Path": Path, "csv": csv, "html": html,
        "json": json, "_TEXT_PREVIEW_LIMIT": 512 * 1024,
        "SPEC": {}, "_ARTIFACT_KINDS": {
            ".png": "image", ".jpg": "image", ".jpeg": "image", ".svg": "SVG",
            ".html": "interactive HTML", ".csv": "CSV table", ".pdf": "PDF",
            ".pdb": "structure", ".ent": "structure", ".cif": "structure",
            ".mmcif": "structure",
        },
    }
    exec(compile(module, str(NOTEBOOK), "exec"), namespace)
    assert wanted <= namespace.keys()
    return namespace


def test_colab_notebook_has_valid_code_cells_and_gpu_metadata() -> None:
    notebook = _notebook()
    introduction = notebook["cells"][0]["source"]

    assert notebook["nbformat"] == 4
    assert notebook["metadata"]["accelerator"] == "GPU"
    assert len(notebook["cells"]) == 4
    assert "[GitHub](https://github.com/t-0hmura/mlmm_toolkit)" in introduction
    assert "ML/MM reaction paths on the full solvated protein" in introduction
    assert "**1 Input → 2 Viewer → 3 Options → 4 Results**" in introduction
    assert (
        "[ChemRxiv](https://chemrxiv.org/doi/full/"
        "10.26434/chemrxiv-2025-jft1k)"
    ) in introduction
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            ast.parse(cell["source"])


def test_colab_setup_is_pinned_to_matching_release_and_one_backend() -> None:
    setup = _notebook()["cells"][1]["source"]

    assert 'mlmm_ref = "v0.3.3"' in setup
    # The release notebook installs the pinned wheel from PyPI, the same way a
    # normal user does, so the version guard below compares the version actually
    # resolved by pip against the requested tag.
    assert "pip('mlmm-toolkit' + ('[dft]' if install_dft else '') + '==' + mlmm_ref.lstrip('v'))" in setup
    assert "git clone" not in setup
    assert "installed_version != mlmm_ref[1:]" in setup
    assert "anywidget" not in setup
    assert "version('mlmm-toolkit')" in setup
    assert "Restart the Colab runtime first" in setup
    assert "mace-torch>=0.3.8" in setup
    assert "HF_TOKEN" in setup
    assert "install_dft = False" in setup
    assert "INSTALL_DFT = install_dft" in setup
    assert "_dft_packages = {'pyscf': 'pyscf', 'gpu4pyscf': 'gpu4pyscf-cuda12x'}" in setup
    assert "importlib.util.find_spec(module)" in setup
    assert "_dft_imports = ('pyscf', 'basis_set_exchange', 'gpu4pyscf.dft')" in setup
    assert "for _module in _dft_imports: importlib.import_module(_module)" in setup
    assert "_cupy.cuda.runtime.getDeviceCount()" in setup
    assert "DFT packages installed but failed their import/GPU check" in setup
    assert "DFT support installed: PySCF %s · GPU4PySCF %s" in setup
    assert "py3Dmol" not in setup
    assert "pip('ipywidgets','matplotlib')" in setup
    assert "[%%d/5] %%s".replace("%%", "%") in setup
    assert "time.monotonic()" in setup
    assert ".pysisyphusrc" in setup
    assert setup.index(".pysisyphusrc") < setup.index("pip('mlmm-toolkit'")
    assert "nglview" not in setup
    assert "first run ~5 min" in setup
    assert "first run ~5-10 min" not in setup


def test_colab_gui_is_mlmm_native_and_tracks_structure_contracts() -> None:
    app = _notebook()["cells"][2]["source"]
    whole = NOTEBOOK.read_text(encoding="utf-8").lower()

    assert isinstance(app, str)
    assert "pdb2reaction" not in whole
    assert "p2r" not in whole
    assert ".pdb,.ent,.cif,.mmcif,.parm7" in app
    assert "'.pdb,.ent,.cif,.mmcif,.parm7,.xyz,.gjf,.com,.inp,.csv'" in app
    assert "prepare_input_structure" in app
    assert "load_pdb_atom_metadata" in app
    assert "Path(_runtime_path('viewer_input.pdb'))" in app
    assert "callback_ns = 'mlmm_gui'" in app
    assert "_co.register_callback('mlmm_gui.on_click', on_click)" in app
    assert "target.invokeFunction(cfg.callback+'.'+suffix,args,{})" in app
    assert "CHAIN:RESNAME:RESSEQ" in app
    assert "String(item.sourceIndex)" in app
    assert "ML/MM compute commands require a matching Amber parm7" in app
    assert "elif parm:" in app
    assert "S['parm'] = parm" in app


def test_colab_gui_uses_public_cli_and_backend_neutral_results() -> None:
    app = _notebook()["cells"][2]["source"]

    assert "cmd += ['--backend-model', mdl]" in app
    assert "model_config.yaml" not in app
    assert "cmd.append('--tsopt')" in app
    assert "cmd.append('--thermo')" in app
    assert "[n, 'True']" not in app
    assert "ps.get('mlip')" in app
    assert "ps.get('gibbs_mlip')" in app
    assert "ps.get('uma')" not in app
    assert "ps.get('gibbs_uma')" not in app


def test_colab_gui_guards_and_reads_the_edited_output_path() -> None:
    app = _notebook()["cells"][2]["source"]

    assert "def _effective_out_dir(argv=None):" in app
    assert "flags = ('-o', '--out', '--out-dir', '--output', '--out-prefix')" in app
    assert "def _trj2fig_output_targets(argv):" in app
    assert "command.make_parser(sub_ctx).parse_args" in app
    assert "if sub == 'mm-parm':" in app
    assert "if sub == 'oniom-import':" in app
    assert "'extract', 'define-layer', 'fix-altloc'" in app
    assert "def _matches_output_scope(path, scope):" in app
    assert "def _snapshot_output_scope(scope):" in app
    assert "def _output_scope_collision(scope):" in app
    assert "scope = _output_scope(a)" in app
    assert "target = scope['target']; out = scope['root']" in app
    assert "out = scope['root']" in app
    assert "os.path.isdir(out)" in app
    assert "_results(out)" in app
    assert "reuse non-empty out dir" in app


def test_colab_gui_keeps_responsive_release_layout() -> None:
    app = _notebook()["cells"][2]["source"]

    assert "@media (max-width: 600px)" in app
    assert "max_width='100%'" in app
    assert "flex:0 0 auto; min-width:0" in app
    assert "layout=W.Layout(width='260px', max_width='100%')" in app
    assert "_MOLSTAR_VERSION = '5.6.1'" in app
    assert "molstar@%s/build/viewer/molstar.js" in app
    assert "layoutShowSequence:cfg.showSequence" in app
    assert "layoutShowControls:true" in app
    assert "collapseRightPanel:true" in app
    assert "viewportFocusBehavior" not in app
    assert "Mol* owns representation, colour, camera, measurement, and screenshot controls." in app
    assert "view_controls = W.HBox([" in app
    assert "cb_water," in app
    assert "last_pick_row = W.HBox([last_pick_html]," in app
    assert "btn_clear_pick" not in app
    assert "sel_lang" not in app
    assert "No structure loaded" in app
    assert "('Scan 1', 'scan')" in app
    assert "style={'button_width': '74px', 'description_width': '0px'}" in app
    assert "No results yet." in app
    assert "frame_slider.disabled = False" in app
    assert "trajectory_box.layout.display = 'none'" in app
    # Colab renders ipywidgets' Tab and Accordion as empty blocks. The tab
    # buttons keep all panes mounted, preserving upload queues and WebGL state.
    assert "_TAB_PAGES = [('1 Input', input_box), ('2 Viewer', viewer_box)," in app
    assert "('3 Options', options_box), ('4 Results', results_box)]" in app
    assert "_tab_body = W.VBox([page for _label, page in _TAB_PAGES])" in app
    assert "_tab_body.children = [_TAB_PAGES[i][1]]" not in app
    assert "_pane.layout.display = '' if _j == i else 'none'" in app
    assert "Options (optional)" not in app
    assert "def _tab_go(i):" in app
    assert "W.Tab(" not in app
    assert "def _collapsible(title, child, on_open=None):" in app
    assert "W.Accordion(" not in app
    # Browser-native details/summary opens without a Python callback.
    assert "def _info_markup(tip, revision=0):" in app
    assert "def _info_control(tip, target=None):" in app
    assert "def _set_info_text(control, tip):" in app
    assert "def _close_info_target(target):" in app
    assert "def _hdr(content, tip):" in app
    assert "def _flag_row(widget, tip, info_target=None):" in app
    assert "'<summary aria-label=\"More information: %s\" title=\"%s\">&#9432;</summary>'" in app
    assert "'<details class=\"rxinfo-details\" data-revision=\"%d\">'" in app
    assert "'<div class=\"rxhelp-panel\" role=\"note\"><small>%s</small>" in app
    assert "_INFO_CONTROLS = weakref.WeakSet()" in app
    assert "description='Show information', icon='info-circle'" not in app
    assert "_rx_info_button" not in app
    assert "rxinfo-popover" not in app
    assert "📥 needs" not in app
    assert "run_log_fold = _collapsible('Run log', logbox)" in app
    assert "command_editor = W.VBox([" in app
    assert "W.HTML('<b>Command line</b>')" in app
    assert "command_editor = _collapsible('Command line'" not in app
    assert "viewer_toolbar = W.HBox(" in app
    assert "[view_input, pick_action, last_pick_info, view_controls, viewer_more]" in app
    assert "viewer_box = W.VBox([workflow_box, viewer_toolbar, view_input_note, selection_box])" in app
    assert "'<div class=\"rxmolstar-embed\">'" in app
    assert "workflow_contract_row = W.HBox([subreq, outputs_html]" in app
    assert "workflow_controls.add_class('rxworkflow-controls')" in app
    assert "grid-template-columns:minmax(280px,320px) minmax(150px,1fr) 240px" in app
    assert "height:clamp(700px,calc(100dvh - 24px),920px); overflow:hidden;" in app
    assert ".rxapp-main { flex:1 1 auto !important; min-height:0; overflow:hidden; }" in app
    assert ".rxpages { flex:1 1 auto !important; min-height:0; overflow:hidden; }" in app
    assert "overscroll-behavior:contain; scrollbar-gutter:stable;" in app
    assert "flex:0 1 clamp(420px,calc(133.333dvh - 600px),640px) !important;" in app
    assert "max-width:clamp(600px,calc(250dvh - 1320px),1000px);" in app
    assert ".rxpath-panel svg, .rxpath-panel img, .rxpath-panel canvas {" in app
    assert "traj_out = W.Output(layout={'width': '100%', 'min_width': '0'})" in app
    assert "plot_out = W.Output(layout={'width': '100%', 'min_width': '0'})" in app
    assert "'flex': '1 1 440px'" not in app
    assert ".rxcommand-dock {" in app
    assert "rootbox = W.VBox([header, app, cmdline_box])" in app
    assert "rootbox = W.VBox([header, app, W.HTML('<hr" not in app
    assert 'role="tooltip"' not in app
    assert "rxworkspace" in app
    assert "rxviewer" in app
    assert "rxinspector" in app
    assert ".rxapp .widget-button .fa, .rxapp .widget-upload .fa { display:none !important; }" in app
    assert "ngl_acc" not in app and "nglview" not in app
    assert "Appended to the command" not in app
    assert "Every remaining option" not in app
    assert "def _ingest_saved_files(" in app
    assert "input_file_rows" in app and "description='Remove file'" in app
    assert "upl = W.FileUpload(accept=_acc, multiple=True, description='Upload files'" in app
    assert "_drop = W.VBox([_drop_prompt, upl, _drop_epoch]," in app
    assert "align_items='center', justify_content='center'" in app
    assert "event.preventDefault();event.stopPropagation();" in app
    assert "box.addEventListener('drop',function(event)" in app
    assert "reader.readAsDataURL(file);" in app
    assert "callback,[out,item.id,item.generation],{})" in app
    assert "_DROP_STATE = {'generation': 0, 'seen': set()}" in app
    assert "def _claim_drop_batch(batch, generation):" in app
    assert "def _bump_drop_generation():" in app
    assert "queue.push({id:opaqueId(),generation:epoch(box),files:selected})" in app
    assert "if(busy||!queue.length)return;" in app
    assert "payload.batch===item.id" in app
    assert "if clear:\n        _bump_drop_generation()" in app
    assert app.count("mlmm_gui.on_drop") == 2
    assert "document.querySelectorAll('.rxapp .rxdrop')" in app
    assert "_dnd_out.register_callback('mlmm_gui.on_drop', _rxgui_drop)" in app
    assert "_accept_upload_pairs(pairs, 'drag & drop')" in app
    assert "def _delete_owned_uploads(paths):" in app
    assert "anywidget" not in app and "_HAS_DROP_WIDGET" not in app
    assert "description='Move earlier'" in app and "description='Move later'" in app
    assert "tooltip='Move earlier'" in app and "tooltip='Move later'" in app
    assert "command_footer = W.HBox(" in app
    assert "cmdline_box.add_class('rxcommand-dock')" in app
    assert "for _label, _page in _TAB_PAGES:" in app
    assert "_page.add_class('rxpage')" in app
    assert "_tab_body.add_class('rxpages')" in app
    assert "app.add_class('rxapp-main')" in app
    assert "def _advanced_coverage(" in app and "adv_extra" not in app
    assert "def _advanced_options(sub):" in app
    assert "adv_acc = _collapsible('Advanced flags', adv_box)" in app
    assert "every CLI option accounted for" not in app
    assert "radius_applies = (sub == 'extract' or" in app
    assert "adv_radius.disabled = not radius_applies" in app
    assert "_set_flag_visible(adv_radius, radius_applies)" in app
    assert "_set_flag_visible(adv_dftfb, _dftfb_applicable)" in app


def test_colab_viewer_persists_exact_atom_and_residue_context() -> None:
    app = _notebook()["cells"][2]["source"]

    for marker in (
        "'_last_pick': None", "'_pick_history': []", "def _remember_pick",
        "_VIEWER_GENERATION", "_MOLSTAR_VERSION = '5.6.1'",
        "layoutShowSequence:cfg.showSequence", "layoutShowControls:true",
        "collapseRightPanel:true", "layoutShowLog:false",
        "layoutShowLeftPanel:false", "layoutShowRemoteState:false",
        "await molstar.Viewer.create", "color:'element-symbol'",
        "alpha:0.55", "structure-component-static-water",
        "await configureWater()", "tryCreateComponentStatic(",
        "'mmcif' if fmt in ('cif', 'mmcif')",
        "ignoreStartupEmpty", "loci.kind==='empty-loci'",
        "clickQueue=clickQueue.then(()=>handleClick(event))",
        "await invoke('clear_highlights',[cfg.generation])",
        "cfg.generation,exact",
        "exact atom", "set current pick", "view_input", "_view_mapping_ok",
        "last_pick_info", "Generated file preview", "Download current run (.zip)",
        "results_box.add_class('rxresults')", "overflow-x:auto",
        "colab_run.log", "energy unavailable", "Command was cancelled",
        "Command failed", "_frame_link = W.jslink", "linked structure + energy",
        "ax.axvline(xs[i]", "artifact_fold._rx_set_open",
        "display(HTML(_molstar_iframe(fr[i], 'xyz', show_sequence=False)))",
        "display(HTML(_molstar_iframe(source, fmt, show_sequence=(fmt != 'xyz'))))",
    ):
        assert marker in app
    assert "py3Dmol" not in app
    assert "v.addSphere(" not in app
    assert "v.addBox(" not in app
    assert "_pick_box_edges" not in app
    assert "'greenCarbon'" not in app
    assert "viewportFocusBehavior" not in app
    assert "tryCreateComponentFromExpression" not in app
    assert "lociSelects.select" not in app
    assert "canvas.setProps" not in app
    assert "managers.interactivity.setProps" not in app
    assert "_viewer_seed_picks" not in app
    assert "skipInitialClick" not in app
    assert "style_pick" not in app
    assert "('Measure (dist/angle/dihedral)'," not in app
    assert "measure_panel" not in app
    assert "freeze_acc = _collapsible('Freezing', freeze_panel)" in app
    assert "_PICK_HINT.get(pick_action.value, 'Choose a click action.')" in app
    assert "_PICK_HINT[pick_action.value]" not in app
    assert "pick_action.value = 'scanB'" in app
    assert "pick_action.value = 'freezeB'" in app
    assert "target.invokeFunction(cfg.callback+'.'+suffix,args,{})" in app
    assert "String(item.sourceIndex)" in app
    assert "cfg.generation,exact" in app
    assert "if exact in (False, 0, 'false', 'False')" in app
    assert "Mol* focused this residue" in app
    click_completion = app[
        app.index("def on_click("):
        app.index("try:\n    from google.colab import output as _co")
    ]
    assert "if not live_marked:" in click_completion
    assert "selected in (False, 0, 'false', 'False')" not in click_completion
    assert "current_generation != _VIEWER_GENERATION['value']" in click_completion
    assert app.count("mlmm_gui.clear_highlights") == 1
    assert "def _resolve_click_meta(" in app
    assert "'viewer_index': viewer_index" in app
    assert "def _load_view_structure(path):" in app
    assert "def _invalidate_last_run(" in app
    assert "def _artifact_kind(" in app

    contract = _viewer_contract()
    metadata = [
        {"index": 0, "chain": "A", "resname": "LIG", "resseq": 10,
         "icode": "", "name": "N", "xyz": (0, 0, 0)},
        {"index": 1, "chain": "A", "resname": "LIG", "resseq": 10,
         "icode": "", "name": "C1", "xyz": (1, 0, 0)},
        {"index": 2, "chain": "B", "resname": "LIG", "resseq": 10,
         "icode": "", "name": "N", "xyz": (5, 0, 0)},
        {"index": 3, "chain": "A", "resname": "LIG", "resseq": 10,
         "icode": "A", "name": "N", "xyz": (8, 0, 0)},
    ]
    serial_metadata = [dict(row, serial=100 + i) for i, row in enumerate(metadata)]
    contract["S"].update(_atom_meta=serial_metadata, _view_format="pdb")
    # Stable PDB serial/signature wins over a drifted browser-side source index.
    assert contract["_resolve_click_meta"](
        0, "101", "LIG", "10", "A", "C1", "",
    )["index"] == 1
    contract["S"]["_atom_meta"] = metadata
    assert contract["_resolve_atom_query"]("A:LIG:10:C1") == 1
    assert contract["_resolve_atom_query"]("2") == 1
    contract["S"].update(
        _view_format="xyz",
        _atom_meta=[{"index": 0, "serial": 1, "resname": "MOL",
                     "resseq": 1, "name": "C1", "xyz": (0, 0, 0)}],
    )
    resolved_xyz = contract["_resolve_click_meta"](
        "0", "0", "", "undefined", "", "C", "",
    )
    assert resolved_xyz["name"] == "C1"
    contract["S"].update(center=["LIG"], center_ids=["B:LIG:10"],
                         _primary_atom_meta=metadata)
    assert contract["_center_cli_selectors"]() == [
        "A:LIG:10", "B:LIG:10", "A:LIG:10A",
    ]
    contract["S"].update(center=[], center_ids=["A:LIG:10"])
    with pytest.raises(ValueError, match="insertion-code sibling"):
        contract["_center_cli_selectors"]()
    contract["SPEC"].update({
        "bond-summary": {"n_in": (2, None)},
        "trj2fig": {"n_in": (1, 1)},
    })
    assert contract["_input_count_error"]("bond-summary", 1)
    assert not contract["_input_count_error"]("bond-summary", 2)
    assert contract["_input_count_error"]("trj2fig", 2)
    assert contract["_artifact_kind"]("profile.html") == "interactive HTML"
    assert contract["_artifact_kind"]("plot.jpeg") == "image"
    assert contract["_artifact_kind"]("final.pdb") == "structure"
    assert contract["_artifact_kind"]("final.cif") == "structure"
    assert contract["_artifact_kind"]("final.xyz") == "structure"
    assert contract["_artifact_kind"]("optimization_trj.xyz") is None
    opt = contract["_trajectory_semantics"]("opt", "optimization_trj.xyz")
    path_semantics = contract["_trajectory_semantics"](
        "path-opt", "mep_trj.xyz",
    )
    vibration = contract["_trajectory_semantics"](
        "tsopt", "vib/imag_120i_trj.xyz",
    )
    assert vibration["title"] == "Vibrational-mode animation"
    assert vibration["x"] == "phase frame" and not vibration["extrema"]
    assert contract["_stationary"]([0.0, 2.0, 0.0], opt) == [
        (0, "initial"), (2, "optimized"),
    ]
    assert (1, "peak candidate") in contract["_stationary"](
        [0.0, 2.0, 0.0], path_semantics,
    )


def test_colab_app_executes_atomic_view_and_result_transitions(
    tmp_path: Path, monkeypatch,
) -> None:
    app, calls = _execute_app(monkeypatch, tmp_path)
    drop_generation = app["_DROP_STATE"]["generation"]
    assert app["_claim_drop_batch"]("batch-a", drop_generation) == (True, "")
    assert app["_claim_drop_batch"]("batch-a", drop_generation) == (
        False, "duplicate batch",
    )
    app["b_clear_inputs"].click()
    assert app["_DROP_STATE"]["generation"] == drop_generation + 1
    assert app["_claim_drop_batch"]("late-batch", drop_generation) == (
        False, "stale generation",
    )
    assert app["workspace"].layout.display == "none"
    assert app["selection_help"].layout.display == "none"
    assert app["selection_route"].layout.display == ""
    primary = tmp_path / "primary.pdb"
    secondary = tmp_path / "secondary.pdb"
    primary_text = (
        "HETATM    1  C1  LIG A  10       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  LIG A  10       1.200   0.000   0.000  1.00  0.00           O\n"
        "HETATM    3  C1  COF B  20       4.000   0.000   0.000  1.00  0.00           C\nEND\n"
    )
    secondary_text = (
        "HETATM    1  C1  ALT A  10       0.100   0.000   0.000  1.00  0.00           C\nEND\n"
    )
    primary.write_text(primary_text, encoding="utf-8")
    secondary.write_text(secondary_text, encoding="utf-8")
    calc_file = tmp_path / "calculator.py"
    read_hess = tmp_path / "initial.hess.npy"
    calc_file.write_text("calculator = 1\n", encoding="utf-8")
    read_hess.write_bytes(b"hessian")
    validation_argv = [
        "mlmm", "irc", "-i", str(primary),
        "--calc-file", str(calc_file), "--read-hess", str(read_hess),
    ]
    app["cmd_box"].value = shlex.join(validation_argv)
    assert set(app["_command_input_files"](validation_argv)) == {
        str(primary), str(calc_file), str(read_hess),
    }
    first_fingerprint = app["_validation_fingerprint"](validation_argv)
    read_hess.write_bytes(b"updated hessian")
    assert app["_validation_fingerprint"](validation_argv) != first_fingerprint
    extra_structure = tmp_path / "extra.pdb"
    extra_structure.write_text(secondary_text, encoding="utf-8")
    positional_argv = ["mlmm", "bond-summary", str(primary), str(extra_structure)]
    assert set(app["_command_input_files"](positional_argv)) == {
        str(primary), str(extra_structure),
    }
    scan_spec = tmp_path / "scan.yaml"
    scan_spec.write_text("scan: first\n", encoding="utf-8")
    scan_argv = ["mlmm", "scan", "-i", str(primary), "-s", str(scan_spec)]
    app["cmd_box"].value = shlex.join(scan_argv)
    assert set(app["_command_input_files"](scan_argv)) == {str(primary), str(scan_spec)}
    scan_fingerprint = app["_validation_fingerprint"](scan_argv)
    scan_spec.write_text("scan: second\n", encoding="utf-8")
    assert app["_validation_fingerprint"](scan_argv) != scan_fingerprint
    metadata = {
        str(primary): [
            {"chain": "A", "resname": "LIG", "resseq": 10, "icode": "", "name": "C1"},
            {"chain": "A", "resname": "LIG", "resseq": 10, "icode": "", "name": "O1"},
            {"chain": "B", "resname": "COF", "resseq": 20, "icode": "", "name": "C1"},
        ],
        str(secondary): [
            {"chain": "A", "resname": "ALT", "resseq": 10, "icode": "", "name": "C1"},
        ],
    }

    def load_view(path):
        path = str(path)
        return Path(path).read_text(encoding="utf-8"), [dict(row) for row in metadata[path]], path

    app["_load_view_structure"] = load_view
    app["load_pdb"]([str(primary), str(secondary)], center=["COF"], lcharge={"COF": 1})
    primary_widget = app["center_widget"]
    app["pick_action"].value = "center"
    app["on_click"]("0", "LIG", "10", "A", "C1")
    assert app["_center_cli_selectors"]() == ["A:LIG:10", "B:COF:20"]

    from Bio.PDB import PDBParser
    from mlmm.workflows.extract import resolve_substrate_residues
    structure = PDBParser(QUIET=True).get_structure("primary", primary)
    resolved = resolve_substrate_residues(structure, "B:COF:20,A:LIG:10")
    assert {(residue.get_parent().id, residue.id[1], residue.resname) for residue in resolved} == {
        ("A", 10, "LIG"), ("B", 20, "COF"),
    }

    saved = (list(app["S"]["center"]), list(app["S"]["center_ids"]), dict(app["S"]["lcharge"]))
    app["S"]["scan_atoms"] = [
        {"index": 1, "chain": "A", "resn": "LIG", "resi": "10", "atom": "O1",
         "xyz": (1.2, 0.0, 0.0)},
        None,
    ]
    app["S"]["freeze_atoms"] = [2]
    app["S"]["measure_atoms"] = [
        {"index": 0, "chain": "A", "resn": "LIG", "resi": "10", "atom": "C1",
         "xyz": (0.0, 0.0, 0.0)},
    ]
    calls.clear()
    app["view_input"].value = 1
    assert app["S"]["_view_mapping_ok"] is False
    assert app["center_widget"] is primary_widget
    assert primary_widget.disabled is True
    assert (app["S"]["center"], app["S"]["center_ids"], app["S"]["lcharge"]) == saved
    assert any(
        isinstance(value, str) and 'class="rxmolstar-frame"' in value
        for value in calls
    )
    app["dd_subcmd"].value = "scan"
    def _descendants(widget):
        yield widget
        for child in getattr(widget, "children", ()):
            yield from _descendants(child)
    clear_pair = next(
        widget for widget in _descendants(app["scan_panel"])
        if getattr(widget, "description", "") == "Clear pair"
    )
    assert clear_pair.disabled
    before_scan = json.dumps(app["S"]["scan_atoms"], sort_keys=True)
    clear_pair.click()
    assert json.dumps(app["S"]["scan_atoms"], sort_keys=True) == before_scan
    app["view_input"].value = 0
    assert app["center_widget"].disabled is False

    old_text = app["S"]["_pdb_text"]
    app["_load_view_structure"] = lambda path: (_ for _ in ()).throw(ValueError("view failed")) \
        if str(path) == str(secondary) else load_view(path)
    app["view_input"].value = 1
    assert app["S"]["_view_input_index"] == 0
    assert app["view_input"].value == 0
    assert app["S"]["_pdb_text"] == old_text

    app["S"].update(_last_out_dir="old", _last_subcmd="opt", _last_argv=["old"],
                    _last_files=[str(primary)], _last_manifest={"status": "success"},
                    _last_log="old log")
    app["artifact_choice"].options = [("old", str(primary))]
    app["dl_btn"].disabled = False
    validate_out = tmp_path / "validate-only"
    app["cmd_box"].value = shlex.join([
        "mlmm", "irc", "-i", str(primary), "--calc-file", str(calc_file),
        "--read-hess", str(read_hess), "-o", str(validate_out),
    ])
    assert app["S"]["_last_manifest"] == {} and app["S"]["_last_files"] == []
    assert app["artifact_choice"].disabled and app["dl_btn"].disabled
    app["S"]["_last_log"] = "old log"
    app["_stream"] = lambda argv: (0, "validation transcript")
    app["_do_validate"](None)
    assert app["S"]["_last_log"] == "old log"
    assert app["_RUN_STATE"]["validation_log"] == "validation transcript"
    app["_stream"] = lambda argv: (2, "invalid options")
    app["_do_validate"](None)
    assert app["run_log_fold"].children[1].layout.display == ""

    # One path-position change drives both the Mol* structure and the energy
    # cursor/status from the same slider value.
    trajectory = tmp_path / "linked_path_trj.xyz"
    trajectory.write_text(
        "2\n-1.000000\n"
        "H 0.000 0.000 0.000\nH 0.700 0.000 0.000\n"
        "2\n-0.990000\n"
        "H 1.500 0.000 0.000\nH 2.200 0.000 0.000\n",
        encoding="utf-8",
    )
    calls.clear()
    app["S"]["_last_subcmd"] = "path-opt"
    app["_load_trajectory"](str(trajectory), str(tmp_path))
    assert app["frame_slider"].max == 1 and not app["frame_slider"].disabled
    assert app["trajectory_box"].layout.display == ""
    assert "Frame 1 of 2" in app["frame_state"].value
    assert "ΔE = 0.0 kcal/mol" in app["frame_state"].value
    first_frame = [
        value for value in calls
        if isinstance(value, str) and 'class="rxmolstar-frame"' in value
    ]
    assert first_frame and "H 0.000 0.000 0.000" in first_frame[-1]

    calls.clear()
    app["frame_slider"].value = 1
    assert "Frame 2 of 2" in app["frame_state"].value
    assert "ΔE = 6.3 kcal/mol" in app["frame_state"].value
    second_frame = [
        value for value in calls
        if isinstance(value, str) and 'class="rxmolstar-frame"' in value
    ]
    assert second_frame and "H 1.500 0.000 0.000" in second_frame[-1]
    assert not app["frame_prev"].disabled and app["frame_next"].disabled

    app["_invalidate_last_run"]("Input identity changed; run again.")
    assert app["S"]["_last_manifest"] == {} and app["S"]["_last_files"] == []
    assert app["artifact_choice"].disabled and app["dl_btn"].disabled
    assert "Input identity changed" in app["results_empty"].value


def test_colab_compact_selection_upload_viewer_and_advanced_contracts(
    tmp_path: Path, monkeypatch,
) -> None:
    app, calls = _execute_app(monkeypatch, tmp_path)
    primary = tmp_path / "primary.pdb"
    secondary = tmp_path / "secondary.pdb"
    topology = tmp_path / "system.parm7"
    primary.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  H   ALA A   1       0.000   1.000   0.000  1.00  0.00           H\n"
        "HETATM    3 MG    MG A   2       2.000   0.000   0.000  1.00  0.00          MG\n"
        "HETATM    4  O   WAT A  10       3.000   0.000   0.000  1.00  0.00           O\n"
        "HETATM    5  C1  LIG A   3       4.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    6  O1  LIG A   3       5.200   0.000   0.000  1.00  0.00           O\n"
        "END\n",
        encoding="utf-8",
    )
    secondary.write_text(
        "HETATM    1  C1  ALT B   4       8.000   0.000   0.000  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    topology.write_text("parm", encoding="utf-8")
    metadata = {
        str(primary): [
            {"serial": 1, "chain": "A", "resname": "ALA", "resseq": 1, "icode": "", "name": "N"},
            {"serial": 2, "chain": "A", "resname": "ALA", "resseq": 1, "icode": "", "name": "H"},
            {"serial": 3, "chain": "A", "resname": "MG", "resseq": 2, "icode": "", "name": "MG"},
            {"serial": 4, "chain": "A", "resname": "WAT", "resseq": 10, "icode": "", "name": "O"},
            {"serial": 5, "chain": "A", "resname": "LIG", "resseq": 3, "icode": "", "name": "C1"},
            {"serial": 6, "chain": "A", "resname": "LIG", "resseq": 3, "icode": "", "name": "O1"},
        ],
        str(secondary): [
            {"serial": 1, "chain": "B", "resname": "ALT", "resseq": 4, "icode": "", "name": "C1"},
        ],
    }

    def load_view(path):
        path = str(path)
        return Path(path).read_text(encoding="utf-8"), [dict(row) for row in metadata[path]], path

    app["_load_view_structure"] = load_view
    assert app["_ingest_saved_files"]([str(primary), str(topology)], "test drop")
    assert app["all_mode"].value == "scan"
    assert app["pick_action"].value == "scanA"
    assert app["workspace"].layout.display == ""
    assert app["selection_help"].layout.display == "none"
    assert app["selection_route"].layout.display == "none"
    assert app["_ingest_saved_files"]([str(secondary)], "test drop")
    assert app["all_mode"].value == "mep"
    assert app["S"]["inputs"] == [str(primary), str(secondary)]
    assert app["S"]["parm"] == str(topology)
    assert len(app["input_file_rows"].children) == 3
    assert app["prep_radius"].value == pytest.approx(2.6)
    assert app["adv_radius"].value == pytest.approx(2.6)
    assert app["dd_col"].value == "element"
    assert app["dd_rep"].value == "cartoon"
    assert app["dd_size"].value == 720
    assert app["S"]["viewer_width"] == 720
    assert app["S"]["viewer_height"] == 540

    info = app["last_pick_info"]
    assert '<details class="rxinfo-details"' in info.value
    assert "<summary" in info.value and "&#9432;" in info.value
    assert 'role="note"' in info.value
    assert " open" not in info.value
    second_info = app["_info_control"]("second")
    assert "second" in second_info.value
    assert 'aria-label="More information: second"' in second_info.value
    assert 'title="second"' in second_info.value
    app["_set_info_text"](second_info, "updated")
    assert "updated" in second_info.value and "second" not in second_info.value
    before_close = second_info.value
    app["_close_info"](second_info)
    assert second_info.value != before_close
    assert "<details" in second_info.value and " open" not in second_info.value
    assert app["key_opts_box"].children[1].layout.display == "none"
    assert app["command_editor"].children[1] is app["cmd_box"]
    assert app["command_editor"].layout.display != "none"
    assert app["run_log_fold"].children[1].layout.display == "none"
    assert set(app["_center_values"]()) == {"LIG", "MG"}
    mg_row = app["charge_rows"]["MG"]
    assert mg_row["auto"] and mg_row["use"].disabled and mg_row["val"].disabled
    assert mg_row["val"].value == 2
    assert "MG" not in app["S"]["lcharge"]

    calls.clear()
    app["cb_water"].value = True
    assert app["S"]["show_water"] is True
    frames = [
        value for value in calls
        if isinstance(value, str) and 'class="rxmolstar-frame"' in value
    ]
    assert frames
    document = app["_molstar_document"](
        app["S"]["_pdb_text"], "pdb", show_water=True, interactive=True,
    )
    assert '"showWater":true' in document
    assert "layoutShowControls:true" in document
    assert "layoutShowSequence:cfg.showSequence" in document
    assert "structure-component-static-water" in document
    cif_document = app["_molstar_document"](
        "data_demo\n_entry.id demo\n", "cif", show_sequence=True,
    )
    assert '"format":"mmcif"' in cif_document
    assert '"showSequence":true' in cif_document

    calls.clear()
    app["pick_action"].value = "center"
    app["on_click"]("1", "LIG", "3", "A", "C1", "5", "")
    assert app["S"]["_last_pick"]["index"] == 4
    assert app["S"]["_last_pick"]["viewer_index"] == 1
    assert any(
        isinstance(value, str) and 'class="rxmolstar-frame"' in value
        for value in calls
    )
    assert [pick["index"] for pick in app["S"]["_pick_history"]] == [4]

    calls.clear()
    app["on_click"]("2", "MG", "2", "A", "MG", "3", "", live_marked=True)
    assert [pick["index"] for pick in app["S"]["_pick_history"]] == [4, 2]
    assert not any(
        isinstance(value, str) and 'class="rxmolstar-frame"' in value
        for value in calls
    )
    committed_centers = list(app["S"]["center_ids"])
    app["_clear_highlights_from_browser"](app["_VIEWER_GENERATION"]["value"])
    assert app["S"]["_pick_history"] == []
    assert app["S"]["center_ids"] == committed_centers

    app["on_click"]("2", "MG", "2", "A", "MG", "3", "", live_marked=True)
    assert [pick["index"] for pick in app["S"]["_pick_history"]] == [2]

    close_secondary = app["input_file_rows"].children[1].children[1].children[-1]
    close_secondary.click()
    assert app["S"]["inputs"] == [str(primary)]
    assert app["S"]["center_ids"] == ["A:LIG:3", "A:MG:2"]
    assert app["S"]["parm"] == str(topology)

    all_statuses = {}
    root = app["PRODUCT_CLI"]
    for subcommand in root.list_commands(app["click"].Context(root)):
        options = app["_advanced_options"](subcommand)
        coverage = app["_advanced_coverage"](subcommand)
        assert set(coverage) == {param.name for param in options}
        assert set(coverage.values()) <= {"owned", "generated", "blocked", "rendered"}
        all_statuses.update(coverage)
    assert all_statuses["tr_projection"] == "blocked"
    assert all_statuses["embedcharge"] == "blocked"
    assert all_statuses["embedcharge_cutoff"] == "blocked"
    assert app["_advanced_coverage"]("dft")
    assert app["_advanced_coverage"]("sp")["hessian_calc_mode"] == "rendered"

    app["cb_advsub"].value = True
    app["dd_subcmd"].value = "add-elem-info"
    app["S"]["advanced_overrides"]["add-elem-info"] = {}
    safe_add_elem = app["build_cmd"]()
    assert "-o" in safe_add_elem and "--inplace" not in safe_add_elem
    app["S"]["advanced_overrides"]["add-elem-info"] = {"inplace": True}
    inplace_add_elem = app["build_cmd"]()
    assert "-o" not in inplace_add_elem and "--inplace" in inplace_add_elem
    # Field re-inference remains independent and keeps the safe output.
    app["S"]["advanced_overrides"]["add-elem-info"] = {"overwrite": True}
    field_overwrite = app["build_cmd"]()
    assert "-o" in field_overwrite and "--overwrite" in field_overwrite

    app["dd_subcmd"].value = "all"
    app["all_mode"].value = "mep"
    app["_render_advanced_rows"]()
    editable = [
        param for param in app["_advanced_options"]("all")
        if app["_advanced_status"]("all", param) == "rendered"
        and app["_advanced_semantic_applicable"]("all", param.name)
    ]
    assert len(app["adv_rows_box"].children) == len(editable)

    click = app["click"]
    text_param = next(
        param for param in editable
        if not param.multiple
        and not param.is_bool_flag
        and not isinstance(param.type, (click.Choice, click.types.BoolParamType))
    )
    bool_param = next(
        param for param in editable
        if param.is_bool_flag or isinstance(param.type, click.types.BoolParamType)
    )
    app["S"]["advanced_overrides"]["all"] = {
        text_param.name: "7", bool_param.name: False,
    }
    advanced_argv = app["_advanced_argv"]("all")
    text_flag = app["_advanced_flag"](text_param)
    assert advanced_argv[advanced_argv.index(text_flag) + 1] == "7"
    bool_flags = set(bool_param.opts + bool_param.secondary_opts)
    assert bool_flags.intersection(advanced_argv)

    sp_hessian = next(param for param in app["_advanced_options"]("sp")
                      if param.name == "hessian_calc_mode")
    hessian_value = (list(sp_hessian.type.choices)[0]
                     if isinstance(sp_hessian.type, click.Choice) else "analytical")
    app["S"]["advanced_overrides"]["sp"] = {"hessian_calc_mode": hessian_value}
    assert "--hessian-calc-mode" in app["_advanced_argv"]("sp")

    for utility, flag in (("trj2fig", "--out-json"),
                          ("energy-diagram", "--out-json"),
                          ("bond-summary", "--json")):
        out_json = next(param for param in app["_advanced_options"](utility)
                        if param.name == "out_json")
        assert app["_advanced_status"](utility, out_json) == "rendered"
        app["S"]["advanced_overrides"][utility] = {"out_json": True}
        assert flag in app["_advanced_argv"](utility)

    app["dd_subcmd"].value = "sp"
    assert app["key_opts_box"].layout.display == ""
    assert [value for _label, value in app["pick_action"].options] == [
        "center", "ligand", "freezeA", "freezeB", "freezeatom",
    ]
    assert app["center_panel"].layout.display == ""
    assert app["charge_panel"].layout.display == ""
    assert app["extract_panel"].layout.display == ""
    assert app["adv_radius"]._rx_flag_row.layout.display == "none"
    assert app["adv_dftfb"]._rx_flag_row.layout.display == "none"
    app["dd_subcmd"].value = "freq"
    assert app["key_opts_box"].layout.display == ""
    app["cb_advsub"].value = True
    app["dd_subcmd"].value = "extract"
    assert app["key_opts_box"].layout.display == ""
    assert app["adv_radius"]._rx_flag_row.layout.display == ""
    app["center_widget"].value = ()
    app["S"]["center_ids"] = []
    with pytest.raises(ValueError, match="needs a center residue"):
        app["build_cmd"]()
    app["center_widget"].value = ("LIG",)
    app["dd_subcmd"].value = "sp"
    app["prep_radius"].value = 4.2
    assert app["adv_radius"].value == 4.2
    extract_commands = []

    def fake_extract(command, **_kwargs):
        extract_commands.append(list(command))
        output = Path(command[command.index("-o") + 1])
        output.write_text(primary.read_text(encoding="utf-8"), encoding="utf-8")
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(app["subprocess"], "run", fake_extract)
    app["b_extract"].click()
    assert extract_commands and extract_commands[-1][-2:] == ["-r", "4.2"]
    assert app["dd_subcmd"].value == "sp" and app["S"]["subcmd"] == "sp"
    assert app["S"]["model_pdb"]
    assert app["b_revert"].layout.display == ""
    app["b_revert"].click()
    assert app["S"]["model_pdb"] is None
    assert app["dd_subcmd"].value == "sp" and app["S"]["subcmd"] == "sp"
    app["dd_subcmd"].value = "add-elem-info"
    assert app["center_panel"].layout.display == "none"
    assert app["charge_panel"].layout.display == "none"
    assert app["extract_panel"].layout.display == "none"
    app["dd_subcmd"].value = "all"
    app["all_mode"].value = "scan"
    assert {"scanA", "scanB"} <= {value for _label, value in app["pick_action"].options}
    app["all_mode"].value = "mep"
    assert not {"scanA", "scanB"} & {value for _label, value in app["pick_action"].options}

    app["S"]["advanced_overrides"]["all"] = {}
    app["adv_refine"].value = True
    app["adv_mep"].value = "gsm"
    app["adv_thresh"].value = "gau"
    app["adv_maxcyc"].value = 9
    app["S"]["advanced_overrides"]["all"] = {
        "reject_uphill": False,
        "opt_mode": "grad",
        "pre_opt": False,
        "irc_step_size": "0.05",
        "opt_mode_post": "grad",
        "thresh_post": "baker",
        "hessian_calc_mode": "FiniteDifference",
        "skip_final_freq": True,
    }
    app["w_ts"].value = False
    app["w_th"].value = False
    app["adv_dft"].value = False
    off_argv = app["_advanced_argv"]("all")
    for flag in (
        "--irc-step-size", "--opt-mode-post", "--thresh-post",
        "--hessian-calc-mode", "--skip-final-freq", "--no-reject-uphill",
    ):
        assert flag not in off_argv
    app["all_mode"].value = "tsonly"
    assert app["adv_refine"].layout.display == "none"
    assert app["w_ts"].value and app["w_ts"].disabled
    command = app["build_cmd"]()
    for flag in ("--refine-path", "--mep-mode", "--thresh", "--max-cycles"):
        assert flag not in command
    assert "--no-reject-uphill" in command
    assert command[command.index("--opt-mode") + 1] == "grad"
    assert command[command.index("--irc-step-size") + 1] == "0.05"
    assert command[command.index("--opt-mode-post") + 1] == "grad"
    assert command[command.index("--thresh-post") + 1] == "baker"
    assert command[command.index("--hessian-calc-mode") + 1] == "FiniteDifference"
    assert "--skip-final-freq" in command
    assert "--preopt" not in command and "--no-preopt" not in command
    assert "--tsopt" in command
    app["all_mode"].value = "mep"
    assert not app["w_ts"].value and not app["w_ts"].disabled
    app["w_th"].value = True
    thermo_argv = app["_advanced_argv"]("all")
    for flag in (
        "--irc-step-size", "--opt-mode-post", "--thresh-post",
        "--hessian-calc-mode",
    ):
        assert flag in thermo_argv
    assert "--no-reject-uphill" in thermo_argv
    assert "--skip-final-freq" not in thermo_argv
    app["w_th"].value = False
    app["S"]["inputs"] = [str(primary), str(secondary)]
    assert "--tsopt" not in app["build_cmd"]()

    result_json = tmp_path / "result.json"
    result_json.write_text('{"energy": -1.25}', encoding="utf-8")
    assert app["_artifact_kind"](str(result_json)) == "JSON"
    preview = app["_text_preview_html"](str(result_json), "JSON")
    assert "&quot;energy&quot;" in preview and "-1.25" in preview
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps({
        "status": "success", "scientific_status": "partial",
        "scientific_status_reasons": ["IRC endpoint mismatch"],
        "segments": [{"index": 1, "barrier_kcal": 8.0, "delta_kcal": -1.0}],
        "post_segments": [{"index": 1, "mlip": {"barrier_kcal": 7.5, "delta_kcal": -1.2}}],
        "rate_limiting_step": {"barrier_kcal": 7.5, "method": "mlip"},
    }), encoding="utf-8")
    summary_html = app["_summary_html"](str(summary))
    assert "Provisional barrier" in summary_html
    assert "IRC endpoint mismatch" in summary_html
    assert "refined MLIP" in summary_html and "⚡" not in summary_html
    calls.clear()
    app["_structure_preview"](str(primary))
    assert any(
        isinstance(value, str) and 'class="rxmolstar-frame"' in value
        and "Mol*" in value
        for value in calls
    )

    # Charge semantics and controls follow the selected ML/MM workflow.
    app["dd_subcmd"].value = "all"
    assert app["w_q"].description == "system charge (-q)"
    assert "full-system charge override" in app["charge_info"].value
    app["dd_subcmd"].value = "sp"
    assert app["w_q"].description == "ML charge (-q)"
    assert "ML region" in app["charge_info"].value
    app["cb_advsub"].value = True
    app["dd_subcmd"].value = "oniom-export"
    assert app["w_q"].description == "QM charge (-q)"
    assert "QM region" in app["charge_info"].value

    # mm-parm exposes its owned ligand-charge editor and gates AmberTools.
    app["S"]["inputs"] = [str(primary)]
    app["dd_subcmd"].value = "mm-parm"
    assert app["charge_panel"].layout.display == ""
    assert app["center_panel"].layout.display == "none"
    app["charge_rows"]["LIG"]["use"].value = True
    app["charge_rows"]["LIG"]["val"].value = 1
    app["_ambertools_available"] = lambda: False
    with pytest.raises(ValueError, match="needs tleap"):
        app["build_cmd"]()
    app["_ambertools_available"] = lambda: True
    mm_parm_command = app["build_cmd"]()
    assert mm_parm_command[1] == "mm-parm" and mm_parm_command[-2:] == ["-l", "LIG:1"]

    # Utility files route directly to their workflow; scan3d CSV needs no parm7.
    oniom = tmp_path / "job.gjf"
    oniom.write_text("# oniom\n", encoding="utf-8")
    assert app["_ingest_saved_files"]([str(oniom)], "test utility")
    oniom_command = app["build_cmd"]()
    assert oniom_command[1] == "oniom-import" and oniom_command[oniom_command.index("-i") + 1] == str(oniom)
    surface = tmp_path / "surface.csv"
    surface.write_text("d1,d2,d3,energy\n1,1,1,0\n", encoding="utf-8")
    assert app["_ingest_saved_files"]([str(surface)], "test csv")
    csv_command = app["build_cmd"]()
    assert csv_command[:2] == ["mlmm", "scan3d"]
    assert "--csv" in csv_command and "--parm" not in csv_command and "-i" not in csv_command

    app["S"]["_pre_extract"] = {"model_pdb": None}
    app["b_revert"].layout.display = ""
    app["_clear_structure_bound_state"]()
    assert app["S"]["_pre_extract"] is None
    assert app["b_revert"].layout.display == "none"


def test_colab_prepared_model_upload_keeps_full_system_inputs(
    tmp_path: Path, monkeypatch,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)
    full = tmp_path / "full.pdb"
    parm = tmp_path / "full.parm7"
    full.write_text("END\n", encoding="utf-8")
    parm.write_text("parm", encoding="utf-8")
    app["S"].update(inputs=[str(full)], parm=str(parm), mode="pdb")
    payload = b"HETATM    1  C1  LIG A  10       0.000   0.000   0.000  1.00  0.00           C\nEND\n"
    app["model_upl"].value = ({"name": "model.pdb", "type": "chemical/x-pdb",
                                "size": len(payload), "content": memoryview(payload),
                                "last_modified": datetime.datetime.now(datetime.timezone.utc)},)
    assert app["S"]["inputs"] == [str(full)] and app["S"]["parm"] == str(parm)
    assert app["S"]["model_pdb"] and Path(app["S"]["model_pdb"]).name == "model.pdb"
    assert app["model_upl"].value == ()
    first_model = Path(app["S"]["model_pdb"])
    assert first_model.exists()
    app["model_clear"].click()
    assert app["S"]["model_pdb"] is None
    assert not first_model.exists()
    app["model_upl"].value = ({"name": "model.pdb", "type": "chemical/x-pdb",
                                "size": len(payload), "content": memoryview(payload),
                                "last_modified": datetime.datetime.now(datetime.timezone.utc)},)
    assert app["S"]["model_pdb"] and Path(app["S"]["model_pdb"]).name.startswith("model")
    assert app["model_upl"].value == ()


def test_colab_adversarial_session_upload_and_view_state(
    tmp_path: Path, monkeypatch,
) -> None:
    app, calls = _execute_app(monkeypatch, tmp_path)
    primary = tmp_path / "primary.pdb"
    product = tmp_path / "product.pdb"
    incompatible = tmp_path / "incompatible.pdb"
    topology = tmp_path / "system.parm7"
    primary.write_text(
        "HETATM    1  C1  LIG A  10       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  LIG A  10       1.200   0.000   0.000  1.00  0.00           O\nEND\n",
        encoding="utf-8",
    )
    product.write_text(
        "HETATM    1  C1  LIG A  10      10.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  LIG A  10      12.400   0.000   0.000  1.00  0.00           O\nEND\n",
        encoding="utf-8",
    )
    incompatible.write_text(
        "HETATM    1  C1  ALT B  20       4.000   0.000   0.000  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    topology.write_text("%VERSION VERSION_STAMP = V0001.000\n", encoding="utf-8")
    metadata = {
        str(primary): [
            {"serial": 1, "chain": "A", "resname": "LIG", "resseq": 10, "icode": "", "name": "C1"},
            {"serial": 2, "chain": "A", "resname": "LIG", "resseq": 10, "icode": "", "name": "O1"},
        ],
        str(product): [
            {"serial": 1, "chain": "A", "resname": "LIG", "resseq": 10, "icode": "", "name": "C1"},
            {"serial": 2, "chain": "A", "resname": "LIG", "resseq": 10, "icode": "", "name": "O1"},
        ],
        str(incompatible): [
            {"serial": 1, "chain": "B", "resname": "ALT", "resseq": 20, "icode": "", "name": "C1"},
        ],
    }
    real_loader = app["_load_view_structure"]

    def load_view(path):
        path = str(path)
        return Path(path).read_text(encoding="utf-8"), [dict(row) for row in metadata[path]], path

    app["_load_view_structure"] = load_view
    assert app["b_run"].disabled and app["b_validate"].disabled
    app["cmd_box"].value = "mlmm --version"
    assert not app["b_run"].disabled and app["b_validate"].disabled
    app["cmd_box"].value = "# incomplete"
    assert app["b_run"].disabled

    advanced_row = next(row for row in app["adv_rows_box"].children if hasattr(row, "_rx_search"))
    advanced_info = advanced_row.children[1]
    assert '<details class="rxinfo-details"' in advanced_info.value
    assert 'role="note"' in advanced_info.value
    app["_render_advanced_rows"]()
    assert " open" not in advanced_info.value

    # Upload order is arbitrary: a queued topology survives the first PDB.
    assert app["_ingest_saved_files"]([str(topology)], "topology first")
    assert app["_ingest_saved_files"]([str(primary)], "structure second")
    assert app["S"]["inputs"] == [str(primary)] and app["S"]["parm"] == str(topology)
    assert "single-structure input" in app["input_order_note"].value

    app["load_pdb"]([str(primary), str(product)], parm=str(topology), center=["LIG"])
    assert "reaction order shown above" in app["input_order_note"].value
    app["S"]["scan_atoms"] = [
        {"index": 0, "chain": "A", "resn": "LIG", "resi": "10", "atom": "C1", "xyz": (0.0, 0.0, 0.0)},
        {"index": 1, "chain": "A", "resn": "LIG", "resi": "10", "atom": "O1", "xyz": (1.2, 0.0, 0.0)},
    ]
    app["view_input"].value = 1
    assert app["S"]["scan_atoms"][0]["xyz"] == pytest.approx((10.0, 0.0, 0.0))
    assert app["S"]["scan_atoms"][1]["xyz"] == pytest.approx((12.4, 0.0, 0.0))
    calls.clear()
    app["render_viewer"]()
    assert any(
        isinstance(value, str) and 'class="rxmolstar-frame"' in value
        for value in calls
    )

    app["load_pdb"](
        [str(primary), str(incompatible)], parm=str(topology), center=["LIG"], lcharge={"LIG": 1},
    )
    app["S"]["scan_atoms"] = [
        {"index": 0, "chain": "A", "resn": "LIG", "resi": "10", "atom": "C1", "xyz": (0.0, 0.0, 0.0)},
        {"index": 1, "chain": "A", "resn": "LIG", "resi": "10", "atom": "O1", "xyz": (1.2, 0.0, 0.0)},
    ]
    primary_owned = json.dumps(app["S"]["scan_atoms"], sort_keys=True)
    app["view_input"].value = 1
    assert json.dumps(app["S"]["scan_atoms"], sort_keys=True) == primary_owned
    before_selection = (list(app["S"]["center"]), dict(app["S"]["lcharge"]))
    assert app["b_clear"].disabled and app["center_widget"].disabled
    app["_clear_sel"](None)
    assert (app["S"]["center"], app["S"]["lcharge"]) == before_selection
    app["view_input"].value = 0
    app["_clear_sel"](None)
    assert app["S"]["center"] == [] and app["S"]["lcharge"] == {}

    app["S"]["center_ids"] = ["A:LIG:10"]
    app["S"]["_last_manifest"] = {"status": "success"}
    before = json.dumps(app["_session_dict"](), sort_keys=True)
    with pytest.raises(ValueError, match="one object"):
        app["_apply_session"]([])
    assert json.dumps(app["_session_dict"](), sort_keys=True) == before
    assert app["S"]["_last_manifest"] == {"status": "success"}
    wrong_tool = app["_session_dict"]()
    wrong_tool["tool"] = "pdb2reaction"
    with pytest.raises(ValueError, match="belong"):
        app["_apply_session"](wrong_tool)

    saved = app["_session_dict"]()
    saved.update(
        inputs=[str(primary), str(product)], parm=str(topology), mode="pdb",
        subcmd="all", all_mode="mep", tsopt=False, backend="uma", model="uma-m-1p1",
        rep="sticks", color="spectrum",
    )
    app["all_mode"].value = "tsonly"
    app["_ALL_MODE_STATE"]["tsopt_before_tsonly"] = True
    stale_command = "mlmm sp -i stale.pdb --parm stale.parm7 -q 99"
    app["cmd_box"].value = stale_command
    with pytest.raises(ValueError, match="not installed"):
        app["_apply_session"](saved)
    saved["backend"] = app["BACKEND"]
    saved["model"] = app["DEFAULT_MODEL"][app["BACKEND"]]
    assert app["_apply_session"](saved) == []
    assert app["all_mode"].value == "mep" and app["w_ts"].value is False
    assert app["dd_backend"].value == app["BACKEND"]
    assert app["dd_model"].value == app["DEFAULT_MODEL"][app["BACKEND"]]
    assert app["dd_rep"].value == "sticks" and app["dd_col"].value == "spectrum"
    assert app["_auto"]["on"] is True and app["cmd_box"].value != stale_command

    app["all_mode"].value = "mep"
    app["w_ts"].value = False
    app["all_mode"].value = "tsonly"
    ts_only_saved = app["_session_dict"]()
    assert ts_only_saved["tsopt"] is False
    assert app["_apply_session"](ts_only_saved) == []
    app["all_mode"].value = "mep"
    assert app["w_ts"].value is False

    app["w_reuse"].value = True
    app["_invalidate_last_run"]("new system")
    assert app["w_reuse"].value is False

    missing = tmp_path / "missing-dir" / "missing.pdb"
    pending = app["_session_dict"]()
    pending.update(
        backend="mace", model="MACE-OMOL-0", inputs=[str(missing)], parm=str(topology),
        mode="pdb", subcmd="opt", all_mode="mep", center=[], center_ids=[], lcharge={},
    )
    assert app["_apply_session"](pending) == [str(missing)]
    uploaded = tmp_path / "missing.pdb"
    uploaded.write_text(primary.read_text(encoding="utf-8"), encoding="utf-8")
    metadata[str(uploaded)] = [dict(row) for row in metadata[str(primary)]]
    assert app["_ingest_saved_files"]([str(uploaded)], "session replacement")
    assert app["S"]["inputs"] == [str(uploaded)] and app["S"]["parm"] == str(topology)

    bad = tmp_path / "empty.pdb"
    bad.write_text("REMARK no atoms\nEND\n", encoding="utf-8")
    old_inputs = list(app["S"]["inputs"])
    app["_load_view_structure"] = real_loader
    assert not app["_ingest_saved_files"]([str(bad)], "invalid")
    assert app["S"]["inputs"] == old_inputs and "not attached" in app["input_msg"].value

    duplicate = {"index": 0, "chain": "A", "resn": "LIG", "resi": "10", "atom": "C1", "xyz": (0, 0, 0)}
    app["S"]["scan_atoms"] = [dict(duplicate), dict(duplicate)]
    with pytest.raises(ValueError, match="different atoms"):
        app["build_cmd"]()


def test_colab_gui_routes_scientific_options_and_round_trips_sessions() -> None:
    app = _notebook()["cells"][2]["source"]

    # SPEC / FLAG_SUBS are the single source of truth, re-derived against the
    # current mlmm CLI. `all` accepts --mep-mode, while --freeze-atoms also
    # reaches sp and dft in this repository.
    assert "SPEC = {" in app
    assert "SUBREQ = {k: v['req'] for k, v in SPEC.items()}" in app
    assert "'adv_mep':     {'all', 'path-opt', 'path-search'}," in app
    assert "'adv_dmf':     {'all', 'path-opt', 'path-search'}," in app
    assert "cmd += ['--dmf-backend', dmf_backend]" in app
    assert "'mep_mode': FLAG_SUBS['adv_mep']," in app
    assert "'threshold': FLAG_SUBS['adv_thresh']," in app
    assert "'adv_thresh':  {'all', 'opt', 'tsopt', 'scan', 'scan2d', 'scan3d', 'path-opt', 'path-search'}," in app
    assert "'adv_maxcyc':  {'all', 'opt', 'tsopt', 'irc', 'scan', 'path-opt', 'path-search'}," in app
    assert "if 'freeze' in SPEC.get(sub, {}).get('panels', ()) and S['freeze_atoms']:" in app
    # `--tr-projection legacy-active` is deprecated (it warns and must not be used
    # for pass/HOSP transition-state certification), and --embedcharge is
    # experimental. Neither may be offered in the GUI, matching pdb2reaction.
    assert "legacy-active" not in app
    assert "--tr-projection" not in app
    assert "--embedcharge" not in app
    # ONIOM utilities are classified, with their inputs/outputs stated.
    for _sub in ("'define-layer'", "'mm-parm'", "'oniom-export'", "'oniom-import'"):
        assert _sub in app
    # Surfaced key flags + the standalone cluster-model button.
    assert "key_opts_box = _collapsible('Key options', key_opts_content)" in app
    assert "cmd += ['--flatten']" in app
    assert "cmd += ['--max-cycles', str(int(mc))]" in app
    assert "b_extract = W.Button(description='Prepare ML-region model'" in app
    assert "prep_radius = W.FloatText(" in app
    # The wheel ships no examples, so Load example resolves them from the git
    # tag matching the installed release (a source checkout is used when present).
    assert "def _example_file(relpath):" in app
    assert "raw.githubusercontent.com/t-0hmura/mlmm_toolkit" in app
    assert "Run Setup first" not in app
    # -r/--radius is extraction-only and requires a structure workflow.
    assert "if radius_applies and r and r > 0: cmd += ['-r', str(r)]" in app
    assert "adv_thresh.disabled = sub not in TOOL_CAPABILITIES['threshold']" in app
    assert "elif sub == 'dft':" in app
    assert "cmd += ['--func-basis', fb]" in app
    assert "adv_mep.disabled = sub not in TOOL_CAPABILITIES['mep_mode']" in app
    assert "_set_flag_visible(adv_dmf, sub in FLAG_SUBS['adv_dmf'] and adv_mep.value == 'dmf')" in app
    assert "d['all_mode'] = _wv('all_mode', 'mep')" in app
    assert "all_mode.value = d['all_mode']" in app
    assert "def _validate_and_normalize_session(payload):" in app
    assert "_SESSION_APPLY = {'active': False}" in app
    assert "bytes(c).decode('utf-8')" in app


def test_colab_gui_preserves_full_system_and_tracks_current_run_only() -> None:
    app = _notebook()["cells"][2]["source"]

    assert "'model_pdb': None" in app
    assert "cmd += ['--model-pdb', S['model_pdb']]" in app
    assert "S['model_pdb'] = out" in app
    assert "S.update(inputs=[out]" not in app
    assert "'oniom-export': dict(n_in=(0, 1)" in app
    assert "def _clear_structure_bound_state():" in app
    for key in (
        "scan_stages", "scan_axes", "freeze_buf", "freeze_pairs",
        "freeze_atoms", "measure_atoms", "_pdb_text", "_atom_meta",
    ):
        assert key in app
    assert "_INITIAL_MODEL" in app
    assert "_bk0 = BACKEND if BACKEND in MODELS else 'mace'" in app
    assert "DFT_READY" in app and "DMF_READY" in app
    assert "options=['(default)', 'gsm', 'dmf']" not in app
    assert "OUT_JSON_SUBS" in app and "if sub in OUT_JSON_SUBS: cmd += ['--out-json']" in app
    assert "upl = W.FileUpload(accept=_acc, multiple=True, description='Upload files'" in app
    assert "_dnd_out.register_callback('mlmm_gui.on_drop', _rxgui_drop)" in app
    assert "reader.readAsDataURL(file);" in app
    assert "def _ingest_saved_files(" in app
    assert "_reset_file_upload(upl)" in app
    assert "anywidget" not in app and "_HAS_DROP_WIDGET" not in app
    assert "pts.append((k, 'TS'))" not in app
    assert "peak candidate" in app and "minimum candidate" in app
    assert "except KeyboardInterrupt:" in app
    assert "start_new_session=(os.name == 'posix')" in app
    assert "os.killpg(" in app and "signal.SIGTERM" in app and "signal.SIGKILL" in app
    assert "process.terminate()" in app and "process.kill()" in app and "process.wait(" in app
    assert "cancelled" in app
    assert "_MANUAL" in app
    assert "current_output_paths" in app and "key_output_files" in app
    assert "def _select_status_json(current, sub):" in app
    assert "if sub in ('all', 'path-search')" in app
    assert "Use machine-readable claims when present" in app
    assert "if real_run:" in app and "'exit_code': rc" in app
    assert "'status': 'success' if rc == 0" in app
    assert "Invalid command line:" in app
    assert "def _command_input_option_flags(argv):" in app
    assert "def _click_parsed_input_files(argv):" in app
    assert "isinstance(value_type, click.Path)" in app
    assert "flags.intersection({'-s', '--scan-lists'})" in app
    assert "_RUN_STATE['validation_log']" in app
    assert "tryCreateComponentFromExpression" not in app
    assert "sub == 'all' and all_kind == 'scan'" in app
    assert "bond table or JSON on stdout (--json)" in app
    assert "colab_run.json" in app and "zipfile.ZipFile" in app
    assert "shutil.make_archive(" not in app
    assert "W.Button(description='Show information'" not in app
    assert "<details class=\"rxinfo-details\" data-revision=" in app
    assert "submit(event.dataTransfer&&event.dataTransfer.files);" in app
    assert "_tab_body.children = [_TAB_PAGES[i][1]]" not in app
    assert "layout=W.Layout(width='560px')" not in app
    assert "ML-region charge (-q)" in app and "charge verified" in app
    assert "Verify the ML-region charge (-q)" in app
    assert "utility .xyz / .gjf / .com / .inp / .csv" in app
    assert app.count("effective = _normalized_scope_argv(a)") == 2
    assert "a = _force_dry_run(a)" in app
    assert "real_run = not _flag_enabled(effective, '--dry-run', '--no-dry-run')" in app


def test_colab_output_scope_executes_cli_grammar_and_utility_defaults(
    tmp_path: Path, monkeypatch,
) -> None:
    contract = _output_contract()
    scope_for = contract["_output_scope"]
    monkeypatch.chdir(tmp_path)

    from mlmm.cli import cli as product_cli
    import click

    root_context = click.Context(product_cli)
    compute_subcommands = {
        "all", "scan", "opt", "path-opt", "path-search", "tsopt", "freq",
        "irc", "dft", "sp", "scan2d", "scan3d",
    }
    for subcommand in compute_subcommands:
        command = product_cli.get_command(root_context, subcommand)
        output_param = next(
            param for param in command.params
            if {"-o", "--out", "--out-dir", "--output"}.intersection(param.opts)
        )
        assert scope_for(["mlmm", subcommand])["root"] == str(output_param.default)
    assert scope_for(["mlmm", "opt", "-oattached"])["root"] == "attached"
    assert scope_for(["mlmm", "opt", "--OUT-DIR", "Upper"])["root"] == "Upper"
    assert scope_for(["mlmm", "opt", "--OUT-DIR=UpperEq"])["root"] == "UpperEq"
    assert scope_for(["mlmm", "-i", "input.pdb"])["root"] == "result_all"

    for info_argv in (
        ["mlmm", "energy-diagram", "--help"],
        ["mlmm", "all", "--help-advanced"],
        ["mlmm", "--version"],
    ):
        info_scope = scope_for(info_argv)
        assert info_scope["stdout_only"] is True
        assert info_scope["targets"] == []
        assert not contract["_output_scope_collision"](info_scope)

    first = tmp_path / "plots-a" / "same.png"
    second = tmp_path / "plots-b" / "same.png"
    positional = tmp_path / "plots-c" / "report.csv"
    scope = scope_for([
        "mlmm", "trj2fig", str(positional), "-i", "trajectory.xyz",
        "-o", str(first), f"--out={second}",
        "--out-json", "--no-out-json", "--out-json",
    ])
    assert scope["targets"] == [
        str(first.resolve()), str(second.resolve()), str(positional.resolve()),
    ]
    assert scope["direct_current"] is True
    assert str(first.parent / "result.json") in scope["exact_targets"]
    assert str(first.parent / "summary.json") in scope["exact_targets"]
    assert len([p for p in scope["targets"] if Path(p).name == "same.png"]) == 2

    second.parent.mkdir(parents=True)
    second.write_text("stale", encoding="utf-8")
    assert contract["_output_scope_collision"](scope)
    before = contract["_snapshot_output_scope"](scope)
    first.parent.mkdir(parents=True)
    positional.parent.mkdir(parents=True)
    first.write_text("new", encoding="utf-8")
    positional.write_text("new", encoding="utf-8")
    unrelated = first.parent / "unrelated.txt"
    unrelated.write_text("ignore", encoding="utf-8")
    after = contract["_snapshot_output_scope"](scope)
    changed = {path for path, stat in after.items() if before.get(path) != stat}
    assert changed == {str(first.resolve()), str(positional.resolve())}

    default_plot = scope_for(["mlmm", "trj2fig", "-i", "trajectory.xyz"])
    assert default_plot["targets"] == [str((tmp_path / "energy.png").resolve())]
    no_json = scope_for([
        "mlmm", "trj2fig", "-i", "trajectory.xyz",
        "--out-json", "--no-out-json",
    ])
    assert not any(Path(path).name == "result.json" for path in no_json["exact_targets"])

    legacy_false = scope_for([
        "mlmm", "trj2fig", "-i", "trajectory.xyz",
        "--out-json", "False",
    ])
    assert legacy_false["targets"] == [str((tmp_path / "energy.png").resolve())]
    assert not any(Path(path).name == "False" for path in legacy_false["exact_targets"])
    assert not any(Path(path).name == "result.json" for path in legacy_false["exact_targets"])
    legacy_negative_false = scope_for([
        "mlmm", "trj2fig", "-i", "trajectory.xyz",
        "--no-out-json", "False",
    ])
    assert legacy_negative_false["targets"] == [str((tmp_path / "energy.png").resolve())]
    assert str((tmp_path / "result.json").resolve()) in legacy_negative_false["exact_targets"]
    uppercase_inline_false = scope_for([
        "mlmm", "energy-diagram", "-i", "0", "-i", "1",
        "--OUT-JSON=False",
    ])
    assert not any(Path(path).name == "result.json" for path in uppercase_inline_false["exact_targets"])

    normalize = contract["_normalized_scope_argv"]
    enabled = contract["_flag_enabled"]
    assert not enabled(normalize([
        "mlmm", "all", "--dry-run", "--no-dry-run",
    ]), "--dry-run", "--no-dry-run")
    assert not enabled(normalize([
        "mlmm", "all", "--dry-run", "False",
    ]), "--dry-run", "--no-dry-run")
    assert enabled(normalize([
        "mlmm", "all", "--no-dry-run", "False",
    ]), "--dry-run", "--no-dry-run")
    forced = contract["_force_dry_run"]([
        "mlmm", "all", "--no-dry-run", "--",
    ])
    assert forced == [
        "mlmm", "all", "--no-dry-run", "--dry-run", "--",
    ]
    assert enabled(normalize(forced), "--dry-run", "--no-dry-run")

    energy = scope_for(["mlmm", "energy-diagram", "-i", "0", "-i", "1"])
    assert energy["targets"] == [str((tmp_path / "energy_diagram.png").resolve())]
    fake_home = tmp_path / "home"
    fake_home.mkdir()
    monkeypatch.setenv("HOME", str(fake_home))
    tilde_energy = scope_for([
        "mlmm", "energy-diagram", "-i", "0", "-i", "1", "-o", "~/profile",
    ])
    assert tilde_energy["targets"] == [str((tmp_path / "~" / "profile.png").resolve())]
    tilde_trj = scope_for([
        "mlmm", "trj2fig", "-i", "trajectory.xyz", "-o", "~/profile.png",
    ])
    assert tilde_trj["targets"] == [str((fake_home / "profile.png").resolve())]
    add_elem = scope_for(["mlmm", "add-elem-info", "-i", "inputs/enzyme.pdb"])
    assert add_elem["targets"] == [str((tmp_path / "inputs/enzyme_add_elem.pdb").resolve())]

    input_file = tmp_path / "inputs/enzyme.pdb"
    input_file.parent.mkdir(parents=True)
    input_file.write_text("ATOM\n", encoding="utf-8")
    assert not contract["_output_scope_collision"](add_elem)
    add_elem_inplace = scope_for([
        "mlmm", "add-elem-info", "-i", str(input_file), "--inplace",
    ])
    assert add_elem_inplace["targets"] == [str(input_file.resolve())]
    assert contract["_output_scope_collision"](add_elem_inplace)
    fixed = scope_for(["mlmm", "fix-altloc", "-i", str(input_file)])
    assert fixed["targets"] == [str(input_file.with_name("enzyme_clean.pdb").resolve())]
    inplace = scope_for(["mlmm", "fix-altloc", "-i", str(input_file), "--inplace"])
    assert set(inplace["targets"]) == {
        str(input_file.resolve()), str(input_file.with_suffix(".pdb.bak").resolve()),
    }

    layered = scope_for([
        "mlmm", "define-layer", "-i", "inputs/enzyme.cif", "--model-indices", "1",
    ])
    assert layered["targets"] == [
        str((tmp_path / "inputs/enzyme_layered.pdb").resolve())
    ]
    assert str((tmp_path / "inputs/enzyme_layered.cif").resolve()) in layered["exact_targets"]
    layered_pdb = scope_for([
        "mlmm", "define-layer", "-i", "inputs/enzyme.pdb", "--model-indices", "1",
    ])
    assert layered_pdb["targets"] == [
        str((tmp_path / "inputs/enzyme_layered.pdb").resolve())
    ]
    assert str((tmp_path / "inputs/enzyme_layered.cif").resolve()) not in layered_pdb["exact_targets"]
    parm = scope_for(["mlmm", "mm-parm", "-i", "inputs/enzyme.pdb"])
    assert parm["targets"] == [
        str((tmp_path / "enzyme.parm7").resolve()),
        str((tmp_path / "enzyme.rst7").resolve()),
    ]
    parm_h = scope_for([
        "mlmm", "mm-parm", "-i", "inputs/enzyme.pdb", "--add-h",
    ])
    assert str((tmp_path / "enzyme_parm.pdb").resolve()) in parm_h["targets"]
    explicit_prefix = scope_for([
        "mlmm", "mm-parm", "-i", "inputs/enzyme.pdb", "-o", "prepared/system",
    ])
    assert set(explicit_prefix["targets"]) == {
        str((tmp_path / "prepared/system.parm7").resolve()),
        str((tmp_path / "prepared/system.rst7").resolve()),
        str((tmp_path / "prepared/system.pdb").resolve()),
    }
    imported = scope_for(["mlmm", "oniom-import", "-i", "inputs/job.log"])
    assert imported["targets"] == [
        str((tmp_path / "job.xyz").resolve()),
        str((tmp_path / "job_layered.pdb").resolve()),
    ]
    tilde_import = scope_for([
        "mlmm", "oniom-import", "-i", "inputs/job.log", "-o", "~/job",
    ])
    assert tilde_import["targets"] == [
        str((tmp_path / "~" / "job.xyz").resolve()),
        str((tmp_path / "~" / "job_layered.pdb").resolve()),
    ]
    dot_import = scope_for([
        "mlmm", "oniom-import", "-i", "inputs/job.log", "-o", ".",
    ])
    assert dot_import["targets"] == [
        str(tmp_path.with_suffix(".xyz")),
        str(tmp_path.parent / (tmp_path.name + "_layered.pdb")),
    ]
    extract_default = scope_for([
        "mlmm", "extract", "-i", "reactant.pdb", "-c", "LIG",
    ])
    assert extract_default["targets"] == [str((tmp_path / "pocket.pdb").resolve())]
    assert str((tmp_path / "pocket.cif").resolve()) not in extract_default["exact_targets"]
    extract_single_extra_outputs = scope_for([
        "mlmm", "extract", "-i", "reactant.pdb", "-c", "LIG",
        "-o", "first.pdb", "ignored.pdb",
    ])
    assert extract_single_extra_outputs["targets"] == [
        str((tmp_path / "first.pdb").resolve())
    ]
    extract_cif = scope_for([
        "mlmm", "extract", "-i", "reactant.cif", "-c", "LIG",
        "-o", "cluster.ent", "--out-json",
    ])
    assert extract_cif["targets"] == [str((tmp_path / "cluster.ent").resolve())]
    assert str((tmp_path / "cluster.cif").resolve()) in extract_cif["exact_targets"]
    duplicate_input = scope_for([
        "mlmm", "extract", "-i", "same.pdb", "same.pdb", "-c", "LIG",
    ])
    assert duplicate_input["targets"] == [str((tmp_path / "pocket_same.pdb").resolve())]
    mixed = scope_for([
        "mlmm", "extract", "-i", "normal.pdb", "wide.cif", "-c", "LIG",
        "-o", "normal.ent", "wide.ent",
    ])
    assert str((tmp_path / "normal.cif").resolve()) not in mixed["exact_targets"]
    assert str((tmp_path / "wide.cif").resolve()) in mixed["exact_targets"]
    combined_first_regular = scope_for([
        "mlmm", "extract", "-i", "normal.pdb", "wide.cif", "-c", "LIG",
        "-o", "combined.ent",
    ])
    assert str((tmp_path / "combined.cif").resolve()) not in combined_first_regular["exact_targets"]

    (tmp_path / "result").mkdir()
    (tmp_path / "result" / "stale.txt").write_text("stale", encoding="utf-8")
    stdout_scope = scope_for([
        "mlmm", "bond-summary", "-i", "reactant.pdb", "product.pdb",
    ])
    assert stdout_scope["stdout_only"] is True
    assert contract["_snapshot_output_scope"](stdout_scope) == {}
    assert contract["_output_scope_collision"](stdout_scope) is False

    exported = scope_for([
        "mlmm", "oniom-export", "--parm", "system.parm7", "-o", "job",
        "--mode", "g16", "-q", "0",
    ])
    assert exported["targets"] == [str((tmp_path / "job").resolve())]
    exported_orca = scope_for([
        "mlmm", "oniom-export", "--parm", "system.parm7", "-o", "job.inp",
        "--mode", "orca", "-q", "0",
    ])
    assert set(exported_orca["targets"]) == {
        str((tmp_path / "job.inp").resolve()),
        str((tmp_path / "system.ORCAFF.prms").resolve()),
    }
