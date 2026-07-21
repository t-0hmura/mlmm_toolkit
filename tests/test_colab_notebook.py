"""Static contracts for the release-matched Colab GUI notebook."""

from __future__ import annotations

import ast
import json
from pathlib import Path


NOTEBOOK = Path(__file__).parents[1] / "examples" / "mlmm_colab.ipynb"


def _notebook() -> dict:
    return json.loads(NOTEBOOK.read_text(encoding="utf-8"))


def test_colab_notebook_has_valid_code_cells_and_gpu_metadata() -> None:
    notebook = _notebook()

    assert notebook["nbformat"] == 4
    assert notebook["metadata"]["accelerator"] == "GPU"
    assert len(notebook["cells"]) == 4
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            ast.parse(cell["source"])


def test_colab_setup_is_pinned_to_matching_release_and_one_backend() -> None:
    setup = _notebook()["cells"][1]["source"]

    assert 'mlmm_ref = "v0.3.3"' in setup
    assert "checkout','--detach'" in setup
    assert "installed_version != mlmm_ref[1:]" in setup
    assert "version('mlmm-toolkit')" in setup
    assert "Restart the Colab runtime first" in setup
    assert "mace-torch>=0.3.8" in setup
    assert "HF_TOKEN" in setup


def test_colab_gui_is_mlmm_native_and_tracks_structure_contracts() -> None:
    app = _notebook()["cells"][2]["source"]
    whole = NOTEBOOK.read_text(encoding="utf-8").lower()

    assert "pdb2reaction" not in whole
    assert "p2r" not in whole
    assert ".pdb,.ent,.cif,.mmcif,.parm7" in app
    assert ".rst7" not in app
    assert "prepare_input_structure" in app
    assert "load_pdb_atom_metadata" in app
    assert "mlmm_gui_viewer_input.pdb" in app
    assert app.count("mlmm_gui.on_click") == 2
    assert "CHAIN:RESNAME:RESSEQ" in app
    assert "String(atom.index)" in app
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
    assert "flags = ('-o', '--out-dir', '--output')" in app
    assert "out = _effective_out_dir(a)" in app
    assert "os.path.isdir(out)" in app
    assert "_results(out)" in app
    assert "reuse non-empty out dir" in app


def test_colab_gui_keeps_responsive_release_layout() -> None:
    app = _notebook()["cells"][2]["source"]

    assert "@media (max-width: 600px)" in app
    assert "max_width='100%'" in app
    assert "No structure loaded" in app
    # Colab renders ipywidgets' Tab and Accordion as empty blocks, so the tab
    # strip is Buttons + a swapping VBox and every collapsible is Button + VBox.
    assert "_TAB_PAGES = [('1 Input', input_box), ('2 Select', select_box)," in app
    assert "('3 Options', options_box), ('4 Results', results_box)]" in app
    assert "def _tab_go(i):" in app
    assert "W.Tab(" not in app
    assert "def _collapsible(title, child, on_open=None):" in app
    assert "W.Accordion(" not in app
    # Hover help is the HTML `title` global attribute with a ⓘ affordance.
    assert "def _hdr(html, tip):" in app
    assert "&#9432;" in app
    assert "rxworkspace" in app
    assert "rxviewer" in app
    assert "rxinspector" in app


def test_colab_gui_routes_scientific_options_and_round_trips_sessions() -> None:
    app = _notebook()["cells"][2]["source"]

    # SPEC / FLAG_SUBS are the single source of truth, re-derived against the
    # mlmm CLI: mlmm's `all` does NOT accept --mep-mode (p2r's does), and
    # --freeze-atoms reaches sp and dft here.
    assert "SPEC = {" in app
    assert "SUBREQ = {k: v['req'] for k, v in SPEC.items()}" in app
    assert "'adv_mep':     {'path-opt', 'path-search'}," in app
    assert "'mep_mode': FLAG_SUBS['adv_mep']," in app
    assert "'threshold': FLAG_SUBS['adv_thresh']," in app
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
    assert "key_opts_box = W.VBox([" in app
    assert "cmd += ['--flatten']" in app
    assert "cmd += ['--max-cycles', str(int(mc))]" in app
    assert "b_extract = W.Button(description='Extract cluster model'" in app
    # -r/--radius is extraction-only; scan must not receive it.
    assert "if sub in ('all', 'extract') and r and r > 0: cmd += ['-r', str(r)]" in app
    assert "sub in TOOL_CAPABILITIES['threshold']" in app
    assert "elif sub == 'dft':" in app
    assert "cmd += ['--func-basis', fb]" in app
    assert "adv_mep.disabled = sub not in TOOL_CAPABILITIES['mep_mode']" in app
    assert "d['all_mode'] = _wv('all_mode', 'mep')" in app
    assert "all_mode.value = saved_all_mode" in app
    assert "bytes(c).decode('utf-8')" in app
