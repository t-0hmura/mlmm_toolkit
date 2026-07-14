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
    assert "['1 Input', '2 Select', '3 Options', '4 Results']" in app
    assert "rxworkspace" in app
    assert "rxviewer" in app
    assert "rxinspector" in app
