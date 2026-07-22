"""Static contracts for the release-matched Colab GUI notebook."""

from __future__ import annotations

import ast
import glob
import json
import os
import shlex
from pathlib import Path


NOTEBOOK = Path(__file__).parents[1] / "examples" / "mlmm_colab.ipynb"


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
    # The release notebook installs the pinned wheel from PyPI, the same way a
    # normal user does, so the version guard below compares the version actually
    # resolved by pip against the requested tag.
    assert "pip('mlmm-toolkit' + ('[dft]' if install_dft else '') + '==' + mlmm_ref.lstrip('v'))" in setup
    assert "git clone" not in setup
    assert "installed_version != mlmm_ref[1:]" in setup
    assert "version('mlmm-toolkit')" in setup
    assert "Restart the Colab runtime first" in setup
    assert "mace-torch>=0.3.8" in setup
    assert "HF_TOKEN" in setup
    assert "install_dft = False" in setup
    assert "INSTALL_DFT = install_dft" in setup


def test_colab_gui_is_mlmm_native_and_tracks_structure_contracts() -> None:
    app = _notebook()["cells"][2]["source"]
    whole = NOTEBOOK.read_text(encoding="utf-8").lower()

    assert isinstance(app, str)
    assert "pdb2reaction" not in whole
    assert "p2r" not in whole
    assert ".pdb,.ent,.cif,.mmcif,.parm7" in app
    assert "else '.pdb,.ent,.cif,.mmcif,.parm7'" in app
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
    assert "No structure loaded" in app
    assert "('Scan (1 input)', 'scan')" in app
    assert "Run a workflow, then choose <b>Show results</b>." in app
    assert "frame_slider.disabled = False" in app
    assert "trajectory_box.layout.display = 'none'" in app
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
    assert "key_opts_box = W.VBox([" in app
    assert "cmd += ['--flatten']" in app
    assert "cmd += ['--max-cycles', str(int(mc))]" in app
    assert "b_extract = W.Button(description='Prepare ML-region model'" in app
    # The wheel ships no examples, so Load example resolves them from the git
    # tag matching the installed release (a source checkout is used when present).
    assert "def _example_file(relpath):" in app
    assert "raw.githubusercontent.com/t-0hmura/mlmm_toolkit" in app
    assert "Run Setup first" not in app
    # -r/--radius is extraction-only; scan must not receive it.
    assert "if sub in ('all', 'extract') and r and r > 0: cmd += ['-r', str(r)]" in app
    assert "sub in TOOL_CAPABILITIES['threshold']" in app
    assert "elif sub == 'dft':" in app
    assert "cmd += ['--func-basis', fb]" in app
    assert "adv_mep.disabled = sub not in TOOL_CAPABILITIES['mep_mode']" in app
    assert "adv_dmf.layout.display = ('' if sub in FLAG_SUBS['adv_dmf'] and adv_mep.value == 'dmf'" in app
    assert "d['all_mode'] = _wv('all_mode', 'mep')" in app
    assert "all_mode.value = saved_all_mode" in app
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
    assert "out = new Array(fl.length)" in app and "out[index]" in app
    assert "out.push(" not in app
    assert "r.onerror = function()" in app
    assert "pts.append((k, 'TS'))" not in app
    assert "peak candidate" in app and "minimum candidate" in app
    assert "except KeyboardInterrupt:" in app
    assert "start_new_session=(os.name == 'posix')" in app
    assert "os.killpg(" in app and "signal.SIGTERM" in app and "signal.SIGKILL" in app
    assert "p.terminate()" in app and "p.kill()" in app and "p.wait(" in app
    assert "cancelled" in app
    assert "_MANUAL" in app
    assert "current_output_paths" in app and "key_output_files" in app
    assert "def _select_status_json(current, sub):" in app
    assert "if sub in ('all', 'path-search')" in app
    assert "Use machine-readable claims when present" in app
    assert "if real_run:" in app and "'exit_code': rc" in app
    assert "'status': 'success' if rc == 0" in app
    assert "Invalid command line:" in app
    assert "sub == 'all' and all_kind == 'scan'" in app
    assert "bond table or JSON on stdout (--json)" in app
    assert "colab_run.json" in app and "zipfile.ZipFile" in app
    assert "shutil.make_archive(" not in app
    assert 'tabindex="0"' in app and "aria-label=" in app
    assert "layout=W.Layout(width='560px')" not in app
    assert "ML-region charge (-q)" in app and "charge verified" in app
    assert "Verify the ML-region charge (-q)" in app
    assert "ML/MM drag & drop accepts PDB/mmCIF + parm7" in app
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
    assert add_elem["targets"] == [str((tmp_path / "inputs/enzyme.pdb").resolve())]

    input_file = tmp_path / "inputs/enzyme.pdb"
    input_file.parent.mkdir(parents=True)
    input_file.write_text("ATOM\n", encoding="utf-8")
    assert contract["_output_scope_collision"](add_elem)
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
