"""Contracts for the release-matched Colab GUI notebook: source-level
assertions plus cells compiled and executed against real ipywidgets."""

from __future__ import annotations

import ast
import base64
import csv
import datetime
import glob
import gzip
import html
import json
import os
import re
import shlex
import subprocess
import sys
import types
import zipfile
from pathlib import Path

import pytest


NOTEBOOK = Path(__file__).parents[1] / "examples" / "mlmm_colab.ipynb"


def _embedded_document(frame: str) -> str:
    """Decode either the direct or large-payload iframe contract."""
    packed = re.search(r'data-rx-document="([^"]+)"', frame)
    if packed is not None:
        return gzip.decompress(base64.b64decode(packed.group(1))).decode()
    srcdoc = re.search(r'srcdoc="([^"]*)"', frame)
    assert srcdoc is not None
    return html.unescape(srcdoc.group(1))


def _execute_app(monkeypatch, tmp_path: Path, *, parent_header: dict | None = None) -> tuple[dict, list]:
    """Execute the complete app cell with real widgets and captured HTML."""
    rendered: list = []
    import IPython.display as ipd
    monkeypatch.setattr(ipd, "display", lambda *args, **kwargs: rendered.extend(args))
    monkeypatch.setattr(ipd, "clear_output", lambda *args, **kwargs: None)
    monkeypatch.setattr(ipd, "HTML", lambda value: value)
    monkeypatch.setattr(ipd, "Image", lambda *args, **kwargs: (args, kwargs))
    monkeypatch.chdir(tmp_path)
    # A Colab-compatible local runtime may have google.colab installed without
    # the native Colab kernel bridges.  Keep ordinary GUI tests on the local
    # path; tests that intentionally exercise Colab install fake modules first.
    colab_module = sys.modules.get("google.colab")
    if colab_module is None or getattr(colab_module, "__file__", None):
        for module_name in (
            "google.colab.output", "google.colab.files", "google.colab.userdata"
        ):
            module = sys.modules.get(module_name)
            if module is None or getattr(module, "__file__", None):
                monkeypatch.delitem(sys.modules, module_name, raising=False)
        monkeypatch.setitem(sys.modules, "google.colab", None)
    namespace = {"TOOL": "mlmm", "BACKEND": "mace", "REPO_DIR": "unused"}
    if parent_header is not None:
        namespace["get_ipython"] = lambda: types.SimpleNamespace(parent_header=parent_header)
    source = _notebook()["cells"][2]["source"]
    exec(compile(source, str(NOTEBOOK), "exec"), namespace)
    return namespace, rendered


def _notebook() -> dict:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    # nbformat stores `source` either as one string or as a list of lines, and
    # Colab exports the list form.  Normalize here so every contract below can
    # treat a cell's source as plain text.
    for cell in notebook["cells"]:
        source = cell.get("source")
        if isinstance(source, list):
            cell["source"] = "".join(source)
    return notebook


def test_colab_debug2_result_workflow_is_integrated_without_regressions() -> None:
    notebook = _notebook()
    setup = notebook["cells"][1]["source"]
    app = notebook["cells"][2]["source"]

    assert "_plot_probe_code" in setup
    assert "'show_water': True" in app
    assert "Toy system - MEP mode (R->P)" in app
    assert "toy_system/r_toy.pdb" in app
    assert "toy_system/p_toy.pdb" in app
    assert "toy_system/p_toy.parm7" in app
    assert "repeat=True, show_repeat=True" in app
    assert "role = _irc_trajectory_role(path)" in app
    assert "if 'results_dir' in globals(): results_dir.value = str(out)" in app
    assert "if name == 'summary.log': return -1" in app
    assert "trajectory_box, res_out, results_empty, result_context, artifact_fold" in app
    assert "Return every distinct imaginary-mode trajectory for each TS." in app
    assert "one TS imaginary-mode animation" not in app


def _css_rule_has(source: str, selector: str, *declarations: str) -> bool:
    """Match declarations in one selector block without coupling to whitespace."""
    compact = lambda value: re.sub(r"\s+", "", value)
    expected = tuple(compact(declaration) for declaration in declarations)
    rules = re.findall(rf"{re.escape(selector)}\s*\{{([^{{}}]*)\}}", source)
    return any(all(declaration in compact(rule) for declaration in expected) for rule in rules)


def _css_media_body(source: str, condition: str) -> str:
    """Return a balanced @media body so responsive rules are checked in scope."""
    marker = f"@media ({condition})"
    start = source.index(marker)
    opening = source.index("{", start + len(marker))
    depth = 0
    for index in range(opening, len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[opening + 1:index]
    raise AssertionError(f"unterminated CSS media block: {condition}")


def test_installation_backend_guidance_is_explicit() -> None:
    setup = _notebook()["cells"][1]["source"]
    expected = (
        "Recommended backend: UMA (Hugging Face licence acceptance + read token required). "
        "Without a token, use ORB or MACE. Restart the runtime before switching backends. "
        "UMA uses fp32; ORB/MACE use fp64 by default."
    )
    assert expected in setup
    assert "Recommended: UMA" not in setup


def _set_file_upload(widget, name: str, media_type: str, content: bytes) -> None:
    """Deliver one browser upload through the widget's public trait mechanism."""
    metadata = {
        "name": name,
        "type": media_type,
        "size": len(content),
        "last_modified": datetime.datetime.now(datetime.timezone.utc),
    }
    if hasattr(widget, "_counter"):  # ipywidgets 7
        widget.set_trait("value", {name: {"metadata": metadata, "content": content}})
    else:  # ipywidgets 8
        widget.set_trait("value", ({**metadata, "content": memoryview(content)},))


def test_setup_command_fields_do_not_repaint_hidden_results_or_viewer(monkeypatch, tmp_path: Path) -> None:
    """Command-only Setup edits must not recreate Colab's live Mol* iframe."""
    app, _ = _execute_app(monkeypatch, tmp_path)

    class OutputSpy:
        def __init__(self) -> None:
            self.writes = 0
            self._outputs = ("keep",)

        @property
        def outputs(self):
            return self._outputs

        @outputs.setter
        def outputs(self, value) -> None:
            self.writes += 1
            self._outputs = value

    result_spy = OutputSpy()
    app["res_out"] = result_spy
    calls = {"summary": 0, "chips": 0, "output": 0}
    app["_render_summary"] = lambda: calls.__setitem__("summary", calls["summary"] + 1)
    app["_render_chips"] = lambda: calls.__setitem__("chips", calls["chips"] + 1)
    app["_render_output_note"] = lambda: calls.__setitem__("output", calls["output"] + 1)
    app["_VIEWER_MOUNTED"]["value"] = True
    app["_VIEWER_GENERATION"]["value"] = 11
    app["viewer_out"].value = '<iframe class="rxmolstar-frame">stable</iframe>'
    app["S"]["mode"] = "pdb"
    app["S"]["inputs"] = [str(tmp_path / "input.pdb")]
    app["_sync_capability_controls"]()
    assert "selectedresn" in dict(app["pick_action"].options).values()
    assert app["prep_radius"].step == 0.5
    assert app["adv_radius"].step == 0.5
    assert app["prep_radius_control"] is app["prep_radius"]
    assert app["adv_radius_control"] is app["adv_radius"]
    assert not hasattr(app["prep_radius_control"], "_rx_step_buttons")
    app["prep_radius"].value += app["prep_radius"].step
    assert app["prep_radius"].value == pytest.approx(3.1)
    assert app["adv_radius"].value == pytest.approx(3.1)
    app["prep_radius"].value -= app["prep_radius"].step
    assert app["prep_radius"].value == pytest.approx(2.6)
    assert app["prep_radius"].min == 0.0
    assert app["b_done_selected_resn"] is app["b_pick_selected_resn"]
    assert app["b_pick_selected_resn"].description == "Pick force-included residues"
    assert app["b_pick_selected_resn"].icon == "mouse-pointer"
    assert app["b_pick_selected_resn"].button_style == ""

    app["selected_resn"].value = "A:123"
    app["prep_radius"].value = 0.0

    assert result_spy.writes == 0
    assert calls == {"summary": 1, "chips": 1, "output": 0}
    assert app["_VIEWER_GENERATION"]["value"] == 11
    assert "stable" in app["viewer_out"].value
    assert "--selected-resn A:123" in app["cmd_box"].value
    assert "-r 0" in app["cmd_box"].value
    app["selected_resn"].value = ""

    app["S"]["_atom_meta"] = [{
        "index": 0, "serial": 1, "chain": "A", "resname": "HIS",
        "resseq": 123, "icode": "", "name": "CA", "xyz": (0.0, 0.0, 0.0),
    }]
    app["b_pick_selected_resn"].click()
    assert app["pick_action"].value == "selectedresn"
    assert app["b_pick_selected_resn"].description == "Done picking"
    assert app["b_pick_selected_resn"].icon == "check"
    assert app["b_pick_selected_resn"].button_style == "primary"
    assert app["b_pick_selected_resn"].layout.display != "none"
    assert "rxactive-card" in app["center_panel"]._dom_classes
    assert "rxactive-card" not in app["charge_panel"]._dom_classes
    # Mol* sequence residue-loci arrive as non-exact clicks and use the same
    # selected-residue toggle as a 3D atom click.
    app["on_click"]("0", live_marked=True, exact=False)
    assert app["selected_resn"].value == "123"
    assert app["_VIEWER_GENERATION"]["value"] == 11
    app["on_click"]("0", live_marked=True, exact=False)
    assert app["selected_resn"].value == ""
    assert app["_VIEWER_GENERATION"]["value"] == 11
    app["b_pick_selected_resn"].click()
    assert app["pick_action"].value == "center"
    assert app["b_pick_selected_resn"].description == "Pick force-included residues"
    assert app["b_pick_selected_resn"].icon == "mouse-pointer"
    assert app["b_pick_selected_resn"].button_style == ""

    def freeze_picker_button(description):
        stack = [app["freeze_panel"]]
        while stack:
            widget = stack.pop()
            if getattr(widget, "description", None) == description:
                return widget
            stack.extend(getattr(widget, "children", ()))
        raise AssertionError("freeze picker button not found: " + description)

    freeze_picker_button("Pick frozen atoms").click()
    assert app["pick_action"].value == "freezeatom"
    done = freeze_picker_button("Done picking")
    assert done.icon == "check"
    assert done.button_style == "primary"
    done.click()
    assert app["pick_action"].value == "center"
    assert freeze_picker_button("Pick frozen atoms").icon == "mouse-pointer"

    app["S"]["mode"] = "small"
    app["pick_action"].value = "freezeatom"
    app["_sync_active_selection_card"]()
    assert "rxactive-card" not in app["freeze_panel"]._dom_classes
    assert "rxactive-card" in app["system_charge_panel"]._dom_classes
    app["S"]["mode"] = "pdb"
    app["_sync_active_selection_card"]()
    assert "rxactive-card" in app["freeze_panel"]._dom_classes
    assert "rxactive-card" not in app["system_charge_panel"]._dom_classes
    assert "0 件選択" not in _notebook()["cells"][2]["source"]

    manual = app["cmd_box"].value + " --manual-flag"
    app["cmd_box"].value = manual
    assert app["manual_mode_notice"].layout.display == ""
    app["prep_radius"].value = 3.3
    assert app["cmd_box"].value == manual
    app["b_rebuild"].click()
    assert app["manual_mode_notice"].layout.display == "none"
    assert "-r 3.3" in app["cmd_box"].value


def _root_normalized_subcommand_argv(app: dict, subcommand: str, argv: list[str]) -> list[str]:
    """Apply the root CLI's compatibility normalization before direct parsing."""
    root = app["PRODUCT_CLI"]
    root_context = app["click"].Context(root)
    bool_values, bool_toggles, negative_aliases, single_flags = (
        root._resolve_bool_options(root_context, subcommand)
    )
    normalized, _legacy = root._normalize_bool_argv(
        [subcommand, *argv],
        {subcommand: bool_values},
        {subcommand: bool_toggles},
        {subcommand: negative_aliases},
        {subcommand: single_flags},
    )
    return normalized[1:]


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
        "_irc_trajectory_role",
        "_trajectory_result_metadata",
        "_trajectory_segment_ranges",
        "_trajectory_semantics",
        "_trajectory_frame_context",
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
            ".gjf": "text", ".com": "text", ".inp": "text", ".prms": "text",
            ".pdb": "structure", ".ent": "structure", ".cif": "structure",
            ".mmcif": "structure",
        },
    }
    exec(compile(module, str(NOTEBOOK), "exec"), namespace)
    assert wanted <= namespace.keys()
    return namespace


def test_small_view_loader_reads_xyz_ordinary_gaussian_and_oniom(
    tmp_path: Path,
) -> None:
    source = _notebook()["cells"][2]["source"]
    tree = ast.parse(source)
    loader = next(
        node for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "_load_small_view_structure"
    )
    namespace = {
        "Path": Path,
        "_runtime_path": lambda name: str(tmp_path / name),
    }
    exec(
        compile(ast.Module(body=[loader], type_ignores=[]), str(NOTEBOOK), "exec"),
        namespace,
    )

    xyz = tmp_path / "input.xyz"
    xyz.write_text("2\nXYZ\nC 0 0 0\nH 0 0 1\n", encoding="utf-8")
    ordinary = tmp_path / "ordinary.gjf"
    ordinary.write_text(
        "#p hf/3-21g\n\nordinary\n\n0 1\nC 0 0 0\nH 0 0 1\n\n",
        encoding="utf-8",
    )
    oniom = tmp_path / "oniom.com"
    oniom.write_text(
        "# oniom(hf/3-21g:amber)\n\noniom\n\n0 1 0 1 0 1\n"
        "C-CT 0 0 0 0 H\nH-HC 0 0 0 1 L\n\n",
        encoding="utf-8",
    )

    for path in (xyz, ordinary, oniom):
        text, metadata, viewer_path = namespace["_load_small_view_structure"](path)
        assert text.splitlines()[0] == "2"
        assert [atom["element"] for atom in metadata] == ["C", "H"]
        assert [atom["index"] for atom in metadata] == [0, 1]
        assert Path(viewer_path).read_text(encoding="utf-8") == text


def _advanced_widget_values(app: dict, param) -> list:
    """Return browser-widget values that exercise each option serializer."""
    click = app["click"]
    if param.name == "verbose":
        return [0, 3]
    if param.is_bool_flag or isinstance(param.type, click.types.BoolParamType):
        return [True, False]
    if isinstance(param.type, click.Choice):
        choices = list(param.type.choices)
        assert choices, f"{param.name} has an empty Click choice"
        return [choices[0]]
    if isinstance(param.type, (click.types.IntParamType, click.types.FloatParamType)):
        default = param.default
        lower = getattr(param.type, "min", None)
        upper = getattr(param.type, "max", None)
        step = 1 if isinstance(param.type, click.types.IntParamType) else 0.5
        candidate = ((lower + step) if lower is not None else step) if default is None else default + step
        if upper is not None and candidate > upper:
            candidate = default - step
        if lower is not None and candidate < lower:
            candidate = lower
        assert candidate != default, f"{param.name} has no editable numeric value"
        return [candidate]
    nargs = max(1, int(getattr(param, "nargs", 1)))
    groups = 2 if param.multiple else 1
    return [" ".join(str(index + 1) for index in range(nargs * groups))]


def _widget_descendants(widget):
    yield widget
    for child in getattr(widget, "children", ()):
        yield from _widget_descendants(child)


def _widget_with_description(widget, description: str, *, enabled: bool = True):
    matches = [
        child for child in _widget_descendants(widget)
        if getattr(child, "description", "") == description
        and (not enabled or not getattr(child, "disabled", False))
    ]
    assert matches, f"missing widget: {description}"
    return matches[0]


def _assert_advanced_argv_parses(app: dict, subcommand: str, argv: list[str]) -> None:
    """Use Click's option parser without invoking a scientific workflow."""
    click = app["click"]
    command = app["_advanced_command"](subcommand)
    assert command is not None
    parser = command.make_parser(click.Context(command, info_name=subcommand))
    _values, leftovers, _order = parser.parse_args(list(argv))
    assert leftovers == [], (subcommand, argv, leftovers)


def test_colab_notebook_has_valid_code_cells_and_gpu_metadata() -> None:
    notebook = _notebook()
    introduction = notebook["cells"][0]["source"]

    assert notebook["nbformat"] == 4
    assert notebook["metadata"]["accelerator"] == "GPU"
    assert len(notebook["cells"]) == 3
    assert (
        "[t-0hmura/mlmm_toolkit](https://github.com/t-0hmura/mlmm_toolkit)"
    ) in introduction
    assert "ML/MM ONIOM reaction paths on full solvated proteins" in introduction
    assert "**① Input → ② Setup → ③ Options → (Run) → ④ Results**" in introduction
    assert (
        "[DOI: 10.26434/chemrxiv-2025-jft1k]"
        "(https://doi.org/10.26434/chemrxiv-2025-jft1k)"
    ) in introduction
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            ast.parse(cell["source"])


def test_colab_responsive_contract_covers_the_whole_four_step_app() -> None:
    app = _notebook()["cells"][2]["source"]

    for breakpoint in ("940px", "720px", "600px", "480px"):
        assert f"@media (max-width: {breakpoint})" in app or (
            breakpoint == "600px" and "@media (max-width: 600px)" in app
        )
    assert ".rxviewer, .rxinspector { flex:1 1 100%" in app
    assert ".rxpath-grid { grid-template-columns:1fr; }" in app
    assert ".rxviewer-toolbar { flex-wrap:wrap" in app
    assert ".rxcommand-actions { display:grid" in app
    assert ".rxapp .rxfile { flex-wrap:wrap" in app
    assert "height:clamp(280px,52vh,420px)" in app
    assert "height:clamp(340px,72vw,460px)" in app


def test_colab_setup_is_pinned_to_matching_release_and_one_backend() -> None:
    setup = _notebook()["cells"][1]["source"]

    assert 'mlmm_toolkit_version = "v0.3.3"' in setup
    # Only the exact `debug` sentinel switches to the matching adjacent source
    # snapshot; release versions never inspect the ZIP.
    assert "_raw_version = str(mlmm_toolkit_version)" in setup
    assert "_debug_install = _raw_version == 'debug'" in setup
    assert "_raw_version.lower() == 'debug'" not in setup
    # One pinned install covers both cases: the `[dft]` extra is selected inside
    # the same expression, so the requested version cannot differ between the
    # plain and DFT paths.
    assert (
        "_release_spec = 'mlmm-toolkit%s==%s' % ('[dft]' if install_dft else '', _requested_version)"
        in setup
    )
    assert "pip(_release_spec)" in setup
    # The [dft] extra goes through the same quiet `pip` helper as every other
    # install; the streaming `pip_logged` variant was removed.
    assert "pip_logged" not in setup
    assert "install_dft is ticked" in setup
    # Gated UMA sign-in is the last step, so no install phase waits on a prompt.
    assert "Hugging Face sign-in runs at the end of Installation" in setup
    assert setup.index("notebook_login()") > setup.index("_phase_done('version verified')")
    assert "except KeyboardInterrupt:" in setup
    assert "type(_hf_login_exc).__name__ != 'DeviceCodeError'" in setup
    assert "UMA authentication is still pending" in setup
    assert "git clone" not in setup
    assert "installed_version != _requested_version" in setup
    assert "anywidget==0.11.0" in setup
    assert "version('mlmm-toolkit')" in setup
    assert "Restart the Colab runtime first" in setup
    assert "mace-torch==0.3.16" in setup
    assert "HF_TOKEN" in setup
    assert 'install_dft = True  #@param {type:"boolean"}' in setup
    assert "INSTALL_DFT = install_dft" in setup
    assert "installs ambertools/openmm dependencies, mlmm-toolkit, the selected mlip backend" in setup.lower()
    assert "Only the **selected backend** is installed" not in setup
    assert "DFT_SETUP_READY = True" in setup
    assert "coinor-libipopt-dev" in setup
    assert "pip('cyipopt')" in setup
    assert "_missing_dmf" in setup
    assert "_dft_packages = {'pyscf': 'pyscf', 'gpu4pyscf': 'gpu4pyscf-cuda12x'}" in setup
    source_id_match = re.search(r"SOURCE_BUNDLE_ID = '([0-9a-f]{64})'", setup)
    assert source_id_match is not None
    source_marker = NOTEBOOK.parents[1] / ".colab-debug-source"
    if source_marker.is_file():
        assert source_id_match.group(1) == source_marker.read_text(
            encoding="utf-8"
        ).strip()
    assert "pip('-e', './' + REPO_DIR + ('[dft]' if install_dft else ''))" in setup
    for rejection in (
        "has no source-snapshot marker",
        "Notebook and ZIP are from different debug builds",
        "contains a file outside its source directory",
        "contains an unsafe path",
        "contains an unsupported symbolic link",
        "did not extract the expected source snapshot",
    ):
        assert rejection in setup


def test_colab_setup_installs_missing_cyipopt(monkeypatch) -> None:
    import importlib.metadata

    setup = _notebook()["cells"][1]["source"].replace(
        "install_dft = True", "install_dft = False", 1
    )
    calls: list[list[str]] = []
    cyipopt_installed = False
    content_checks = 0

    def fake_run(argv, **_kwargs):
        nonlocal cyipopt_installed
        call = [str(value) for value in argv]
        calls.append(call)
        if "cyipopt" in call:
            cyipopt_installed = True
        return types.SimpleNamespace(stdout="GPU 0", stderr="", returncode=0)

    def fake_find_spec(name):
        if name == "cyipopt" and not cyipopt_installed:
            return None
        return object()

    real_isdir = os.path.isdir

    def fake_isdir(path):
        nonlocal content_checks
        if str(path) == "/content":
            content_checks += 1
            return content_checks > 1
        return real_isdir(path)

    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(importlib.util, "find_spec", fake_find_spec)
    monkeypatch.setattr(
        importlib.metadata,
        "version",
        lambda name: "0.3.3" if name == "mlmm-toolkit" else "test",
    )
    monkeypatch.setattr(os.path, "isdir", fake_isdir)

    exec(compile(setup, str(NOTEBOOK), "exec"), {})

    joined = [" ".join(call) for call in calls]
    apt_update = next(i for i, call in enumerate(joined) if "apt-get update" in call)
    apt_ipopt = next(
        i for i, call in enumerate(joined) if "coinor-libipopt-dev" in call
    )
    pip_cyipopt = next(i for i, call in enumerate(joined) if "pip install -q cyipopt" in call)
    pip_pydmf = next(i for i, call in enumerate(joined) if "pydmf>=1.2" in call)
    assert apt_update < apt_ipopt < pip_cyipopt < pip_pydmf
    assert cyipopt_installed


def test_colab_setup_dft_branch_installs_extra_and_checks_gpu(monkeypatch, capsys) -> None:
    import importlib.metadata

    # `install_dft` already defaults to True; the backend is pinned explicitly so
    # this contract stays about the DFT extra rather than the current default.
    setup = _notebook()["cells"][1]["source"].replace(
        'backend = "uma"', 'backend = "mace"', 1,
    )
    calls: list[list[str]] = []

    def fake_run(argv, **_kwargs):
        calls.append([str(value) for value in argv])
        return types.SimpleNamespace(stdout="GPU 0", stderr="", returncode=0)

    popen_calls: list[list[str]] = []

    class _FakePopen:
        """Records the argv of every subprocess the Setup cell launches."""

        def __init__(self, argv, **_kwargs):
            popen_calls.append([str(value) for value in argv])
            self.stdout = iter([
                "Collecting pyscf\n",
                "Downloading pyscf-2.13.0-cp311.whl (30.2 MB)\n",
                "a quiet line that must stay unechoed\n",
                "Installing collected packages: pyscf\n",
                "Successfully installed pyscf-2.13.0\n",
            ])

        def wait(self):
            return 0

    real_import_module = importlib.import_module
    fake_cupy = types.SimpleNamespace(
        cuda=types.SimpleNamespace(
            runtime=types.SimpleNamespace(getDeviceCount=lambda: 1),
        ),
    )

    def fake_import_module(name):
        if name == "cupy":
            return fake_cupy
        if name in {"pyscf", "basis_set_exchange", "gpu4pyscf.dft"}:
            return types.SimpleNamespace()
        return real_import_module(name)

    versions = {
        "pyscf": "2.11.0",
        "gpu4pyscf-cuda12x": "1.5.2",
        "mlmm-toolkit": "0.3.3",
    }
    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(subprocess, "Popen", _FakePopen)
    monkeypatch.setattr(importlib.util, "find_spec", lambda _name: object())
    monkeypatch.setattr(importlib, "import_module", fake_import_module)
    monkeypatch.setattr(importlib.metadata, "version", lambda name: versions[name])
    namespace: dict = {}
    exec(compile(setup, str(NOTEBOOK), "exec"), namespace)

    installs = [argv for argv in calls if "install" in argv]
    # One pinned install carries the extra, so the DFT branch differs from the
    # plain branch only by the `[dft]` marker on the same requested version.
    assert any("mlmm-toolkit[dft]==0.3.3" in argv for argv in installs)
    assert not any("mlmm-toolkit==0.3.3" in argv for argv in installs)
    assert popen_calls == []          # no streamed pip log, only the announcement
    logged = capsys.readouterr().out
    assert "install_dft is ticked" in logged
    assert "DFT support installed: PySCF 2.11.0" in logged
    assert "Collecting" not in logged
    assert "Downloading MACE-OMOL-0 model weights (fp64)" in setup
    assert namespace["DFT_SETUP_READY"] is True
    assert namespace["INSTALL_DFT"] is True
    assert namespace["BACKEND"] == "mace"
    assert "importlib.util.find_spec(module)" in setup
    assert "_dft_imports = ('pyscf', 'basis_set_exchange', 'gpu4pyscf.dft')" in setup
    assert "for _module in _dft_imports: importlib.import_module(_module)" in setup
    assert "_cupy.cuda.runtime.getDeviceCount()" in setup
    assert "DFT packages installed but failed their import/GPU check" in setup
    assert "DFT support installed: PySCF %s · GPU4PySCF %s" in setup
    assert "py3Dmol" not in setup
    assert "pip('ipywidgets==7.7.2','anywidget==0.11.0','matplotlib')" in setup
    assert "[%%d/6] %%s".replace("%%", "%") in setup
    assert "time.monotonic()" in setup
    assert ".pysisyphusrc" not in setup
    assert "nglview" not in setup
    assert "first run takes several minutes" in setup


@pytest.mark.parametrize(
    ("backend", "use_token", "expected_login"),
    [
        ("orb", False, None),
        ("uma", True, "token"),
        ("uma", False, "notebook"),
    ],
)
def test_colab_setup_operates_orb_and_uma_branches(
    monkeypatch, backend, use_token, expected_login,
) -> None:
    import importlib.metadata

    setup = _notebook()["cells"][1]["source"].replace(
        'backend = "uma"', f'backend = "{backend}"', 1,
    ).replace("install_dft = True", "install_dft = False", 1)
    calls: list[list[str]] = []
    logins: list[tuple] = []

    def fake_run(argv, **_kwargs):
        calls.append([str(value) for value in argv])
        return types.SimpleNamespace(stdout="GPU 0", stderr="", returncode=0)

    fake_hf = types.ModuleType("huggingface_hub")
    fake_hf.login = lambda **kwargs: logins.append(("token", kwargs))
    fake_hf.notebook_login = lambda: logins.append(("notebook", {}))
    monkeypatch.setitem(sys.modules, "huggingface_hub", fake_hf)
    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(importlib.util, "find_spec", lambda _name: object())
    monkeypatch.setattr(
        importlib.metadata,
        "version",
        lambda name: "0.3.3" if name == "mlmm-toolkit" else "test",
    )
    monkeypatch.delenv("HF_TOKEN", raising=False)
    if use_token:
        monkeypatch.setenv("HF_TOKEN", "test-token")

    namespace: dict = {}
    exec(compile(setup, str(NOTEBOOK), "exec"), namespace)

    joined = [" ".join(argv) for argv in calls]
    assert namespace["BACKEND"] == backend
    assert not any("mace-torch" in command for command in joined)
    assert not any("uninstall -y -q fairchem-core" in command for command in joined)
    if backend == "orb":
        assert any("orb-models" in command for command in joined)
        assert logins == []
    else:
        assert not any("orb-models" in command for command in joined)
        assert [entry[0] for entry in logins] == [expected_login]
        if use_token:
            assert logins[0][1]["token"] == "test-token"


def test_colab_setup_handles_cancelled_uma_sign_in(monkeypatch, capsys) -> None:
    import importlib.metadata

    setup = _notebook()["cells"][1]["source"].replace(
        "install_dft = True", "install_dft = False", 1,
    )

    def fake_run(_argv, **_kwargs):
        return types.SimpleNamespace(stdout="GPU 0", stderr="", returncode=0)

    def cancelled_login():
        raise KeyboardInterrupt

    fake_hf = types.ModuleType("huggingface_hub")
    fake_hf.login = lambda **_kwargs: None
    fake_hf.notebook_login = cancelled_login
    monkeypatch.setitem(sys.modules, "huggingface_hub", fake_hf)
    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(importlib.util, "find_spec", lambda _name: object())
    monkeypatch.setattr(
        importlib.metadata,
        "version",
        lambda name: "0.3.3" if name == "mlmm-toolkit" else "test",
    )
    monkeypatch.delenv("HF_TOKEN", raising=False)

    namespace: dict = {}
    exec(compile(setup, str(NOTEBOOK), "exec"), namespace)

    assert namespace["_hf_auth_ready"] is False
    output = capsys.readouterr().out
    assert "Hugging Face sign-in was cancelled" in output
    assert "Installation is complete" in output
    assert "authenticate UMA before running a calculation" in output


def test_colab_setup_explains_unavailable_release(monkeypatch) -> None:
    import importlib.metadata

    setup = _notebook()["cells"][1]["source"].replace(
        "install_dft = True", "install_dft = False", 1,
    )

    def fake_run(argv, **_kwargs):
        command = [str(value) for value in argv]
        failed = any(value == "mlmm-toolkit==0.3.3" for value in command)
        return types.SimpleNamespace(stdout="GPU 0", stderr="", returncode=1 if failed else 0)

    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(importlib.util, "find_spec", lambda _name: object())
    monkeypatch.setattr(importlib.metadata, "version", lambda _name: "0.3.3")

    with pytest.raises(RuntimeError) as error:
        exec(compile(setup, str(NOTEBOOK), "exec"), {})

    message = str(error.value)
    assert "Could not install mlmm-toolkit==0.3.3 from PyPI" in message
    assert "version v0.3.3 may not be published" in message
    assert "mlmm-toolkit-src.zip pair, then enter debug" in message


def test_colab_setup_explains_missing_debug_zip(monkeypatch, tmp_path) -> None:
    setup = _notebook()["cells"][1]["source"].replace(
        'mlmm_toolkit_version = "v0.3.3"', 'mlmm_toolkit_version = "debug"', 1,
    ).replace("install_dft = True", "install_dft = False", 1)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda _argv, **_kwargs: types.SimpleNamespace(
            stdout="GPU 0", stderr="", returncode=0,
        ),
    )

    with pytest.raises(FileNotFoundError) as error:
        exec(compile(setup, str(NOTEBOOK), "exec"), {})

    message = str(error.value)
    assert "Debug mode needs the adjacent source ZIP" in message
    assert "mlmm-toolkit-src.zip" in message
    assert "enter v0.3.3 for a published release install" in message


def test_colab_gui_is_mlmm_native_and_tracks_structure_contracts() -> None:
    app = _notebook()["cells"][2]["source"]
    whole = NOTEBOOK.read_text(encoding="utf-8").lower()

    assert isinstance(app, str)
    assert "other-tool" not in whole
    assert ".pdb,.ent,.cif,.mmcif,.parm7" in app
    assert "'.pdb,.ent,.cif,.mmcif,.parm7,.xyz,.gjf,.com,.inp,.csv'" in app
    assert "prepare_input_structure" in app
    assert "load_pdb_atom_metadata" in app
    assert "Path(_runtime_path('viewer_input.pdb'))" in app
    assert "callback_ns = 'mlmm_gui'" in app
    assert "_co.register_callback('mlmm_gui.on_click', on_click)" in app
    assert "for(let depth=0;depth<8&&scope;depth++)" in app
    assert "const safeArgs=bridge.scope.JSON.parse(JSON.stringify(args))" in app
    assert "bridge.api.invokeFunction(cfg.callback+'.'+suffix,safeArgs,kwargs)" in app
    assert "window.parent.postMessage({type:'rx-kernel-invoke',suffix:suffix,args:args},'*')" in app
    assert "if(window.parent&&window.parent!==window)" in app
    assert "const seenClickEvents=new WeakSet()" in app
    assert "function queueInteractiveClick(event)" in app
    assert "seenClickEvents.has(event)" in app
    assert "seenClickEvents.add(event)" in app
    assert "viewer.plugin.behaviors.interaction.click.subscribe(queueInteractiveClick)" in app
    assert "viewer.plugin.canvas3d.interaction.click.subscribe(queueInteractiveClick)" in app
    assert "function installIframeKernelRelay(scope)" in app
    assert "function wireIframeKernelRelay()" in app
    assert "scope.document.querySelectorAll('iframe.rxmolstar-frame')" in app
    assert "for(var depth=0;depth<6&&scope;depth++)" in app
    assert "frame.contentWindow===event.source" in app
    assert "['on_click','clear_highlights','set_frame'].indexOf(suffix)<0" in app
    assert "'viewer_callback_base': 'mlmm_gui'" in app
    assert "CONFIG.viewer_callback_base+'.'+suffix" in app
    assert "stateNode.style.display='flex'" in app
    assert "else if(loci.kind==='bond-loci')" in app
    assert "SE.Location.create(bond.aStructure,bond.aUnit,bond.aUnit.elements[bond.aIndex])" in app
    assert "CHAIN:RESNAME:RESSEQ" in app
    assert "String(item.sourceIndex)" in app
    assert "ML/MM compute commands require a matching Amber parm7" not in app
    assert "use mlmm all to generate it automatically" in app
    assert "existing parm7 optional" in app
    assert "Structures (.pdb/.cif/.xyz/.gjf)" in app
    assert "elif parm:" in app
    assert "S['parm'] = parm" in app


def test_colab_gui_uses_public_cli_and_backend_neutral_results() -> None:
    app = _notebook()["cells"][2]["source"]

    assert "cmd += ['--backend-model', mdl]" in app
    assert "model_config.yaml" not in app
    assert "cmd.append('--tsopt')" in app
    assert "cmd.append('--thermo')" in app
    assert "[n, 'True']" not in app
    assert "'MLIP_Gibbs': 'gibbs_mlip', 'MLIP': 'mlip'" in app
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
    assert "scope = _preflight_output_scope(a)" in app
    assert "target = scope['target']; out = scope['root']" in app
    assert "out = scope['root']" in app
    assert "os.path.isdir(out)" in app
    assert "_results(out)" in app
    assert "fallback = _next_numbered_output_dir(previous)" in app
    assert "Run not started: refusing to overwrite the existing output:" in app


def test_colab_gui_keeps_responsive_release_layout() -> None:
    app = _notebook()["cells"][2]["source"]

    assert "@media (max-width: 600px)" in app
    assert "max_width='100%'" in app
    assert "flex:0 0 auto; min-width:0" in app
    assert "layout=W.Layout(width='260px', max_width='100%')" in app
    assert "_MOLSTAR_VERSION = '5.11.0'" in app
    assert "molstar@%s/build/viewer/molstar.js" in app
    assert "layoutShowSequence:cfg.showSequence" in app
    assert "layoutShowControls:true" in app
    assert "viewportShowControls:true" in app
    assert "viewportShowSelectionMode:true" in app
    assert "collapseRightPanel:true" in app
    assert "viewportFocusBehavior" not in app
    assert "options=[('Cartoon', 'cartoon'), ('Stick', 'stick')]" in app
    assert "description='Representation'" in app
    assert "view_controls = W.HBox([" in app
    assert "    dd_rep," in app
    assert "cb_water," in app
    assert "view_controls.add_class('rxview-controls')" in app
    assert ".rxview-controls { flex-flow:row nowrap !important;" in app
    assert ".rxview-controls > .rxinfo { flex:0 0 26px !important; }" in app
    assert "last_pick_row = W.HBox([last_pick_html]," in app
    assert "btn_clear_pick" not in app
    assert "sel_lang" not in app
    assert "No structure loaded" in app
    assert "options=[('MEP', 'mep'), ('Scan', 'scan'), ('TS-only', 'tsonly')]" in app
    assert "style={'button_width': '74px', 'description_width': '0px'}" in app
    assert "No results yet." in app
    assert "frame_slider.disabled = last == 0" in app
    assert "trajectory_box.layout.display = 'none'" in app
    # Colab renders ipywidgets' Tab and Accordion as empty blocks. ToggleButtons
    # updates immediately while all panes remain mounted for uploads and WebGL.
    assert "_TAB_PAGES = [('① Input', input_box), ('② Setup', viewer_box)," in app
    assert "('③ Options', options_box), ('④ Results', results_box)]" in app
    assert "_tab_body = W.VBox([page for _label, page in _TAB_PAGES])" in app
    assert "_tab_body.children = [_TAB_PAGES[i][1]]" not in app
    assert "_pane.add_class('rxpage-prewarm')" in app
    assert "_pane.remove_class('rxpage-prewarm')" in app
    assert "elif _j == 1:" in app
    assert "Options (optional)" not in app
    assert "def _tab_go(i):" in app
    assert "def _finish_tab_go(i, request):" in app
    assert "def _defer_results_for_tab(request):" in app
    assert "timer = threading.Timer(0.16" in app
    assert "request != _TAB_NAV['request'] or _TAB_NAV['active'] != 3" in app
    assert "_pane.layout.display = 'block'" in app
    assert "_pane.layout.display = ('flex' if _j in (i, 1) else 'none')" in app
    assert "_normalise_result_label(options[selected_index][0])" in app
    assert "def _defer_tab_finish" not in app
    assert "_finish_tab_go(i, _request)" in app
    assert "_TAB_NAV = {'active': 0, 'syncing': False, 'request': 0}" in app
    assert "_tab_choice = W.ToggleButtons(" not in app
    assert "_tab_strip = W.HBox(_tab_buttons" in app
    assert "_button.add_class('rxtab-button')" in app
    assert "if not IS_COLAB_FRONTEND:" in app
    assert "_button.on_click(lambda _clicked, _index=_index: _tab_go(_index))" in app
    assert "Loading <b>%s</b>…" in app
    assert "function wireTabs()" in app
    assert "setNativeTabLoading(index,true)" in app
    assert "function setNativeTabPane(index)" in app
    assert "setNativeTabChoice(index); setNativeTabPane(index); setNativeTabLoading(index,true);" in app
    assert "function settle(ok)" in app
    assert "setNativeTabChoice(previous)" in app
    assert "Number(result.active)===index" in app
    assert "bridge.invokeFunction(CONFIG.tab_callback,[index],{})" in app
    assert "160-(performance.now()-shownAt)" in app
    assert "'tab_callback': 'mlmm_gui.tab_go'" in app
    assert "register_callback('mlmm_gui.tab_go', _on_colab_tab)" in app
    assert "register_callback('mlmm_gui.select_result', _on_colab_result)" in app
    assert "function wireResultChoices()" in app
    assert "document.__rxResultChoicesWired" in app
    assert "document.addEventListener('change',function(event)" in app
    assert "select.addEventListener('change'" not in app
    assert "document.querySelector(selector)||select" in app
    assert "rxResultBusy" in app
    assert "setNativeResultLoading(label,true)" in app
    assert "var selectedIndex=select.selectedIndex" in app
    assert "bridge.invokeFunction(CONFIG.result_callback,[kind,selectedIndex,generation,label],{})" in app
    assert "generation != _RESULT_SET_GENERATION['value']" in app
    assert "value = options[selected_index][1]" in app
    assert "var nativeFrameTimer=0" in app
    assert "playButtons[0].dataset.rxNativePlay='true'" in app
    assert "playButtons[1].textContent='Ⅱ'" in app
    assert "playButtons[1].title='Pause'" in app
    assert "var repeatButton=playButtons.length>=3?playButtons[playButtons.length-1]:null" in app
    assert "if(!repeatButton||!repeatButton.classList.contains('mod-active'))" in app
    assert "if(!vibration){stopNativePlayback(true);return;}" not in app
    assert "function clearNativeFrameTimer()" in app
    assert app.count("clearNativeFrameTimer();") >= 3
    assert "frame_play.add_class('rxframe-play')" in app
    result_selector_css = re.search(r"\.rxresult-selector\s*\{([^{}]*)\}", app)
    assert result_selector_css is not None
    assert "display:grid!important" not in re.sub(r"\s+", "", result_selector_css.group(1))
    assert "var warningHost = document.querySelector('.rxdrop');" in app
    assert "warningHost.parentElement.insertBefore(warning, warningHost.nextSibling);" in app
    assert "File upload controls did not initialize. Rerun the Launch GUI cell." in app
    assert "Upload controls could not attach." not in app
    assert "try: _tab_status.send_state()" in app
    assert "W.Tab(" not in app
    assert "def _collapsible(title, child, on_open=None):" in app
    assert "W.Accordion(" not in app
    assert "example_fold = _collapsible('Examples'" in app
    assert "_collapsible('Manual ML-region PDB settings', model_upload_box)" in app
    assert "Try an example" not in app
    assert "Use existing ML-region PDB" not in app
    assert "--refine-path (MEP refinement & multi-step reaction detection)" in app
    assert ".widget-checkbox .widget-label-basic { white-space:normal !important;" in app
    assert ".rxflagrow:not(:has(.rxinfo-details[open])) { flex-wrap:nowrap !important; }" in app
    assert ".rxflagrow:has(.rxinfo-details[open]) { flex-wrap:wrap !important; }" in app
    assert ".rxflagrow:not(:has(.rxinfo-details[open])) > .rxinfo" in app
    assert "height:auto !important; align-items:flex-start !important;" in app
    assert "text-overflow:clip !important;" in app
    assert "_input_box_children.append(example_fold)" in app
    # Browser-native details/summary opens without a Python callback.
    assert "def _info_markup(tip, revision=0, rich=False):" in app
    assert "def _info_control(tip, target=None, rich=False):" in app
    assert "def _set_info_text(control, tip):" in app
    assert "def _close_info_target(target):" in app
    assert "def _hdr(content, tip):" in app
    assert "def _flag_row(widget, tip, info_target=None, rich=False):" in app
    assert "control._rx_rich = bool(rich)" in app
    assert "getattr(item, '_rx_rich', False)" in app
    assert ".rxhelp-table-scroll { max-width:100%; overflow-x:auto; }" in app
    assert "'<div class=\"rxhelp-table-scroll\"><table>" in app
    assert "_flag_row(adv_thresh, _THRESH_HELP, rich=True)" in app
    assert "_flag_row(adv_thresh_post, _THRESH_HELP, rich=True)" in app
    assert "if param.name in {'thresh', 'thresh_post'}:" in app
    assert "row = _flag_row(widget, _THRESH_HELP, rich=True)" in app
    assert app.count("rich=True)") >= 7
    assert "'<summary aria-label=\"More information: %s\" title=\"%s\">&#9432;</summary>'" in app
    assert "'<details class=\"rxinfo-details\" data-revision=\"%d\">'" in app
    assert "'<div class=\"rxhelp-panel\" role=\"note\"><small>%s</small>" in app
    assert "_INFO_CONTROLS = weakref.WeakSet()" in app
    assert "description='Show information', icon='info-circle'" not in app
    assert "_rx_info_button" not in app
    assert "rxinfo-popover" not in app
    assert "📥 needs" not in app
    assert "w_show_run_log = W.Checkbox(value=True, description='Show log'" in app
    assert "cmdline_box = W.VBox([command_editor, command_footer, logbox])" in app
    assert "command_editor = W.VBox([" in app
    assert "W.HTML('<b>Command line</b>')" in app
    assert "command_editor = _collapsible('Command line'" not in app
    assert "viewer_more = W.VBox([" in app
    assert "viewer_more.layout.display = 'none'" in app
    assert "_collapsible('Exact atom input'" not in app
    assert "viewer_toolbar = W.HBox(" in app
    assert "[view_input, pick_action, last_pick_info, view_controls]" in app
    assert "[view_input, pick_action, last_pick_info, view_controls, viewer_more]" not in app
    assert "for _control in (pick_action, last_pick_info, viewer_more):" not in app
    assert "viewer_box = W.VBox([workflow_box, viewer_toolbar, view_input_note, selection_box])" in app
    assert "'<div class=\"rxmolstar-embed\">%s</div>' % _document_iframe(" in app
    assert "workflow_contract_row = W.HBox([subreq, outputs_html]" in app
    assert "workflow_controls, workflow_contract_row, subcmd_note])" not in app
    assert "workflow_controls, subcmd_note])" in app
    assert "workflow_controls.add_class('rxworkflow-controls')" in app
    assert "grid-template-columns:minmax(230px,270px) minmax(140px,1fr) 246px" in app
    assert "grid-template-columns:minmax(0,2fr) minmax(260px,1fr);" in app
    assert _css_rule_has(
        app,
        ".rxapp .rxworkspace",
        "overflow-x:clip",
        "overflow-y:visible",
    )
    assert not _css_rule_has(app, ".rxapp .rxworkspace", "max-height:clamp(")
    assert not _css_rule_has(app, ".rxapp .rxworkspace", "scrollbar-gutter:stable")
    assert ".rxapp .rxworkspace::-webkit-scrollbar" not in app
    assert _css_rule_has(app, ".rxviewer", "position:static")
    assert ".rxviewer-page { min-height:0 !important; row-gap:6px !important; align-items:stretch !important; }" in app
    assert "padding:0 0 2px; box-sizing:border-box;" in app
    assert ".rxviewer-page > * { width:100% !important; max-width:100% !important;" in app
    assert ".rxworkspace > * { width:100% !important; max-width:100% !important;" in app
    assert ".rxviewer-toolbar > .widget-dropdown { flex:1 1 15rem !important;" in app
    assert "layout=W.Layout(width='64px', flex='0 0 64px')" in app
    assert ".rxviewer > .rxpick-footer { margin-top:clamp(8px,.7vw,12px) !important; }" in app
    assert "@media (min-width: 821px) and (max-height: 900px)" in app
    assert ".rxapp-main { flex:0 0 auto !important; min-height:0; overflow:visible; }" in app
    assert ".rxpages { flex:0 0 auto !important; min-height:0; overflow:visible; position:relative; }" in app
    assert ".rxpage-prewarm {" in app
    assert "left:-200vw !important" in app
    assert ".rxtabs .rxtab-button { flex:1 1 0 !important; width:auto !important;" in app
    assert ".rxpage > * { flex:0 0 auto !important; }" not in app
    assert "height:auto; max-height:none; min-height:0; overflow:visible;" in app
    assert "height:calc(clamp(460px,50vw,660px) + 52px); max-height:712px;" in app
    assert "height:clamp(460px,50vw,660px) !important; max-height:660px !important;" in app
    assert "overflow-y:auto !important; overflow-x:hidden !important;" in app
    assert "scrollbar-gutter:auto; scrollbar-width:auto; scrollbar-color:#2563eb #dbeafe;" in app
    assert ".rxinspector::-webkit-scrollbar { width:14px; }" in app
    assert ".rxinspector::-webkit-scrollbar-track { background:#dbeafe; border:1px solid #bfdbfe;" in app
    assert ".rxinspector::-webkit-scrollbar-thumb { min-height:52px; background:#2563eb;" in app
    assert ".rxinspector::-webkit-scrollbar-thumb:hover { background:#1d4ed8; }" in app
    assert "scrollbar-color:#9CA3AF #EEF1F4" not in app
    assert ".rxscan-pair-controls, .rxscan-add-actions {" in app
    assert "overflow-x:hidden !important; overflow-y:visible !important;" in app
    assert "border:0; border-radius:0; background:transparent; outline:none; box-shadow:none;" in app
    assert "padding-right:7px; box-sizing:border-box;" in app
    assert "align-items:stretch !important;" in app
    assert ".rxinspector > * {" in app
    assert "margin-left:0 !important; margin-right:0 !important; box-sizing:border-box !important;" in app
    assert ".rxinspector > .widget-html > .widget-html-content," in app
    assert ".rxviewer .rxmolstar-embed { width:100%; margin-inline:0; border:0;" in app
    assert "background:#fff; outline:1px solid #E0DDD4; outline-offset:-1px;" in app
    assert "box-shadow:inset 0 0 8px rgba(51,43,31,.055)" in app
    assert ".rxviewer > .widget-html, .rxviewer > .rxpick-footer {" in app
    assert ".rxpick-footer .widget-html-content, .rxpick-footer [role=\"status\"] {" in app
    assert "last_pick_row.add_class('rxpick-footer')" in app
    assert "aspect-ratio:auto !important; border:0 !important;" in app
    assert "box-shadow:0 0 5px rgba(16,24,40,0.065)" in app
    assert "border-right:2px solid #bfdbfe" not in app
    assert re.search(r"box-shadow\s*:\s*0\s+(?:[1-9]\d*|-\d+)px", app) is None
    assert "width:100% !important; min-width:0; max-width:none;" in app
    assert ".rxinspector > .rxcard { width:100%; min-width:0; max-width:100%;" in app
    mobile_css = _css_media_body(app, "max-width: 940px")
    assert _css_rule_has(mobile_css, ".rxapp .rxworkspace", "display:flex")
    assert _css_rule_has(mobile_css, ".rxinspector", "overflow-y:visible")
    assert _css_rule_has(mobile_css, ".rxinspector", "height:auto")
    assert _css_rule_has(mobile_css, ".rxinspector", "border:0")
    assert _css_rule_has(mobile_css, ".rxinspector", "box-shadow:none")
    assert "setup_layout" not in app
    assert "setup_folds" not in app
    assert "selection_inspector = W.VBox([" in app
    assert "center_panel, charge_panel, system_charge_panel, scan_panel, freeze_panel," in app
    assert "extract_panel])" in app
    assert "rxradius-step-buttons" not in app
    assert "rxradius-control" not in app
    assert "appearance:textfield" not in app
    assert ".rxcenter-panel { overflow-x:clip !important; overflow-y:visible !important; }" in app
    assert ".rxcenter-panel > *, .rxcenter-panel .widget-hbox" in app
    assert ".rxcenter-panel::-webkit-scrollbar, .rxcenter-panel > *::-webkit-scrollbar" in app
    assert "background:#fff !important" in app
    assert "background:linear-gradient(145deg,#f8fbff 0%,#edf5ff 100%)" not in app
    assert "0 0 0 1px rgba(37,99,235,.12),0 0 16px rgba(30,64,175,.09)" in app
    assert ".rxfold-toggle::before { content:'▶'" in app
    assert ".rxfold-toggle.rxfold-open::before { content:'▼'" in app
    assert "_btn.description = ('Hide ' if _st['open'] else 'Show ') + title" in app
    assert "_btn.description = ('▾ ' if _st['open'] else '▸ ')" not in app
    assert "inset 4px 0 0 var(--rx-blue)" not in app
    assert ".rxtabs .widget-toggle-button" not in app
    assert ".rxpages-loading { opacity:.18 !important; pointer-events:none !important; }" in app
    assert ".rxpath-panel svg, .rxpath-panel img, .rxpath-panel canvas {" in app
    assert "traj_out = W.HTML(layout={'width': '100%', 'min_width': '0'})" in app
    assert "plot_out = W.HTML(layout={'width': '100%', 'min_width': '0'})" in app
    assert ".rxresults-actions { flex:0 0 auto !important; min-height:38px; }" not in app
    assert "results_actions.add_class('rxresults-actions')" in app
    assert "'flex': '1 1 440px'" not in app
    assert ".rxcommand-dock {" in app
    assert "width:100% !important; max-width:100% !important; min-width:0;" in app
    assert ".rxapp .rxcommand-editor {" in app
    assert "padding:8px 10px !important; row-gap:8px !important; align-items:stretch !important;" in app
    assert ".rxcommand-editor > * { width:100% !important; max-width:100% !important;" in app
    assert ".rxapp .rxcommand-editor { padding:6px 8px !important; row-gap:6px !important; }" in app
    assert "layout=W.Layout(width='99%'" not in app
    assert "border-top:1px solid var(--rx-line); padding-top:13px" in app
    assert "b_validate = W.Button(description='Validate', button_style='info'" in app
    assert "b_validate.add_class('rxvalidate')" not in app
    assert "b_run = W.Button(description='Run', button_style='danger'" in app
    assert "b_run.add_class('rxexecute')" not in app
    assert (
        "    if real_run: _open_run_log()\n"
        "    _RUN_EXECUTION['argv'] = list(a)\n"
        "    _set_running(True)\n"
        "    run_log_path = None\n"
        "    try:"
    ) in app
    assert "rootbox = W.VBox([header, manual_mode_notice, plotly_preload_out, app, cmdline_box])" in app
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
    assert "class _DropUpload(anywidget.AnyWidget):" in app
    assert "upl = _DropUpload(accept=_acc, formats=_drop_formats, formats_detail=_drop_formats_detail)" in app
    assert "upl = W.FileUpload(accept=_acc, multiple=True, description='Upload files'" in app
    assert "_drop = W.VBox(_drop_children," in app
    assert "align_items='center', justify_content='center'" in app
    assert "el.addEventListener('drop', onDrop)" in app
    assert "file.arrayBuffer()" in app
    assert "model.send({" in app
    assert "_DROP_STATE = {'generation': 0, 'seen': set()}" in app
    assert "def _claim_drop_batch(batch, generation):" in app
    assert "def _bump_drop_generation():" in app
    assert "queue.push({id: opaqueId(), files: selected, generation})" in app
    assert "if (busy || !queue.length || signal.aborted) return;" in app
    assert "message.batch !== active" in app
    assert "if clear:\n        _bump_drop_generation()" in app
    assert "mlmm_gui.on_drop" not in app
    # Hosted Colab uses the native bridge; a Colab local runtime needs the
    # standard FileUpload fallback because google.colab is unavailable there.
    assert "_UPLOAD_MODE = ('colab' if IN_COLAB else" in app
    assert "_colab_frontend = bool(get_ipython().parent_header.get('metadata', {}).get('colab'))" in app
    assert "IS_COLAB_FRONTEND = IN_COLAB or _colab_frontend" in app
    assert "'basic' if IS_COLAB_FRONTEND else" in app
    assert "if _HAS_DROP_WIDGET and (IN_COLAB or not IS_COLAB_FRONTEND):" in app
    assert "if _HAS_DROP_WIDGET and not IS_COLAB_FRONTEND:" not in app
    assert "host.scrollTop = host.scrollHeight" in app
    assert 'host.dataset.rxFollowLog = "smart"' in app
    assert 'style.textContent = `' in app
    assert 'overflow-y:scroll' in app
    assert 'el.replaceChildren(style, shell)' in app
    assert 'rail.className = "rxlog-scrollbar"' in app
    assert 'thumb.className = "rxlog-scroll-thumb"' in app
    assert 'host.addEventListener("scroll", updateThumb' in app
    assert 'thumb.addEventListener("pointermove"' in app
    assert "host.scrollHeight - host.scrollTop - host.clientHeight" in app
    assert "settleScroll(shouldFollow)" in app
    assert "if(!pinned)return" in app
    assert "artifact_fold = _collapsible('Generated file preview', artifact_box, on_open=_render_artifact)" in app
    assert "artifact_fold._rx_set_open(bool(visuals and not result_views and not energy_options))" in app
    assert "artifact_fold._rx_body.layout.display == 'none'" in app
    assert "tick0:0,dtick:cfg.xTickStep" in app
    assert "showticklabels:true,ticks:'',ticklen:0" in app
    assert "_drop_children = ([upl] if _UPLOAD_MODE == 'anywidget'" in app
    assert "_cwm.register_callback('mlmm_gui.upload_files', _on_colab_upload)" in app
    assert "_cwm.register_callback('mlmm_gui.load_example', _on_colab_example)" in app
    assert "var uploadPending = 0;" in app
    assert "if(!wire(CONFIG.zones[i])) uploadPending += 1;" in app
    assert "if not IS_CLUSTER and 'model_upl' in globals():" in app
    assert "if(!wireTabs())pending" not in app
    assert "if(!wireCancel())pending" not in app
    assert "if(warning)warning.remove();" in app
    assert "def _set_operation_loading(label, active):" in app
    assert "Loading <b>%s</b>…</div>" in app
    assert "ex_btn.add_class('rxoperation-trigger')" in app
    assert "def _load_example_impl(_):" in app
    assert "try: return _load_example_impl(_)" in app
    assert "_set_operation_loading('files', True)" in app
    assert "_set_operation_loading('session', True)" in app
    assert "_set_operation_loading('ML-region PDB', True)" in app
    assert "busy_label = {'model': 'ML-region PDB', 'session': 'session'}.get(role, 'files')" in app
    assert "function setNativeOperationLoading(label,active)" in app
    assert "'example_callback': 'mlmm_gui.load_example'" in app
    assert "bridge.invokeFunction(CONFIG.example_callback,[],{})" in app
    assert "setNativeOperationLoading('example',false)" in app
    assert "host.dataset.rxOperationBusy==='true'" in app
    assert "wireOperationTriggers();" in app
    assert "Atom identifiers match input 1" not in app
    # Incremental logs use AnyWidget in hosted Colab and ordinary Jupyter,
    # while Colab local runtimes retain the standard widget path.
    assert "if _HAS_DROP_WIDGET and (IN_COLAB or not IS_COLAB_FRONTEND):" in app
    # callback_ns is local to _molstar_document; the module-level block needs a literal.
    assert "callback_ns" not in app.split("display(rootbox)")[1]
    assert "bridge.invokeFunction(CONFIG.callback, [spec.role, payload], {})" in app
    assert "reader.readAsDataURL(file)" in app
    assert "def _decode_colab_batch(files):" in app
    assert "base64.b64decode(entry['b64'], validate=True)" in app
    for _role in ("'input'", "'model'", "'session'"):
        assert "role=%s" % _role in app or "role == %s" % _role in app
    assert "_accept_upload_pairs(pairs, 'browser upload')" in app
    assert "_accept_upload_pairs(pairs, 'browser upload')" in app
    assert "def _delete_owned_uploads(paths):" in app
    assert "_HAS_DROP_WIDGET" in app
    assert "if _UPLOAD_MODE == 'anywidget': upl.on_msg(_on_drop_upload)" in app
    assert "description='Move earlier', icon='arrow-up'" in app
    assert "description='Move later', icon='arrow-down'" in app
    assert "description='Remove file', icon='times'" in app
    assert "tooltip='Move earlier'" in app and "tooltip='Move later'" in app
    assert "command_footer = W.VBox([command_actions, w_show_run_log]" in app
    assert "cmdline_box.add_class('rxcommand-dock')" in app
    assert "for _label, _page in _TAB_PAGES:" in app
    assert "_page.add_class('rxpage')" in app
    assert "_tab_body.add_class('rxpages')" in app
    assert "app.add_class('rxapp-main')" in app
    assert "def _advanced_coverage(" in app and "adv_extra" not in app
    assert "def _advanced_options(sub):" in app
    assert "adv_acc = _collapsible('All flags', adv_box)" in app
    assert "row._rx_tier = 'advanced' if bool(getattr(param, 'hidden', False)) else 'key'" in app
    assert "every CLI option accounted for" not in app
    assert "radius_applies = (sub == 'extract' or" in app
    assert "adv_radius.disabled = not radius_applies" in app
    assert "_set_flag_visible(adv_radius, radius_applies)" in app
    assert "_set_flag_visible(adv_dftfb, _dftfb_applicable)" in app


def test_colab_viewer_persists_exact_atom_and_residue_context() -> None:
    app = _notebook()["cells"][2]["source"]

    for marker in (
        "'_last_pick': None", "'_pick_history': []", "def _remember_pick",
        "_VIEWER_GENERATION", "_MOLSTAR_VERSION = '5.11.0'",
        "layoutShowSequence:cfg.showSequence", "layoutShowControls:true",
        "collapseRightPanel:true", "layoutShowLog:false",
        "layoutShowLeftPanel:false", "layoutShowRemoteState:false",
        "await molstar.Viewer.create", "params.values.type.name!=='cartoon'",
        "alpha:0.4", "structure-component-static-water",
        "await configureWater()", "tryCreateComponentStatic(",
        "'mmcif' if fmt in ('cif', 'mmcif')",
        "ignoreStartupEmpty", "loci.kind==='empty-loci'",
        "clickQueue=clickQueue.then(()=>handleClick(event))",
        "await invoke('clear_highlights',[cfg.generation])",
        "cfg.generation,exact",
        "exact atom", "set current pick", "view_input", "_view_mapping_ok",
        "last_pick_info", "Generated file preview", "Download current run (.zip)",
        "Results directory", "Load results",
        "results_box.add_class('rxresults')", "overflow-x:auto",
        "run.log", "energy unavailable", "No completed result",
        "Command failed", "frame_slider = W.IntSlider(", "channel='trajectory'",
        "host.on('plotly_click'", "Plotly.restyle", "artifact_fold._rx_set_open",
        "message.type!=='rx-set-frame'", "update.to(model).update",
        "show_sequence=False, channel='trajectory',",
        "generation=generation, frame_count=len(frames)",
        "source = ''.join(frames)",
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
    assert "freeze_acc = freeze_panel" in app
    assert "_PICK_HINT.get(pick_action.value, 'Choose a click action.')" in app
    assert "_PICK_HINT[pick_action.value]" not in app
    assert "pick_action.value = 'scanB'" in app
    assert "pick_action.value = 'freezeB'" in app
    assert "bridge.api.invokeFunction(cfg.callback+'.'+suffix,safeArgs,kwargs)" in app
    assert "String(item.sourceIndex)" in app
    assert "cfg.generation,exact" in app
    assert "if exact in (False, 0, 'false', 'False')" in app
    assert "Mol* focused this residue" in app
    click_completion = app[
        app.index("def on_click("):
        app.index("_co.register_callback('mlmm_gui.on_click', on_click)")
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
    assert contract["_artifact_kind"]("mep_trj.xyz") == "MEP profile"
    assert contract["_artifact_kind"]("final_geometries_trj.xyz") == "MEP profile"
    assert contract["_artifact_kind"]("finished_irc_trj.xyz") == "IRC trajectory"
    assert contract["_artifact_kind"]("mode_0001_-100.00cm-1_trj.xyz") == "vibrational mode"
    assert contract["_artifact_kind"]("imag_-421.25cm-1.pdb") == "frequency structure"
    assert contract["_artifact_kind"]("job.gjf") == "text"
    assert contract["_artifact_kind"]("job.com") == "text"
    assert contract["_artifact_kind"]("job.inp") == "text"
    assert contract["_artifact_kind"]("forcefield.prms") == "text"
    assert contract["_artifact_kind"]("optimization_trj.xyz") is None
    opt = contract["_trajectory_semantics"]("opt", "optimization_trj.xyz")
    path_semantics = contract["_trajectory_semantics"](
        "path-opt", "mep_trj.xyz",
    )
    vibration = contract["_trajectory_semantics"](
        "tsopt", "vib/imag_120i_trj.xyz",
    )
    assert vibration["title"] == "Vibrational-mode trajectory"
    assert vibration["x"] == "phase frame" and not vibration["extrema"]
    assert contract["_stationary"]([0.0, 2.0, 0.0], opt) == [
        (0, "initial"), (2, "optimized"),
    ]
    # A trajectory plot labels R and P only. An
    # energy-only extremum is neither a certified transition state nor a
    # certified intermediate; the energy diagram is where the profile is read.
    assert contract["_stationary"]([0.0, 2.0, 0.0], path_semantics) == [
        (0, "R"), (1, "TS candidate"), (2, "P"),
    ]
    assert contract["_stationary"]([2.0, 0.0, 2.0], path_semantics) == [
        (0, "R"), (2, "P"),
    ]


def test_colab_results_bind_irc_truth_and_skip_bridge_extrema(tmp_path: Path) -> None:
    contract = _viewer_contract()
    irc_dir = tmp_path / "irc"
    irc_dir.mkdir()
    irc_path = irc_dir / "finished_irc_trj.xyz"
    irc_path.write_text("", encoding="utf-8")
    (irc_dir / "result.json").write_text(
        json.dumps({
            "n_frames_forward": 2,
            "n_frames_backward": 2,
            "n_frames_total": 5,
            "forward_converged": True,
            "backward_converged": False,
            "scientific_status": "partial",
            "scientific_status_reasons": ["backward IRC did not converge"],
        }),
        encoding="utf-8",
    )
    irc = contract["_trajectory_semantics"]("irc", str(irc_path), n_frames=5)
    assert irc["title"] == "Combined IRC trajectory"
    assert irc["ts_index"] == 2
    assert "partial" in irc["trajectory_status"]
    assert "forward ✓" in irc["trajectory_status"]
    assert "backward ✗" in irc["trajectory_status"]
    mismatch = contract["_trajectory_semantics"](
        "irc", str(irc_path), n_frames=4,
    )
    assert "ts_index" not in mismatch
    assert mismatch["metadata_warning"] == "IRC frame metadata mismatch"
    forward = contract["_trajectory_semantics"](
        "all", str(irc_dir / "forward_irc_trj.xyz"), n_frames=2,
    )
    backward = contract["_trajectory_semantics"](
        "all", str(irc_dir / "backward_irc_trj.xyz"), n_frames=2,
    )
    assert (forward["title"], forward["start"], forward["end"]) == (
        "Forward IRC branch", "near TS", "forward endpoint",
    )
    assert (backward["title"], backward["start"], backward["end"]) == (
        "Backward IRC branch", "near TS", "backward endpoint",
    )
    assert forward["x"] == "Forward IRC step"
    assert backward["x"] == "Backward IRC step"

    mep_path = tmp_path / "mep_trj.xyz"
    mep_path.write_text("", encoding="utf-8")
    (tmp_path / "summary.json").write_text(
        json.dumps({
            "segments": [
                {"index": 1, "kind": "seg", "frame_ranges": [[0, 2]]},
                {"index": 2, "kind": "bridge", "frame_ranges": [[2, 5]]},
                {"index": 3, "kind": "seg", "frame_ranges": [[5, 7]]},
            ]
        }),
        encoding="utf-8",
    )
    path = contract["_trajectory_semantics"](
        "path-search", str(mep_path), n_frames=7,
    )
    assert path["bridge_ranges"] == [(2, 5)]
    assert contract["_trajectory_frame_context"](1, path)["label"] == "seg 01"
    assert contract["_trajectory_frame_context"](3, path)["label"] == "bridge"
    assert contract["_trajectory_frame_context"](6, path)["label"] == "seg 03"
    irc_context = contract["_trajectory_frame_context"](2, irc)
    assert "backward IRC did not converge" in " ".join(irc_context["details"])
    points = contract["_stationary"](
        [0.0, 0.5, 0.0, 3.0, 0.0, 2.0, 0.0], path,
    )
    assert (3, "peak candidate") not in points
    assert (5, "peak candidate") not in points
    assert points == [(0, "R"), (6, "P")]


def test_colab_drop_protocol_appends_and_rejects_incomplete_batches(
    tmp_path: Path, monkeypatch,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)

    class FakeDropWidget:
        def __init__(self, generation):
            self.generation = generation
            self.messages = []

        def send(self, message):
            self.messages.append(dict(message))

    generation = app["_DROP_STATE"]["generation"]
    widget = FakeDropWidget(generation)
    pdb = (
        b"HETATM    1  C1  LIG A  10       0.000   0.000   0.000"
        b"  1.00  0.00           C\nEND\n"
    )

    def upload(batch, name, payload=pdb, *, size=None, buffers=None):
        content = {
            "event": "upload",
            "batch": batch,
            "generation": widget.generation,
            "files": [{"name": name, "size": len(payload) if size is None else size}],
        }
        app["_on_drop_upload"](
            widget, content, [payload] if buffers is None else buffers,
        )

    upload("batch-1", "structure.pdb")
    upload("batch-2", "structure.pdb")
    assert [Path(path).name for path in app["S"]["inputs"]] == [
        "structure.pdb",
        "structure_2.pdb",
    ]
    assert [message["ok"] for message in widget.messages[-2:]] == [True, True]

    before = list(app["S"]["inputs"])
    upload("bad-size", "broken.pdb", size=len(pdb) + 1)
    upload("bad-count", "missing.pdb", buffers=[])
    assert app["S"]["inputs"] == before
    assert [message["ok"] for message in widget.messages[-2:]] == [False, False]
    assert "existing files were kept" in app["input_msg"].value

    old_generation = widget.generation
    app["_bump_drop_generation"]()
    widget.generation = app["_DROP_STATE"]["generation"]
    stale = {
        "event": "upload",
        "batch": "stale",
        "generation": old_generation,
        "files": [{"name": "late.pdb", "size": len(pdb)}],
    }
    app["_on_drop_upload"](widget, stale, [pdb])
    assert widget.messages[-1]["ok"] is False
    assert app["S"]["inputs"] == before


def test_colab_local_runtime_uses_standard_widgets(monkeypatch, tmp_path: Path) -> None:
    app, _ = _execute_app(
        monkeypatch, tmp_path, parent_header={"metadata": {"colab": {"test": True}}}
    )

    assert app["IN_COLAB"] is False
    assert app["IS_COLAB_FRONTEND"] is True
    assert app["_UPLOAD_MODE"] == "basic"
    assert isinstance(app["upl"], app["W"].FileUpload)
    assert app["_RUN_LOG_INCREMENTAL"] is False


def test_colab_molstar_cluster_bonds_and_representation_state(
    monkeypatch, tmp_path: Path,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)

    def pdb_row(
        serial: int, atom: str, resid: int, x: float, y: float, z: float,
        element: str, record: str = "ATOM",
    ) -> str:
        return (
            f"{record:<6}{serial:5d} {atom:>4s} LIG A{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n"
        )

    def backbone(residues: tuple[int, ...]) -> str:
        serial = 0
        rows = []
        for resid in residues:
            for atom, element in (("N", "N"), ("CA", "C"), ("C", "C")):
                serial += 1
                rows.append(pdb_row(serial, atom, resid, float(serial * 3), 0.0, 0.0, element))
        return "".join(rows) + "END\n"

    intact = backbone((1, 2, 3, 10, 20, 21))
    scattered = backbone((1, 10, 20))
    assert app["_molstar_default_preset"](intact, "pdb") == "polymer-and-ligand"
    assert app["_molstar_default_preset"](scattered, "pdb") == "ball-and-stick"
    assert app["_molstar_default_preset"]("data_demo\n", "mmcif") == "polymer-and-ligand"

    cluster = (
        pdb_row(1, "C1", 1, 0.0, 0.0, 0.0, "C", "HETATM")
        + pdb_row(2, "O1", 1, 1.2, 0.0, 0.0, "O", "HETATM")
        + pdb_row(3, "H1", 1, 5.0, 0.0, 0.0, "H", "HETATM")
        + pdb_row(4, "H2", 1, 5.6, 0.0, 0.0, "H", "HETATM")
        + "END\n"
    )
    prepared = app["_molstar_viewer_source"](cluster, "pdb")
    assert cluster.endswith("END\n")
    assert [line for line in prepared.splitlines() if line.startswith(("ATOM", "HETATM"))] == [
        line for line in cluster.splitlines() if line.startswith(("ATOM", "HETATM"))
    ]
    assert "CONECT    1    2" in prepared
    assert "CONECT    2    1" in prepared
    assert "CONECT    3    4" not in prepared and "CONECT    4    3" not in prepared
    assert prepared.splitlines()[-1] == "END"
    existing = cluster.replace("END\n", "CONECT    1    2\nCONECT    2    1\nEND\n")
    assert app["_molstar_viewer_source"](existing, "pdb") == existing
    assert app["_molstar_viewer_source"](intact, "pdb") == intact
    assert app["_molstar_viewer_source"](cluster, "mmcif") == cluster
    assert app["_molstar_viewer_source"](cluster, "pdb", source_format="mmcif") == cluster
    multimodel = "MODEL        1\n" + cluster + "MODEL        2\n" + cluster + "ENDMDL\n"
    assert app["_molstar_viewer_source"](multimodel, "pdb") == multimodel
    duplicate_serial = (
        pdb_row(1, "C1", 1, 0.0, 0.0, 0.0, "C", "HETATM")
        + pdb_row(1, "O1", 1, 1.2, 0.0, 0.0, "O", "HETATM")
        + "END\n"
    )
    assert app["_molstar_viewer_source"](duplicate_serial, "pdb") == duplicate_serial
    digit_hydrogen = pdb_row(5, "1HG", 1, 0.0, 0.0, 0.0, "H", "HETATM").rstrip("\n")
    digit_hydrogen = digit_hydrogen[:12] + "1HG " + digit_hydrogen[16:76] + "  "
    assert app["_pdb_element_for_conect"](digit_hydrogen, {"H": 1, "Hg": 80}) == "H"

    assert dict(app["dd_rep"].options) == {"Cartoon": "cartoon", "Stick": "stick"}
    assert app["dd_rep"].description == "Representation"
    assert app["dd_rep"] in app["view_controls"].children
    assert app["dd_rep"].layout.display != "none"
    app["_apply_default_rep"](scattered, "pdb")
    assert app["dd_rep"].value == "stick" and app["S"]["rep"] == "stick"
    assert app["_rep_user_set"]["value"] is False
    app["dd_rep"].value = "cartoon"
    assert app["_rep_user_set"]["value"] is True
    app["_apply_default_rep"](scattered, "pdb")
    assert app["dd_rep"].value == "cartoon"

    app["_molstar_iframe"] = lambda *args, **kwargs: '<iframe class="rxmolstar-frame">stable</iframe>'
    app["S"].update(_pdb_text=cluster, _view_format="pdb")
    app["render_viewer"]()
    initial_generation = app["_VIEWER_GENERATION"]["value"]
    initial_frame = app["viewer_out"].value
    assert app["_VIEWER_CONTENT_KEY"]["value"][-1] == "cartoon"
    app["render_viewer"]()
    assert app["_VIEWER_GENERATION"]["value"] == initial_generation
    app["dd_rep"].value = "stick"
    assert app["_VIEWER_GENERATION"]["value"] == initial_generation + 1
    assert app["viewer_out"].value == initial_frame
    assert app["_VIEWER_CONTENT_KEY"]["value"][-1] == "stick"
    assert '"representationPreset":"ball-and-stick"' in html.unescape(
        app["viewer_signal_out"].value
    )

    legacy_payload = app["_session_dict"]()
    legacy_payload.pop("rep_user_set")
    for legacy in ("sticks", "ball+stick", "spheres", "line"):
        candidate = dict(legacy_payload, rep=legacy)
        normalized = app["_validate_and_normalize_session"](candidate)
        assert normalized["rep"] == "stick" and normalized["rep_user_set"] is True
    old_default = app["_validate_and_normalize_session"](
        dict(legacy_payload, rep="cartoon")
    )
    assert old_default["rep_user_set"] is False
    app["_apply_session"](dict(legacy_payload, rep="spheres"))
    assert app["dd_rep"].value == "stick"
    assert app["_rep_user_set"]["value"] is True

    app["_queue_change"](clear=True)
    assert app["_rep_user_set"]["value"] is False
    assert app["dd_rep"].value == "cartoon" and app["S"]["rep"] == "cartoon"

    import builtins
    real_import = builtins.__import__

    def without_ase(name, *args, **kwargs):
        if name == "ase.data":
            raise ImportError("ASE unavailable for fallback test")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_ase)
    assert app["_molstar_viewer_source"](cluster, "pdb") == cluster


def test_colab_molstar_origin_format_and_legacy_session_boundaries(
    monkeypatch, tmp_path: Path,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)

    def pdb_row(serial: int, atom: str, x: float, element: str) -> str:
        return (
            f"HETATM{serial:5d} {atom:>4s} LIG A   1    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          {element:>2s}\n"
        )

    viewer_pdb = pdb_row(1, "C1", 0.0, "C") + pdb_row(2, "O1", 1.2, "O") + "END\n"
    metadata = [
        {"serial": 1, "chain": "A", "resname": "LIG", "resseq": 1,
         "icode": "", "name": "C1"},
        {"serial": 2, "chain": "A", "resname": "LIG", "resseq": 1,
         "icode": "", "name": "O1"},
    ]
    origin = tmp_path / "origin.cif"
    origin.write_text("data_origin\n", encoding="utf-8")
    app["_load_view_structure"] = lambda _path: (
        viewer_pdb, [dict(row) for row in metadata], str(origin),
    )
    captured = {}

    def capture_iframe(source, fmt="pdb", **kwargs):
        captured.update(source=source, fmt=fmt, kwargs=dict(kwargs))
        return '<iframe class="rxmolstar-frame">stable</iframe>'

    app["_molstar_iframe"] = capture_iframe
    app["S"].update(mode="mmcif", inputs=[str(origin)], _view_input_index=0)
    app["_reset_representation_choice"]()
    assert app["build_selection"]()
    assert app["S"]["_view_format"] == "pdb"
    assert app["S"]["_view_source_format"] == "mmcif"
    assert captured["fmt"] == "pdb" and captured["kwargs"]["source_format"] == "mmcif"
    assert "CONECT" not in captured["source"]
    assert app["_VIEWER_CONTENT_KEY"]["value"][3] == "mmcif"
    app["_inline_molstar_assets"] = lambda: (None, None)
    document = app["_molstar_document"](
        viewer_pdb, "pdb", source_format="mmcif", representation="stick",
    )
    update = html.unescape(app["_molstar_update_script"](
        viewer_pdb, "pdb", 9, source_format="mmcif", representation="stick",
    ))
    assert "CONECT" not in document and "CONECT" not in update

    legacy_auto = app["_session_dict"]()
    legacy_auto.pop("rep_user_set")
    legacy_auto["rep"] = "cartoon"
    normalized_auto = app["_validate_and_normalize_session"](legacy_auto)
    assert normalized_auto["rep_user_set"] is False
    assert app["_apply_session"](legacy_auto) == []
    assert app["_rep_user_set"]["value"] is False
    assert app["dd_rep"].value == "stick"

    legacy_manual = app["_session_dict"]()
    legacy_manual.pop("rep_user_set")
    legacy_manual["rep"] = "sticks"
    assert app["_apply_session"](legacy_manual) == []
    assert app["_rep_user_set"]["value"] is True
    assert app["dd_rep"].value == "stick"

    explicit_manual = app["_session_dict"]()
    explicit_manual.update(rep="cartoon", rep_user_set=True)
    assert app["_apply_session"](explicit_manual) == []
    assert app["_rep_user_set"]["value"] is True
    assert app["dd_rep"].value == "cartoon"
    normalized_explicit_auto = app["_validate_and_normalize_session"](
        dict(explicit_manual, rep="spheres", rep_user_set=False),
    )
    assert normalized_explicit_auto["rep"] == "stick"
    assert normalized_explicit_auto["rep_user_set"] is False
    with pytest.raises(ValueError, match="rep_user_set"):
        app["_validate_and_normalize_session"](
            dict(explicit_manual, rep_user_set="false"),
        )

    app["_queue_change"](path=str(origin), remove=True)
    assert app["S"]["inputs"] == []
    assert app["_rep_user_set"]["value"] is False
    assert app["dd_rep"].value == "cartoon"
    assert app["S"]["_view_source_format"] == "pdb"


def test_colab_clear_files_clears_example_status_contract() -> None:
    app_source = _notebook()["cells"][2]["source"]
    queue_source = app_source.split("def _queue_change(", 1)[1]
    clear_branch = queue_source.split("    if clear:\n", 1)[1].split(
        "    elif path in paths:\n", 1,
    )[0]
    assert "example_status = globals().get('example_msg')" in clear_branch
    assert "if example_status is not None: example_status.value = ''" in clear_branch


def test_colab_app_executes_atomic_view_and_result_transitions(
    tmp_path: Path, monkeypatch,
) -> None:
    app, calls = _execute_app(monkeypatch, tmp_path)
    drop_generation = app["_DROP_STATE"]["generation"]
    app["example_msg"].value = "⭐ Toy system ready"
    assert app["_claim_drop_batch"]("batch-a", drop_generation) == (True, "")
    assert app["_claim_drop_batch"]("batch-a", drop_generation) == (
        False, "duplicate batch",
    )
    app["b_clear_inputs"].click()
    assert app["_DROP_STATE"]["generation"] == drop_generation + 1
    assert app["example_msg"].value == ""
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
    structure_signal = html.unescape(app["viewer_signal_out"].value)
    assert '"type":"rx-load-structure"' in structure_signal
    assert '"generation":%d' % app["_VIEWER_GENERATION"]["value"] in structure_signal
    app["dd_subcmd"].value = "scan"
    def _descendants(widget):
        yield widget
        for child in getattr(widget, "children", ()):
            yield from _descendants(child)
    clear_pair = next(
        widget for widget in _descendants(app["scan_panel"])
        if getattr(widget, "tooltip", "") == "Clear the atom pair being prepared."
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
    worker = app["_RUN_EXECUTION"].get("thread")
    if worker is not None:
        worker.join(timeout=5)
    assert app["S"]["_last_log"] == "old log"
    assert app["_RUN_STATE"]["validation_log"] == "validation transcript"
    app["_stream"] = lambda argv: (2, "invalid options")
    app["_do_validate"](None)
    worker = app["_RUN_EXECUTION"].get("thread")
    if worker is not None:
        worker.join(timeout=5)
    assert app["w_show_run_log"].value is True
    assert app["logbox"].layout.display == ""

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
    app["traj_out"].value = '<div class="rxmolstar-frame">stale viewer</div>'
    app["_load_trajectory"](str(trajectory), str(tmp_path))
    assert app["frame_slider"].max == 1 and not app["frame_slider"].disabled
    assert app["trajectory_box"].layout.display == ""
    assert "Frame 1 of 2" in app["frame_state"].value
    assert "ΔE = 0.0 kcal/mol" in app["frame_state"].value
    # The Results viewer shows this run only: one frame, no stale content.
    trajectory_viewer = app["traj_out"].value
    assert "stale viewer" not in trajectory_viewer
    assert trajectory_viewer.count('class="rxmolstar-frame"') == 1
    trajectory_document = _embedded_document(trajectory_viewer)
    assert "H 0.000 0.000 0.000" in trajectory_document
    assert "H 1.500 0.000 0.000" in trajectory_document
    viewer_generation = app["_TRAJ"]["generation"]

    calls.clear()
    app["frame_slider"].value = 1
    assert "Frame 2 of 2" in app["frame_state"].value
    assert "ΔE = 6.3 kcal/mol" in app["frame_state"].value
    # Moving the slider drives the mounted viewer; it must not rebuild it.
    assert app["traj_out"].value == trajectory_viewer
    replacement_viewers = [
        value for value in calls
        if isinstance(value, str) and 'class="rxmolstar-frame"' in value
    ]
    assert replacement_viewers == []
    frame_signal = app["traj_signal_out"].value
    signal_srcdoc = re.search(r'srcdoc="([^"]*)"', frame_signal)
    assert signal_srcdoc is not None
    assert html.unescape(signal_srcdoc.group(1)).count('"type":"rx-set-frame"') == 1
    assert app["_TRAJ"]["generation"] == viewer_generation
    assert app["_TRAJ"]["frame_message"] == {
        "type": "rx-set-frame", "generation": viewer_generation, "index": 1,
    }
    assert not app["frame_prev"].disabled and app["frame_next"].disabled

    app["_invalidate_last_run"]("Input identity changed; run again.")
    assert app["S"]["_last_manifest"] == {} and app["S"]["_last_files"] == []
    assert app["artifact_choice"].disabled and app["dl_btn"].disabled
    assert "Input identity changed" in app["results_empty"].value


def test_colab_run_fails_closed_when_implicit_validation_fails(
    tmp_path: Path, monkeypatch,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)
    streamed = []
    app["_sync_run"] = lambda: None
    app["_argv"] = lambda: [
        "mlmm", "opt", "-i", "missing.pdb", "--parm", "missing.parm7",
        "-o", "result",
    ]
    app["_stream"] = lambda argv: (streamed.append(list(argv)) or (2, "invalid"))

    app["_do_run"](None)

    assert len(streamed) == 1
    assert "--dry-run" in streamed[0]
    assert app["_RUN_STATE"]["validated_fingerprint"] is None
    assert "validation failed" in app["run_status"].value


def test_colab_manifest_hashes_inplace_input_before_execution(
    tmp_path: Path, monkeypatch,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)
    structure = tmp_path / "structure.pdb"
    structure.write_text("ORIGINAL\n", encoding="utf-8")
    original_digest = app["_sha256"](str(structure))
    argv = ["mlmm", "add-elem-info", "-i", str(structure), "--inplace"]
    app["_sync_run"] = lambda: None
    app["_argv"] = lambda: list(argv)
    app["_output_scope"] = lambda _argv: {
        "target": str(structure), "root": str(tmp_path), "shallow": True,
        "stdout_only": False, "direct_current": True,
    }
    app["_output_scope_collision"] = lambda _scope: False
    app["_snapshot_output_scope"] = lambda _scope: {}
    app["_begin_results_attempt"] = lambda *_args: None
    app["_results"] = lambda *_args: None

    def mutate_input(_argv):
        structure.write_text("MUTATED\n", encoding="utf-8")
        return 0, "updated in place"

    app["_stream"] = mutate_input
    app["_do_run"](None)

    assert app["S"]["_last_manifest"]["inputs"] == [{
        "path": str(structure), "sha256": original_digest,
    }]
    assert app["_sha256"](str(structure)) != original_digest


def test_colab_xyz_reference_inputs_are_transactional_and_stay_paired(
    tmp_path: Path, monkeypatch,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)
    xyz_a = tmp_path / "reactant.xyz"
    xyz_b = tmp_path / "product.xyz"
    ref_a = tmp_path / "reactant_ref.pdb"
    ref_b = tmp_path / "product_ref.pdb"
    topology = tmp_path / "system.parm7"
    xyz_a.write_text(
        "2\nreactant geometry\n"
        "C 0.000 0.000 0.000\nO 1.200 0.000 0.000\n",
        encoding="utf-8",
    )
    xyz_b.write_text(
        "2\nproduct geometry\n"
        "C 0.100 0.000 0.000\nO 1.300 0.000 0.000\n",
        encoding="utf-8",
    )
    ref_template = (
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  1.00           C\n"
        "HETATM    2  O1  LIG A   1       1.200   0.000   0.000  1.00  1.00           O\n"
        "END\n"
    )
    ref_a.write_text(ref_template, encoding="utf-8")
    ref_b.write_text(ref_template, encoding="utf-8")
    topology.write_text("%VERSION VERSION_STAMP = V0001.000\n", encoding="utf-8")

    assert app["_ingest_saved_files"]([str(xyz_a)], "XYZ first")
    assert app["S"]["mode"] == "xyz"
    assert app["S"]["inputs"] == [str(xyz_a)]
    assert app["S"]["ref_pdbs"] == []
    assert app["_ingest_saved_files"](
        [str(ref_a), str(topology)], "reference and topology"
    )
    assert app["S"]["ref_pdbs"] == [str(ref_a)]
    assert app["S"]["parm"] == str(topology)
    assert app["S"]["_pdb_text"]

    app["set_subcmd"]("path-opt")
    assert app["_ingest_saved_files"]([str(xyz_b), str(ref_b)], "second endpoint")
    assert app["S"]["subcmd"] == "path-opt"
    assert app["S"]["inputs"] == [str(xyz_a), str(xyz_b)]
    assert app["S"]["ref_pdbs"] == [str(ref_a), str(ref_b)]

    app["_queue_change"](path=str(xyz_b), delta=-1)
    assert app["S"]["inputs"] == [str(xyz_b), str(xyz_a)]
    assert app["S"]["ref_pdbs"] == [str(ref_b), str(ref_a)]
    app["w_q"].value = 0
    app["w_charge_ok"].value = True
    command = app["build_cmd"]()
    assert command[1] == "path-opt"
    assert command[command.index("-i") + 1 : command.index("-o")] == [
        str(xyz_b), str(xyz_a), "-b", app["BACKEND"],
    ]
    emitted_refs = [
        command[index + 1]
        for index, token in enumerate(command)
        if token == "--ref-pdb"
    ]
    assert emitted_refs == [str(ref_b), str(ref_a)]
    assert command[command.index("--parm") + 1] == str(topology)
    assert command[command.index("-q") + 1] == "0"

    app["_queue_change"](path=str(xyz_b), remove=True)
    assert app["S"]["inputs"] == [str(xyz_a)]
    assert app["S"]["ref_pdbs"] == [str(ref_a)]

    malformed = tmp_path / "bad.xyz"
    malformed.write_text("2\nbad\nC 0 0 0\n", encoding="utf-8")
    before = (
        list(app["S"]["inputs"]),
        list(app["S"]["ref_pdbs"]),
        app["S"]["subcmd"],
    )
    assert not app["_ingest_saved_files"]([str(malformed)], "bad XYZ")
    assert (
        app["S"]["inputs"],
        app["S"]["ref_pdbs"],
        app["S"]["subcmd"],
    ) == before


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
    uploaded_parm_command = app["build_cmd"]()
    assert uploaded_parm_command[uploaded_parm_command.index("--parm") + 1] == str(topology)
    app["S"]["parm"] = None
    generated_parm_command = app["build_cmd"]()
    assert "--parm" not in generated_parm_command
    app["S"]["parm"] = str(topology)
    assert len(app["input_file_rows"].children) == 3
    assert app["prep_radius"].value == pytest.approx(2.6)
    assert app["adv_radius"].value == pytest.approx(2.6)
    assert app["dd_col"].value == "element"
    assert app["dd_rep"].value == "stick"
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
    assert app["w_show_run_log"].value is True
    assert app["logbox"].layout.display == ""
    app["w_show_run_log"].value = False
    assert app["logbox"].layout.display == "none"
    app["w_show_run_log"].value = True
    assert app["_input_box_children"][1] is app["_drop"]
    assert app["_input_box_children"][2] is app["input_msg"]
    assert app["_input_box_children"][-1] is app["example_fold"]
    assert app["example_fold"]._rx_body.layout.display == "none"
    assert app["example_fold"]._rx_button.description == "Show Examples"
    assert "rxfold-toggle" in app["example_fold"]._rx_button._dom_classes
    app["example_fold"]._rx_button.click()
    assert app["example_fold"]._rx_body.layout.display == ""
    assert app["example_fold"]._rx_button.description == "Hide Examples"
    assert "rxfold-open" in app["example_fold"]._rx_button._dom_classes
    assert "Click sets:" in app["last_pick_html"].value
    assert "#E0DDD4" in app["last_pick_html"].value
    assert "#F3F2EE" in app["last_pick_html"].value
    assert "#332B1F" in app["last_pick_html"].value
    assert "#B66D1D" in app["last_pick_html"].value
    assert "rxcard" in app["depth_box"]._dom_classes
    assert set(app["_center_values"]()) == {"LIG", "MG"}
    mg_row = app["charge_rows"]["MG"]
    assert mg_row["auto"] and mg_row["use"].disabled and mg_row["val"].disabled
    assert mg_row["val"].value == 2
    assert "MG" not in app["S"]["lcharge"]

    app["center_widget"].value = ("MG",)
    with pytest.raises(ValueError, match="Confirm each ligand charge"):
        app["build_cmd"]()
    app["charge_rows"]["LIG"]["use"].value = True
    # -l is a defensible charge source, so the CLI derives the ML-region charge
    # and the GUI emits no -q until the user opts into an explicit override.
    derived_command = app["build_cmd"]()
    assert "-q" not in derived_command
    assert derived_command[derived_command.index("-l") + 1] == "LIG:0,MG:2"
    # Without -l, the visible system charge is emitted directly.
    app["center_widget"].value = ()
    app["charge_rows"]["LIG"]["use"].value = False
    app["charge_rows"]["MG"]["use"].value = False
    direct_charge_command = app["build_cmd"]()
    assert direct_charge_command[direct_charge_command.index("-q") + 1] == "0"
    app["charge_rows"]["MG"]["use"].value = True
    app["charge_rows"]["LIG"]["use"].value = True
    app["center_widget"].value = ("MG",)
    app["w_q"].value = -1
    app["w_charge_ok"].value = True
    mg_row["val"].value = 3  # disabled widgets are still mutable from Python
    assert app["w_charge_ok"].value is True
    charged_command = app["build_cmd"]()
    assert charged_command[charged_command.index("-q") + 1] == "-1"
    ligand_values = dict(
        item.split(":", 1)
        for item in charged_command[charged_command.index("-l") + 1].split(",")
    )
    assert ligand_values == {"LIG": "0", "MG": "2"}
    stage_scope = app["_charge_scope_fingerprint"]("sp")
    assert app["_charge_scope_fingerprint"]("opt") == stage_scope
    assert app["_charge_scope_fingerprint"]("all") != stage_scope
    mg_row["val"].value = 2
    app["charge_rows"]["LIG"]["use"].value = False
    app["center_widget"].value = ()
    app["w_charge_ok"].value = False

    assert app["S"]["show_water"] is True
    app["cb_water"].value = False
    app["cb_water"].value = True
    assert app["S"]["show_water"] is True
    document = app["_molstar_document"](
        app["S"]["_pdb_text"], "pdb", show_water=True, interactive=True,
    )
    assert '"showWater":true' in document
    assert "layoutShowControls:true" in document
    assert "layoutShowSequence:cfg.showSequence" in document
    assert "structureRequestKey=JSON.stringify" in document
    assert "if(requestKey===structureRequestKey)return;" in document
    assert "structure-component-static-water" in document
    assert "WebGL is unavailable in this browser session." in document
    partial_pdb = "ATOM      1   CB ALA A   1      0.000   0.000   0.000\n"
    intact_pdb = "".join(
        f"ATOM  {serial:5d} {atom:>4s} ALA A{resid:4d}    0.000\n"
        for serial, (resid, atom) in enumerate(
            ((resid, atom) for resid in range(1, 4) for atom in ("N", "CA", "C")),
            start=1,
        )
    )
    assert app["_molstar_default_preset"](partial_pdb, "pdb") == "ball-and-stick"
    assert app["_molstar_default_preset"](intact_pdb, "pdb") == "polymer-and-ligand"
    partial_document = app["_molstar_document"](partial_pdb, "pdb")
    assert '"representationPreset":"ball-and-stick"' in partial_document
    assert "representationPreset:'empty'" in partial_document
    assert "component,{type:'ball-and-stick'}" in partial_document
    intact_document = app["_molstar_document"](intact_pdb, "pdb")
    assert '"representationPreset":"polymer-and-ligand"' in intact_document
    partial_update = html.unescape(app["_molstar_update_script"](
        partial_pdb, "pdb", 7, show_water=False,
    ))
    assert '"representationPreset":"ball-and-stick"' in partial_update
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
    assert [pick["index"] for pick in app["S"]["_pick_history"]] == [4]

    previous_signal = app["viewer_signal_out"].value
    app["on_click"]("2", "MG", "2", "A", "MG", "3", "", live_marked=True)
    assert [pick["index"] for pick in app["S"]["_pick_history"]] == [4, 2]
    assert app["viewer_signal_out"].value == previous_signal
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
    assert "tr_projection" not in all_statuses
    assert all_statuses["verbose"] == "rendered"
    assert all_statuses["embedcharge"] == "rendered"
    assert all_statuses["embedcharge_cutoff"] == "rendered"
    assert all_statuses["detect_layer"] == "rendered"
    assert app["_advanced_coverage"]("dft")
    assert app["_advanced_coverage"]("sp")["hessian_calc_mode"] == "rendered"
    print_every_param = next(
        param for param in app["_advanced_options"]("all")
        if param.name == "print_every"
    )
    assert app["_advanced_status"]("all", print_every_param) == "rendered"
    app["S"]["advanced_overrides"]["all"] = {"print_every": 7}
    print_every_argv = app["_advanced_argv"]("all")
    assert print_every_argv[print_every_argv.index("--print-every") + 1] == "7"
    verbose_param = next(
        param for param in app["_advanced_options"]("all")
        if param.name == "verbose"
    )
    app["S"]["advanced_overrides"]["all"] = {"verbose": "7"}
    verbose_row = app["_advanced_widget"]("all", verbose_param)
    assert verbose_row.children[0].value == app["_ADVANCED_DEFAULT"]
    assert "verbose" not in app["S"]["advanced_overrides"]["all"]

    app["dd_subcmd"].options = app["_sub_options"](app["SUBS"])
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
    ]
    rendered_rows = [child for child in app["adv_rows_box"].children
                     if "rxflagrow" in getattr(child, "_dom_classes", ())]
    assert len(rendered_rows) == len(editable)
    assert "Standard flags" in app["adv_rows_box"].children[0].value
    assert any("Advanced flags" in getattr(child, "value", "")
               for child in app["adv_rows_box"].children)

    click = app["click"]
    text_param = next(
        param for param in editable
        if param.name != "verbose"
        and not param.multiple
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

    app["dd_subcmd"].value = "energy-diagram"
    app["S"]["advanced_overrides"]["energy-diagram"] = {}
    app["refresh"]()
    assert not app["_ACTION_STATE"]["auto_ready"]
    assert app["b_run"].disabled
    assert "Advanced flags" in app["subcmd_note"].value
    for invalid_values in (
        {"input_values": "0"},
        {"input_values": "0 oops"},
        {"input_values": "0 1", "label_x": "R"},
    ):
        app["S"]["advanced_overrides"]["energy-diagram"] = invalid_values
        invalid_command = app["build_cmd"]()
        assert not app["_utility_autofill_complete"](invalid_command)
        app["refresh"]()
        assert app["b_run"].disabled
    app["S"]["advanced_overrides"]["energy-diagram"] = {
        "input_values": "0 12.5 4.3",
        "output_path": str(tmp_path / "profile.png"),
    }
    energy_command = app["build_cmd"]()
    assert energy_command.count("--input") == 3
    assert app["_utility_autofill_complete"](energy_command)
    app["refresh"]()
    assert app["_ACTION_STATE"]["auto_ready"]
    assert not app["b_run"].disabled

    atom_a = {"index": 0, "chain": "A", "resn": "ALA", "resi": "1", "atom": "N"}
    atom_b = {"index": 1, "chain": "A", "resn": "ALA", "resi": "1", "atom": "H"}
    app["S"]["freeze_pairs"] = [{"a": atom_a, "b": atom_b, "t": None}]
    app["dd_subcmd"].value = "opt"
    opt_picks = {value for _label, value in app["pick_action"].options}
    assert {"freezeA", "freezeB", "freezeatom"} <= opt_picks
    assert "Distance restraints" in " ".join(
        getattr(child, "value", "") for child in app["freeze_panel"].children
    )

    app["dd_subcmd"].value = "sp"
    assert app["key_opts_box"].layout.display == ""
    assert [value for _label, value in app["pick_action"].options] == [
        "center", "ligand", "selectedresn", "freezeatom",
    ]
    assert "Distance restraints" not in " ".join(
        getattr(child, "value", "") for child in app["freeze_panel"].children
    )
    assert len(app["S"]["freeze_pairs"]) == 1
    app["_render_summary"]()
    assert "freeze:" not in app["summary_html"].value
    app["center_widget"].value = ()
    app["S"]["center_ids"] = []
    app["charge_rows"]["LIG"]["use"].value = False
    app["w_q"].value = 0
    app["w_charge_ok"].value = True
    assert "--dist-freeze" not in app["build_cmd"]()
    app["dd_subcmd"].value = "opt"
    assert "Distance restraints" in " ".join(
        getattr(child, "value", "") for child in app["freeze_panel"].children
    )
    app["_render_summary"]()
    assert "freeze:" in app["summary_html"].value
    app["dd_subcmd"].value = "sp"
    assert app["center_panel"].layout.display == ""
    assert app["charge_panel"].layout.display == ""
    assert app["extract_panel"].layout.display == ""
    assert app["adv_radius"]._rx_flag_row.layout.display == "none"
    assert app["adv_dftfb"]._rx_flag_row.layout.display == "none"
    app["dd_subcmd"].value = "freq"
    assert app["key_opts_box"].layout.display == ""
    app["dd_subcmd"].options = app["_sub_options"](app["SUBS"])
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
    app["charge_rows"]["LIG"]["use"].value = True
    app["w_charge_ok"].value = True
    extract_commands = []

    def fake_extract(command, **_kwargs):
        extract_commands.append(list(command))
        output = Path(command[command.index("-o") + 1])
        output.write_text(primary.read_text(encoding="utf-8"), encoding="utf-8")
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(app["subprocess"], "run", fake_extract)
    app["b_extract"].click()
    app["_EXTRACT_TASK"]["thread"].join(timeout=5)
    assert extract_commands
    assert extract_commands[-1][extract_commands[-1].index("-r") + 1] == "4.2"
    assert extract_commands[-1][extract_commands[-1].index("-l") + 1] == "LIG:0,MG:2"
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
    app["w_th"].value = False
    app["adv_dft"].value = False
    app["w_ts"].value = False
    assert not app["w_ts"].value and not app["w_ts"].disabled
    # Advanced controls remain discoverable and retain explicit user values
    # even while their optional stage is disabled.  Re-enabling the stage must
    # not silently discard those values.
    off_argv = app["_advanced_argv"]("all")
    for flag in (
        "--irc-step-size", "--opt-mode-post", "--thresh-post",
        "--hessian-calc-mode", "--skip-final-freq", "--no-reject-uphill",
    ):
        assert flag in off_argv
    app["all_mode"].value = "tsonly"
    assert app["adv_refine"].layout.display == "none"
    assert app["w_ts"].value and not app["w_ts"].disabled
    app["w_charge_ok"].value = True
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
    assert "--no-preopt" in command and "--preopt" not in command
    assert "--tsopt" in command
    app["all_mode"].value = "mep"
    assert not app["w_ts"].value and not app["w_ts"].disabled
    app["w_th"].value = True
    assert app["w_ts"].value and not app["w_ts"].disabled
    app["w_ts"].value = False
    assert not app["w_ts"].value and not app["w_ts"].disabled
    assert not app["w_th"].value
    app["w_th"].value = True
    assert app["w_ts"].value and not app["w_ts"].disabled
    thermo_argv = app["_advanced_argv"]("all")
    for flag in (
        "--irc-step-size", "--opt-mode-post", "--thresh-post",
        "--hessian-calc-mode",
    ):
        assert flag in thermo_argv
    assert "--no-reject-uphill" in thermo_argv
    assert "--skip-final-freq" in thermo_argv
    app["S"]["inputs"] = [str(primary), str(secondary)]
    app["w_charge_ok"].value = False
    app["w_charge_ok"].value = True
    thermo_cmd = app["build_cmd"]()
    assert thermo_cmd.count("--tsopt") == 1
    assert thermo_cmd.count("--thermo") == 1
    app["w_th"].value = False
    assert app["w_ts"].value and not app["w_ts"].disabled
    app["w_ts"].value = False
    assert "--tsopt" not in app["build_cmd"]()
    app["S"]["tsopt"] = False
    app["S"]["thermo"] = True
    with pytest.raises(ValueError, match="require TS optimization"):
        app["build_cmd"]()
    app["S"]["thermo"] = False

    result_json = tmp_path / "result.json"
    result_json.write_text('{"energy": -1.25}', encoding="utf-8")
    assert app["_artifact_kind"](str(result_json)) == "JSON"
    assert app["_artifact_kind"]("job.gjf") == "text"
    assert app["_artifact_kind"]("job.com") == "text"
    assert app["_artifact_kind"]("job.inp") == "text"
    assert app["_artifact_kind"]("forcefield.prms") == "text"
    preview = app["_text_preview_html"](str(result_json), "JSON")
    assert "&quot;energy&quot;" in preview and "-1.25" in preview
    assert "DejaVu Sans Mono" in preview and "Liberation Mono" in preview
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps({
        "status": "success", "scientific_status": "partial",
        "scientific_status_reasons": ["IRC endpoint mismatch"],
        "segments": [{"index": 1, "barrier_kcal": 8.0, "delta_kcal": -1.0}],
        "post_segments": [{"index": 1, "mlip": {"barrier_kcal": 7.5, "delta_kcal": -1.2}}],
        "rate_limiting_step": {"barrier_kcal": 7.5, "method": "mlip"},
    }), encoding="utf-8")
    summary_html = app["_summary_html"](str(summary))
    assert "Highest local barrier: 7.5" in summary_html
    assert "IRC endpoint mismatch" in summary_html
    assert "MLIP" in summary_html and "partial result" not in summary_html
    assert ">partial</span>" not in summary_html
    ts_only_summary = tmp_path / "ts_only_summary.json"
    ts_only_summary.write_text(json.dumps({
        "status": "success",
        "scientific_status": "success",
        "pipeline_mode": "tsopt-only",
        "n_images": 5,
        "mlip_backend": "mace",
        "mlip_model": "MACE-OMOL-0",
        "mlip_model_label": "MACE-OMOL-0",
        "endpoint_assignment": {"chemical_direction_known": False},
        "segments": [{
            "index": 1, "kind": "tsopt",
            "barrier_kcal": 8.0, "delta_kcal": -1.0,
        }],
        "post_segments": [{
            "index": 1,
            "dft": {"barrier_kcal": 7.8, "delta_kcal": -1.2},
        }],
        "rate_limiting_step": {
            "segment": 1, "barrier_kcal": 7.8, "method": "DFT",
        },
    }), encoding="utf-8")
    ts_only_html = app["_summary_html"](str(ts_only_summary))
    assert "raw MEP" not in ts_only_html
    assert "IRC frames: 5" not in ts_only_html
    assert "backend/model:" not in ts_only_html
    assert "MLIP model: <b>MACE-OMOL-0</b>" in ts_only_html
    assert "model-region DFT" in ts_only_html
    assert "energy order" in ts_only_html
    extra_artifacts = []
    for name in ("a.txt", "b.txt", "c.txt", "z.txt"):
        artifact = tmp_path / name
        artifact.write_text(name, encoding="utf-8")
        extra_artifacts.append(str(artifact))
    app["S"].update(
        _last_out_dir=str(tmp_path),
        _last_subcmd="all",
        _last_files=[str(summary), *extra_artifacts],
        _last_manifest={"status": "success", "exit_code": 0},
    )
    context_html = app["_result_context_html"](str(tmp_path))
    assert "IRC endpoint mismatch" not in context_html
    assert "<b>all</b>" in context_html
    assert context_html.count("<code>") == len(extra_artifacts) + 1
    calls.clear()
    app["_structure_preview"](str(primary))
    assert any(
        isinstance(value, str) and 'class="rxmolstar-frame"' in value
        and "Mol*" in value
        for value in calls
    )

    # The Setup page keeps one explicit system-charge contract across workflows.
    app["dd_subcmd"].value = "all"
    assert app["w_q"].description == "system charge (-q)"
    assert "Ligand charges calculate the ML-region charge automatically" in app["charge_info"].value
    assert "overwrite system charge" in app["charge_info"].value
    app["dd_subcmd"].value = "sp"
    assert app["w_q"].description == "system charge (-q)"
    app["dd_subcmd"].options = app["_sub_options"](app["SUBS"])
    app["dd_subcmd"].value = "oniom-export"
    assert app["w_q"].description == "system charge (-q)"

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
    assert mm_parm_command[1] == "mm-parm"
    assert mm_parm_command[-2:] == ["-l", "LIG:1,MG:2"]

    # ONIOM export chooses an output extension that matches its selected engine.
    app["S"].update(inputs=[str(primary)], parm=str(topology), model_pdb=None)
    app["dd_subcmd"].value = "oniom-export"
    app["S"]["advanced_overrides"]["oniom-export"] = {"mode": "orca"}
    app["w_q"].value = 0
    app["w_charge_ok"].value = False
    app["w_charge_ok"].value = True
    orca_export = app["build_cmd"]()
    assert Path(orca_export[orca_export.index("-o") + 1]).suffix == ".inp"
    assert orca_export[orca_export.index("--mode") + 1] == "orca"
    app["S"]["advanced_overrides"]["oniom-export"] = {"mode": "g16"}
    assert app["w_charge_ok"].value is True
    gaussian_export = app["build_cmd"]()
    assert Path(gaussian_export[gaussian_export.index("-o") + 1]).suffix == ".gjf"

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


def test_colab_exercises_every_workflow_and_advanced_flag_widget(
    tmp_path: Path, monkeypatch,
) -> None:
    """Operate every workflow selector and every editable live Click option."""
    app, _ = _execute_app(monkeypatch, tmp_path)
    for subcommand in (
        "all", "opt", "tsopt", "irc", "freq", "sp",
        "scan", "scan2d", "scan3d", "path-opt", "path-search", "dft",
    ):
        param = next(
            option
            for option in app["_advanced_options"](subcommand)
            if option.name == "embedcharge"
        )
        assert app["_advanced_status"](subcommand, param) == "rendered"
        help_text = param.help.lower()
        assert "experimental" in help_text
        if subcommand == "dft":
            assert "pyscf qm hamiltonian" in help_text
        else:
            assert "xtb" in help_text and "delta correction" in help_text
            assert "computationally expensive" in help_text
        assert "experimental" in app["_advanced_widget"](
            subcommand, param
        )._rx_search
    app["dd_subcmd"].options = app["_sub_options"](app["SUBS"])
    option_values = {
        item[1] if isinstance(item, tuple) else item
        for item in app["dd_subcmd"].options
    }
    assert option_values == set(app["SUBS"])

    for subcommand in sorted(option_values):
        app["dd_subcmd"].value = subcommand
        assert app["S"]["subcmd"] == subcommand
        assert app["subreq"].value and app["outputs_html"].value
        app["_render_advanced_rows"]()

    root = app["PRODUCT_CLI"]
    root_context = app["click"].Context(root)
    default_dropdowns: set[str] = set()
    for subcommand in root.list_commands(root_context):
        rendered = {
            param.name for param in app["_advanced_options"](subcommand)
            if app["_advanced_status"](subcommand, param) == "rendered"
        }
        exercised: set[str] = set()
        scenarios = [None]
        if subcommand == "all":
            scenarios = [
                ("mep", False, False, False),
                ("scan", False, False, False),
                ("tsonly", True, False, False),
                ("mep", True, True, True),
            ]
        for scenario in scenarios:
            if scenario is not None:
                mode, tsopt, thermo, dft = scenario
                app["all_mode"].value = mode
                app["w_ts"].value = tsopt
                app["w_th"].value = thermo
                app["adv_dft"].value = dft
            for param in app["_advanced_options"](subcommand):
                if app["_advanced_status"](subcommand, param) != "rendered":
                    continue
                app["S"].setdefault("advanced_overrides", {})[subcommand] = {}
                row = app["_advanced_widget"](subcommand, param)
                widget = row.children[0]
                if isinstance(widget, app["W"].Dropdown):
                    assert widget.index == 0, param.name
                    assert widget.value == app["_ADVANCED_DEFAULT"], param.name
                    assert widget.label.startswith("default:"), param.name
                    default_dropdowns.add(param.name)
                for value in _advanced_widget_values(app, param):
                    widget.value = value
                    assert (
                        app["S"]["advanced_overrides"][subcommand][param.name]
                        == widget.value
                    )
                    argv = app["_advanced_argv"](subcommand)
                    _assert_advanced_argv_parses(app, subcommand, argv)
                # Nullable numeric controls cannot accept an empty string.
                # Reset through the one authoritative override API; the next
                # render then recreates the control at its CLI default.
                app["_set_advanced_override"](subcommand, param.name, None)
                assert (
                    param.name
                    not in app["S"]["advanced_overrides"][subcommand]
                )
                exercised.add(param.name)
        assert exercised == rendered, (
            subcommand,
            sorted(rendered - exercised),
        )
    assert {"verbose", "hessian_calc_mode"} <= default_dropdowns


def test_colab_operates_every_workflow_through_validate_run_and_results(
    tmp_path: Path, monkeypatch,
) -> None:
    """Build, parse, validate, run, and render every GUI workflow state."""
    app, _ = _execute_app(monkeypatch, tmp_path)
    pdb_a = tmp_path / "reactant.pdb"
    pdb_b = tmp_path / "product.pdb"
    parm = tmp_path / "system.parm7"
    trajectory = tmp_path / "path.xyz"
    oniom_input = tmp_path / "job.gjf"
    pdb_text = (
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  LIG A   1       1.200   0.000   0.000  1.00  0.00           O\n"
        "HETATM    3  N1  LIG A   1       0.000   1.200   0.000  1.00  0.00           N\nEND\n"
    )
    pdb_a.write_text(pdb_text, encoding="utf-8")
    pdb_b.write_text(pdb_text.replace("1.200", "1.300"), encoding="utf-8")
    parm.write_text("minimal topology fixture\n", encoding="utf-8")
    trajectory.write_text(
        "2\n0.0\nH 0 0 0\nH 0 0 1\n2\n1.0\nH 0 0 0\nH 0 0 1.1\n",
        encoding="utf-8",
    )
    oniom_input.write_text("# oniom\n\nmatrix\n\n", encoding="utf-8")

    app["DFT_READY"] = True
    if "dft" not in app["SUBS"]:
        app["SUBS"].append("dft")
    app["dd_subcmd"].options = app["_sub_options"](app["SUBS"])
    app["_ambertools_available"] = lambda: True

    streamed: list[list[str]] = []
    results: list[tuple[str, str]] = []
    app["_stream"] = lambda argv: (
        streamed.append(list(argv)) or (0, "operation-matrix success")
    )
    app["_output_scope"] = lambda argv: {
        "target": str(tmp_path / ("out-" + argv[1])),
        "root": str(tmp_path / ("out-" + argv[1])),
        "shallow": False,
        "stdout_only": False,
        "direct_current": False,
    }
    app["_output_scope_collision"] = lambda _scope: False
    app["_snapshot_output_scope"] = lambda _scope: {}
    app["_begin_results_attempt"] = lambda root, sub: results.append((root, sub))
    app["_results"] = lambda root: results.append((root, "rendered"))

    atoms = [
        {"index": index, "chain": "A", "resn": "LIG", "resi": "1",
         "atom": name, "xyz": xyz}
        for index, (name, xyz) in enumerate((
            ("C1", (0.0, 0.0, 0.0)),
            ("O1", (1.2, 0.0, 0.0)),
            ("N1", (0.0, 1.2, 0.0)),
        ))
    ]
    primary_atom_meta = [
        {
            "index": atom["index"], "chain": atom["chain"],
            "resname": atom["resn"], "resseq": atom["resi"], "icode": "",
            "atom": atom["atom"], "xyz": atom["xyz"],
        }
        for atom in atoms
    ]
    cases = [
        ("all", "mep"), ("all", "scan"), ("all", "tsonly"),
        ("opt", None), ("sp", None), ("tsopt", None), ("freq", None),
        ("irc", None), ("dft", None), ("scan", None), ("scan2d", None),
        ("scan3d", None), ("path-opt", None), ("path-search", None),
        ("extract", None), ("fix-altloc", None), ("add-elem-info", None),
        ("energy-diagram", None), ("bond-summary", None), ("trj2fig", None),
        ("define-layer", None), ("mm-parm", None),
        ("oniom-export", None), ("oniom-import", None),
    ]

    for index, (subcommand, all_kind) in enumerate(cases):
        streamed.clear()
        app["S"].update(
            mode="pdb", inputs=[str(pdb_a)], ref_pdbs=[], parm=str(parm),
            model_pdb=None, center=[], center_ids=[], lcharge={},
            scan_atoms=[None, None], scan_stages=[], scan_axes=[],
            freeze_buf=[None, None], freeze_pairs=[], freeze_atoms=[],
            advanced_overrides={}, _session_file_identities={},
            _installed_backend="mace", backend="mace",
            _primary_atom_meta=primary_atom_meta,
        )
        app["w_ts"].value = False
        app["w_th"].value = False
        app["adv_dft"].value = False
        app["w_out"].value = "matrix-%02d-%s" % (index, subcommand)
        app["dd_subcmd"].value = subcommand
        app["all_mode"].value = all_kind or "mep"

        if subcommand in {"all", "path-opt", "path-search", "bond-summary"}:
            app["S"]["inputs"] = [str(pdb_a), str(pdb_b)]
        if subcommand == "all" and all_kind in {"scan", "tsonly"}:
            app["S"]["inputs"] = [str(pdb_a)]
        if subcommand == "all" and all_kind == "scan":
            bond = {"a": atoms[0], "b": atoms[1], "t": 1.5}
            app["S"]["scan_atoms"] = [atoms[0], atoms[1]]
            app["S"]["scan_stages"] = [[bond]]
        if subcommand == "scan":
            bond = {"a": atoms[0], "b": atoms[1], "t": 1.5}
            app["S"]["scan_atoms"] = [atoms[0], atoms[1]]
            app["S"]["scan_stages"] = [[bond]]
        if subcommand in {"scan2d", "scan3d"}:
            axes = [
                {"a": atoms[0], "b": atoms[1], "lo": 1.0, "hi": 2.0},
                {"a": atoms[0], "b": atoms[2], "lo": 1.0, "hi": 2.0},
                {"a": atoms[1], "b": atoms[2], "lo": 1.0, "hi": 2.0},
            ]
            app["S"]["scan_axes"] = axes[:2 if subcommand == "scan2d" else 3]
        if subcommand == "extract":
            app["S"]["center_ids"] = ["A:LIG:1"]
        elif subcommand == "energy-diagram":
            app["S"]["inputs"] = []
            app["S"]["mode"] = "utility"
            app["S"]["advanced_overrides"][subcommand] = {
                "input_values": "0 5", "output_path": str(tmp_path / "energy.png"),
            }
        elif subcommand == "trj2fig":
            app["S"]["inputs"] = [str(trajectory)]
            app["S"]["mode"] = "xyz"
        elif subcommand == "define-layer":
            app["S"]["advanced_overrides"][subcommand] = {
                "model_indices_str": "1,2",
            }
        elif subcommand == "oniom-export":
            app["S"]["model_pdb"] = str(pdb_a)
            app["S"]["advanced_overrides"][subcommand] = {"mode": "orca"}
        elif subcommand == "oniom-import":
            app["S"]["inputs"] = [str(oniom_input)]
            app["S"]["parm"] = None
            app["S"]["mode"] = "utility"

        if subcommand in app["COMPUTE"] or subcommand == "oniom-export":
            app["w_q"].value = 0
            app["w_charge_ok"].value = False
            app["w_charge_ok"].value = True
        command = app["build_cmd"]()
        assert command[:2] == ["mlmm", subcommand]
        click_command = app["PRODUCT_CLI"].get_command(
            app["click"].Context(app["PRODUCT_CLI"]), subcommand,
        )
        parser_argv = _root_normalized_subcommand_argv(
            app, subcommand, command[2:],
        )
        with click_command.make_context(
            subcommand, parser_argv, resilient_parsing=False,
        ) as context:
            expected_extra = (
                app["S"]["inputs"][1:]
                if subcommand in {"all", "path-search"} else []
            )
            assert list(context.args) == expected_extra, (
                subcommand, command, context.args,
            )

        app["_auto"]["on"] = False
        app["cmd_box"].value = shlex.join(command)
        if subcommand in app["COMPUTE"]:
            assert not app["b_validate"].disabled
            app["b_validate"].click()
            assert streamed and "--dry-run" in streamed[0]
        else:
            before_validate = len(streamed)
            app["b_validate"].click()
            assert len(streamed) == before_validate
        assert not app["b_run"].disabled, (subcommand, app["run_status"].value)
        app["b_run"].click()
        assert streamed[-1] == command
        assert app["S"]["_last_manifest"]["status"] == "success"
        assert app["S"]["_last_manifest"]["subcommand"] == subcommand
        assert results[-1][1] == "rendered"


def test_colab_results_keep_missing_energies_unknown(
    tmp_path: Path, monkeypatch,
) -> None:
    app, rendered = _execute_app(monkeypatch, tmp_path)
    missing_first = tmp_path / "missing_first_trj.xyz"
    missing_first.write_text(
        "2\nenergy unavailable\nH 0 0 0\nH 0.7 0 0\n"
        "2\nenergy=-0.990000 unit=hartree\nH 0.2 0 0\nH 0.9 0 0\n",
        encoding="utf-8",
    )
    app["S"]["_last_subcmd"] = "path-opt"
    rendered.clear()
    app["_load_trajectory"](str(missing_first), str(tmp_path))
    assert app["_TRAJ"]["energies"] == [None, -0.99]
    assert app["_rel_kcal"]() is None
    assert "profile not re-referenced" in app["frame_state"].value
    assert app["plot_out"].value == ""
    assert app["energy_panel"].layout.display == "none"
    assert "rxstructure-only" in app["path_grid"]._dom_classes

    missing_middle = tmp_path / "missing_middle_trj.xyz"
    missing_middle.write_text(
        "2\nenergy=-1.000000 unit=hartree\nH 0 0 0\nH 0.7 0 0\n"
        "2\nenergy unavailable\nH 0.2 0 0\nH 0.9 0 0\n",
        encoding="utf-8",
    )
    rendered.clear()
    app["_load_trajectory"](str(missing_middle), str(tmp_path))
    assert app["_rel_kcal"]() == [0.0, None]
    app["frame_slider"].value = 1
    assert "ΔE unavailable" in app["frame_state"].value
    assert app["plot_out"].value.count('class="rxenergy-frame"') == 1
    assert "stale energy view" not in app["plot_out"].value


def test_colab_defers_active_browser_input_deletion(
    tmp_path: Path, monkeypatch,
) -> None:
    app, _rendered = _execute_app(monkeypatch, tmp_path)
    source = _notebook()["cells"][2]["source"]
    # One run path: the unreachable async twin was removed.
    run_path = source[source.index("def _do_run_sync"):source.index("def _do_run(")]
    assert "_RUN_EXECUTION['argv'] = list(a)\n    _set_running(True)" in run_path
    sync_start = source.index("def _do_run_sync")
    sync_prefix = source[sync_start:source.index("    try:\n", sync_start)]
    assert "_RUN_EXECUTION['argv']" not in sync_prefix

    detached = tmp_path / "detached.pdb"
    detached.write_text("END\n", encoding="utf-8")
    app["_remember_uploaded_paths"]([str(detached)])
    app["S"]["inputs"] = [str(detached)]
    app["_RUN_EXECUTION"]["argv"] = [
        app["CLI"], "opt", "-i", str(detached), "--out-dir", "result",
    ]
    app["_set_running"](True)
    app["S"]["inputs"] = []
    app["_delete_owned_uploads"]([str(detached)])
    assert detached.exists()
    assert str(detached) in app["_RUN_EXECUTION"]["deferred_deletes"]
    app["_set_running"](False)
    assert not detached.exists()

    reattached = tmp_path / "reattached.pdb"
    reattached.write_text("END\n", encoding="utf-8")
    app["_remember_uploaded_paths"]([str(reattached)])
    app["_RUN_EXECUTION"]["argv"] = [app["CLI"], "opt", "-i", str(reattached)]
    app["_set_running"](True)
    app["_delete_owned_uploads"]([str(reattached)])
    app["S"]["inputs"] = [str(reattached)]
    app["_set_running"](False)
    assert reattached.exists()
    assert str(reattached) in app["S"]["_uploaded_paths"]
    app["S"]["inputs"] = []
    app["_delete_owned_uploads"]([str(reattached)])
    assert not reattached.exists()


def test_colab_operates_scientific_selectors_and_remaining_buttons(
    tmp_path: Path, monkeypatch,
) -> None:
    app, _rendered = _execute_app(monkeypatch, tmp_path)
    first = tmp_path / "first.pdb"
    second = tmp_path / "second.pdb"
    parm = tmp_path / "system.parm7"
    pdb_text = (
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  LIG A   1       1.200   0.000   0.000  1.00  0.00           O\n"
        "HETATM    3  N1  LIG A   1       0.000   1.200   0.000  1.00  0.00           N\n"
        "HETATM    4  O   HOH A   2       3.000   0.000   0.000  1.00  0.00           O\n"
        "END\n"
    )
    first.write_text(pdb_text, encoding="utf-8")
    second.write_text(pdb_text.replace("3.000", "3.100"), encoding="utf-8")
    parm.write_text("parm", encoding="utf-8")
    app["load_pdb"]([str(first), str(second)], parm=str(parm))
    assert app["viewer_out"].__class__.__name__ == "HTML"
    assert 'class="rxmolstar-frame"' in app["viewer_out"].value
    notebook_text = NOTEBOOK.read_text(encoding="utf-8")
    assert "viewer_out = W.HTML(" in notebook_text
    assert "viewer_signal_out = W.HTML(" in notebook_text
    assert "with viewer_out:" not in notebook_text
    assert "output.clear_output()" not in notebook_text
    assert "with output: clear_output()" not in notebook_text
    assert "with res_out: clear_output()" not in notebook_text
    assert "output.outputs = ()" in notebook_text
    setup_page = app["_TAB_PAGES"][1][1]
    viewer_generation = app["_VIEWER_GENERATION"]["value"]
    for index, (_label, expected_page) in enumerate(app["_TAB_PAGES"]):
        app["_tab_go"](index)
        assert expected_page.layout.display == "flex"
        assert "rxpage-prewarm" not in expected_page._dom_classes
        for page_index, (_name, page) in enumerate(app["_TAB_PAGES"]):
            if page is expected_page:
                continue
            if page_index == 1:
                assert page.layout.display == "flex"
                assert "rxpage-prewarm" in page._dom_classes
            else:
                assert page.layout.display == "none"
        assert f"Step {index + 1} of 4" in app["_tab_status"].value
        assert app["_tab_choice"].value == index
        assert app["_tab_loading"].value == ""
    assert app["_VIEWER_GENERATION"]["value"] == viewer_generation
    assert len(app["_tab_buttons"]) == 4
    assert all(button.layout.flex == "1 1 0" for button in app["_tab_buttons"])
    app["_tab_buttons"][1].click()
    assert app["_tab_choice"].value == 1
    app["_tab_go"](1)
    assert "rxpage-prewarm" not in setup_page._dom_classes
    assert 'class="rxmolstar-frame"' in app["viewer_out"].value

    app["_tab_go"](2)
    assert app["_tab_choice"].value == 2
    assert app["_TAB_NAV"]["active"] == 2
    assert app["_tab_loading"].value == ""
    assert "rxpages-loading" not in app["_tab_body"]._dom_classes

    real_example_impl = app["_load_example_impl"]
    example_calls = []
    app["_load_example_impl"] = lambda _button: example_calls.append("loaded")
    assert app["_on_colab_example"]() == {"ok": True}
    assert example_calls == ["loaded"]
    assert not app["ex_btn"].disabled
    assert app["_tab_loading"].value == ""
    assert "rxpages-loading" not in app["_tab_body"]._dom_classes
    app["_load_example_impl"] = real_example_impl

    loading_seen = {}
    real_results = app["_results"]
    def _fake_results(root):
        loading_seen["value"] = app["_tab_loading"].value
        app["S"]["_results_presented_dir"] = os.path.abspath(root)
    app["_results"] = _fake_results
    app["S"]["_last_out_dir"] = str(tmp_path)
    app["S"]["_last_manifest"] = {"status": "success"}
    app["S"]["_results_presented_dir"] = None
    app["res_btn"].disabled = False
    app["_tab_go"](3)
    assert "Loading <b>④ Results</b>…" in loading_seen["value"]
    assert app["_tab_loading"].value == ""
    app["_results"] = real_results

    launch_disabled_before = {
        name: app[name].disabled for name in ("b_validate", "b_run")
    }
    app["_set_running"](True)
    assert app["b_validate"].disabled
    assert app["b_run"].disabled
    assert not app["dd_subcmd"].disabled
    assert not app["w_out"].disabled
    assert not app["upl"].disabled
    assert app["_view_is_editable"]()
    assert app["_on_colab_upload"]("input", [])["failed"] == []
    assert not _widget_with_description(app["input_file_rows"], "Remove file").disabled
    app["_set_running"](False)
    assert {
        name: app[name].disabled for name in ("b_validate", "b_run")
    } == launch_disabled_before
    assert not app["dd_subcmd"].disabled
    assert not app["w_out"].disabled
    assert not app["upl"].disabled
    assert not _widget_with_description(app["input_file_rows"], "Remove file").disabled

    app["dd_subcmd"].options = app["_sub_options"](app["SUBS"])
    app["dd_subcmd"].value = "energy-diagram"
    for control in (
        app["pick_action"], app["last_pick_info"], app["last_pick_row"],
    ):
        assert control.layout.display == "none"
    assert app["viewer_more"].layout.display == "none"
    assert app["viewer_more"] not in app["viewer_toolbar"].children
    app["dd_subcmd"].value = "scan"
    assert app["pick_action"].layout.display == ""
    assert app["last_pick_row"].layout.display == ""
    assert app["viewer_more"].layout.display == "none"

    _widget_with_description(app["input_file_rows"], "Move later").click()
    assert app["S"]["inputs"] == [str(second), str(first)]
    assert "Mol* update bridge" in app["viewer_signal_out"].value
    assert "rx-load-structure" in app["viewer_signal_out"].value
    _widget_with_description(app["input_file_rows"], "Move earlier").click()
    assert app["S"]["inputs"] == [str(first), str(second)]
    app["load_pdb"]([str(first)], parm=str(parm))

    def pick(index: int) -> None:
        app["on_click"](str(index), live_marked=True)

    def scan_pair(first_index: int, second_index: int) -> None:
        _widget_with_description(app["scan_panel"], "① Pick atom A").click()
        assert app["pick_action"].value == "scanA"
        pick(first_index)
        _widget_with_description(app["scan_panel"], "② Pick atom B").click()
        assert app["pick_action"].value == "scanB"
        pick(second_index)

    app["dd_subcmd"].value = "scan"
    scan_pair(0, 1)
    target = _widget_with_description(app["scan_panel"], "③ Set target Å")
    assert target.step == pytest.approx(0.5)
    target.value = 1.8
    _widget_with_description(app["scan_panel"], "Add sequential stage").click()
    assert len(app["S"]["scan_stages"]) == 1
    scan_pair(0, 1)
    _widget_with_description(app["scan_panel"], "Add to stage").click()
    assert len(app["S"]["scan_stages"][0]) == 1
    assert "already in the current stage" in app["S"]["_last_pick_message"]
    scan_pair(1, 0)
    _widget_with_description(app["scan_panel"], "Add to stage").click()
    assert len(app["S"]["scan_stages"][0]) == 1
    scan_pair(0, 2)
    _widget_with_description(app["scan_panel"], "Add to stage").click()
    assert len(app["S"]["scan_stages"][0]) == 2
    staged_controls = list(_widget_descendants(app["scan_panel"]))
    assert not any(str(getattr(widget, "tooltip", "") or "").startswith("Move coordinate")
                   for widget in staged_controls)
    assert not any(getattr(widget, "description", "") == "Use picked pair"
                   for widget in staged_controls)
    coordinate_removals = [
        widget for widget in staged_controls
        if getattr(widget, "tooltip", "") == "Remove coordinate"
    ]
    assert len(coordinate_removals) == 2
    coordinate_removals[-1].click()
    assert len(app["S"]["scan_stages"][0]) == 1
    assert not any(getattr(widget, "tooltip", "") == "Remove coordinate"
                   for widget in _widget_descendants(app["scan_panel"]))
    next(
        widget for widget in _widget_descendants(app["scan_panel"])
        if getattr(widget, "tooltip", "") == "Remove stage"
    ).click()
    assert app["S"]["scan_stages"] == []
    next(
        widget for widget in _widget_descendants(app["scan_panel"])
        if getattr(widget, "tooltip", "") == "Clear the atom pair being prepared."
    ).click()
    assert app["S"]["scan_atoms"] == [None, None]

    for subcommand, pairs in (
        ("scan2d", ((0, 1), (0, 2))),
        ("scan3d", ((0, 1), (0, 2), (1, 2))),
    ):
        app["dd_subcmd"].value = subcommand
        for pair_index, (left, right) in enumerate(pairs):
            scan_pair(left, right)
            low = _widget_with_description(app["scan_panel"], "③ Set low Å")
            high = _widget_with_description(app["scan_panel"], "High Å")
            low.value = 1.0 + pair_index * 0.1
            high.value = 2.0 + pair_index * 0.1
            _widget_with_description(app["scan_panel"], "4  Add axis").click()
        assert len(app["S"]["scan_axes"]) == len(pairs)
        _widget_with_description(app["scan_panel"], "×").click()
        assert len(app["S"]["scan_axes"]) == len(pairs) - 1
        app["b_clear_scan"].click()
        assert app["S"]["scan_axes"] == []

    app["dd_subcmd"].value = "opt"
    app["pick_action"].value = "freezeA"
    pick(0)
    app["pick_action"].value = "freezeB"
    pick(1)
    _widget_with_description(app["freeze_panel"], "add freeze pair").click()
    assert len(app["S"]["freeze_pairs"]) == 1
    _widget_with_description(app["freeze_panel"], "clear pairs").click()
    assert app["S"]["freeze_pairs"] == []
    app["pick_action"].value = "freezeatom"
    pick(2)
    assert app["S"]["freeze_atoms"] == [3]
    assert app["chips_box"].layout.display == "flex"
    assert app["chips_box"].layout.min_height == "38px"
    _widget_with_description(app["chips_box"], "⚓ 3 ✕").click()
    assert app["S"]["freeze_atoms"] == []
    pick(2)
    assert app["S"]["freeze_atoms"] == [3]
    _widget_with_description(app["freeze_panel"], "Clear").click()
    assert app["S"]["freeze_atoms"] == []

    app["pick_action"].value = "center"
    app["exact_atom"].value = "1"
    app["exact_atom_btn"].click()
    assert app["S"]["_last_pick"]["index"] == 0
    exact_center_chip = _widget_with_description(app["chips_box"], "A:LIG:1 ✕")
    assert exact_center_chip.button_style == "info"
    app["pick_action"].value = "selectedresn"
    pick(0)
    selected_chip = _widget_with_description(app["chips_box"], "🔒 1 ✕")
    assert "force-included --selected-resn" in app["summary_html"].value
    selected_chip.click()
    assert app["selected_resn"].value == ""
    app["center_widget"].value = ("LIG",)
    app["S"]["center_ids"] = ["A:LIG:1"]
    app["selected_resn"].value = "A:123"
    app["S"]["freeze_atoms"] = [3]
    app["_render_chips"]()
    app["_render_summary"]()
    assert "center · <code>LIG,A:LIG:1</code>" in app["summary_html"].value
    assert "exact <code>" not in app["summary_html"].value
    assert "force-included --selected-resn" in app["summary_html"].value
    app["dd_subcmd"].value = "all"
    app["all_mode"].value = "mep"
    expected_chips = [widget.description for widget in app["chips_box"].children]
    app["all_mode"].value = "scan"
    assert app["chips_box"].layout.display == "flex"
    assert [widget.description for widget in app["chips_box"].children] == expected_chips
    for description in ("LIG ✕", "A:LIG:1 ✕", "🔒 A:123 ✕", "⚓ 3 ✕"):
        _widget_with_description(app["chips_box"], description).click()
    assert app["center_widget"].value == ()
    assert app["selected_resn"].value == ""
    assert app["S"]["center_ids"] == [] and app["S"]["freeze_atoms"] == []
    app["S"]["center_ids"] = ["A:LIG:1"]
    app["S"]["freeze_atoms"] = [2]
    app["b_clear"].click()
    assert app["S"]["center_ids"] == []
    assert app["S"]["freeze_atoms"] == [2]

    app["dd_subcmd"].value = "all"
    app["adv_search"].value = "cycles"
    assert "total" in app["adv_count"].value
    app["adv_search"].value = ""
    folds = [
        widget for widget in _widget_descendants(app["rootbox"])
        if "rxfold" in getattr(widget, "_dom_classes", ())
    ]
    assert folds
    for fold in folds:
        fold._rx_set_open(False)
        fold._rx_button.click()
        assert fold._rx_body.layout.display == ""
        fold._rx_button.click()
        assert fold._rx_body.layout.display == "none"

    clipboard: list[str] = []
    downloads: list[str] = []
    google = types.ModuleType("google")
    colab = types.ModuleType("google.colab")
    output = types.ModuleType("google.colab.output")
    files = types.ModuleType("google.colab.files")
    output.eval_js = lambda script: clipboard.append(str(script))
    files.download = lambda path: downloads.append(str(path))
    colab.output = output
    colab.files = files
    google.colab = colab
    monkeypatch.setitem(sys.modules, "google", google)
    monkeypatch.setitem(sys.modules, "google.colab", colab)
    monkeypatch.setitem(sys.modules, "google.colab.output", output)
    monkeypatch.setitem(sys.modules, "google.colab.files", files)

    app["dd_subcmd"].value = "fix-altloc"
    app["b_rebuild"].click()
    app["b_copy"].click()
    assert clipboard and json.dumps(app["cmd_box"].value) in clipboard[-1]
    saved_out = app["w_out"].value
    app["b_save"].click()
    session_path = Path(downloads[-1])
    session_bytes = session_path.read_bytes()
    app["w_out"].value = "changed-after-save"
    _set_file_upload(app["up_sess"], "session.json", "application/json", session_bytes)
    assert app["w_out"].value == saved_out
    assert not app["up_sess"].value

    trajectory_a = tmp_path / "path_a_trj.xyz"
    trajectory_b = tmp_path / "path_b_trj.xyz"
    trajectory_text = (
        "2\nenergy=-1.000000 unit=hartree\nH 0 0 0\nH 0.7 0 0\n"
        "2\nenergy=-0.990000 unit=hartree\nH 0.2 0 0\nH 0.9 0 0\n"
    )
    trajectory_a.write_text(trajectory_text, encoding="utf-8")
    trajectory_b.write_text(
        trajectory_text.replace("-0.990000", "-0.980000"), encoding="utf-8",
    )
    artifact_a = tmp_path / "artifact_a.pdb"
    artifact_b = tmp_path / "artifact_b.pdb"
    artifact_a.write_text(pdb_text, encoding="utf-8")
    artifact_b.write_text(pdb_text.replace("1.200", "1.300"), encoding="utf-8")
    current_files = [trajectory_a, trajectory_b, artifact_a, artifact_b]
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps({
        "mlmm_toolkit_version": "0.3.3",
        "status": "success",
        "command": "mlmm path-search -i input.pdb",
        "current_output_paths": [path.name for path in current_files],
    }), encoding="utf-8")
    app["results_dir"].value = str(tmp_path)
    app["res_btn"].click()
    assert app["res_btn"].description == "Load results"
    assert app["S"]["_last_out_dir"] == str(tmp_path.resolve())
    app["S"]["_last_log"] = "button coverage log"
    # The selector names the scientifically selected result rather than using
    # a generic "Primary result" label.
    assert [label for label, _ in app["traj_choice"].options] == [
        "Trajectory · path_a_trj.xyz"
    ]
    app["traj_choice"].value = app["traj_choice"].options[0][1]
    app["frame_next"].click()
    assert app["frame_slider"].value == 1
    app["frame_prev"].click()
    assert app["frame_slider"].value == 0
    assert app["frame_play"].disabled is False
    assert app["frame_prev"].description == "‹"
    assert app["frame_next"].description == "›"
    assert app["frame_controls"].children == (
        app["frame_prev"], app["frame_next"], app["frame_play"], app["frame_slider"],
        app["playback_interval_signal"],
    )
    notebook_text = NOTEBOOK.read_text(encoding="utf-8")
    assert "grid-template-columns:40px 40px 116px minmax(0,1fr)" in notebook_text
    assert "sliderObserver.observe(slider,{attributes:true,attributeFilter:['aria-valuenow']})" in notebook_text
    assert "queueFrame(index); broadcast(index); invoke(index);" in notebook_text
    app["artifact_choice"].value = str(artifact_b)
    app["artifact_fold"]._rx_set_open(False)
    app["artifact_fold"]._rx_button.click()
    assert app["artifact_fold"]._rx_body.layout.display == ""
    complete_log = "begin\n" + ("x" * 250100) + "\nend\n"
    app["_RUN_LOG_BUFFER"]["text"] = complete_log
    assert app["_visible_run_log_text"]() == complete_log
    persisted = Path(app["_persist_run_log"](tmp_path / "logged", complete_log))
    assert persisted == tmp_path / "logged" / "run.log"
    assert persisted.read_text(encoding="utf-8") == complete_log
    app["dl_btn"].click()
    archive_path = Path(downloads[-1])
    with zipfile.ZipFile(archive_path) as archive:
        members = set(archive.namelist())
        assert archive.read("run.log").decode() == "button coverage log"
    assert {"colab_run.json", "run.log"} <= members
    assert {path.name for path in current_files} <= members

    app["_example_file"] = lambda relpath: str(NOTEBOOK.parent / relpath)
    for choice in app["ex_choice"].options:
        app["ex_choice"].value = choice
        app["ex_btn"].click()
        assert app["S"]["inputs"], choice
        assert app["S"]["parm"] and Path(app["S"]["parm"]).is_file()
        assert "⚠️" not in app["example_msg"].value
        assert app["S"]["subcmd"] == "all"
        if choice.startswith("Toy system"):
            assert app["all_mode"].value == "mep"
            assert app["S"]["tsopt"] is True
            assert app["S"]["thermo"] is True
            command = app["build_cmd"]()
            assert command[:2] == ["mlmm", "all"]
            assert command.count("--tsopt") == 1
            assert command.count("--thermo") == 1


def test_colab_stop_child_reaps_after_disappeared_process_group(
    tmp_path: Path, monkeypatch,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)

    class FakeProcess:
        pid = 4242
        returncode = None

        def __init__(self):
            self.waited = False

        def poll(self):
            return self.returncode

        def wait(self, timeout=None):
            self.waited = True
            self.returncode = -15
            return self.returncode

    process = FakeProcess()
    monkeypatch.setattr(app["os"], "getpgid", lambda pid: pid)
    monkeypatch.setattr(
        app["os"], "killpg",
        lambda *_args: (_ for _ in ()).throw(ProcessLookupError()),
    )

    app["_stop_child"](process)

    assert process.waited


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
    _set_file_upload(app["model_upl"], "model.pdb", "chemical/x-pdb", payload)
    assert app["S"]["inputs"] == [str(full)] and app["S"]["parm"] == str(parm)
    assert app["S"]["model_pdb"] and Path(app["S"]["model_pdb"]).name == "model.pdb"
    assert not app["model_upl"].value
    first_model = Path(app["S"]["model_pdb"])
    assert first_model.exists()
    app["model_clear"].click()
    assert app["S"]["model_pdb"] is None
    assert not first_model.exists()
    _set_file_upload(app["model_upl"], "model.pdb", "chemical/x-pdb", payload)
    assert app["S"]["model_pdb"] and Path(app["S"]["model_pdb"]).name.startswith("model")
    assert not app["model_upl"].value
    second_model = Path(app["S"]["model_pdb"])
    relative_alias = os.path.relpath(second_model, tmp_path)
    app["_clear_structure_bound_state"](preserve_model=relative_alias)
    assert second_model.exists()
    assert str(second_model.resolve()) in app["S"]["_uploaded_paths"]
    app["_delete_owned_uploads"]([relative_alias])
    assert not second_model.exists()


def test_colab_topology_replacement_removes_only_detached_owned_copy(
    tmp_path: Path, monkeypatch,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)
    first = app["_save_upload"]("first.parm7", memoryview(b"first topology\n"))
    app["_remember_uploaded_paths"]([first])
    app["S"]["parm"] = first
    assert Path(first).exists()

    assert app["_accept_upload_pairs"](
        [("replacement.parm7", memoryview(b"replacement topology\n"))],
        "test replacement",
    )

    replacement = app["S"]["parm"]
    assert replacement and Path(replacement).exists()
    assert not Path(first).exists()
    assert os.path.abspath(replacement) in app["S"]["_uploaded_paths"]


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

    app["_render_advanced_rows"]()
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
    generation = app["_VIEWER_GENERATION"]["value"]
    signal = app["viewer_signal_out"].value
    app["render_viewer"]()
    update = html.unescape(app["viewer_signal_out"].value)
    assert "rx-load-structure" in update
    assert str(app["_VIEWER_GENERATION"]["value"]) in update
    assert app["_VIEWER_GENERATION"]["value"] == generation
    assert app["viewer_signal_out"].value == signal
    app["center_widget"].value = ()
    assert app["_VIEWER_GENERATION"]["value"] == generation
    assert app["viewer_signal_out"].value == signal
    app["center_widget"].value = ("LIG",)
    assert app["_VIEWER_GENERATION"]["value"] == generation
    assert app["viewer_signal_out"].value == signal

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
    wrong_tool["tool"] = "other-tool"
    with pytest.raises(ValueError, match="belong"):
        app["_apply_session"](wrong_tool)

    saved = app["_session_dict"]()
    saved.update(
        inputs=[str(primary), str(product)], parm=str(topology), mode="pdb",
        subcmd="all", all_mode="mep", tsopt=False, thermo=False,
        backend="uma", model="uma-m-1p1",
        rep="sticks", rep_user_set=True, color="spectrum",
    )
    saved["file_identities"]["inputs"] = [
        app["_file_identity"](str(primary)),
        app["_file_identity"](str(product)),
    ]
    saved["file_identities"]["parm"] = app["_file_identity"](str(topology))
    saved["primary_atom_signatures"] = [
        list(item) for item in app["_atom_signatures"](metadata[str(primary)])
    ]
    saved["advanced_overrides"] = {"all": {"verbose": 3}}
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
    assert app["dd_rep"].value == "stick" and app["dd_col"].value == "spectrum"
    assert app["_rep_user_set"]["value"] is True
    assert app["S"]["advanced_overrides"]["all"]["verbose"] == 3
    assert app["_auto"]["on"] is True and app["cmd_box"].value != stale_command
    bad_verbose = app["_session_dict"]()
    bad_verbose["advanced_overrides"] = {"all": {"verbose": "3"}}
    before_bad_verbose = json.dumps(app["_session_dict"](), sort_keys=True)
    with pytest.raises(ValueError, match="integer from 0 to 3"):
        app["_apply_session"](bad_verbose)
    assert json.dumps(app["_session_dict"](), sort_keys=True) == before_bad_verbose
    untrusted_confirmation = app["_session_dict"]()
    untrusted_confirmation["charge_explicit"] = True
    assert app["_apply_session"](untrusted_confirmation) == []
    assert app["w_charge_ok"].value is False
    assert app["S"]["charge_explicit"] is False
    assert app["S"]["_charge_scope"] is None

    dependent = app["_session_dict"]()
    dependent["tsopt"] = False
    dependent["thermo"] = True
    dependent["advanced"]["dft"] = False
    assert app["_apply_session"](dependent) == []
    assert app["w_ts"].value is True and app["w_ts"].disabled is False
    assert app["_session_dict"]()["tsopt"] is True
    app["w_ts"].value = False
    assert app["w_th"].value is False

    app["all_mode"].value = "mep"
    app["w_ts"].value = False
    app["all_mode"].value = "tsonly"
    ts_only_saved = app["_session_dict"]()
    assert ts_only_saved["tsopt"] is False
    assert app["_apply_session"](ts_only_saved) == []
    app["all_mode"].value = "mep"
    assert app["w_ts"].value is False

    legacy = app["_session_dict"]()
    legacy["schema_version"] = 1
    before_legacy = json.dumps(app["_session_dict"](), sort_keys=True)
    with pytest.raises(ValueError, match="Save settings again"):
        app["_apply_session"](legacy)
    assert json.dumps(app["_session_dict"](), sort_keys=True) == before_legacy

    tampered_signature = app["_session_dict"]()
    tampered_signature["primary_atom_signatures"][0][4] = "WRONG"
    before_tamper = json.dumps(app["_session_dict"](), sort_keys=True)
    with pytest.raises(ValueError, match="atom identity/order"):
        app["_apply_session"](tampered_signature)
    assert json.dumps(app["_session_dict"](), sort_keys=True) == before_tamper

    # The reuse opt-in is gone: an occupied output can never be silently
    # overwritten, so no per-run opt-in survives into the next run.
    app["_invalidate_last_run"]("new system")
    assert "w_reuse" not in app
    occupied = tmp_path / "occupied-out"
    occupied.mkdir()
    (occupied / "leftover.txt").write_text("x", encoding="utf-8")
    occupied_argv = ["mlmm", "opt", "-i", str(primary), "-o", str(occupied)]
    assert app["_output_scope_collision"](app["_output_scope"](occupied_argv))
    diverted = app["_preflight_output_scope"](occupied_argv)
    assert diverted is not None
    assert Path(diverted["root"]).resolve() != occupied.resolve()
    assert not app["_output_scope_collision"](diverted)

    missing = tmp_path / "missing-dir" / "missing.pdb"
    pending = app["_session_dict"]()
    pending.update(
        backend="mace", model="MACE-OMOL-0", inputs=[str(missing)], parm=str(topology),
        mode="pdb", subcmd="opt", all_mode="mep", center=[], center_ids=[], lcharge={},
    )
    pending["file_identities"]["inputs"] = [
        app["_file_identity"](str(primary)),
    ]
    pending["file_identities"]["parm"] = app["_file_identity"](str(topology))
    pending["primary_atom_signatures"] = [
        list(item) for item in app["_atom_signatures"](metadata[str(primary)])
    ]
    assert app["_apply_session"](pending) == [str(missing)]
    uploaded = tmp_path / "missing.pdb"
    uploaded.write_text(incompatible.read_text(encoding="utf-8"), encoding="utf-8")
    metadata[str(uploaded)] = [dict(row) for row in metadata[str(incompatible)]]
    with pytest.raises(ValueError, match="SHA-256 identity"):
        app["_ingest_saved_files"]([str(uploaded)], "wrong session re-upload")
    assert app["S"]["inputs"] == [str(missing)]
    uploaded.write_text(primary.read_text(encoding="utf-8"), encoding="utf-8")
    metadata[str(uploaded)] = [dict(row) for row in metadata[str(primary)]]
    assert app["_ingest_saved_files"]([str(uploaded)], "session replacement")
    assert app["S"]["inputs"] == [str(uploaded)] and app["S"]["parm"] == str(topology)

    missing_topology = tmp_path / "missing-dir" / "system.parm7"
    pending_topology = app["_session_dict"]()
    pending_topology["parm"] = str(missing_topology)
    pending_topology["file_identities"]["parm"] = app["_file_identity"](
        str(topology)
    )
    assert app["_apply_session"](pending_topology) == [str(missing_topology)]
    wrong_topology = tmp_path / "wrong" / "system.parm7"
    wrong_topology.parent.mkdir()
    wrong_topology.write_text("different topology\n", encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256 identity"):
        app["_ingest_saved_files"]([str(wrong_topology)], "wrong topology")
    assert app["S"]["parm"] == str(missing_topology)
    assert app["_ingest_saved_files"]([str(topology)], "matching topology")
    assert app["S"]["parm"] == str(topology)

    bad = tmp_path / "empty.pdb"
    bad.write_text("REMARK no atoms\nEND\n", encoding="utf-8")
    old_inputs = list(app["S"]["inputs"])
    app["_load_view_structure"] = real_loader
    assert not app["_ingest_saved_files"]([str(bad)], "invalid")
    assert app["S"]["inputs"] == old_inputs and "not attached" in app["input_msg"].value

    duplicate = {"index": 0, "chain": "A", "resn": "LIG", "resi": "10", "atom": "C1", "xyz": (0, 0, 0)}
    app["dd_subcmd"].value = "scan"
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
    assert "def _advanced_options(sub):" in app
    assert "if 'freeze' in SPEC.get(sub, {}).get('panels', ()) and S['freeze_atoms']:" in app
    # The public TR-projection option and legacy treatment are removed.
    assert "legacy-active" not in app
    assert "--tr-projection" not in app
    # ONIOM utilities are classified, with their inputs/outputs stated.
    for _sub in ("'define-layer'", "'mm-parm'", "'oniom-export'", "'oniom-import'"):
        assert _sub in app
    # Surfaced key flags + the standalone cluster-model button.
    assert "key_opts_box = _collapsible('Key options', key_opts_content)" in app
    assert "key_opts_content = W.VBox([\n    _flag_row(adv_mep" in app
    assert "W.HBox([_flag_row(adv_mep" not in app
    assert "adv_prec = W.Dropdown" in app
    assert "width='420px', max_width='92%', min_width='0'" in app
    assert "cmd += ['--flatten']" in app
    assert "_ADVANCED_DEFAULT" in app
    assert "b_extract = W.Button(description='Extract ML region & use it'" in app
    assert "prep_radius = W.BoundedFloatText(" in app
    # The wheel ships no examples, so Load example resolves them from the git
    # tag matching the installed release (a source checkout is used when present).
    assert "def _example_file(relpath):" in app
    assert "raw.githubusercontent.com/t-0hmura/mlmm_toolkit" in app
    assert "Run Setup first" not in app
    # -r/--radius is extraction-only and requires a structure workflow.
    assert "if radius_applies: cmd += ['-r', str(r)]" in app
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
    assert "<b>Drop an ML-region PDB here</b> or choose it below" in app
    assert "Used only as <code>--model-pdb</code>; full-system inputs remain unchanged." in app
    assert "model_upload_box.add_class('rxmodel-upload')" in app
    assert "_model_actions.add_class('rxmodel-actions')" in app
    assert "model_upload_msg = W.HTML()" in app
    assert "model_clear = W.Button(description='Clear'" in app
    assert "var actionHost = spec.role === 'model' ? zone.querySelector('.rxmodel-actions') : null;" in app
    assert "if(actionHost) actionHost.insertBefore(button, actionHost.firstChild);" in app
    assert ".rxapp .rxmodel-upload {" in app and "min-height:116px" in app
    assert "margin-top:8px !important" in app
    assert ".rxapp .rxmodel-state:has(.widget-html-content:empty) { display:none !important; }" in app
    assert "if(files.length!==1)" in app
    assert "name.endsWith('.pdb')||name.endsWith('.ent')" in app
    assert "input.dispatchEvent(new Event('change',{bubbles:true}))" in app
    assert "_accept_model_pairs(pairs); accepted = True" in app
    assert "selector='.rxmodel-upload', role='model'" in app
    assert "'oniom-export': dict(n_in=(0, 1)" in app
    assert "def _clear_structure_bound_state(preserve_model=None):" in app
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
    assert "class _DropUpload(anywidget.AnyWidget):" in app
    assert "upl = _DropUpload(accept=_acc, formats=_drop_formats, formats_detail=_drop_formats_detail)" in app
    assert "mlmm_gui.on_drop" not in app
    assert "def _ingest_saved_files(" in app
    assert "_reset_file_upload(upl)" in app
    assert "if _UPLOAD_MODE == 'anywidget': upl.on_msg(_on_drop_upload)" in app
    assert "pts.append((k, 'TS'))" not in app
    assert "points[k] = 'TS candidate'" not in app
    assert "points[k] = 'minimum candidate'" not in app
    assert "'extrema': True" not in app
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
    assert "submit(event.dataTransfer ? event.dataTransfer.files : []);" in app
    assert "_tab_body.children = [_TAB_PAGES[i][1]]" not in app
    assert "layout=W.Layout(width='560px')" not in app
    assert "system charge (-q)" in app and "overwrite system charge" in app
    assert ("No <code>-l</code> charge source is active, so <code>-q</code> "
            "is used directly.") in app
    assert "No ligand-charge source is available" not in app
    assert ".pdb / .cif / .xyz / .gjf + optional .parm7" in app
    assert "utility .gjf / .com / .inp / .csv" not in app
    assert app.count("effective = _normalized_scope_argv(a)") == 2
    assert "dry_argv = _force_dry_run(list(a))" in app
    assert "if not _validate_command(a): return" in app
    assert "real_run = not _flag_enabled(effective, '--dry-run', '--no-dry-run')" in app
    assert "def _utility_autofill_complete(argv):" in app
    assert (
        app.index("input_records = (")
        < app.index("rc, transcript = _stream(a)")
        < app.index("'inputs': input_records")
    )


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
        output_params = [
            param for param in command.params
            if {"-o", "--out", "--out-dir", "--output"}.intersection(param.opts)
        ]
        assert len(output_params) == 1, (
            subcommand,
            [(param.name, tuple(param.opts)) for param in command.params],
        )
        output_param = output_params[0]
        assert scope_for(["mlmm", subcommand])["root"] == str(output_param.default)
    assert scope_for(["mlmm", "opt", "-oattached"])["root"] == "attached"
    assert scope_for(["mlmm", "opt", "--out-dir", "Upper"])["root"] == "Upper"
    assert scope_for(["mlmm", "opt", "--out-dir=UpperEq"])["root"] == "UpperEq"
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


def test_colab_release_state_and_linked_results_regressions(
    tmp_path: Path, monkeypatch,
) -> None:
    app, rendered = _execute_app(monkeypatch, tmp_path)

    notebook_source = _notebook()["cells"][2]["source"]
    assert "rxpath-panel-head" in notebook_source
    assert "xanchor:item.index===0?'left'" in notebook_source
    assert "shapes.push({type:'line',xref:'x',yref:'paper',x0:1" not in notebook_source
    assert "the profile, controls, and molecular structure stay synchronized" not in notebook_source
    assert "never" in [value for _, value in app["adv_thresh"].options]
    assert app["adv_dftfb"].placeholder == "wb97m-v/def2-tzvpd"
    app["S"]["mode"] = "xyz"
    assert app["_aspec"]({
        "index": 3, "chain": "", "resn": "MOL", "resi": "1",
        "icode": "", "atom": "C4",
    }) == 4
    app["S"]["mode"] = "pdb"
    assert app["_aspec"]({
        "index": 4, "chain": "", "resn": "LIG", "resi": "12A",
        "icode": "A", "atom": "C1",
    }) == 5

    primary = tmp_path / "primary.pdb"
    product = tmp_path / "product.pdb"
    pdb_text = (
        "HETATM    1  C1  LIG A  10       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  LIG A  10       1.200   0.000   0.000  1.00  0.00           O\n"
        "END\n"
    )
    primary.write_text(pdb_text, encoding="utf-8")
    product.write_text(pdb_text, encoding="utf-8")
    topology = tmp_path / "system.parm7"
    topology.write_text("%VERSION VERSION_STAMP = V0001.000\n", encoding="utf-8")

    app["load_pdb"](
        [str(primary), str(product)], parm=str(topology), center=["LIG"],
    )
    app["charge_rows"]["LIG"]["use"].value = True
    app["charge_rows"]["LIG"]["val"].value = 1
    app["w_charge_ok"].value = True
    app["view_input"].value = 1
    app["view_input"].value = 0
    assert app["w_charge_ok"].value is True
    assert app["S"]["lcharge"] == {"LIG": 1.0}
    app["S"]["scan_preset"] = "[('OLD 1 A','OLD 1 B',1.4)]"
    app["S"]["scan_atoms"] = [None, None]
    app["dd_subcmd"].value = "scan"
    app["_render_scan_panel"]()
    app["_render_summary"]()
    scan_panel_text = "\n".join(
        str(getattr(widget, "value", ""))
        for widget in _widget_descendants(app["scan_panel"])
    )
    assert "<b>A</b>: <code>OLD 1 A</code>" in scan_panel_text
    assert "<b>B</b>: <code>OLD 1 B</code>" in scan_panel_text
    assert "Loaded example · Clear to redefine." in scan_panel_text
    assert "OLD 1 A — OLD 1 B → 1.4 Å" in app["summary_html"].value
    assert "[(&#x27;OLD 1 A&#x27;" not in app["summary_html"].value
    app["pick_action"].value = "scanA"
    app["on_click"]("0", "LIG", "10", "A", "C1", "1", "", live_marked=True)
    assert app["S"]["scan_preset"] == ""
    assert app["scan_literals"]() == []

    center_file = tmp_path / "centers.pdb"
    center_file.write_text("MODEL        1\nENDMDL\n", encoding="utf-8")
    center_argv = [
        "mlmm", "all", "-i", str(primary), "--parm", str(topology),
        "-c", str(center_file),
    ]
    assert str(center_file.resolve()) in app["_command_input_files"](center_argv)
    first_center_fingerprint = app["_validation_fingerprint"](center_argv)
    center_file.write_text("MODEL        2\nENDMDL\n", encoding="utf-8")
    assert app["_validation_fingerprint"](center_argv) != first_center_fingerprint

    claimed_root = tmp_path / "claimed"
    claimed_root.mkdir()
    mode_path = claimed_root / "mode_0001_-100.00cm-1_trj.xyz"
    mode_path.write_text(
        "2\nenergy=-1 unit=hartree\nH 0 0 0\nH 0.7 0 0\n",
        encoding="utf-8",
    )
    result_path = claimed_root / "result.json"
    result_path.write_text(json.dumps({
        "files": {"frequencies_txt": "frequencies_cm-1.txt",
                  "mode_files": [mode_path.name]},
    }), encoding="utf-8")
    assert mode_path.resolve() in {
        Path(path).resolve() for path in app["_structured_current_paths"](
            str(claimed_root), [str(result_path), str(mode_path)],
        )
    }
    old_mode = claimed_root / "old_mode_trj.xyz"
    old_mode.write_text(mode_path.read_text(encoding="utf-8"), encoding="utf-8")
    new_mode = claimed_root / "new_mode_trj.xyz"
    new_mode.write_text(
        mode_path.read_text(encoding="utf-8").replace("energy=-1", "energy=-0.9"),
        encoding="utf-8",
    )
    result_path.write_text(json.dumps({
        "files": {"mode_files": [old_mode.name, new_mode.name]},
    }), encoding="utf-8")
    reused_paths = {
        Path(path).resolve() for path in app["_structured_current_paths"](
            str(claimed_root), [str(result_path), str(new_mode)],
        )
    }
    assert new_mode.resolve() in reused_paths
    assert old_mode.resolve() not in reused_paths

    aggregate_root = tmp_path / "aggregate"
    aggregate_root.mkdir()
    aggregate_summary = aggregate_root / "summary.json"
    aggregate_summary.write_text(json.dumps({
        "current_output_paths": ["summary.json"],
    }), encoding="utf-8")
    aggregate_vib = aggregate_root / "segments" / "seg_02" / "ts" / "vib"
    aggregate_vib.mkdir(parents=True)
    aggregate_imag = aggregate_vib / "imag_01_-250.00cm-1_trj.xyz"
    aggregate_imag.write_text(mode_path.read_text(encoding="utf-8"), encoding="utf-8")
    aggregate_noise = aggregate_vib / "unrelated.xyz"
    aggregate_noise.write_text(mode_path.read_text(encoding="utf-8"), encoding="utf-8")
    aggregate_mep_dir = aggregate_root / "_work" / "path_search"
    aggregate_mep_dir.mkdir(parents=True)
    aggregate_mep = aggregate_mep_dir / "mep_seg_02_trj.xyz"
    aggregate_mep.write_text(mode_path.read_text(encoding="utf-8"), encoding="utf-8")
    aggregate_paths = {
        Path(path).resolve() for path in app["_structured_current_paths"](
            str(aggregate_root),
            [str(aggregate_summary), str(aggregate_imag), str(aggregate_noise),
             str(aggregate_mep)],
        )
    }
    assert aggregate_imag.resolve() in aggregate_paths
    assert aggregate_mep.resolve() in aggregate_paths
    assert aggregate_noise.resolve() not in aggregate_paths

    nested_ts = claimed_root / "segments" / "seg_01" / "ts"
    nested_vib = nested_ts / "vib"
    nested_vib.mkdir(parents=True)
    nested_final = nested_ts / "final_geometry.xyz"
    nested_final.write_text(mode_path.read_text(encoding="utf-8"), encoding="utf-8")
    nested_mode = nested_vib / "imag_-421.25cm-1_trj.xyz"
    nested_mode.write_text(mode_path.read_text(encoding="utf-8"), encoding="utf-8")
    nested_result = nested_ts / "result.json"
    nested_result.write_text(json.dumps({
        "files": {"final_geometry": "final_geometry.xyz",
                  "mode_files": ["vib/imag_-421.25cm-1_trj.xyz"]},
    }), encoding="utf-8")
    nested_claims = {
        Path(path).resolve() for path in app["_structured_current_paths"](
            str(claimed_root),
            [str(nested_result), str(nested_final), str(nested_mode)],
        )
    }
    assert nested_final.resolve() in nested_claims
    assert nested_mode.resolve() in nested_claims
    assert '"mode_files": _mode_output_files' in (
        NOTEBOOK.parents[1] / "mlmm" / "workflows" / "freq.py"
    ).read_text(encoding="utf-8")

    good = tmp_path / "good_trj.xyz"
    good.write_text(
        "2\nenergy=-1.000000 unit=hartree\nH 0 0 0\nH 0.7 0 0\n"
        "2\nenergy=-0.990000 unit=hartree\nH 0.2 0 0\nH 0.9 0 0\n",
        encoding="utf-8",
    )
    bad = tmp_path / "bad_trj.xyz"
    bad.write_text("not an xyz trajectory\n", encoding="utf-8")
    app["_result_pick_guard"]["active"] = True
    try:
        app["traj_choice"].options = [
            ("bad", str(bad)), ("good", str(good)),
        ]
        app["traj_choice"].value = str(bad)
    finally:
        app["_result_pick_guard"]["active"] = False
    app["_load_trajectory"](str(bad), str(tmp_path))
    assert app["trajectory_box"].layout.display == ""
    assert app["trajectory_content"].layout.display == "none"

    rendered.clear()
    app["plot_out"].value = '<div class="rxenergy-frame">stale energy view</div>'
    app["_load_trajectory"](str(good), str(tmp_path))
    energy_view = app["plot_out"].value
    assert energy_view.count('class="rxenergy-frame"') == 1
    assert "stale energy view" not in energy_view
    assert "plotly_click" in html.unescape(energy_view)
    generation = app["_TRAJ"]["generation"]
    app["_set_frame_from_browser"](generation, 1)
    assert app["frame_slider"].value == 1
    app["_set_frame_from_browser"](generation - 1, 0)
    assert app["frame_slider"].value == 1

    semantics = {
        "start": "R", "end": "P", "extrema": True,
        "bridge_ranges": [(2, 3)],
    }
    assert app["_stationary"]([0.0, 4.0, 10.0, 4.0, 0.0], semantics) == [
        (0, "R"), (4, "P"),
    ]

    leaf = tmp_path / "leaf-summary.json"
    leaf.write_text(json.dumps({
        "status": "completed", "scientific_status": "success",
        "n_modes": 6, "n_imaginary": 1,
    }), encoding="utf-8")
    leaf_html = app["_summary_html"](str(leaf))
    assert "modes: <b>6</b>" in leaf_html
    assert "images:" not in leaf_html and "None" not in leaf_html

    app["trajectory_box"].layout.display = ""
    app["artifact_fold"].layout.display = ""
    app["_begin_results_attempt"](str(tmp_path / "new-result"), "scan")
    assert app["trajectory_box"].layout.display == "none"
    assert app["artifact_fold"].layout.display == "none"
    assert "Run in progress" in app["results_empty"].value

    app["all_mode"].value = "tsonly"
    assert app["_ALL_MODE_STATE"]["user"] is True
    app["_queue_change"](clear=True)
    app["load_pdb"]([str(primary)], parm=str(topology))
    assert app["all_mode"].value == "scan"
    app["load_pdb"]([str(product)], parm=str(topology), append=True)
    assert app["all_mode"].value == "mep"

    first_xyz = tmp_path / "first.xyz"
    second_xyz = tmp_path / "second.xyz"
    xyz_text = "2\nenergy=-1 unit=hartree\nH 0 0 0\nH 0.7 0 0\n"
    first_xyz.write_text(xyz_text, encoding="utf-8")
    second_xyz.write_text(xyz_text, encoding="utf-8")
    app["_queue_change"](clear=True)
    app["_preflight_structures"] = lambda _paths: None
    app["_preflight_xyz_files"] = lambda _paths: None
    app["_preflight_xyz_reference_pairs"] = lambda _xyz, _refs: None
    app["build_selection"] = lambda: True
    assert app["_ingest_saved_files"](
        [str(first_xyz), str(second_xyz), str(primary), str(product)],
        "endpoint references",
    )
    assert app["S"]["subcmd"] == "path-search"
    assert app["dd_subcmd"].value == "path-search"

    app["S"]["_session_file_identities"] = {
        "parm": {"sha256": "0" * 64, "size": 1, "name": "missing.parm7"},
    }
    app["S"]["_session_primary_atom_signatures"] = [("A", "LIG", "10", "", "C1")]
    app["_queue_change"](clear=True)
    assert app["S"]["_session_file_identities"] == {}
    assert app["S"]["_session_primary_atom_signatures"] == []

    surface = tmp_path / "surface.csv"
    surface.write_text("d1,d2,d3,energy\n1,1,1,0\n", encoding="utf-8")
    assert app["_ingest_saved_files"]([str(surface)], "plot-only")
    app["_sync_capability_controls"]()
    assert app["S"]["mode"] == "utility"
    for name in ("scan_panel", "freeze_acc", "center_panel", "charge_panel", "extract_panel"):
        assert app[name].layout.display == "none"


def test_document_iframes_use_one_script_independent_srcdoc_payload(monkeypatch, tmp_path) -> None:
    """Widget refreshes cannot strand an iframe behind an inert loader script."""
    app, _rendered = _execute_app(monkeypatch, tmp_path)

    documents = [
        "<!doctype html><html><body>rx-small-unique<script>void 0;</script></body></html>",
        "<!doctype html><html><body>rx-large-unique" + ("<p>&lt;payload&gt;</p>" * 60_000)
        + "</body></html>",
    ]
    for upload_mode in ("colab", "basic"):
        app["_UPLOAD_MODE"] = upload_mode
        for document in documents:
            marker = "rx-large-unique" if "rx-large-unique" in document else "rx-small-unique"
            frame = app["_document_iframe"](
                document, 'class="rxprobe" data-rx-channel="viewer"', "width:100%;",
            )
            assert frame.count("<iframe") == 1 and frame.count("</iframe>") == 1
            assert 'srcdoc="' in frame
            assert 'class="rxprobe"' in frame and 'data-rx-channel="viewer"' in frame
            assert "URL.createObjectURL" not in frame
            assert "MutationObserver" not in frame
            assert "<script>(function(){var self=document.currentScript" not in frame
            encoded = re.search(r'srcdoc="([^"]*)"', frame)
            assert encoded is not None
            srcdoc = html.unescape(encoded.group(1))
            if len(document) > app["_DOC_FRAME_INLINE_LIMIT"]:
                packed = re.search(r'data-rx-document="([^"]+)"', frame)
                assert packed is not None and marker not in srcdoc
                assert "DecompressionStream('gzip')" in srcdoc
                assert gzip.decompress(base64.b64decode(packed.group(1))).decode() == document
            else:
                assert frame.count(marker) == 1 and srcdoc == document

        viewer = app["_molstar_iframe"]("HETATM\n", "pdb", interactive=True, generation=3)
        assert 'srcdoc="' in viewer
        assert "frame._rxObjectUrl" not in viewer
        assert 'data-rx-generation="3"' in viewer and 'class="rxmolstar-frame"' in viewer
        energy = app["_energy_plot_iframe"]([0.0, 1.0], {}, 4)
        assert 'srcdoc="' in energy
        assert "frame._rxObjectUrl" not in energy
        assert 'data-rx-channel="trajectory"' in energy and 'data-rx-generation="4"' in energy

    pdf = app["_binary_iframe"]("application/pdf", "AAAA", 'title="PDF result"', "width:100%;")
    assert 'src="data:application/pdf;base64,AAAA"' in pdf
    assert "atob(" in pdf and 'frame.getAttribute("src")' in pdf


def test_plotly_documents_use_one_stable_payload(monkeypatch, tmp_path) -> None:
    """Large Plotly HTML stays renderable after a widget value refresh."""
    app, _rendered = _execute_app(monkeypatch, tmp_path)
    document = "<main>RX_PLOTLY_UNIQUE" + ("x" * 512_000) + "</main>"

    for upload_mode in ("colab", "basic"):
        app["_UPLOAD_MODE"] = upload_mode
        frame = app["_plotly_document_iframe"](
            document, 'class="rxscan-pes-frame"', "width:100%;",
        )
        assert 'srcdoc="' in frame
        assert "URL.createObjectURL" not in frame
        assert "MutationObserver" not in frame
        assert frame.count("RX_PLOTLY_UNIQUE") == 1
        assert frame.count("<iframe") == 1 and frame.count("</iframe>") == 1
        encoded = re.search(r'srcdoc="([^"]*)"', frame)
        assert encoded is not None and html.unescape(encoded.group(1)) == document


def test_colab_uma_login_accepts_a_colab_secret(monkeypatch) -> None:
    """A stored Colab secret signs in without showing the token prompt."""
    import importlib.metadata

    setup = _notebook()["cells"][1]["source"].replace(
        'backend = "uma"', 'backend = "uma"', 1,
    ).replace("install_dft = True", "install_dft = False", 1)
    logins: list[tuple] = []
    fake_hf = types.ModuleType("huggingface_hub")
    fake_hf.login = lambda **kwargs: logins.append(("token", kwargs))
    fake_hf.notebook_login = lambda: logins.append(("notebook", {}))
    fake_hf.get_token = lambda: None
    userdata = types.ModuleType("google.colab.userdata")
    userdata.get = lambda name: "secret-token" if name == "HF_TOKEN" else None
    colab = types.ModuleType("google.colab")
    colab.userdata = userdata
    google = types.ModuleType("google")
    google.colab = colab
    monkeypatch.setitem(sys.modules, "huggingface_hub", fake_hf)
    monkeypatch.setitem(sys.modules, "google", google)
    monkeypatch.setitem(sys.modules, "google.colab", colab)
    monkeypatch.setitem(sys.modules, "google.colab.userdata", userdata)
    monkeypatch.setattr(
        subprocess, "run",
        lambda argv, **_kwargs: types.SimpleNamespace(stdout="GPU 0", stderr="", returncode=0),
    )
    monkeypatch.setattr(importlib.util, "find_spec", lambda _name: object())
    monkeypatch.setattr(
        importlib.metadata, "version",
        lambda name: "0.3.3" if name == "mlmm-toolkit" else "test",
    )
    monkeypatch.delenv("HF_TOKEN", raising=False)

    exec(compile(setup, str(NOTEBOOK), "exec"), {})

    assert [entry[0] for entry in logins] == ["token"]
    assert logins[0][1]["token"] == "secret-token"


def test_extract_panel_is_shown_only_for_workflows_that_need_it() -> None:
    """`all` extracts internally and `extract` IS the extraction command, so the
    preparation panel is hidden for both, with no extra notice in its place."""
    app = _notebook()["cells"][2]["source"]

    assert "extract_panel])" in app
    assert "extract_hint" not in app
    assert "_extract_panel.layout.display = '' if _prep_active else 'none'" in app
    # The button says what it does.
    assert "description='Extract ML region & use it'" in app


def test_option_controls_report_their_effective_default() -> None:
    """Every Options control names the default it falls back to, in one format.

    Dropdowns included: a bare "(default)" entry told the user nothing about
    what omitting the flag actually does. The label reads `show_default` first,
    because an option whose real default lives in a config block is declared
    `default=None` so an explicit value stays distinguishable from an omission.
    """
    app = _notebook()["cells"][2]["source"]

    assert "return 'default: %s' % value" in app
    assert "shown = getattr(param, 'show_default', None)" in app
    assert "if isinstance(shown, str) and shown: value = shown" in app
    assert "elif default is None: value = 'None'" in app
    # Hand-built dropdowns carry the same format as the generated ones.
    for label in ("('default: gsm', '(default)')", "('default: gpu', '(default)')",
                  "('default: gau', '(default)')", "('default: baker', '(default)')",
                  "('default: per backend', 'auto')",
                  "('default: 2', _ADVANCED_DEFAULT)"):
        assert label in app, label
    # One format only.
    assert "CLI default \u00b7" not in app


def test_advanced_dropdowns_select_their_default_label() -> None:
    source = _notebook()["cells"][2]["source"]
    wanted = {"_cli_default_label", "_set_advanced_override", "_advanced_widget"}
    functions = [
        node for node in ast.parse(source).body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    import click
    import ipywidgets as W

    state = {"advanced_overrides": {"all": {}}}
    namespace = {
        "S": state, "W": W, "click": click,
        "_ADVANCED_DEFAULT": "(default)",
        "_advanced_flag": lambda param: param.opts[0],
        "refresh": lambda: None,
        "_flag_row": lambda widget, _help: types.SimpleNamespace(children=(widget,)),
    }
    exec(compile(ast.Module(body=functions, type_ignores=[]),
                 str(NOTEBOOK), "exec"), namespace)

    params = [
        click.Option(["--verbose"], type=int, default=2, show_default="2"),
        click.Option(["--feature/--no-feature"], default=False,
                     show_default="off"),
        click.Option(["--mode"], type=click.Choice(["a", "b"]),
                     default="a", show_default="a"),
    ]
    for param in params:
        widget = namespace["_advanced_widget"]("all", param).children[0]
        assert widget.index == 0
        assert widget.value == "(default)"
        assert widget.label.startswith("default:")
        widget.value = widget.options[1][1]
        assert state["advanced_overrides"]["all"][param.name] == widget.value
        widget.value = "(default)"
        assert param.name not in state["advanced_overrides"]["all"]


def test_the_gui_keeps_one_run_path() -> None:
    """The cell used to carry an unreachable async twin of the whole execution
    path (`_start_async_task`, `_do_validate_async`, `_do_run_guarded`,
    `_do_run_async`, `_stream_async`, `_stop_async_process`,
    `_validate_command_async`, `_async_task_done`): nothing called any of them,
    so a reader had two implementations to reason about and only one ever ran.
    One path only, so the live behaviour is the readable one.
    """
    app = _notebook()["cells"][2]["source"]

    for gone in ("_stream_async", "_stop_async_process", "_validate_command_async",
                 "_do_validate_async", "_async_task_done", "_start_async_task",
                 "_do_run_async", "_do_run_guarded"):
        assert gone not in app, gone
    # The live path, and the loop reference `_dispatch_ui` still needs.
    assert "def _stream(cmd):" in app
    assert "_RUN_LOG_VIEW_LIMIT" not in app
    assert "[earlier output trimmed from this view]" not in app
    assert "def _persist_run_log(out, transcript):" in app
    assert "if sub in COMPUTE: run_log_path = _persist_run_log(out, transcript)" in app
    assert "def _do_run_sync(" in app
    assert "subprocess.Popen(cmd" in app
    assert "_tensor_warning_filter = ('ignore:using cupy as the tensor contraction '" in app
    assert "'PYTHONWARNINGS': ','.join(" in app
    assert "env=_child_env" in app
    assert "asyncio.create_subprocess_exec" not in app
    assert "_UI_ASYNC_LOOP.call_soon_threadsafe(callback)" in app


def test_colab_poll_repairs_a_finished_but_stale_frontend() -> None:
    source = _notebook()["cells"][2]["source"]
    tree = ast.parse(source)
    poll_node = next(
        node for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and node.name == "_colab_poll_run"
    )
    events = []
    namespace = {
        "os": os,
        "_RUN_EXECUTION": {
            "task": None, "thread": None, "process": None,
            "result_attempt": True,
        },
        "_ACTION_STATE": {"running": False},
        "_RUN_STATE": {"text": "✓ done", "tone": "ok", "kind": "done"},
        "S": {"_last_manifest": {"status": "success"},
              "_last_out_dir": "result",
              "_results_presented_dir": os.path.abspath("result")},
        "_flush_run_log": lambda force=False: events.append(("flush", force)),
        "_set_running": lambda value, on_loop=False: events.append(
            ("running", value, on_loop)),
        "_set_run_status": lambda text, tone, kind, on_loop=False: events.append(
            ("status", text, tone, kind, on_loop)),
        "_results": lambda out: events.append(("results", out)),
        "_finish_cancelled_attempt": lambda out: events.append(
            ("cancelled_results", out)),
        "_tab_go": lambda index: events.append(("tab", index)),
        "_publish_result_widget_state": lambda: events.append(("publish_results",)),
        "_publish_run_widget_state": lambda: events.append(("publish_run",)),
    }
    exec(compile(ast.Module(body=[poll_node], type_ignores=[]),
                 str(NOTEBOOK), "exec"), namespace)

    assert "bridge.invokeFunction(CONFIG.poll_callback,[active],{})" in source
    assert namespace["_colab_poll_run"](False) == {"running": False}
    assert events == [("flush", True)]

    events.clear()
    assert namespace["_colab_poll_run"](True) == {"running": False}
    assert events == [
        ("flush", True),
        ("running", False, True),
        ("status", "✓ done", "ok", "done", True),
        ("tab", 3),
        ("publish_run",),
    ]

    events.clear()
    namespace["_RUN_EXECUTION"]["result_attempt"] = False
    namespace["S"]["_last_manifest"]["status"] = "success"
    namespace["S"]["_results_presented_dir"] = os.path.abspath("result")
    namespace["_RUN_STATE"].update(
        text="■ cancelled", tone="warn", kind="cancelled")
    assert namespace["_colab_poll_run"](True) == {"running": False}
    assert events == [
        ("flush", True),
        ("running", False, True),
        ("status", "■ cancelled", "warn", "cancelled", True),
        ("publish_run",),
    ]

    # A scientifically partial result is still a terminal, displayable run.
    events.clear()
    namespace["_RUN_EXECUTION"]["result_attempt"] = True
    namespace["_RUN_STATE"].update(text="✓ done", tone="ok", kind="done")
    namespace["S"]["_last_manifest"]["status"] = "partial"
    namespace["S"]["_results_presented_dir"] = None
    assert namespace["_colab_poll_run"](False) == {"running": False}
    assert ("tab", 3) in events

    events.clear()
    namespace["S"]["_last_manifest"]["status"] = "cancelled"
    namespace["S"]["_results_presented_dir"] = None
    namespace["_RUN_STATE"].update(
        text="■ cancelled", tone="warn", kind="cancelled")
    assert namespace["_colab_poll_run"](False) == {"running": False}
    assert events == [
        ("flush", True),
        ("running", False, True),
        ("status", "■ cancelled", "warn", "cancelled", True),
        ("cancelled_results", "result"),
        ("publish_run",),
    ]


def test_results_render_only_after_results_tab_owns_the_view(
    monkeypatch, tmp_path: Path,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)
    result_root = tmp_path / "result"
    result_root.mkdir()
    app["S"].update(
        _last_manifest={"status": "success"},
        _last_files=[str(result_root / "summary.log")],
        _last_out_dir=str(result_root),
        _results_presented_dir=None,
    )
    events = []
    app["_results"] = lambda out: events.append(
        ("render", app["_TAB_NAV"]["active"], out))
    app["_publish_result_widget_state"] = lambda: events.append(
        ("publish", app["_TAB_NAV"]["active"]))
    display_events = {index: [] for index in range(4)}
    for index, (_label, pane) in enumerate(app["_TAB_PAGES"]):
        pane.layout.observe(
            lambda change, index=index: display_events[index].append(change["new"]),
            names="display",
        )

    app["_tab_go"](3)

    assert events == [
        ("render", 3, str(result_root)),
        ("publish", 3),
    ]
    assert [pane.layout.display for _label, pane in app["_TAB_PAGES"]] == [
        "none", "flex", "none", "flex",
    ]
    assert all(values[-2:] == ["block", "flex" if index in (1, 3) else "none"]
               for index, values in display_events.items())


def test_cancelled_run_stays_on_the_active_tab_without_result_rendering(
    tmp_path: Path, monkeypatch,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)
    input_path = tmp_path / "input.pdb"
    input_path.write_text("END\n", encoding="utf-8")
    partial = tmp_path / "result" / "partial_mep_trj.xyz"
    partial.parent.mkdir()
    partial.write_text("1\npartial\nH 0 0 0\n", encoding="utf-8")
    scope = {
        "target": str(partial.parent), "root": str(partial.parent),
        "shallow": False, "stdout_only": False, "direct_current": True,
    }
    snapshots = iter([{}, {str(partial): (1, 1)}])
    app["_auto"]["on"] = False
    app["_preflight_output_scope"] = lambda _argv: scope
    app["_snapshot_output_scope"] = lambda _scope: next(snapshots)
    app["_matches_output_scope"] = lambda _path, _scope: True
    app["_stream"] = lambda _argv: (130, "cancelled")
    rendered = []
    tab_calls = []
    real_results = app["_results"]
    app["_results"] = lambda out: rendered.append(out)
    published_results = []
    app["_publish_result_widget_state"] = lambda: published_results.append(True)
    real_tab_go = app["_tab_go"]
    real_tab_go(2)
    app["_tab_go"] = lambda index: (tab_calls.append(index), real_tab_go(index))[-1]

    app["_do_run_sync"](None, [app["CLI"], "fix-altloc", "-i", str(input_path)])

    assert app["S"]["_last_manifest"]["status"] == "cancelled"
    assert app["_RUN_EXECUTION"]["result_attempt"] is True
    assert app["S"]["_last_files"] == [str(partial)]
    assert app["S"]["_results_presented_dir"] == os.path.abspath(partial.parent)
    assert app["_TAB_NAV"]["active"] == 2
    assert tab_calls == [] and rendered == [] and published_results == []
    assert app["trajectory_box"].layout.display == "none"
    assert app["artifact_fold"].layout.display == "none"
    assert "No completed result" in app["results_empty"].value
    assert app["dl_btn"].disabled
    app["S"]["_results_presented_dir"] = None
    real_results(str(partial.parent))
    assert rendered == []
    assert app["trajectory_box"].layout.display == "none"
    assert not app["traj_choice"].options
    real_tab_go(3)
    assert rendered == []
    assert "No completed result" in app["results_empty"].value


def test_cancel_keeps_polling_until_worker_publishes_terminal_state() -> None:
    source = _notebook()["cells"][2]["source"]
    tree = ast.parse(source)
    cancel_node = next(
        node for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and node.name == "_cancel_run"
    )

    class LiveThread:
        @staticmethod
        def is_alive() -> bool:
            return True

    class CancelEvent:
        def set(self) -> None:
            events.append(("cancel",))

    process = object()
    events = []
    namespace = {
        "_RUN_EXECUTION": {
            "task": None,
            "thread": LiveThread(),
            "process": process,
            "cancel": CancelEvent(),
        },
        "_set_run_status": lambda text, tone, kind: events.append(
            ("status", text, tone, kind)
        ),
        "_stop_child": lambda child: events.append(("stop", child)),
        "_set_running": lambda *args: events.append(("running",) + args),
        "_publish_run_widget_state": lambda: events.append(("publish",)),
        "os": __import__("os"),
        "signal": __import__("signal"),
    }
    exec(compile(ast.Module(body=[cancel_node], type_ignores=[]),
                 str(NOTEBOOK), "exec"), namespace)

    namespace["_cancel_run"](None)
    assert events == [
        ("cancel",),
        ("status", "■ cancelling…", "warn", "cancelling"),
        ("stop", process),
    ]


def test_results_replaces_trajectory_with_exact_stationary_model_set(
    tmp_path: Path, monkeypatch,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)
    assert app["_normalise_result_label"]("MLIP\u00a0ΔG") == "MLIP ΔG"
    assert app["_normalise_result_label"](
        "Energy\u00a0profile  &\tTrajectory"
    ) == "Energy profile & Trajectory"

    trajectory = tmp_path / "irc_trj.xyz"
    trajectory.write_text("".join(
        "2\nTRJ_FRAME_%d energy=%.6f unit=hartree\nH %d 0 0\nH %d.7 0 0\n"
        % (index, -1.0 + index * 0.001, index, index)
        for index in range(5)
    ), encoding="utf-8")
    app["S"]["_last_subcmd"] = "all"
    app["_load_trajectory"](str(trajectory), str(tmp_path))
    assert len(app["_TRAJ"]["frames"]) == 5
    assert "TRJ_FRAME_4" in _embedded_document(app["traj_out"].value)

    structures = []
    for index, label in enumerate(("R", "TS1", "P")):
        path = tmp_path / (label + ".xyz")
        path.write_text(
            "2\n%s_ONLY\nH %d 0 0\nH %d.7 0 0\n" % (label, index, index),
            encoding="utf-8",
        )
        structures.append(str(path))
    view = {
        "structures": structures,
        "labels": ["R", "TS1", "P"],
        "energies_kcal": [0.0, 12.0, -1.0],
    }

    assert app["_load_energy_structures"](view, "mlip_g", "MLIP ΔG")
    stationary_viewer = _embedded_document(app["traj_out"].value)
    assert app["_TRAJ"]["mode"] == "levels"
    assert len(app["_TRAJ"]["frames"]) == 3
    assert app["_TRAJ"]["state_labels"] == ["R", "TS1", "P"]
    assert all(marker in stationary_viewer for marker in ("R_ONLY", "TS1_ONLY", "P_ONLY"))
    assert "TRJ_FRAME_" not in stationary_viewer
    assert app["_TRAJ"]["structure_message"] == {}

    app["_load_trajectory"](str(trajectory), str(tmp_path))
    assert len(app["_TRAJ"]["frames"]) == 5
    assert "R_ONLY" not in _embedded_document(app["traj_out"].value)
    assert app["_load_energy_structures"](view, "mlip_g", "MLIP ΔG")
    assert len(app["_TRAJ"]["frames"]) == 3
    assert "TRJ_FRAME_" not in _embedded_document(app["traj_out"].value)

    # Reloading a completed run with a static energy image must also replace
    # a previously selected imaginary-mode trajectory with the primary MEP.
    stale_mode = tmp_path / "segments" / "seg_01" / "ts" / "vib" / "imag_-550.0cm-1_trj.xyz"
    stale_mode.parent.mkdir(parents=True)
    stale_mode.write_text(
        "2\nSTALE_MODE_0\nH 0 0 0\nH 0.7 0 0\n"
        "2\nSTALE_MODE_1\nH 0.1 0 0\nH 0.8 0 0\n",
        encoding="utf-8",
    )
    mep = tmp_path / "mep_trj.xyz"
    mep.write_text(trajectory.read_text(encoding="utf-8"), encoding="utf-8")
    energy_image = tmp_path / "energy_diagram_G_MLIP.png"
    energy_image.write_bytes(b"\x89PNG\r\n\x1a\n")
    app["S"].update({
        "_last_out_dir": str(tmp_path),
        "_last_subcmd": "all",
        "_last_files": [str(stale_mode), str(mep), str(energy_image)],
        "_last_manifest": {"status": "success", "exit_code": 0},
    })
    app["_load_trajectory"](str(stale_mode), str(tmp_path))
    assert "STALE_MODE_1" in _embedded_document(app["traj_out"].value)
    app["_results"](str(tmp_path))
    energy_token = app["energy_choice"].value
    assert app["_ENERGY"]["views"][energy_token]["mode"] == "image"
    assert app["traj_choice"].value == str(mep)
    assert app["_TRAJ"]["path"] == str(mep)
    assert len(app["_TRAJ"]["frames"]) == 5
    assert "MEP profile" in app["trajectory_intro"].value
    assert "Vibrational mode" not in app["trajectory_intro"].value
    assert "STALE_MODE_" not in _embedded_document(app["traj_out"].value)

    source = _notebook()["cells"][2]["source"]
    assert "selected_energy_view.get('mode') == 'image'" in source
    assert "_load_trajectory(primary_path, out)" in source
    assert "white-space:normal" in source and "overflow-wrap:anywhere" in source
    assert "flex_flow='row wrap'" in source
    profile = html.unescape(app["_energy_plot_document"](
        [0.0, 10.0, 5.0],
        {"x": "Path image", "start": "R", "end": "P", "extrema": True}, 17,
    ))
    assert '"xRange":[0.5,3.5]' in profile
    assert '"yRange":[-1.2,11.2]' in profile
    assert '"xTickStep":1' in profile
    assert "range:cfg.xRange,autorange:false" in profile
    assert "range:cfg.yRange,autorange:false" in profile
    assert profile.count("if(!cfg.linked)return;") == 2
    medium = app["_energy_plot_document"](
        [float(index) for index in range(22)],
        {"x": "Path image", "start": "R", "end": "P", "extrema": True}, 18,
    )
    many = app["_energy_plot_document"](
        [float(index) for index in range(51)],
        {"x": "Intrinsic reaction coordinate", "start": "R", "end": "P", "extrema": True}, 19,
    )
    assert '"xTickStep":5' in medium
    assert '"xTickStep":10' in many
    assert "showticklabels:true,ticks:'',ticklen:0,tickfont:{size:14}" in profile
    assert "font:{size:18,color:'#253047'}" in profile
    assert "const below=/endpoint$/i.test" in profile
    assert "yshift:below?-18:16,yanchor:below?'top':'bottom'" in profile
    irc_profile = html.unescape(app["_energy_plot_document"](
        [0.0, 10.0, 5.0],
        {"x": "IRC point", "start": "forward endpoint",
         "end": "backward endpoint", "ts_index": 1, "ts_label": "TS"}, 20,
    ))
    assert '"yRange":[-2.2,12.2]' in irc_profile


def test_results_route_single_structures_modes_and_scan_grids(
    monkeypatch, tmp_path: Path,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)
    assert "Structure Viewer" in app["structure_panel_title"].value
    assert "Distance and angle measurements update" not in app["structure_panel_title"].value
    assert "Distance and angle measurements update" in app["structure_panel_info"]._rx_tip
    assert "rxstructure-panel-head" in app["structure_panel_head"]._dom_classes
    xyz = (
        "2\nenergy=-1.000000 unit=hartree\nH 0 0 0\nH 0.7 0 0\n"
        "2\nenergy=-0.990000 unit=hartree\nH 0.1 0 0\nH 0.8 0 0\n"
    )
    final = tmp_path / "final_geometry.xyz"
    final.write_text(
        "2\nenergy=-1.23456789 unit=hartree\nH 0 0 0\nH 0.7 0 0\n",
        encoding="utf-8",
    )
    app["S"].update(
        _last_out_dir=str(tmp_path), _last_subcmd="opt",
        _last_files=[str(final)], _last_manifest={"status": "success", "exit_code": 0},
    )
    app["_results"](str(tmp_path))
    assert app["traj_choice"].options[0][0] == "Final optimized structure"
    assert len(app["_TRAJ"]["frames"]) == 1
    assert app["frame_play"].disabled is True
    assert app["frame_controls"].layout.display == "none"
    assert app["energy_panel"].layout.display == "none"
    assert "rxstructure-only" in app["path_grid"]._dom_classes
    assert "Energy = -1.23456789 Ha" in app["frame_state"].value

    optimization = tmp_path / "optimization_all_trj.xyz"
    optimization.write_text(xyz, encoding="utf-8")
    app["S"].update(
        _last_subcmd="opt", _last_files=[str(final), str(optimization)],
        _last_manifest={"status": "success", "exit_code": 0},
    )
    app["_results"](str(tmp_path))
    assert [label for label, _ in app["traj_choice"].options] == [
        "Optimization trajectory", "Final optimized structure",
    ]
    assert len(app["_TRAJ"]["frames"]) == 2
    assert app["frame_controls"].layout.display == ""
    assert app["energy_panel"].layout.display == ""
    assert "rxstructure-only" not in app["path_grid"]._dom_classes
    assert "Optimization progress" in app["trajectory_intro"].value
    assert "rxenergy-frame" in app["plot_out"].value
    assert [label for label, _ in app["_result_view_candidates"](
        [str(final), str(optimization)], str(tmp_path), "tsopt"
    )] == ["TS-refinement trajectory", "Refined transition-state structure"]

    vib_dir = tmp_path / "vib"
    vib_dir.mkdir()
    imaginary = vib_dir / "imag_-550.15cm-1_trj.xyz"
    imaginary.write_text(xyz, encoding="utf-8")
    assert [label for label, _ in app["_result_view_candidates"](
        [str(final), str(optimization), str(imaginary)], str(tmp_path), "tsopt"
    )] == [
        "TS-refinement trajectory", "Refined transition-state structure",
        "Imaginary mode · −550.15 cm⁻¹",
    ]

    segment_vib = tmp_path / "segments" / "seg_02" / "ts" / "vib"
    segment_vib.mkdir(parents=True)
    segment_imaginary = segment_vib / "imag_-421.25cm-1_trj.xyz"
    segment_imaginary.write_text(xyz, encoding="utf-8")
    weaker_segment_imaginary = segment_vib / "imag_-120.00cm-1_trj.xyz"
    weaker_segment_imaginary.write_text(xyz, encoding="utf-8")
    fallback_vib = tmp_path / "segments" / "seg_03" / "ts" / "vib"
    fallback_vib.mkdir(parents=True)
    fallback_imaginary = fallback_vib / "mode_1_-333.50cm-1_trj.xyz"
    fallback_imaginary.write_text(xyz, encoding="utf-8")
    segment_irc_dir = tmp_path / "segments" / "seg_02" / "irc"
    segment_irc_dir.mkdir(parents=True)
    segment_irc = segment_irc_dir / "finished_irc_trj.xyz"
    segment_irc.write_text(xyz, encoding="utf-8")
    segment_mep = tmp_path / "mep_seg_02_trj.xyz"
    segment_mep.write_text(xyz, encoding="utf-8")
    mep = tmp_path / "mep_trj.xyz"
    mep.write_text(xyz, encoding="utf-8")
    all_views = app["_result_view_candidates"](
        [str(mep), str(segment_imaginary), str(weaker_segment_imaginary),
         str(fallback_imaginary), str(segment_irc), str(segment_mep)],
        str(tmp_path), "all"
    )
    assert [label for label, _ in all_views] == [
        "Reaction-path trajectory",
        "Imaginary mode · TS2 · −421.25 cm⁻¹",
        "Imaginary mode · TS2 · −120.00 cm⁻¹",
        "Imaginary mode · TS3 · −333.50 cm⁻¹",
    ]

    app["S"].update(
        _last_subcmd="all",
        _last_files=[str(mep), str(segment_imaginary), str(fallback_imaginary)],
        _last_manifest={"status": "success", "exit_code": 0},
    )
    app["_results"](str(tmp_path))
    app["traj_choice"].value = str(segment_imaginary)
    assert app["_TRAJ"]["path"] == str(segment_imaginary)
    assert len(app["_TRAJ"]["frames"]) == 2
    assert "Vibrational mode" in app["trajectory_intro"].value
    assert app["energy_panel"].layout.display == "none"
    assert "rxstructure-only" in app["path_grid"]._dom_classes
    assert app["plot_out"].value == ""

    mode = tmp_path / "mode_0001_-100.00cm-1_trj.xyz"
    mode.write_text(xyz, encoding="utf-8")
    app["S"].update(
        _last_subcmd="freq", _last_files=[str(mode)],
        _last_manifest={"status": "success", "exit_code": 0},
    )
    app["_results"](str(tmp_path))
    assert app["traj_choice"].options[0][0].startswith("Mode 1 of 1 · −100.00 cm⁻¹")
    assert app["result_selector_row"].layout.display == ""
    assert len(app["_TRAJ"]["frames"]) == 2
    assert app["frame_play"].disabled is False
    assert "Vibrational mode" in app["trajectory_intro"].value
    assert app["energy_panel"].layout.display == "none"
    assert "rxstructure-only" in app["path_grid"]._dom_classes
    assert app["plot_out"].value == ""

    grid = tmp_path / "grid"
    grid.mkdir()
    point_a = grid / "point_i100_j150.xyz"
    point_b = grid / "point_i110_j160.xyz"
    point_a.write_text("2\npoint A\nH 0 0 0\nH 0.7 0 0\n", encoding="utf-8")
    point_b.write_text("2\npoint B\nH 0 0 0\nH 0.8 0 0\n", encoding="utf-8")
    surface = tmp_path / "surface.csv"
    surface.write_text(
        "i,j,d1_A,d2_A,target_d1_A,target_d2_A,energy_kcal,is_preopt\n"
        "0,0,1.01,1.49,1.00,1.50,0.0,False\n"
        "1,0,1.09,1.61,1.10,1.60,2.5,False\n", encoding="utf-8",
    )
    landscape = tmp_path / "scan2d_landscape.html"
    landscape.write_text(
        "<!doctype html><html><body><div class='plotly-graph-div'>original landscape</div></body></html>",
        encoding="utf-8",
    )
    files = [surface, landscape, point_a, point_b]
    app["S"].update(
        _last_subcmd="scan2d", _last_files=[str(path) for path in files],
        _last_manifest={"status": "success", "exit_code": 0},
    )
    app["_results"](str(tmp_path))
    assert app["_TRAJ"]["mode"] == "grid"
    assert len(app["_TRAJ"]["frames"]) == 2
    assert app["_SCAN_GRID"]["view"]["points"][0]["coords"] == [1.0, 1.5]
    assert app["_SCAN_GRID"]["view"]["points"][0]["plot_coords"] == [1.01, 1.49]
    assert "d1 1.00 Å ↔ d2 1.50 Å" in app["frame_state"].value
    assert "ΔE = 0.0 kcal/mol" in app["frame_state"].value
    assert "?" not in app["frame_state"].value
    assert "original landscape" in app["plot_out"].value
    assert "Selected structure" not in app["plot_out"].value
    assert "max-width:100%" in app["plot_out"].value
    assert "Plotly.relayout(plot,{autosize:true,height:height})" in app["plot_out"].value
    linked_plot = html.unescape(app["plot_out"].value)
    assert "Plotly.restyle(graph,{hoverinfo:'skip'},surfaces)" in linked_plot
    assert "Plotly.addTraces(graph,[pointTrace()])" in linked_plot
    assert "marker:{size:1.8" in linked_plot
    assert "item.plot_coords||item.coords" in linked_plot
    assert "rx-set-frame" not in app["plot_out"].value
    assert app["frame_controls"].layout.display == "none"
    assert app["scan_axis_controls"].layout.display == ""
    assert app["scan_axis_rows"][0].layout.display == ""
    assert app["scan_axis_rows"][1].layout.display == ""
    assert app["scan_axis_rows"][2].layout.display == "none"
    app["scan_axis_sliders"][0].value = 1.10
    # The requested tuple (1.10, 1.50) is absent from this partial grid.  Do not
    # substitute a nearest point and do not move any other coordinate.
    assert "Grid point 1 of 2" in app["frame_state"].value
    assert app["frame_slider"].value == 0
    assert app["scan_axis_sliders"][0].value == 1.10
    assert app["scan_axis_sliders"][1].value == 1.50
    app["_set_frame_from_browser"](app["_TRAJ"]["generation"], 0)
    assert app["frame_slider"].value == 0
    assert app["scan_axis_sliders"][0].value == 1.00
    assert app["scan_axis_sliders"][1].value == 1.50
    # Colab can deliver a stale selection index after Results options shrink.
    # Ignore that stale state: selecting a different valid result would silently
    # desynchronize the plot, structure, and coordinate controls.
    app["_result_pick_guard"]["active"] = True
    try:
        app["traj_choice"].options = [("Only result", str(point_a))]
        choice_before = app["traj_choice"].index
        app["traj_choice"].set_state({"index": 999})
    finally:
        app["_result_pick_guard"]["active"] = False
    assert app["traj_choice"].index == choice_before
    generation_before = app["_TRAJ"]["generation"]
    frame_before = app["frame_slider"].value
    axis_index_before = app["scan_axis_sliders"][0].index
    axis_value_before = app["scan_axis_sliders"][0].value
    app["scan_axis_sliders"][0].set_state({"index": 999})
    assert app["scan_axis_sliders"][0].index == axis_index_before
    assert app["scan_axis_sliders"][0].value == axis_value_before
    assert app["frame_slider"].value == frame_before
    assert app["_TRAJ"]["generation"] == generation_before

    # Colab's slider bridge can serialize an integer index with harmless
    # floating-point noise.  Normalize that before traitlets' Int type check,
    # but ignore genuinely fractional or out-of-range stale messages.
    noisy_slider = app["_SafeResultSelectionSlider"](
        options=[(str(i), i) for i in range(40)], value=0
    )
    noisy_slider.set_state({"index": 34.99999999999999})
    assert noisy_slider.index == 35
    noisy_slider.set_state({"index": 34.5})
    assert noisy_slider.index == 35
    noisy_slider.set_state({"index": 999.0})
    assert noisy_slider.index == 35

    noisy_dropdown = app["_SafeResultDropdown"](
        options=[("zero", 0), ("one", 1)], value=0
    )
    noisy_dropdown.set_state({"index": 1.0000000000000002})
    assert noisy_dropdown.index == 1
    noisy_dropdown.set_state({"index": 0.5})
    assert noisy_dropdown.index == 1

    app["_configure_scan_axis_controls"]({
        "dims": 3,
        "points": [{"coords": (1.0, 1.5, 2.0)}, {"coords": (1.1, 1.6, 2.1)}],
    })
    assert [row.layout.display for row in app["scan_axis_rows"]] == ["", "", ""]

    # A one-point grid keeps each real coordinate selected in a disabled,
    # non-operable control.  A second inert option prevents Colab/noUiSlider
    # from rejecting a collapsed one-index frontend range.
    app["_configure_scan_axis_controls"]({
        "dims": 3,
        "points": [{"coords": (1.23, 1.50, 2.00)}],
    })
    assert [tuple(slider.options) for slider in app["scan_axis_sliders"]] == [
        (("1.23 Å", 1.23), ("", None)),
        (("1.50 Å", 1.5), ("", None)),
        (("2.00 Å", 2.0), ("", None)),
    ]
    assert all(slider.disabled for slider in app["scan_axis_sliders"])
    generation_before = app["_TRAJ"]["generation"]
    frame_before = app["frame_slider"].value
    app["scan_axis_sliders"][0].set_state({"index": 999})
    assert app["scan_axis_sliders"][0].index == 0
    assert app["scan_axis_sliders"][0].value == 1.23
    assert app["frame_slider"].value == frame_before
    assert app["_TRAJ"]["generation"] == generation_before


def test_scan_grid_manifest_links_structures_and_html_is_a_plot_only_fallback(
    monkeypatch, tmp_path: Path,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)
    grid = tmp_path / "grid"
    grid.mkdir()
    point_a = grid / "point_i100_j150_grid_000_000.xyz"
    point_b = grid / "point_i1001_j1602_grid_001_000.xyz"
    for path, distance in ((point_a, 0.7), (point_b, 0.8)):
        path.write_text(
            f"2\nlinked point\nH 0 0 0\nH {distance} 0 0\n",
            encoding="utf-8",
        )
    surface = tmp_path / "surface.csv"
    surface.write_text(
        "i,j,d1_A,d2_A,energy_kcal,is_preopt\n"
        "0,0,1.000,1.500,0.0,False\n"
        "1,0,1.001,1.602,2.5,False\n",
        encoding="utf-8",
    )
    landscape = tmp_path / "scan2d_landscape.html"
    landscape.write_text(
        "<!doctype html><html><head></head><body>manifest landscape"
        "<div class='plotly-graph-div'></div></body></html>",
        encoding="utf-8",
    )
    result = tmp_path / "result.json"
    result.write_text(
        json.dumps(
            {
                "grid_points": [
                    {
                        "index": [0, 0],
                        "distances_angstrom": [1.0, 1.5],
                        "targets_angstrom": [1.0, 1.5],
                        "geometry_file": str(point_a.relative_to(tmp_path)),
                    },
                    {
                        "index": [1, 0],
                        "distances_angstrom": [1.001, 1.602],
                        "targets_angstrom": [1.001, 1.602],
                        "geometry_file": str(point_b.relative_to(tmp_path)),
                    },
                ]
            }
        ),
        encoding="utf-8",
    )
    app["S"].update(
        _last_out_dir=str(tmp_path),
        _last_subcmd="scan2d",
        _last_files=[str(surface), str(landscape), str(result)],
        _last_manifest={"status": "success", "exit_code": 0},
    )
    app["_results"](str(tmp_path))
    assert len(app["_TRAJ"]["frames"]) == 2
    assert "Computed grid points" in app["plot_out"].value
    assert "manifest landscape" in app["plot_out"].value
    assert "Interactive 2D PES" in app["energy_panel_title"].value
    assert "Plotly.relayout(graph,{'title.text':'','margin.t':20})" in html.unescape(app["plot_out"].value)
    assert "rxplot-only" not in app["path_grid"]._dom_classes

    plot_only = tmp_path / "plot-only"
    plot_only.mkdir()
    density = plot_only / "scan3d_density.html"
    density.write_text(
        "<!doctype html><html><head></head><body>plot-only density"
        "<div class='plotly-graph-div'></div></body></html>",
        encoding="utf-8",
    )
    app["S"].update(
        _last_out_dir=str(plot_only),
        _last_subcmd="scan3d",
        _last_files=[str(density)],
        _last_manifest={"status": "success", "exit_code": 0},
    )
    app["_results"](str(plot_only))
    assert "plot-only density" in app["plot_out"].value
    assert app["structure_panel"].layout.display == "none"
    assert app["frame_controls"].layout.display == "none"
    assert "rxplot-only" in app["path_grid"]._dom_classes
    assert "Structure linking is unavailable" in app["frame_state"].value


def test_scan3d_axes_plot_selection_and_structure_stay_bidirectionally_linked(
    monkeypatch, tmp_path: Path,
) -> None:
    """Every 3-D grid axis and the browser plot select the same structure."""
    app, _ = _execute_app(monkeypatch, tmp_path)
    grid = tmp_path / "grid"
    grid.mkdir()
    points = []
    files = []
    for i, d1 in enumerate((1.0, 1.1)):
        for j, d2 in enumerate((1.5, 1.6)):
            for k, d3 in enumerate((2.0, 2.1)):
                path = grid / f"point_{i}_{j}_{k}.xyz"
                path.write_text(
                    f"2\ngrid {i} {j} {k}\nH 0 0 0\nH {0.7 + .01*(i+j+k)} 0 0\n",
                    encoding="utf-8",
                )
                files.append(path)
                points.append({
                    "index": [i, j, k],
                    "distances_angstrom": [d1, d2, d3],
                    "targets_angstrom": [d1, d2, d3],
                    "energy_kcal": float(i + 2*j + 3*k),
                    "geometry_file": str(path.relative_to(tmp_path)),
                })
    density = tmp_path / "scan3d_density.html"
    density.write_text(
        "<!doctype html><html><body><div class='plotly-graph-div'>density</div></body></html>",
        encoding="utf-8",
    )
    result = tmp_path / "result.json"
    result.write_text(json.dumps({"grid_points": points}), encoding="utf-8")
    app["S"].update(
        _last_out_dir=str(tmp_path), _last_subcmd="scan3d",
        _last_files=[str(density), str(result), *map(str, files)],
        _last_manifest={"status": "success", "exit_code": 0},
    )
    app["_results"](str(tmp_path))

    assert app["_TRAJ"]["mode"] == "grid"
    assert len(app["_TRAJ"]["frames"]) == 8
    assert [row.layout.display for row in app["scan_axis_rows"]] == ["", "", ""]
    assert all(len(slider.options) == 2 for slider in app["scan_axis_sliders"])
    generation = app["_TRAJ"]["generation"]
    app["scan_axis_sliders"][0].value = 1.1
    assert [slider.value for slider in app["scan_axis_sliders"]] == [1.1, 1.5, 2.0]
    app["scan_axis_sliders"][1].value = 1.6
    assert [slider.value for slider in app["scan_axis_sliders"]] == [1.1, 1.6, 2.0]
    app["scan_axis_sliders"][2].value = 2.1
    assert "Grid point 8 of 8" in app["frame_state"].value
    assert app["frame_slider"].value == 7

    # A Plotly click reaches this registered browser callback.  The PES itself
    # remains unchanged while the sliders and molecular structure follow it.
    app["_set_frame_from_browser"](generation, 2)
    assert app["frame_slider"].value == 2
    assert [slider.value for slider in app["scan_axis_sliders"]] == [1.0, 1.6, 2.0]
    assert "grid 0 1 0" in app["_TRAJ"]["frames"][2]
    view = app["_scan_grid_view"](
        [str(density), str(result), *map(str, files)], str(tmp_path), "scan3d")
    document = app["_scan_grid_document"](view, generation)
    assert "graph.on('plotly_click'" in document
    assert "Plotly.relayout(graph,{'title.text':'','margin.t':20})" in document
    assert "invoke(index)" in document
    assert "nearest(point)" not in document
    assert "rx-set-frame" not in document
    assert "Plotly.restyle(graph,{hoverinfo:'skip'},surfaces)" in document
    assert "Plotly.addTraces(graph,[pointTrace()])" in document
    assert "marker:{size:1.8" in document


def test_results_playback_uses_fifty_ms_for_mep_irc_and_imaginary_modes(
    monkeypatch, tmp_path: Path,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)

    assert app["_PROFILE_PLAYBACK_INTERVAL_MS"] == 50
    assert app["frame_play"].interval == 850
    app["_sync_playback_interval"](token="mep")
    assert app["frame_play"].interval == 50
    app["_sync_playback_interval"](path="finished_irc_trj.xyz")
    assert app["frame_play"].interval == 50
    app["_sync_playback_interval"](path="mode_0001_-100.00cm-1_trj.xyz")
    assert app["frame_play"].interval == 50
    app["_sync_playback_interval"](token="vibration")
    assert app["frame_play"].interval == 50
    app["_sync_playback_interval"](token="mlip_g")
    assert app["frame_play"].interval == 850
    app["_sync_playback_interval"](path="scan_trj.xyz")
    assert app["frame_play"].interval == 850

    irc_plot = tmp_path / "irc_plot_all.png"
    irc_plot.write_bytes(b"\x89PNG\r\n\x1a\n")
    segment_irc = tmp_path / "segments" / "seg_01" / "irc" / "finished_irc_trj.xyz"
    segment_irc.parent.mkdir(parents=True)
    segment_irc.write_text("1\nenergy=-1.0\nH 0 0 0\n", encoding="utf-8")
    aggregate = app["_aggregate_irc_trajectory"](
        [str(irc_plot), str(segment_irc)], str(tmp_path))
    app["_ENERGY"]["aggregate_irc"] = aggregate
    options, _ = app["_energy_diagram_options"](
        [str(irc_plot), str(segment_irc)],
        [str(segment_irc), aggregate], str(tmp_path), "all")
    assert ("IRC", "energy:irc") in options
    assert app["_ENERGY"]["views"]["energy:irc"]["path"] == aggregate
    app["S"].update(
        _last_subcmd="all",
        _last_out_dir=str(tmp_path),
        _last_files=[str(irc_plot), str(segment_irc)],
        _last_manifest={"status": "success", "exit_code": 0},
    )
    app["_results"](str(tmp_path))
    assert app["_TRAJ"]["path"] == aggregate
    assert 'class="rxenergy-frame"' in html.unescape(app["plot_out"].value)
    assert "data:image/png;base64," not in app["plot_out"].value

    segment_irc_2 = tmp_path / "segments" / "seg_02" / "irc" / "finished_irc_trj.xyz"
    segment_irc_2.parent.mkdir(parents=True)
    segment_irc_2.write_text("1\nenergy=-0.9\nH 0 0 0\n", encoding="utf-8")
    aggregate = app["_aggregate_irc_trajectory"](
        [str(segment_irc), str(segment_irc_2)], str(tmp_path))
    app["_ENERGY"]["aggregate_irc"] = aggregate
    options, _ = app["_energy_diagram_options"](
        [str(segment_irc), str(segment_irc_2)],
        [str(segment_irc), str(segment_irc_2), aggregate], str(tmp_path), "all")
    assert options.count(("IRC", "energy:irc")) == 1

    standalone_irc = tmp_path / "finished_irc_trj.xyz"
    forward_irc = tmp_path / "forward_irc_trj.xyz"
    backward_irc = tmp_path / "backward_irc_trj.xyz"
    for path in (standalone_irc, forward_irc, backward_irc):
        path.write_text("1\nIRC\nH 0 0 0\n", encoding="utf-8")
    app["_energy_diagram_options"](
        [str(forward_irc), str(backward_irc), str(standalone_irc)],
        [str(forward_irc), str(backward_irc), str(standalone_irc)],
        str(tmp_path), "irc",
    )
    assert app["_ENERGY"]["views"]["energy:irc"]["path"] == str(standalone_irc)

    source = _notebook()["cells"][2]["source"]
    assert "_frame_play_jslink = W.jslink" in source
    assert "_frame_play_link" not in source
    assert "W.link((frame_play, 'value'), (frame_slider, 'value'))" not in source
    assert "if _UPLOAD_MODE != 'colab': _send_trajectory_frame(i)" in source
    assert "const modelUiPoll=setInterval(readModelUiFrame,250)" not in source
    assert "if(!event.isTrusted||!event.target.closest('button'))return" in source
    assert "if(!ready||!viewer||applyingFrame" in source
    assert "if getattr(frame_play, '_playing', False):\n        return" in source
    assert "intervalNode.dataset.rxPlaybackInterval" in source
    assert "var interval=50;" not in source
    assert "shape:'linear'" in source
    assert "shape:'spline'" not in source
    assert "vibration?850:50" not in source
    app["frame_slider"].disabled = False
    app["frame_slider"].min = 0
    app["frame_slider"].max = 2
    app["frame_slider"].value = 0
    app["_TRAJ"]["generation"] = 7
    app["frame_play"]._playing = True
    app["_set_frame_from_browser"](7, 2)
    assert app["frame_slider"].value == 0
    app["frame_play"]._playing = False

    # Generic path-opt output names still inherit reaction-path semantics and
    # therefore the same 50 ms playback speed as explicitly named MEP files.
    path_opt = tmp_path / "final_geometries_trj.xyz"
    path_opt.write_text(
        "1\n-1.0\nH 0 0 0\n1\n-0.9\nH 0.1 0 0\n", encoding="utf-8")
    app["S"]["_last_subcmd"] = "path-opt"
    app["_load_trajectory"](str(path_opt), str(tmp_path))
    assert app["_TRAJ"]["semantics"]["title"].startswith("Reaction-path")
    assert app["frame_play"].interval == 50


def test_results_aggregate_segment_data_without_segment_controls(
    monkeypatch, tmp_path: Path,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)

    structures = {}
    for segment in (1, 2):
        structures[segment] = {}
        for state in ("R", "TS", "P"):
            path = tmp_path / "segments" / f"seg_{segment:02d}" / f"{state}.xyz"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"1\nsegment {segment} {state}\nH {segment} 0 0\n", encoding="utf-8")
            structures[segment][state] = str(path)
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps({
        "energy_diagrams": [{
            "name": "energy_diagram_DFT_all",
            "labels": ["R", "TS1", "P1", "R2", "TS2", "P"],
            "energies_kcal": [0.0, 6.0, -1.0, -1.0, 7.0, -2.0],
            "structures": [
                structures[1]["R"], structures[1]["TS"], structures[1]["P"],
                structures[2]["R"], structures[2]["TS"], structures[2]["P"],
            ],
        }],
        "post_segments": [
            {
                "index": segment,
                "dft": {
                    "labels": ["R", f"TS{segment}", "P"],
                    "energies_kcal": [0.0, 5.0 + segment, -1.0],
                    "structures": structures[segment],
                },
            }
            for segment in (1, 2)
        ]
    }), encoding="utf-8")
    current = [str(summary), *[path for group in structures.values() for path in group.values()]]

    options, _ = app["_energy_diagram_options"](current, [], str(tmp_path), "all")
    assert options.count(("DFT//MLIP ΔE", "energy:dft_e")) == 1
    assert all("Segment" not in label for label, _token in options)
    assert app["_ENERGY"]["views"]["energy:dft_e"]["structures"] == [
        structures[1]["R"], structures[1]["TS"], structures[1]["P"],
        structures[2]["R"], structures[2]["TS"], structures[2]["P"],
    ]
    nested_result = tmp_path / "segments" / "seg_02" / "result.json"
    nested_result.write_text(json.dumps({
        "energy_diagrams": [{
            "name": "DFT", "labels": ["R", "TS2", "P"],
            "energies_kcal": [0.0, 7.0, -1.0], "structures": structures[2],
        }]
    }), encoding="utf-8")
    nested_options, _ = app["_energy_diagram_options"](
        [str(nested_result), *structures[2].values()], [], str(tmp_path), "all")
    assert not nested_options

    irc_dir = tmp_path / "segments" / "seg_02" / "irc"
    irc_dir.mkdir(parents=True, exist_ok=True)
    finished = irc_dir / "seg02_finished_irc_trj.xyz"
    forward = irc_dir / "seg02_forward_irc_trj.xyz"
    for path in (finished, forward):
        path.write_text("1\nIRC\nH 0 0 0\n", encoding="utf-8")
    irc_plot = irc_dir / "irc_plot.png"
    irc_plot.write_bytes(b"\x89PNG\r\n\x1a\n")
    aggregate = app["_aggregate_irc_trajectory"](
        [str(irc_plot), str(finished), str(forward)], str(tmp_path))
    app["_ENERGY"]["aggregate_irc"] = aggregate
    irc_options, _ = app["_energy_diagram_options"](
        [str(irc_plot), str(finished), str(forward)],
        [str(finished), str(forward), aggregate], str(tmp_path), "all")
    assert irc_options == [("IRC", "energy:irc")]
    labels = [label for label, _path in app["_result_view_candidates"](
        [str(finished), str(forward), aggregate], str(tmp_path), "all")]
    assert labels == ["IRC trajectory"]

    (irc_dir / "result.json").write_text(json.dumps({
        "n_frames_forward": 0, "n_frames_backward": 1, "n_frames_total": 2,
        "forward_requested": False, "backward_requested": True,
        "forward_converged": None, "backward_converged": True,
        "scientific_status": "success",
    }), encoding="utf-8")
    semantics = app["_trajectory_semantics"]("irc", str(finished), n_frames=2)
    assert semantics["title"] == "Backward IRC trajectory"
    assert semantics["start"] == "input TS"
    assert "forward" not in semantics["trajectory_status"]
    assert app["_stationary"]([0.0], {"start": "start", "end": "end"}) == [(0, "only frame")]

    range_dir = tmp_path / "range_check"
    range_dir.mkdir()
    range_path = range_dir / "mep_trj.xyz"
    range_path.write_text("1\na\nH 0 0 0\n", encoding="utf-8")
    (range_dir / "summary.json").write_text(json.dumps({
        "segments": [
            {"index": 1, "kind": "seg", "frame_ranges": [[0, 2]]},
            {"index": 2, "kind": "seg", "frame_ranges": [[1, 3]]},
        ]
    }), encoding="utf-8")
    range_metadata = app["_trajectory_segment_ranges"](str(range_path), n_frames=3)
    assert "gap or overlap" in range_metadata["metadata_warning"]

    image_only = tmp_path / "preview.png"
    image_only.write_bytes(b"\x89PNG\r\n\x1a\n")
    app["_TRAJ"].update(frames=["old"], path="old.xyz")
    app["S"].update(
        _last_out_dir=str(tmp_path), _last_subcmd="all", _last_files=[str(image_only)],
        _last_manifest={"status": "success", "exit_code": 0},
    )
    app["_results"](str(tmp_path))
    assert app["_TRAJ"]["frames"] == []
    assert app["_TRAJ"]["path"] is None
    assert app["trajectory_box"].layout.display == "none"

    recovered_root = tmp_path / "recovered"
    recovered_irc = recovered_root / "segments" / "seg_02" / "irc" / "finished_irc_trj.xyz"
    recovered_irc.parent.mkdir(parents=True)
    recovered_irc.write_text("1\nIRC\nH 0 0 0\n", encoding="utf-8")
    stale = recovered_root / "unclaimed_old.xyz"
    stale.write_text("1\nstale\nH 0 0 0\n", encoding="utf-8")
    (recovered_root / "summary.json").write_text(json.dumps({
        "mlmm_toolkit_version": "0.3.3", "status": "success",
        "current_output_paths": ["summary.json"],
        "key_output_files": {
            "seg_02": {"files": ["irc/finished_irc_trj.xyz", "missing.xyz"]}
        },
    }), encoding="utf-8")
    assert app["_recover_existing_results"](str(recovered_root), replace=True)
    assert str(recovered_irc.resolve()) in app["S"]["_last_files"]
    assert str(stale.resolve()) not in app["S"]["_last_files"]
    assert any(path.endswith("missing.xyz") for path in
               app["S"]["_last_manifest"]["rejected_claims"])


def test_run_log_does_not_probe_or_announce_model_weight_cache() -> None:
    app = _notebook()["cells"][2]["source"]
    assert "def _model_weights_cached" not in app
    assert "def _model_download_line" not in app
    assert "Downloading %s model weights" not in app
    assert "transcript = [command_line]" in app


def test_results_playback_preview_reset_and_cli_output_defaults(
    monkeypatch, tmp_path: Path,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)

    # Image previews are stateful HTML, not background Output displays that can
    # escape from the Results card into the outer Colab cell.
    png = tmp_path / "preview.png"
    png.write_bytes(b"\x89PNG\r\n\x1a\npreview")
    app["S"].update(_last_out_dir=str(tmp_path), _last_files=[str(png)])
    app["artifact_choice"].options = [("image", str(png))]
    app["artifact_choice"].value = str(png)
    app["artifact_fold"].layout.display = ""
    app["artifact_fold"]._rx_set_open(True)
    app["_render_artifact"]()
    assert "data:image/png;base64," in app["artifact_out"].value
    assert 'class="rxartifact-image"' in app["artifact_out"].value
    notebook_source = NOTEBOOK.read_text(encoding="utf-8")
    assert ".rxartifact-image, .rxpath-panel img.rxartifact-image" in notebook_source
    assert "width:50% !important; max-width:50% !important" in notebook_source
    assert "overflow-x:hidden!important" in notebook_source

    standalone = tmp_path / "surface.html"
    standalone.write_text(
        "<html><head></head><body><h1>Energy Landscape with 2D PES Scan</h1>"
        "<div class='js-plotly-plot'></div></body></html>",
        encoding="utf-8",
    )
    preview = html.unescape(app["_artifact_preview_html"](str(standalone), str(tmp_path)))
    assert "rx-artifact-fit" in preview
    assert "Plotly.Plots.resize" in preview
    assert "aspect-ratio:16/10" in preview
    assert "min-height:320px;max-height:520px" in preview
    assert "ResizeObserver" in preview
    assert "Energy Landscape with 2D PES Scan" in preview
    assert 'srcdoc="' in preview
    assert "MutationObserver" not in preview

    # Important trajectory outputs remain inspectable in Generated file preview.
    path_profile = tmp_path / "final_geometries_trj.xyz"
    path_profile.write_text(
        "2\nR energy=-1.0\nH 0 0 0\nH 0.7 0 0\n"
        "2\nP energy=-0.9\nH 0 0 0\nH 0.8 0 0\n",
        encoding="utf-8",
    )
    app["S"]["_last_subcmd"] = "path-opt"
    profile_preview = html.unescape(
        app["_artifact_preview_html"](str(path_profile), str(tmp_path))
    )
    assert 'class="rxenergy-frame"' in profile_preview
    assert "Interactive energy profile" in profile_preview

    mode_path = tmp_path / "mode_0001_-100.00cm-1_trj.xyz"
    mode_path.write_text(path_profile.read_text(encoding="utf-8"), encoding="utf-8")
    mode_preview = html.unescape(
        app["_artifact_preview_html"](str(mode_path), str(tmp_path))
    )
    assert "2 frames" in mode_preview
    assert "vibrational mode" in mode_preview
    assert 'data-rx-channel="artifact"' in mode_preview

    # Starting another workflow replaces every previous Results fragment.
    app["res_out"].value = "old all summary"
    app["_begin_results_attempt"]("./result_scan/", "scan")
    assert app["res_out"].value == ""
    assert app["artifact_out"].value == ""

    defaults = {
        "all": "./result_all/", "opt": "./result_opt/",
        "tsopt": "./result_tsopt/", "freq": "./result_freq/",
        "irc": "./result_irc/", "scan": "./result_scan/",
        "scan2d": "./result_scan2d/", "scan3d": "./result_scan3d/",
        "path-opt": "./result_path_opt/", "path-search": "./result_path_search/",
        "sp": "./result_sp/", "dft": "./result_dft/",
    }
    # Exercise the same dropdown transition used by the browser.  Checking the
    # table helper directly misses observer-order bugs that can leave the output
    # directory one workflow behind.
    for subcommand, expected in defaults.items():
        command = app["_advanced_command"](subcommand)
        cli_out = next(param.default for param in command.params
                       if param.name == "out_dir")
        assert os.path.normpath(str(cli_out)) == os.path.normpath(expected)
        if subcommand not in app["SUBS"]:
            continue
        app["set_subcmd"](subcommand)
        assert app["dd_subcmd"].value == subcommand
        assert app["S"]["out_dir"] == expected
        assert app["w_out"].value == expected


def test_notebook_default_controls_defer_to_live_cli(
    monkeypatch, tmp_path: Path,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)

    # Every generic Advanced control starts in the omission state; Click/YAML
    # remains the source of truth until the user chooses an override.
    for subcommand in app["COMPUTE"]:
        assert app["_advanced_argv"](subcommand) == []

    # The visible shared threshold control must name each command's effective
    # Click default rather than retaining one global label.
    expected_thresholds = {
        "all": "default: gau", "opt": "default: gau",
        "path-opt": "default: gau", "path-search": "default: gau",
        "scan": "default: gau", "tsopt": "default: baker",
        "scan2d": "default: baker", "scan3d": "default: baker",
    }
    for subcommand, expected in expected_thresholds.items():
        app["set_subcmd"](subcommand)
        assert app["adv_thresh"].options[0][0] == expected
        assert app["adv_thresh"].value == "(default)"

    # freq inherits the CLI/YAML export policy: ten value-sorted modes.
    freq_params = {param.name: param for param in
                   app["_advanced_command"]("freq").params}
    assert freq_params["max_write"].default == 10
    assert freq_params["sort"].default == "value"
    all_params = {param.name: param for param in
                  app["_advanced_command"]("all").params}
    assert app["_cli_default_label"](all_params["dft_func_basis"]) == (
        "default: wb97m-v/def2-tzvpd"
    )


def test_freq_results_expose_ten_written_modes_in_top_selector(
    monkeypatch, tmp_path: Path,
) -> None:
    app, _ = _execute_app(monkeypatch, tmp_path)
    trajectory = (
        "2\nmode frame 1\nH 0 0 0\nH 0.7 0 0\n"
        "2\nmode frame 2\nH 0.1 0 0\nH 0.8 0 0\n"
    )
    frequencies = [-500.0, -120.0, 25.0, 80.0, 150.0, 300.0, 450.0, 700.0, 1000.0, 1500.0]
    modes = []
    for index, frequency in enumerate(frequencies):
        path = tmp_path / f"mode_{index:04d}_{frequency:+.2f}cm-1_trj.xyz"
        path.write_text(trajectory, encoding="utf-8")
        modes.append(path)
    reference = tmp_path / "final_geometry.xyz"
    reference.write_text(trajectory.split("2\nmode frame 2", 1)[0], encoding="utf-8")
    app["S"].update(
        _last_out_dir=str(tmp_path), _last_subcmd="freq",
        _last_files=[str(reference), *[str(path) for path in modes]],
        _last_manifest={"status": "success", "exit_code": 0},
    )
    app["_results"](str(tmp_path))

    labels = [label for label, _ in app["traj_choice"].options]
    assert len(labels) == 11
    assert labels[0] == "Frequency reference structure"
    assert labels[1] == "Mode 1 of 10 · −500.00 cm⁻¹ · imaginary"
    assert labels[-1] == "Mode 10 of 10 · 1500.00 cm⁻¹"
    assert app["traj_choice"].description == "Mode"
    assert app["result_selector_row"].layout.display == ""
    assert "data-rx-result-generation" in app["result_selector_meta"].value
    assert app["trajectory_head"].children[1] is app["result_loading"]
    assert app["trajectory_head"].children[2] is app["result_selector_row"]

    app["traj_choice"].value = str(modes[-1])
    assert app["_TRAJ"]["path"] == str(modes[-1])
    assert len(app["_TRAJ"]["frames"]) == 2


def test_colab_setup_exposes_linked_nonnegative_radius_and_selected_resn() -> None:
    app = _notebook()["cells"][2]["source"]
    setup = _notebook()["cells"][1]["source"]
    assert 'backend = "uma"' in setup
    assert 'mlmm_toolkit_version = "v0.3.3"' in setup
    assert 'prep_radius = W.BoundedFloatText(' in app
    assert 'value=2.6, min=0.0, max=1000000.0, step=0.5' in app
    assert 'adv_radius = W.BoundedFloatText(value=2.6, min=0.0, max=1000000.0, step=0.5' in app
    assert "widget.add_class('rxhalf-step')" in app
    assert "target.add_class('rxhalf-step')" in app
    assert "target_edit.add_class('rxhalf-step')" in app
    assert "document.addEventListener('keydown'" in app
    assert "document.addEventListener('pointerdown'" in app
    assert "input.setAttribute('value',input.value)" in app
    assert "input.removeAttribute('min')" in app
    assert "value=Math.max(minimum,value)" in app
    assert "value=Math.min(maximum,value)" in app
    assert "input.setAttribute('min',input.dataset.rxHalfMin)" in app
    assert "input.setAttribute('max',input.dataset.rxHalfMax)" in app
    assert "document.addEventListener('input',settle,true)" in app
    assert "document.addEventListener('change',settle,true)" in app
    assert "wireHalfStepSpinners();" in app
    assert "if radius_applies: cmd += ['-r', str(r)]" in app
    assert "selected_resn = W.Text" in app
    assert "cmd += ['--selected-resn', selected]" in app
    assert "b_done_selected_resn = b_pick_selected_resn" in app
    assert "Done picking" in app
    assert "Add or remove residue IDs in --selected-resn by clicking the sequence or 3D viewer." in app
    assert "those clicks update this field" in app
    assert "At <code>-r 0</code>" in app
    assert "without radius-based expansion" in app
    assert "0 is allowed" not in app
    assert "description=('%s charge (auto)' if auto else 'set %s charge') % rn" in app
    assert "description='', style={'description_width': '0px'}" in app
    assert "placeholder='123, A:456'" in app
    assert "b_clear_selected_resn" not in app
    assert "center_actions = W.HBox([b_pick_selected_resn, b_clear_center]" in app
    assert "center_actions.add_class('rxcenter-actions')" in app
    assert "b_pick_selected_resn.add_class('rxpickresn-button')" in app
    assert "b_clear_center.add_class('rxclear-center')" in app
    assert ".rxcenter-actions > .rxpickresn-button" in app
    assert "padding-inline:6px !important" in app
    assert "flex='1 1 174px'" in app
    assert "flex='0 0 68px'" in app
    assert "optional_stage_note = W.HTML(" in app
    assert app.count("DFT is hidden (rerun Installation with install_dft enabled)") == 1
    assert "]), optional_stage_note])" in app
    assert "if not DFT_READY: _missing_features.append" not in app
    assert "tooltip='Clear the center, ligand charges, and force-included residues.'" in app
    assert "Click a residue in the sequence panel or any atom in the 3D viewer" in app
    assert "if panel is freeze_panel and S.get('mode') == 'small':" in app
    assert "system_charge_panel.add_class('rxactive-card')" in app
    assert "atom_info = _info_control(" in app
    assert "Click a selected atom again, or its × chip, to remove it." in app
    assert "atom_entry = W.Text(" not in app
    assert "'selected_resn': {'all', 'extract'}" in app
    assert "center_widget = W.SelectMultiple" in app
    assert "rows=min(5, max(2, len(het)))" in app
    assert "description='Save settings'" in app
    assert "Settings only; input files are not embedded." not in app


def test_scan_editor_uses_inspector_width_and_content_sized_stage_cards() -> None:
    app = _notebook()["cells"][2]["source"]
    assert "container-type:inline-size" in app
    assert "@container (max-width: 400px)" in app
    assert ".rxscan-pair-controls > :last-child { grid-column:1 / -1; }" in app
    assert ".rxscan-add-actions { grid-template-columns:minmax(0,.82fr) minmax(0,1.18fr) !important; gap:4px !important; }" in app
    assert "@container (max-width: 280px)" in app
    assert ".rxscan-pair-controls > :last-child { grid-column:1 / -1; }" in app
    assert ".rxscan-add-actions { grid-template-columns:1fr !important; }" in app
    assert ".rxscan-bond > * { flex:0 1 auto !important;" in app
    assert "label.layout = W.Layout(width='auto', flex='1 1 90px', min_width='70px')" in app
    assert "coordinate_children = [label, target_edit]" in app
    assert "if len(stage) > 1:" in app
    assert "bond_remove.add_class('rxscan-coordinate-remove')" in app
    assert "Move coordinate up" not in app
    assert "label.layout = W.Layout(flex='1 1 230px'" not in app
    assert "description='① Pick atom A'" in app
    assert "description='② Pick atom B'" in app
    assert "description='Add to stage'" in app
    assert "description='Add sequential stage'" in app
    assert "tooltip='Validate the command shown above without starting a calculation.'" in app
    assert "tooltip='Execute the command shown above.'" in app
    assert "tooltip='Rebuild the command from the current settings.'" in app
    assert "tooltip='Stop the current validation or run.'" in app
    assert "tooltip='Download the current settings as session.json.'" in app
    assert "_FRIENDLY_OPTION_HELP = {" in app
    assert "'thermo': 'Run PHVA vibrational and QRRHO thermochemistry analyses" in app
    assert 'display(HTML("""<style>\nhtml, body {' not in app


def test_cli_flag_introspection_does_not_leak_banner_output() -> None:
    app = _notebook()["cells"][2]["source"]
    assert app.count("with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):") >= 3
    assert "_ADV_COMMAND_CACHE[sub] = PRODUCT_CLI.get_command" in app
    assert "_ADV_BOOL_CACHE[sub] = PRODUCT_CLI._resolve_bool_options" in app
    assert "with command.make_context(" in app


def test_mep_commands_use_algorithm_specific_cycle_flags() -> None:
    import click

    from mlmm.cli import cli as root_cli

    app = _notebook()["cells"][2]["source"]
    assert "adv_maxcyc" not in app
    assert "'--max-cycles'" not in app
    assert "--max-cycles (0=None)" not in app

    context = click.Context(root_cli)

    def opts_of(sub: str) -> set[str]:
        command = root_cli.get_command(context, sub)
        return {opt for param in command.params for opt in param.opts}

    for sub in ("all", "path-opt", "path-search"):
        available = opts_of(sub)
        assert "--max-cycles-gsm" in available, sub
        assert "--max-cycles-dmf" in available, sub
        assert "--max-cycles" not in available, sub
    for sub in ("opt", "tsopt", "irc"):
        assert "--max-cycles" in opts_of(sub), sub


def test_plot_export_installation_does_not_launch_a_render_probe() -> None:
    setup = _notebook()["cells"][1]["source"]
    assert "subprocess.run(['plotly_get_chrome','-y']" in setup
    assert "['apt-cache','show','libasound2t64']" in setup
    for library in (
        "libnss3", "libatk-bridge2.0-0", "libcups2", "libxcomposite1",
        "libxdamage1", "libxfixes3", "libxrandr2", "libgbm1",
        "libxkbcommon0", "libpango-1.0-0", "libcairo2",
    ):
        assert library in setup
    assert "Chrome installation for Plotly export failed" in setup
    assert "write_image" in setup
    assert "timeout=90" in setup
    assert "_plot_probe_code" in setup
    assert "PNG export verified" in setup
    assert "def _run_installer(command, label, check=True):" in setup
    assert "completed = subprocess.run(command)" in setup
    assert "Still installing" not in setup
    assert "worker.join(timeout=max(1, interval))" not in setup
    assert "_ChoreographerChromium.find_browser" not in setup


def test_colab_setup_uses_the_hashed_published_pdbfixer_release() -> None:
    setup = _notebook()["cells"][1]["source"]
    tree = ast.parse(setup)
    pip_calls = [
        node for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "pip"
    ]
    pdbfixer_call = next(
        node for node in pip_calls
        if any(
            isinstance(arg, ast.Constant)
            and "pdbfixer-1.12.0.tar.gz" in str(arg.value)
            for arg in node.args
        )
    )
    args = [ast.literal_eval(arg) for arg in pdbfixer_call.args]
    assert args[:3] == ["--timeout=60", "--retries=8", "openmm[cuda12]"]
    assert args[-1] == "ninja"
    assert args[3].startswith("https://files.pythonhosted.org/packages/")
    assert args[3].endswith(
        "pdbfixer-1.12.0.tar.gz#sha256="
        "82e54074b33eb00270e9a439c8633fcf82cd3ab7a36561b01425a260c0fc4143"
    )
    assert "github.com/openmm/pdbfixer/archive" not in setup
    assert setup.count("--timeout=60") == 1
    assert setup.count("--retries=8") == 1
    assert "subprocess.run(command, timeout=" not in setup
