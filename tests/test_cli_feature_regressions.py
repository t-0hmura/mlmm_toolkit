"""CLI regressions for recently added/removed options."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest
from click.testing import CliRunner

from mlmm.cli import cli as root_cli

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm CLI requires Python >= 3.11",
)


def test_freeze_links_removed_from_help_outputs() -> None:
    runner = CliRunner()
    for command_name in ("opt", "tsopt", "freq", "irc"):
        short = runner.invoke(root_cli, [command_name, "--help"])
        assert short.exit_code == 0, short.output
        assert "--freeze-links" not in short.output

        advanced = runner.invoke(root_cli, [command_name, "--help-advanced"])
        assert advanced.exit_code == 0, advanced.output
        assert "--freeze-links" not in advanced.output


def test_freeze_links_now_fails_as_unknown_option() -> None:
    runner = CliRunner()
    for command_name in ("opt", "tsopt", "freq", "irc"):
        result = runner.invoke(root_cli, [command_name, "--freeze-links"])
        assert result.exit_code != 0
        assert "No such option" in result.output


def test_path_search_help_shows_refine_mode() -> None:
    runner = CliRunner()
    result = runner.invoke(root_cli, ["path-search", "--help"])
    assert result.exit_code == 0, result.output
    assert "--refine-mode" in result.output


def test_path_search_dry_run_uses_prepared_layer_source(tmp_path: Path) -> None:
    repo = Path(__file__).resolve().parents[1]
    smoke = repo / "tests" / "smoke"
    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "path-search",
            "-i",
            str(smoke / "r_complex_layered.pdb"),
            "-i",
            str(smoke / "p_complex_layered.pdb"),
            "--parm",
            str(smoke / "p_complex.parm7"),
            "-q",
            "-1",
            "-m",
            "1",
            "--dry-run",
            "--out-dir",
            str(tmp_path / "path-search"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "[dry-run] Validation complete. Path search execution was skipped." in result.output


def test_path_opt_help_shows_fix_ends() -> None:
    runner = CliRunner()
    result = runner.invoke(root_cli, ["path-opt", "--help"])
    assert result.exit_code == 0, result.output
    assert "--fix-ends" in result.output


def test_irc_help_names_the_ml_region_charge() -> None:
    result = CliRunner().invoke(root_cli, ["irc", "--help"])

    assert result.exit_code == 0, result.output
    assert "Net charge of the ML region/model system" in result.output
    assert "overrides calc.model_charge from YAML" in result.output


def test_scan2d_declares_scan_lists_required() -> None:
    """Reject an incomplete scan2d invocation during Click parsing."""
    command = root_cli.get_command(None, "scan2d")
    param = next(p for p in command.params if p.name == "scan_list_raw")
    assert param.required is True


@pytest.mark.parametrize("repeated", [False, True])
def test_scan_accepts_grouped_and_repeated_stages(
    tmp_path: Path,
    repeated: bool,
) -> None:
    repo = Path(__file__).resolve().parents[1]
    smoke = repo / "tests" / "smoke"
    first = "[(1,2,1.8)]"
    second = "[(1,2,2.0)]"
    stages = (
        ["--scan-lists", first, "--scan-lists", second]
        if repeated
        else ["--scan-lists", first, second]
    )
    result = CliRunner().invoke(
        root_cli,
        [
            "scan", "-i", str(smoke / "r_complex_layered.pdb"),
            "--parm", str(smoke / "p_complex.parm7"),
            "-q", "-1", "-m", "1", *stages,
            "--dry-run", "--out-dir", str(tmp_path / f"scan-{repeated}"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Received 2 stage(s)" in result.output


@pytest.mark.parametrize("repeated", [False, True])
def test_all_accepts_grouped_and_repeated_inputs(
    tmp_path: Path,
    repeated: bool,
) -> None:
    repo = Path(__file__).resolve().parents[1]
    smoke = repo / "tests" / "smoke"
    reactant = str(smoke / "r_complex_layered.pdb")
    product = str(smoke / "p_complex_layered.pdb")
    inputs = (
        ["-i", reactant, "-i", product]
        if repeated
        else ["-i", reactant, product]
    )
    result = CliRunner().invoke(
        root_cli,
        [
            "all", *inputs, "--parm", str(smoke / "p_complex.parm7"),
            "-q", "-1", "-m", "1", "--detect-layer", "--dry-run",
            "--out-dir", str(tmp_path / f"all-{repeated}"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "inputs=2" in result.output


def test_all_rejects_trailing_orphan_equal_to_consumed_input(tmp_path: Path) -> None:
    repo = Path(__file__).resolve().parents[1]
    smoke = repo / "tests" / "smoke"
    reactant = str(smoke / "r_complex_layered.pdb")
    result = CliRunner().invoke(
        root_cli,
        [
            "all", "-i", reactant,
            "--parm", str(smoke / "p_complex.parm7"),
            "-q", "-1", "-m", "1",
            "--scan-lists", "[(1,2,1.8)]", "--detect-layer", "--dry-run",
            "--out-dir", str(tmp_path / "duplicate-orphan"), reactant,
        ],
    )

    assert result.exit_code == 2
    assert f"Unexpected extra argument: {reactant}" in result.output


@pytest.mark.parametrize("repeated", [False, True])
def test_all_accepts_grouped_and_repeated_scan_stages(
    tmp_path: Path,
    repeated: bool,
) -> None:
    repo = Path(__file__).resolve().parents[1]
    smoke = repo / "tests" / "smoke"
    first = "[(1,2,1.8)]"
    second = "[(1,2,2.0)]"
    stages = (
        ["--scan-lists", first, "--scan-lists", second]
        if repeated
        else ["--scan-lists", first, second]
    )
    result = CliRunner().invoke(
        root_cli,
        [
            "all", "-i", str(smoke / "r_complex_layered.pdb"),
            "--parm", str(smoke / "p_complex.parm7"),
            "-q", "-1", "-m", "1", *stages, "--detect-layer", "--dry-run",
            "--out-dir", str(tmp_path / f"all-scan-{repeated}"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "scan=yes" in result.output


def test_all_rejects_scan_lists_with_multiple_inputs(tmp_path: Path) -> None:
    repo = Path(__file__).resolve().parents[1]
    smoke = repo / "tests" / "smoke"
    result = CliRunner().invoke(
        root_cli,
        [
            "all",
            "-i",
            str(smoke / "r_complex_layered.pdb"),
            str(smoke / "p_complex_layered.pdb"),
            "--parm",
            str(smoke / "p_complex.parm7"),
            "-q",
            "-1",
            "-m",
            "1",
            "--scan-lists",
            "[(1,2,1.8)]",
            "--detect-layer",
            "--dry-run",
            "--out-dir",
            str(tmp_path / "invalid-mixed-mode"),
        ],
    )

    assert result.exit_code != 0
    assert "--scan-lists requires exactly one input structure" in result.output


def test_energy_diagram_rejects_unknown_options() -> None:
    result = CliRunner().invoke(
        root_cli,
        ["energy-diagram", "-i", "0", "-i", "1", "--bogus", "123"],
    )
    assert result.exit_code == 2
    assert "No such option: --bogus" in result.output

    orphan = CliRunner().invoke(
        root_cli,
        ["energy-diagram", "-i", "0", "1", "--label-y", "E", "orphan"],
    )
    assert orphan.exit_code == 2
    assert "Unexpected extra argument: orphan" in orphan.output


def test_energy_diagram_preserves_equals_attached_and_grouped_values(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    from mlmm.io import energy_diagram

    captured: dict = {}

    class _DummyFigure:
        def write_image(self, path: str, scale: int = 2) -> None:
            captured["path"] = path
            captured["scale"] = scale
            Path(path).write_text("image", encoding="utf-8")

    def _build(**kwargs):
        captured.update(kwargs)
        return _DummyFigure()

    monkeypatch.setattr(energy_diagram, "build_energy_diagram", _build)
    output = tmp_path / "mixed.png"
    result = CliRunner().invoke(
        root_cli,
        [
            "energy-diagram", "--input=0", "-i1", "2",
            "--label-x=R", "--label-x", "TS", "P", "-o", str(output),
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert captured["energies"] == [0.0, 1.0, 2.0]
    assert captured["labels"] == ["R", "TS", "P"]
    assert output.exists()


@pytest.mark.parametrize(
    ("command_name", "base_args"),
    [
        ("extract", ["-i", "r_complex_layered.pdb", "-c", "1"]),
        (
            "scan",
            [
                "-i",
                "r_complex_layered.pdb",
                "--parm",
                "p_complex.parm7",
                "-s",
                "[(0, 1, 1.5)]",
            ],
        ),
        (
            "path-search",
            [
                "-i",
                "r_complex_layered.pdb",
                "p_complex_layered.pdb",
                "--parm",
                "p_complex.parm7",
            ],
        ),
        (
            "all",
            [
                "-i",
                "r_complex_layered.pdb",
                "p_complex_layered.pdb",
                "--parm",
                "p_complex.parm7",
            ],
        ),
    ],
)
@pytest.mark.parametrize("unknown", ["--bogus-option", "-Z"])
def test_legacy_multi_value_commands_reject_unknown_options(
    command_name: str,
    base_args: list[str],
    unknown: str,
) -> None:
    """Compatibility parsing must not turn option typos into no-ops."""
    smoke = Path(__file__).resolve().parent / "smoke"
    argv = [
        command_name,
        *(str(smoke / value) if value.endswith((".pdb", ".parm7")) else value for value in base_args),
        unknown,
    ]
    result = CliRunner().invoke(root_cli, argv)

    assert result.exit_code == 2, result.output
    assert f"No such option: {unknown}" in result.output


@pytest.mark.parametrize(
    ("command_name", "base_args", "separator"),
    [
        ("extract", ["-i", "r_complex_layered.pdb", "-c", "1"], ["--out-json"]),
        (
            "scan",
            ["-i", "r_complex_layered.pdb", "--parm", "p_complex.parm7",
             "-s", "[(0, 1, 1.5)]"],
            ["--dry-run"],
        ),
        (
            "path-search",
            ["-i", "r_complex_layered.pdb", "p_complex_layered.pdb",
             "--parm", "p_complex.parm7"],
            ["--dry-run"],
        ),
        (
            "all",
            ["-i", "r_complex_layered.pdb", "p_complex_layered.pdb",
             "--parm", "p_complex.parm7"],
            ["--dry-run"],
        ),
    ],
)
def test_legacy_multi_value_commands_reject_unclaimed_bare_arguments(
    command_name: str,
    base_args: list[str],
    separator: list[str],
) -> None:
    smoke = Path(__file__).resolve().parent / "smoke"
    argv = [
        command_name,
        *(str(smoke / value) if value.endswith((".pdb", ".parm7")) else value for value in base_args),
        *separator,
        "orphan",
    ]
    result = CliRunner().invoke(root_cli, argv)

    assert result.exit_code == 2, result.output
    assert "Unexpected extra argument: orphan" in result.output


def test_path_search_accepts_repeated_inputs_and_reference_templates() -> None:
    smoke = Path(__file__).resolve().parent / "smoke"
    reactant = smoke / "r_complex_layered.pdb"
    product = smoke / "p_complex_layered.pdb"
    parm = smoke / "p_complex.parm7"
    result = CliRunner().invoke(
        root_cli,
        [
            "path-search",
            "-i", str(reactant),
            "-i", str(product),
            "--parm", str(parm),
            "--ref-pdb", str(reactant),
            "--ref-pdb", str(product),
            "-q", "-1",
            "-m", "1",
            "--dry-run",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Validation complete" in result.output


def test_all_exposes_irc_step_size_and_tsopt_exposes_charge_guard() -> None:
    """Keep the composite IRC controls and the shared parity escape hatch aligned."""
    runner = CliRunner()

    all_help = runner.invoke(root_cli, ["all", "--help-advanced"])
    assert all_help.exit_code == 0, all_help.output
    assert "--irc-step-size" in all_help.output
    assert "--mep-mode" in all_help.output
    assert "--dmf-backend" in all_help.output

    basic_all_help = runner.invoke(root_cli, ["all", "--help"])
    assert basic_all_help.exit_code == 0, basic_all_help.output
    assert "--mep-mode" not in basic_all_help.output
    assert "--dmf-backend" not in basic_all_help.output

    tsopt_help = runner.invoke(root_cli, ["tsopt", "--help-advanced"])
    assert tsopt_help.exit_code == 0, tsopt_help.output
    assert "--allow-charge-mult-mismatch" in tsopt_help.output


def test_all_forwards_irc_step_size_to_child(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The all-level override must reach the child as IRC's ``--step-size``."""
    from mlmm.workflows import all as all_workflow

    captured: list[str] = []

    class _StopHere(RuntimeError):
        pass

    def _capture(_name, _command, args, **_kwargs):
        captured.extend(args)
        raise _StopHere

    monkeypatch.setattr(all_workflow, "_run_cli_main", _capture)
    template = all_workflow._ResolvedCalculatorTemplate.from_mapping({})

    with pytest.raises(_StopHere):
        all_workflow._irc_and_match(
            seg_idx=1,
            seg_dir=tmp_path,
            ref_pdb_for_seg=tmp_path / "ts.pdb",
            seg_pocket_pdb=tmp_path / "model.pdb",
            g_ts=object(),
            q_int=0,
            spin=1,
            resolved_calc_template=template,
            real_parm7=tmp_path / "system.parm7",
            model_pdb=tmp_path / "model.pdb",
            irc_step_size=0.05,
        )

    idx = captured.index("--step-size")
    assert captured[idx + 1] == "0.05"


def test_scan3d_csv_mode_runs_without_scan_inputs(tmp_path: Path) -> None:
    csv_path = tmp_path / "surface.csv"
    out_dir = tmp_path / "scan3d_out"

    rows = []
    for i in (0.0, 1.0):
        for j in (0.0, 1.0):
            for k in (0.0, 1.0):
                rows.append(
                    {
                        "d1_A": i,
                        "d2_A": j,
                        "d3_A": k,
                        "energy_kcal": i + j + k,
                        "d1_label": "d1",
                        "d2_label": "d2",
                        "d3_label": "d3",
                    }
                )
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "scan3d",
            "--csv",
            str(csv_path),
            "--out-dir",
            str(out_dir),
            "--out-json",
        ],
    )
    assert result.exit_code == 0, result.output
    assert (out_dir / "scan3d_density.html").exists()
    payload = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    assert payload["mlip_backend"] is None
    assert payload["mlip_model"] is None
    assert "backend" not in payload


@pytest.mark.parametrize("extra_args", [["-q", "0"], ["-m", "2"]])
def test_trj2fig_recompute_branch_triggers_on_charge_or_multiplicity(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    extra_args: list[str],
) -> None:
    from mlmm.io import trj2fig as trj2fig_mod

    xyz_path = tmp_path / "traj.xyz"
    xyz_path.write_text(
        "1\ncomment\nH 0.0 0.0 0.0\n1\ncomment\nH 0.0 0.0 0.1\n",
        encoding="utf-8",
    )
    out_csv = tmp_path / "energy.csv"
    called = {"value": False}

    def _fake_recompute(
        _traj: Path,
        _charge: int | None,
        _mult: int | None,
        **_kwargs,
    ):
        called["value"] = True
        return [0.0, 0.001]

    monkeypatch.setattr(trj2fig_mod, "recompute_energies", _fake_recompute)

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        ["trj2fig", "-i", str(xyz_path), "-o", str(out_csv), *extra_args],
    )
    assert result.exit_code == 0, result.output
    assert called["value"] is True
    assert out_csv.exists()


def test_trj2fig_json_records_selected_backend_provenance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from mlmm.io import trj2fig as trj2fig_mod

    xyz_path = tmp_path / "traj.xyz"
    xyz_path.write_text(
        "1\nno-energy-needed\nH 0.0 0.0 0.0\n",
        encoding="utf-8",
    )
    out_csv = tmp_path / "energy.csv"
    seen: dict = {}

    def _fake_recompute(path, charge, multiplicity, **kwargs):
        seen.update(
            path=path,
            charge=charge,
            multiplicity=multiplicity,
            **kwargs,
        )
        return [-0.5]

    monkeypatch.setattr(trj2fig_mod, "recompute_energies", _fake_recompute)
    result = CliRunner().invoke(
        root_cli,
        [
            "trj2fig",
            "-i",
            str(xyz_path),
            "-o",
            str(out_csv),
            "-q",
            "-1",
            "-m",
            "2",
            "-b",
            "orb",
            "--backend-model",
            "orb-test-model",
            "--precision",
            "fp64",
            "--out-json",
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    payload = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    assert payload["energy_source"] == "mlip_recomputed"
    assert payload["energy_provenance"] == ["mlip-recomputed"]
    assert payload["energy_unit"] == "hartree"
    assert payload["mlip_backend"] == "orb"
    assert payload["mlip_model"] == "orb-test-model"
    assert payload["mlip_precision"] == "fp64"
    assert payload["charge"] == -1
    assert payload["multiplicity"] == 2
    assert seen["backend"] == "orb"
    assert seen["backend_model"] == "orb-test-model"
    assert seen["precision"] == "fp64"


def test_trj2fig_comment_json_does_not_claim_calculator_provenance(
    tmp_path: Path,
) -> None:
    xyz_path = tmp_path / "traj.xyz"
    xyz_path.write_text(
        "1\n-0.500000\nH 0.0 0.0 0.0\n",
        encoding="utf-8",
    )
    out_csv = tmp_path / "energy.csv"

    result = CliRunner().invoke(
        root_cli,
        [
            "trj2fig",
            "-i",
            str(xyz_path),
            "-o",
            str(out_csv),
            "-b",
            "mace",
            "--out-json",
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    payload = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    assert payload["energy_source"] == "trajectory_comment"
    assert payload["energy_provenance"] == ["bare-assumed-Ha"]
    assert payload["energy_unit"] == "hartree"
    assert payload["mlip_backend"] is None
    assert payload["mlip_model"] is None
    assert payload["mlip_precision"] is None
    assert payload["charge"] is None
    assert payload["multiplicity"] is None


def test_trj2fig_rejects_missing_frame_zero_energy_without_outputs(
    tmp_path: Path,
) -> None:
    xyz_path = tmp_path / "missing-frame-zero.xyz"
    xyz_path.write_text(
        "1\nframe 0\nH 0.0 0.0 0.0\n"
        "1\nE=-0.500000 Ha\nH 0.0 0.0 0.1\n",
        encoding="utf-8",
    )
    out_csv = tmp_path / "must-not-exist.csv"
    result = CliRunner().invoke(
        root_cli,
        [
            "trj2fig", "-i", str(xyz_path), "-o", str(out_csv),
            "--out-json",
        ],
    )
    assert result.exit_code != 0
    assert "frame 1" in (result.output + str(result.exception))
    assert not out_csv.exists()
    assert not (tmp_path / "result.json").exists()


def test_trj2fig_reports_out_of_range_reference_without_traceback(
    tmp_path: Path,
) -> None:
    xyz_path = tmp_path / "traj.xyz"
    xyz_path.write_text(
        "1\n-0.5\nH 0 0 0\n1\n-0.4\nH 0 0 0.1\n",
        encoding="utf-8",
    )
    out_csv = tmp_path / "energy.csv"

    result = CliRunner().invoke(
        root_cli,
        ["trj2fig", "-i", str(xyz_path), "-o", str(out_csv), "-r", "9"],
    )

    assert result.exit_code != 0
    assert "Reference index 9 out of range" in result.output
    assert "Traceback" not in result.output
    assert not out_csv.exists()


def test_trj2fig_json_preserves_same_named_outputs(tmp_path: Path) -> None:
    xyz_path = tmp_path / "traj.xyz"
    xyz_path.write_text(
        "1\n-0.500000\nH 0.0 0.0 0.0\n",
        encoding="utf-8",
    )
    first = tmp_path / "plots-a" / "same.csv"
    second = tmp_path / "plots-b" / "same.csv"
    first.parent.mkdir()
    second.parent.mkdir()

    result = CliRunner().invoke(
        root_cli,
        [
            "trj2fig", "-i", str(xyz_path), "-o", str(first), str(second),
            "--out-json",
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert first.exists() and second.exists()
    payload = json.loads((first.parent / "result.json").read_text(encoding="utf-8"))
    assert payload["output_files"] == [str(first), str(second)]
    assert payload["files"] == {"same.csv": str(second)}


def test_trj2fig_module_entrypoint_executes(tmp_path: Path) -> None:
    xyz_path = tmp_path / "traj.xyz"
    xyz_path.write_text(
        "1\nenergy=-1.000000\nH 0 0 0\n"
        "1\nenergy=-0.900000\nH 0 0 0.1\n",
        encoding="utf-8",
    )
    out_csv = tmp_path / "energy.csv"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "mlmm.io.trj2fig",
            "-i",
            str(xyz_path),
            "-o",
            str(out_csv),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert out_csv.is_file()
    assert "energy_hartree" in out_csv.read_text(encoding="utf-8")


def test_coord_type_dlc_falls_back_to_cart_under_lbfgs() -> None:
    """`--coord-type dlc` is only meaningful with Hessian-based optimization;
    under `--opt-mode grad` (L-BFGS) it must fall back to Cartesian, while
    `--opt-mode hess` keeps DLC."""
    repo = Path(__file__).resolve().parents[1]
    in_pdb = repo / "examples" / "toy_system" / "p_complex_layered.pdb"
    parm = repo / "examples" / "toy_system" / "p_complex.parm7"
    if not (in_pdb.exists() and parm.exists()):
        pytest.skip("toy_system example inputs not present")
    runner = CliRunner()
    base = [
        "opt", "-i", str(in_pdb), "--parm", str(parm), "-q", "0",
        "--detect-layer", "--coord-type", "dlc", "--dry-run",
    ]
    grad = runner.invoke(root_cli, base + ["--opt-mode", "grad"])
    assert grad.exit_code == 0, grad.output
    assert "falling back to cart" in grad.output

    hess = runner.invoke(root_cli, base + ["--opt-mode", "hess"])
    assert hess.exit_code == 0, hess.output
    assert "falling back to cart" not in hess.output


@pytest.mark.parametrize(
    ("extra", "expected"),
    [([], "legacy-active"), (["--tr-projection", "constrained"], "constrained")],
)
def test_opt_tr_projection_cli_overrides_yaml(
    tmp_path: Path, extra: list[str], expected: str,
) -> None:
    repo = Path(__file__).resolve().parents[1]
    in_pdb = repo / "examples" / "toy_system" / "p_complex_layered.pdb"
    parm = repo / "examples" / "toy_system" / "p_complex.parm7"
    if not (in_pdb.exists() and parm.exists()):
        pytest.skip("toy_system example inputs not present")
    config = tmp_path / "projection.yaml"
    config.write_text(
        "geom:\n  tr_projection: legacy-active\n",
        encoding="utf-8",
    )
    result = CliRunner().invoke(
        root_cli,
        [
            "opt", "-i", str(in_pdb), "--parm", str(parm), "-q", "0",
            "--detect-layer", "--config", str(config), "--dry-run", "-v", "3",
            *extra,
        ],
    )

    assert result.exit_code == 0, result.output
    assert f"tr_projection: {expected}" in result.output


def test_verbose_is_a_per_subcommand_option() -> None:
    """`-v/--verbose LEVEL` is injected into every subcommand (including the
    compat-parsing `extract`); it is no longer a root-group option, so a
    root-placed `-v` is rejected."""
    runner = CliRunner()
    for name in (
        "opt", "tsopt", "freq", "irc", "sp", "scan",
        "path-search", "dft", "all", "extract",
    ):
        res = runner.invoke(root_cli, [name, "--help"])
        assert res.exit_code == 0, res.output
        assert "-v, --verbose" in res.output, f"{name} --help is missing -v"

    root = runner.invoke(root_cli, ["-v", "2", "opt", "--help"])
    assert root.exit_code != 0
    assert "No such option" in root.output

    # The level is an IntRange(0, 3); 0/1/2/3 are accepted (default 2) and
    # out-of-range values are rejected, for the injected commands, including the
    # parser-wrapper `extract`.
    for cmd in (["opt", "-v", "4"], ["extract", "-v", "4"]):
        bad = runner.invoke(root_cli, cmd)
        assert bad.exit_code != 0, cmd
        assert "is not in the range" in bad.output, cmd
