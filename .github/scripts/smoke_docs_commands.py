#!/usr/bin/env python3
"""Lightweight smoke tests for commands embedded in docs markdown files."""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

from click.testing import CliRunner

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_ROOT = REPO_ROOT / "docs"
TOOL_NAME = "mlmm"
CLI_MODULE = "mlmm"
DOCS_SMOKE_COMMAND_TIMEOUT_SEC = float(os.environ.get("DOCS_SMOKE_COMMAND_TIMEOUT_SEC", "120"))

sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from mlmm.cli import cli as root_cli  # noqa: E402

from docs_command_contract import (  # noqa: E402
    extract_docs_commands,
    subcommand_from_tokens as _subcommand_from_tokens,
    validate_option_names,
)


_ALL_ONLY_PATH_EXTS = {".pdb", ".xyz", ".gjf", ".yaml", ".yml", ".json"}


def _prepare_fixture_files(tmp: Path) -> dict[str, Path]:
    pdb_text = (
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  C2  LIG A   1       1.400   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )
    r_pdb = tmp / "R.pdb"
    p_pdb = tmp / "P.pdb"
    xyz = tmp / "input.xyz"
    gjf = tmp / "input.gjf"
    cfg = tmp / "config.yaml"
    parm7 = tmp / "fixture.parm7"
    out_dir = tmp / "result_all"

    r_pdb.write_text(pdb_text, encoding="utf-8")
    p_pdb.write_text(pdb_text, encoding="utf-8")
    xyz.write_text("1\n\nC 0.0 0.0 0.0\n", encoding="utf-8")
    gjf.write_text("%chk=test\n#p hf/3-21g\n\nTitle\n\n0 1\nC 0.0 0.0 0.0\n\n", encoding="utf-8")
    cfg.write_text("extract:\n  radius: 2.6\n", encoding="utf-8")
    # all --dry-run now performs the real topology atom-count/order check.
    # Build a minimal, parameterized two-atom Amber topology so staged scans
    # can name a real pair without requiring AmberTools executables. PDB bond
    # inference has no force-field BondType, so omit that inferred bond: the
    # dry-run contract under test is atom count/order, not bonded parameters.
    import parmed as pmd
    from parmed.topologyobjects import AtomType

    structure = pmd.load_file(str(r_pdb))
    atom_type = AtomType("C", 1, 12.011, 6)
    atom_type.set_lj_params(0.1, 1.7)
    for atom in structure.atoms:
        atom.atom_type = atom_type
        atom.type = "C"
    structure.bonds.clear()
    pmd.amber.AmberParm.from_structure(structure).save(str(parm7))
    out_dir.mkdir(parents=True, exist_ok=True)

    return {
        "r_pdb": r_pdb,
        "p_pdb": p_pdb,
        "xyz": xyz,
        "gjf": gjf,
        "config": cfg,
        "parm7": parm7,
        "out_dir": out_dir,
    }


def _sanitize_all_args(args: list[str], fixture: dict[str, Path]) -> list[str]:
    out: list[str] = []
    staged_scan = any(tok in {"-s", "--scan-lists"} for tok in args)
    fixture_inputs = [str(fixture["r_pdb"])]
    if not staged_scan:
        fixture_inputs.append(str(fixture["p_pdb"]))
    saw_input = False
    saw_dry_run = False
    saw_center = False
    saw_parm = False
    i = 0
    while i < len(args):
        tok = args[i]
        if tok in {"--version", "-h", "--help"}:
            i += 1
            continue
        if tok in {"-c", "--center"}:
            # Pin the center to the fixture's LIG residue: documented examples
            # may name real cofactors that the synthetic smoke PDB cannot hold,
            # and the dry-run pre-check only needs a resolvable center to test
            # that the command parses and plans.
            saw_center = True
            out.extend([tok, "LIG"])
            i += 2
            continue
        if tok in {"-l", "--ligand-charge", "-q", "--charge"}:
            # Skip explicit charge inputs tied to the example's real residues;
            # the extractor derives a charge consistent with the LIG fixture.
            i += 2
            continue
        if tok in {"-s", "--scan-lists"}:
            out.append(tok)
            i += 1
            n_values = 0
            while i < len(args) and not args[i].startswith("-"):
                out.append(f"[(1,2,{1.5 + 0.1 * n_values:.1f})]")
                n_values += 1
                i += 1
            continue
        if tok == "--parm":
            saw_parm = True
            out.extend([tok, str(fixture["parm7"])])
            i += 2
            continue
        if tok in {"-i", "--input"}:
            saw_input = True
            out.extend([tok, *fixture_inputs])
            i += 1
            while i < len(args) and not args[i].startswith("-"):
                i += 1
            continue
        if tok == "--config":
            out.extend([tok, str(fixture["config"])])
            i += 2
            continue
        if tok == "--out-dir":
            out.extend([tok, str(fixture["out_dir"])])
            i += 2
            continue
        if tok == "--dry-run":
            saw_dry_run = True
            out.append(tok)
            i += 1
            continue
        if tok == "--no-dry-run":
            saw_dry_run = True
            out.append("--dry-run")
            i += 1
            continue
        if (not tok.startswith("-")) and Path(tok).suffix.lower() in _ALL_ONLY_PATH_EXTS:
            ext = Path(tok).suffix.lower()
            if ext == ".pdb":
                out.append(str(fixture["r_pdb"]))
            elif ext == ".xyz":
                out.append(str(fixture["xyz"]))
            elif ext == ".gjf":
                out.append(str(fixture["gjf"]))
            else:
                out.append(str(fixture["config"]))
            i += 1
            continue

        out.append(tok)
        i += 1

    if not saw_input:
        out.extend(["-i", *fixture_inputs])
    if not saw_center:
        out.extend(["-c", "LIG"])
    if not saw_parm:
        out.extend(["--parm", str(fixture["parm7"])])
    if not saw_dry_run:
        out.append("--dry-run")
    if "--out-dir" not in out:
        out.extend(["--out-dir", str(fixture["out_dir"])])
    return out


def _run_help_smoke(commands: list[str]) -> None:
    runner = CliRunner()
    subcommands = sorted({_subcommand_from_tokens(shlex.split(cmd)) for cmd in commands})
    for subcmd in subcommands:
        result = runner.invoke(root_cli, [subcmd, "--help"], catch_exceptions=False)
        if result.exit_code != 0:
            raise RuntimeError(
                f"[help-smoke] failed for '{TOOL_NAME} {subcmd} --help':\n{result.output}"
            )
    print(f"[help-smoke] validated {len(subcommands)} subcommands from docs.")


def _run_all_dry_run_smoke(commands: list[str]) -> None:
    try:
        probe = subprocess.run(
            [sys.executable, "-m", CLI_MODULE, "all", "--help"],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            timeout=DOCS_SMOKE_COMMAND_TIMEOUT_SEC,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            "[dry-run-smoke] timeout while probing availability for "
            f"'{TOOL_NAME} all --help' ({DOCS_SMOKE_COMMAND_TIMEOUT_SEC:g}s)."
        ) from exc
    probe_output = f"{probe.stdout}\n{probe.stderr}"
    if "Command 'all' is unavailable" in probe_output or "Missing dependency:" in probe_output:
        raise RuntimeError(
            "[dry-run-smoke] required 'all' command is unavailable in this environment."
        )

    all_cmds: set[str] = set()
    for cmd in commands:
        tokens = shlex.split(cmd)
        if not tokens or tokens[0] != TOOL_NAME:
            continue
        if len(tokens) >= 2 and tokens[1] == "all":
            all_cmds.add(cmd)
            continue
        if len(tokens) >= 2 and tokens[1].startswith("-"):
            if any(tok in {"-i", "--input"} for tok in tokens[1:]):
                all_cmds.add(cmd)
    all_cmds = sorted(all_cmds)
    if not all_cmds:
        raise RuntimeError("No 'all' command examples found in docs.")

    with tempfile.TemporaryDirectory(prefix=f"{TOOL_NAME}_docs_smoke_") as tmpdir:
        fixture = _prepare_fixture_files(Path(tmpdir))
        # Many EN/JA pages intentionally repeat the same canonical invocation.
        # Sanitize first and execute each distinct CLI contract once.
        cases: dict[tuple[str, ...], str] = {}
        for raw in all_cmds:
            tokens = shlex.split(raw)
            if not tokens or tokens[0] != TOOL_NAME:
                continue
            args = tokens[1:]
            if not args or args[0].startswith("-"):
                args = ["all", *args]
            if args[0] != "all":
                continue
            dry_args = _sanitize_all_args(args, fixture)
            cases.setdefault(tuple(dry_args), raw)

        for dry_args_tuple, raw in cases.items():
            dry_args = list(dry_args_tuple)
            try:
                completed = subprocess.run(
                    [sys.executable, "-m", CLI_MODULE, *dry_args],
                    cwd=REPO_ROOT,
                    text=True,
                    capture_output=True,
                    timeout=DOCS_SMOKE_COMMAND_TIMEOUT_SEC,
                )
            except subprocess.TimeoutExpired as exc:
                raise RuntimeError(
                    f"[dry-run-smoke] timeout for docs command ({DOCS_SMOKE_COMMAND_TIMEOUT_SEC:g}s):\n"
                    f"  {raw}\n"
                    f"sanitized args:\n"
                    f"  {TOOL_NAME} {' '.join(dry_args)}"
                ) from exc
            if completed.returncode != 0:
                raise RuntimeError(
                    f"[dry-run-smoke] failed for docs command:\n  {raw}\n"
                    f"sanitized args:\n  {TOOL_NAME} {' '.join(dry_args)}\n\n"
                    f"stdout:\n{completed.stdout}\n\nstderr:\n{completed.stderr}"
                )
    print(
        f"[dry-run-smoke] validated {len(all_cmds)} docs examples "
        f"through {len(cases)} distinct plans "
        f"(timeout={DOCS_SMOKE_COMMAND_TIMEOUT_SEC:g}s)."
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()

    authored = extract_docs_commands()
    if not authored:
        raise RuntimeError("No commands were extracted from docs markdown code fences.")

    # Static validation retains EVERY authored command (including bracket-bearing
    # and data-literal examples); execution eligibility is classified separately.
    errors = validate_option_names(authored, root_cli)
    if errors:
        raise RuntimeError(
            "[option-smoke] docs option validation failed:\n" + "\n".join(errors)
        )
    print(f"[option-smoke] validated option names in {len(authored)} docs examples.")

    _run_help_smoke([cmd.text for cmd in authored])
    _run_all_dry_run_smoke([cmd.text for cmd in authored if cmd.executable])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
