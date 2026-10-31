#!/usr/bin/env python3
"""Make the advertised example shell scripts a checked path contract.

The README points users at working scripts under ``examples/``. This checker
proves each advertised script exists, is syntactically valid shell (``bash -n``),
and that every ``mlmm`` invocation in it uses only real live-CLI option names.
Scientific execution stays in the dedicated GPU smoke lanes; here we validate
authored syntax only. It also checks the full-system example structures that the
Colab notebook loads (atom counts) and the options of their run scripts.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from docs_command_contract import (  # noqa: E402
    REPO_ROOT,
    extract_shell_commands,
    load_root_cli,
    public_shell_examples,
    validate_option_names,
)

# Floor guarding against a silently-empty extraction. The current advertised
# scripts contribute 49 invocations (toy_system 32, methyltransferase 15,
# BezA 2).
MIN_INVOCATIONS = 40
# Atoms per endpoint structure of the full-system examples.
FULL_SYSTEM_ATOMS = {"comt": 3420, "beza": 9215}


def _check_full_system_examples() -> list[str]:
    """The COMT and BezA endpoints keep their atom counts, BezA's endpoint run
    refines the path, and COMT's run script selects its centre and charges."""
    examples = REPO_ROOT / "examples"
    errors: list[str] = []
    for name, expected_atoms in FULL_SYSTEM_ATOMS.items():
        for endpoint in ("1.R.pdb", "3.P.pdb"):
            path = examples / name / endpoint
            lines = path.read_text(encoding="utf-8").splitlines()
            atoms = sum(line.startswith(("ATOM  ", "HETATM")) for line in lines)
            if atoms != expected_atoms:
                errors.append(
                    f"examples/{name}/{endpoint}: {atoms} atoms; expected {expected_atoms}"
                )
    beza_parts = (examples / "beza" / "run.sh").read_text(encoding="utf-8").split("mlmm all", 2)
    if len(beza_parts) < 2 or "--refine-path" not in beza_parts[1]:
        errors.append("examples/beza/run.sh: the first `mlmm all` command must pass --refine-path")
    comt_script = (examples / "comt" / "run.sh").read_text(encoding="utf-8")
    for fragment in ("-c 'CAT,SAM,MG'", "-l 'CAT:-1,SAM:0,MG:2' -r 4.0"):
        if fragment not in comt_script:
            errors.append(f"examples/comt/run.sh: missing {fragment!r}")
    return errors


def main() -> int:
    argparse.ArgumentParser(description=__doc__).parse_args()

    errors: list[str] = _check_full_system_examples()
    scripts = public_shell_examples()

    missing = [s for s in scripts if not s.exists()]
    for s in missing:
        errors.append(f"missing advertised example script: {s.relative_to(REPO_ROOT)}")
    if errors:
        print("[example-scripts] FAIL:")
        for e in errors:
            print(f"- {e}")
        return 1

    for s in scripts:
        text = s.read_text(encoding="utf-8")
        lines = text.splitlines()
        rel = s.relative_to(REPO_ROOT)
        if not lines or lines[0] != "#!/usr/bin/env bash":
            errors.append(f"{rel}: first line must be '#!/usr/bin/env bash'")
        if "set -euo pipefail" not in lines[:5]:
            errors.append(f"{rel}: missing 'set -euo pipefail' near the top")
        if not s.stat().st_mode & 0o111:
            errors.append(f"{rel}: script is not executable")
        if "BASH_SOURCE[0]" not in text:
            errors.append(f"{rel}: fixture paths are not anchored to the script directory")
        completed = subprocess.run(
            ["bash", "-n", str(s)], text=True, capture_output=True
        )
        if completed.returncode != 0:
            errors.append(
                f"bash -n failed for {s.relative_to(REPO_ROOT)}:\n{completed.stderr.strip()}"
            )

    stepwise = REPO_ROOT / "examples" / "methyltransferase" / "run_stepwise.sh"
    stepwise_text = stepwise.read_text(encoding="utf-8")
    if (
        "${ENERGIES:-}" not in stepwise_text
        or "${ENDPOINT_LABELS:-}" not in stepwise_text
        or "freq_reac" in stepwise_text
        or "freq_prod" in stepwise_text
    ):
        errors.append(
            "examples/methyltransferase/run_stepwise.sh: endpoint identity and energies must be explicit"
        )

    commands = extract_shell_commands(scripts)
    per_script: dict[str, int] = {}
    for cmd in commands:
        per_script[cmd.rel] = per_script.get(cmd.rel, 0) + 1
    for s in scripts:
        rel = str(s.relative_to(REPO_ROOT))
        if per_script.get(rel, 0) < 1:
            errors.append(f"no mlmm invocations discovered in {rel}")
    if len(commands) < MIN_INVOCATIONS:
        errors.append(
            f"discovered {len(commands)} invocations; expected >= {MIN_INVOCATIONS}"
        )

    root_cli = load_root_cli()
    errors.extend(validate_option_names(commands, root_cli))

    if errors:
        print("[example-scripts] FAIL:")
        for e in errors:
            print(f"- {e}")
        return 1

    print(
        f"[example-scripts] OK: {len(scripts)} scripts, "
        f"{len(commands)} invocations validated."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
