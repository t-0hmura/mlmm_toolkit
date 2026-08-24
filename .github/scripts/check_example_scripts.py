#!/usr/bin/env python3
"""Make the advertised example shell scripts a checked path contract.

The README points users at working scripts under ``examples/``. This checker
proves each advertised script exists, is syntactically valid shell (``bash -n``),
and that every ``mlmm`` invocation in it uses only real live-CLI option names.
Scientific execution stays in the dedicated GPU smoke lanes; here we validate
authored syntax only.
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


def main() -> int:
    argparse.ArgumentParser(description=__doc__).parse_args()

    errors: list[str] = []
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
