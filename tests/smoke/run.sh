#!/usr/bin/env bash
# mlmm smoke tests — GPU required.
#
# This script assumes the calling environment already has:
#   - a Python with `mlmm` installed and importable
#   - AmberTools (antechamber / parmchk2 / tleap) on PATH
#   - CUDA available
#   - a writable scratch copy of this directory (artefacts land in `testNN*`)
# It does NOT activate conda, load modules, or contain HPC scheduler
# directives. Copy the fixtures to scratch, then run from that copy:
#
#   cp -a tests/smoke /path/to/scratch/mlmm-smoke
#   cd /path/to/scratch/mlmm-smoke
#   bash run.sh
#
# If you need an HPC scheduler wrapper, keep that in
# your own out-of-tree submission script and have it invoke this file as
# the body.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel 2>/dev/null || true)"
if [[ -n "$REPO_ROOT" ]]; then
  echo "ERROR: copy tests/smoke to writable scratch and run that copy; do not run the repository script in place." >&2
  exit 2
fi

# Deterministic run gate: pin PYTHONHASHSEED so Set / dict iteration order
# does not leak hash randomisation into the produced output. Combined with the
# UMA backend's deterministic-algorithms wiring this removes the remaining
# I/O-side non-determinism between two consecutive runs.
export PYTHONHASHSEED=0
# Reduce CUDA allocator fragmentation across the 40+ stage processes.
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

python - <<'PY'
from importlib.metadata import version
from pathlib import Path
import subprocess
import sys

import mlmm
import torch

installed = version("mlmm-toolkit")
if mlmm.__version__ != installed:
    raise SystemExit(
        "[smoke] BLOCKED: distribution/module version mismatch: "
        f"metadata={installed}, module={mlmm.__version__}"
    )
cli_version = subprocess.run(
    [sys.executable, "-m", "mlmm", "--version"],
    check=True,
    capture_output=True,
    text=True,
).stdout.strip().split()[-1]
if cli_version != installed:
    raise SystemExit(
        "[smoke] BLOCKED: distribution/CLI version mismatch: "
        f"metadata={installed}, cli={cli_version}"
    )
if not torch.cuda.is_available():
    raise SystemExit("[smoke] BLOCKED: CUDA is required by this lane")
print(
    "[smoke] package mlmm "
    f"version={installed} source={Path(mlmm.__file__).resolve()}"
)
PY
mlmm() { python -m mlmm "$@"; }
python assert_tr_cuda_parity.py

# Clean only artifacts authored by this harness. The digit-qualified glob must
# not be widened to `test*`, which would also match the repository's tests/.
rm -rf -- test[0-9]* pocket_r.pdb r_complex_layered.pdb r_complex_layered.cif r_complex_elem.pdb r_complex_fixalt.pdb

MLMM_COMPLEX_FREEZE_ATOMS="1,32"

# test1: extract
mlmm extract -i r_complex.pdb -c PRE -r 5.0 --no-exclude-backbone --ligand-charge 'PRE:0' -o pocket_r.pdb > test1.out 2>&1

# test2: define-layer
mlmm define-layer -i r_complex.pdb --model-pdb pocket_r.pdb --radius-freeze 8.0 -o r_complex_layered.pdb > test2.out 2>&1

# test3: mm-parm
mlmm mm-parm -i r_complex.pdb --ligand-charge 'PRE:0' > test3.out 2>&1

# test4: opt (grad)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode grad --max-cycles 5 --thresh gau_loose --dump --out-dir test4 > test4.out 2>&1

# test5: opt (hess)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --max-cycles 3 --thresh gau_loose --out-dir test5 > test5.out 2>&1

# test6: opt (hess, microiter)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --microiter --max-cycles 2 --thresh gau_loose --out-dir test6 > test6.out 2>&1

# test7: tsopt (grad / dimer)
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode grad --max-cycles 100 --thresh gau --out-dir test7 > test7.out 2>&1

# test8: tsopt (hess / rsirfo)
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --max-cycles 5 --thresh gau --out-dir test8 > test8.out 2>&1

# test9: freq
mlmm freq -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --out-dir test9 > test9.out 2>&1

# test10: irc
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-cycles 3 --out-dir test10 > test10.out 2>&1

# test11: dft (hf/sto-3g, cpu — gpu4pyscf may not be available in all envs)
mlmm dft -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --func-basis 'hf/sto-3g' --grid-level 0 --conv-tol 1e-5 --max-cycle 40 --engine cpu --out-dir test11 > test11.out 2>&1

# test12: scan (1D)
mlmm scan -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',3.5),('PRE 8 C1','PRE 8 C8',1.5)]" --max-step-size 2.0 --max-cycles 3 --no-preopt --no-endopt --out-dir test12 > test12.out 2>&1

# test13: scan2d
mlmm scan2d -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',1.4,1.8),('PRE 8 C1','PRE 8 C8',3.2,3.6)]" --max-step-size 0.4 --relax-max-cycles 100 --thresh gau_loose --out-dir test13 > test13.out 2>&1

# test14: scan3d
mlmm scan3d -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',1.4,1.8),('PRE 8 C1','PRE 8 C8',3.2,3.6),('PRE 8 C1','PRE 8 C7',1.4,1.8)]" --max-step-size 0.4 --relax-max-cycles 100 --thresh gau_loose --out-dir test14 > test14.out 2>&1

# test15: path-opt (gsm)
mlmm path-opt -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-nodes 5 --max-cycles 5 --no-preopt --no-climb --out-dir test15 > test15.out 2>&1

# test16: path-opt (dmf)
mlmm path-opt -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --mep-mode dmf --max-cycles 3 --no-preopt --out-dir test16 > test16.out 2>&1

# test17: path-search
mlmm path-search -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-cycles 5 --out-dir test17 > test17.out 2>&1

# test18: all (no tsopt/thermo/dft)
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge 'PRE:0' -q -1 -m 1 --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test18 > test18.out 2>&1

# test19: required positive MEP -> TSopt -> IRC -> thermo -> DFT handoff.
# The lane uses flattening to obtain a first-order saddle and deterministic
# execution to keep the release comparison reproducible. Saddle certification
# requires exactly one imaginary mode; its magnitude is not a pass/fail gate.
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 4.0 --ligand-charge 'PRE:0' -q -1 -m 1 --deterministic --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --tsopt --thermo --dft --flatten --irc-never-stop --tsopt-max-cycles 2000 --dft-func-basis 'hf/sto-3g' --dft-grid-level 0 --dft-conv-tol 1e-5 --dft-max-cycle 40 --dft-engine cpu --out-dir test19 > test19.out 2>&1
python assert_release_result.py all test19 --require-thermo --require-dft >> test19.out 2>&1

# test20: all (--parm + --model-pdb override, reuse test19 outputs)
mapfile -t test19_parms < <(find test19/mm_parm -maxdepth 1 -type f -name '*.parm7' -print)
if [[ "${#test19_parms[@]}" -ne 1 ]]; then
  echo "[smoke] FAIL test20: expected exactly one reusable test19 parm7, found ${#test19_parms[@]}" >&2
  exit 1
fi
mlmm all -i r_complex.pdb p_complex.pdb --parm "${test19_parms[0]}" --model-pdb test19/ml_region.pdb -q -1 -m 1 --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test20 > test20.out 2>&1

# test21: tsopt (radius-hessian 0.0)
mlmm tsopt -i p_complex.pdb --parm p_complex.parm7 --model-pdb pocket_r.pdb --no-detect-layer -q -1 -m 1 --opt-mode grad --max-cycles 5 --radius-hessian 0.0 --active-dof-mode ml-only --thresh gau_loose --out-dir test21 > test21.out 2>&1

# test22: tsopt (radius-hessian 3.6)
mlmm tsopt -i p_complex.pdb --parm p_complex.parm7 --model-pdb pocket_r.pdb --no-detect-layer -q -1 -m 1 --opt-mode grad --max-cycles 5 --radius-hessian 3.6 --active-dof-mode ml-only --thresh gau_loose --out-dir test22 > test22.out 2>&1
python - <<'PY'
import re
from pathlib import Path

def initial_active_count(name: str) -> int:
    text = Path(name).read_text(encoding="utf-8")
    match = re.search(r"\[tsopt\] H_act=\d+ active_atoms=(\d+)", text)
    if match is None:
        raise SystemExit(f"[smoke] FAIL: {name} lacks the initial Hessian coverage record")
    return int(match.group(1))

ml_only = initial_active_count("test21.out")
expanded = initial_active_count("test22.out")
if expanded <= ml_only:
    raise SystemExit(
        "[smoke] FAIL: --radius-hessian 3.6 did not expand the Dimer Hessian "
        f"coverage ({expanded} <= {ml_only})"
    )
print(f"[smoke] PASS: radius-hessian coverage expanded {ml_only} -> {expanded} atoms")
PY

# test23: opt --dry-run
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode grad --dry-run --out-dir test23 > test23.out 2>&1

# test24: tsopt --dry-run
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --dry-run --out-dir test24 > test24.out 2>&1

# test25: freq --dry-run
mlmm freq -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --dry-run --out-dir test25 > test25.out 2>&1

# test26: scan --dry-run
mlmm scan -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',3.5),('PRE 8 C1','PRE 8 C8',1.5)]" --dry-run --out-dir test26 > test26.out 2>&1

# test27: dft --dry-run
mlmm dft -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --dry-run --out-dir test27 > test27.out 2>&1

# test28: path-search --dry-run
mlmm path-search -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --dry-run --out-dir test28 > test28.out 2>&1

# test29: irc --dry-run
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --dry-run --out-dir test29 > test29.out 2>&1

# test30: add-elem-info
mlmm add-elem-info -i r_complex.pdb -o r_complex_elem.pdb > test30.out 2>&1

# test31: trj2fig
mlmm trj2fig -i test4/optimization_trj.xyz -o test31.png > test31.out 2>&1

# test32: energy-diagram
mlmm energy-diagram -i "[0, 12.5, 4.3, 18.7, -1.2]" -o test32.png > test32.out 2>&1

# test33: oniom-export
mlmm oniom-export --parm p_complex.parm7 -i r_complex_layered.pdb --model-pdb pocket_r.pdb -q -1 -m 1 -o test33.gjf > test33.out 2>&1

# --- Bond-summary, fix-altloc, oniom-import ---

# test34: bond-summary (two layered PDBs)
mlmm bond-summary -i r_complex_layered.pdb p_complex_layered.pdb > test34.out 2>&1

# test35: fix-altloc
mlmm fix-altloc -i r_complex.pdb -o r_complex_fixalt.pdb > test35.out 2>&1

# test36: oniom-import (Gaussian input → layered PDB)
mlmm oniom-import -i test33.gjf -o test36 > test36.out 2>&1

# --- refine-path ---

# test37: all (--refine-path)
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge 'PRE:0' -q -1 -m 1 --refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test37 > test37.out 2>&1

# test38: retired electronic embedding fails before calculation
if mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode grad --max-cycles 3 --thresh gau_loose --embedcharge --embedcharge-cutoff 6.0 --out-dir test38 > test38.out 2>&1; then
  echo "[smoke] FAIL test38: --embedcharge was accepted" >&2
  exit 1
fi
grep -Fq "Electronic embedding is unavailable in v0.3.3" test38.out

# --- Opt-in TS and IRC methods ---

# test39: tsopt --opt-mode trim (Helgaker trust-region image-min; non-microiter)
mlmm tsopt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode trim --no-microiter --max-cycles 5 --thresh gau_loose --out-dir test39 > test39.out 2>&1

# test40: tsopt --opt-mode rsprfo (Banerjee P-RFO; non-microiter)
mlmm tsopt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode rsprfo --no-microiter --max-cycles 5 --thresh gau_loose --out-dir test40 > test40.out 2>&1

# test42: irc --irc-pos-def (PSD-Hessian convergence guard)
mlmm irc -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-cycles 3 --irc-pos-def --out-dir test42 > test42.out 2>&1

# test43: opt --print-every 3 (diagnostic output throttle)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --max-cycles 5 --thresh gau_loose --print-every 3 --out-dir test43 > test43.out 2>&1

# --- Determinism gate ---

# test44: `all` pipeline determinism GATE (`--deterministic`, ONIOM end-to-end).
# Runs the full pipeline twice with identical inputs / args + `--deterministic`
# and REQUIRES the two runs to be bit-identical. Default (non-deterministic) GPU
# runs carry ~ULP scatter/atomic non-determinism and are not asserted here;
# `--deterministic` enables torch deterministic algorithms and MUST be
# bit-reproducible, so any drift is a real regression and fails the smoke.
#
# The MM parm (antechamber AM1-BCC ligand charges) is a NON-deterministic INPUT,
# not part of the compute this gate exercises: sqm's AM1 SCF for this ligand is
# poorly convergent, so its early-stop point — and thus the charges — vary
# run-to-run. Regenerating it in both runs would make the gate test antechamber,
# not `--deterministic` compute reproducibility. So build the parm once (run a)
# and REUSE it via `--parm` in run b, giving both runs identical MM charges;
# `--deterministic` then yields bit-identical geometry/MEP output. (For full
# end-to-end reproducibility across separate invocations, pass a fixed `--parm`.)
det_args="-i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --deterministic"
mlmm all $det_args --out-dir test44_a > test44_a.out 2>&1
mapfile -t test44_parms < <(find test44_a/mm_parm -maxdepth 1 -type f -name '*.parm7' -print)
if [[ "${#test44_parms[@]}" -ne 1 ]]; then
  echo "[smoke] FAIL test44: expected exactly one reusable parm7, found ${#test44_parms[@]}" >&2
  exit 1
fi
mlmm all $det_args --parm "${test44_parms[0]}" --out-dir test44_b > test44_b.out 2>&1
# Run b is given the parm via `--parm`, so by design it never runs the mm-parm
# stage and never writes `mm_parm/`. That directory is the gate's fixed INPUT,
# not its output, so comparing it would fail on the very workaround above.
find test44_a -type f \( -name "*.pdb" -o -name "*.xyz" \) -not -path '*/mm_parm/*' -printf '%P\n' | LC_ALL=C sort > test44_a.manifest
find test44_b -type f \( -name "*.pdb" -o -name "*.xyz" \) -not -path '*/mm_parm/*' -printf '%P\n' | LC_ALL=C sort > test44_b.manifest
if ! cmp -s test44_a.manifest test44_b.manifest; then
  echo "[smoke] FAIL test44: deterministic runs produced different file manifests" > test44.out
  comm -3 test44_a.manifest test44_b.manifest >> test44.out
  cat test44.out
  exit 1
fi
total=$(wc -l < test44_a.manifest)
if [ "$total" -eq 0 ]; then
  echo "[smoke] FAIL test44: deterministic gate found no PDB/XYZ artifacts" > test44.out
  cat test44.out
  exit 1
fi
drifted=0
while IFS= read -r rel; do
  if ! cmp -s "test44_a/$rel" "test44_b/$rel"; then
    drifted=$((drifted + 1))
    echo "DRIFT: $rel" >> test44.out
  fi
done < test44_a.manifest
echo "[det_check] compared $total PDB/XYZ files; $drifted differ" >> test44.out
if [ "$drifted" -ne 0 ]; then
  echo "[smoke] FAIL test44: --deterministic runs differ" >> test44.out
  cat test44.out
  exit 1
fi

# --- --coord-type CLI plumbing (throttled, fast) ---

# test45: `all --coord-type cart` — explicit cart (== default), verifies CLI plumbing.
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --coord-type cart --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test45 > test45.out 2>&1

# test46: `all --coord-type dlc` — DLC propagated to the child opt / path-opt
# stages this run enables (it passes --no-tsopt).
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --coord-type dlc --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test46 > test46.out 2>&1

# test47: `sp` (single-point ONIOM) — energy + forces.
mlmm sp -i r_complex_layered.pdb --real-parm7 p_complex.parm7 -q -1 -m 1 --out-dir test47 > test47.out 2>&1

# test48: `sp --hess` — energy + forces + ONIOM Hessian (default FiniteDifference).
mlmm sp -i r_complex_layered.pdb --real-parm7 p_complex.parm7 -q -1 -m 1 --hess --out-dir test48 > test48.out 2>&1

# --- Full-pipeline release-gate runs ---
# test49 is the untrottled one: the canonical `all` flow with default
# convergence thresholds and production-realistic optimizer cycle budgets, so it
# takes substantially longer (~30-90 min for ONIOM) than the throttled tests
# above and is the "does the pipeline actually finish on a real input" gate.
# test50 stays capped, for the reason spelled out at its own comment below.

# test49: full `all` cart — default thresh, no max-cycles cap.
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --out-dir test49 > test49.out 2>&1

# test50: `all` dlc — verifies the DLC code-path lights up end-to-end.
# Capped at max-cycles 5 + thresh gau_loose + --no-tsopt/thermo/dft because
# DLC GSM on this 122-atom complex needs hundreds of cycles to converge with
# default `gau` thresh (3h+ on consumer GPU) and the post-stages (TS / IRC /
# DFT) depend on a converged HEI from the MEP — silently broken structure
# handoff otherwise. test49 (cart) keeps the no-cap default-behaviour check.
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --coord-type dlc --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test50 > test50.out 2>&1

# --- Per-stage internal-coordinate code-path verification ---
# Each test is scoped at a 2-3 cycle cap (plus gau_loose where the stage
# needs it) so it exercises the
# coordinate paths without requiring convergence. Frequency analysis remains
# Cartesian because its PHVA contract consumes a Cartesian Hessian directly.

# test50a: opt --coord-type dlc
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --coord-type dlc --max-cycles 3 --thresh gau_loose --out-dir test50a_opt_dlc > test50a_opt_dlc.out 2>&1

# test50b: opt --opt-mode hess --coord-type dlc (microiter+DLC regression: ML internals, MM cart twin)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --coord-type dlc --max-cycles 3 --thresh gau_loose --out-dir test50b_opt_hess_dlc > test50b_opt_hess_dlc.out 2>&1

# test50c: opt --coord-type dlc with explicit frozen atoms
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --coord-type dlc --freeze-atoms "$MLMM_COMPLEX_FREEZE_ATOMS" --max-cycles 3 --thresh gau_loose --out-dir test50c_opt_freeze_dlc > test50c_opt_freeze_dlc.out 2>&1
python check_frozen_atoms.py r_complex_layered.pdb test50c_opt_freeze_dlc/final_geometry.pdb "$MLMM_COMPLEX_FREEZE_ATOMS" test50c >> test50c_opt_freeze_dlc.out 2>&1

# test50d: Cartesian scan with explicit frozen atoms.
# Small non-reactive target: this checks coordinate integrity, not chemistry.
mlmm scan -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --freeze-atoms "$MLMM_COMPLEX_FREEZE_ATOMS" --scan-lists "[('PRE 8 O1\'','PRE 8 C3',1.5)]" --max-step-size 0.1 --max-cycles 3 --no-endopt --out-dir test50d_scan_freeze_cart > test50d_scan_freeze_cart.out 2>&1
python check_frozen_atoms.py r_complex_layered.pdb test50d_scan_freeze_cart/stage_01/result.pdb "$MLMM_COMPLEX_FREEZE_ATOMS" test50d >> test50d_scan_freeze_cart.out 2>&1
if grep -q "Covalent-bond changes (start vs final): Yes" test50d_scan_freeze_cart.out; then
  echo "[bond-check] test50d: unexpected covalent-bond changes in non-reactive Cartesian freeze scan" >> test50d_scan_freeze_cart.out
  exit 1
fi

# test50e: Hessian TS microiteration with DLC and frozen atoms. This is the
# partial-Cartesian-Hessian -> internal-coordinate handoff regression.
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --coord-type dlc --freeze-atoms "$MLMM_COMPLEX_FREEZE_ATOMS" --microiter --max-cycles 2 --thresh gau_loose --out-dir test50e_ts_hess_dlc > test50e_ts_hess_dlc.out 2>&1

# test50g: opt --coord-type redund
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --coord-type redund --max-cycles 3 --thresh gau_loose --out-dir test50g_opt_redund > test50g_opt_redund.out 2>&1

# test50k: opt --coord-type tric
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --coord-type tric --max-cycles 3 --thresh gau_loose --out-dir test50k_opt_tric > test50k_opt_tric.out 2>&1

# --- Multi-mode flag code-path verify (single-stage) ---

# test50m: opt --precision fp64 (UMA backend, alternate precision)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --precision fp64 --max-cycles 3 --thresh gau_loose --out-dir test50m_opt_fp64 > test50m_opt_fp64.out 2>&1

# test50n: opt --precision fp32 (explicit fp32 dispatch; default is fp32, this pins the explicit path alongside test50m fp64)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --precision fp32 --max-cycles 3 --thresh gau_loose --out-dir test50n_opt_fp32 > test50n_opt_fp32.out 2>&1

# test50p: opt --mm-backend openmm (alternate MM backend; analytical Hessian path → FD)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --mm-backend openmm --max-cycles 3 --thresh gau_loose --out-dir test50p_opt_openmm > test50p_opt_openmm.out 2>&1

# test50q: opt --link-atom-method fixed (legacy 1.09/1.01 Å placement)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --link-atom-method fixed --max-cycles 3 --thresh gau_loose --out-dir test50q_opt_linkfixed > test50q_opt_linkfixed.out 2>&1

# --- Non-default MLIP backend, full pipeline ---

# test51: full `all` with `--backend orb` — exercises the non-default MLIP backend.
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --backend orb --out-dir test51 > test51.out 2>&1
python assert_release_result.py provenance test51 --expected-backend orb --expected-model orb_v3_conservative_omol --expected-precision fp64 >> test51.out 2>&1

# ---- Subcommand-specific regression coverage ----
# test52: opt --mm-only (MM-only minimization; skips the MLIP component entirely)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --mm-only --opt-mode grad --max-cycles 3 --thresh gau_loose --out-dir test52_opt_mmonly > test52_opt_mmonly.out 2>&1

# test53: freq --active-dof-mode ml-only (alternate PHVA active-DOF subspace)
mlmm freq -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --active-dof-mode ml-only --max-write 5 --out-dir test53_freq_mlonly > test53_freq_mlonly.out 2>&1

# test54: freq --hessian-calc-mode Analytical (workflow analytical path)
mlmm freq -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --hessian-calc-mode Analytical --max-write 5 --out-dir test54_freq_anahess > test54_freq_anahess.out 2>&1

# test55: irc --hessian-calc-mode analytical (IRC initial Hessian path)
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --hessian-calc-mode analytical --max-cycles 2 --out-dir test55_irc_anahess > test55_irc_anahess.out 2>&1

# test56: irc --mm-backend openmm (MM Hessian via OpenMM finite-difference)
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --mm-backend openmm --max-cycles 2 --out-dir test56_irc_openmm > test56_irc_openmm.out 2>&1

# test57: irc --freeze-atoms (DOF-reduction / reduced-Hessian projection path)
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --freeze-atoms 1,2,3 --max-cycles 2 --out-dir test57_irc_freeze > test57_irc_freeze.out 2>&1

# test58: DFT also rejects retired electronic embedding before SCF setup
if mlmm dft -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --func-basis 'hf/sto-3g' --grid-level 0 --conv-tol 1e-5 --max-cycle 40 --engine cpu --embedcharge --embedcharge-cutoff 8.0 --out-dir test58_dft_embed > test58_dft_embed.out 2>&1; then
  echo "[smoke] FAIL test58: DFT accepted --embedcharge" >&2
  exit 1
fi
grep -Fq "Electronic embedding is unavailable in v0.3.3" test58_dft_embed.out

# test59: dft --link-atom-method fixed (legacy 1.09/1.01 Å link-atom placement)
mlmm dft -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --func-basis 'hf/sto-3g' --grid-level 0 --conv-tol 1e-5 --max-cycle 40 --engine cpu --link-atom-method fixed --out-dir test59_dft_linkfixed > test59_dft_linkfixed.out 2>&1

# test60: path-search --mep-mode dmf (Direct Max Flux vs GrowingString)
mlmm path-search -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --mep-mode dmf --max-cycles 3 --no-preopt --out-dir test60_psdmf > test60_psdmf.out 2>&1

# test61: all --scan-lists (single-PDB scan->path mode of `all`, distinct from the multi-PDB MEP branch)
mlmm all -i r_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',3.5),('PRE 8 C1','PRE 8 C8',1.5)]" --no-refine-path --max-cycles 3 --thresh gau_loose --no-tsopt --no-thermo --no-dft --out-dir test61_all_scan > test61_all_scan.out 2>&1

# --- refine-path opt-in (recursive path_search) extra coverage ---
# The `all` default is now single-pass path-opt; exercise the recursive
# path_search opt-in (`--refine-path`) in scan->path mode too (test37 already
# covers the multi-input endpoint MEP with --refine-path).

# test63: extract MULTI-INPUT via space-separated '-i a.pdb b.pdb' (one flag, two paths).
# Regression guard: a single -i with several space-separated paths must NOT drop the 2nd input.
# A single -o yields one multi-MODEL PDB, so both endpoints must appear (-> exactly 2 MODEL records).
mlmm extract -i r_complex.pdb p_complex.pdb -c PRE -r 5.0 --no-exclude-backbone --ligand-charge 'PRE:0' -o pocket_multi.pdb > test63_multi_extract.out 2>&1
n_models=$(grep -c '^MODEL' pocket_multi.pdb 2>/dev/null)
if [ "${n_models:-0}" -ne 2 ]; then
  echo "[extract-multi] test63: space-separated '-i a b' yielded ${n_models:-0} MODEL records (expected 2); the 2nd input was dropped" >> test63_multi_extract.out
  exit 1
fi

# test62: all --scan-lists --refine-path (single-PDB scan -> recursive path_search)
mlmm all -i r_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',3.5),('PRE 8 C1','PRE 8 C8',1.5)]" --refine-path --max-cycles 3 --thresh gau_loose --no-tsopt --no-thermo --no-dft --out-dir test62_rp_scan > test62_rp_scan.out 2>&1

# test64: --backend-model routing — a non-default model must reach the resolved
# runtime header. Dry-run avoids downloading the alternate model.
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --backend-model uma-m-1p1 --dry-run --out-dir test64_backend_model > test64_backend_model.out 2>&1
grep -Eq '^\[backend\] uma \(uma-m-1p1, fp32\)$' test64_backend_model.out || { echo "[smoke] FAIL test64: non-default backend model missing from resolved runtime summary" >> test64_backend_model.out; exit 1; }

# Build an mmCIF equivalent of the layered fixture while exercising identifiers
# that cannot be represented in fixed-column PDB.
python - <<'PY'
from Bio.PDB import MMCIFIO, PDBParser

structure = PDBParser(QUIET=True).get_structure("smoke", "r_complex_layered.pdb")
chains = list(structure.get_chains())
assert len(chains) == 1
chain = chains[0]
chain.id = "LONG_CHAIN"
ligands = [res for res in chain if res.get_resname().strip() == "PRE"]
assert len(ligands) == 1
het, _resseq, icode = ligands[0].id
ligands[0].id = (het, 10001, icode)
writer = MMCIFIO()
writer.set_structure(structure)
writer.save("r_complex_layered.cif")
PY

# test65: a real ML/MM optimization crosses the mmCIF bridge and restores the
# original long chain and five-digit residue identifier in its public output.
mlmm opt -i r_complex_layered.cif --parm p_complex.parm7 -q -1 -m 1 --max-cycles 1 --thresh gau_loose --out-dir test65_opt_cif > test65_opt_cif.out 2>&1
test -s test65_opt_cif/final_geometry.pdb || { echo "[smoke] FAIL test65: final PDB missing" >> test65_opt_cif.out; exit 1; }
test -s test65_opt_cif/final_geometry.cif || { echo "[smoke] FAIL test65: final CIF missing" >> test65_opt_cif.out; exit 1; }
grep -q 'LONG_CHAIN' test65_opt_cif/final_geometry.cif || { echo "[smoke] FAIL test65: auth chain was not restored" >> test65_opt_cif.out; exit 1; }
grep -q '10001' test65_opt_cif/final_geometry.cif || { echo "[smoke] FAIL test65: auth residue number was not restored" >> test65_opt_cif.out; exit 1; }

# test66: exact chain/residue-name/residue-number selection remains stable
# after normalization, with both internal PDB and identifier-preserving CIF.
mlmm extract -i r_complex_layered.cif -c 'LONG_CHAIN:PRE:10001' -r 0.1 --no-add-linkh -o test66_model_from_cif.pdb -v 0 > test66_extract_cif.out 2>&1
test -s test66_model_from_cif.pdb || { echo "[smoke] FAIL test66: extracted PDB missing" >> test66_extract_cif.out; exit 1; }
test -s test66_model_from_cif.cif || { echo "[smoke] FAIL test66: extracted CIF missing" >> test66_extract_cif.out; exit 1; }

# test67/68: a partial Hessian dumped by freq must be consumed unchanged by
# IRC, including its active-DOF metadata. Never-stop is verified at runtime.
# The dump must use IRC's own active-DOF basis (ML + MovableMM, i.e. freq's
# default). IRC has no --active-dof-mode/--hess-cutoff flag and always analyzes
# every movable atom, so an ml-only dump would be rejected as a basis mismatch;
# the partial nature is still exercised because --freeze-atoms keeps it < full.
mlmm freq -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --freeze-atoms 1,2,3 --max-write 1 --dump-hess test67_freq/hessian.npz --out-json --out-dir test67_freq > test67_freq.out 2>&1
test -s test67_freq/hessian.npz || { echo "[smoke] FAIL test67: dumped Hessian missing" >> test67_freq.out; exit 1; }
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --freeze-atoms 1,2,3 --read-hess test67_freq/hessian.npz --never-stop --config never_stop_config.yaml --max-cycles 2 --out-json --out-dir test68_irc_handoff > test68_irc_handoff.out 2>&1
python assert_release_result.py irc-handoff test68_irc_handoff --hessian-file test67_freq/hessian.npz >> test68_irc_handoff.out 2>&1

# A same-size Hessian from a different geometry must be rejected before IRC.
python - <<'PY'
from pathlib import Path

lines = Path("p_complex_layered.pdb").read_text(encoding="utf-8").splitlines(True)
for index, line in enumerate(lines):
    if line.startswith(("ATOM  ", "HETATM")):
        x = float(line[30:38]) + 0.100
        lines[index] = line[:30] + f"{x:8.3f}" + line[38:]
        break
Path("test68_wrong_geometry.pdb").write_text("".join(lines), encoding="utf-8")
PY
rc=0
mlmm irc -i test68_wrong_geometry.pdb --parm p_complex.parm7 -q -1 -m 1 --freeze-atoms 1,2,3 --read-hess test67_freq/hessian.npz --max-cycles 1 --out-dir test68_wrong > test68_wrong.out 2>&1 || rc=$?
if [ "$rc" -eq 0 ] || ! grep -Eq 'coordinates do not match|PES identity does not match' test68_wrong.out; then
  echo "[smoke] FAIL test68: stale same-size Hessian was not rejected" >> test68_wrong.out
  exit 1
fi

# test69: standalone --ref-mode is an actual Cartesian mode vector. The
# all-workflow path-tangent handoff is exercised by required-positive test19.
python - <<'PY'
from pathlib import Path
import numpy as np

def read(path):
    coords = []
    layers = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            layers.append(float(line[60:66]))
    return np.asarray(coords), np.asarray(layers)

reactant, _ = read("r_complex_layered.pdb")
product, layers = read("p_complex_layered.pdb")
mode = (reactant - product).reshape(-1)
active = np.repeat(np.isclose(layers, 0.0, atol=1.0), 3)
if np.linalg.norm(mode[active]) <= 1.0e-8:
    raise SystemExit("reference path tangent is zero in the active ML Hessian space")
np.savetxt("test69_reference_mode.txt", mode)
PY
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --ref-mode test69_reference_mode.txt --max-cycles 2 --thresh gau_loose --out-json --out-dir test69_ref_mode > test69_ref_mode.out 2>&1
python assert_release_result.py tsopt-reference test69_ref_mode >> test69_ref_mode.out 2>&1

# test70: YAML backend-model/precision settings reach a real calculation while
# an explicit CLI max-cycles value retains precedence.
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --config runtime_config.yaml --max-cycles 1 --out-json --out-dir test70_config > test70_config.out 2>&1
python assert_release_result.py opt-config test70_config --expected-model uma-s-1p2 --expected-max-cycles 1 --expected-precision fp64 --expected-link-atom-method fixed --expected-thresh gau_loose >> test70_config.out 2>&1

# test71: a user ASE calculator with only energy/forces supports the ML-region
# finite-difference Hessian path.
mlmm sp -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --calc-file harmonic_calc.py --hess --hessian-calc-mode FiniteDifference --out-json --out-dir test71_custom > test71_custom.out 2>&1
python assert_release_result.py sp-hessian test71_custom >> test71_custom.out 2>&1

# test72: an explicit analytical Hessian request cannot silently fall back when
# predictor workers remove the autograd model.
rc=0
mlmm sp -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --hess --hessian-calc-mode Analytical --workers 2 --out-dir test72_workers > test72_workers.out 2>&1 || rc=$?
if [ "$rc" -ne 1 ] || ! grep -Fq "Analytical Hessian cannot be combined with workers>1: the parallel predictor exposes no autograd model. Use workers=1 or select hessian_calc_mode='FiniteDifference'." test72_workers.out || [ -e test72_workers/hessian.npy ]; then
  echo "[smoke] FAIL test72: Analytical + workers>1 was not rejected exactly" >> test72_workers.out
  exit 1
fi

# test73: ORCA QM/MM export/import round-trip preserves coordinates and all
# three layer classes without requiring an ORCA installation.
python - <<'PY'
from pathlib import Path

lines = Path("r_complex_layered.pdb").read_text(encoding="utf-8").splitlines(True)
changed = False
for index in range(len(lines) - 1, -1, -1):
    line = lines[index]
    if line.startswith(("ATOM  ", "HETATM")) and abs(float(line[60:66]) - 10.0) <= 1.0:
        lines[index] = line[:60] + f"{20.0:6.2f}" + line[66:]
        changed = True
        break
if not changed:
    raise SystemExit("could not create a frozen-MM layer for ORCA round-trip")
Path("test73_three_layer.pdb").write_text("".join(lines), encoding="utf-8")
PY
mlmm oniom-export --parm p_complex.parm7 -i test73_three_layer.pdb --model-pdb pocket_r.pdb -q -1 -m 1 --mode orca --no-convert-orcaff -o test73_orca.inp > test73_orca_export.out 2>&1
for token in '! QMMM' 'QMAtoms {' 'ActiveAtoms {' 'Charge_Total -1' '* xyz -1 1'; do
  grep -Fq "$token" test73_orca.inp || { echo "[smoke] FAIL test73: ORCA input missing $token" >> test73_orca_export.out; exit 1; }
done
mlmm oniom-import -i test73_orca.inp --ref-pdb test73_three_layer.pdb -o test73_orca_import > test73_orca_import.out 2>&1
test -s test73_orca_import.xyz || { echo "[smoke] FAIL test73: restored XYZ missing" >> test73_orca_import.out; exit 1; }
test -s test73_orca_import_layered.pdb || { echo "[smoke] FAIL test73: restored layered PDB missing" >> test73_orca_import.out; exit 1; }
sed -n '2p' test73_orca_import.xyz | grep -Fq 'mode=orca' || { echo "[smoke] FAIL test73: restored XYZ omits ORCA provenance" >> test73_orca_import.out; exit 1; }
sed -n '2p' test73_orca_import.xyz | grep -Fq 'q=-1 m=1' || { echo "[smoke] FAIL test73: restored XYZ charge/multiplicity mismatch" >> test73_orca_import.out; exit 1; }
python assert_orca_roundtrip.py test73_three_layer.pdb test73_orca_import_layered.pdb >> test73_orca_import.out 2>&1

# test74: force a known higher-order candidate through the actual flatten
# branch. The checker requires n_imag>1 before flattening, an executed RS-I-RFO
# flatten iteration, and no increase in saddle order.
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --no-microiter --flatten --config flatten_branch_config.yaml --max-cycles 50 --thresh gau_loose --out-json --out-dir test74_flatten > test74_flatten.out 2>&1
python assert_flatten_branch.py test74_flatten.out test74_flatten/result.json >> test74_flatten.out 2>&1

# Numerical analytical-vs-FD agreement for every backend installed in the
# default strict environment. MACE/AIMNet2 use this same required wrapper in
# their dependency-isolated cluster environments.
bash run_backend_hessian.sh uma orb > backend_hessian.out 2>&1

echo "[smoke] PASS: required GPU, ML/MM, Hessian-handoff, and structure-I/O lane completed with zero skips."
