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

# Clean only artifacts authored by this harness. The positive-lane input PDBs
# share the test73 prefix and are static fixtures, not run output.
for artifact in test[0-9]*; do
  case "$artifact" in
    test73_r_complex.pdb|test73_p_complex.pdb) ;;
    *) rm -rf -- "$artifact" ;;
  esac
done
rm -rf -- pocket_r.pdb r_complex_layered.pdb r_complex_layered.cif r_complex_elem.pdb r_complex_fixalt.pdb
for fixture in test73_r_complex.pdb test73_p_complex.pdb; do
  test -s "$fixture" || { echo "[smoke] BLOCKED: required fixture missing after cleanup: $fixture" >&2; exit 1; }
done

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
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode grad --max-cycles 100 --thresh gau --out-json --out-dir test7 > test7.out 2>&1
python assert_release_result.py tsopt-optimizer test7 --expected-mode grad --expected-optimizer dimer >> test7.out 2>&1

# test8: tsopt (hess / rsprfo, microiteration default)
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --max-cycles 5 --thresh gau --out-dir test8 > test8.out 2>&1

# test9: freq
mlmm freq -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --out-dir test9 > test9.out 2>&1

# test10: irc
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-cycles 3 --out-dir test10 > test10.out 2>&1

# test11: dft (hf/sto-3g, cpu — gpu4pyscf may not be available in all envs)
mlmm dft -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --func-basis 'hf/sto-3g' --grid-level 0 --conv-tol 1e-5 --max-cycle 40 --engine cpu --out-dir test11 > test11.out 2>&1

# test12: scan (1D)
mlmm scan -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',3.5),('PRE 8 C1','PRE 8 C8',1.5)]" --max-step-size 2.0 --max-cycles 3 --no-preopt --no-endopt --out-json --out-dir test12 > test12.out 2>&1
python assert_release_result.py scan-optimizer test12 --expected-mode grad --expected-optimizer lbfgs >> test12.out 2>&1

# test13: scan2d
mlmm scan2d -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',1.4,1.8),('PRE 8 C1','PRE 8 C8',3.2,3.6)]" --max-step-size 0.4 --relax-max-cycles 100 --thresh gau_loose --out-dir test13 > test13.out 2>&1

# test14: scan3d
mlmm scan3d -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',1.4,1.8),('PRE 8 C1','PRE 8 C8',3.2,3.6),('PRE 8 C1','PRE 8 C7',1.4,1.8)]" --max-step-size 0.4 --relax-max-cycles 100 --thresh gau_loose --out-dir test14 > test14.out 2>&1

# test15: path-opt (gsm)
mlmm path-opt -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-nodes 5 --max-cycles-gsm 5 --thresh-gsm gau_loose --no-preopt --no-climb --out-dir test15 > test15.out 2>&1
grep -Fq "====== Growing String optimization ======" test15.out || { echo "[smoke] FAIL: path-opt GSM section heading missing" >&2; exit 1; }

# test16: path-opt (dmf)
mlmm path-opt -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --mep-mode dmf --max-cycles-dmf 3 --thresh-dmf middle --no-preopt --out-dir test16 > test16.out 2>&1

# test17: path-search
mlmm path-search -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-cycles-gsm 5 --out-dir test17 > test17.out 2>&1
grep -Fq "====== [seg_000_refine] GSM ======" test17.out || { echo "[smoke] FAIL: tagged recursive GSM section heading missing" >&2; exit 1; }

# test18: all (no tsopt/thermo/dft)
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge 'PRE:0' -q -1 -m 1 --no-refine-path --max-cycles-gsm 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test18 > test18.out 2>&1

# test19: tsopt (radius-hessian 0.0)
mlmm tsopt -i p_complex.pdb --parm p_complex.parm7 --model-pdb pocket_r.pdb -q -1 -m 1 --opt-mode grad --max-cycles 5 --radius-hessian 0.0 --active-dof-mode ml-only --thresh gau_loose --out-dir test19 > test19.out 2>&1

# test20: tsopt (radius-hessian 3.6)
mlmm tsopt -i p_complex.pdb --parm p_complex.parm7 --model-pdb pocket_r.pdb -q -1 -m 1 --opt-mode grad --max-cycles 5 --radius-hessian 3.6 --active-dof-mode ml-only --thresh gau_loose --out-dir test20 > test20.out 2>&1
python - <<'PY'
import re
from pathlib import Path

def initial_active_count(name: str) -> int:
    text = Path(name).read_text(encoding="utf-8")
    match = re.search(r"\[tsopt\] H_act=\d+ active_atoms=(\d+)", text)
    if match is None:
        raise SystemExit(f"[smoke] FAIL: {name} lacks the initial Hessian coverage record")
    return int(match.group(1))

ml_only = initial_active_count("test19.out")
expanded = initial_active_count("test20.out")
if expanded <= ml_only:
    raise SystemExit(
        "[smoke] FAIL: --radius-hessian 3.6 did not expand the Dimer Hessian "
        f"coverage ({expanded} <= {ml_only})"
    )
print(f"[smoke] PASS: radius-hessian coverage expanded {ml_only} -> {expanded} atoms")
PY

# test21: opt --dry-run
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode grad --dry-run --out-dir test21 > test21.out 2>&1

# test22: tsopt --dry-run
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --dry-run --out-dir test22 > test22.out 2>&1

# test23: freq --dry-run
mlmm freq -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --dry-run --out-dir test23 > test23.out 2>&1

# test24: scan --dry-run
mlmm scan -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',3.5),('PRE 8 C1','PRE 8 C8',1.5)]" --dry-run --out-dir test24 > test24.out 2>&1

# test25: dft --dry-run
mlmm dft -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --dry-run --out-dir test25 > test25.out 2>&1

# test26: path-search --dry-run
mlmm path-search -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --dry-run --out-dir test26 > test26.out 2>&1

# test27: irc --dry-run
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --dry-run --out-dir test27 > test27.out 2>&1

# test28: add-elem-info
mlmm add-elem-info -i r_complex.pdb -o r_complex_elem.pdb > test28.out 2>&1

# test29: trj2fig
mlmm trj2fig -i test4/optimization_trj.xyz -o test29.png > test29.out 2>&1

# test30: energy-diagram
mlmm energy-diagram -i "[0, 12.5, 4.3, 18.7, -1.2]" -o test30.png > test30.out 2>&1

# test31: oniom-export
python - <<'PY'
import parmed as pmd

parm = pmd.load_file("p_complex.parm7")
parm.cmaps[:] = []
parm.save("p_complex_nocmap.parm7", overwrite=True)
PY
mlmm oniom-export --parm p_complex_nocmap.parm7 -i r_complex_layered.pdb --model-pdb pocket_r.pdb -q -1 -m 1 -o test31.gjf > test31.out 2>&1
if mlmm oniom-export --parm p_complex.parm7 -i r_complex_layered.pdb --model-pdb pocket_r.pdb -q -1 -m 1 -o test31_cmap.gjf > test31_cmap_g16.out 2>&1; then
  echo "[smoke] FAIL test31: Gaussian export accepted a CMAP topology" >> test31_cmap_g16.out
  exit 1
fi
grep -Fq "CMAP" test31_cmap_g16.out
test ! -e test31_cmap.gjf
if mlmm oniom-export --parm p_complex.parm7 -i r_complex_layered.pdb --model-pdb pocket_r.pdb -q -1 -m 1 --mode orca --no-convert-orcaff -o test31_cmap.inp > test31_cmap_orca.out 2>&1; then
  echo "[smoke] FAIL test31: ORCA export accepted a CMAP topology" >> test31_cmap_orca.out
  exit 1
fi
grep -Fq "CMAP" test31_cmap_orca.out
test ! -e test31_cmap.inp

# --- Bond-summary, fix-altloc, oniom-import ---

# test32: bond-summary (two layered PDBs)
mlmm bond-summary -i r_complex_layered.pdb p_complex_layered.pdb > test32.out 2>&1

# test33: fix-altloc
mlmm fix-altloc -i r_complex.pdb -o r_complex_fixalt.pdb > test33.out 2>&1

# test34: oniom-import (Gaussian input → layered PDB)
mlmm oniom-import -i test31.gjf -o test34 > test34.out 2>&1

# --- refine-path ---

# test35: all (--refine-path)
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge 'PRE:0' -q -1 -m 1 --refine-path --max-cycles-gsm 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test35 > test35.out 2>&1

# test36: experimental MLIP/MM embedding remains accepted by the CLI.
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode grad --max-cycles 3 --thresh gau_loose --embedcharge --embedcharge-cutoff 6.0 --dry-run --out-dir test36 > test36.out 2>&1

# --- Opt-in TS and IRC methods ---

# test37: tsopt --opt-mode trim (Helgaker trust-region image-min; non-microiter)
mlmm tsopt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode trim --no-microiter --max-cycles 5 --thresh gau_loose --out-dir test37 > test37.out 2>&1

# test38: explicit RS-I-RFO (non-microiter); default hess/RS-P-RFO is covered by test8
mlmm tsopt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode rsirfo --no-microiter --max-cycles 5 --thresh gau_loose --out-dir test38 > test38.out 2>&1

# test39: irc --irc-pos-def (PSD-Hessian convergence guard)
mlmm irc -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-cycles 3 --irc-pos-def --out-dir test39 > test39.out 2>&1

# test40: opt --print-every 3 (diagnostic output throttle)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --max-cycles 5 --thresh gau_loose --print-every 3 --out-dir test40 > test40.out 2>&1

# --- Determinism gate ---

# test41: fixed-stack `all --deterministic` artifact comparison.
# On this smoke input, software stack, and reused MM topology, exact artifact
# drift is a regression for the tested stack. Topology generation is outside
# the flag's scope, so run b reuses run a's parm7.
det_args="-i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --no-refine-path --max-cycles-gsm 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --deterministic"
mlmm all $det_args --out-dir test41_a > test41_a.out 2>&1
mapfile -t test41_parms < <(find test41_a/mm_parm -maxdepth 1 -type f -name '*.parm7' -print)
if [[ "${#test41_parms[@]}" -ne 1 ]]; then
  echo "[smoke] FAIL test41: expected exactly one reusable parm7, found ${#test41_parms[@]}" >&2
  exit 1
fi
mlmm all $det_args --parm "${test41_parms[0]}" --out-dir test41_b > test41_b.out 2>&1
# Run b is given the parm via `--parm`, so by design it never runs the mm-parm
# stage and never writes `mm_parm/`. That directory is the gate's fixed INPUT,
# not its output, so comparing it would fail on the very workaround above.
find test41_a -type f \( -name "*.pdb" -o -name "*.xyz" \) -not -path '*/mm_parm/*' -printf '%P\n' | LC_ALL=C sort > test41_a.manifest
find test41_b -type f \( -name "*.pdb" -o -name "*.xyz" \) -not -path '*/mm_parm/*' -printf '%P\n' | LC_ALL=C sort > test41_b.manifest
if ! cmp -s test41_a.manifest test41_b.manifest; then
  echo "[smoke] FAIL test41: deterministic runs produced different file manifests" > test41.out
  comm -3 test41_a.manifest test41_b.manifest >> test41.out
  cat test41.out
  exit 1
fi
total=$(wc -l < test41_a.manifest)
if [ "$total" -eq 0 ]; then
  echo "[smoke] FAIL test41: deterministic gate found no PDB/XYZ artifacts" > test41.out
  cat test41.out
  exit 1
fi
drifted=0
while IFS= read -r rel; do
  if ! cmp -s "test41_a/$rel" "test41_b/$rel"; then
    drifted=$((drifted + 1))
    echo "DRIFT: $rel" >> test41.out
  fi
done < test41_a.manifest
echo "[det_check] compared $total PDB/XYZ files; $drifted differ" >> test41.out
if [ "$drifted" -ne 0 ]; then
  echo "[smoke] FAIL test41: --deterministic runs differ" >> test41.out
  cat test41.out
  exit 1
fi

# --- --coord-type CLI plumbing (throttled, fast) ---

# test42: `all --coord-type cart` — explicit cart (== default), verifies CLI plumbing.
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --coord-type cart --no-refine-path --max-cycles-gsm 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test42 > test42.out 2>&1

# test43: `all --coord-type dlc` — DLC propagated to the child opt / path-opt
# stages this run enables (it passes --no-tsopt).
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --coord-type dlc --no-refine-path --max-cycles-gsm 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test43 > test43.out 2>&1

# test44: `sp` (single-point ONIOM) — energy + forces.
mlmm sp -i r_complex_layered.pdb --real-parm7 p_complex.parm7 -q -1 -m 1 --out-dir test44 > test44.out 2>&1

# test45: `sp --hess` — energy + forces + ONIOM Hessian (default FiniteDifference).
mlmm sp -i r_complex_layered.pdb --real-parm7 p_complex.parm7 -q -1 -m 1 --hess --out-dir test45 > test45.out 2>&1

# --- Full-pipeline release-gate runs ---
# test46 exercises the unthrottled default `all` flow without a cycle cap.
# test47 is capped to cover the DLC path without requiring downstream
# convergence.

# test46: full `all` cart — default thresh, no max-cycles cap.
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --out-dir test46 > test46.out 2>&1

# test47: `all` dlc — verifies the DLC code-path lights up end-to-end.
# Capped at max-cycles 5 + thresh gau_loose + --no-tsopt/thermo/dft so this
# lane exercises DLC setup and trajectory handling without requiring a
# converged HEI. test46 keeps the no-cap default-behaviour check.
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --coord-type dlc --no-refine-path --max-cycles-gsm 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test47 > test47.out 2>&1

# --- Per-stage internal-coordinate code-path verification ---
# Each test is scoped at a 2-3 cycle cap (plus gau_loose where the stage
# needs it) so it exercises the
# coordinate paths without requiring convergence. Frequency analysis remains
# Cartesian because its PHVA contract consumes a Cartesian Hessian directly.

# test47a: opt --coord-type dlc
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --coord-type dlc --max-cycles 3 --thresh gau_loose --out-dir test47a_opt_dlc > test47a_opt_dlc.out 2>&1

# test47b: opt --opt-mode hess --coord-type dlc (microiter+DLC regression: ML internals, MM cart twin)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --coord-type dlc --max-cycles 3 --thresh gau_loose --out-dir test47b_opt_hess_dlc > test47b_opt_hess_dlc.out 2>&1

# test47c: opt --coord-type dlc with explicit frozen atoms
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --coord-type dlc --freeze-atoms "$MLMM_COMPLEX_FREEZE_ATOMS" --max-cycles 3 --thresh gau_loose --out-dir test47c_opt_freeze_dlc > test47c_opt_freeze_dlc.out 2>&1
python check_frozen_atoms.py r_complex_layered.pdb test47c_opt_freeze_dlc/final_geometry.pdb "$MLMM_COMPLEX_FREEZE_ATOMS" test47c >> test47c_opt_freeze_dlc.out 2>&1

# test47e: Hessian TS microiteration with DLC and frozen atoms. This is the
# partial-Cartesian-Hessian -> internal-coordinate handoff regression.
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --coord-type dlc --freeze-atoms "$MLMM_COMPLEX_FREEZE_ATOMS" --microiter --max-cycles 2 --thresh gau_loose --out-dir test47e_ts_hess_dlc > test47e_ts_hess_dlc.out 2>&1

# test47g: opt --coord-type redund
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --coord-type redund --max-cycles 3 --thresh gau_loose --out-dir test47g_opt_redund > test47g_opt_redund.out 2>&1

# test47k: opt --coord-type tric
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --coord-type tric --max-cycles 3 --thresh gau_loose --out-dir test47k_opt_tric > test47k_opt_tric.out 2>&1

# --- Multi-mode flag code-path verify (single-stage) ---

# test47m: opt --precision fp64 (UMA backend, alternate precision)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --precision fp64 --max-cycles 3 --thresh gau_loose --out-dir test47m_opt_fp64 > test47m_opt_fp64.out 2>&1

# test47n: opt --precision fp32 (explicit UMA fp32 dispatch alongside test47m fp64)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --precision fp32 --max-cycles 3 --thresh gau_loose --out-dir test47n_opt_fp32 > test47n_opt_fp32.out 2>&1

# test47p: opt --mm-backend openmm (alternate MM backend; analytical Hessian path → FD)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --mm-backend openmm --max-cycles 3 --thresh gau_loose --out-dir test47p_opt_openmm > test47p_opt_openmm.out 2>&1

# test47q: opt --link-atom-method fixed (legacy 1.09/1.01 Å placement)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --link-atom-method fixed --max-cycles 3 --thresh gau_loose --out-dir test47q_opt_linkfixed > test47q_opt_linkfixed.out 2>&1

# --- Non-default MLIP backend, full pipeline ---

# test48: full `all` with `--backend orb` — exercises the non-default MLIP backend.
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --backend orb --out-dir test48 > test48.out 2>&1
python assert_release_result.py provenance test48 --expected-backend orb --expected-model orb_v3_conservative_omol --expected-precision fp64 >> test48.out 2>&1

# ---- Subcommand-specific regression coverage ----
# test49: opt --mm-only (MM-only minimization; skips the MLIP component entirely)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --mm-only --opt-mode grad --max-cycles 3 --thresh gau_loose --out-dir test49_opt_mmonly > test49_opt_mmonly.out 2>&1

# test50: freq --active-dof-mode ml-only (alternate PHVA active-DOF subspace)
mlmm freq -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --active-dof-mode ml-only --max-write 5 --out-dir test50_freq_mlonly > test50_freq_mlonly.out 2>&1

# test51: freq --hessian-calc-mode Analytical (workflow analytical path)
mlmm freq -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --hessian-calc-mode Analytical --max-write 5 --out-dir test51_freq_anahess > test51_freq_anahess.out 2>&1

# test52: irc --hessian-calc-mode analytical (IRC initial Hessian path)
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --hessian-calc-mode analytical --max-cycles 2 --out-dir test52_irc_anahess > test52_irc_anahess.out 2>&1

# test53: irc --mm-backend openmm (MM Hessian via OpenMM finite-difference)
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --mm-backend openmm --max-cycles 2 --out-dir test53_irc_openmm > test53_irc_openmm.out 2>&1

# test54: irc --freeze-atoms (DOF-reduction / reduced-Hessian projection path)
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --freeze-atoms 1,2,3 --max-cycles 2 --out-dir test54_irc_freeze > test54_irc_freeze.out 2>&1

# test55: experimental DFT/MM embedding remains accepted by the CLI.
mlmm dft -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --func-basis 'hf/sto-3g' --grid-level 0 --conv-tol 1e-5 --max-cycle 40 --engine cpu --embedcharge --embedcharge-cutoff 8.0 --dry-run --out-dir test55_dft_embed > test55_dft_embed.out 2>&1

# test56: dft --link-atom-method fixed (legacy 1.09/1.01 Å link-atom placement)
mlmm dft -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --func-basis 'hf/sto-3g' --grid-level 0 --conv-tol 1e-5 --max-cycle 40 --engine cpu --link-atom-method fixed --out-dir test56_dft_linkfixed > test56_dft_linkfixed.out 2>&1

# test57: path-search --mep-mode dmf (Direct Max Flux vs GrowingString)
mlmm path-search -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --mep-mode dmf --max-cycles-dmf 3 --no-preopt --out-dir test57_psdmf > test57_psdmf.out 2>&1

# test58: all --scan-lists (single-PDB scan->path mode of `all`, distinct from the multi-PDB MEP branch)
mlmm all -i r_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',3.5),('PRE 8 C1','PRE 8 C8',1.5)]" --no-refine-path --max-cycles-gsm 3 --thresh gau_loose --no-tsopt --no-thermo --no-dft --out-dir test58_all_scan > test58_all_scan.out 2>&1

# --- refine-path opt-in (recursive path_search) extra coverage ---
# The `all` default is now single-pass path-opt; exercise the recursive
# path_search opt-in (`--refine-path`) in scan->path mode too (test35 already
# covers the multi-input endpoint MEP with --refine-path).

# test59: extract MULTI-INPUT via space-separated '-i a.pdb b.pdb' (one flag, two paths).
# Regression guard: a single -i with several space-separated paths must NOT drop the 2nd input.
# A single -o yields one multi-MODEL PDB, so both endpoints must appear (-> exactly 2 MODEL records).
mlmm extract -i r_complex.pdb p_complex.pdb -c PRE -r 5.0 --no-exclude-backbone --ligand-charge 'PRE:0' -o pocket_multi.pdb > test59_multi_extract.out 2>&1
n_models=$(grep -c '^MODEL' pocket_multi.pdb 2>/dev/null)
if [ "${n_models:-0}" -ne 2 ]; then
  echo "[extract-multi] test59: space-separated '-i a b' yielded ${n_models:-0} MODEL records (expected 2); the 2nd input was dropped" >> test59_multi_extract.out
  exit 1
fi

# test60: all --scan-lists --refine-path (single-PDB scan -> recursive path_search)
mlmm all -i r_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',3.5),('PRE 8 C1','PRE 8 C8',1.5)]" --refine-path --max-cycles-gsm 3 --thresh gau_loose --no-tsopt --no-thermo --no-dft --out-dir test60_rp_scan > test60_rp_scan.out 2>&1

# test61: --backend-model routing — a non-default model must reach the resolved
# runtime header. Dry-run avoids downloading the alternate model.
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --backend-model uma-m-1p1 --dry-run --out-dir test61_backend_model > test61_backend_model.out 2>&1
grep -Eq '^\[backend\] UMA \(UMA-M-1\.1 \(OMol\), fp32\)$' test61_backend_model.out || { echo "[smoke] FAIL test61: non-default backend model missing from resolved runtime summary" >> test61_backend_model.out; exit 1; }

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

# test62: a real ML/MM optimization crosses the mmCIF bridge and restores the
# original long chain and five-digit residue identifier in its public output.
mlmm opt -i r_complex_layered.cif --parm p_complex.parm7 -q -1 -m 1 --max-cycles 1 --thresh gau_loose --out-dir test62_opt_cif > test62_opt_cif.out 2>&1
test -s test62_opt_cif/final_geometry.pdb || { echo "[smoke] FAIL test62: final PDB missing" >> test62_opt_cif.out; exit 1; }
test -s test62_opt_cif/final_geometry.cif || { echo "[smoke] FAIL test62: final CIF missing" >> test62_opt_cif.out; exit 1; }
grep -q 'LONG_CHAIN' test62_opt_cif/final_geometry.cif || { echo "[smoke] FAIL test62: auth chain was not restored" >> test62_opt_cif.out; exit 1; }
grep -q '10001' test62_opt_cif/final_geometry.cif || { echo "[smoke] FAIL test62: auth residue number was not restored" >> test62_opt_cif.out; exit 1; }

# test63: exact chain/residue-name/residue-number selection remains stable
# after normalization, with both internal PDB and identifier-preserving CIF.
mlmm extract -i r_complex_layered.cif -c 'LONG_CHAIN:PRE:10001' -r 0.1 --no-add-linkh -o test63_model_from_cif.pdb -v 0 > test63_extract_cif.out 2>&1
test -s test63_model_from_cif.pdb || { echo "[smoke] FAIL test63: extracted PDB missing" >> test63_extract_cif.out; exit 1; }
test -s test63_model_from_cif.cif || { echo "[smoke] FAIL test63: extracted CIF missing" >> test63_extract_cif.out; exit 1; }

# test64: dump a partial Hessian with its active-DOF metadata.
# test65: consume that Hessian unchanged in IRC and verify never-stop at runtime.
# The dump must use IRC's own active-DOF basis (ML + MovableMM, i.e. freq's
# default). IRC has no --active-dof-mode/--hess-cutoff flag and always analyzes
# every movable atom, so an ml-only dump would be rejected as a basis mismatch;
# the partial nature is still exercised because --freeze-atoms keeps it < full.
mlmm freq -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --freeze-atoms 1,2,3 --max-write 1 --dump-hess test64_freq/hessian.npz --out-json --out-dir test64_freq > test64_freq.out 2>&1
test -s test64_freq/hessian.npz || { echo "[smoke] FAIL test64: dumped Hessian missing" >> test64_freq.out; exit 1; }
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --freeze-atoms 1,2,3 --read-hess test64_freq/hessian.npz --never-stop --config never_stop_config.yaml --max-cycles 2 --out-json --out-dir test65_irc_handoff > test65_irc_handoff.out 2>&1
python assert_release_result.py irc-handoff test65_irc_handoff --hessian-file test64_freq/hessian.npz >> test65_irc_handoff.out 2>&1
if grep -Fq '[irc] IRC stopped after only a few frames' test65_irc_handoff.out; then
  echo '[smoke] FAIL test65: cycle-cap completion was reported as early IRC termination' >> test65_irc_handoff.out
  exit 1
fi

# A same-size Hessian from a different geometry must be rejected before IRC.
python - <<'PY'
from pathlib import Path

lines = Path("p_complex_layered.pdb").read_text(encoding="utf-8").splitlines(True)
for index, line in enumerate(lines):
    if line.startswith(("ATOM  ", "HETATM")):
        x = float(line[30:38]) + 0.100
        lines[index] = line[:30] + f"{x:8.3f}" + line[38:]
        break
Path("test65_wrong_geometry.pdb").write_text("".join(lines), encoding="utf-8")
PY
rc=0
mlmm irc -i test65_wrong_geometry.pdb --parm p_complex.parm7 -q -1 -m 1 --freeze-atoms 1,2,3 --read-hess test64_freq/hessian.npz --max-cycles 1 --out-dir test65_wrong > test65_wrong.out 2>&1 || rc=$?
if [ "$rc" -eq 0 ] || ! grep -Eq 'coordinates do not match|PES identity does not match' test65_wrong.out; then
  echo "[smoke] FAIL test65: stale same-size Hessian was not rejected" >> test65_wrong.out
  exit 1
fi

# test66: standalone --ref-mode is an actual Cartesian mode vector. The
# all-workflow path-tangent handoff is exercised by required-positive test73.
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
np.savetxt("test66_reference_mode.txt", mode)
PY
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --ref-mode test66_reference_mode.txt --max-cycles 2 --thresh gau_loose --out-json --out-dir test66_ref_mode > test66_ref_mode.out 2>&1
python assert_release_result.py tsopt-reference test66_ref_mode >> test66_ref_mode.out 2>&1

# test67: YAML backend-model/precision settings reach a real calculation while
# an explicit CLI max-cycles value retains precedence.
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --config runtime_config.yaml --max-cycles 1 --out-json --out-dir test67_config > test67_config.out 2>&1
python assert_release_result.py opt-config test67_config --expected-model uma-s-1p2 --expected-max-cycles 1 --expected-precision fp64 --expected-link-atom-method fixed --expected-thresh gau_loose >> test67_config.out 2>&1

# test68: a user ASE calculator with only energy/forces supports the ML-region
# finite-difference Hessian path.
mlmm sp -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --calc-file harmonic_calc.py --hess --hessian-calc-mode FiniteDifference --out-json --out-dir test68_custom > test68_custom.out 2>&1
python assert_release_result.py sp-hessian test68_custom >> test68_custom.out 2>&1

# test69: an explicit analytical Hessian request cannot silently fall back when
# predictor workers remove the autograd model.
rc=0
mlmm sp -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --hess --hessian-calc-mode Analytical --workers 2 --out-dir test69_workers > test69_workers.out 2>&1 || rc=$?
if [ "$rc" -ne 1 ] || ! grep -Fq "Analytical Hessian cannot be combined with workers>1: the parallel predictor exposes no autograd model. Use workers=1 or select hessian_calc_mode='FiniteDifference'." test69_workers.out || [ -e test69_workers/hessian.npy ]; then
  echo "[smoke] FAIL test69: Analytical + workers>1 was not rejected exactly" >> test69_workers.out
  exit 1
fi

# test70: ORCA QM/MM export/import round-trip preserves coordinates and all
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
Path("test70_three_layer.pdb").write_text("".join(lines), encoding="utf-8")
PY
mlmm oniom-export --parm p_complex_nocmap.parm7 -i test70_three_layer.pdb --model-pdb pocket_r.pdb -q -1 -m 1 --mode orca --no-convert-orcaff -o test70_orca.inp > test70_orca_export.out 2>&1
for token in '! QMMM' 'QMAtoms {' 'ActiveAtoms {' 'Charge_Total -1' '* xyz -1 1'; do
  grep -Fq "$token" test70_orca.inp || { echo "[smoke] FAIL test70: ORCA input missing $token" >> test70_orca_export.out; exit 1; }
done
mlmm oniom-import -i test70_orca.inp --ref-pdb test70_three_layer.pdb -o test70_orca_import > test70_orca_import.out 2>&1
test -s test70_orca_import.xyz || { echo "[smoke] FAIL test70: restored XYZ missing" >> test70_orca_import.out; exit 1; }
test -s test70_orca_import_layered.pdb || { echo "[smoke] FAIL test70: restored layered PDB missing" >> test70_orca_import.out; exit 1; }
sed -n '2p' test70_orca_import.xyz | grep -Fq 'mode=orca' || { echo "[smoke] FAIL test70: restored XYZ omits ORCA provenance" >> test70_orca_import.out; exit 1; }
sed -n '2p' test70_orca_import.xyz | grep -Fq 'q=-1 m=1' || { echo "[smoke] FAIL test70: restored XYZ charge/multiplicity mismatch" >> test70_orca_import.out; exit 1; }
python assert_orca_roundtrip.py test70_three_layer.pdb test70_orca_import_layered.pdb >> test70_orca_import.out 2>&1

# test71: force a known higher-order candidate through the actual flatten
# branch. The checker requires n_imag>1 before flattening, an executed RS-P-RFO
# flatten iteration, and no increase in saddle order.
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --no-microiter --flatten --config flatten_branch_config.yaml --thresh gau_loose --out-json --out-dir test71_flatten > test71_flatten.out 2>&1
python assert_flatten_branch.py test71_flatten.out test71_flatten/result.json >> test71_flatten.out 2>&1

# test72: numerical analytical-vs-FD agreement for every backend installed in the
# default strict environment. MACE/AIMNet2 use this same required wrapper in
# their dependency-isolated cluster environments.
bash run_backend_hessian.sh uma orb > test72_backend_hessian.out 2>&1

# test73: required positive MEP -> TSopt -> IRC -> thermo -> DFT handoff.
# Endpoint and GSM thresholds are pinned independently for this positive lane.
# The long lane runs last with its dependent manual-topology reuse check.
mlmm all \
    -i test73_r_complex.pdb test73_p_complex.pdb \
    -c PRE \
    -r 4.0 \
    --ligand-charge 'PRE:0' \
    -q -1 \
    -m 1 \
    --deterministic \
    --no-refine-path \
    --thresh gau \
    --thresh-gsm gau \
    --thresh-post baker \
    --tsopt \
    --thermo \
    --dft \
    --flatten \
    --irc-never-stop \
    --irc-max-cycles 3 \
    --dft-func-basis 'hf/sto-3g' \
    --dft-grid-level 0 \
    --dft-conv-tol 1e-5 \
    --dft-max-cycle 40 \
    --dft-engine cpu \
    --out-dir test73 \
    > test73.out 2>&1
python assert_release_result.py all test73 --require-thermo --require-dft >> test73.out 2>&1

# test74: all (manual --parm + --model-pdb override, reuse test73 outputs)
mapfile -t test73_parms < <(find test73/mm_parm -maxdepth 1 -type f -name '*.parm7' -print)
if [[ "${#test73_parms[@]}" -ne 1 ]]; then
  echo "[smoke] FAIL test74: expected exactly one reusable test73 parm7, found ${#test73_parms[@]}" >&2
  exit 1
fi
mlmm all -i test73/layered/test73_r_complex_layered.pdb test73/layered/test73_p_complex_layered.pdb --parm "${test73_parms[0]}" --model-pdb test73/ml_region.pdb -q -1 -m 1 --no-refine-path --max-cycles-gsm 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test74 > test74.out 2>&1

echo "[smoke] PASS: required GPU, ML/MM, Hessian-handoff, and structure-I/O lane completed with zero skips."
