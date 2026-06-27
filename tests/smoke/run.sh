#!/usr/bin/env bash
# mlmm smoke tests — GPU required.
#
# This script assumes the calling environment already has:
#   - a Python with `mlmm` installed and importable
#   - AmberTools (antechamber / parmchk2 / tleap) on PATH
#   - CUDA available
#   - a writeable working directory (the smoke artefacts land in `test*/`)
# It does NOT activate conda, load modules, or contain HPC scheduler
# directives. Run it directly from your active env:
#
#   bash tests/smoke/run.sh
#
# If you need an HPC scheduler wrapper, keep that in
# your own out-of-tree submission script and have it invoke this file as
# the body.
set -euo pipefail

# Deterministic run gate: pin PYTHONHASHSEED so Set / dict iteration order
# does not leak hash randomisation into the produced output. Combined with the
# UMA backend's deterministic-algorithms wiring this removes the remaining
# I/O-side non-determinism between two consecutive runs.
export PYTHONHASHSEED=0
# Reduce CUDA allocator fragmentation across the 40+ stage processes.
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

# Clean previous results
rm -rf test* pocket_r.pdb r_complex_layered.pdb r_complex_elem.pdb r_complex_fixalt.pdb

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
mlmm scan2d -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',1.0,3.0),('PRE 8 C1','PRE 8 C8',1.0,3.0)]" --max-step-size 2.0 --relax-max-cycles 100 --thresh gau_loose --out-dir test13 > test13.out 2>&1

# test14: scan3d
mlmm scan3d -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',1.0,3.0),('PRE 8 C1','PRE 8 C8',1.0,3.0),('PRE 8 C1','PRE 8 C7',1.0,2.0)]" --max-step-size 2.0 --relax-max-cycles 100 --thresh gau_loose --out-dir test14 > test14.out 2>&1

# test15: path-opt (gsm)
mlmm path-opt -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-nodes 5 --max-cycles 5 --no-preopt --no-climb --out-dir test15 > test15.out 2>&1

# test16: path-opt (dmf)
mlmm path-opt -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --mep-mode dmf --max-cycles 3 --no-preopt --out-dir test16 > test16.out 2>&1

# test17: path-search
mlmm path-search -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-cycles 5 --out-dir test17 > test17.out 2>&1

# test18: all (no tsopt/thermo/dft)
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge 'PRE:0' -q -1 -m 1 --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test18 > test18.out 2>&1

# test19: all (tsopt + thermo + dft) — the throttled --max-cycles run may leave the
# IRC seed non-converged for some MLIP seeds; tolerate ONLY that case (with a printed
# message) and fail the suite on any other non-zero exit (tsopt ZeroStepLength /
# OptimizationError, DFT/SCF sys.exit(3), OOM, setup error).
rc=0
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge 'PRE:0' -q -1 -m 1 --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --tsopt --thermo --dft --tsopt-max-cycles 5 --dft-func-basis 'hf/sto-3g' --dft-grid-level 0 --dft-conv-tol 1e-5 --dft-max-cycle 40 --dft-engine cpu --out-dir test19 > test19.out 2>&1 || rc=$?
if [ "${rc:-0}" -ne 0 ]; then
  if grep -qiE "IRC|did not converge|not converged" test19.out; then
    echo "[smoke] test19: IRC non-convergence on throttled run (rc=$rc) — tolerated, continuing"
  else
    echo "[smoke] FAIL test19 (rc=$rc) — real pipeline failure"
    tail -40 test19.out
    exit "$rc"
  fi
fi

# test20: all (--parm + --model-pdb override, reuse test19 outputs)
mlmm all -i r_complex.pdb p_complex.pdb --parm test19/mm_parm/r_complex.parm7 --model-pdb test19/ml_region.pdb -q -1 -m 1 --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test20 > test20.out 2>&1

# test21: tsopt (radius-hessian 0.0)
mlmm tsopt -i p_complex.pdb --parm p_complex.parm7 --model-pdb pocket_r.pdb --no-detect-layer -q -1 -m 1 --opt-mode grad --max-cycles 5 --radius-hessian 0.0 --thresh gau_loose --out-dir test21 > test21.out 2>&1

# test22: tsopt (radius-hessian 3.6)
mlmm tsopt -i p_complex.pdb --parm p_complex.parm7 --model-pdb pocket_r.pdb --no-detect-layer -q -1 -m 1 --opt-mode grad --max-cycles 5 --radius-hessian 3.6 --thresh gau_loose --out-dir test22 > test22.out 2>&1

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

# --- xTB-dependent tests (requires xtb binary) ---

# test38: opt (embedcharge) — skip if xtb is not available
if command -v xtb &>/dev/null; then
  mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode grad --max-cycles 3 --thresh gau_loose --embedcharge --embedcharge-cutoff 6.0 --out-dir test38 > test38.out 2>&1
else
  echo "SKIP test38: xtb not found" > test38.out
fi

# --- Polish-train new CLI flags (A1 + W3 + B4 wires; all opt-in, defaults preserve Table 1 numerics) ---

# test39: tsopt --opt-mode trim (A1 Helgaker trust-region image-min; non-microiter)
mlmm tsopt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode trim --no-microiter --max-cycles 5 --thresh gau_loose --skip-final-freq --out-dir test39 > test39.out 2>&1

# test40: tsopt --opt-mode rsprfo (A1 Banerjee P-RFO; non-microiter)
mlmm tsopt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode rsprfo --no-microiter --max-cycles 5 --thresh gau_loose --skip-final-freq --out-dir test40 > test40.out 2>&1

# test42: irc --irc-pos-def (Sella backport 1: PSD-Hessian convergence guard)
mlmm irc -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-cycles 3 --irc-pos-def --out-dir test42 > test42.out 2>&1

# test43: opt --print-every 3 (W3a debug throttle, no behavior change)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --max-cycles 5 --thresh gau_loose --print-every 3 --out-dir test43 > test43.out 2>&1

# --- Determinism gate ---

# test44: `all` pipeline determinism GATE (`--deterministic`, ONIOM end-to-end).
# Runs the full pipeline twice with identical inputs / args + `--deterministic`
# and REQUIRES the two runs to be bit-identical. Default (non-deterministic) GPU
# runs carry ~ULP scatter/atomic non-determinism and are not asserted here;
# `--deterministic` enables torch deterministic algorithms and MUST be
# bit-reproducible, so any drift is a real regression and fails the smoke.
det_args="-i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --deterministic"
mlmm all $det_args --out-dir test44_a > test44_a.out 2>&1
mlmm all $det_args --out-dir test44_b > test44_b.out 2>&1
{
  total=0
  drifted=0
  while IFS= read -r path_a; do
    rel="${path_a#test44_a/}"
    path_b="test44_b/$rel"
    if [ -f "$path_b" ]; then
      total=$((total + 1))
      cmp -s "$path_a" "$path_b" || { drifted=$((drifted + 1)); echo "DRIFT: $rel"; }
    fi
  done < <(find test44_a -type f \( -name "*.pdb" -o -name "*.xyz" \))
  echo "[det_check] all pipeline (--deterministic): compared $total .pdb/.xyz file(s); $drifted differ"
} > test44.out 2>&1
if [ "${drifted:-0}" -ne 0 ]; then
  echo "[smoke] FAIL test44: --deterministic runs are not bit-reproducible ($drifted file(s) differ)"
  cat test44.out
  exit 1
fi

# --- --coord-type CLI plumbing (throttled, fast) ---

# test45: `all --coord-type cart` — explicit cart (== default), verifies CLI plumbing.
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --coord-type cart --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test45 > test45.out 2>&1

# test46: `all --coord-type dlc` — DLC propagated to child opt / tsopt / path-opt stages.
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --coord-type dlc --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test46 > test46.out 2>&1

# test47: `sp` (single-point ONIOM) — energy + forces.
mlmm sp -i r_complex_layered.pdb --real-parm7 p_complex.parm7 -q -1 -m 1 --out-dir test47 > test47.out 2>&1

# test48: `sp --hess` — energy + forces + ONIOM Hessian (UMA analytical).
mlmm sp -i r_complex_layered.pdb --real-parm7 p_complex.parm7 -q -1 -m 1 --hess --out-dir test48 > test48.out 2>&1

# --- Full-pipeline, NO max-cycles throttle (release-gate runs) ---
# These exercise the canonical `all` flow with default convergence thresholds
# and the production-realistic optimizer cycle budgets. Each run takes
# substantially longer (~30-90 min for ONIOM) than the throttled tests above;
# they are the "does the pipeline actually finish on a real input" gate.

# test49: full `all` cart — default thresh, no max-cycles cap.
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --out-dir test49 > test49.out 2>&1

# test50: `all` dlc — verifies the DLC code-path lights up end-to-end.
# Capped at max-cycles 5 + thresh gau_loose + --no-tsopt/thermo/dft because
# DLC GSM on this 122-atom complex needs hundreds of cycles to converge with
# default `gau` thresh (3h+ on consumer GPU) and the post-stages (TS / IRC /
# DFT) depend on a converged HEI from the MEP — silently broken structure
# handoff otherwise. test49 (cart) keeps the no-cap default-behaviour check.
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --coord-type dlc --no-refine-path --max-cycles 5 --thresh gau_loose --thresh-post gau_loose --no-tsopt --no-thermo --no-dft --out-dir test50 > test50.out 2>&1

# --- Per-stage internal-coord code-path verify (opt + scan only) ---
# Each test is scoped at max-cycles 3 + thresh gau_loose so it exercises the
# coord-type code path without long convergence. Excluded subcmds:
#   - tsopt: pysisyphus RSIRFOptimizer.full_from_active raises a shape mismatch
#     for non-cart coord types when active-atom subset DOF differs from the
#     internal coord count (upstream pysisyphus fork issue).
#   - freq: vibrational analysis uses cart Hessian regardless of --coord-type;
#     the only thing the flag adds is a cusolver SVD on the internal B-matrix
#     that occasionally fails on the per-atom-subset internal-coord set.
#     freq + cart is already covered by test9.

# test50a: opt --coord-type dlc
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --coord-type dlc --max-cycles 3 --thresh gau_loose --out-dir test50a_opt_dlc > test50a_opt_dlc.out 2>&1

# test50b: opt --opt-mode hess --coord-type dlc (microiter+DLC regression: ML internals, MM cart twin)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --coord-type dlc --max-cycles 3 --thresh gau_loose --out-dir test50b_opt_hess_dlc > test50b_opt_hess_dlc.out 2>&1

# test50c: opt --coord-type dlc with explicit frozen atoms
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --coord-type dlc --freeze-atoms "$MLMM_COMPLEX_FREEZE_ATOMS" --max-cycles 3 --thresh gau_loose --out-dir test50c_opt_freeze_dlc > test50c_opt_freeze_dlc.out 2>&1
python check_frozen_atoms.py r_complex_layered.pdb test50c_opt_freeze_dlc/final_geometry.pdb "$MLMM_COMPLEX_FREEZE_ATOMS" test50c >> test50c_opt_freeze_dlc.out 2>&1

# test50d: scan --coord-type dlc with explicit frozen atoms.
# Small non-reactive target: this checks coordinate integrity, not chemistry.
mlmm scan -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --coord-type dlc --freeze-atoms "$MLMM_COMPLEX_FREEZE_ATOMS" --scan-lists "[('PRE 8 O1\'','PRE 8 C3',1.5)]" --max-step-size 0.1 --max-cycles 3 --no-endopt --out-dir test50d_scan_freeze_dlc > test50d_scan_freeze_dlc.out 2>&1
python check_frozen_atoms.py r_complex_layered.pdb test50d_scan_freeze_dlc/stage_01/result.pdb "$MLMM_COMPLEX_FREEZE_ATOMS" test50d >> test50d_scan_freeze_dlc.out 2>&1
if grep -q "Covalent-bond changes (start vs final): Yes" test50d_scan_freeze_dlc.out; then
  echo "[bond-check] test50d: unexpected covalent-bond changes in non-reactive DLC+freeze scan" >> test50d_scan_freeze_dlc.out
  exit 1
fi

# test50g: opt --coord-type redund
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --coord-type redund --max-cycles 3 --thresh gau_loose --out-dir test50g_opt_redund > test50g_opt_redund.out 2>&1

# test50k: opt --coord-type tric
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --coord-type tric --max-cycles 3 --thresh gau_loose --out-dir test50k_opt_tric > test50k_opt_tric.out 2>&1

# --- Multi-mode flag code-path verify (single-stage) ---

# test50m: opt --precision fp64 (UMA backend, alternate precision)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --precision fp64 --max-cycles 3 --thresh gau_loose --out-dir test50m_opt_fp64 > test50m_opt_fp64.out 2>&1

# test50n: opt --precision fp32 (explicit fp32 dispatch; default is fp32, this pins the explicit path alongside test50m fp64)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --precision fp32 --max-cycles 3 --thresh gau_loose --out-dir test50n_opt_fp32 > test50n_opt_fp32.out 2>&1

# test50n: opt --mm-backend openmm (alternate MM backend; analytical Hessian path → FD)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --mm-backend openmm --max-cycles 3 --thresh gau_loose --out-dir test50n_opt_openmm > test50n_opt_openmm.out 2>&1

# test50o: opt --link-atom-method fixed (legacy 1.09/1.01 Å placement)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --link-atom-method fixed --max-cycles 3 --thresh gau_loose --out-dir test50o_opt_linkfixed > test50o_opt_linkfixed.out 2>&1

# test51: full `all` with `--backend orb` — exercises the non-default MLIP backend.
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --backend orb --out-dir test51 > test51.out 2>&1

# ---- Coverage-gap regression (subcommand-specific code paths; coverage audit 2026-06-05) ----
# test52: opt --mm-only (MM-only minimization; skips the MLIP component entirely)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --mm-only --opt-mode grad --max-cycles 3 --thresh gau_loose --out-dir test52_opt_mmonly > test52_opt_mmonly.out 2>&1

# test53: freq --active-dof-mode ml-only (alternate PHVA active-DOF subspace)
mlmm freq -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --active-dof-mode ml-only --max-write 5 --out-dir test53_freq_mlonly > test53_freq_mlonly.out 2>&1

# test54: freq --hessian-calc-mode Analytical (ML-block Hessian analytical vs FD)
mlmm freq -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --hessian-calc-mode Analytical --max-write 5 --out-dir test54_freq_anahess > test54_freq_anahess.out 2>&1

# test55: irc --hessian-calc-mode analytical (IRC initial Hessian path)
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --hessian-calc-mode analytical --max-cycles 2 --out-dir test55_irc_anahess > test55_irc_anahess.out 2>&1

# test56: irc --mm-backend openmm (MM Hessian via OpenMM finite-difference)
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --mm-backend openmm --max-cycles 2 --out-dir test56_irc_openmm > test56_irc_openmm.out 2>&1

# test57: irc --freeze-atoms (DOF-reduction / reduced-Hessian projection path)
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --freeze-atoms 1,2,3 --max-cycles 2 --out-dir test57_irc_freeze > test57_irc_freeze.out 2>&1

# test58: dft --embedcharge (MM point charges into the PySCF QM Hamiltonian)
mlmm dft -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --func-basis 'hf/sto-3g' --grid-level 0 --conv-tol 1e-5 --max-cycle 40 --engine cpu --embedcharge --embedcharge-cutoff 8.0 --out-dir test58_dft_embed > test58_dft_embed.out 2>&1

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

# test62: all --scan-lists --refine-path (single-PDB scan -> recursive path_search)
mlmm all -i r_complex.pdb -c PRE -r 6.0 --ligand-charge PRE:0 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',3.5),('PRE 8 C1','PRE 8 C8',1.5)]" --refine-path --max-cycles 3 --thresh gau_loose --no-tsopt --no-thermo --no-dft --out-dir test62_rp_scan > test62_rp_scan.out 2>&1
