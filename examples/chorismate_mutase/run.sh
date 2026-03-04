#!/usr/bin/env bash
# Chorismate mutase — Claisen rearrangement
# Full 10-stage ML/MM workflow using the mlmm CLI.
#
# Prerequisites:
#   parm/complex.pdb, parm/complex.parm7, parm/complex.rst7, parm/ml_region.pdb
#
# Intermediate geometries are pre-computed in coord/ so that each step
# can also be run independently.

set -euo pipefail

# (0) Define layers
mlmm define-layer -i parm/complex.pdb --model-pdb parm/ml_region.pdb --radius-freeze 10.0 -o complex_layered.pdb

# (i) Initial geometry optimization
mlmm opt -i complex_layered.pdb --parm parm/complex.parm7 -q 0 -m 1 --opt-mode grad --out-dir dump/opt1

# (ii) 1-D bond-length scan (atom indices 5686 and 5696, 0-based)
mlmm scan -i coord/1_opt_final_geometry.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q 0 -m 1 --scan-lists "[(5686,5696,4.0)]" --max-step-size 0.05 --out-dir dump/scan

# (iii) Reactant / product optimization
mlmm opt -i coord/2_scan_reac.pdb --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q 0 -m 1 --opt-mode grad --out-dir dump/opt2
mlmm opt -i coord/2_scan_prod.pdb --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q 0 -m 1 --opt-mode grad --out-dir dump/opt3

# (iv) Growing String Method
mlmm path-opt -i coord/3_opt_final_geometry.xyz coord/4_opt_final_geometry.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q 0 -m 1 --max-nodes 20 --climb True --out-dir dump/gs

# (v) Transition state refinement
mlmm tsopt -i coord/5_gs_peak.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q 0 -m 1 --opt-mode grad --out-dir dump/tsopt

# (vi) IRC propagation
mlmm irc -i coord/6_tsopt_final_geometry.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q 0 -m 1 --out-dir dump/irc

# (vii) Endpoint relaxation
mlmm opt -i coord/7_irc_forward_last.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q 0 -m 1 --opt-mode grad --out-dir dump/opt4
mlmm opt -i coord/7_irc_backward_last.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q 0 -m 1 --opt-mode grad --out-dir dump/opt5

# (viii) Vibrational analysis & thermochemistry
mlmm freq -i coord/6_tsopt_final_geometry.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q 0 -m 1 --out-dir dump/freq_ts
mlmm freq -i coord/8_opt_final_geometry.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q 0 -m 1 --out-dir dump/freq_prod
mlmm freq -i coord/9_opt_final_geometry.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q 0 -m 1 --out-dir dump/freq_reac

echo "Chorismate mutase workflow complete."
