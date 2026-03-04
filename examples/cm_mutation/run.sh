#!/usr/bin/env bash
# Chorismate mutase Arg90Ala mutant — barrier comparison
# Steps (iii)-(viii): reuses scan results from the wild-type system.
#
# Prerequisites:
#   parm/complex.pdb, parm/complex.parm7, parm/complex.rst7, parm/ml_region.pdb
#   coord/2_scan_reac.pdb, coord/2_scan_prod.pdb  (from wild-type scan)

set -euo pipefail

# (iii) Reactant / product optimization
mlmm opt -i coord/2_scan_reac.pdb --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q -1 -m 1 --opt-mode grad --out-dir dump/opt2
mlmm opt -i coord/2_scan_prod.pdb --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q -1 -m 1 --opt-mode grad --out-dir dump/opt3

# (iv) Growing String Method
mlmm path-opt -i coord/3_opt_final_geometry.xyz coord/4_opt_final_geometry.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q -1 -m 1 --max-nodes 20 --climb True --out-dir dump/gs

# (v) Transition state refinement
mlmm tsopt -i coord/5_gs_peak.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q -1 -m 1 --opt-mode grad --out-dir dump/tsopt

# (vi) IRC propagation
mlmm irc -i coord/6_tsopt_final_geometry.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q -1 -m 1 --out-dir dump/irc

# (vii) Endpoint relaxation
mlmm opt -i coord/7_irc_forward_last.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q -1 -m 1 --opt-mode grad --out-dir dump/opt4
mlmm opt -i coord/7_irc_backward_last.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q -1 -m 1 --opt-mode grad --out-dir dump/opt5

# (viii) Vibrational analysis & thermochemistry
mlmm freq -i coord/6_tsopt_final_geometry.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q -1 -m 1 --out-dir dump/freq_ts
mlmm freq -i coord/8_opt_final_geometry.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q -1 -m 1 --out-dir dump/freq_prod
mlmm freq -i coord/9_opt_final_geometry.xyz --parm parm/complex.parm7 --model-pdb parm/ml_region.pdb -q -1 -m 1 --out-dir dump/freq_reac

echo "CM mutation workflow complete."
