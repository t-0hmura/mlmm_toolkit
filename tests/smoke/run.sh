#!/usr/bin/env bash
set -euo pipefail
# mlmm smoke tests — GPU required

# test1: extract
mlmm extract -i r_complex.pdb -c PRE -r 5.0 --exclude-backbone False --ligand-charge 'PRE:0' -o pocket_r.pdb --verbose > test1.out 2>&1

# test2: define-layer
mlmm define-layer -i r_complex.pdb --model-pdb pocket_r.pdb --radius-freeze 8.0 -o r_complex_layered.pdb --verbose > test2.out 2>&1

# test3: mm-parm
mlmm mm-parm -i r_complex.pdb --ligand-charge 'PRE:0' --verbose > test3.out 2>&1

# test4: opt (grad)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode grad --max-cycles 5 --thresh gau_loose --dump --out-dir test4 --verbose > test4.out 2>&1

# test5: opt (hess)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --max-cycles 3 --thresh gau_loose --out-dir test5 --verbose > test5.out 2>&1

# test6: opt (hess, microiter)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --microiter --max-cycles 2 --thresh gau_loose --out-dir test6 --verbose > test6.out 2>&1

# test7: tsopt (grad / dimer)
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode grad --max-cycles 100 --thresh gau --out-dir test7 --verbose > test7.out 2>&1

# test8: tsopt (hess / rsirfo)
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --max-cycles 5 --thresh gau --out-dir test8 --verbose > test8.out 2>&1

# test9: freq
mlmm freq -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --out-dir test9 --verbose > test9.out 2>&1

# test10: irc
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-cycles 3 --out-dir test10 --verbose > test10.out 2>&1

# test11: dft (hf/sto-3g)
mlmm dft -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --func-basis 'hf/sto-3g' --grid-level 0 --conv-tol 1e-5 --max-cycle 40 --out-dir test11 --verbose > test11.out 2>&1

# test12: scan (1D)
mlmm scan -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',2.0)]" --max-step-size 2.0 --max-cycles 3 --no-preopt --no-endopt --out-dir test12 --verbose > test12.out 2>&1

# test13: scan2d
mlmm scan2d -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',1.0,3.0),('PRE 8 C1','PRE 8 C8',1.0,3.0)]" --max-step-size 2.0 --relax-max-cycles 100 --thresh gau_loose --out-dir test13 --verbose > test13.out 2>&1

# test14: scan3d
mlmm scan3d -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',1.0,3.0),('PRE 8 C1','PRE 8 C8',1.0,3.0),('PRE 8 C1','PRE 8 C7',1.0,2.0)]" --max-step-size 2.0 --relax-max-cycles 100 --thresh gau_loose --out-dir test14 --verbose > test14.out 2>&1

# test15: path-opt (gsm)
mlmm path-opt -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-nodes 5 --max-cycles 5 --no-preopt --no-climb --out-dir test15 --verbose > test15.out 2>&1

# test16: path-opt (dmf)
mlmm path-opt -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --mep-mode dmf --max-cycles 3 --no-preopt --out-dir test16 --verbose > test16.out 2>&1

# test17: path-search
mlmm path-search -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-cycles 5 --out-dir test17 --verbose > test17.out 2>&1

# test18: all (no tsopt/thermo/dft)
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge 'PRE:0' -q -1 -m 1 --max-cycles 5 --thresh gau_loose --no-tsopt --no-thermo --no-dft --out-dir test18 --verbose > test18.out 2>&1

# test19: all (tsopt + thermo + dft)
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge 'PRE:0' -q -1 -m 1 --max-cycles 5 --thresh gau_loose --tsopt --thermo --dft --tsopt-max-cycles 5 --dft-func-basis 'hf/sto-3g' --dft-grid-level 0 --dft-conv-tol 1e-5 --dft-max-cycle 40 --out-dir test19 --verbose > test19.out 2>&1

# test20: all (--parm + --model-pdb override, reuse test19 outputs)
mlmm all -i r_complex.pdb p_complex.pdb --parm test19/mm_parm/r_complex.parm7 --model-pdb test19/ml_region.pdb -q -1 -m 1 --max-cycles 5 --thresh gau_loose --no-tsopt --no-thermo --no-dft --out-dir test20 --verbose > test20.out 2>&1

# test21: tsopt (radius-hessian 0.0)
mlmm tsopt -i p_complex.pdb --parm p_complex.parm7 --model-pdb pocket_r.pdb --no-detect-layer -q -1 -m 1 --opt-mode grad --max-cycles 5 --radius-hessian 0.0 --thresh gau_loose --out-dir test21 --verbose > test21.out 2>&1

# test22: tsopt (radius-hessian 3.6)
mlmm tsopt -i p_complex.pdb --parm p_complex.parm7 --model-pdb pocket_r.pdb --no-detect-layer -q -1 -m 1 --opt-mode grad --max-cycles 5 --radius-hessian 3.6 --thresh gau_loose --out-dir test22 --verbose > test22.out 2>&1

# test23: opt --dry-run
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode grad --dry-run --out-dir test23 --verbose > test23.out 2>&1

# test24: tsopt --dry-run
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --dry-run --out-dir test24 --verbose > test24.out 2>&1

# test25: freq --dry-run
mlmm freq -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --dry-run --out-dir test25 --verbose > test25.out 2>&1

# test26: scan --dry-run
mlmm scan -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',2.0)]" --dry-run --out-dir test26 --verbose > test26.out 2>&1

# test27: dft --dry-run
mlmm dft -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --dry-run --out-dir test27 --verbose > test27.out 2>&1

# test28: path-search --dry-run
mlmm path-search -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --dry-run --out-dir test28 --verbose > test28.out 2>&1

# test29: irc --dry-run
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --dry-run --out-dir test29 --verbose > test29.out 2>&1

# test30: add-elem-info
mlmm add-elem-info -i r_complex.pdb -o r_complex_elem.pdb --verbose > test30.out 2>&1

# test31: trj2fig
mlmm trj2fig -i test4/optimization_trj.xyz -o test31.png --verbose > test31.out 2>&1

# test32: energy-diagram
mlmm energy-diagram -i "[0, 12.5, 4.3, 18.7, -1.2]" -o test32.png --verbose > test32.out 2>&1

# test33: oniom-export
mlmm oniom-export --parm p_complex.parm7 -i r_complex_layered.pdb --model-pdb pocket_r.pdb -q -1 -m 1 -o test33.gjf --verbose > test33.out 2>&1

# --- xTB-dependent tests (requires xtb binary) ---

# test34: opt (embedcharge)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode grad --max-cycles 3 --thresh gau_loose --embedcharge --embedcharge-cutoff 6.0 --out-dir test34 --verbose > test34.out 2>&1
