# test1: extract
mlmm extract -i r_complex.pdb -c PRE -r 5.0 --exclude-backbone False --ligand-charge 'PRE:0' -o pocket_r_smoke.pdb > test1.out 2>&1

# test2: define-layer
mlmm define-layer -i r_complex.pdb --model-pdb pocket_r_smoke.pdb --radius-freeze 8.0 -o r_complex_layered_smoke.pdb > test2.out 2>&1

# test3: path-opt (gsm)
mlmm path-opt -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-nodes 5 --max-cycles 10 --preopt False --climb False --out-dir test3 > test3.out 2>&1

# test4: path-opt (dmf)
mlmm path-opt -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --mep-mode dmf --max-cycles 5 --preopt False --out-dir test4 > test4.out 2>&1

# test5: path-search
mlmm path-search -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-cycles 10 --out-dir test5 > test5.out 2>&1

# test6: opt (grad / lbfgs)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode grad --max-cycles 10 --out-dir test6 --dump True > test6.out 2>&1

# test7: opt (hess / rfo)
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --max-cycles 5 --out-dir test7 --dump True > test7.out 2>&1

# test8: tsopt (grad / dimer)
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode grad --max-cycles 10 --out-dir test8 --dump True > test8.out 2>&1

# test9: tsopt (hess / rsirfo)
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --max-cycles 5 --out-dir test9 --dump True > test9.out 2>&1

# test10: freq
mlmm freq -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --out-dir test10 > test10.out 2>&1

# test11: irc
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --max-cycles 5 --out-dir test11 > test11.out 2>&1

# test12: dft (lightweight: hf/sto-3g)
mlmm dft -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --func-basis 'hf/sto-3g' --grid-level 0 --conv-tol 1e-5 --max-cycle 40 --out-dir test12 > test12.out 2>&1

# test13: scan (1D)
mlmm scan -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',2.0)]" --max-step-size 2.0 --max-cycles 5 --preopt False --endopt False --out-dir test13 > test13.out 2>&1

# test14: scan2d
mlmm scan2d -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',1.0,3.0),('PRE 8 C1','PRE 8 C8',1.0,3.0)]" --max-step-size 2.0 --out-dir test14 > test14.out 2>&1

# test15: scan3d
mlmm scan3d -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',1.0,3.0),('PRE 8 C1','PRE 8 C8',1.0,3.0),('PRE 8 C1','PRE 8 C7',1.0,2.0)]" --max-step-size 2.0 --out-dir test15 > test15.out 2>&1

# test16: all (minimal: no tsopt/thermo/dft)
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge 'PRE:0' -q -1 -m 1 --tsopt False --thermo False --dft False --out-dir test16 > test16.out 2>&1

# test17: all (full: tsopt + thermo + dft)
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 --ligand-charge 'PRE:0' -q -1 -m 1 --tsopt True --thermo True --dft True --tsopt-mode grad --tsopt-max-cycles 5 --dft-func-basis 'hf/sto-3g' --dft-grid-level 0 --dft-conv-tol 1e-5 --dft-max-cycle 40 --out-dir test17 > test17.out 2>&1

# test18: tsopt radius-hessian 0.0
mlmm tsopt -i p_complex.pdb --parm p_complex.parm7 --model-pdb test17/ml_region.pdb --no-detect-layer -q -1 -m 1 --opt-mode grad --max-cycles 10 --radius-hessian 0.0 --out-dir test18 > test18.out 2>&1

# test19: tsopt radius-hessian 3.6
mlmm tsopt -i p_complex.pdb --parm p_complex.parm7 --model-pdb test17/ml_region.pdb --no-detect-layer -q -1 -m 1 --opt-mode grad --max-cycles 10 --radius-hessian 3.6 --out-dir test19 > test19.out 2>&1

# test20-26: --dry-run validation
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode grad --dry-run --out-dir test20 > test20.out 2>&1
mlmm tsopt -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --dry-run --out-dir test21 > test21.out 2>&1
mlmm freq -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --dry-run --out-dir test22 > test22.out 2>&1
mlmm scan -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --scan-lists "[('PRE 8 O1\'','PRE 8 C3',2.0)]" --dry-run --out-dir test23 > test23.out 2>&1
mlmm dft -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --dry-run --out-dir test24 > test24.out 2>&1
mlmm path-search -i r_complex_layered.pdb p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --dry-run --out-dir test25 > test25.out 2>&1
mlmm irc -i p_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --dry-run --out-dir test26 > test26.out 2>&1

# test27: --microiter
mlmm opt -i r_complex_layered.pdb --parm p_complex.parm7 -q -1 -m 1 --opt-mode hess --microiter --max-cycles 2 --out-dir test27 > test27.out 2>&1

# test28: add-elem-info
cp r_complex.pdb /tmp/test_add_elem.pdb
mlmm add-elem-info -i /tmp/test_add_elem.pdb -o /tmp/test_add_elem_out.pdb > test28.out 2>&1

# test29: trj2fig (uses test6 trajectory if available)
if [ -d test6 ] && [ -f test6/opt_trj.xyz ]; then
  mlmm trj2fig -i test6/opt_trj.xyz --ref-pdb r_complex_layered.pdb -o test29.png > test29.out 2>&1
fi

# test30: energy-diagram
mlmm energy-diagram -i 0 12.5 4.3 18.7 -1.2 -o test30.png > test30.out 2>&1
