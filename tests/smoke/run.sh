#!/usr/bin/env bash
# ============================================================
#  mlmm_toolkit GPU smoke tests
#
#  Requires:  conda activate mlmm   (GPU + CUDA environment)
#  Usage:     bash run.sh           (local)
#             qsub run.sh           (PBS)
#
#  All test outputs go into test*/ and test*.out inside this
#  directory.  A successful run prints the final banner and
#  exits 0.
# ============================================================
#PBS -N mlmm_smoke
#PBS -q default
#PBS -l nodes=1:ppn=32:gpus=1,mem=120GB,walltime=72:00:00
#PBS -o /dev/null
#PBS -e /dev/null

set -euo pipefail

# Run from this script directory so relative fixture paths resolve.
cd "$(dirname "$0")"

hostname
test "${PBS_O_WORKDIR:-}" && cd "${PBS_O_WORKDIR}" || true

. /home/apps/Modules/init/profile.sh
module load cuda/12.9

source /home/tohmura/miniconda3/etc/profile.d/conda.sh
conda activate mlmm

# Clean up previous runs
rm -rf test[0-9]* test_dryrun_* test_config_* test_microiter test_trj2fig* test_add_elem* test_energy_diagram*

# ---------- fixtures ----------
Q=-1
M=1
R_FULL="r_complex.pdb"
P_FULL="p_complex.pdb"
R_PDB="r_complex_layered.pdb"
P_PDB="p_complex_layered.pdb"
PARM7="p_complex.parm7"

echo "================================================================"
echo " mlmm_toolkit smoke tests"
echo " $(date)"
echo "================================================================"

# ------------------------------------------------------------------
# 1. CLI availability: --version, --help, --help-advanced
# ------------------------------------------------------------------
mlmm --version > test1.out 2>&1
for cmd in all extract define-layer path-search path-opt opt tsopt freq \
           irc dft scan scan2d scan3d mm-parm trj2fig add-elem-info \
           oniom-export oniom-import fix-altloc energy-diagram; do
  mlmm ${cmd} --help >> test1.out 2>&1
done
for cmd in opt tsopt freq scan dft extract path-search path-opt irc all; do
  mlmm ${cmd} --help-advanced >> test1.out 2>&1
done

# ------------------------------------------------------------------
# 2. extract
# ------------------------------------------------------------------
mlmm extract \
  -i "${R_FULL}" \
  -c PRE \
  -r 5.0 \
  --exclude-backbone False \
  --ligand-charge "PRE:0" \
  -o pocket_r_smoke.pdb \
  > test2.out 2>&1

# ------------------------------------------------------------------
# 3. define-layer
# ------------------------------------------------------------------
mlmm define-layer \
  -i "${R_FULL}" \
  --model-pdb pocket_r_smoke.pdb \
  --radius-partial-hessian 3.6 \
  --radius-freeze 8.0 \
  -o r_complex_layered_smoke.pdb \
  > test3.out 2>&1

# ------------------------------------------------------------------
# 4. path-opt (gsm)
# ------------------------------------------------------------------
rm -rf test4
mlmm path-opt \
  -i "${R_PDB}" "${P_PDB}" \
  --parm "${PARM7}" \
  -q ${Q} -m ${M} \
  --max-nodes 5 \
  --max-cycles 10 \
  --preopt False \
  --climb False \
  --out-dir test4 \
  > test4.out 2>&1

# ------------------------------------------------------------------
# 5. path-opt (dmf)
# ------------------------------------------------------------------
rm -rf test5
mlmm path-opt \
  -i "${R_PDB}" "${P_PDB}" \
  --parm "${PARM7}" \
  -q ${Q} -m ${M} \
  --mep-mode dmf \
  --max-cycles 5 \
  --preopt False \
  --out-dir test5 \
  > test5.out 2>&1

# ------------------------------------------------------------------
# 6. path-search
# ------------------------------------------------------------------
rm -rf test6
mlmm path-search \
  -i "${R_PDB}" "${P_PDB}" \
  --parm "${PARM7}" \
  -q ${Q} -m ${M} \
  --max-cycles 10 \
  --out-dir test6 \
  > test6.out 2>&1

# ------------------------------------------------------------------
# 7. opt (grad / lbfgs)
# ------------------------------------------------------------------
rm -rf test7
mlmm opt \
  -i "${R_PDB}" \
  --parm "${PARM7}" \
  -q ${Q} -m ${M} \
  --opt-mode grad \
  --max-cycles 10 \
  --out-dir test7 \
  --dump True \
  > test7.out 2>&1

# ------------------------------------------------------------------
# 8. opt (hess / rfo)
# ------------------------------------------------------------------
rm -rf test8
mlmm opt \
  -i "${R_PDB}" \
  --parm "${PARM7}" \
  -q ${Q} -m ${M} \
  --opt-mode hess \
  --max-cycles 5 \
  --out-dir test8 \
  --dump True \
  > test8.out 2>&1

# ------------------------------------------------------------------
# 9. tsopt (grad / dimer)
# ------------------------------------------------------------------
rm -rf test9
mlmm tsopt \
  -i "${P_PDB}" \
  --parm "${PARM7}" \
  -q ${Q} -m ${M} \
  --opt-mode grad \
  --max-cycles 10 \
  --out-dir test9 \
  --dump True \
  > test9.out 2>&1

# ------------------------------------------------------------------
# 10. tsopt (hess / rsirfo)
# ------------------------------------------------------------------
rm -rf test10
mlmm tsopt \
  -i "${P_PDB}" \
  --parm "${PARM7}" \
  -q ${Q} -m ${M} \
  --opt-mode hess \
  --max-cycles 5 \
  --out-dir test10 \
  --dump True \
  > test10.out 2>&1

# ------------------------------------------------------------------
# 11. freq
# ------------------------------------------------------------------
rm -rf test11
mlmm freq \
  -i "${R_PDB}" \
  --parm "${PARM7}" \
  -q ${Q} -m ${M} \
  --out-dir test11 \
  > test11.out 2>&1

# ------------------------------------------------------------------
# 12. irc
# ------------------------------------------------------------------
rm -rf test12
mlmm irc \
  -i "${P_PDB}" \
  --parm "${PARM7}" \
  -q ${Q} -m ${M} \
  --max-cycles 5 \
  --out-dir test12 \
  > test12.out 2>&1

# ------------------------------------------------------------------
# 13. dft (lightweight: hf/sto-3g)
# ------------------------------------------------------------------
rm -rf test13
mlmm dft \
  -i "${R_PDB}" \
  --parm "${PARM7}" \
  -q ${Q} -m ${M} \
  --func-basis "hf/sto-3g" \
  --grid-level 0 \
  --conv-tol 1e-5 \
  --max-cycle 40 \
  --out-dir test13 \
  > test13.out 2>&1

# ------------------------------------------------------------------
# 15. scan (1D)
# ------------------------------------------------------------------
rm -rf test15
mlmm scan \
  -i "${R_PDB}" \
  --parm "${PARM7}" \
  -q ${Q} -m ${M} \
  --scan-lists "[('PRE 8 O1\\'','PRE 8 C3',2.0)]" \
  --max-step-size 2.0 \
  --max-cycles 5 \
  --preopt False \
  --endopt False \
  --out-dir test15 \
  > test15.out 2>&1

# ------------------------------------------------------------------
# 16. scan2d
# ------------------------------------------------------------------
rm -rf test16
mlmm scan2d \
  -i "${R_PDB}" \
  --parm "${PARM7}" \
  -q ${Q} -m ${M} \
  --scan-lists "[('PRE 8 O1\\'','PRE 8 C3',1.0,3.0),('PRE 8 C1','PRE 8 C8',1.0,3.0)]" \
  --max-step-size 2.0 \
  --out-dir test16 \
  > test16.out 2>&1

# ------------------------------------------------------------------
# 17. scan3d
# ------------------------------------------------------------------
rm -rf test17
mlmm scan3d \
  -i "${R_PDB}" \
  --parm "${PARM7}" \
  -q ${Q} -m ${M} \
  --scan-lists "[('PRE 8 O1\\'','PRE 8 C3',1.0,3.0),('PRE 8 C1','PRE 8 C8',1.0,3.0),('PRE 8 C1','PRE 8 C7',1.0,2.0)]" \
  --max-step-size 2.0 \
  --out-dir test17 \
  > test17.out 2>&1

# ------------------------------------------------------------------
# 18. all (minimal: no tsopt/thermo/dft)
# ------------------------------------------------------------------
rm -rf test18
mlmm all \
  -i "${R_FULL}" "${P_FULL}" \
  -c PRE \
  -r 6.0 \
  --ligand-charge "PRE:0" \
  -q ${Q} -m ${M} \
  --tsopt False \
  --thermo False \
  --dft False \
  --out-dir test18 \
  > test18.out 2>&1

# ------------------------------------------------------------------
# 19. all (full: tsopt + thermo + dft)
# ------------------------------------------------------------------
rm -rf test19
mlmm all \
  -i "${R_FULL}" "${P_FULL}" \
  -c PRE \
  -r 6.0 \
  --ligand-charge "PRE:0" \
  -q ${Q} -m ${M} \
  --tsopt True \
  --thermo True \
  --dft True \
  --tsopt-mode grad \
  --tsopt-max-cycles 5 \
  --dft-func-basis "hf/sto-3g" \
  --dft-grid-level 0 \
  --dft-conv-tol 1e-5 \
  --dft-max-cycle 40 \
  --out-dir test19 \
  > test19.out 2>&1

# ------------------------------------------------------------------
# 20-21. tsopt radius-hessian behavior
# ------------------------------------------------------------------
if [ ! -f test19/ml_region.pdb ]; then
  echo "[error] test19/ml_region.pdb not found (required for radius-hessian tests)." > test20.out
  exit 2
fi

rm -rf test20
mlmm tsopt \
  -i p_complex.pdb \
  --parm p_complex.parm7 \
  --model-pdb test19/ml_region.pdb \
  --no-detect-layer \
  -q ${Q} -m ${M} \
  --opt-mode grad \
  --max-cycles 10 \
  --radius-hessian 0.0 \
  --out-dir test20 \
  > test20.out 2>&1

rm -rf test21
mlmm tsopt \
  -i p_complex.pdb \
  --parm p_complex.parm7 \
  --model-pdb test19/ml_region.pdb \
  --no-detect-layer \
  -q ${Q} -m ${M} \
  --opt-mode grad \
  --max-cycles 10 \
  --radius-hessian 3.6 \
  --out-dir test21 \
  > test21.out 2>&1

# ------------------------------------------------------------------
# --dry-run validation (fast, no GPU computation)
# ------------------------------------------------------------------
mlmm opt \
  -i "${R_PDB}" --parm "${PARM7}" -q ${Q} -m ${M} \
  --opt-mode grad --dry-run --out-dir test_dryrun_opt \
  > test_dryrun_opt.out 2>&1

mlmm tsopt \
  -i "${P_PDB}" --parm "${PARM7}" -q ${Q} -m ${M} \
  --opt-mode hess --dry-run --out-dir test_dryrun_tsopt \
  > test_dryrun_tsopt.out 2>&1

mlmm freq \
  -i "${R_PDB}" --parm "${PARM7}" -q ${Q} -m ${M} \
  --dry-run --out-dir test_dryrun_freq \
  > test_dryrun_freq.out 2>&1

mlmm scan \
  -i "${R_PDB}" --parm "${PARM7}" -q ${Q} -m ${M} \
  --scan-lists "[('PRE 8 O1\\'','PRE 8 C3',2.0)]" \
  --dry-run --out-dir test_dryrun_scan \
  > test_dryrun_scan.out 2>&1

mlmm dft \
  -i "${R_PDB}" --parm "${PARM7}" -q ${Q} -m ${M} \
  --dry-run --out-dir test_dryrun_dft \
  > test_dryrun_dft.out 2>&1

mlmm path-search \
  -i "${R_PDB}" "${P_PDB}" --parm "${PARM7}" -q ${Q} -m ${M} \
  --dry-run --out-dir test_dryrun_ps \
  > test_dryrun_ps.out 2>&1

mlmm irc \
  -i "${P_PDB}" --parm "${PARM7}" -q ${Q} -m ${M} \
  --dry-run --out-dir test_dryrun_irc \
  > test_dryrun_irc.out 2>&1

# ------------------------------------------------------------------
# --config YAML test
# ------------------------------------------------------------------
cat > /tmp/mlmm_smoke_config.yaml << 'YAMLEOF'
opt:
  thresh: gau_loose
  max_cycles: 5
YAMLEOF

rm -rf test_config_opt
mlmm opt \
  -i "${R_PDB}" --parm "${PARM7}" -q ${Q} -m ${M} \
  --opt-mode grad \
  --config /tmp/mlmm_smoke_config.yaml \
  --dry-run --out-dir test_config_opt \
  > test_config_opt.out 2>&1

# ------------------------------------------------------------------
# --microiter test
# ------------------------------------------------------------------
rm -rf test_microiter
mlmm opt \
  -i "${R_PDB}" --parm "${PARM7}" -q ${Q} -m ${M} \
  --opt-mode hess --microiter --max-cycles 2 \
  --out-dir test_microiter \
  > test_microiter.out 2>&1

# ------------------------------------------------------------------
# add-elem-info
# ------------------------------------------------------------------
cp "${R_FULL}" /tmp/test_add_elem.pdb
mlmm add-elem-info \
  -i /tmp/test_add_elem.pdb \
  -o /tmp/test_add_elem_out.pdb \
  > test_add_elem.out 2>&1

# ------------------------------------------------------------------
# trj2fig (uses test7 trajectory if available)
# ------------------------------------------------------------------
if [ -d test7 ] && [ -f test7/opt_trj.xyz ]; then
  mlmm trj2fig \
    -i test7/opt_trj.xyz \
    --ref-pdb "${R_PDB}" \
    -o test_trj2fig.png \
    > test_trj2fig.out 2>&1
fi

# ------------------------------------------------------------------
# energy-diagram
# ------------------------------------------------------------------
mlmm energy-diagram \
  -i 0 12.5 4.3 18.7 -1.2 \
  -o test_energy_diagram.png \
  > test_energy_diagram.out 2>&1

echo "================================================================"
echo " mlmm_toolkit smoke tests completed successfully"
echo " $(date)"
echo "================================================================"
