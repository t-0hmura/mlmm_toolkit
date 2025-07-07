#!/bin/sh
#PBS -N uma_test
#PBS -q default
#PBS -l nodes=1:ppn=16:gpus=1,mem=60GB,walltime=72:00:00

hostname
set -e
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR

. /home/apps/Modules/init/profile.sh
module load cuda/12.8

source /home/tohmura/miniconda3/etc/profile.d/conda.sh
conda activate mlmm
touch ~/.pysisyphusrc

python3 3-1_opt_ase.py
python3 3-2a_opt_pysis.py
mlmm 3-2b_opt_pysis.yaml
python3 3-3_test_core.py
