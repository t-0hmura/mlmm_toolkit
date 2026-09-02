#!/usr/bin/env bash
set -euo pipefail

# COMT catechol O-methyltransferase: full-system endpoint workflow.
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
work_dir=${1:-"$PWD/mlmm_comt_output"}
mkdir -p -- "$(dirname -- "$work_dir")"
mkdir -- "$work_dir"  # refuse to overwrite an earlier run
cd -- "$work_dir"

mlmm all -i "$script_dir/1.R.pdb" "$script_dir/3.P.pdb" \
  -c 'CAT,SAM,MG' -l 'CAT:-1,SAM:0,MG:2' -r 4.0 \
  --tsopt --thermo --out-dir result > run.log 2>&1
