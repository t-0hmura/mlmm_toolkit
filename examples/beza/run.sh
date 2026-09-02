#!/usr/bin/env bash
set -euo pipefail

# BezA GPP C6-methyltransferase: endpoint-MEP and staged-scan workflows.
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
work_dir=${1:-"$PWD/mlmm_beza_output"}
mkdir -p -- "$(dirname -- "$work_dir")"
mkdir -- "$work_dir"  # refuse to overwrite an earlier run
cd -- "$work_dir"

mlmm all -i "$script_dir/1.R.pdb" "$script_dir/3.P.pdb" \
  -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
  --refine-path --tsopt --thermo --out-dir result_mep > mep.log 2>&1

mlmm all -i "$script_dir/1.R.pdb" \
  -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
  --scan-lists '[("CS1 SAM 320","C7 GPP 321",1.50),("CS1 SAM 320","SD SAM 320",3.30)]' \
               '[("C7 GPP 321","H11 GPP 321",2.90),("OE2 GLU 186","H11 GPP 321",1.00)]' \
  --tsopt --thermo --out-dir result_scan > scan.log 2>&1
