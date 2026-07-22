#!/usr/bin/env bash
set -euo pipefail

# Methyltransferase — automated end-to-end workflow via `mlmm all`.
# Single PDB + scan mode: extract → scan → MEP → TS → IRC → freq

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
work_dir=${1:-"$PWD/mlmm_methyltransferase_all"}
mkdir -p -- "$(dirname -- "$work_dir")"
mkdir -- "$work_dir"  # refuse to overwrite an earlier run
cd -- "$work_dir"

mlmm all -i "$script_dir/complex.pdb" -c "SAM,PHN" -l "SAM:1,PHN:-1" \
  -s "[('SAM 359 CS1','PHN 360 C8',1.3)]" \
  --tsopt --thermo -o result_all > all.log 2>&1
