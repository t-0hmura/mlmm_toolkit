#!/usr/bin/env bash
# Toy system — quick smoke test for the ML/MM calculator.
#
# Tests all three interfaces:
#   1. ASE Calculator (MLMMASECalculator)
#   2. PySisyphus Calculator (mlmm)
#   3. MLMMCore API (compute → energy, forces, hessian)
#   4. CLI optimization (mlmm opt)

set -euo pipefail

python3 opt_ase.py
python3 opt_pysisyphus.py
python3 test_core.py

mlmm opt -i structure.pdb --parm complex.parm7 --model-pdb ml_region.pdb -q -1 -m 1 --opt-mode grad --max-cycles 10 --out-dir dump/opt_cli

echo "Toy system tests complete."
