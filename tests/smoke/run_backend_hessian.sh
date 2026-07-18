#!/usr/bin/env bash
# Strict numerical analytical-vs-FD Hessian lane for explicitly named backends.
set -euo pipefail

if [[ "$#" -eq 0 ]]; then
  echo "usage: bash run_backend_hessian.sh BACKEND [BACKEND ...]" >&2
  exit 2
fi
python - <<'PY'
from importlib.metadata import version

import mlmm
import torch

installed = version("mlmm-toolkit")
if mlmm.__version__ != installed:
    raise SystemExit(
        "distribution/module version mismatch: "
        f"metadata={installed}, module={mlmm.__version__}"
    )
if not torch.cuda.is_available():
    raise SystemExit("CUDA is required for the backend Hessian lane")
PY

for backend in "$@"; do
  case "$backend" in
    uma|orb|mace|aimnet2) ;;
    *) echo "unknown backend: $backend" >&2; exit 2 ;;
  esac
  python backend_analytical_hessian.py "$backend"
done
