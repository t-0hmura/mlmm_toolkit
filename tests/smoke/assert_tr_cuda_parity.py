#!/usr/bin/env python3
"""Require CPU/CUDA parity for the constrained two-anchor TR projector."""

import torch

from pysisyphus.tr_projection import active_tr_basis


coords = torch.tensor(
    [
        [0.0, 0.0, 0.0],
        [1.2, 0.0, 0.0],
        [0.1, 1.1, 0.0],
        [0.2, 0.3, 1.3],
        [1.0, 0.8, 0.7],
    ],
    dtype=torch.float32,
)
masses = torch.tensor([12.0, 1.0, 16.0, 14.0, 32.0], dtype=torch.float32)
active = [2, 3, 4]

if not torch.cuda.is_available():
    raise SystemExit("CUDA is required for the TR projector parity assertion")

cpu_basis, cpu_info = active_tr_basis(coords, masses, active)
gpu_basis, gpu_info = active_tr_basis(coords.cuda(), masses.cuda(), active)
if cpu_info.effective_rank != 1 or gpu_info.effective_rank != 1:
    raise SystemExit(
        "unexpected constrained TR rank: "
        f"cpu={cpu_info.effective_rank}, gpu={gpu_info.effective_rank}"
    )
torch.testing.assert_close(
    cpu_basis @ cpu_basis.T,
    (gpu_basis @ gpu_basis.T).cpu(),
    atol=1.0e-12,
    rtol=1.0e-12,
)
print("[smoke] PASS: constrained TR projector CPU/CUDA parity")
