# Reproducibility and determinism

MLIP inference on a GPU is **not bit-reproducible by default**: parallel
reductions (atomic adds, scatter operations) accumulate in a
hardware-scheduling-dependent order, so two runs with identical inputs differ
at the floating-point ULP level. For `mlmm` the practical size of this drift is
**~1e-7 Å in coordinates** and below 1e-7 a.u. in energies — far below any
chemically meaningful threshold. Results are *scientifically* reproducible; they
are not *bit*-identical.

If you need bit-identical output (e.g. golden-file regression tests, exact
re-runs for an audit), use the `--deterministic` flag.

## `--deterministic`

`--deterministic` is accepted by every compute subcommand
(`opt`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`,
`path-search`, `all`, `sp`). It turns on `torch.use_deterministic_algorithms`
and an `index_reduce_` shim so that the GPU run (the ML-region high-level
inference of the ONIOM calculation) is bit-reproducible.

```bash
mlmm opt -i complex.pdb --parm enzyme.parm7 -q 0 --deterministic
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -q -1 --deterministic
```

- It is **process-global**: setting it on `all` propagates to every internal
  stage; you do not pass it per stage.
- It is **slower**: deterministic scatter/reduce kernels have lower throughput
  than the default ones. Use it only when you need exact reproducibility.
- It **fails loudly**: if the current PyTorch build cannot provide a
  deterministic kernel for an operation in your run, the command raises rather
  than silently producing non-reproducible output.
- The environment variable `MLMM_STRICT_DETERMINISTIC=1` is the equivalent
  entry point for CI or the direct Python API.

### Verified behavior by backend

| ML backend | `--deterministic` |
|---|---|
| `uma` | bit-identical energy **and** forces |
| `orb` | bit-identical energy **and** forces |
| `mace` | bit-identical energy **and** forces |
| `aimnet2` | **not supported — rejected** (see below) |

The MM low-level layer (analytical `hessian_ff` force field) is deterministic
on CPU; the non-determinism that `--deterministic` addresses is the GPU ML-region
inference.

## Precision and reproducibility

Running in `--precision fp64` *reduces* the default drift but does **not** make a
GPU run bit-identical — the reduction-order non-determinism is independent of
precision. Only `--deterministic` gives bit-exactness.

`--precision fp64` and the (internal, always-on) fp64 Hessian (`H_double`) are
independent knobs; passing `--precision fp64` additionally forces the Hessian to
fp64 so the optimizer linear algebra cannot silently run in a lower precision
than the model.

For *which* precision to choose by GPU class — `fp64` on HPC datacenter GPUs
(H100 / H200 / A100) for deterministic-grade low-noise results, `fp32` (default)
on consumer cards where `fp64` is markedly slower — see
[Device & HPC Setup → Precision by GPU class](device-hpc.md#precision-by-gpu-class).

## AIMNet2 limitations

AIMNet2 supports neither route to high reproducibility, and both combinations
are rejected with a clear error rather than run misleadingly:

- **`--precision fp64`** — AIMNet2's model inputs are cast to float32 upstream,
  so an "fp64" run would not actually be fp64.
- **`--deterministic`** — AIMNet2 computes forces through a custom CUDA kernel
  that lies outside `torch.use_deterministic_algorithms` control, so its forces
  are not bit-reproducible (energy is). PyTorch's deterministic mode neither
  detects nor controls the custom op, so the limitation is reported explicitly.

For bit-reproducible runs use `uma`, `orb`, or `mace`.
