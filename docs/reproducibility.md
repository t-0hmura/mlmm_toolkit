# Reproducibility and determinism

MLIP inference on a GPU is not guaranteed to be bit-reproducible by default:
parallel reductions can accumulate in a hardware- and software-dependent
order. Assess numerical sensitivity on the target backend, model, hardware,
and software stack.

When exact comparison is required, use `--deterministic` to request
deterministic algorithms, then verify the produced artifacts on the complete
target stack.

## `--deterministic`

`--deterministic` is accepted by every compute subcommand
(`opt`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`,
`path-search`, `all`, `sp`). It turns on `torch.use_deterministic_algorithms`
and an `index_reduce_` shim for operations controlled by mlmm-toolkit. PyTorch
raises if a selected operation lacks a deterministic implementation.

```bash
mlmm opt -i complex.pdb --parm enzyme.parm7 -q 0 --deterministic
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -q -1 --deterministic
```

- It is **process-global**: setting it on `all` propagates to every internal
  stage; you do not pass it per stage.
- It can change performance because deterministic and default kernels may use
  different implementations.
- For PyTorch operations under its control, an unavailable deterministic
  implementation raises instead of silently using a known nondeterministic
  implementation. Backend SDK and custom operations still require separate
  verification.
- The environment variable `MLMM_STRICT_DETERMINISTIC=1` is the equivalent
  entry point for CI or the direct Python API.

### Backend support

| ML backend | `--deterministic` |
|---|---|
| `uma` | deterministic mode accepted; verify the installed model/SDK on the target system |
| `orb` | deterministic mode accepted; verify the installed model/SDK on the target system |
| `mace` | deterministic mode accepted; verify the installed model/SDK on the target system |
| `aimnet2` | **not supported — rejected** (see below) |
| `custom` (`--calc-file`) | **not supported — rejected** because the supplied calculator is outside mlmm-toolkit's control |

The default `hessian_ff` MM low-level layer runs on CPU; OpenMM can instead use
its selected device, including CUDA. Exact end-to-end comparison still requires
fixed inputs (including the topology), software versions, hardware, and backend
configuration.

## Precision and reproducibility

Running in `--precision fp64` changes numerical precision but does not by itself
guarantee bit-identical GPU execution. `--deterministic` requests deterministic
algorithms; confirm exact reproducibility for the complete target stack.

`--precision fp64` and Hessian storage precision (`H_double`, default `true`) are
independent knobs; supported fp32 configurations can set `H_double: false`, while
passing `--precision fp64` forces the Hessian to
fp64 so the optimizer linear algebra cannot silently run in a lower precision
than the model.

Precision defaults are backend-specific (UMA fp32; ORB and MACE fp64). For
selection guidance and performance tradeoffs, see
[Device & HPC Setup → Backend precision defaults](device-hpc.md#backend-precision-defaults).

## AIMNet2 limitations

AIMNet2 does not support these features:

- **`--precision fp64`** — AIMNet2's model inputs are cast to float32 upstream,
  so an "fp64" run would not actually be fp64.
- **`--deterministic`** — AIMNet2 computes forces through a custom CUDA kernel
  that lies outside `torch.use_deterministic_algorithms` control, so the flag
  cannot enforce exact force repeatability. PyTorch's deterministic mode neither
  detects nor controls the custom op, so the limitation is reported explicitly.

UMA, Orb, and MACE accept the flag; verify exact repeatability for the installed
backend/model/SDK and target stack.
