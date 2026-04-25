# UMA backend (uma.md)

UMA (**U**niversal **M**odel for **A**toms, Meta FAIR) is the default
backend for `mlmm-toolkit`. It covers the broadest element / chemistry
range of the four backends and is the reference for the published
benchmark numbers.

## Install

UMA pulls in via `fairchem-core`, which is a **core dependency**, so
you don't need an extras flag:

```bash
pip install mlmm                              # fairchem-core comes along
```

Confirm:

```bash
python -c "import fairchem; print('fairchem :', fairchem.__version__)"
python -c "from mlmm.backends import create_calculator; create_calculator(backend='uma', charge=0, spin=1)"
```

## HuggingFace authentication (required)

UMA model weights are gated on HuggingFace and need an authenticated
download:

```bash
pip install huggingface_hub[cli]
huggingface-cli login        # paste a Read token from huggingface.co/settings/tokens
```

The token is cached in `~/.cache/huggingface/`. Once it's there, future
runs (and PBS jobs) pick it up automatically.

If you hit `huggingface_hub.errors.GatedRepoError` or
`401 Client Error: Unauthorized`, re-run `huggingface-cli login` and
make sure the token has access to the UMA model repos
(`facebook/UMA-S-1.1`, `facebook/UMA-M-1.1`).

## CLI usage

`uma` is the default — `mlmm all -i ...` uses UMA-s-1.1 unless
overridden:

```bash
mlmm all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -b uma                       # explicit, identical to default
```

Available models (set via `--model` or via `mlmm.defaults.MLMM_CALC_KW`).
Two equivalent notations are common:

| config string (`--model`) | paper notation | HuggingFace repo | Notes |
|---|---|---|---|
| `uma-s-1p1` (default) | UMA-S-1.1 / UMA-s-1.1 | `facebook/UMA-S-1.1` | Smaller / faster, sufficient for most workflows |
| `uma-m-1p1` | UMA-M-1.1 / UMA-m-1.1 | `facebook/UMA-M-1.1` | Larger, slightly more accurate, ~3× slower |

`p` is the dot replacement used by fairchem-core's config parser
(`1p1` ↔ `1.1`). Pass the config string (`uma-s-1p1`) on the CLI.

Inspect the full default kwarg dict:

```bash
python -c "import mlmm.defaults as d; print(d.MLMM_CALC_KW)"
```

## Backend-specific flags

UMA accepts these calculator kwargs (canonical list in
`backends/__init__.py:_BACKEND_ACCEPTED_KEYS['uma']`):

| Key | Purpose |
|---|---|
| `charge`, `spin` | Total charge and spin multiplicity |
| `device` | `'cuda'`, `'cpu'`, or `'auto'` |
| `model` | `'uma-s-1.1'` or `'uma-m-1.1'` |
| `task_name` | `'omol'` (default — organic molecules + 1st-row metals) |
| `freeze_atoms` | Indices of atoms held fixed (link-atom parents, frozen residues) |
| `hessian_calc_mode` | `'FiniteDifference'` (default) or `'Analytical'` |
| `return_partial_hessian`, `hessian_double` | Memory / numerical-precision toggles |
| `workers`, `workers_per_node` | Multi-GPU inference (uses Ray) |
| `max_neigh`, `radius` | Neighbor-list cutoffs |

These are passed through `--config` YAML or the appropriate flag of the
subcommand; see `mlmm-cli/SKILL.md`.

## Multi-GPU inference (advanced)

Under heavy MEP search load you can shard inference across multiple GPUs.
Configure via a YAML config (`--config`) — there is no `--calc-kwargs`
CLI flag:

```yaml
# multi_worker.yaml
calc:
  workers: 4
  workers_per_node: 4
```

```bash
mlmm all -i ... -b uma --config multi_worker.yaml
```

This spawns a Ray worker pool. Limitations:

- **Single node only** (no cross-node sharding via `mlmm-toolkit`).
- All workers must see the same GPUs (e.g. `CUDA_VISIBLE_DEVICES=0,1,2,3`).
- Adds overhead for small graphs — disable for systems below ~100 atoms.
- **Silent Hessian downgrade**: when `workers > 1`, analytical Hessian
  is silently disabled in favor of finite-difference. See
  `mlmm/docs/uma-pysis.md` for the full caveat.

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `e3nn` install conflict | UMA's `fairchem-core` pin clashes with `mace-torch`. Use a separate env for MACE (see `mace.md`). |
| `uma-m-1.1` runs out of VRAM during freq | Switch `hessian_calc_mode` to `'FiniteDifference'`, or use `uma-s-1.1`. |
| First call is slow (10–30 s) | One-time model download + JIT compile. The cache lives at `~/.cache/huggingface/hub/`. |
| Multi-worker run crashes with `Ray actor died` | Mismatched CUDA versions across processes; fall back to single-worker. |

## See also

- `env-cuda.md` — torch + CUDA setup (UMA needs CUDA-enabled torch).
- `core.md` — `mlmm-toolkit` itself.
- `mace.md` — alternate backend, requires a **separate** env.
- `mlmm-cli/tsopt.md`, `mlmm-cli/freq.md` — `--hessian-calc-mode` choices.