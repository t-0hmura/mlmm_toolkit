# UMA backend (uma.md)

UMA (**U**niversal **M**odel for **A**toms, Meta FAIR) is the default
backend for `mlmm-toolkit`. It covers the broadest element / chemistry
range of the four backends and is the reference for the published
benchmark numbers.

## Install

UMA pulls in via `fairchem-core`, which is a **core dependency** of
`mlmm-toolkit`, so you don't need an extras flag:

```bash
pip install mlmm-toolkit                      # fairchem-core comes along
```

Confirm:

```bash
python -c "import fairchem; print('fairchem :', fairchem.__version__)"
python -c "import mlmm.defaults as d; print('default backend:', d.MLMM_CALC_KW['backend'])"
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

Available models (set via `mlmm.defaults.MLMM_CALC_KW["uma_model"]` or
the `calc.uma_model` YAML key):

| config string (`uma_model`) | paper notation | HuggingFace repo | Notes |
|---|---|---|---|
| `uma-s-1p1` (default) | UMA-S-1.1 / UMA-s-1.1 | `facebook/UMA-S-1.1` | Smaller / faster, sufficient for most workflows |
| `uma-m-1p1` | UMA-M-1.1 / UMA-m-1.1 | `facebook/UMA-M-1.1` | Larger, slightly more accurate, ~3× slower |

`p` is the dot replacement used by fairchem-core's config parser
(`1p1` ↔ `1.1`).

Inspect the full default kwarg dict:

```bash
python -c "import mlmm.defaults as d; print(d.MLMM_CALC_KW)"
```

## Backend-specific kwargs (`MLMM_CALC_KW`)

UMA-relevant entries in `mlmm.defaults.MLMM_CALC_KW` (set via `--config`
YAML under the `calc:` block, or via the appropriate CLI flag):

| Key | Purpose |
|---|---|
| `backend` | `'uma'` (default) / `'orb'` / `'mace'` / `'aimnet2'` |
| `uma_model` | `'uma-s-1p1'` (default) or `'uma-m-1p1'` |
| `uma_task_name` | `'omol'` (default — organic molecules + 1st-row metals) |
| `ml_device` | `'auto'` (default), `'cuda'`, or `'cpu'` |
| `ml_cuda_idx` | GPU ordinal when `ml_device='cuda'` |
| `hessian_calc_mode` | `'FiniteDifference'` (default) or `'Analytical'` |
| `H_double` | Upcast Hessian to FP64 (default `True`) |
| `out_hess_torch` | Return torch tensor Hessian (default `True`) |
| `freeze_atoms` | 1-based indices of atoms held fixed (link-atom parents, frozen residues) |
| `return_partial_hessian` | Skip frozen-atom Hessian rows for memory (default `True`) |

The MM side of `mlmm-toolkit` is configured separately (`mm_backend`,
`mm_threads`, `mm_device`, …) — see `mlmm-cli/SKILL.md` for the full
list. `mlmm-toolkit` is single-GPU on the ML side; there is no Ray /
multi-worker path.

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `e3nn` install conflict | UMA's `fairchem-core` pin clashes with `mace-torch`. Use a separate env for MACE (see `mace.md`). |
| `uma-m-1p1` runs out of VRAM during freq | Switch `hessian_calc_mode` to `'FiniteDifference'`, or use `uma-s-1p1`. |
| First call is slow (10–30 s) | One-time model download + JIT compile. The cache lives at `~/.cache/huggingface/hub/`. |
| `GatedRepoError` / `401 Unauthorized` | HuggingFace token missing or lacks access to the gated UMA repo — re-run `huggingface-cli login`. |

## See also

- `env-cuda.md` — torch + CUDA setup (UMA needs CUDA-enabled torch).
- `core.md` — `mlmm-toolkit` itself.
- `mace.md` — alternate backend, requires a **separate** env.
- `mlmm-cli/tsopt.md`, `mlmm-cli/freq.md` — `--hessian-calc-mode` choices.