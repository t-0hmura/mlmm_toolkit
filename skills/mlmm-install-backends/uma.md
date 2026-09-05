# UMA backend (uma.md)

UMA (**U**niversal **M**odel for **A**toms, Meta FAIR) is the default
backend for `mlmm-toolkit`. It covers the broadest element / chemistry
range of the four bundled backends.

## Install

UMA pulls in via `fairchem-core`, which is a **core dependency** of
`mlmm-toolkit`, so you don't need an extras flag:

```bash
pip install mlmm-toolkit                      # fairchem-core comes along
```

Confirm:

```bash
python -c "import fairchem; print('fairchem :', fairchem.__version__)"
python -c "import mlmm.core.defaults as d; print('default backend:', d.MLMM_CALC_KW['backend'])"
```

## HuggingFace authentication (required)

UMA model weights are gated on HuggingFace and need an authenticated
download:

```bash
hf auth login               # paste a Read token from huggingface.co/settings/tokens
```

The token is cached in `~/.cache/huggingface/`. Once it's there, future
runs (and PBS jobs) pick it up automatically.

If you hit `huggingface_hub.errors.GatedRepoError` or
`401 Client Error: Unauthorized`, re-run `hf auth login` and
make sure the token has access to the gated `facebook/UMA` repository
(model variants are selected by `uma_model` config string, not by
separate repos).

## CLI usage

`uma` is the default — `mlmm all -i ...` uses UMA-s-1.2 unless
overridden:

```bash
mlmm all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -b uma                       # explicit, identical to default

# pick a non-default model variant from the CLI (no YAML needed):
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -b uma --backend-model uma-m-1p1
```

Available models (select with the `--backend-model` CLI flag, or via
`mlmm.core.defaults.MLMM_CALC_KW["uma_model"]` / the `calc.uma_model` YAML key):

| config string (`uma_model`) | paper notation | HuggingFace repo | Notes |
|---|---|---|---|
| `uma-s-1p2` (default) | UMA-s-1.2 | `facebook/UMA` | Current small-model default |
| `uma-s-1p1` | UMA-s-1.1 | `facebook/UMA` | Previous small model; explicit opt-in |
| `uma-m-1p1` | UMA-m-1.1 | `facebook/UMA` | Larger model; benchmark accuracy and cost on the target system |

`p` is the dot replacement used by fairchem-core's config parser
(`1p1` ↔ `1.1`).

Inspect the full default kwarg dict:

```bash
python -c "import mlmm.core.defaults as d; print(d.MLMM_CALC_KW)"
```

## Backend-specific kwargs (`MLMM_CALC_KW`)

UMA-relevant entries in `mlmm.core.defaults.MLMM_CALC_KW` (set via `--config`
YAML under the `calc:` block, or via the appropriate CLI flag):

| Key | Purpose |
|---|---|
| `backend` | `'uma'` (default) / `'orb'` / `'mace'` / `'aimnet2'` |
| `uma_model` | `'uma-s-1p2'` (default), `'uma-s-1p1'`, or `'uma-m-1p1'` |
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
list. UMA supports `--workers > 1` with `fairchem-core[extras]` and finite-difference Hessians; analytical Hessians require one worker.

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `e3nn` install conflict | UMA's `fairchem-core` pin clashes with `mace-torch`. Use a separate env for MACE (see `mace.md`). |
| Frequency calculation runs out of VRAM | Compare Hessian modes and compatible model sizes on a representative pilot, or move Hessian assembly to CPU. |
| First call is slower than later calls | One-time model download + JIT compile. The cache lives at `~/.cache/huggingface/hub/`. |
| `GatedRepoError` / `401 Unauthorized` | HuggingFace token missing or lacks access to the gated UMA repo — re-run `hf auth login`. |

## See also

- `env-cuda.md` — torch + CUDA setup (UMA needs CUDA-enabled torch).
- `core.md` — `mlmm-toolkit` itself.
- `mace.md` — alternate backend, requires a **separate** env.
- `mlmm-cli/tsopt.md`, `mlmm-cli/freq.md` — `--hessian-calc-mode` choices.
