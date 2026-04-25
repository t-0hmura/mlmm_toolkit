# `mlmm dft`

## Purpose

Single-point DFT energy on an arbitrary geometry, via PySCF (CPU) or
GPU4PySCF (CUDA, x86_64). Use as a post-MLIP refinement on R / TS / P
geometries from `irc` / `tsopt`, or as a standalone DFT driver on any
input.

## Synopsis

```bash
mlmm dft -i geom.{pdb,xyz,gjf} \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--func-basis 'wb97m-v/def2-svp'] \
    [--engine gpu|cpu] \
    [-o ./result_dft/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | `.pdb` / `.xyz` / `.gjf` |
| `-q` / `-l` / `-m` | — | — | Charge / spin (required for `.xyz` without `--ref-pdb`) |
| `--ref-pdb` | path | none | Reference PDB so `-l` works on `.xyz` input |
| `--func-basis` | str | `wb97m-v/def2-tzvpd` | `'FUNC/BASIS'` |
| `--engine` | str | `gpu` | `gpu` (GPU4PySCF) or `cpu` (PySCF) |
| `--config` | path | none | YAML config file |
| `-o, --out-dir` | path | `./result_dft/` | Output directory |
| `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

### Default DFT//MLIP on a TS

```bash
mlmm dft -i seg_01/ts.pdb \
    -l 'SAM:1,GPP:-3' \
    --func-basis 'wb97m-v/def2-tzvpd' \
    --engine gpu
```

### Lighter basis for benchmark scans

```bash
mlmm dft -i seg_01/ts.pdb -l 'SAM:1,GPP:-3' \
    --func-basis 'wb97m-v/def2-svp' \
    --engine gpu \
    -o result_dft_svp
```

### CPU PySCF (aarch64 / no GPU)

```bash
mlmm dft -i ts.xyz -q 0 -m 1 \
    --func-basis 'wb97m-v/def2-svp' \
    --engine cpu \
    -o result_dft_cpu
```

## Output

```
result_dft/
├── result.json
├── result.yaml             # PySCF-style detail dump
├── input_geometry.xyz       # echoed input
└── dft.log
```

`result.json` keys:

```python
import json
d = json.load(open("result_dft/result.json"))
print(d["energy_hartree"])
print(d["method"])           # e.g. "wb97m-v/def2-tzvpd"
print(d["engine"])           # "gpu" / "cpu"
print(d["status"])
```

`result.yaml` carries the full PySCF / GPU4PySCF runtime info: basis
expansion, grid_level, SCF iterations, Mulliken / Loewdin / IAO charges,
spin densities. Useful for debugging convergence problems.

## Engine choice

| `--engine` | When | Cost |
|---|---|---|
| `gpu` | x86_64 + CUDA + > 100 atoms | ~1–10 h on RTX-class GPU per single point |
| `cpu` | aarch64, no GPU, or < 100 atoms | ~10–100× slower than GPU |

aarch64 (`uname -m`) **forces CPU** because `gpu4pyscf-cuda12x` ships
x86_64 wheels only.

## Common errors

| Symptom | Fix |
|---|---|
| `OSError: libcusolver.so.11 not found` | `mlmm-install-backends/env-cuda.md` (LD_LIBRARY_PATH order) |
| `cupy ... invalid device ordinal` | `unset CUDA_VISIBLE_DEVICES` |
| `RuntimeError: CUDA out of memory` | Lower `grid_level`, switch to `def2-svp`, or `--engine cpu` |
| aarch64 `--engine gpu` falls back silently | Expected; check log says `engine: cpu` |

## Caveats

- `mlmm dft` runs only **single points**, not optimization.
  For DFT-level geometry refinement, run `mlmm tsopt -b dft`
  (DFT calculator selectable via `-b`) or use a separate code.
- `--func-basis` follows PySCF naming; cross-check with
  `python -c "from pyscf import gto; print(gto.basis._BASIS_DEFAULT)"`.
- xTB-ALPB (`--solvent`) **does not stack** with PySCF's PCM; pick one.

## See also

- `mlmm-install-backends/dft.md` — install + aarch64 handling.
- `tsopt.md`, `irc.md` — produce the geometries you DFT-refine.
- `mlmm-workflows-output/SKILL.md` — DFT//MLIP recipe.
- Defaults: `import mlmm.defaults as d; print(d.GEOM_KW_DEFAULT)`
