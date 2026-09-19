# `mlmm dft`

## Purpose

Single-point DFT energy on an arbitrary geometry, via PySCF (CPU) or
GPU4PySCF (CUDA, x86_64). Use as a post-MLIP single-point energy evaluation on R / TS / P
geometries from `irc` / `tsopt`, or as a standalone DFT driver on any
input.

## Synopsis

```bash
mlmm dft -i geom.{pdb,xyz} --parm real.parm7 \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--func-basis 'wb97m-v/def2-svp'] \
    [--engine gpu|cpu] \
    [--embedcharge --embedcharge-cutoff ANGSTROM] \
    [--dft-nprocs INT --dft-mem SIZE] \
    [-o ./result_dft/]
```


## ML/MM-aware flags (mlmm-toolkit specific)

Beyond the common flags below, **`mlmm-toolkit` requires an Amber
topology** and supports layer-aware selection. Most subcommands accept:

| flag | purpose |
|---|---|
| `--parm FILE` | Amber `parm7` topology of the whole enzyme — **required** |
| `--model-pdb FILE` | Explicit ML-region PDB; takes precedence over B-factor ML membership |
| `--detect-layer` | Automatically read B-factor layers; explicit ML membership retains valid movable/frozen MM layers. Enabled by default. |
| `--model-indices` | Explicit ML atom indices used when `--model-pdb` is omitted; takes precedence over B-factor ML membership |
| `--ref-pdb FILE` | Full-enzyme PDB used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `-q, --charge` | **ML-region** charge (not whole-system) |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | `.pdb` / `.xyz` (XYZ requires `--ref-pdb`) |
| `-q` / `-l` / `-m` | — | — | ML-region charge / ligand-charge mapping / multiplicity (XYZ input always needs `--ref-pdb`) |
| `--ref-pdb` | path | none | Reference PDB so `-l` works on `.xyz` input |
| `--func-basis` | str | `wb97m-v/def2-svp` | `'FUNC/BASIS'` |
| `--engine` | choice {gpu,cpu} | `gpu` | `gpu` (GPU4PySCF) or `cpu` (PySCF) |
| `--lowmem/--no-lowmem` | bool | `True` | `gpu4pyscf.dft.rks_lowmem.RKS` on closed-shell GPU across calculator workflows, including electrostatic embedding; open-shell GPU and CPU use standard direct-JK RKS/UKS. `--no-lowmem` enables density fitting by default. |
| `--embedcharge` / `--embedcharge-cutoff` | toggle / Å | off / `12.0` | Embed selected MM point charges in the PySCF Hamiltonian. |
| `--dft-nprocs` / `--dft-mem` | int / size | auto / auto | Override scheduler/environment-derived thread count and memory limit. |
| `--config` | path | none | YAML config file |
| `-o, --out-dir` | path | `./result_dft/` | Output directory |
| `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

### Default DFT//MLIP/MM on a TS

```bash
mlmm dft -i seg_01/ts.pdb --parm real.parm7 \
    -l 'SAM:1,GPP:-3' \
    --func-basis 'wb97m-v/def2-tzvpd' \
    --engine gpu
```

### Lighter basis for benchmark scans

```bash
mlmm dft -i seg_01/ts.pdb --parm real.parm7 -l 'SAM:1,GPP:-3' \
    --func-basis 'wb97m-v/def2-svp' \
    --engine gpu \
    -o result_dft_svp
```

### CPU PySCF (aarch64 / no GPU)

```bash
mlmm dft -i ts.xyz --ref-pdb real.pdb --parm real.parm7 -q 0 -m 1 \
    --func-basis 'wb97m-v/def2-svp' \
    --engine cpu \
    -o result_dft_cpu
```

## Output

```
result_dft/
├── result.yaml                 # ML(dft)/MM energies + per-atom charges/spin densities
├── result.json                 # only when --out-json
├── summary.json                # mirror of result.json (only when --out-json)
├── ml_region_without_linkH.xyz # exact ML selection before generated link-H
├── ml_region_with_linkH.xyz    # ML region + generated link-H (PySCF input snapshot)
├── ml_region_without_linkH.pdb # PDB input with --convert-files; topology-bearing companion
└── ml_region_with_linkH.pdb    # PDB input with --convert-files; generated link-H as HL/LKH
```

`result.json` keys (when `--out-json`):

```python
import json
d = json.load(open("result_dft/result.json"))
print(d["energy_hartree"])     # ML-region DFT energy
print(d["xc_functional"])      # e.g. "wb97m-v"
print(d["basis_set"])          # e.g. "def2-tzvpd"
print(d["engine"])             # "gpu4pyscf(rks_lowmem)" / "gpu4pyscf" / "pyscf(cpu)"
print(d["used_lowmem"])        # True when rks_lowmem.RKS was used
print(d["converged"])
```

`result.yaml` records the effective functional/basis, grid and convergence
settings, engine/low-memory state, recombined energy, and population analyses.

## Engine choice

| `--engine` | When | Cost |
|---|---|---|
| `gpu` | x86_64 + a supported CUDA/GPU4PySCF stack | Pilot the target system |
| `cpu` | aarch64, no supported GPU stack, or explicit CPU execution | Pilot the target system |

aarch64 (`uname -m`) **requires `--engine cpu` explicitly**:
`gpu4pyscf-cuda13x` ships x86_64 wheels only, so `--engine gpu` (the
default) raises `ClickException` on aarch64 rather than silently
falling back.

## Common errors

| Symptom | Fix |
|---|---|
| `OSError: libcusolver.so.11 not found` | `mlmm-install-backends/env-cuda.md` (LD_LIBRARY_PATH order) |
| `cupy ... invalid device ordinal` | Keep scheduler-provided `CUDA_VISIBLE_DEVICES`; use a valid local ordinal (usually 0 for a one-GPU allocation). |
| `RuntimeError: CUDA out of memory` | Try the same method on CPU or a larger-memory GPU. A smaller basis/grid is a different method and must be labeled and revalidated. |
| aarch64 `--engine gpu` raises `ClickException` ("GPU backend failed...") | `gpu4pyscf-cuda13x` is x86_64 only; re-submit with `--engine cpu` |

## Caveats

- `mlmm dft` remains an energy/population-analysis **single-point** command.
  As an optional high-level calculator backend, `sp`, `opt`, `tsopt`, `irc`,
  `freq`, scan/path workflows, and `all` also accept `-b dft`; those iterative
  DFT/MM paths reuse the last converged density/orbitals.
- `--func-basis` follows PySCF naming; cross-check with
  `python -c "from pyscf import gto; print(gto.basis._BASIS_DEFAULT)"`.

## See also

- `mlmm-install-backends/dft.md` — install + aarch64 handling.
- `tsopt.md`, `irc.md` — produce the geometries used for DFT single points.
- `mlmm-workflows-output/SKILL.md` — DFT//MLIP/MM recipe.
- Defaults: `import mlmm.core.defaults as d; print(d.GEOM_KW_DEFAULT)`
