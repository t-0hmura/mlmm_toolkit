# `mlmm oniom-export`

## Purpose

Export an `mlmm-toolkit` system (parm7 + layer-encoded PDB / XYZ) as a
Gaussian g16 ONIOM input (`.gjf`/`.com`) or an ORCA ONIOM input
(`.inp`). Useful for hand-comparing to a reference DFT/MM
calculation, re-running with a third-party engine, or feeding into a
pipeline that expects a Gaussian-style input.

The reverse direction is `oniom-import.md`.

## Synopsis

```bash
mlmm oniom-export --parm enzyme.parm7 [-i complex.pdb] \
    [--model-pdb model.pdb] \
    -o oniom.gjf \
    [--mode g16|orca] \
    [--method 'wB97X-D/def2-svp'] \
    [-q 0] [-m 1] \
    [--near 6.0] [--nproc 8] [--mem 16GB]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `--parm` | path | required | Amber `parm7` topology |
| `-i, --input` | path | none | Coordinate file (`.pdb` or `.xyz`); atom order must match parm7 |
| `--model-pdb` | path | none | PDB defining QM-region atoms (B-factor 0 atoms used otherwise) |
| `-o, --output` | path | required | Output path. Suffix `.gjf` / `.com` → g16; `.inp` → ORCA (when `--mode` omitted) |
| `--mode` | choice | inferred | `g16` or `orca`; falls back to `-o` suffix |
| `--method` | str | mode-dependent | QM method + basis set, e.g. `'wB97X-D/def2-svp'` |
| `-q, --charge` | int | required | Charge of the QM region |
| `-m, --multiplicity` | int | `1` | Multiplicity of the QM region |
| `--near` | float | `6.0` | Distance cutoff (Å) for movable/active atoms |
| `--nproc` | int | `8` | Processor count (g16) |
| `--mem` | str | `16GB` | Memory allocation (g16) |
| `--total-charge` / `--total-mult` | int | none | ORCA `Charge_Total` / `Mult_Total` for the full QM+MM system |
| `--orcaff` | path | derived | Path to `ORCAFF.prms` (ORCA mode) |
| `--convert-orcaff / --no-convert-orcaff` | flag | `--convert-orcaff` | Auto-run `orca_mm -convff -AMBER` if ORCAFF.prms missing |
| `--help-advanced` | flag | — | Reveal advanced flags |

There is **no `--engine` flag** (use `--mode`), **no `--high-method` /
`--low-method`** (single `--method` covers the QM region; MM is
parm7-driven), and **no `--two-layer / --three-layer`**.

## Examples

### g16 ONIOM

```bash
mlmm oniom-export --parm enzyme.parm7 -i complex_layered.pdb \
    -o complex_oniom.gjf \
    --method 'wB97X-D/def2-svp' \
    -q 0 -m 1 \
    --near 6.0 --nproc 16 --mem 32GB
```

### ORCA ONIOM

```bash
mlmm oniom-export --parm enzyme.parm7 -i complex_layered.pdb \
    -o complex_oniom.inp \
    --method 'B3LYP def2-SVP' \
    -q 0 -m 1 \
    --total-charge 0 --total-mult 1
```

The mode is inferred from the `.inp` suffix; `--mode orca` makes it
explicit. The exporter calls `orca_mm -convff -AMBER` to convert the
parm7 to `ORCAFF.prms` if the prms file is missing (controlled by
`--convert-orcaff`).

## Output

A single text input file the chosen engine can run directly:

- **g16 `.gjf`/`.com`**: standard Gaussian ONIOM input (with QM + MM
  partitioning by layer).
- **ORCA `.inp`** (+ `<parm7_stem>.ORCAFF.prms`): ORCA-style multi-layer
  block referencing the converted MM force field.

## Caveats

- The QM-region selection comes from `--model-pdb` if provided, else
  from the input PDB's B-factor 0 atoms (`--detect-layer`-style).
- Multiplicity defaults to `1`; specify explicitly for radicals.
- For ORCA, `--orcaff` is required when the prms file does not yet
  exist and `--convert-orcaff` is disabled.
- Round-trip (`oniom-export` → `oniom-import`) is lossy in atom
  numbering / chain / residue rebuilding; pass `--ref-pdb` on import
  to recover names.

## See also

- `oniom-import.md` — reverse direction.
- `../mlmm-structure-io/gjf.md` — Gaussian gjf format reference.
- `../mlmm-structure-io/pdb.md` § "B-factor layer encoding" — what the
  QM/MM layer assignment is read from.
- `../mlmm-structure-io/parm7.md` — parm7 contents needed by the exporter.