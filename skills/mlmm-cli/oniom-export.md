# `mlmm oniom-export`

## Purpose

Export an `mlmm-toolkit` system (parm7 + layer-encoded PDB) as a
Gaussian g16 ONIOM input (`.gjf`/`.com`) or an ORCA ONIOM input
(`.inp`). Useful for manually comparing against a reference DFT/MM
calculation, re-running with a third-party engine, or feeding into a
pipeline that expects a Gaussian-style input.

The reverse direction is `oniom-import.md`.

Both modes require a CMAP-free `parm7`. The exporter fails before writing if
CMAP terms are present because neither target MM representation applies them
faithfully; runtime mlmm calculations may still use CMAP in both MM layers.

## Synopsis

```bash
mlmm oniom-export --parm enzyme.parm7 -i complex_layered.pdb \
    [--model-pdb model.pdb] \
    -o oniom.gjf \
    [--mode g16|orca] \
    [--method 'wB97X-D/def2-svp'] \
    -q 0 [-m 1] \
    [--nproc 8] [--mem 16GB]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `--parm` | path | required | Amber `parm7` topology |
| `-i, --input` | path | required | MLMM layered PDB; atom order must match parm7 and B-factors define movable/frozen atoms |
| `--model-pdb` | path | none | PDB defining QM-region atoms (B-factor 0 atoms used otherwise) |
| `-o, --output` | path | required | Output path. Suffix `.gjf` / `.com` → g16; `.inp` → ORCA (when `--mode` omitted) |
| `--mode` | choice | inferred | `g16` or `orca`; falls back to `-o` suffix |
| `--method` | str | mode-dependent | QM method + basis set, e.g. `'wB97X-D/def2-svp'` |
| `-q, --charge` | int | required | Charge of the QM region |
| `-m, --multiplicity` | int | `1` | Multiplicity of the QM region |
| `--nproc` | int | `8` | Processor count (g16 nprocshared / ORCA %pal nprocs) |
| `--mem` | str | `16GB` | Memory allocation (g16) |
| `--total-charge` / `--total-mult` | int | none | ORCA `Charge_Total` / `Mult_Total` for the full QM+MM system |
| `--orcaff` | path | derived | Path to `ORCAFF.prms` (ORCA mode) |
| `--convert-orcaff / --no-convert-orcaff` | flag | `--convert-orcaff` | Auto-run `orca_mm -convff -AMBER` if ORCAFF.prms missing |
| `--link-atom-method` | choice | `scaled` | Link-H placement: `scaled` (g-factor) or `fixed` (1.09/1.01 Å) |
| `--element-check / --no-element-check` | flag | `--element-check` | Validate element symbols before export |
| `--help-advanced` | flag | — | Reveal advanced flags |

Use `--mode` to select the engine. A single `--method` covers the QM
region; MM is parm7-driven.

## Examples

### g16 ONIOM

```bash
mlmm oniom-export --parm enzyme.parm7 -i complex_layered.pdb \
    -o complex_oniom.gjf \
    --method 'wB97X-D/def2-svp' \
    -q 0 -m 1 \
    --nproc 16 --mem 32GB
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
explicit. When the parameter file is missing and `--convert-orcaff` is
enabled, the exporter attempts `orca_mm -convff -AMBER`. If conversion is
disabled or `orca_mm` is unavailable, it still writes the `.inp`, reports the
unresolved parameter path, and leaves conversion as a required step before
running ORCA.

## Output

A target-format text input file:

- **g16 `.gjf`/`.com`**: standard Gaussian ONIOM input (with QM + MM
  partitioning by layer). PDB input adds an atom-order identity marker.
- **ORCA `.inp`**: ORCA-style multi-layer block referencing
  `<parm7_stem>.ORCAFF.prms`; that parameter file is reused or conditionally
  generated when conversion is available. PDB input adds the same identity
  marker as a comment.

## Caveats

- The QM-region selection comes from `--model-pdb` if provided, else
  from the input PDB's layer B-factors (B-factor=0 marks the ML/QM region).
- Multiplicity defaults to `1`; specify explicitly for radicals.
- A missing ORCAFF file does not prevent input export, but the resulting
  `.inp` cannot run until the reported parameter path is populated.
- Pass the original `--ref-pdb` on import to recover atom numbering, chains,
  residues, and atom names. The embedded digest verifies its atom order;
  mismatches fail even when the legacy-order override is requested.

## See also

- `oniom-import.md` — reverse direction.
- `../mlmm-structure-io/gjf.md` — Gaussian gjf format reference.
- `../mlmm-structure-io/pdb.md` § "B-factor layer encoding" — what the
  QM/MM layer assignment is read from.
- `../mlmm-structure-io/parm7.md` — parm7 contents needed by the exporter.
