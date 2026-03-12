# ONIOM Import (`oniom-import`)

## Overview

> **Summary:** Import a Gaussian/ORCA ONIOM input file and reconstruct an XYZ plus a B-factor layered PDB.

## Minimal example

```bash
mlmm oniom-import -i ts_guess.inp -o ts_guess_imported
```

## Output checklist
- `<out_prefix>.xyz` exists and atom count matches the source ONIOM input.
- `<out_prefix>_layered.pdb` exists and B-factor values encode layers (`0/10/20`).
- Log lines report parsed mode, atom count, and QM/movable/frozen counts.

## Common examples

```bash
# Let mode be inferred from suffix (.gjf/.com -> g16, .inp -> orca)
mlmm oniom-import -i model.gjf -o model_imported

# Force mode explicitly
mlmm oniom-import -i model.inp --mode orca -o model_imported

# Keep atom/residue metadata from a reference PDB (atom count must match)
mlmm oniom-import -i model.inp --ref-pdb complex_layered.pdb -o model_imported
```

## Usage

```bash
mlmm oniom-import -i INPUT.[gjf|com|inp] [--mode g16|orca] \
  [-o OUT_PREFIX] [--ref-pdb REF.pdb]
```

## Workflow
1. Resolve import mode from `--mode` or input suffix.
2. Parse coordinates and layer-related atom sets:
 - Gaussian mode: ONIOM coordinate rows (`H`/`L` layer markers).
 - ORCA mode: `%qmmm` (`QMAtoms`/`ActiveAtoms`) + `* xyz` block.
3. Write `<out_prefix>.xyz`.
4. Write `<out_prefix>_layered.pdb`:
 - Without `--ref-pdb`: generic atom/residue labels are generated.
 - With `--ref-pdb`: original naming/residue metadata is preserved while coordinates/B-factors are updated.

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input FILE` | Input ONIOM file (`.gjf`/`.com` for g16, `.inp` for ORCA). | Required |
| `--mode [g16|orca]` | Import mode. If omitted, inferred from input suffix. | Inferred |
| `-o, --out-prefix PATH` | Output prefix. | Input stem in current directory |
| `--ref-pdb FILE` | Reference PDB to preserve atom naming/residue metadata (atom count must match). | _None_ |

## Outputs

```text
<out_prefix>.xyz
<out_prefix>_layered.pdb
```

## See Also

- [oniom_export](oniom_export.md) -- Export to Gaussian ONIOM / ORCA QM/MM
- [define_layer](define_layer.md) -- Layer assignment conventions
- [troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
