# `oniom-import`

Import a Gaussian/ORCA ONIOM input file and reconstruct an XYZ plus a B-factor layered PDB. The XYZ comment line carries the QM-region charge and multiplicity (`q=<charge> m=<multiplicity>`), and the layered PDB encodes ML/Movable-MM/Frozen-MM layers in its B-factor column. Import mode is taken from `--mode` or inferred from the input suffix.

## When to use

- Bring an externally prepared Gaussian/ORCA ONIOM input back into the ML/MM toolchain as an XYZ + layered PDB pair.
- Recover atom/residue metadata from a reference PDB when reconstructing a layered structure (atom count must match).

## Quick examples

```bash
mlmm oniom-import -i ts_guess.inp -o ts_guess_imported
```

```bash
# Let mode be inferred from suffix (.gjf/.com -> g16, .inp -> orca)
mlmm oniom-import -i model.gjf -o model_imported

# Force mode explicitly
mlmm oniom-import -i model.inp --mode orca -o model_imported

# Keep atom/residue metadata from a reference PDB (atom count must match)
mlmm oniom-import -i model.inp --ref-pdb complex_layered.pdb -o model_imported
```

## Inputs

Command form:

```bash
mlmm oniom-import -i INPUT.[gjf|com|inp] [--mode g16|orca] \
  [-o OUT_PREFIX] [--ref-pdb REF.pdb]
```

| Input | Required | Notes |
| --- | --- | --- |
| `-i, --input` | yes | Input ONIOM file (`.gjf`/`.com` for g16, `.inp` for ORCA). |
| `--mode` | no | Import mode (`g16`/`orca`). If omitted, inferred from input suffix. |
| `-o, --out-prefix` | no | Output prefix. Defaults to the input stem in the current directory. |
| `--ref-pdb` | no | Reference PDB to preserve atom naming/residue metadata (atom count must match). |

The full flag list (including advanced options) is in the generated [command reference](reference/commands/index.md).

## Workflow

1. Resolve import mode from `--mode` or input suffix.
2. Parse coordinates and layer-related atom sets:
 - Gaussian mode: ONIOM coordinate rows (`H`/`L` layer markers).
 - ORCA mode: `%qmmm` (`QMAtoms`/`ActiveAtoms`) + `* xyz` block.
3. Parse the QM-region charge and multiplicity (Gaussian: the ONIOM charge/multiplicity line — the six-integer line; ORCA: the `* xyz <charge> <mult>` header) and emit them as `q=<charge> m=<multiplicity>` in the XYZ comment line.
4. Write `<out_prefix>.xyz`.
5. Write `<out_prefix>_layered.pdb`:
 - Without `--ref-pdb`: generic atom/residue labels are generated.
 - With `--ref-pdb`: original naming/residue metadata is preserved while coordinates/B-factors are updated.

## Outputs

```text
<out_prefix>.xyz
<out_prefix>_layered.pdb
```

- `<out_prefix>.xyz` exists and atom count matches the source ONIOM input.
- `<out_prefix>_layered.pdb` exists and B-factor values encode layers (`0/10/20`).
- Log lines report parsed mode, atom count, and QM/movable/frozen counts.

## See Also

- [oniom_export](oniom-export.md) — Export to Gaussian ONIOM / ORCA QM/MM
- [define_layer](define-layer.md) — Layer assignment conventions
- [troubleshooting](troubleshooting.md) — Detailed troubleshooting guide
