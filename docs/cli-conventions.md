# CLI Conventions

This page documents the conventions used across all `mlmm` commands. Understanding these conventions helps you write correct commands and avoid common errors.

---

## Boolean Options

All boolean CLI options require an explicit value -- you cannot use flag-style `--tsopt` alone:

```bash
# Correct (any of these work)
--tsopt True --thermo True --dft False
--tsopt true --thermo TRUE --dft false   # case-insensitive
--tsopt 1 --thermo yes --dft 0           # 1/0, yes/no also accepted

# Wrong (will not work)
--tsopt         # flag-style (no value) is not supported
```

The CLI accepts `True`, `true`, `TRUE`, `1`, `yes`, `Yes`, `y`, `t` for truthy values, and `False`, `false`, `FALSE`, `0`, `no`, `No`, `n`, `f` for falsy values.

Common boolean options:
- `--tsopt`, `--thermo`, `--dft` -- enable post-processing stages
- `--freeze-links` -- freeze link-hydrogen parents (default: `True`)
- `--dump` -- write trajectory files
- `--preopt`, `--endopt` -- pre/post optimization toggles
- `--climb` -- enable climbing image in MEP search
- `--convert-files` -- generate PDB/GJF companion files

---

## ML/MM Required Options

Most subcommands (except `all`, `extract`, `mm-parm`, and `define-layer`) require two additional options for ML/MM calculations:

```bash
--real-parm7 real.parm7    # Amber parm7 topology file for the full (real) system
--model-pdb model.pdb      # PDB file defining the ML (model) region atoms
```

The `all` command generates these files automatically during the workflow (via `mm-parm` and `define-layer`). When running subcommands individually, you must provide them explicitly.

```bash
# Example: running path-search individually
mlmm path-search -i R.pdb P.pdb --real-parm7 real.parm7 --model-pdb model.pdb -q 0 -m 1
```

---

## B-factor Layer Encoding

mlmm_toolkit uses the PDB B-factor column to encode the 3-layer ML/MM partitioning:

| Layer | B-factor | Meaning |
|-------|----------|---------|
| ML | 0.0 | UMA energy/force/Hessian |
| Movable-MM | 10.0 | MM atoms allowed to move |
| Frozen | 20.0 | Coordinates fixed |

The `define-layer` subcommand writes these B-factors into the PDB. You can inspect layer assignments in any molecular viewer that supports B-factor coloring.

A tolerance of 1.0 is used when reading B-factors, so values near 0/10/20 are mapped to ML/Movable/Frozen.

---

## Residue Selectors

Residue selectors identify which residues to use as substrates or extraction centers.

### By residue name
```bash
-c 'SAM,GPP'          # Select all residues named SAM or GPP
-c 'LIG'              # Select all residues named LIG
```

### By residue ID
```bash
-c '123,456'          # Residues 123 and 456
-c 'A:123,B:456'      # Chain A residue 123, Chain B residue 456
-c '123A'             # Residue 123 with insertion code A
-c 'A:123A'           # Chain A, residue 123, insertion code A
```

### By PDB file
```bash
-c substrate.pdb      # Use coordinates from a separate PDB to locate substrates
```

```{note}
When selecting by residue name, if multiple residues share the same name, **all** matches are included and a warning is logged.
```

---

## Charge Specification

### Per-residue mapping (recommended)
```bash
--ligand-charge 'SAM:1,GPP:-3'    # SAM has charge +1, GPP has charge -3
--ligand-charge 'LIG:-2'          # LIG has charge -2
```

### Total charge override
```bash
-q 0                              # Force total system charge to 0
-q -1                             # Force total system charge to -1
```

### Charge resolution order
1. `-q/--charge` (explicit CLI override) -- highest priority
2. Pocket extraction (sums amino acids, ions, `--ligand-charge`)
3. `--ligand-charge` as fallback (when extraction skipped)
4. `.gjf` template metadata
5. Default: none (unresolved charge aborts; provide `-q` or `.gjf` charge metadata, or use PDB `--ligand-charge`)

```{note}
`--ligand-charge` derivation is only applied for PDB inputs (including XYZ/GJF inputs when `--ref-pdb` is supplied) and only when charge is **not yet resolved**. If a `.gjf` template already provides a charge value before `--ligand-charge` is evaluated, the template charge takes precedence and `--ligand-charge` will not override it.
```

```{tip}
Always provide `--ligand-charge` for non-standard residues (substrates, cofactors, unusual ligands) to ensure correct charge propagation.
```

---

## Ligand Charge Format

The `--ligand-charge` option supports two formats:

### Mapping format (recommended)
```bash
--ligand-charge 'SAM:1,GPP:-3'     # Per-residue name mapping
--ligand-charge 'LIG:-2'           # Single residue mapping
```

### Integer format
```bash
--ligand-charge -3                  # Total charge for all unknown residues
```

In the mapping format, residue names are matched case-insensitively. Unmapped non-standard residues default to charge 0.

---

## Spin Multiplicity

```bash
-m 1      # Singlet (default)
-m 2      # Doublet
-m 3      # Triplet
```

```{note}
In the `all` command, use `-m/--mult`. In other subcommands, use `-m/--multiplicity`.
```

---

## Atom Selectors

Atom selectors identify specific atoms for scans and restraints. They can be:

### Integer index (1-based by default)
```bash
--scan-lists '[(1, 5, 2.0)]'      # Atoms 1 and 5, target distance 2.0 A
```

### PDB-style selector string
```bash
--scan-lists '[("TYR,285,CA", "MMT,309,C10", 2.20)]'
```

Selector fields can be separated by:
- Space: `'TYR 285 CA'`
- Comma: `'TYR,285,CA'`
- Slash: `'TYR/285/CA'`
- Backtick: `` 'TYR`285`CA' ``
- Backslash: `'TYR\285\CA'`

The three tokens (residue name, residue number, atom name) can appear in any order -- the parser uses a fallback heuristic if the order is non-standard.

---

## Input File Requirements

### PDB files
- Must contain **hydrogen atoms** (use `reduce`, `pdb2pqr`, or Open Babel to add them)
- Must have **element symbols** in columns 77-78 (use `mlmm add-elem-info` if missing)
- Multiple PDBs must have **identical atoms in the same order** (only coordinates may differ)

### XYZ and GJF files
- Can be used when pocket extraction is skipped (omit `-c/--center`)
- `.gjf` files can provide charge/spin defaults from embedded metadata

### Amber topology files
- `--real-parm7`: Amber parm7 file containing force field parameters for the full system
- The parm7 must match the atom ordering of the input PDB exactly

---

## YAML Configuration

Advanced settings can be passed via `--args-yaml`:

```bash
mlmm all -i R.pdb P.pdb -c 'LIG' --args-yaml config.yaml
```

YAML values take **highest precedence**:
```
defaults -> CLI options -> YAML (wins)
```

See [YAML Reference](yaml-reference.md) for all available options.

---

## Output Directory

Use `--out-dir` to specify where results are saved:

```bash
--out-dir ./my_results/    # Custom output directory
```

Default output directories:
- `all`: `./result_all/`
- `extract`: current directory or specified `-o`
- `mm-parm`: current directory or specified `--out-prefix`
- `define-layer`: current directory or specified `-o`
- `opt`: `./result_opt/`
- `tsopt`: `./result_tsopt/`
- `path-opt`: `./result_path_opt/`
- `path-search`: `./result_path_search/`
- `scan`: `./result_scan/`
- `freq`: `./result_freq/`
- `irc`: `./result_irc/`
- `dft`: `./result_dft/`

---

## See Also

- [Getting Started](getting-started.md) -- Installation and first run
- [Troubleshooting](troubleshooting.md) -- Common errors and fixes
- [YAML Reference](yaml-reference.md) -- Complete configuration options
