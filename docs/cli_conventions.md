# CLI Conventions

This page documents the conventions used across all `mlmm` commands. Understanding these conventions helps you write correct commands and avoid common errors.

---

## Boolean Options

Boolean options are normalized at the root CLI.
Both notations are accepted:

```bash
# Recommended
--tsopt --thermo --no-dft

# Also accepted
--tsopt True --thermo yes --dft 0
```

For options that are defined only as `--flag`, the root CLI also accepts `--no-flag` and `--flag False` as compatibility aliases.
`extract` and `fix-altloc` are parser-wrapper subcommands (argparse backend), but the same bool normalization still applies at the root CLI.
For these parser wrappers, root normalization discovers bool options from the argparse definitions, so `--flag/--no-flag` and `--flag True/False` stay in sync without a separate hand-maintained map.

Common boolean options:
- `--tsopt`, `--thermo`, `--dft` -- enable post-processing stages
- `--dump` -- write trajectory files
- `--preopt`, `--endopt` -- pre/post optimization toggles
- `--climb` -- enable climbing image in MEP search

---

## Progressive Help (`all`)

`mlmm all` uses two help levels:

```bash
mlmm all --help # core options only
mlmm all --help-advanced # full option list
```

`scan`, `scan2d`, `scan3d`, the calculation commands (`opt`, `path-opt`, `path-search`, `tsopt`, `freq`, `irc`, `dft`), and selected utility commands (`mm-parm`, `define-layer`, `add-elem-info`, `trj2fig`, `energy-diagram`, `oniom-export`) now follow the same progressive-help pattern (`--help` core, `--help-advanced` full). `extract` and `fix-altloc` also support progressive help (`--help` core, `--help-advanced` full parser options).

---

## Init Template

Generate a starter YAML and run a parse-only check:

```bash
mlmm init --out mlmm_all.config.yaml
mlmm all --config mlmm_all.config.yaml --dry-run
```

---

## ML/MM Required Options

Most subcommands (except `all`, `extract`, `mm-parm`, and `define-layer`) require two additional options for ML/MM calculations:

```bash
--real-parm7 real.parm7 # Amber parm7 topology file for the full (real) system
--model-pdb model.pdb # PDB file defining the ML (model) region atoms
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
-c 'SAM,GPP' # Select all residues named SAM or GPP
-c 'LIG' # Select all residues named LIG
```

### By residue ID
```bash
-c '123,456' # Residues 123 and 456
-c 'A:123,B:456' # Chain A residue 123, Chain B residue 456
-c '123A' # Residue 123 with insertion code A
-c 'A:123A' # Chain A, residue 123, insertion code A
```

### By PDB file
```bash
-c substrate.pdb # Use coordinates from a separate PDB to locate substrates
```

```{note}
When selecting by residue name, if multiple residues share the same name, **all** matches are included and a warning is logged.
```

---

## Charge Specification

For PDB inputs, `--ligand-charge` lets you specify charges only for non-standard residues (substrates, cofactors, metal ions). The total system charge is then **automatically derived** by summing standard amino-acid charges, ion charges, and your ligand charges -- no need to manually count atoms across the entire complex. This is especially useful for large enzyme-substrate systems where the total charge is not obvious.

### Per-residue mapping (recommended)
```bash
--ligand-charge 'SAM:1,GPP:-3' # SAM has charge +1, GPP has charge -3
--ligand-charge 'LIG:-2' # LIG has charge -2
```

### Total charge override
```bash
-q 0 # Force total system charge to 0
-q -1 # Force total system charge to -1
```

### Charge resolution order
1. `-q/--charge` (explicit CLI override) -- highest priority.
2. Pocket extraction charge summary (sum of amino acids, ions, and `--ligand-charge`).
3. `--ligand-charge` fallback when extraction is skipped.
4. Default: none (abort if unresolved).

Calculation subcommands (`scan`, `scan2d`, `scan3d`, `opt`, `path-opt`, `path-search`, `tsopt`, `freq`, `irc`, `dft`, `oniom-export`) still require explicit `-q/--charge`.

```{tip}
Always provide `--ligand-charge` for non-standard residues (substrates, cofactors, unusual ligands) to ensure correct charge propagation.
```

---

## Ligand Charge Format

The `--ligand-charge` option supports two formats:

### Mapping format (recommended)
```bash
--ligand-charge 'SAM:1,GPP:-3' # Per-residue name mapping
--ligand-charge 'LIG:-2' # Single residue mapping
```

### Integer format
```bash
--ligand-charge -3 # Total charge for all unknown residues
```

In the mapping format, residue names are matched case-insensitively. Unmapped non-standard residues default to charge 0.

---

## Spin Multiplicity

```bash
-m 1 # Singlet (default)
-m 2 # Doublet
-m 3 # Triplet
```

```{note}
Use `-m/--multiplicity` consistently in `all` and the calculation subcommands. `mm-parm` uses `--ligand-mult` for residue multiplicity metadata, which is a separate option.
```

---

## Atom Selectors

Atom selectors identify specific atoms for scans and restraints. They can be:

### Integer index (1-based by default)
```bash
--scan-lists '[(1, 5, 2.0)]' # Atoms 1 and 5, target distance 2.0 A
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

Advanced settings can be passed via layered YAML inputs:

```bash
```

Precedence:
```
defaults < config < CLI options < override-yaml
```

See [YAML Reference](yaml_reference.md) for all available options.

---

## Output Directory

Use `--out-dir` to specify where results are saved:

```bash
--out-dir ./my_results/ # Custom output directory
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

- [Getting Started](getting_started.md) -- Installation and first run
- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Common errors and fixes
- [YAML Reference](yaml_reference.md) -- Complete configuration options
