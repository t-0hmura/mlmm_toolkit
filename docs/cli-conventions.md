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
All subcommands (including `extract` and `fix-altloc`) use Click as their CLI backend.

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

## Show Config

View the current configuration (useful for verifying YAML overrides):

```bash
mlmm opt -i input.pdb --parm real.parm7 -q -1 --show-config --dry-run
```

---

## ML/MM Required Options

Most subcommands (except `all`, `extract`, `mm-parm`, and `define-layer`) require two additional options for ML/MM calculations:

```bash
--parm real.parm7 # Amber parm7 topology file for the full (real) system
--model-pdb model.pdb # PDB file defining the ML (model) region atoms
```

The `all` command generates these files automatically during the workflow (via `mm-parm` and `define-layer`). When running subcommands individually, you must provide them explicitly.

```bash
# Example: running path-search individually
mlmm path-search -i R.pdb P.pdb --parm real.parm7 --model-pdb model.pdb -q 0 -m 1
```

---

## B-factor Layer Encoding

mlmm-toolkit uses the PDB B-factor column to encode the 3-layer ML/MM partitioning:

| Layer | B-factor | Meaning |
|-------|----------|---------|
| ML | 0.0 | MLIP energy/force/Hessian |
| Movable-MM | 10.0 | MM atoms allowed to move |
| Frozen | 20.0 | Coordinates fixed |

The `define-layer` subcommand writes these B-factors into the PDB. You can inspect layer assignments in any molecular viewer that supports B-factor coloring.

A tolerance of 1.0 is used when reading B-factors, so values near 0/10/20 are mapped to ML/Movable/Frozen.

### Layer definition methods

1. **`define-layer` subcommand** (recommended):
   ```bash
   mlmm define-layer -i system.pdb --model-pdb ml_region.pdb -o labeled.pdb
   ```

2. **Distance cutoffs** (YAML/CLI):
   ```yaml
   calc:
     hess_cutoff: 3.6       # Distance cutoff for Hessian-target MM atoms
     movable_cutoff: 8.0    # Distance cutoff for Movable-MM (beyond = Frozen)
   ```

3. **Read from B-factors**:
   ```yaml
   calc:
     use_bfactor_layers: true   # Read layers from input PDB B-factors
   ```

4. **Explicit index specification** (YAML):
   ```yaml
   calc:
     hess_mm_atoms: [100, 101, 102, ...]
     movable_mm_atoms: [200, 201, 202, ...]
     frozen_mm_atoms: [300, 301, 302, ...]
   ```

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

For PDB inputs, `--ligand-charge` lets you specify charges only for non-standard residues (substrates, cofactors, metal ions). The net system charge is then **automatically derived** by summing standard amino-acid charges, ion charges, and your ligand charges -- no need to manually count atoms across the entire complex. This is especially useful for large enzyme-substrate systems where the net system charge is not obvious.

### Per-residue mapping (recommended)
```bash
-l 'SAM:1,GPP:-3' # SAM has charge +1, GPP has charge -3
-l 'LIG:-2' # LIG has charge -2
```

### Total charge override
```bash
-q 0 # Force net system charge to 0
-q -1 # Force net system charge to -1
```

### Charge resolution order
1. `-q/--charge` (explicit CLI override) -- highest priority.
2. ML-region determination charge summary (sum of amino acids, ions, and `--ligand-charge`).
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
-l 'SAM:1,GPP:-3' # Per-residue name mapping
-l 'SAM=1,GPP=-3' # Same meaning (= separator)
-l 'LIG:-2' # Single residue mapping
```

Both colon (`:`) and equals (`=`) separators are accepted.

### Integer format
```bash
-l -3 # Total charge for all unknown residues
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
--scan-lists '[(1, 5, 2.0)]' # Atoms 1 and 5, target distance 2.0 Å
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

### XYZ files
- Can be used when ML-region determination is skipped (omit `-c/--center`)

### Amber topology files
- `--parm`: Amber parm7 file containing force field parameters for the full system
- The parm7 must match the atom ordering of the input PDB exactly

---

## Backend Selection

All computation subcommands (`opt`, `tsopt`, `freq`, `irc`, `dft`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, `all`) accept:

| Option | Description | Default |
|--------|-------------|---------|
| `-b, --backend` | MLIP backend for the ML region: `uma`, `orb`, `mace`, `aimnet2` | `uma` |
| `--embedcharge/--no-embedcharge` | Enable xTB point-charge embedding correction | `--no-embedcharge` |

Alternative backends are installed via optional dependency groups:

```bash
pip install "mlmm-toolkit[orb]"       # ORB backend
pip install "mlmm-toolkit[aimnet]"   # AIMNet2 backend
pip install --no-deps mace-torch      # MACE backend
```

---

## YAML Configuration

Advanced settings can be passed via layered YAML inputs.

Precedence:
```
defaults < config < CLI options < override-yaml
```

See [YAML Reference](yaml-reference.md) for all available options.

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
- `scan2d`: `./result_scan2d/`
- `scan3d`: `./result_scan3d/`
- `freq`: `./result_freq/`
- `irc`: `./result_irc/`
- `dft`: `./result_dft/`

---

## See Also

- [Getting Started](getting-started.md) -- Installation and first run
- [Concepts & Workflow](concepts.md) -- ML/MM 3-layer system, ONIOM decomposition overview
- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Common errors and fixes
- [YAML Reference](yaml-reference.md) -- Complete configuration options
