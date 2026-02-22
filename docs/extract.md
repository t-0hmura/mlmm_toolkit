# extract

## Overview

The `extract` subcommand extracts an **active-site binding pocket** (cluster model) from a protein-ligand complex PDB. It identifies residues near specified substrate(s), applies biochemically aware truncation (backbone/side-chain capping with safeguards), and can append link hydrogens for cut bonds. The extracted pocket serves as the starting structure for downstream ML/MM calculations.

This is typically the **first step** in an mlmm_toolkit workflow, producing a smaller, computationally tractable model from a full protein-ligand complex.

---

## Usage

```bash
mlmm extract -i INPUT.pdb [INPUT2.pdb ...] -c <substrate_spec> \
    [-o OUTPUT.pdb ...] [-r <A>] [--radius-het2het <A>] \
    [--include-H2O/--no-include-H2O] [--exclude-backbone/--no-exclude-backbone] \
    [--add-linkH/--no-add-linkH] [--selected-resn "CHAIN:RES" ...] \
    [--ligand-charge <number|"RES:Q,...">] [--verbose/--no-verbose]
```

---

## Examples

```bash
# Minimal (ID-based substrate) with explicit total ligand charge
mlmm extract -i complex.pdb -c A:123 -o pocket.pdb --ligand-charge -3

# Substrate provided as a PDB; per-resname charge mapping
mlmm extract -i complex.pdb -c substrate.pdb -o pocket.pdb \
    --ligand-charge "GPP:-3,MMT:-1"

# Name-based substrate selection including all matches
mlmm extract -i complex.pdb -c "GPP,MMT" -o pocket.pdb --ligand-charge -4

# Multi-structure to single multi-MODEL output with hetero-hetero proximity
mlmm extract -i complex1.pdb complex2.pdb -c A:123 \
    -o pocket_multi.pdb --radius-het2het 2.6 --ligand-charge -3 --verbose
```

---

## CLI Options

| Option | Default | Description |
|--------|---------|-------------|
| `-i, --input PATH...` | (required) | One or more input PDB files (full protein-ligand complex). |
| `-c, --center TEXT` | (required) | Substrate specification: residue names (`'SAM,GPP'`), residue IDs (`'A:123,B:456'`), or a PDB file path. |
| `-o, --output PATH...` | auto | Output PDB path(s). One path produces a multi-MODEL PDB; N paths (matching N inputs) produce per-structure PDBs. Omitted: `pocket_<source>.pdb`. |
| `-r, --radius FLOAT` | 2.6 | Cutoff distance (A) for including residues around the substrate. |
| `--radius-het2het FLOAT` | 0 (off) | Independent cutoff for hetero-atom to hetero-atom proximity. |
| `--include-H2O/--no-include-H2O` | true | Include water molecules (HOH/WAT/TIP3/SOL) in the pocket. |
| `--exclude-backbone/--no-exclude-backbone` | true | Exclude backbone atoms from non-substrate amino acids. |
| `--add-linkH/--no-add-linkH` | false | Add link hydrogens (HL in residue LKH) to cap severed bonds. |
| `--selected-resn TEXT...` | (none) | Force-include residues by chain:ID (e.g., `'A:123,B:456'`). |
| `--ligand-charge TEXT` | (none) | Charge specification: mapping (`'SAM:1,GPP:-3'`) or total integer. |
| `--verbose/--no-verbose` | true | Enable verbose logging of residue selection and atom counts. |

---

## Substrate Specification (`-c/--center`)

The `--center` option accepts three forms:

### By PDB file
```bash
-c substrate.pdb      # Exact coordinate match on the first input; IDs propagated to others
```

### By residue ID
```bash
-c '123,456'          # Residues 123 and 456
-c 'A:123,B:456'      # Chain A residue 123, Chain B residue 456
-c '123A'             # Residue 123 with insertion code A
```

### By residue name
```bash
-c 'GPP,MMT'          # All residues named GPP or MMT (case-insensitive)
```

If multiple residues share the same name, **all** matches are included and a WARNING is logged.

---

## Residue Inclusion Logic

1. **Substrate residues** are always included.
2. **Standard cutoff** (`--radius`, default 2.6 A):
   - With `--exclude-backbone` (default): amino-acid residues qualify only if a **non-backbone** atom is within the cutoff.
   - With `--no-exclude-backbone`: any atom within the cutoff qualifies the residue.
3. **Hetero-hetero proximity** (`--radius-het2het`): adds residues if a substrate hetero atom (non-C/H) is within the cutoff of a protein hetero atom.
4. **Waters** are included by default (`--include-H2O`).
5. **`--selected-resn`** force-includes residues regardless of distance.
6. **Disulfide safeguard**: if a selected CYS/CYX forms an SG-SG contact <= 2.5 A, both partners are included.
7. **Proline safeguard**: if a selected PRO is not N-terminal, the immediately preceding amino acid is included.

---

## Truncation (Capping)

- **Isolated residues**: keep pure side-chain only (remove backbone atoms N, CA, C, O, OXT and associated H).
  - PRO/HYP retain N, CA, HA, H to keep the ring intact.
- **Continuous peptide stretches**: keep internal backbone; only terminal caps are removed.
- With `--exclude-backbone` (default): delete main-chain atoms on all non-substrate amino acids (with PRO/HYP exceptions).

---

## Link Hydrogens (`--add-linkH`)

When enabled (`--add-linkH`), link hydrogens are added at **1.09 A** along the cut-bond vector (carbon-only caps). They are written as **HETATM** records with:
- Atom name: `HL`
- Residue name: `LKH`
- Chain: `L`

In multi-structure mode, link-H targets and ordering are enforced to be identical across models; coordinates remain model-specific.

---

## Charge Summary

The extraction process computes a charge summary:
- **Amino acids**: nominal integer charges from a built-in dictionary (standard + modified residues).
- **Ions**: built-in charges for common ions (ZN, MG, FE2, etc.). Waters are 0.
- **Unknown residues**: 0 unless `--ligand-charge` is given.
  - `--ligand-charge <number>`: total charge distributed across unknown substrate residues.
  - `--ligand-charge "RES1:Q1,RES2:Q2"`: per-residue-name charges; other unknowns remain 0.

---

## Outputs

| Output | Description |
|--------|-------------|
| `pocket.pdb` | Default single-input pocket (or custom path via `-o`). |
| `pocket_<source>.pdb` | Default per-input pocket when multiple inputs and `-o` is omitted. |
| Multi-MODEL PDB | When one `-o` path is given for multiple inputs. |
| Charge summary | Logged via INFO (amino acids, ions, substrates, total). |

The link-H block (when added) follows a TER record as contiguous HETATM records.

---

## Multi-structure Mode

- Accepts multiple input PDBs (same atom count; ordering assumed identical).
- Each structure is selected independently; the **union** of selected residues is applied to all.
- Disulfides, PRO-adjacency, and backbone-contact neighbor augmentation are also unioned.
- Outputs: one `-o` path produces a multi-MODEL PDB; N paths produce N single-model PDBs.
- Atom counts (raw vs after truncation) are logged per model.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [Getting Started](getting_started.md) -- Installation and first run
- [Concepts](concepts.md) -- Full system vs. pocket (cluster model)
- [CLI Conventions](cli_conventions.md) -- Residue selectors and charge specification
- [mm-parm](mm_parm.md) -- Generate Amber topology from the extracted pocket
- [define-layer](define_layer.md) -- Assign 3-layer ML/MM partitioning
