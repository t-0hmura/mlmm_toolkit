# `extract`

## Overview

> **Summary:** Extract a cluster model (active-site pocket) from a protein-ligand PDB. Specify substrates with `-c` by residue name, residue ID, or a PDB path. Link hydrogens are added to cap cut bonds when `--add-linkH` is enabled. Use `--ligand-charge` for non-standard residue charges.

### At a glance
- **Input:** One or more complex PDBs with consistent atom ordering (ensemble mode supported).
- **Substrate selection (`-c`):** residue IDs (`A:123A`), residue names (`GPP,SAM`), or a substrate PDB that matches the complex coordinates.
- **Selection logic:** distance cutoff (`--radius`) plus optional hetero-hetero proximity (`--radius-het2het`) and peptide/disulfide/PRO safeguards.
- **Truncation & capping:** trims residues/segments and optionally adds link hydrogens (`--add-linkH`).
- **Charges:** unknown residues default to 0 unless `--ligand-charge` supplies a total charge or per-resname mapping.

`mlmm extract` creates an active-site pocket (cluster model) from a protein-ligand PDB. It selects residues near the substrate, truncates the model according to backbone/side-chain rules, optionally caps severed bonds with link hydrogens, and can process single structures or ensembles.

This is typically the **first step** in an mlmm_toolkit workflow, producing a smaller, computationally tractable model from a full protein-ligand complex.

If you run into misclassification (e.g., unusual residue/atom naming), see the appendix below on naming requirements and the internal reference lists.

## Minimal example

```bash
mlmm extract -i complex.pdb -c A:123 -o pocket.pdb -l -3
```

## Output checklist

- `pocket.pdb` (or custom path via `-o`)
- Charge summary logged via INFO (amino acids, ions, substrates, total)

## Common examples

1. Minimal (ID-based substrate) with explicit total ligand charge.

```bash
mlmm extract -i complex.pdb -c A:123 -o pocket.pdb -l -3
```

2. Substrate provided as a PDB; per-resname charge mapping.

```bash
mlmm extract -i complex.pdb -c substrate.pdb -o pocket.pdb \
 -l "GPP:-3,MMT:-1"
```

3. Multi-structure to single multi-MODEL output with hetero-hetero proximity.

```bash
mlmm extract -i complex1.pdb complex2.pdb -c A:123 \
 -o pocket_multi.pdb --radius-het2het 2.6 -l -3 --verbose
```

## Usage
```bash
mlmm extract -i COMPLEX.pdb [COMPLEX2.pdb...]
 -c SUBSTRATE_SPEC
 [-o POCKET.pdb [POCKET2.pdb...]]
 [--radius Å] [--radius-het2het Å]
 [--include-H2O/--no-include-H2O]
 [--exclude-backbone/--no-exclude-backbone]
 [--add-linkH/--no-add-linkH]
 [--selected-resn LIST]
 [-l, --ligand-charge MAP_OR_NUMBER]
 [--verbose/--no-verbose]
```

### Examples
```bash
# Minimal (ID-based substrate) with explicit total ligand charge
mlmm extract -i complex.pdb -c '123' -o pocket.pdb -l -3

# Substrate provided as a PDB; per-resname charge mapping (others remain 0)
mlmm extract -i complex.pdb -c substrate.pdb -o pocket.pdb -l 'GPP:-3,SAM:1'

# Name-based substrate selection including all matches (WARNING is logged)
mlmm extract -i complex.pdb -c 'GPP,SAM' -o pocket.pdb -l 'GPP:-3,SAM:1'

# Multi-structure to single multi-MODEL output with hetero-hetero proximity enabled
mlmm extract -i complex1.pdb complex2.pdb -c 'GPP,SAM' -o pocket_multi.pdb --radius-het2het 2.6 -l 'GPP:-3,SAM:1'

# Multi-structure to multiple outputs with hetero-hetero proximity enabled
mlmm extract -i complex1.pdb complex2.pdb -c 'GPP,SAM' -o pocket1.pdb pocket2.pdb --radius-het2het 2.6 -l 'GPP:-3,SAM:1'
```

## Workflow

### Residue inclusion
- Always include the substrate residues from `-c/--center`.
- **Standard cutoff (`--radius`, default 2.6 Å):**
 - When `--no-exclude-backbone`, any atom within the cutoff qualifies a residue.
 - When `--exclude-backbone`, amino-acid residues must contact the substrate with a **non-backbone** atom (not N/H*/CA/HA*/C/O). Non-amino acids use any atom.
- **Independent hetero-hetero cutoff (`--radius-het2het`):** adds residues when a substrate hetero atom (non C/H) lies within the specified Å of a protein hetero atom. With backbone exclusion enabled the protein atom must be non-backbone.
- **Water handling:** HOH/WAT/H2O/DOD/TIP/TIP3/SOL are included by default (`--include-H2O`).
- **Forced inclusion:** `--selected-resn` accepts IDs with chains/insertion codes (e.g., `A:123A`).
- **Neighbor safeguards:**
 - When backbone exclusion is off and a residue contacts the substrate with a backbone atom, auto-include the peptide-adjacent N/C neighbors (C-N <= 1.9 Å). Termini keep caps (N/H* or C/O/OXT).
 - Disulfide bonds (SG-SG <= 2.5 Å) bring both cysteines.
 - Non-terminal PRO residues always pull in the N-side amino acid; CA is preserved even if backbone atoms are removed, and when `--exclude-backbone`, the neighbor's C/O/OXT remain to maintain the peptide bond.

### Truncation/capping
- Isolated residues retain only side-chain atoms; amino-acid backbone atoms (N, CA, C, O, OXT plus N/CA hydrogens) are removed except for PRO/HYP safeguards.
- Continuous peptide stretches keep internal backbone atoms; only terminal caps (N/H* or C/O/OXT) are removed. TER awareness prevents capping across chain breaks.
- With `--exclude-backbone`, main-chain atoms on all **non-substrate** amino acids are stripped (subject to PRO/HYP safeguards and PRO neighbor retention).
- Non-amino-acid residues never lose atoms named like backbone (N/CA/HA/H/H1/H2/H3).

### Link hydrogens (`--add-linkH`)
- Adds carbon-only link hydrogens at 1.09 Å along severed bond vectors (CB-CA, CA-N, CA-C; PRO/HYP use CA-C only).
- Inserted after a `TER` as contiguous `HETATM` records named `HL` in residue `LKH` (chain `L`). Serial numbers continue from the main block.
- In multi-structure mode the same bonds are capped across all models; coordinates remain model-specific.

### Charge summary (`--ligand-charge`)
- Amino acids and common ions draw charges from internal dictionaries; waters are zero.
- Unknown residues default to 0 unless `--ligand-charge` supplies either a total charge (distributed across unknown substrate residues, or all unknowns when no unknown substrate) or a per-resname mapping like `GPP:-3,SAM:1`.
- Summaries (protein/ligand/ion/total) are logged for the first input when verbose mode is enabled.

### Substrate specification (`-c/--center`)
- PDB path: the coordinates must match the first input exactly (tolerance 1e-3 Å); residue IDs propagate to other structures.
- Residue IDs: `'123,124'`, `'A:123,B:456'`, `'123A'`, `'A:123A'` (insertion codes supported).
- Residue names: comma-separated list (case insensitive). If multiple residues share a name, **all** matches are included and a warning is logged.

### Multi-structure ensembles
- Accepts multiple input PDBs (identical atom ordering is validated at the head/tail of each file). Each structure is processed independently and the **union** of selected residues is applied to every model so that outputs remain consistent.
- Output policy:
 - No `-o`, multiple inputs -> per-file `pocket_<original_basename>.pdb`.
 - One `-o` path -> single multi-MODEL PDB.
 - N outputs where N == number of inputs -> N individual PDBs.
- Diagnostics echo raw vs. kept atom counts per model along with residue IDs.

## CLI options

> **Note:** Default values shown are used when the option is not specified.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | One or more protein-ligand PDB files (identical atom ordering required). | Required |
| `-c, --center SPEC` | Substrate specification (PDB path, residue IDs, or residue names). | Required |
| `-o, --output PATH...` | Pocket PDB output(s). One path => multi-MODEL, N paths => per input. | Auto (`pocket.pdb` or `pocket_<input>.pdb`) |
| `-r, --radius FLOAT` | Atom-atom distance cutoff (Å) for inclusion. | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero-hetero cutoff (Å, non C/H). | `0.0` (internally 0.001 Å when zero) |
| `--include-H2O/--no-include-H2O` | Include HOH/WAT/H2O/DOD/TIP/TIP3/SOL waters. | `True` |
| `--exclude-backbone/--no-exclude-backbone` | Remove backbone atoms on non-substrate amino acids (PRO/HYP safeguards). | `True` |
| `--add-linkH/--no-add-linkH` | Add carbon-only link hydrogens at 1.09 Å along severed bonds. | `False` |
| `--selected-resn TEXT` | Force-include residues (IDs with optional chains/insertion codes). | `""` |
| `-l, --ligand-charge TEXT` | Total charge or per-resname mapping (e.g., `GPP:-3,SAM:1`). | _None_ |
| `-v, --verbose/--no-verbose` | Emit INFO-level logging (`True`) or keep warnings only (`False`). | `True` |

## Outputs
```text
<output>.pdb # Pocket PDB(s) with optional link hydrogens after a TER record
 # Single input -> pocket.pdb by default
 # Multiple inputs without -o -> pocket_<original_basename>.pdb per structure
 # One -o path with multiple inputs -> single multi-MODEL PDB
 # Output directories are not created automatically; ensure they exist
```
- Charge summary (protein/ligand/ion/total) is logged for model #1 when verbose mode is enabled.
- Programmatic use (`extract_api`) returns `{"outputs": [...], "counts": [...], "charge_summary": {...}}`.

## Appendix: PDB naming requirements and reference lists

This appendix is mainly for debugging cases where `extract` misclassifies residues due to **non-standard residue/atom naming**. If your inputs follow standard PDB conventions, you can usually skip it.

```{important}
For `extract` to work correctly, **residue names and atom names in the input PDB must conform to standard PDB naming conventions**. The tool relies on internal dictionaries to recognize amino acids, ions, water molecules, and backbone atoms. Non-standard naming will cause residues to be misclassified or charges to be incorrectly assigned.
```

The following internal constants define the recognized names:

### `AMINO_ACIDS`

A dictionary mapping residue names to their nominal integer charges. Membership in this dictionary determines whether a residue is treated as an amino acid for backbone handling, truncation, and charge calculation.

**Standard 20 amino acids** (charges reflect physiological pH):
- Neutral: `ALA`, `ASN`, `CYS`, `GLN`, `GLY`, `HIS`, `ILE`, `LEU`, `MET`, `PHE`, `PRO`, `SER`, `THR`, `TRP`, `TYR`, `VAL`
- Positive (+1): `ARG`, `LYS`
- Negative (-1): `ASP`, `GLU`

**Canonical extras:**
- `SEC` (selenocysteine, 0), `PYL` (pyrrolysine, +1)

**Protonation/tautomer variants** (Amber/CHARMM style):
- `HIP` (+1, fully protonated His), `HID` (0, Nd-protonated His), `HIE` (0, Ne-protonated His)
- `ASH` (0, neutral Asp), `GLH` (0, neutral Glu), `LYN` (0, neutral Lys), `ARN` (0, neutral Arg)
- `TYM` (-1, deprotonated Tyr phenolate)

**Phosphorylated residues:**
- Dianionic (-2): `SEP`, `TPO`, `PTR`
- Monoanionic (-1): `S1P`, `T1P`, `Y1P`
- Phospho-His (phosaa19SB): `H1D` (0), `H2D` (-1), `H1E` (0), `H2E` (-1)

**Cysteine variants:**
- `CYX` (0, disulfide), `CSO` (0, sulfenic acid), `CSD` (-1, sulfinic acid), `CSX` (0, generic derivative)
- `OCS` (-1, cysteic acid), `CYM` (-1, deprotonated Cys)

**Lysine variants / carboxylation:**
- `MLY` (+1), `LLP` (+1), `KCX` (-1, Nz-carboxylic acid)

**D-amino acids** (19 residues):
- `DAL`, `DAR`, `DSG`, `DAS`, `DCY`, `DGN`, `DGL`, `DHI`, `DIL`, `DLE`, `DLY`, `MED`, `DPN`, `DPR`, `DSN`, `DTH`, `DTR`, `DTY`, `DVA`

**Other modified residues:**
- `CGU` (-2, gamma-carboxy-glutamate), `CGA` (-1), `PCA` (0, pyroglutamate), `MSE` (0, selenomethionine), `OMT` (0, methionine sulfone), `HYP` (0, hydroxyproline)
- Various others: `ASA`, `CIR`, `FOR`, `MVA`, `IIL`, `AIB`, `HTN`, `SAR`, `NMC`, `PFF`, `NFA`, `ALY`, `AZF`, `CNX`, `CYF`

**N-terminal variants** (prefix `N`): `NALA` (+1), `NARG` (+2), `NASP` (0), `NGLU` (0), `NLYS` (+2), etc., plus `ACE` (0), `NTER` (+1, generic)

**C-terminal variants** (prefix `C`): `CALA` (-1), `CARG` (0), `CASP` (-2), `CGLU` (-2), `CLYS` (0), etc., plus `NHE` (0), `NME` (0), `CTER` (-1, generic)

### `BACKBONE_ATOMS`

A set of atom names considered backbone atoms for amino acids. These are used when `--exclude-backbone` to determine which atoms to remove from non-substrate residues:

```
N, C, O, CA, OXT, H, H1, H2, H3, HN, HA, HA2, HA3
```

### `ION`

A dictionary mapping ion residue names to their formal charges. Recognized ions are automatically assigned correct charges in the charge summary.

| Charge | Residue Names |
|--------|---------------|
| +1 | `LI`, `NA`, `K`, `RB`, `CS`, `TL`, `AG`, `CU1`, `Ag`, `K+`, `Na+`, `NH4`, `H3O+`, `HE+`, `HZ+`, `Tl` |
| +2 | `MG`, `CA`, `SR`, `BA`, `MN`, `FE2`, `CO`, `NI`, `CU`, `ZN`, `CD`, `HG`, `PB`, `Be`, `PD`, `PT`, `Sn`, `Ra`, `YB2`, `V2+` |
| +3 | `FE`, `AU3`, `AL`, `GA`, `IN`, `CE`, `Ce`, `CR`, `Cr`, `Dy`, `EU`, `EU3`, `Er`, `GD3`, `LA`, `LU`, `Nd`, `PR`, `SM`, `Sm`, `TB`, `Tm`, `Y`, `Pu` |
| +4 | `U4+`, `Th`, `Hf`, `Zr` |
| -1 | `F`, `CL`, `BR`, `I`, `Cl-`, `IOD` |

### `WATER_RES`

A set of residue names recognized as water molecules. Waters are included by default (`--include-H2O`) and assigned zero charge:

```
HOH, WAT, H2O, DOD, TIP, TIP3, SOL
```

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- `--radius` defaults to 2.6 Å; `0` is nudged to 0.001 Å to avoid empty selections. `--radius-het2het` is off by default (also nudged to 0.001 Å when zero is provided).
- Waters can be excluded with `--no-include-H2O`.
- Backbone trimming plus capping respect chain breaks and PRO/HYP safeguards as outlined above; non-amino residues never lose backbone-like atom names.
- Link hydrogens are inserted only on carbon cuts and reuse identical bonding patterns across models in ensemble mode.
- INFO logs summarize residue selection, truncation counts, and charge breakdowns.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [Getting Started](getting_started.md) -- Installation and first run
- [Concepts](concepts.md) -- Full system vs. pocket (cluster model)
- [CLI Conventions](cli_conventions.md) -- Residue selectors and charge specification
- [mm-parm](mm_parm.md) -- Generate Amber topology from the extracted pocket
- [define-layer](define_layer.md) -- Assign 3-layer ML/MM partitioning
