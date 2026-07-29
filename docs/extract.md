# `extract`

`mlmm extract` carves an active-site pocket from a protein–ligand PDB to define the ML region (and the surrounding MM environment for downstream stages). In `mlmm all`, this selection stage is managed internally. For reusable manual files, first run `mm-parm --out-prefix system`, then extract from its topology-matched `system.pdb`. The command selects residues near the substrate, truncates the model according to backbone / side-chain rules, optionally caps severed bonds with link hydrogens, and accepts either a single structure or a multi-structure ensemble (details under [Multi-structure ensembles](#multi-structure-ensembles)). Pick a substrate-selection mode (residue IDs, a substrate PDB, or residue names) depending on what you have on hand — see [Input syntax](#input-syntax) for the exact grammar.

If you run into misclassification (e.g. unusual residue / atom naming), see the appendix below on naming requirements and the internal reference lists.

## Examples

Minimal run with an ID-based substrate and an explicit total ligand charge:

```bash
# Minimal (ID-based substrate) with explicit total ligand charge
mlmm extract -i complex.pdb -c A:123 -o pocket.pdb -l -3
```

Substrate supplied as a PDB, with a per-resname charge mapping:

```bash
# Substrate provided as a PDB; per-resname charge mapping (others remain 0)
mlmm extract -i complex.pdb -c substrate.pdb -o pocket.pdb -l "GPP:-3,MMT:-1"
# name-based selection includes all matches (WARNING logged): -c 'GPP,SAM'
```

Multi-structure ensemble collapsed into one multi-MODEL output, using hetero-hetero proximity:

```bash
# Multi-structure → single multi-MODEL output with hetero-hetero proximity
mlmm extract -i complex1.pdb complex2.pdb -c A:123 \
    -o pocket_multi.pdb --radius-het2het 2.6 -l -3 --verbose 3
```

## Workflow

### Residue inclusion

- Always include the substrate residues from `-c/--center`.
- **Standard cutoff (`--radius`, default 2.6 Å)**: with `--no-exclude-backbone` (default), any atom within the cutoff qualifies a residue. With `--exclude-backbone`, amino-acid residues must contact the substrate with a **non-backbone** atom (not N / H* / CA / HA* / C / O / OXT). Non-amino acids always use any atom.
- **Independent hetero-hetero cutoff (`--radius-het2het`)**: adds residues when a substrate hetero atom (non C / H) lies within the specified Å of a protein hetero atom. With backbone exclusion enabled, the protein atom must be non-backbone.
- **Water handling**: HOH / WAT / H2O / DOD / TIP / TIP3 / SOL are included by default (`--include-h2o`).
- **Forced inclusion**: `--selected-resn` accepts residue IDs with chains / insertion codes (e.g. `A:123A`).
- **Neighbor safeguards**:
  - When backbone exclusion is off and a residue contacts the substrate with a backbone atom, the peptide-adjacent N / C neighbors (C–N ≤ 1.9 Å) are auto-included; termini keep caps (N/H* or C/O/OXT).
  - Disulfide bonds (SG–SG ≤ 2.5 Å) bring both cysteines.
  - Non-terminal PRO residues always pull in the N-side amino acid; CA is preserved even when backbone atoms are removed, and under `--exclude-backbone` the neighbor's C / O / OXT remain to maintain the peptide bond.

### Truncation and capping

- Isolated residues retain only side-chain atoms; amino-acid backbone atoms (N, CA, C, O, OXT plus N/CA hydrogens) are removed except for PRO / HYP safeguards.
- Continuous peptide stretches keep internal backbone atoms; only terminal caps (N/H* or C/O/OXT) are removed. TER awareness prevents capping across chain breaks.
- With `--exclude-backbone`, main-chain atoms on all **non-substrate** amino acids are stripped (subject to PRO / HYP safeguards and PRO neighbor retention).
- Non-amino-acid residues never lose atoms named like backbone (N / CA / HA / H / H1 / H2 / H3).

### Link hydrogens (`--add-linkh`)

- Carbon-only link hydrogens are placed at 1.09 Å along severed bond vectors (CB–CA, CA–N, CA–C; PRO / HYP use CA–C only).
- Inserted after a `TER` as contiguous `HETATM` records named `HL` in residue `LKH` (chain `L`). Serial numbers continue from the main block.
- In multi-structure mode the same bonds are capped across all models; coordinates remain model-specific.

### Charge summary (`--ligand-charge`)

Amino acids and common ions draw charges from internal dictionaries; waters are zero. Unknown residues default to 0 unless `--ligand-charge` supplies either a total charge (distributed across unknown substrate residues, or all unknowns when no unknown substrate) or a per-resname mapping like `GPP:-3,SAM:1`. Summaries (protein / ligand / ion / total) are logged for the first input when verbose mode is enabled.

### Multi-structure ensembles

`extract` accepts multiple input PDBs and compares the complete ordered atom-identity sequence across every file. Each structure is processed independently and the **union** of selected residues is applied to every model so outputs stay consistent.

| Output policy | Layout |
|---|---|
| No `-o`, multiple inputs | per-file `pocket_<original_basename>.pdb` |
| One `-o` path | single multi-MODEL PDB |
| N outputs matching N inputs | N individual PDBs |

Diagnostics echo raw vs. kept atom counts per model along with residue IDs.

## Outputs

```text
<output>.pdb        # Pocket PDB(s) with optional link hydrogens after a TER record.
                    # See the "Multi-structure ensembles" output-policy table above
                    # for the -o / multi-input naming rules.
                    # Parent directories of -o are created automatically if missing.
```

Programmatic use (`extract_api`) returns `{"outputs": [...], "counts": [...], "charge_summary": {...}}` (the verbose-mode charge summary is described under Workflow > Charge summary).

## CLI options

Command form:

```bash
mlmm extract -i COMPLEX.pdb [COMPLEX2.pdb ...]
    -c SUBSTRATE_SPEC
    [-o POCKET.pdb [POCKET2.pdb ...]]
    [--radius Å] [--radius-het2het Å]
    [--include-h2o / --no-include-h2o]
    [--exclude-backbone / --no-exclude-backbone]
    [--add-linkh / --no-add-linkh]
    [--selected-resn LIST]
    [-l, --ligand-charge MAP_OR_NUMBER]
    [-v LEVEL]
```

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | One or more protein–ligand PDB files (identical atom ordering required). | Required |
| `-c, --center SPEC` | Substrate specification (PDB path, residue IDs, or residue names). | Required |
| `-o, --output PATH...` | Pocket PDB output(s). One path ⇒ multi-MODEL; N paths ⇒ per input. | Auto (`pocket.pdb` or `pocket_<input>.pdb`) |
| `-r, --radius FLOAT` | Atom-atom distance cutoff (Å) for inclusion. | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero-hetero cutoff (Å, non C / H). | `0.0` |
| `--include-h2o / --no-include-h2o` | Include HOH / WAT / H2O / DOD / TIP / TIP3 / SOL waters. | `True` |
| `--exclude-backbone / --no-exclude-backbone` | Remove backbone atoms on non-substrate amino acids (PRO / HYP safeguards). | `False` |
| `--add-linkh / --no-add-linkh` | Add carbon-only link hydrogens at 1.09 Å along severed bonds (distance-based). Not needed for an mlmm `--model-pdb` (the ML/MM calculator caps the boundary from the `--parm` topology); use for standalone pocket models only. | `False` |
| `--selected-resn TEXT` | Force-include residues (IDs with optional chains / insertion codes). | `""` |
| `--modified-residue TEXT` | Comma-separated residue names (with optional charge) to treat as amino acids for backbone truncation and charge assignment (e.g. `HD1,HD2,HD3` or `HD1:0,SEP:-2`). Useful for modified amino acids with non-standard names. | `""` |
| `-l, --ligand-charge TEXT` | Total charge or per-resname mapping (e.g. `GPP:-3,SAM:1`). | _None_ |

### Input syntax

Substrate specification (`-c/--center`):

- **PDB path**: coordinates must match the first input exactly (tolerance 1e-3 Å); residue IDs propagate to other structures.
- **Residue IDs**: `'123,124'`, `'A:123,B:456'`, `'123A'`, `'A:123A'` (insertion codes supported).
- **Residue names**: comma-separated, case-insensitive. If multiple residues share a name, **all** matches are included and a warning is logged.

```{tip}
If the extracted pocket is too small, calculated energies and barriers may be unreliable — increasing the extraction radius (e.g. `-r 4.0` or higher) improves accuracy by including more of the protein environment.
```

## Notes

### Systems with non-standard residues (MCPB, etc.)

When metal-coordinating amino-acid parameters are generated by tools such as Amber's `MCPB.py` (Metal Center Parameter Builder), the coordinating residues are assigned non-standard names (e.g. `HD1`, `HE1`, `CM1`, `AP1`). These are not in `extract`'s internal `AMINO_ACIDS` dictionary, so **backbone truncation and link-hydrogen capping will not be applied correctly**, and a warning is emitted:

```text
[extract] WARNING: Residue HD1 83 may be an amino acid (has N, CA, C, O)
but is not recognized as a standard residue name.
Backbone truncation was not applied.
Consider preparing the pocket model manually.
```

```{tip}
Register these names with `--modified-residue` (e.g. `--modified-residue HD1,HE1,CM1,AP1`); append `:charge` to set an integer charge, e.g. `--modified-residue HD1:0,SEP:-2` (charge defaults to `0`). `extract` then treats them as amino acids and applies backbone truncation, link-hydrogen capping, and charge assignment automatically, and the warning above is suppressed.
```

```{important}
If `--modified-residue` cannot cover your case (e.g. unusual backbone topology), **construct the pocket model manually**:

1. Select residues around the active site and determine truncation points.
2. Add a link hydrogen on the parent atom (the atom that remains) of each severed covalent bond.
3. Use residue name `LKH` (chain `L`) and atom name `HL` for the link hydrogen.
4. Place it at **1.09 Å** along the original bond direction.
```

## Appendix: PDB naming requirements and reference lists

This appendix exists mainly for debugging cases where `extract` misclassifies residues due to **non-standard residue or atom naming**. If your inputs follow standard PDB conventions, you can usually skip it.

```{important}
For `extract` to work correctly, **residue and atom names in the input PDB must conform to standard PDB naming conventions**. The tool relies on internal dictionaries to recognize amino acids, ions, water molecules, and backbone atoms. Non-standard naming will cause residues to be misclassified or charges to be incorrectly assigned.
```

### `AMINO_ACIDS`

A dictionary mapping residue names to their nominal integer charges. Membership determines whether a residue is treated as an amino acid for backbone handling, truncation, and charge calculation.

**Standard 20** (charges reflect physiological pH):

- Neutral: `ALA`, `ASN`, `CYS`, `GLN`, `GLY`, `HIS`, `ILE`, `LEU`, `MET`, `PHE`, `PRO`, `SER`, `THR`, `TRP`, `TYR`, `VAL`
- Positive (+1): `ARG`, `LYS`
- Negative (−1): `ASP`, `GLU`

**Canonical extras:** `SEC` (selenocysteine, 0), `PYL` (pyrrolysine, 0).

**Protonation / tautomer variants** (Amber / CHARMM): `HIP` (+1, fully protonated His), `HID` (0, Nδ-protonated His), `HIE` (0, Nε-protonated His), `ASH` (0, neutral Asp), `GLH` (0, neutral Glu), `LYN` (0, neutral Lys), `ARN` (0, neutral Arg), `TYM` (−1, deprotonated Tyr phenolate).

**Phosphorylated:** dianionic (−2) `SEP`, `TPO`, `PTR`; monoanionic (−1) `S1P`, `T1P`, `Y1P`; phospho-His (phosaa19SB) `H1D` (0), `H2D` (−1), `H1E` (0), `H2E` (−1).

**Cysteine variants:** `CYX` (0, disulfide), `CSO` (0, sulfenic acid), `CSD` (−1, sulfinic acid), `CSX` (0, generic), `OCS` (−1, cysteic acid), `CYM` (−1, deprotonated Cys).

**Lysine variants / carboxylation:** `MLY` (+1), `LLP` (0), `DLY` (+1), `KCX` (−1, Nz-carboxylic acid).

**D-amino acids** (19): `DAL`, `DAR`, `DSG`, `DAS`, `DCY`, `DGN`, `DGL`, `DHI`, `DIL`, `DLE`, `DLY`, `MED`, `DPN`, `DPR`, `DSN`, `DTH`, `DTR`, `DTY`, `DVA`.

**Other modified:** `CGU` (−2, γ-carboxy-glutamate), `CGA` (−1), `PCA` (0, pyroglutamate), `MSE` (0, selenomethionine), `OMT` (0, methionine sulfone), `HYP` (0, hydroxyproline); also `ASA`, `CIR`, `FOR`, `MVA`, `IIL`, `AIB`, `HTN`, `SAR`, `NMC`, `PFF`, `NFA`, `ALY`, `AZF`, `CNX`, `CYF`.

**N-terminal variants** (`N` prefix): `NALA` (+1), `NARG` (+2), `NASP` (0), `NGLU` (0), `NLYS` (+2), … plus `ACE` (0), `NTER` (+1, generic).
**C-terminal variants** (`C` prefix): `CALA` (−1), `CARG` (0), `CASP` (−2), `CGLU` (−2), `CLYS` (0), … plus `NHE` (0), `NME` (0), `CTER` (−1, generic).

### `BACKBONE_ATOMS`

Atom names treated as backbone for amino acids; under `--exclude-backbone` these are removed from non-substrate residues:

```
N, C, O, CA, OXT, H, H1, H2, H3, HN, HA, HA2, HA3
```

### `ION`

Recognized ion residue names with formal charges:

| Charge | Residue names |
|---|---|
| +1 | `LI`, `NA`, `K`, `RB`, `CS`, `TL`, `AG`, `CU1`, `K+`, `NA+`, `NH4`, `H3O+`, `HE+`, `HZ+` |
| +2 | `MG`, `CA`, `SR`, `BA`, `MN`, `FE2`, `CO`, `NI`, `CU`, `ZN`, `CD`, `HG`, `PB`, `BE`, `PD`, `PT`, `SN`, `RA`, `YB2`, `V2+` |
| +3 | `FE`, `AU3`, `AL`, `GA`, `IN`, `CE`, `CR`, `DY`, `EU`, `EU3`, `ER`, `GD3`, `LA`, `LU`, `ND`, `PR`, `SM`, `TB`, `TM`, `Y`, `PU` |
| +4 | `U4+`, `TH`, `HF`, `ZR` |
| −1 | `F`, `CL`, `BR`, `I`, `CL-`, `IOD` |

### `WATER_RES`

Recognized water residue names (included by default with `--include-h2o`, assigned zero charge):

```
HOH, WAT, H2O, DOD, TIP, TIP3, SOL
```

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide
- [Getting Started](getting-started.md) — Installation and first run
- [Concepts](concepts.md) — Full system vs. ML region
- [CLI Conventions](cli-conventions.md) — Residue selectors and charge specification
- [mm-parm](mm-parm.md) — Generate the full-system Amber topology and the
  topology-matched PDB to use for standalone extraction/layering
- [define-layer](define-layer.md) — Assign 3-layer ML/MM partitioning
