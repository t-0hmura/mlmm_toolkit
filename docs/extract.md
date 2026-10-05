# `extract` (cut out the ML region)

## Overview

`extract` cuts the residues around a substrate out of a protein–ligand PDB/mmCIF file, cuts the main chain by fixed rules, and counts the charge of the cut-out region. In mlmm-toolkit this region becomes the ML region; the rest of the protein stays in the calculation as MM.

### What it is for

* **Defining the ML region**: write the atom selection that `define-layer` and the calculation commands take as `--model-pdb`.
* **Checking the ML region before a run**: see which residues fall in the region, how many atoms it has, and its charge.
* **Cutting several states the same way**: give reactant and product structures with the same atom order in one run, and every output gets the same residues.
* **Handling non-standard residues**: register residue names from MCPB.py or similar tools as amino acids with `--modified-residue`.

`all` runs `extract` for you when you pass `-c`. To choose how large the ML region should be, see [Building the ML region and layers](model-setup.md).

---

## Examples

### 1. Select by residue ID with a total ligand charge

Give the substrate as chain:name:number and its total charge as one number.

```bash
mlmm extract -i complex.pdb -c 'A:GPP:301' -o pocket.pdb -l -3 --out-json
```

The console prints `[extract] Atoms after truncation: N`; the region has N atoms. The charge is on the line `[extract] Total active site model charge`. `result.json` has the same numbers in `n_atoms_extracted` and `total_charge`. Open `pocket.pdb` in a viewer and check that the residues of the reaction are in the region.

### 2. Substrate given as a PDB file

Pass a PDB file of the substrate as the center, with a charge for each residue name.

```bash
mlmm extract -i complex.pdb -c substrate.pdb -o pocket.pdb -l 'GPP:-3,SAM:1'
```

The substrate file must have the same coordinates as the complex (within 0.001 Å).

### 3. Select by residue name

Name the residues; every residue with that name is a center.

```bash
mlmm extract -i complex.pdb -c 'GPP,SAM' -o pocket.pdb -l 'GPP:-3,SAM:1'
```

### 4. Several structures in one run

List the reactant and product after one `-i`; both get the same residues, written as one multi-MODEL PDB.

```bash
mlmm extract -i complex_R.pdb complex_P.pdb -c 'A:GPP:301,A:SAM:302' \
    -o pocket_multi.pdb -l 'GPP:-3,SAM:1'
```

(extract-modified-residue)=
### 5. Non-standard residues (`--modified-residue`)

Tools such as Amber's MCPB.py give metal-coordinating residues non-standard names (`HD1`, `HE1`, `CM1`, `AP1`). `extract` does not know these names, so it does not cut their main chain, and it prints:

```text
[extract] WARNING: Residue HD1 83 may be an amino acid (has N, CA, C, O) but is not recognized as a standard residue name. Backbone truncation was not applied. Consider preparing the active site model manually.
```

Register the names as amino acids with their charges as `NAME:charge`.

```bash
mlmm extract -i complex.pdb -c 'A:SUB:301' -o pocket.pdb \
    --modified-residue 'HD1:0,HE1:0,CM1:0,AP1:0'
```

---

## How it works

1. **Centers**: `-c` lists the substrate, cofactors, and metals. Each entry is a [residue selector](cli-conventions.md#residue-selectors); the recommended form is `A:TYR:44` (chain:name:number), and `A:SAM`, a name such as `SAM`, a number, or a PDB/mmCIF file of the substrate also work. `--selected-resn` adds residues in the same forms without starting a distance search.
2. **Neighbors**: a residue joins the region when one of its atoms lies within `-r` (default 2.6 Å) of a center atom, or within `--radius-het2het` when both atoms are other than C and H; either cutoff is enough. Waters count unless `--no-include-h2o` is given, and with `--exclude-backbone` contacts through main-chain atoms of amino acids do not count. Three kinds of residues are then added: the disulfide partner of a selected cysteine (S–S ≤ 2.5 Å), the N-side neighbor of a selected proline, and, without `--exclude-backbone`, the two residues peptide-bonded to an amino acid whose main-chain atom touches a center.
3. **Main-chain cuts**: a run of consecutive amino acids keeps its internal main chain and is cut at both ends so that each end stops at CA; a residue that comes in alone keeps only its side chain. Amino acids in `-c` follow the same rule; the residues peptide-bonded to them lie within `-r` and join, so they keep all their atoms. Prolines keep their ring. With `--exclude-backbone`, amino acids lose all main-chain atoms, except between amino acids in `-c` that are peptide-bonded to each other. Waters and non-amino-acid residues are never cut.
4. **Charge**: amino acids and ions take their charges from built-in tables, waters are 0, and other residues are 0 unless `-l` gives them a charge (see below).
5. **Cap hydrogens (only with `--add-linkh`)**: where a cut leaves CA or CB without its bonded partner (CB–CA, CA–N, CA–C; only CA–C for proline), a hydrogen is placed 1.09 Å from that carbon along the old bond. The caps are written after a `TER` record as `HETATM` atoms `HL` in residue `LKH`, chain `L`. A `--model-pdb` file must not contain them: the ML/MM calculator adds its own link hydrogens on the `parm7` bonds that cross the ML/MM boundary, so leave `--add-linkh` off unless you want a capped pocket on its own.

To decide the boundary yourself and check it, see {ref}`How to construct a reliable model.pdb <model-pdb-selection>`.

### Charge summary

`-l` takes either a mapping such as `'GPP:-3,SAM:1'` or one number. Unknown residues are those that the appendix does not list as amino acids, ions, or waters. A number is split evenly over the unknown residues in `-c`, or over all unknown residues when `-c` has none; `-l -3` over two residues gives −1.5 each, and the total stays −3. With a mapping, unknown residues that are not listed stay 0. At the default verbosity the console prints the protein, ligand, and ion charges and then `Total active site model charge`. With several inputs, the summary is for the first one.

### Several structures

With several inputs, each structure selects its residues, and the union of the selections is applied to every structure, so all outputs have the same atoms. Each output keeps its own coordinates. The console prints `[extract:multi] Atoms after truncation (model k): N` for each model.

---

## Output files

```text
./
├─ pocket.pdb    # the cut-out region (cap hydrogens after a TER record only with --add-linkh)
├─ pocket.cif    # mmCIF input, or PDB input too large for the PDB columns
├─ result.json   # with --out-json, next to the first output file
└─ summary.json  # copy of result.json; read result.json (with --out-json)
```

| Inputs | `-o` | Output |
| --- | --- | --- |
| One | not given | `pocket.pdb` |
| Several | not given | `pocket_<input name>.pdb` for each input |
| Several | one path | one multi-MODEL PDB |
| Several | one path per input | one PDB per input |

Any other number of `-o` paths stops with an error, and so does an output path that is the input file itself. Missing parent directories are created. `result.json` holds the atom counts (`n_atoms_raw`, `n_atoms_extracted`, `n_link_hydrogens`), the charges (`total_charge`, `protein_charge`, `ligand_total_charge`, `ion_total_charge`), and the settings used; see [JSON Output Reference](json-output.md). mmCIF input, and PDB input too large for the PDB columns, also get `.cif` files that keep the original identifiers (see {ref}`mmCIF input <mmcif-input>`).

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path(s) | (required) | Protein–ligand PDB/mmCIF files. List several after one `-i`, or repeat `-i`; they must have the same atoms in the same order |
| `-c, --center` | text | (required) | Center residues or a PDB/mmCIF file of the substrate (e.g. `'A:TYR:44,A:SAM:301'`) |
| `-o, --output` | path(s) | see [Output files](#output-files) | Output PDB path(s) |
| `-r, --radius` | float | `2.6` | Distance cutoff (Å) around center atoms. `0` adds no neighbors by distance (see [Notes](#notes)) |
| `--radius-het2het` | float | `0` (off) | Second cutoff (Å) between atoms other than C and H |
| `--selected-resn` | text | `""` | Residues to add without a distance search, in the same forms as `-c` |
| `--include-h2o/--no-include-h2o` | flag | `True` | Include waters (HOH, WAT, H2O, DOD, TIP, TIP3, SOL) |
| `--exclude-backbone/--no-exclude-backbone` | flag | `False` | Remove main-chain atoms from amino acids, except between peptide-bonded `-c` residues |
| `--add-linkh/--no-add-linkh` | flag | `False` | Add cap hydrogens where a cut leaves CA or CB without its partner. Leave off for a `--model-pdb` file |
| `--modified-residue` | text | `""` | Residue names to treat as amino acids, as `NAME:charge` (bare `NAME` only for names in the built-in table) |
| `-l, --ligand-charge` | text | `None` | Total charge, or charge per residue name (e.g. `'GPP:-3,SAM:1'`) |
| `--out-json/--no-out-json` | flag | `False` | Write `result.json` and `summary.json` |

See the [generated CLI reference](reference/commands/extract.md) for every option.

---

## Notes

* **`-r 0`** adds no neighbors by distance (the cutoff is evaluated as 0.001 Å): the region is built from the `-c` and `--selected-resn` residues, plus the disulfide partners and the N-side neighbor of a proline that step 2 adds. The same holds for `--radius-het2het 0`.
* **Region size**: check for your system that the result does not change when the ML region grows; a larger `-r` costs more and does not always improve accuracy. See [Make the model larger](model-setup.md#make-the-model-larger).
* **To build a reusable `--model-pdb` by hand**, extract from the PDB that `mm-parm` writes, so that the atoms match the `parm7` ([mm-parm example 4](mm-parm.md#examples)).
* **Names match everywhere**: a name such as `TYR` selects every TYR in every chain, with a warning when there is more than one.
* **`TYR:44` means chain TYR**: with two fields the first is always the chain, and the second is a number (`TYR:44`) or a name (`A:SAM`), so write `A:TYR:44`. In a PDB with an empty chain column, such as the bundled examples, use the name or the number alone.
* **One form per list**: a list that mixes names and numbers, such as `'SAM,44'`, stops with an error.
* **Boundary warnings**: `extract` warns when a bond between nonmetal atoms is cut anywhere other than the CA and CB cuts of step 5. Check the boundary, the charge, and the multiplicity before a calculation; how to choose the boundary is in [Building the ML region and layers](model-setup.md).
* **Same atoms in every input**: inputs with different atom counts or order stop with `[multi] Atom count mismatch` or `[multi] Atom order mismatch`; see [Troubleshooting](troubleshooting.md).
* **Multi-MODEL input** uses only the first MODEL and prints a warning.
* **Alternate locations (altLoc)**: `extract` keeps one altLoc per residue, the one with the highest mean occupancy; use [fix-altloc](fix-altloc.md) when you need the cleaned file itself.
* **`--modified-residue`**: a bare `NAME` is accepted only for names already in the built-in table (see the appendix) and keeps the table's charge, such as −2 for `SEP`; any other bare name stops with an error that asks for `NAME:charge`. `NAME:charge` also overrides a built-in charge for this run (`LYS:0` for a neutral lysine). When `--modified-residue` is not enough, select the ML atoms by hand.
* **Built-in residue names** follow Amber/CHARMM. If a PDB residue shares a name with a different chemical component, give the intended charge with `--modified-residue NAME:charge`.

---

## See also

* [Building the ML region and layers](model-setup.md) — make the ML region smaller or larger, set the MM layers, and freeze atoms
* [all](all.md) — the full workflow; runs `extract` with `-c`
* [mm-parm](mm-parm.md) — build the Amber topology and the matching PDB to extract from
* [define-layer](define-layer.md) — assign the ML, Movable-MM, and Frozen-MM layers from the extracted region
* [fix-altloc](fix-altloc.md) — write a PDB with one alternate location per residue
* [add-elem-info](add-elem-info.md) — fill missing element columns before extraction
* [Common options and selectors](cli-conventions.md) — residue selectors and charge
* [Troubleshooting](troubleshooting.md) — extraction errors
* [Glossary](glossary.md) — ML region, link atom

## Appendix: PDB naming requirements and reference lists

Use this appendix when `extract` classifies a residue or assigns a charge wrongly because of non-standard residue or atom names. With standard PDB names you can skip it.

```{important}
`extract` recognizes amino acids, ions, waters, and main-chain atoms by their PDB residue and atom names. Inputs must follow the standard PDB chemical-component names; other names can misclassify residues or give wrong charges.
```

### Amino acids

Residue names treated as amino acids, with their nominal charges. Only these residues get main-chain cuts and amino-acid charges.

**Standard 20** (charges reflect physiological pH):

- Neutral: `ALA`, `ASN`, `CYS`, `GLN`, `GLY`, `HIS`, `ILE`, `LEU`, `MET`, `PHE`, `PRO`, `SER`, `THR`, `TRP`, `TYR`, `VAL`
- Positive (+1): `ARG`, `LYS`
- Negative (−1): `ASP`, `GLU`

**Canonical extras:** `SEC` (selenocysteine, 0), `PYL` (pyrrolysine, 0).

**Protonation / tautomer variants** (Amber / CHARMM style): `HIP` (+1, fully protonated His), `HID` (0, Nδ-protonated His), `HIE` (0, Nε-protonated His), `ASH` (0, neutral Asp), `GLH` (0, neutral Glu), `LYN` (0, neutral Lys), `ARN` (0, neutral Arg), `TYM` (−1, deprotonated Tyr phenolate).

**Phosphorylated:** dianionic (−2) `SEP`, `TPO`, `PTR`; monoanionic (−1) `S1P`, `T1P`, `Y1P`; phospho-His (phosaa19SB) `H1D` (0), `H2D` (−1), `H1E` (0), `H2E` (−1).

**Cysteine variants:** `CYX` (0, disulfide), `CSO` (0, sulfenic acid), `CSD` (−1, sulfinic acid), `CSX` (0, generic), `OCS` (−1, cysteic acid), `CYM` (−1, deprotonated Cys).

**Lysine variants / carboxylation:** `MLY` (+1), `LLP` (0), `KCX` (−1, Nz-carboxylic acid).

**D-amino acids** (19): `DAL`, `DAR`, `DSG`, `DAS`, `DCY`, `DGN`, `DGL`, `DHI`, `DIL`, `DLE`, `DLY`, `MED`, `DPN`, `DPR`, `DSN`, `DTH`, `DTR`, `DTY`, `DVA`.

**Other modified:** `CGU` (−2, γ-carboxy-glutamate), `CGA` (−1), `PCA` (0, pyroglutamate), `MSE` (0, selenomethionine), `OMT` (0, methionine sulfone), `HYP` (0, hydroxyproline); also `ASA`, `CIR`, `FOR`, `MVA`, `IIL`, `AIB`, `HTN`, `SAR`, `NMC`, `PFF`, `NFA`, `ALY`, `AZF`, `CNX`, `CYF` (all 0).

**N-terminal variants** (`N` prefix): `NALA` (+1), `NARG` (+2), `NASP` (0), `NGLU` (0), `NLYS` (+2), … plus `ACE` (0), `NTER` (+1, generic).
**C-terminal variants** (`C` prefix): `CALA` (−1), `CARG` (0), `CASP` (−2), `CGLU` (−2), `CLYS` (0), … plus `NHE` (0), `NME` (0), `CTER` (−1, generic).

The `N`- and `C`-prefixed Amber names are read as the standard residue (`NALA` → `ALA`); the terminal charge is counted only while the model keeps the N-terminal H1–H3 (H2 and H3 for proline) or OXT.

### Main-chain atoms

Atom names treated as the main chain of an amino acid; under `--exclude-backbone` they are removed, except between peptide-bonded amino acids in `-c`:

```
N, C, O, CA, OXT, H, H1, H2, H3, HN, HA, HA2, HA3
```

### Ions

Ion residue names and their formal charges:

| Charge | Residue names |
|---|---|
| +1 | `LI`, `NA`, `K`, `RB`, `CS`, `TL`, `AG`, `CU1`, `K+`, `NA+`, `NH4`, `H3O+`, `H3O`, `HE+`, `HZ+` |
| +2 | `MG`, `CA`, `SR`, `BA`, `MN`, `FE2`, `CO`, `NI`, `CU`, `ZN`, `CD`, `HG`, `PB`, `BE`, `PD`, `PT`, `SN`, `RA`, `YB2`, `V2+` |
| +3 | `FE`, `AU3`, `AL`, `GA`, `IN`, `CE`, `CR`, `DY`, `EU`, `EU3`, `ER`, `GD3`, `LA`, `LU`, `ND`, `PR`, `SM`, `TB`, `TM`, `Y`, `PU` |
| +4 | `U4+`, `TH`, `HF`, `ZR` |
| −1 | `F`, `CL`, `BR`, `I`, `CL-`, `IOD` |

### Waters

Water residue names (included by default with `--include-h2o`, charge 0):

```
HOH, WAT, H2O, DOD, TIP, TIP3, SOL
```
