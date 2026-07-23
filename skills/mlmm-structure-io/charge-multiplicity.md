# Charge and multiplicity (charge-multiplicity.md)

Every run needs a total charge and a multiplicity, and a wrong value can
silently produce a chemically wrong trajectory.

**For PDB/mmCIF input, give the charge with `-l 'RES:Q'` and let
`mlmm-toolkit` derive the ML-region total.** Name only unknown/non-standard
ligand residues; standard amino acids and recognized ions come from internal
tables, and waters and link atoms are neutral. Recheck the reported breakdown
whenever the ML region, residue naming, protonation state, or oxidation state
changes.

Use `-q INTEGER` when no residue metadata is available or when deliberately
overriding the derived ML-region charge. Set multiplicity explicitly whenever
the electronic state is not a verified singlet.

## Multiplicity (`-m`)

| Default | Use case |
|---|---|
| **1 (singlet, closed shell)** | Use only when the modeled electron count and electronic state are known to be closed-shell. Do not infer singlet merely because the structure is biological or metal-bound. |
| 2 (doublet) | Radical species, unpaired-electron transition states (e.g. radical SAM enzymes, Fe(III) low-spin) |
| 3 (triplet) | O₂, some carbenes, Ni(II) (d⁸) high-spin tetrahedral / weak-field octahedral |
| 4 (quartet) | Co²⁺ (d⁷) high-spin, Cr³⁺ / V²⁺ (d³) |
| 5 (quintet) | Mn(III), Fe(II) high-spin |
| 6 (sextet) | Mn(II) high-spin, S=5/2 ferric |

> These are examples, not a spin-state calculator. For metals, radicals,
> antiferromagnetically coupled centers, or uncertain protonation/oxidation
> states, derive charge and multiplicity from the modeled mechanism and primary
> literature.

## Charge (`-q`, or summed via `-l 'RES:Q'`)

1. **Per-residue mapping (for PDB/mmCIF)** — pass
   `-l 'RES1:Q1,RES2:Q2,...'` and let `mlmm-toolkit` sum amino-acid, recognized
   ion, and unknown-ligand charges over the ML region.
2. **Direct total / override** — pass `-q INTEGER` when the input has no
   residue metadata or to deliberately replace the derived charge. In
   `mlmm all`, the explicit charge override wins and the workflow reports the
   value it would otherwise have derived.

The residue and ion tables are internal:

```bash
python -c "from mlmm.core.residue_data import AMINO_ACIDS, ION; print(dict(AMINO_ACIDS)); print(dict(ION))"
```

For unknown/non-standard ligand residues, supply `-l`. Recognized monatomic
ions use the internal `ION` table and must not be repeated in `-l`; a mapping
does not override a recognized ion. To represent a different oxidation state,
use the appropriate residue name in the model or provide the verified total
with `-q`.

## Lookup workflow for an unfamiliar substrate

When you don't know a ligand's formal charge:

### Step 1 — check the primary paper

Most enzyme-mechanism papers state the charge state of the substrate
explicitly in the Methods. The PDB summary page links to the
reference; check there first.

### Step 2 — PubChem / ChEBI

For small-molecule ligands:

- **PubChem** (https://pubchem.ncbi.nlm.nih.gov) — search by ligand
  3-letter code or name. The "Computed Properties" panel lists
  `Formal Charge`.
- **ChEBI** (https://www.ebi.ac.uk/chebi) — biological compound focus,
  often has the charge state used in published mechanisms.
- **PDB Ligand summary** (e.g. `https://www.rcsb.org/ligand/SAM`) —
  shows the canonical SMILES and charge in the deposited model.

### Step 3 — derive from the SMILES

Given a SMILES (e.g. from PubChem), compute the formal charge:

```python
from rdkit import Chem
mol = Chem.MolFromSmiles("CC(=O)[O-]")     # acetate
print(sum(a.GetFormalCharge() for a in mol.GetAtoms()))    # → -1
```

### Step 4 — protonation state at physiological pH

Many ligands have multiple protonation states. Common rule of thumb:

| Group | At pH 7 | Typical formal charge contribution |
|---|---|---|
| Carboxylate (`-COO⁻`) | deprotonated | −1 each |
| Phosphate, monoester | mostly `-OPO₃²⁻` | −2 |
| Phosphate, diester | mostly `-OPO₂⁻` | −1 |
| Triphosphate (e.g. ATP) | fully deprotonated | −4 |
| Sulfonium (e.g. SAM cofactor) | quaternary | +1 |
| Lysine / Arginine side chain | protonated | +1 |
| Aspartate / Glutamate side chain | deprotonated | −1 |
| Histidine | mostly neutral, possibly +1 | 0 or +1 |

Check the literature for the cluster you are modeling — biological
mechanisms sometimes invoke an unusual protonation state.

### Step 5 — sanity-check the total

After summing residue + ligand + metal charges, sanity-check by:

- Letting `mlmm-toolkit` echo the charge it parsed:

  ```bash
  mlmm extract -i complex.pdb -c '...' -l '...' -o cluster.pdb --out-json
  python -c "import json; print(json.load(open('result.json'))['total_charge'])"
  ```

  `--out-json` writes `result.json` next to the output PDB; it records
  `total_charge` plus a per-source breakdown (`protein_charge`,
  `ligand_total_charge`, `ion_total_charge`). The `--verbose` INFO logs
  (on by default) also print the same charge summary.

- Or run a tiny optimization and read `summary.json`:

  ```bash
  mlmm opt -i cluster.pdb --parm real.parm7 -q ... -m 1 -o /tmp/check --out-json
  python -c "import json; print(json.load(open('/tmp/check/result.json'))['charge'])"
  ```

## Source policy for an unknown value

Use authoritative sources in this order: the mechanism's primary paper or
deposited structure documentation, then PubChem/ChEBI/RCSB CCD. Cite the source
and state the modeled protonation and oxidation state. If the sources do not
determine one unambiguous state, ask rather than defaulting a metal or radical
model to `-q 0 -m 1`.

## Quick-reference ligand charges (commonly seen)

Always confirm against the relevant mechanism.

| Ligand | Resname (PDB) | Charge at pH 7 |
|---|---|---|
| Methionine sulfonium (SAM) | `SAM` | +1 |
| Adenosylhomocysteine | `SAH` | 0 |
| Geranyl pyrophosphate | `GPP` | −3 |
| ATP | `ATP` | −4 |
| ADP | `ADP` | −3 |
| GTP | `GTP` | −4 |
| NADH | `NAI`/`NDH` | −2 |
| NAD⁺ | `NAD` | −1 |
| FAD | `FAD` | −2 |
| Pyridoxal phosphate (PLP) | `PLP` | −2 |
| Heme (Fe(III) protoporphyrin) | `HEM` | +1 (with Fe³⁺ + porphyrin²⁻) |
| Phosphate ion (free) | `PO4` | −2 to −3 |

### Monatomic ions are summed from the internal `ION` table

`-l` applies only to residues that are in none of `AMINO_ACIDS`, `ION`, or the
water set. A token such as `-l 'MG:3'` or `-l 'FE:2'` therefore matches no
unknown residue, emits a warning, and is ignored; it does not change the
built-in charge.

| Ion | Resname (PDB) | `ION` value |
|---|---|---|
| Mg²⁺ | `MG` | +2 |
| Zn²⁺ | `ZN` | +2 |
| Mn²⁺ | `MN` | +2 |
| Fe³⁺ | `FE` | +3 |
| Fe²⁺ | `FE2` | +2 |
| Cu²⁺ / Cu⁺ | `CU` / `CU1` | +2 / +1 |
| Na⁺ / K⁺ | `NA` / `K` | +1 |
| Cl⁻ | `CL` | −1 |

When a deposited resname does not represent the intended oxidation state,
correct the model's residue naming or pass the verified ML-region total with
`-q`.

## Multiplicity for metals (look-up shortcuts)

| Metal | Common high-spin S | Multiplicity (2S+1) |
|---|---|---|
| Mn²⁺ (d⁵) | 5/2 | 6 |
| Fe²⁺ (d⁶), high-spin tetrahedral / weak field | 2 | 5 |
| Fe²⁺ (d⁶), low-spin octahedral / strong field | 0 | 1 |
| Fe³⁺ (d⁵) high-spin | 5/2 | 6 |
| Co²⁺ (d⁷) high-spin | 3/2 | 4 |
| Cu²⁺ (d⁹) | 1/2 | 2 |
| Zn²⁺ (d¹⁰) | 0 | 1 |

## See also

- `pdb.md` — `-l 'RES:Q'` syntax and where it parses from.
- `xyz.md` — XYZ has no header, so `-q`/`-m` must be on the CLI.
- `gjf.md` — gjf encodes charge / spin in the header.
- `mlmm-cli/extract.md` — the subcommand that consumes
  `-l` and `-q` first.
