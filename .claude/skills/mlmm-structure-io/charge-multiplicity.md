# Charge and multiplicity (charge-multiplicity.md)

Sometimes, failure of "the optimizer ran but the chemistry is wrong"
traces back to a wrong total charge or multiplicity. `mlmm-toolkit`
needs both as integers; getting them right is **non-negotiable** for
meaningful energies.

This file is the workflow for figuring them out when they aren't given.

## Multiplicity (`-m`)

| Default | Use case |
|---|---|
| **1 (singlet, closed shell)** | The default for almost every organic / biological / metal-coordination system whose textbook description is closed-shell. |
| 2 (doublet) | Radical species, unpaired-electron transition states (e.g. radical SAM enzymes, Fe(III) low-spin) |
| 3 (triplet) | O₂, some carbenes, ferromagnetic Mn(IV)–Mn(IV) couples in low-symmetry environments |
| 4 (quartet) | Mn(II) / Fe(III) high-spin in some geometries |
| 5 (quintet) | Mn(III), Fe(II) high-spin |
| 6 (sextet) | Mn(II) high-spin, S=5/2 ferric |

> **Default to `-m 1` unless you have a positive reason to disagree.**
> If the system contains a known paramagnetic metal, look up the
> oxidation state and use the high-spin/low-spin assignment from the
> primary literature for that enzyme.

## Charge (`-q`, or summed via `-l 'RES:Q'`)

Two routes:

1. **Direct total** — pass `-q INTEGER` if you already know the total
   charge of the cluster.
2. **Per-residue mapping** — pass `-l 'RES1:Q1,RES2:Q2,...'` and let
   `mlmm-toolkit` sum amino-acid + ligand charges.

The amino-acid table is internal:

```bash
python -c "from mlmm.extract import AMINO_ACIDS as a; print(a)"  # Dict[str, int] of canonical residue charges
```

(Or read `mlmm/extract.py` directly if `dir()` shows other
relevant attributes.)

For non-standard residues / ligands / metals, you must supply `-l`.

## Lookup workflow for an unfamiliar substrate

When you don't know a ligand's formal charge:

### Step 1 — check the primary paper

Most enzyme-mechanism papers state the charge state of the substrate
explicitly in the Methods. The PDB Bank summary page links to the
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
  mlmm extract -i complex.pdb -c '...' -l '...' -o cluster.pdb \
      --show-config --dry-run
  ```

  The output prints the resolved `charge` and `spin`.

- Or run a tiny optimization and read `summary.json`:

  ```bash
  mlmm opt -i cluster.pdb -q ... -m 1 -o /tmp/check
  python -c "import json; print(json.load(open('/tmp/check/result.json'))['charge'])"
  ```

## Permission to web-search

When the agent does not know a charge / multiplicity:

- **Confirm with the user before running a web search.** Many users
  prefer to point to the relevant paper themselves.
- If web search is allowed, prefer authoritative sources in this order:
  primary research paper → PubChem / ChEBI → general databases →
  general web. Cite the source in the agent's output.
- If neither user input nor a clean web source is available, **stop
  and ask** — do not silently default to `-q 0 -m 1` for a metal
  cluster.

## Quick-reference ligand charges (commonly seen)

These are starting points, not gospel — always confirm against the
relevant mechanism.

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
| Mg²⁺ | `MG` | +2 |
| Zn²⁺ | `ZN` | +2 |
| Mn²⁺ / Mn³⁺ | `MN` | +2 / +3 |
| Fe²⁺ / Fe³⁺ | `FE` / `FE2` / `FE3` | +2 / +3 |
| Heme (Fe(III) protoporphyrin) | `HEM` | +1 (with Fe³⁺ + porphyrin²⁻) |
| Phosphate ion (free) | `PO4` | −2 to −3 |

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