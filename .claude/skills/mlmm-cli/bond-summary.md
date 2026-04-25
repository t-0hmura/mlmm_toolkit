# `mlmm bond-summary`

## Purpose

Detect bond changes between two structures. Reports formed / broken
bonds using the same 1.3× covalent-radius cutoff that `path-search`
uses for recursive segmentation and that `irc` reports as
`bond_changes`.

Use as a sanity check on R vs. P, or to understand how a recursive
`path-search` decided where to split a multi-step mechanism.

## Synopsis

```bash
mlmm bond-summary -i reactant.pdb product.pdb [--threshold 1.3]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path×2 | required | Two structures with identical atom ordering |
| `--threshold` | float | 1.3 | Covalent-radius multiplier for bond cutoff |
| `--csv` | path | none | Optionally write a CSV of pairs + distances |

## Examples

```bash
mlmm bond-summary -i 1.R.pdb 3.P.pdb
```

Output (text on stdout):

```
FORMED:
  C(SAM 320) — C(GPP 321 7)   3.17 Å → 1.68 Å
  O(GLU 186) — H(GPP 321 H11) 2.23 Å → 0.98 Å
BROKEN:
  S(SAM 320) — C(SAM 320 CS1) 1.80 Å → 3.43 Å
  C(GPP 321 7) — H(GPP 321 H11) 1.10 Å → 2.32 Å
```

## Output (when `--csv`)

```
atom1, atom2, dist_R, dist_P, change
"CS1 SAM 320", "C7 GPP 321", 3.170, 1.678, "FORMED"
...
```

## Caveats

- Atom ordering must match between the two inputs. If they don't,
  `mlmm extract` first to canonicalize.
- The 1.3× threshold is empirical. Metal–ligand bonds often hover
  around the cutoff; raise to 1.5 or 1.6 to be more permissive.
- Bond-change detection is geometry-based, not physics-based — it
  doesn't tell you whether the bond is covalent vs. ionic vs.
  hydrogen-bonded.

## See also

- `irc.md` — same algorithm, applied to IRC endpoints.
- `path-search.md` — uses bond changes for recursive segmentation.
- `mlmm-workflows-output/SKILL.md` — interpreting
  `summary.json["segments"][i]["bond_changes"]`.
