# `mlmm bond-summary`

## Purpose

Detect bond changes between consecutive structures. Reports formed /
broken bonds using the same algorithm `path-search` invokes for
recursive segmentation and that `irc` reports under `bond_changes`.

Use it as a sanity check on R vs P, or to understand how a recursive
`path-search` decided where to split a multi-step mechanism.

## Synopsis

```bash
mlmm bond-summary -i a.pdb -i b.pdb [-i c.pdb ...] \
    [--bond-factor 1.2] [--device cpu] [--one-based|--zero-based]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required (≥2) | Repeat `-i` for each structure (XYZ / PDB / GJF). Atom ordering must be identical across all inputs. |
| `--device` | str | `cpu` | Compute device for distance calculations |
| `--bond-factor` | float | `1.2` | Covalent-radius multiplier for bond cutoff |
| `--one-based / --zero-based` | flag | `--one-based` | Atom-index numbering convention in the report |

(Internal `margin_fraction` of 0.05 further shrinks the threshold; see
`bond_changes.py`.)

## Examples

```bash
# R vs P
mlmm bond-summary -i 1.R.pdb -i 3.P.pdb

# Multi-frame check across an MEP
mlmm bond-summary -i frame_01.xyz -i frame_05.xyz -i frame_10.xyz
```

## Output

Text on stdout, e.g.:

```
Pair 1 -> 2:
  Bond formed (2):
    CS1 SAM 320 — C7 GPP 321 :  3.17 Å -> 1.68 Å
    OE2 GLU 186 — H11 GPP 321 :  2.23 Å -> 0.98 Å
  Bond broken (2):
    S SAM 320 — CS1 SAM 320 :  1.80 Å -> 3.43 Å
    C7 GPP 321 — H11 GPP 321 :  1.10 Å -> 2.32 Å
```

## Caveats

- Atom ordering must match across all inputs. If they don't, run
  `mlmm extract` first to canonicalize.
- The default `1.2` × covalent-radius cutoff (with internal margin
  fraction `0.05`) is geometry-only — it does not classify covalent
  vs ionic vs hydrogen-bonded; metal–ligand interactions may hover
  near the cutoff. Raise `--bond-factor` to 1.5 / 1.6 for permissive
  detection.
- Bond-change blocks are also embedded in `summary.json`'s per-segment
  output as a list of `{"Bond formed (k)": [...], "Bond broken (k)": [...]}`
  dicts (one per consecutive pair). See `../mlmm-workflows-output/SKILL.md`.

## See also

- `irc.md` — same algorithm applied to IRC endpoints.
- `path-search.md` — uses bond changes for recursive segmentation.
- `../mlmm-workflows-output/SKILL.md` — interpreting
  `summary.json["segments"][i]["bond_changes"]`.
- Defaults: `import mlmm.defaults as d; print(d.BOND_KW)`