# `mlmm all`: Single structure + scan

Give one full-system reactant and the coordinates to drive with `-s`; `all`
runs the scan stages in order, searches the MEP through the stage ends, and,
with `--tsopt`, optimizes each TS candidate and runs IRC. It succeeded when the
console prints `Scientific status: success` under the last
`====== Pipeline summary ======` (and `[Imaginary modes] n=1 (...)` for each
TS).

## When to use

You have only the reactant and can write the chemistry as a sequence of scans
of distances, angles, or dihedrals, for example "first move the methyl from S
of SAM to C7 of GPP, then move H11 of GPP onto OE2 of Glu186" in a
methyltransferase. The start and the stage ends become the inputs of the MEP
search: single-pass `path-opt` by default, or the recursive `path-search` with
`--refine-path`, which can add intermediates it finds. `--mep-mode dmf` selects
DMF in either route.

## Minimal run

```bash
mlmm all --parm7 enzyme.parm7 -i 1.R.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --scan-lists \
        '[("CS1 SAM 320","GPP 321 C7",1.60)]' \
        '[("GPP`321/H11","GLU`186/OE2",0.90)]' \
    --tsopt --thermo \
    -o result_scan
```

Give `-s` once and list every literal after it. Each literal is one stage.
Stages run in order; the final geometry of stage k is the input geometry of
stage k+1. Check the input first with `--dry-run`; it has passed when the
console ends with `[Dry run] --dry-run completed. Input command is valid.`

## Writing --scan-lists

Each literal is a list of target tuples: distance `(i, j, target_Å)`, angle
`(i, j, k, target_deg)`, or dihedral `(i, j, k, l, target_deg)`. `all` takes
target values only; ranges and YAML or JSON spec files belong to the standalone
`scan`. A four-element tuple is an angle target here and a distance range in
`scan`.

An atom is either a 1-based atom number in the full input (`--scan-zero-based`
for 0-based) or a selector in double quotes. A three-field selector gives the
atom name, residue name, and residue number in any order, separated by spaces,
commas, colons, slashes, backticks, or backslashes (`"CS1 SAM 320"`,
`"SAM,320,CS1"`). To name the chain, use the four-field form
`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` in this order (`"A:SAM:320:CS1"`). In a PDB
with an empty chain column (the bundled examples), use three fields; `_` does
not mean an empty chain.

Several tuples in one literal move together in one stage; to drive them one
after another, put them in separate literals. Which to choose for your
reaction: [Staged vs concerted scan](../mlmm-overview/ts-strategy.md#5-staged-vs-concerted-scan).

## Judge success

Read the console, `summary.json`, and the endpoints as in
[all.md](all.md#judge-success). For the scan itself:

- **Stages**: open `_work/scan/stage_NN/scan_trj.xyz` and check that the coordinates change as intended. Each stage prints `[stage 1] Covalent-bond changes (start vs final): Yes` or `No`.
- **Scan record**: `all` runs the scan with `--out-json`, so `_work/scan/result.json` holds the stages; the top-level `summary.json` has no scan stages.
- **Stage ends**: `_work/scan/stage_NN/result.*` are restrained structures, not minima or TS until an unrestrained optimization, or a TS optimization and IRC, confirms them.
- **MEP**: open `mep_trj.pdb` and the TS candidate `_work/path_opt/hei_seg_01.pdb`, and check that `energy_diagram_MEP.png` shows a clear barrier.

```python
import json
d = json.load(open("result_scan/_work/scan/result.json"))
for stage in d["stages"]:
    print(stage["index"], stage["converged"], stage["bond_changes"], stage["target_distances_angstrom"])
```

## Pitfalls and recovery

- **A stage reaches an unexpected geometry.** The restraint was not strong enough, or the stage relaxed into a side product. Inspect the trajectory, tighten the target, or split a complex stage into two simpler ones; do not assume the side product is valid.
- **Python literal error.** Wrap each stage in single quotes outside and double quotes inside; backticks survive bash inside the outer single quotes.
- **Atom not found or matched twice.** Atom names must match those in the input PDB (case is ignored); editing tools such as PyMOL and Maestro sometimes rename `CB` to `CB1`. If a three-field selector matches more than one atom, the run stops; add the chain with `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` or use the atom number.
- **Several `-i` inputs.** `-s` takes exactly one structure; with two or more, the run stops with an error ([all-endpoint-mep.md](all-endpoint-mep.md)).
- **`-s` with `--tsopt`.** This is the scan mode with TS optimization, not TS-only mode.
- **More segments than expected** (`--refine-path` only). Bond-change splitting proposed another candidate intermediate; check it and the neighbouring TS and IRC. The default `path-opt` adds no segments.
- **Walltime.** One stage can take longer than the MEP search; time a pilot stage and budget from it.

## Next step

- [all.md](all.md): mode choice, success criteria, resume, outputs.
- [scan.md](scan.md): `scan`, `scan2d`, and `scan3d` on their own.
- [path.md](path.md): the MEP search after the scans.
- Defaults: `python -c "import mlmm.core.defaults as d; print(d.SEARCH_KW, d.STOPT_KW)"`.
