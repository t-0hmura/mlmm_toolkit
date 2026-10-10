# `mlmm all`: Endpoint mode

Give two or more full-system structures in reaction order; `all` finds the MEP
between each neighbouring pair and, with `--tsopt`, optimizes each TS candidate
and runs IRC. It succeeded when the console prints `[Imaginary modes] n=1 (...)`
for each TS and `Scientific status: success` under the last
`====== Pipeline summary ======`.

## When to use

You have two or more structures in reaction order (reactant, optional
intermediates, product) with the same atoms in the same order, typically R and
P (sometimes IM) from a published QM or QM/MM study. By default `all` runs one
single-pass `path-opt` per neighbouring pair, so the structures you pass are
taken as the steps. With `--refine-path`, the recursive `path-search` splits a
pair further where bonds change, so you need to give only the steps you already
know.

## Minimal run

```bash
mlmm all --parm7 enzyme.parm7 -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo -o result_mep
```

Add `--dft` (and `--func-basis 'wb97m-v/def2-svp'`) for DFT single points of
the ML region on R, TS, and P. Leave out `--parm7` to build the topology from
the first input with AmberTools. For a known multistep mechanism, give each
intermediate:

```bash
mlmm all --parm7 enzyme.parm7 -i 1.R.pdb 2.IM.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo -o result_mep_3pt
```

## Same atoms in the same order

Every input, and the parm7, needs the same number of atoms, the same element
sequence, and the same residue assignments. `all` builds the ML region from the
first input and applies it to every input. The extraction checks the series and
stops with `[multi] Atom count mismatch between input #1 and input #2: ...` or
`[multi] Atom order mismatch between input #1 and input #2.`, but it does not
map or repair a mismatched series. To check a series before a long job, cut all
inputs in one `extract` run:

```bash
mlmm extract -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -o pocket_R.pdb pocket_P.pdb
```

If they do not match, regenerate every structure with the same protonation
tool and settings (for MD snapshots, the same trajectory and topology), or
start from one structure with `-s` ([all-scan-list.md](all-scan-list.md)).
Rules for states and variants: [model-setup](../mlmm-model-setup/SKILL.md#same-atoms-across-states-and-variants).

## Judge success

Read the console, `summary.json`, and the endpoints as in
[all.md](all.md#judge-success). For this mode also check:

- **MEP**: open `mep_trj.pdb` (or `.xyz`) and `energy_diagram_MEP.png`; the TS candidates are `_work/path_opt/hei_seg_NN.*` (`_work/path_search/` with `--refine-path`).
- **Segments**: `summary.json["segments"]` lists `index`, `kind`, `barrier_kcal`, `delta_kcal`, and `bond_changes` for each segment. The default gives one segment per neighbouring pair. With `--refine-path`, compare `n_segments_reactive` with the number of inputs minus one; `n_segments` also counts bridge segments.
- **R/TS/P**: with `--tsopt`, `segments/seg_NN/{reactant,ts,product}.*` are written after IRC and the endpoint optimizations (RFO by default, L-BFGS with `--opt-mode-post grad`).

## Pitfalls and recovery

- **More bond changes than the inputs imply.** The reaction encoded by the inputs and the path the optimizer found differ. Check which bonds changed with `mlmm bond-summary -i 1.R.pdb 3.P.pdb`, then supply the intermediate yourself or rerun the standalone `path-search` with `--refine-mode minima`. `all` has no `--refine-mode`; with `--refine-path`, set `search.refine_mode` in the `--config` YAML.
- **More reactive segments than input pairs** (`--refine-path` only). This is a candidate decomposition, not proof that the hidden intermediates are real; validate each IM and its TS and IRC.
- **n_imag ≥ 2.** A higher-order saddle or unresolved soft modes, not a validated TS. Compare a Hessian-based optimizer with Dimer (`--opt-mode-post grad`) on the same seed and backend, or try `--flatten`, then rerun the frequency analysis and check the IRC connects the intended states. See [Wrong n_imag](../mlmm-overview/ts-strategy.md#3-wrong-n_imag-after-ts-optimization).
- **Different atoms or order across inputs.** Compare the ordered atom identities, then regenerate the series as above; equal PDB line counts are not enough.
- **GSM or DMF.** GSM is the default. `--mep-mode dmf` selects DMF in either route; set `--dmf-backend cpu` when the GPU implementation runs out of memory. Inspect and validate either MEP.
- **`-s` with several inputs.** It stops with an error; `-s` takes exactly one structure ([all-scan-list.md](all-scan-list.md)).

## Next step

- [all.md](all.md): mode choice, success criteria, resume, outputs.
- [path.md](path.md): what `path-opt` and `path-search` do.
- [bond-summary](utilities.md#bond-summary): what bond-change detection reports.
- [Reading outputs](../mlmm-overview/outputs.md#summaryjson): multi-segment results.
