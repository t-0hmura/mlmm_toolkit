# `mlmm all` — TS-only mode

## When to use

You already have a **TS candidate** (typically from another QM code, an
older `mlmm-toolkit` run, or a manual guess) and want to run only the
TS validation stages — `tsopt`, then IRC after saddle validation, plus `freq`
with `--thermo` and DFT with `--dft` —
without an MEP search.

## Synopsis

```bash
mlmm all --parm enzyme.parm7 -i ts_candidate.xyz --ref-pdb enzyme_layered.pdb \
    -q -1 -m 1 -b uma \
    --tsopt --thermo \
    [--dft --dft-func-basis 'wb97m-v/def2-svp'] \
    -o result_ts_only
```

Or with a PDB that carries residue / charge info:

```bash
mlmm all --parm enzyme.parm7 -i ts_candidate.pdb \
    -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -o result_ts_only
```

## How it differs from the other two modes

`mlmm all` falls into TS-only mode when:

- exactly **one** `-i` input is given,
- `--tsopt` is enabled,
- **no** `--scan-lists` is provided.

After system preparation, the orchestrator skips the MEP and starts at
`tsopt`. There is **no explicit "force TS-only" flag** — the
mode is selected purely from the input shape. TS-only mode requires
`--tsopt`; passing `--no-tsopt` with a single input raises a
validation error.

IRC starts only after numerical TS convergence, completed terminal PHVA,
and selection of a valid negative root. A converged higher-order result can
continue only as warning-labelled diagnostic IRC and is not a certified
first-order TS. `--skip-final-freq` retains the TS structure but leaves the
reaction direction unverified, so `all` stops before IRC.

For finer control, check the TS result before running the downstream commands:

```bash
mlmm tsopt -i ts.xyz --parm enzyme.parm7 --ref-pdb enzyme_layered.pdb -q -1 -m 1 -o result_tsopt -b uma
mlmm irc   -i result_tsopt/final_geometry.xyz --parm enzyme.parm7 --ref-pdb enzyme_layered.pdb -q -1 -m 1 -o result_irc -b uma
mlmm freq  -i result_tsopt/final_geometry.xyz --parm enzyme.parm7 --ref-pdb enzyme_layered.pdb -q -1 -m 1 -o result_freq -b uma
```

## Pipeline collapses to

```
ts_candidate.{xyz,pdb,cif,mmcif}
       │
       ▼
   [tsopt]            (Dimer or Hessian TS optimizer; default RS-P-RFO)
       │  optimization_status=converged
       │  hessian_status=completed, valid negative root
       │  first_order OR warning-labelled higher_order diagnostic
       ▼
   [irc]              (forward + backward; RFO endpoint refinement by default, via --opt-mode-post hess)
       │
       ▼
   [freq]             (with --thermo)
       │
       ▼
   [dft]              (with --dft)
```

MEP search is skipped; model preparation follows the supplied inputs/options.
The TS child outputs are kept when written;
IRC-derived entries below appear only when the IRC gate passes:

```
result_ts_only/
├── summary.json
├── summary.log
└── segments/
    └── seg_01/
        ├── e1.pdb         chemically unassigned IRC endpoint 1
        ├── ts.pdb         optimized TS
        ├── e2.pdb         chemically unassigned IRC endpoint 2
        ├── ts/            final_geometry.{xyz,pdb}, result.json (requested by all)
        ├── irc/           forward_irc_trj.xyz, backward_irc_trj.xyz, finished_irc_trj.xyz
        ├── freq/          frequencies_cm-1.txt, thermoanalysis.yaml
        ├── structures/    nested copies + raw IRC endpoints ({endpoint_1_irc,ts,endpoint_2_irc}.{xyz,pdb})
        └── (dft/)
```

## Output keys

```python
import json
from pathlib import Path
d = json.load(open("result_ts_only/summary.json"))
seg = d["segments"][0]

# n_imaginary and IRC endpoint energies are not on the summary segment.
# When the TS child reaches its result writer:
ts = json.load(open("result_ts_only/segments/seg_01/ts/result.json"))
print(ts["n_imaginary_modes"])         # should be 1

# The IRC child exists only when IRC was started:
irc_path = Path("result_ts_only/segments/seg_01/irc/result.json")
if irc_path.exists():
    irc = json.load(open(irc_path))
    print(seg["barrier_from_endpoint_1_kcal"])
    print(seg["barrier_from_endpoint_2_kcal"])
    print(seg["bond_changes"])         # bonds broken/formed along the IRC
    print(irc["energy_first_hartree"], irc["energy_ts_hartree"], irc["energy_last_hartree"])
```

The child IRC result reports directional first/last endpoints only. TS-only
mode preserves them as `E1`/`E2`; inspect the structures before attaching
chemical R/P identity.

If `n_imaginary_modes != 1`, the geometry is **not a true first-order
saddle**; see "Distinctive failure modes" below.

`all` preserves the TS child result before deciding whether to continue. It
stops before IRC for numerical non-convergence, zero imaginary modes,
failed/skipped PHVA, or no valid negative root. A numerically converged
higher-order stationary point may continue only as warning-labelled diagnostic
IRC and remains uncertified.

## Distinctive failure modes

| Symptom | Likely cause | Fix |
|---|---|---|
| `tsopt.optimization_status == "not_converged"` | Numerical optimizer did not converge; the terminal structure is retained and PHVA is skipped | Inspect the stop reason, then retry from a better seed or with an appropriate optimizer/coordinate setting |
| `tsopt.n_imaginary_modes == 0` | Geometry collapsed to a minimum during refinement | TS guess was not a real saddle; re-do `path-search` instead |
| `tsopt.n_imaginary_modes >= 2` | Higher-order saddle or unresolved constrained mode; first-order certification failed | Inspect the modes, tighten convergence/frozen-boundary setup, then flatten or reoptimize from a better TS seed. Certification requires exactly one imaginary mode plus the intended displacement and IRC connectivity. |
| `irc.bond_changes == {}` (no bonds change) | TS connects two essentially identical wells (numerical ringing) | Verify the imaginary mode visualization in `freq/`; this is sometimes a non-physical TS |

## When *not* to use TS-only mode

- You do not yet have a TS candidate. Run `path-search` (or the
  full `all` in endpoint-MEP / scan-list mode) instead.
- You have a candidate but suspect the connectivity is wrong (i.e.
  you're not sure whether your "TS" sits between the right reactant
  and product). Use `path-search` to discover the connectivity.

## Caveats

- `--tsopt` is mandatory in TS-only mode; `--no-tsopt` with
  a single PDB triggers a validation error.
- For an XYZ TS candidate, supply `--ref-pdb` for topology and B-factor
  layers, plus `-q` and `-m` because XYZ has no charge or spin metadata.
- Inspect `segments/seg_01/{e1,e2}.pdb` to determine which chemical states the
  IRC reached. IRC direction and endpoint energy do not assign R/P identity.

## See also

- `all.md` — base orientation.
- `tsopt.md`, `irc.md`, `freq.md`, `dft.md` — the underlying
  subcommands (which you can also run standalone if you want
  fine-grained control).
- `mlmm-workflows-output/SKILL.md` — IRC interpretation
  and bond-change conventions.

## ML/MM-aware flags (mlmm-toolkit specific)

In addition to the common flags below,
**`mlmm-toolkit` requires an Amber topology** and supports layer-aware
selection. Most subcommands accept:

| flag | purpose |
|---|---|
| `--parm FILE` | Amber `parm7` topology of the whole enzyme — optional; when omitted, `mm_parm` generates a parm7 from the input PDB |
| `--model-pdb FILE` | Explicit ML-region PDB; takes precedence over extraction- or B-factor-derived ML membership |
| `--detect-layer` | Automatically read valid B-factor MM sublayers; without explicit or extraction-derived ML membership, B-factors also define ML membership. Enabled by default. |
| `--ref-pdb FILE` | Full-enzyme PDB used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `-q, --charge` | Override the net ML-region/model charge (highest priority) |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.

## Mutant-vs-WT barrier comparison (preserve the WT ML region)

Compare barriers formed within each system, then compare those barriers. Do
not subtract mutant and WT absolute energies when their compositions differ.

1. Build and parameterize each complete structure independently.
2. Define chemically corresponding ML and movable regions. Transfer layer
   labels only for atoms with an unambiguous correspondence, and assign every
   added or deleted atom explicitly.
3. Determine charge and multiplicity independently for each system.
4. Use matched backend/method, force field, convergence, and thermochemistry
   settings.
5. Validate each TS with exactly one imaginary mode and inspect the
   displacement. Inspect both IRC endpoints before assigning chemical R/P
   labels.
