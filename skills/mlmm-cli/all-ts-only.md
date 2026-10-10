# `mlmm all`: TS-only mode

Give one full-system TS candidate with `--tsopt` and no `-s`; `all` optimizes
the TS, runs IRC, and optimizes both IRC ends (`--thermo` and `--dft` add R/TS/P
frequencies and DFT). It succeeded when the console prints
`[Imaginary modes] n=1 (...)` and `Scientific status: success` under the last
`====== Pipeline summary ======`; then check that the IRC ends are the intended
R and P.

## When to use, and when not

Use it when you already have a TS candidate (from another QM code, an earlier
run such as `result_all/_work/path_opt/hei_seg_01.pdb`, or a manual guess) and
want only the validation stages, without an MEP search.

Without a TS candidate, use [Endpoint mode](all-endpoint-mep.md)
or [Scan-list mode](all-scan-list.md), or `path-search`
([path.md](path.md)). If you suspect the candidate does not sit between the
right reactant and product, find the connectivity with `path-search` first.

## Minimal run

```bash
mlmm all --parm7 enzyme.parm7 -i ts_candidate.xyz --ref-pdb enzyme_layered.pdb \
    -q -1 -m 1 -b uma \
    --tsopt --thermo \
    -o result_ts_only
```

With a PDB that carries the layers in its B-factors, `-l` derives the charge:

```bash
mlmm all --parm7 enzyme.parm7 -i ts_candidate.pdb \
    -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -o result_ts_only
```

Add `--dft` (and `--func-basis 'wb97m-v/def2-svp'`) for DFT single points of
the ML region on R, TS, and P. With `-c`, the ML region is cut out around the
given residues; without it, it comes from the B-factor layers or `--model-pdb`.

## How the mode is chosen

`all` runs TS-only mode for exactly one `-i` input with `--tsopt` and no `-s`;
there is no flag that forces it, and `summary.log` shows `Pipeline mode` as
`TS-only`. After the preparation, the MEP search is skipped and the run starts
at the TS optimization. One input without `-s` or `--tsopt` stops with
`BadParameter`. One input with both `-s` and `--tsopt` runs Scan-list mode.

IRC starts only when the TS optimization converged, its final Hessian was
computed, and n_imag ≥ 1. With n_imag ≥ 2, IRC runs with a warning as a
diagnostic, not as a first-order TS. `--skip-final-freq` keeps the TS but leaves
n_imag unknown, so the run stops before IRC.

## Judge success

- **TS**: a successful TS optimization gives one imaginary mode along the reaction coordinate. `post_segments[0].tsopt.n_imaginary_modes` should be 1 and `.imaginary_frequencies_cm` gives its wavenumber; `segments/seg_01/ts/result.json` has the same `n_imaginary_modes`. Play `segments/seg_01/ts/vib/imag_*_trj.xyz` to see that the mode moves the bonds that form or break.
- **Stopped before IRC**: `summary.json` has `pipeline_stop` with `stage` `before_irc` and the reason; the TS files stay in `segments/seg_01/ts/` with a copy in `structures/ts.*`. n_imag is computed after a `--stop-plateau` stop but not at the cycle limit.
- **Status**: `scientific_status` is `success` only when every requested stage converged and n_imag = 1; otherwise read `scientific_status_reasons`.
- **Endpoints**: open `segments/seg_01/irc/finished_irc_trj.xyz` and `segments/seg_01/reactant.*` and `product.*`, and read `segments[0].bond_changes`. Even if the IRC does not converge, the result is usable when the endpoint optimizations reach the intended R and P.
- **R and P names**: with no MEP, the higher-energy IRC end is named the reactant (on an exact tie, the left end). The names and the barrier follow this energy order, not a known chemical direction; `post_segments[0].endpoint_assignment` records the rule as `policy` `higher_energy_endpoint_as_reactant` with `chemical_direction_known: false`. The barrier from P is `barrier_kcal − delta_kcal`. Compare both ends with the intended states before reporting a forward barrier.
- **Energies**: `post_segments[0].mlip.barrier_kcal` and `.delta_kcal` (same values in `segments[0]`); `gibbs_mlip` (`--thermo`) and `dft` (`--dft`) carry the same keys. `irc/result.json` gives the raw IRC energies.

```python
import json
d = json.load(open("result_ts_only/summary.json"))
seg, post = d["segments"][0], d["post_segments"][0]
print(d["scientific_status"], d.get("scientific_status_reasons"), d.get("pipeline_stop"))
print(post["tsopt"].get("n_imaginary_modes"), post["tsopt"].get("imaginary_frequencies_cm"))
# the rest exists only when IRC ran
print(seg["barrier_kcal"], seg["delta_kcal"], seg["bond_changes"])
print(post["endpoint_assignment"], post["mlip"]["energies_kcal"])
irc = json.load(open("result_ts_only/segments/seg_01/irc/result.json"))
print(irc["energy_first_hartree"], irc["energy_ts_hartree"], irc["energy_last_hartree"])
```

## Pitfalls and recovery

- **TS optimization not converged.** The last structure is kept and no final Hessian is computed. Read the stop reason in `summary.log` and `segments/seg_01/ts/`, then retry from a better seed or with another optimizer setting (`--opt-mode-post grad` for Dimer).
- **n_imag = 0.** The geometry fell to a minimum; the candidate was not a saddle. The run stops before IRC and is not `success`. Start from a better seed, such as the HEI of an MEP or the top of a scan.
- **n_imag ≥ 2.** The result is `partial`; IRC follows one mode only as a diagnostic. Inspect every mode, check the frozen boundary, then tighten convergence with `--thresh-post gau_tight` or, if the extra mode persists, re-optimize with `--flatten`. A first-order TS needs exactly one imaginary mode along the intended displacement and an IRC that connects the intended states. See [Wrong n_imag](../mlmm-overview/ts-strategy.md#3-wrong-n_imag-after-ts-optimization).
- **`bond_changes` is `(no covalent changes detected)`, or an end is not the intended state.** The TS may connect two nearly identical wells or other minima. Watch the imaginary mode and the IRC before trusting the TS.
- **`--no-tsopt` with one input.** It stops with `BadParameter`; TS-only mode needs `--tsopt`.
- **XYZ candidate.** Give `--ref-pdb` for the topology and the B-factor layers, and `-q` and `-m`, because XYZ carries no charge or multiplicity.
- **R and P labels.** See R and P names above; inspect `reactant.*` and `product.*` to tell which chemical states the IRC reached.
- **TS still not found.** See [When the TS does not come out](../mlmm-overview/ts-strategy.md#7-when-the-ts-does-not-come-out).

## Run the stages yourself

For finer control, check the TS before running the next commands:

```bash
mlmm tsopt -i ts.xyz --parm7 enzyme.parm7 --ref-pdb enzyme_layered.pdb -q -1 -m 1 --out-json -o result_tsopt -b uma
mlmm irc   -i result_tsopt/final_geometry.xyz --parm7 enzyme.parm7 --ref-pdb enzyme_layered.pdb -q -1 -m 1 --out-json -o result_irc -b uma
mlmm opt   -i result_irc/forward_first.xyz --parm7 enzyme.parm7 --ref-pdb enzyme_layered.pdb -q -1 -m 1 --out-json -o result_end_forward -b uma
mlmm opt   -i result_irc/backward_last.xyz --parm7 enzyme.parm7 --ref-pdb enzyme_layered.pdb -q -1 -m 1 --out-json -o result_end_backward -b uma
mlmm freq  -i result_tsopt/final_geometry.xyz --parm7 enzyme.parm7 --ref-pdb enzyme_layered.pdb -q -1 -m 1 --out-json -o result_freq -b uma
```

## Outputs

```text
ts_candidate (PDB, mmCIF, or XYZ with --ref-pdb)
  └─ tsopt (RS-P-RFO by default; Dimer with --opt-mode-post grad)
       └─ converged, final Hessian, n_imag ≥ 1
            └─ irc (forward and backward) → endpoint optimization (RFO by default)
                 ├─ freq (--thermo)
                 └─ dft (--dft)
```

`summary.json` and `summary.log` sit at the top of `--out-dir` with
`ml_region.pdb`, `mm_parm/` (without `--parm7`), and `layered/`. Cite
`segments/seg_01/reactant.*`, `ts.*`, and `product.*` (in the input format).
`seg_01/` also has `ts/` (`final_geometry.*`, `vib/imag_*_trj.xyz`,
`result.json`), `irc/` (`{forward,backward,finished}_irc_trj.xyz`,
`result.json`), `structures/` (the raw IRC ends `reactant_irc.*` and
`product_irc.*`, and `ts.*`), `freq/{R,TS,P}/` with `frequencies_cm-1.txt` and
`thermoanalysis.yaml` (`--thermo`), `dft/{R,TS,P}/result.yaml` (`--dft`), and
the energy diagrams. The IRC-derived files appear only when IRC ran. There are
no MEP files and no `_work/path_opt/`.

## Next step

- [all.md](all.md): mode choice, success criteria, resume.
- [tsopt.md](tsopt.md), [irc.md](irc.md), [freq.md](freq.md), [dft.md](dft.md): each stage on its own.
- [Reading outputs](../mlmm-overview/outputs.md#bond-changes): IRC ends and bond changes.
