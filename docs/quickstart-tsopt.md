# Quickstart: `mlmm all --tsopt` (TS-only mode)

TS-only mode checks one transition-state (TS) candidate without a minimum energy path (MEP) search. `mlmm all --tsopt` optimizes the TS on the full ML/MM system, follows the intrinsic reaction coordinate (IRC) in both directions, and optimizes the two endpoints, the reactant (R) and the product (P). `--thermo` adds vibrational analysis and thermochemistry, and `--dft` adds DFT single points of the ML region on R, TS, and P. If you already have a TS candidate, give it to this mode to start the TS search directly.

---

## What it is for

* **Refining a candidate from a scan or an MEP**: optimize the top of a scan or the highest-energy image (HEI) of an MEP into a TS.
* **Checking a candidate made another way**: confirm that a structure from another program, or one built by hand, is a TS (n_imag = 1) that connects the intended R and P.
* **Checking an MLIP TS before DFT**: confirm a TS from the machine-learning interatomic potential (MLIP) before you refine it with [DFT/MM](dft-backend.md).

## Minimal command

Pass one full-system TS candidate with `--tsopt`. The bundled example has no TS candidate, so the command below uses the HEI from the [`all` quickstart](quickstart-all.md) run. For your own reaction, pass your own candidate.

```bash
mlmm all -i result_all/_work/path_opt/hei_seg_01.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo -o ./result_ts_only
```

The run succeeded when the `====== Pipeline summary ======` block near the end of the console shows `Scientific status: success`, and the TS optimization prints `[Imaginary modes] n=1 (...)` when the TS has one imaginary mode. `summary.json` holds the status in `scientific_status`.

To compute the TS on exactly the system of an earlier run, pass its topology with `--parm7` and its ML region with `--model-pdb`, and leave out `-c`:

```bash
mlmm all -i result_all/_work/path_opt/hei_seg_01.pdb \
    --parm7 result_all/mm_parm/1.R.parm7 --model-pdb result_all/ml_region.pdb \
    -l 'SAM:1,GPP:-3' --tsopt --thermo -o ./result_ts_only
```

### (Optional) Add DFT single-points

`--dft` adds DFT single points of the ML region on R, TS, and P, and `--func-basis` sets the functional and basis.

```bash
mlmm all -i result_all/_work/path_opt/hei_seg_01.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft --func-basis 'wb97m-v/def2-tzvpd' \
    -o ./result_ts_only
```

To optimize the TS itself with DFT/MM (`-b dft`), and for the DFT extra and GPU memory, see [Refine an MLIP TS with DFT](dft-backend.md).

## Before you run

* **Input**: one full-system TS candidate as PDB or mmCIF, or as XYZ with `--ref-pdb` (a PDB with the same atoms). With `-c`, the ML region is cut out around the given residues. Without `-c`, the ML region comes from `--model-pdb` or from the B-factor layers of the input; the HEI PDB of `all` carries those layers.
* **Charge and multiplicity**: `-q` is the charge of the ML region, not of the whole system. Without `-q`, `all` derives it from the residues in the ML region and `-l`; see {ref}`Charge specification <charge-specification>`. The multiplicity comes from `-m`, then YAML `calc.model_mult`, then 1.
* **When TS-only mode runs**: one input, `--tsopt`, and no `--scan-lists`. Two or more inputs run the MEP search, and one input with `--scan-lists` runs a scan.

## Expected output

A successful run writes:

```text
result_ts_only/
├── summary.log                     # Run summary
├── summary.json                    # Results, with scientific_status
├── ml_region.pdb                   # ML region (reusable with --model-pdb)
├── mm_parm/                        # Amber topology built from the candidate (not with --parm7)
├── layered/                        # hei_seg_01_layered.pdb, the candidate with the three layers in the B-factors (with -c)
└── segments/
    └── seg_01/
        ├── reactant.pdb            # R/TS/P structures, in the format of the input
        ├── ts.pdb
        ├── product.pdb
        ├── energy_diagram_MLIP.png # R–TS–P ML/MM energy diagram (energy_diagram_G_MLIP.png with --thermo)
        ├── ts/
        │   ├── final_geometry.{xyz,pdb}
        │   └── vib/imag_*_trj.xyz  # Animation of each imaginary mode
        ├── irc/
        │   └── {forward,backward,finished}_irc_trj.xyz
        ├── freq/{R,TS,P}/          # --thermo
        │   ├── frequencies_cm-1.txt
        │   └── thermoanalysis.yaml
        └── dft/{R,TS,P}/           # --dft
            └── result.yaml
```

## Checking the result

1. **Completion**: `scientific_status` is `success` when every requested stage converged; otherwise it is `partial` or `failed`, with the [reasons](json-output.md#execution-and-requested-stage-completion) in `scientific_status_reasons`. Two checks are left for you: that the imaginary mode moves the bonds that form or break, and that the endpoints are the intended R and P.
2. **TS mode**: a successful TS optimization gives one imaginary mode along the reaction coordinate. The console then prints `[microiter] Converged!` and then the imaginary wavenumber in cm⁻¹, for example `[Imaginary modes] n=1 ([-593.1])`. Open `segments/seg_01/ts/vib/imag_*_trj.xyz` in a viewer and check that the mode moves the bonds that form or break.
3. **Endpoints**: open `segments/seg_01/irc/finished_irc_trj.xyz` and the R/TS/P structures (`reactant.pdb`, `ts.pdb`, `product.pdb`), and read `segments[0].bond_changes`. The endpoints should be the intended R and P. Even if the IRC does not converge, the result is usable when the endpoint optimizations reach the intended R and P.
4. **Endpoint frequencies**: with `--thermo`, `segments/seg_01/freq/{R,TS,P}/frequencies_cm-1.txt` lists every frequency with its sign. R and P should have no imaginary mode (no value below −5.00 cm⁻¹).
5. **Energies**: `post_segments[0].mlip.barrier_kcal` is ΔE‡ (TS − R) and `.delta_kcal` is ΔE (P − R) in kcal/mol, from the ML/MM energies of the optimized TS and endpoints. Since there is no MEP, `segments[0].barrier_kcal` and `.delta_kcal` hold the same values. With `--thermo`, `post_segments[0].gibbs_mlip.barrier_kcal` and `.delta_kcal` give ΔG‡ and ΔG; with `--dft`, `post_segments[0].dft.barrier_kcal` and `.delta_kcal` give the DFT values.

For how `all` judges each stage, see [Reading the run status](all.md#reading-the-run-status).

| Result | What to try |
|---|---|
| n_imag = 0 | Start from a better candidate, such as the HEI of an MEP or the top of a scan; TS-only mode has no path to guide it. |
| n_imag ≥ 2 | Watch every imaginary mode. Re-optimize with `--flatten`, or tighten convergence with `all --thresh-post gau_tight` (stricter than the default [`baker`](tsopt.md#how-it-works)) or `tsopt --thresh gau_tight`. |
| `bond_changes` is empty, or an endpoint is not the intended one | Check the TS mode and the IRC; the path may connect other minima. |
| R or P keeps an imaginary mode | Check the endpoint geometry and the mode, and tighten the endpoint optimization with `--thresh-post gau_tight`. |

If n_imag is not 1 or the IRC endpoints are not the intended ones, see {ref}`Check the TS <mechanism-check-ts>` and {ref}`When the TS search fails <ts-search-fails>` in [Tips for studying reaction mechanisms](mechanism-tips.md).

## Notes

* **When IRC runs**: `all` goes on to IRC only when the TS optimization converged, its final Hessian finished, and n_imag ≥ 1. With n_imag ≥ 2 the IRC follows one imaginary mode as a diagnostic; it does not make the structure a first-order saddle point.
* **Extra imaginary modes**: `--flatten` displaces the structure along the extra imaginary modes and optimizes again, for up to 50 rounds; see {ref}`When --flatten is on <flatten-precedence-caveat>`.
* **Hessian mode**: keep the default `--hessian-calc-mode FiniteDifference`. Set `--hessian-calc-mode Analytical` only after comparing its speed, memory use, and results with the default on a representative structure of your system.
* **R and P labels**: without an MEP the direction of the reaction is unknown, so TS-only mode labels the higher-energy IRC endpoint R and the lower one P, and records this rule in `endpoint_assignment` in `summary.json`. The labels are not the chemical direction; the barrier from P is `barrier_kcal − delta_kcal`.
* **`tsopt` and `freq` on their own**: for `--opt-mode`, `--max-cycles`, `--no-microiter`, `--hessian-cutoff`, and the other Hessian options, run [`tsopt`](tsopt.md) on its own with the `--parm7` and `--model-pdb` of the run. It writes the final geometry to `final_geometry.{xyz,pdb}` and the imaginary-mode animations to `vib/`. For the full frequency list and thermochemistry of that structure, run [`freq`](freq.md) on `final_geometry.pdb` with the same `--parm7`, `--model-pdb`, `-q`, and `-m`. `mlmm all --help-advanced` lists every option of `all`.

## Next steps

- [Refine an MLIP TS with DFT](dft-backend.md): refine and check the TS with DFT/MM
- [Tips for studying reaction mechanisms](mechanism-tips.md): check the TS, and what to try when the TS search fails
- [`tsopt`](tsopt.md), [`irc`](irc.md), [`freq`](freq.md): run each stage on its own
- [Quickstart: `mlmm all`](quickstart-all.md): build an MEP from R and P
- [Quickstart: scan](quickstart-scan.md): build a path from one structure
- [`all`](all.md), [`dft`](dft.md): full option references
- [Troubleshooting](troubleshooting.md): find an error message or symptom and its fix
