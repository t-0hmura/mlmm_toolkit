# BezA endpoint and scan example

This example computes the two-step reaction of the geranyl pyrophosphate (GPP)
C6-methyltransferase BezA on the 9,215-atom full system: methyl transfer from
SAM to GPP, then proton abstraction from GPP by Glu186 (Glu170 in the original
study). The mechanism was reported by Tsutsumi et al.,
*Angew. Chem. Int. Ed.* **2022**, 61, e202111217
([DOI: 10.1002/anie.202111217](https://doi.org/10.1002/anie.202111217)).

## Files

- `1.R.pdb`: reactant.
- `2.IM.pdb`: carbocation intermediate. `run.sh` does not use it; it is there
  for comparison with the result, or as the middle structure of a
  multi-structure `-i` input.
- `3.P.pdb`: product.
- `run.sh`: the two runs below.

## Run

AmberTools is required; a GPU and a job scheduler are strongly recommended for
this full system. Give a new output directory (the script stops if it already
exists):

```bash
bash examples/beza/run.sh /path/to/mlmm_beza_output
```

`run.sh` runs `mlmm all` twice, both with `--tsopt --thermo`:

1. The MEP between `1.R.pdb` and `3.P.pdb`, with `--refine-path`.
2. A two-stage distance scan from `1.R.pdb`: stage 1 moves the methyl group of
   SAM onto C7 of GPP, and stage 2 moves H11 of GPP onto OE2 of Glu186.

`mlmm all` builds the Amber topology and the ML/MM layers from the input
structures.

## Main outputs

The output directory holds `result_mep/` with its console log `mep.log`, and
`result_scan/` with `scan.log`. Each result directory contains:

- `segments/seg_NN/{reactant,ts,product}.pdb`: the optimized R, TS, and P of
  each step;
- `energy_diagram_MEP.png`;
- `summary.log` / `summary.json`.

`--refine-path` splits the path where bonds change, so in `result_mep/` the
methyl transfer and the proton abstraction are reported as separate segments.

## Next steps

- [Quickstart: `mlmm all`](../../docs/quickstart-all.md)
- [Quickstart: `mlmm all --scan-lists`](../../docs/quickstart-scan.md)
- [`all`](../../docs/all.md)
- [Output Directory Layout](../../docs/output-layout.md)
