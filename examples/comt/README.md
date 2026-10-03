# COMT endpoint example

This example computes the S<sub>N</sub>2 methyl transfer from SAM to a
catecholate oxygen in catechol O-methyltransferase (COMT), on the 3,420-atom
full system.

## Files

- `1.R.pdb`: reactant.
- `3.P.pdb`: product.
- `run.sh`: the run below.

## Run

AmberTools is required; a GPU and a job scheduler are recommended for this full
system. Give a new output directory (the script stops if it already exists):

```bash
bash examples/comt/run.sh /path/to/mlmm_comt_output
```

`run.sh` runs `mlmm all` on `1.R.pdb` and `3.P.pdb` with `--tsopt --thermo`.
`mlmm all` builds the Amber topology and the ML/MM layers from the input
structures: the residues within 4.0 Å of CAT, SAM, and Mg form the ML region
(`-c 'CAT,SAM,MG' -r 4.0`), and the rest of the enzyme is MM. The Mg<sup>2+</sup>
charge is recognized automatically; `MG:2` in `-l` only restates it.

## Main outputs

The output directory holds `result/` and its console log `run.log`. `result/`
has the same layout as in the [BezA example](../beza/README.md):
`segments/seg_NN/{reactant,ts,product}.pdb`, `energy_diagram_MEP.png`, and
`summary.log` / `summary.json`.

## Next steps

- [Quickstart: `mlmm all`](../../docs/quickstart-all.md)
- [`all`](../../docs/all.md)
- [Building the ML region and layers](../../docs/model-setup.md): how to choose
  the ML region
