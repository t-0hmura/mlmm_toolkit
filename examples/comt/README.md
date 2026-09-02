# COMT endpoint example

This example uses the 3,420-atom reactant and product structures of catechol
O-methyltransferase (COMT). COMT catalyzes an S<sub>N</sub>2 methyl transfer from
SAM to a catecholate oxygen.

Run the endpoint workflow in a new output directory:

```bash
bash examples/comt/run.sh /path/to/mlmm_comt_output
```

The script uses `mlmm all` to build the Amber topology and ML/MM layers from
the input structures. The 4.0 Å region around CAT, SAM, and Mg is assigned to
the ML potential, while the rest of the enzyme remains in the MM environment;
the Mg<sup>2+</sup> charge is recognized automatically. AmberTools is required;
a backend-compatible GPU environment and scheduled execution are recommended
for this full-system example.
