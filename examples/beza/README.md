# BezA endpoint and scan example

This example uses the 9,215-atom full-system structures for geranyl
pyrophosphate (GPP) C6-methyltransferase BezA:

- `1.R.pdb`: reactant;
- `2.IM.pdb`: carbocation intermediate;
- `3.P.pdb`: product.

The reaction has two chemical steps: methyl transfer from SAM to GPP, followed
by proton abstraction from GPP by Glu186 (Glu170 in the original study). The
coordinates accompany the mechanism reported by Tsutsumi et al.,
*Angew. Chem. Int. Ed.* **2022**, 61, e202111217
([DOI: 10.1002/anie.202111217](https://doi.org/10.1002/anie.202111217)).

Run both the endpoint-MEP and staged-scan workflows in a new output directory:

```bash
bash examples/beza/run.sh /path/to/mlmm_beza_output
```

The script uses `mlmm all`, which builds the Amber topology and ML/MM layers
from the input structures. Its endpoint workflow enables `--refine-path`, so
the MEP is recursively divided where bonding changes and the methyl-transfer
and proton-abstraction steps are reported separately. AmberTools is required;
a backend-compatible GPU environment and scheduled execution are strongly
recommended for this full-system example. `2.IM.pdb` is included as the known
intermediate for inspection or an explicitly guided multi-structure run; the
default endpoint workflow searches from `1.R.pdb` to `3.P.pdb` without using
the intermediate as an input.
