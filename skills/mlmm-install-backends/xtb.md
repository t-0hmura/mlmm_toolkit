# Experimental xTB point-charge correction

For MLIP/MM commands, `--embedcharge` adds this experimental correction:

```text
ΔE = E_xTB(ML + MM point charges) - E_xTB(ML)
```

Forces and Hessians use the corresponding difference. It is disabled by
default and computationally expensive. In `mlmm dft`, the same flag instead
adds Amber MM point charges directly to the PySCF QM Hamiltonian.

## Install

The MLIP/MM correction calls the standalone `xtb` executable:

```bash
conda install -c conda-forge xtb
xtb --version
```

The executable must remain on `PATH` in batch jobs. Configure the correction
under `calc` as listed in `docs/yaml-reference.md`; use `--help-advanced` for
the corresponding CLI options.
