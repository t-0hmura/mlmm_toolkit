# xTB point-charge correction

For MLIP/MM commands, `--embedcharge` adds this correction:

```text
ΔE = E_xTB(ML + MM point charges) - E_xTB(ML)
```

Forces and Hessians use the corresponding difference. It is disabled by
default and computationally expensive because each correction evaluates xTB
with and without the MM point charges. Keep the ML region to roughly 200–300
atoms or fewer as a practical guideline, then benchmark the actual system,
point-charge count, and hardware. This is not a hard atom limit.

In `mlmm dft`, the same flag instead adds Amber MM point charges directly to
the PySCF QM Hamiltonian; that path does not use the xTB correction.

## Install

The MLIP/MM correction calls the standalone `xtb` executable:

```bash
conda install -c conda-forge xtb
xtb --version
```

The executable must remain on `PATH` in batch jobs. Configure the correction
under `calc` as listed in `docs/yaml-reference.md`; use `--help-advanced` for
the corresponding CLI options.
