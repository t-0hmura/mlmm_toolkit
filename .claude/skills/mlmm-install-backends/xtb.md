# xTB / ALPB solvent layer (xtb.md)

`mlmm` can apply an xTB-ALPB **implicit solvent correction** to
any backend-computed energy. The correction is a separate semi-empirical
SCF call that adds a continuum-solvation term; it is optional and
turned off by default.

## When to use

| Use it when | Skip it when |
|---|---|
| You need a quick, approximate solvent correction for a bare cluster | Fully embedded QM/MM or explicit waters (use `mlmm_toolkit`) |
| Comparing barrier heights between gas-phase and solvated approximation | Reporting absolute solvent free energies (xTB-ALPB is empirical) |
| Substrate is in a public-friendly solvent (water, DMSO, methanol, …) | Solvent isn't in the ALPB parameter set (unusual organics) |

## Install

Two routes; pick whichever your env can use.

**Route A — `xtb-python` (recommended, pip-installable):**

```bash
pip install xtb-python
```

Verify:

```bash
python -c "import xtb; print('xtb-python:', xtb.__version__)"
```

**Route B — system `xtb` binary (xtb 6.7+):**

If your site already has `xtb` as a module or in `$PATH`:

```bash
xtb --version
which xtb
```

`mlmm` calls the binary via subprocess when it can't find the
Python bindings. Make sure `xtb` is on `$PATH` for any PBS / SLURM job
that needs it.

## CLI usage

The solvent correction wraps any other calculator transparently:

```bash
mlmm all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    --solvent water        # ALPB water on top of MLIP
```

Available solvents (pass to `--solvent`): `water`, `methanol`, `ethanol`,
`acetone`, `acetonitrile`, `dmso`, `dmf`, `chloroform`, `dichloromethane`,
`hexane`, `benzene`, `toluene`, `thf`. The exact list comes from xTB's
ALPB parameter set; check `xtb --help` for the most current options.

To turn off: simply omit `--solvent` or pass `--solvent none`.

## How it composes with MLIP / DFT

The corrected energy at each step is:

```
E_total = E_MLIP_or_DFT (in vacuo) + ΔE_ALPB (xTB)
```

The ALPB term is computed from the same atomic positions as the
backbone calculator, so it's geometry-consistent but does not feed back
into the gradient unless you also enable `--solvent-gradient` (off by
default — adds noise on small clusters).

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `xtb-python` import fails on aarch64 | Wheel not yet published for ARM. Use Route B (system binary). |
| `xtb` binary version too old (< 6.7) | ALPB introduced in 6.4 but parameter set updated repeatedly; upgrade if results disagree with documentation. |
| Different barrier vs literature | Could be `--solvent` mismatch or the literature used a different solvation model (CPCM / SMD). State the model in any comparison. |

## See also

- `core.md` — `mlmm` install (xTB layer is plumbed via the
  bundled `backends/xtb_alpb_correction.py` and `backends/solvent.py`).
- `dft.md` — note that xTB-ALPB does **not** stack with PySCF's own
  PCM/COSMO; pick one.
- `mlmm-cli/SKILL.md` — `--solvent` is accepted by `all`,
  `tsopt`, `freq`, `irc`, and `dft`.