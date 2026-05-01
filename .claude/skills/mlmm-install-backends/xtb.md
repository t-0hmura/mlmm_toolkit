# xTB point-charge embedding correction (xtb.md)

`mlmm-toolkit` can apply an **xTB point-charge embedding correction** on
top of the MLIP energy and forces, replacing the bare ML/MM coupling at
the boundary by an electrostatic embedding obtained from an xTB
single-point with the surrounding MM-region atoms inserted as point
charges. This is exposed through the `--embedcharge` /
`--embedcharge-cutoff` flags on every ML/MM-evaluating subcommand
(`opt`, `tsopt`, `freq`, `irc`, `scan`, `path-opt`, `path-search`,
`all`, …).

## When to use it

| Use it for | Don't use it for |
|---|---|
| Polar active sites where the bare ML/MM electrostatics underestimate the ligand–environment coupling | Reproducing absolute solvation free energies (the correction is geared toward TS / barrier shifts, not bulk thermodynamics) |
| Sanity-checking that an MLIP barrier doesn't change drastically when the surrounding charges are felt with a different embedding | Cluster models without an explicit MM region (use `pdb2reaction` instead) |
| Systems where the MLIP backend was trained on isolated organics and may miss long-range polarization | Systems already well-described by the ML/MM ONIOM coupling alone — `--embedcharge` adds wall-clock cost per energy/force call |

## Install

xTB is shipped as a stand-alone binary (Fortran), not as a PyPI wheel.
Install via conda-forge:

```bash
conda install -c conda-forge xtb
```

Verify:

```bash
which xtb
xtb --version
```

`mlmm-toolkit` calls the `xtb` binary through `subprocess` from
`mlmm/xtb_embedcharge_correction.py`, so the binary must be on `$PATH`
inside the conda env that runs `mlmm`. There is no Python `xtb`
package on PyPI — `pip install xtb` does **not** provide the binary
used by mlmm.

## CLI usage

The flag is a Boolean toggle plus a cutoff radius:

```bash
mlmm tsopt -i ts_guess.pdb -q 0 -m 1 \
    --embedcharge \
    --embedcharge-cutoff 12.0
```

`--embedcharge-cutoff` (default 12.0 Å) controls how far around each
ML-region atom point charges are gathered from the MM region. Atoms
outside the cutoff contribute nothing to the embedding correction;
charges inside are passed through to xTB as the point-charge field.

## Configuration

The corresponding YAML key in the `calc` block is `embedcharge: true`
plus `embedcharge_cutoff: 12.0`. See `mlmm-cli/SKILL.md` and
`mlmm.defaults.MLMM_CALC_KW` for the full set of related keys (e.g.,
`embedcharge_step` for the finite-difference step used inside the
correction).

## Verify the install

```bash
mlmm tsopt --help | grep embedcharge   # should list --embedcharge / --no-embedcharge
xtb --version                          # binary on PATH
```

If the second command fails, install xtb via conda-forge (above). If
the first command does not list the flag, your `mlmm-toolkit` install
is older than the one that introduced `--embedcharge`; upgrade with
`pip install --upgrade mlmm-toolkit`.

## Live source

```bash
python -c "from mlmm import xtb_embedcharge_correction as x; print([n for n in dir(x) if not n.startswith('_')])"
```