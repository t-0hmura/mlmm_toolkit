# Quickstart: `mlmm all`

## Goal

Run the end-to-end workflow once from two full PDB structures.

## Minimal command

```bash
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --out-dir ./result_all
```

If you want post-processing in the same run:

```bash
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft --out-dir ./result_all
```

## What to check

- `result_all/summary.log`
- `result_all/summary.json`
- `result_all/path_opt/mep.pdb` (or segment outputs under `result_all/path_opt/seg_*/`)

## Tips

- Use `--dry-run` (shown in `--help-advanced`) to validate parsing and execution plan without running heavy stages.
- `mlmm all --help` shows core options; `mlmm all --help-advanced` shows the full list including `--dry-run`.
- To use a different MLIP backend, add `-b orb` (or `mace`, `aimnet2`). Default is `uma`.
- Add `--embedcharge` to enable xTB point-charge embedding for MM-to-ML environmental corrections.

## Next step

- Single-structure scan route: [Quickstart: `mlmm scan` with `-s` (YAML spec)](quickstart-scan-spec.md)
- TS validation route: [Quickstart: `mlmm tsopt` -> `mlmm freq`](quickstart-tsopt-freq.md)
- Full option reference: [all](all.md)
