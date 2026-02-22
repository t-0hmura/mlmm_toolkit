# Quickstart: `mlmm all`

## Goal

Run the end-to-end workflow once from two full PDB structures.

## Minimal command

```bash
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --out-dir ./result_all
```

If you want post-processing in the same run:

```bash
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
  --tsopt --thermo --dft --out-dir ./result_all
```

## What to check

- `result_all/summary.log`
- `result_all/summary.yaml`
- `result_all/path_search/mep.pdb` (or segment outputs under `result_all/path_search/seg_*/`)

## Tips

- Use `--dry-run` first to validate parsing and execution plan without running heavy stages.
- `mlmm all --help` shows core options; `mlmm all --help-advanced` shows the full list.

## Next step

- Single-structure scan route: [Quickstart: `mlmm scan` with `--spec`](quickstart_scan_spec.md)
- TS validation route: [Quickstart: `mlmm tsopt` -> `mlmm freq`](quickstart_tsopt_freq.md)
- Full option reference: [all](all.md)
