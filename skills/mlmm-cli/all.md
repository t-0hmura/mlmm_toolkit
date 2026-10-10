# `mlmm all`

`all` builds the ML region, the Amber topology, and the layers of the full
system, runs the MEP search or a staged scan, and, when asked, TS optimization,
IRC, frequencies, and DFT in one command. It succeeded when the console prints
`[Imaginary modes] n=1 (...)` for each TS and `Scientific status: success`
under the last `====== Pipeline summary ======`.

## When to use

Use `all` when one job should produce R, TS, P (and IM) structures and barrier
candidates for one or more steps. The MEP stage runs single-pass `path-opt` by
default; `--refine-path` runs the recursive `path-search` instead. `--thermo`
and `--dft` add frequencies with thermochemistry and DFT single points on R,
TS, and P; both need `--tsopt`. To inspect each stage before the next, run the
stages one by one ([overview](../mlmm-overview/SKILL.md#run-stage-by-stage)).

## Pick the mode

| Input | Mode (`Pipeline mode` in `summary.log`) | Page |
|---|---|---|
| Two or more structures in reaction order | Endpoint mode (`MEP`) | [all-endpoint-mep.md](all-endpoint-mep.md) |
| One structure with `-s` | Scan-list mode (`Scan`) | [all-scan-list.md](all-scan-list.md) |
| One structure with `--tsopt` and no `-s` | TS-only mode (`TS-only`) | [all-ts-only.md](all-ts-only.md) |

One structure without `-s` or `--tsopt` stops with `BadParameter` ("Provide at
least two structures with -i/--input in reaction order, or use one structure
with --scan-lists or --tsopt."). `-s` with two or more structures also stops
("--scan-lists requires exactly one input structure"). One structure with both
`-s` and `--tsopt` runs Scan-list mode.

## Minimal run

```bash
mlmm all -i <inputs> [-c <centers>] [-l 'RES:Q,...'] [-s '...'] \
    [--tsopt] [--thermo] [--dft] [-b uma|orb|mace|aimnet2|dft] [-o result_all/]
```

`-i` takes full-system structures. Without `--parm7`, `all` builds the topology
with AmberTools; without `-c`, the ML region comes from the B-factor layers of
the input or from `--model-pdb`. Add `--dry-run` first: it runs the preparation
and the charge and electron-parity checks in a temporary directory, prints the
plan, and skips the calculations. Each mode page has a complete command.

## Judge success

- **Console**: `[Imaginary modes] n=1 (...)` per TS, and `Execution status:` and `Scientific status:` under the last `====== Pipeline summary ======`. When the result is not `success`, `RESULT WARNING:` lines give the reasons.
- **`summary.json`**: `scientific_status` is `success`, `partial`, or `failed`, with `scientific_status_reasons`. `success` means every requested stage converged. A TS with n_imag ≥ 2 gives `partial`; n_imag = 0 stops before IRC and is never `success`.
- **TS**: a successful TS optimization gives one imaginary mode along the reaction coordinate. Read `post_segments[].tsopt.n_imaginary_modes` and `.imaginary_frequencies_cm`, and play `segments/seg_NN/ts/vib/imag_*_trj.xyz` to see that the mode moves the bonds that form or break. n_imag comes from the final Hessian of the ML and movable MM atoms (PHVA); only the rigid motions that leave the frozen atoms in place are projected out, so with a frozen MM layer usually none are.
- **Endpoints**: whether they are the intended R and P is for you to check. Compare `segments/seg_NN/reactant.*` and `product.*`, and the bond changes in section [2] of `summary.log`, with the intended states. Even if the IRC does not converge, the result is usable when the endpoint optimizations reach the intended R and P.
- **Barriers**: `segments[].barrier_kcal` is the barrier on the MEP before TS optimization (TS − R in TS-only mode). After `--tsopt`, read `post_segments[].mlip.barrier_kcal`; `gibbs_mlip` (`--thermo`), `dft` (`--dft`), and `gibbs_dft_mlip` (both) carry the same keys. `rate_limiting_step` is the highest local barrier at the highest method available for every segment, not a microkinetic assignment.

```python
import json
d = json.load(open("result_all/summary.json"))
print(d["execution_status"], d["scientific_status"], d.get("scientific_status_reasons"))
print(d["mlmm_toolkit_version"], d["charge"], d["spin"], d["rate_limiting_step"])
for seg in d["segments"]:
    print(seg["index"], seg["kind"], seg["barrier_kcal"], seg["delta_kcal"], seg["bond_changes"])
for post in d.get("post_segments", []):   # match to segments by "index"
    print(post["index"], post.get("tsopt", {}).get("n_imaginary_modes"))
```

Each `segments` record has `index`, `tag`, `kind`, `barrier_kcal`,
`delta_kcal`, and `bond_changes`, and no stage sub-objects; requested
post-processing is in `post_segments`. Key details:
[summary.json](../mlmm-overview/outputs.md#summaryjson).

## Resume a failed segment

Repeat the original command with the same inputs, extraction, layer, path, and
calculator options and the same `--out-dir`, add `--resume-segment N`, and
change only post-processing options (for example `--tsopt-max-cycles`). The
saved artifacts are checked, earlier segments are kept, and segment N onward,
the summary, and the diagrams are written again. `--resume-segment` needs
`--tsopt`, `--thermo`, or `--dft`, cannot be combined with `--dry-run`, and
stops with an error when the saved inputs, ML region, topology, layered
structures, or MEP do not match the command.

## Pitfalls and recovery

- **`-s` quoting.** Give `-s` once and list every literal after it; each literal is one stage. Each literal is a Python literal: single quotes outside, double quotes inside. Most quoting trouble comes from mixing the two.
- **Cycle limits.** `--max-cycles-gsm` and `--dmf-max-iterations` bound only the selected MEP stage; scan, TS, IRC, freq, and DFT keep their own limits (for example `--tsopt-max-cycles`, `--irc-max-cycles`).
- **Status is not `success`.** Read `scientific_status_reasons`, then the matching block of `summary.log`. The TS and IRC stages (and endpoint and path optimization) write `result.json`; `freq/` and `dft/` do not, so rerun those standalone with `--out-json` when you need one.
- **A segment directory exists but the stage failed.** `segments/seg_NN/` can hold partial files from this run after a later stage failed; check `summary.json` and the stage `result.json`, not the directory. A failed segment can also leave MEP scratch in `_work/path_opt/` (`_work/path_search/` with `--refine-path`).
- **TS stops before IRC.** IRC starts only when the TS optimization converged, its final Hessian was computed, and n_imag ≥ 1. With n_imag ≥ 2, IRC runs with a warning along the mode closest to the MEP direction; it is a diagnostic, not a first-order TS. With n_imag = 0, a cycle limit (no final Hessian), a plateau stop (`--stop-plateau`, the Hessian still gives n_imag), `--skip-final-freq`, or a failed Hessian, the run stops before IRC, keeps the TS files in `segments/seg_NN/ts/`, and does not post-process later segments. Next: [Wrong n_imag](../mlmm-overview/ts-strategy.md#3-wrong-n_imag-after-ts-optimization) and [When the TS does not come out](../mlmm-overview/ts-strategy.md#6-when-the-ts-does-not-come-out).
- **Endpoint optimization trouble.** A non-converged endpoint gives `partial` and keeps `segments/seg_NN/endpoint_opt/`; an endpoint error is written to `endpoint_opt/failure.json`, and that segment skips freq and DFT.
- **`--dft` with `-b dft`.** The run stops at startup. Run `mlmm sp -b dft` as a separate job.
- **No AmberTools.** Without `--parm7`, the run stops when AmberTools is not found; pass a topology built elsewhere with `--parm7`.
- **No ML region.** Without `-c`, `--no-detect-layer` and no `--model-pdb` stop with an error.

## Outputs

Cite `segments/seg_NN/reactant.*`, `ts.*`, and `product.*` (written with
`--tsopt`). The top of `--out-dir` has `summary.json`, `summary.log`,
`mep_trj.xyz` (and `.pdb`), `energy_diagram_MEP.png`,
`energy_diagram_*_all.png`, and the reusable setup: `ml_region.pdb` (for
`--model-pdb`), `mm_parm/` (for `--parm7`), `layered/`, and
`ml_region_{without,with}_linkH.xyz`. Each `seg_NN/` has `ts/`, `irc/`,
`endpoint_opt/`, `freq/{R,TS,P}/` (`--thermo`), and `dft/{R,TS,P}/` (`--dft`).
With `--thermo`, `thermoanalysis.yaml` is kept even under `--no-dump`.
`_work/` holds intermediate files, including the TS candidates (HEI): `pockets/`
(with `-c`), `scan/` (with `-s`), and `path_opt/` (`path_search/` with
`--refine-path`) with `hei_seg_NN.*`; keep it while you use them. Full tree:
[Output tree](../mlmm-overview/outputs.md#output-tree).

## Next step

- Mode pages: [all-endpoint-mep.md](all-endpoint-mep.md), [all-scan-list.md](all-scan-list.md), [all-ts-only.md](all-ts-only.md).
- The stages it runs: [extract.md](extract.md), [mm-parm.md](mm-parm.md), [define-layer.md](define-layer.md), [path.md](path.md), [tsopt.md](tsopt.md), [irc.md](irc.md), [freq.md](freq.md), [dft.md](dft.md).
- [Reading outputs](../mlmm-overview/outputs.md): `summary.json` keys and R/TS/P paths.
- Defaults (`OUT_DIR_ALL` and the per-stage `*_KW`): [Where flags and defaults live](SKILL.md#where-flags-and-defaults-live).
