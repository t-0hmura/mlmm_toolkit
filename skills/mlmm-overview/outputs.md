# mlmm outputs

How to read what `mlmm all` writes: the output tree, the ML-region files, `summary.json`, the R/TS/P structures, bond changes, failed runs, and the energy diagrams.

## Output tree

```text
result_all/
├─ summary.json                  # Machine-readable results of every stage
├─ summary.log                   # Text summary and the output tree
├─ mep_trj.xyz, mep_trj.pdb      # MEP over all segments
├─ mep_plot.png, energy_diagram_MEP.png
├─ energy_diagram_*_all.png      # R → TS → P diagrams over all segments
├─ irc_plot_all.png              # IRC profiles over all segments (--tsopt)
├─ ml_region.pdb                 # ML region (reusable with --model-pdb)
├─ ml_region_without_linkH.xyz   # ML region without and with link hydrogens
├─ ml_region_with_linkH.xyz      #   (.pdb versions too for PDB input)
├─ mm_parm/                      # parm7 and rst7 (reusable with --parm7)
├─ layered/                      # Full structures with the layers in the B-factors
├─ segments/
│  └─ seg_NN/                    # One step: seg_01, seg_02, ...
│     ├─ reactant.*, ts.*, product.*   # Optimized R, TS, and P (--tsopt)
│     ├─ structures/             # XYZ and PDB of R, TS, P, and the raw IRC ends
│     ├─ energy_diagram_*.png    # R → TS → P diagrams of this step
│     ├─ ts/                     # TS optimization; vib/imag_*_trj.xyz shows each imaginary mode
│     ├─ irc/                    # IRC trajectories and irc_plot.png
│     ├─ endpoint_opt/           # Kept with --dump or when an endpoint did not converge
│     ├─ freq/{R,TS,P}/          # Frequencies and thermochemistry (--thermo)
│     └─ dft/{R,TS,P}/           # DFT single points (--dft)
└─ _work/                        # Intermediate files, including the TS candidates (HEI)
   ├─ pockets/                   # Extracted models (with -c)
   ├─ scan/                      # Staged scan (with -s)
   └─ path_opt/                  # MEP search and hei_seg_NN.* (path_search/ with --refine-path)
```

Report `segments/seg_NN/reactant.*`, `ts.*`, and `product.*`. They follow the input format: `.pdb` for PDB and mmCIF input, with a `.cif` that keeps the original identifiers for mmCIF input and for PDB input too large for the PDB columns, and `.xyz` for XYZ input. TS-only mode has no MEP files and no `_work/path_opt/`; R, TS, and P go to `segments/seg_01/`.

## ML-region files

`--model-pdb` and `--model-indices` select real ML atoms only. Link hydrogens are derived from parm7 bonds crossing that selection, not from a distance cutoff. `all` and `dft` write `ml_region_{without,with}_linkH.xyz`, and PDB input also gets matching `.pdb` files.

Standalone ML/MM commands need `--parm7`; `all` builds it with `mm-parm` when it is omitted. The layers come from `--model-pdb` or `--model-indices`, or from the B-factors 0/10/20 with `--detect-layer`. `dft` computes the DFT energy of the ML region with its link hydrogens only, not of the whole enzyme; the MM energies of the full system and of the ML region come from the parm7 force field.

## summary.json

Status:

- `execution_status` is `completed` or `failed`: whether the requested stages ran. `scientific_status` is `success`, `partial`, or `failed`: whether they met their numerical criteria. `scientific_status_reasons` lists what is missing or unusable and is omitted on success.
- In `all`, a TS with n_imag ≥ 2 gives `partial`, and a TS with n_imag = 0 stops the run before the IRC, so `success` means n_imag = 1 for every TS.
- `expected_item_ids` and `observed_item_ids` list the expected and observed stages; compare them before accepting the run. `stage_outcomes` has one record per stage and `point_outcomes` one per scan point. IRC records say how the integration stopped, with no success or failure verdict per direction.
- `pipeline_stop` appears only when the run stopped early and names the stage and the reason.

Run record:

- `command` is the full invocation and `mlmm_toolkit_version` the version that wrote the file. `pipeline_mode` is `path-opt`, `path-search` (with `--refine-path`), or `tsopt-only`.
- `config` holds the effective settings after the CLI, YAML, and defaults are merged. `mep_mode` names GSM or DMF, and `ts_opt_mode` and `endpoint_opt_mode` the post-processing optimizers. `path_opt_mode` is the single-structure optimizer used to pre-optimize the endpoints, not the MEP algorithm.
- `charge` and `spin` are the ML-region charge and multiplicity. `freeze_atoms` lists the 0-based indices frozen with `--freeze-atoms` or YAML `geom.freeze_atoms`, when there are any. `environment` is `{device, gpu_name, gpu_vram_gb, cuda_version, cpu, n_cpus, ram_gb}`.
- `mlip_backend` names the backend. `mlip_model` is the exact model or checkpoint, `filename:factory` for a custom calculator, and `FUNCTIONAL/BASIS` for `-b dft`. `mlip_precision` is the effective `fp32` or `fp64`, and null for DFT and custom calculators.
- `references` lists the methods actually used by the run, as `{method, citation, doi}` records. The same set is printed at the end of `summary.log` and of the console output, just before the elapsed time.

Segments:

- `n_segments` counts the MEP segments and `n_segments_reactive` those that are not bridges. Check the chemistry before treating a segment as an elementary step.
- Each `segments[]` record has `index`, `tag` (a kink segment, which has no covalent bond change, has `kink` in its tag), `kind` (`seg`, `bridge` for a short connecting path, or `tsopt` in TS-only mode), `converged`, `barrier_kcal`, `delta_kcal`, and `bond_changes`. It holds no structures and no stage records: the structures are files under `segments/seg_NN/`, and the stage results are in `post_segments[]`.
- Each `post_segments[]` record covers one post-processed segment: `tag` and `post_dir`; `tsopt`, with `n_imaginary_modes`, `imaginary_frequencies_cm`, `optimization_status`, and `n_opt_cycles` against `max_cycles`; `irc`, the stop diagnostics of each direction, and `irc_plot` and `irc_traj`; `endpoint_assignment`, how the IRC ends were named R and P; `endpoint_opt`, with `optimization_status`, `n_opt_cycles`, `max_cycles`, and any `stop_reason` for `reactant` and `product`; `thermo_symmetry`, the point group and symmetry number of R, TS, and P when found; and `mep_barrier_kcal` and `mep_delta_kcal`, the MEP values of the segment.
- `mlip`, `gibbs_mlip`, `dft`, and `gibbs_dft_mlip` in `post_segments[]` give R, TS, and P at one level each, with `energies_kcal`, `barrier_kcal`, `delta_kcal`, and a `structures` map keyed R/TS/P: ML/MM energies (`--tsopt`), ML/MM Gibbs energies (`--thermo`), DFT energies of the ML region (`--dft`), and DFT//ML/MM Gibbs energies (`--dft` with `--thermo`).

Barriers:

- n_imag of each TS is `post_segments[].tsopt.n_imaginary_modes`.
- With `--tsopt`, the barrier of a segment is `post_segments[].mlip.barrier_kcal`, the ML/MM energy of the optimized TS minus R; the other levels use the same key under `gibbs_mlip`, `dft`, and `gibbs_dft_mlip`. `segments[].barrier_kcal` is the barrier on the MEP before TS optimization, and in TS-only mode it is TS − R. The barrier from P is `barrier_kcal − delta_kcal`.
- `rate_limiting_step` is the highest local barrier among the reactive segments, as `{segment, barrier_kcal, method}`, plus `mep_barrier_kcal` outside TS-only mode, and null when there is no reactive segment. It takes the highest method available for every segment, in the order `DFT//MLIP/MM_Gibbs`, `DFT`, `MLIP_Gibbs`, `MLIP`, `MEP`; with `-b dft`, `DFT/MM_Gibbs` and `DFT/MM` take the place of `MLIP_Gibbs` and `MLIP`. It is not a microkinetic assignment of the rate-limiting step.
- `overall_reaction_energy_kcal` is R → P over the whole run, and `overall_reaction_energy_method` its method, no higher than that of `rate_limiting_step`.

Files and diagrams:

- `key_output_files` maps each root file to a description, and each `seg_NN` to `{description, files}` with paths relative to that segment directory. `current_output_paths` lists the paths, relative to `--out-dir`, written by this run. Both cover only the current run, so trust them over the existence of a directory.
- `energy_diagrams` lists each diagram with its name, labels, `energies_kcal`, and `image` path.
- `all` passes `--out-json` to `tsopt`, `irc`, and the endpoint optimizations, so `ts/`, `irc/`, and `endpoint_opt/` have a `result.json`. It does not pass it to `freq` or `dft`; rerun those alone with `--out-json` when you need one.

## Oriented R/TS/P paths

```text
segments/seg_NN/
├─ reactant.*, ts.*, product.*   # Endpoint-optimized R and P and the optimized TS; read these
└─ structures/
   ├─ reactant.{xyz,pdb}         # Same optimized R
   ├─ ts.{xyz,pdb}               # Same optimized TS
   ├─ product.{xyz,pdb}          # Same optimized P
   ├─ reactant_irc.{xyz,pdb}     # Raw IRC end named R
   └─ product_irc.{xyz,pdb}      # Raw IRC end named P
```

Read from `segments/seg_NN/` downstream. Use `structures/reactant_irc.*` and `product_irc.*` only to see where the IRC end and the optimized end differ.

In MEP modes the IRC ends are named R and P by matching them to the ends of the MEP segment, by bond topology first and RMSD second; `endpoint_assignment.method` records which one decided.

In TS-only mode there is no MEP, so the higher-energy IRC end is named the reactant and the other the product, with the left end as the reactant on a tie. The names, the file names, `barrier_kcal`, and `delta_kcal` follow this energy order, not a known chemical direction. `endpoint_assignment.policy` is `higher_energy_endpoint_as_reactant` and `chemical_direction_known` is false. Inspect the structures to identify the chemical states; the barrier from P is `barrier_kcal − delta_kcal`.

## Reading keys in Python

```python
import json

d = json.load(open("result_all/summary.json"))
print(d["execution_status"], d["scientific_status"], d.get("scientific_status_reasons"))

# MEP barriers, before TS optimization
for seg in d["segments"]:
    print(f"seg_{seg['index']:02d} ({seg['kind']}): "
          f"MEP barrier = {seg['barrier_kcal']:.1f} kcal/mol, "
          f"ΔE = {seg['delta_kcal']:.1f} kcal/mol")

# n_imag and the barriers after TS optimization
for ps in d.get("post_segments", []):
    n_imag = (ps.get("tsopt") or {}).get("n_imaginary_modes")
    if n_imag != 1:
        print(f"WARNING: {ps['tag']} has n_imag = {n_imag}")
    for level in ("mlip", "gibbs_mlip", "dft", "gibbs_dft_mlip"):
        block = ps.get(level)
        if block:
            print(ps["tag"], level, block["barrier_kcal"], block["delta_kcal"])

# Highest local barrier (a dict, not an int)
rls = d.get("rate_limiting_step")
if rls is not None:
    print(f"highest local barrier: seg_{rls['segment']:02d}, "
          f"{rls['barrier_kcal']:.1f} kcal/mol ({rls['method']})")

# Output files of this run, whatever the shape of each entry
for key, value in (d.get("key_output_files") or {}).items():
    if isinstance(value, str):
        print(key, value)
    elif isinstance(value, dict):
        for rel in value.get("files", []):
            print(f"segments/{key}/{rel}")
```

## Bond changes

`segments[].bond_changes` is a multi-line string. Each atom is its element and 1-based index, the lower index first:

```text
Bond formed (2):
  - C7-C12 : 3.170 Å --> 1.680 Å
  - H38-O45 : 2.410 Å --> 1.020 Å
Bond broken (2):
  - S11-C12 : 1.810 Å --> 3.050 Å
  - C7-H38 : 1.100 Å --> 2.940 Å
```

- An empty list prints as `Bond formed: None` or `Bond broken: None`. A segment with no change gives `(no covalent changes detected)`, a bridge segment gives `""`, and a failed analysis gives `(bond-change analysis unavailable)`.
- In MEP modes the changes are between the first and last images of the MEP segment; in TS-only mode they are between the optimized R and P. In MEP modes, compare the optimized `reactant.*` and `product.*` with the intended R and P yourself.
- A pair counts as bonded within 1.20 times the sum of the covalent radii, less a 5% margin.
- One elementary step usually has 1 to 4 entries in all. More than 8 in one segment suggests that the segmentation failed; inspect the geometries before trusting the barrier.
- The `result.json` of the `irc` command uses another shape, an object `{formed: [...], broken: [...]}` from the first to the last frame; do not confuse the two.

## When a run fails

When `execution_status` is `failed` or `scientific_status` is not `success`, look at:

1. `summary.log`: its header gives both statuses, and an early stop is shown as `Pipeline stop`. On the console, `RESULT WARNING:` lines after `====== Pipeline summary ======` give the reasons.
2. `segments/seg_NN/{ts,irc,endpoint_opt}/result.json`: the status of each stage. `endpoint_opt/failure.json` is written when an endpoint optimization could not run. In an `all` run, `freq/` and `dft/` have no `result.json`.
3. The terminal or scheduler stderr, for tracebacks that are not in the JSON. A run that fails while its options or inputs are checked can stop before writing any JSON, so treat a nonzero exit code as a failure.

Partial outputs are kept. The MEP intermediates are under `_work/path_opt/` or `_work/path_search/`, and `segments/seg_NN/` can hold files of the current run even when a later stage failed. Trust the stage outcomes and `current_output_paths`, not the existence of a directory.

## Energy diagrams

`mlmm all` writes these diagrams at the output root when the energies exist and the image export succeeds:

- `energy_diagram_MEP.png`: MEP energies over all segments, without thermochemistry.
- `energy_diagram_MLIP_all.png`: R → TS → P with ML/MM energies (`--tsopt`).
- `energy_diagram_G_MLIP_all.png`: ML/MM Gibbs energies with QRRHO thermochemistry (`--thermo`).
- `energy_diagram_DFT_all.png`: DFT energies of the ML region on the ML/MM geometries (`--dft`).
- `energy_diagram_G_DFT_plus_MLIP_all.png`: DFT energies plus the ML/MM thermal correction (`--dft` and `--thermo`).

Each `segments/seg_NN/` has the same diagrams for its own step, without `_all`, and `irc/irc_plot.png` shows its IRC. Energies are in kcal/mol relative to the first state. When the PNG cannot be written, the console prints a `NOTE`, and the values stay in `summary.json["energy_diagrams"]`. To draw a diagram from the numbers of several runs, use `mlmm energy-diagram` ([cli/utilities.md](../mlmm-cli/utilities.md)).

## See also

- [cli/all.md](../mlmm-cli/all.md) and the three mode pages: [all-endpoint-mep.md](../mlmm-cli/all-endpoint-mep.md), [all-scan-list.md](../mlmm-cli/all-scan-list.md), [all-ts-only.md](../mlmm-cli/all-ts-only.md).
- [cli/tsopt.md](../mlmm-cli/tsopt.md), [cli/freq.md](../mlmm-cli/freq.md), [cli/irc.md](../mlmm-cli/irc.md), [cli/dft.md](../mlmm-cli/dft.md): the `result.json` of each stage.
- [cli/utilities.md](../mlmm-cli/utilities.md): `bond-summary`, the same bond-change analysis on its own.
- [mlmm-model-setup](../mlmm-model-setup/SKILL.md): the input formats.
