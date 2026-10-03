---
name: mlmm-cli
description: "Per-subcommand guidance for mlmm-toolkit's 22 CLI subcommands: when to use each, a minimal run, how to judge success, pitfalls and recovery, and the next step. SKILL.md is a one-line input-to-output cheatsheet plus the shared ML/MM conventions; heavy commands have their own file (`all.md`, `all-*.md`, `tsopt.md`, `irc.md`, `freq.md`, and others), and `path.md`, `scan.md`, `oniom.md`, and `utilities.md` group related commands. Full flag lists come from `--help-advanced` and the generated reference. TRIGGER on questions about a specific subcommand or shell invocation. SKIP for install, HPC, output-parsing, or structure-format questions, and for choosing the ML region or layers (mlmm-model-setup)."
---

# mlmm CLI

One line per subcommand: what goes in, what comes out, and which file to read.

## Cheatsheet

| Command | In → out | File |
|---|---|---|
| `all` | full-system structures → R/TS/P of each step, `summary.json`, energy diagrams | [all.md](all.md); modes in [all-endpoint-mep.md](all-endpoint-mep.md), [all-scan-list.md](all-scan-list.md), [all-ts-only.md](all-ts-only.md) |
| `extract` | full PDB/mmCIF + `-c` residues → active-site model `pocket.pdb` | [extract.md](extract.md) |
| `mm-parm` | full PDB → Amber `<prefix>.parm7` and `.rst7` | [mm-parm.md](mm-parm.md) |
| `define-layer` | full PDB + ML region → `<input>_layered.pdb` with B-factor layers 0/10/20 | [define-layer.md](define-layer.md) |
| `opt` | one structure → minimum `final_geometry.{xyz,pdb}` | [opt.md](opt.md) |
| `tsopt` | TS candidate → `final_geometry.*`, imaginary-mode animations `vib/imag_*` | [tsopt.md](tsopt.md) |
| `irc` | TS → `finished_irc_trj.xyz`, branch ends `forward_first.*` and `backward_last.*` | [irc.md](irc.md) |
| `freq` | structure → `frequencies_cm-1.txt`, mode animations, thermochemistry | [freq.md](freq.md) |
| `dft` | structure → DFT single point of the ML region in `result.yaml` | [dft.md](dft.md) |
| `path-opt` | two endpoints → one MEP, HEI `hei.pdb` (TS candidate) | [path.md](path.md) |
| `path-search` | two or more structures → MEP split where bonds change, `hei_seg_NN.*`, `summary.json` | [path.md](path.md) |
| `scan` | one structure + `-s` stages → restrained scan, `stage_NN/result.*` | [scan.md](scan.md) |
| `scan2d` | one structure + two coordinates → `surface.csv`, `scan2d_map.png` | [scan.md](scan.md) |
| `scan3d` | one structure + three coordinates → `surface.csv`, `scan3d_density.html` | [scan.md](scan.md) |
| `oniom-export` | full PDB + parm7 + ML region → Gaussian or ORCA ONIOM input | [oniom.md](oniom.md) |
| `oniom-import` | Gaussian or ORCA ONIOM input → `<prefix>.xyz`, `<prefix>_layered.pdb` | [oniom.md](oniom.md) |
| `sp` | structure → ML/MM energy and `forces.npy` (`hessian.npy` with `--hess`) | [utilities.md](utilities.md) |
| `fix-altloc` | PDB with alternate locations → `<input>_clean.pdb` | [utilities.md](utilities.md) |
| `add-elem-info` | PDB with blank element columns → `<input>_add_elem.pdb` | [utilities.md](utilities.md) |
| `bond-summary` | two or more structures → bond changes on stdout | [utilities.md](utilities.md) |
| `trj2fig` | XYZ trajectory with energies → `energy.png` or CSV | [utilities.md](utilities.md) |
| `energy-diagram` | state energies → `energy_diagram.png` | [utilities.md](utilities.md) |

## Pipeline at a glance

`all` prepares the system (extract → mm-parm → define-layer), runs single-pass `path-opt` by default, and runs recursive `path-search` with `--refine-path`. `--tsopt` adds TS optimization, IRC, and endpoint optimization; `--thermo` and `--dft` add frequencies with thermochemistry and DFT single points. Each stage is also its own subcommand. The stage diagram is in [Pipeline at a glance](../mlmm-overview/SKILL.md#pipeline-at-a-glance).

## Shared ML/MM conventions

`-i` takes the full system, not a cut-out model: PDB, mmCIF, or XYZ with `--ref-pdb` (a PDB with the same atoms). Every structure of one run, and the parm7, has the same atoms in the same order. Gaussian and ORCA ONIOM inputs go through `oniom-import`.

The calculation commands other than `all` need the full-system Amber topology (`--parm7`) and an ML region; `all` builds both. The ML region comes from `--model-pdb`, `--model-indices`, or the B-factor layers read by `--detect-layer` (on by default); which one wins, and the 0/10/20 encoding, are in [ML region and layers](../mlmm-structure-io/SKILL.md#ml-region-and-layers). `--link-atom-method` places the link atoms: `scaled` (g-factor, the default) or `fixed` (1.09/1.01 Å).

`-q` is the charge of the ML region, not of the whole system. `-l 'SAM:1,GPP:-3'` gives the charges of non-standard residues, and the ML-region charge is derived from them. The charge comes from explicit `-q`, then the `-l` derivation, then `calc.model_charge` in the `--config` YAML; otherwise the run stops with an error. `-m` is the ML-region multiplicity, otherwise `calc.model_mult`, otherwise 1.

`-b` selects the ML-region backend: `uma` (default), `orb`, `mace`, `aimnet2`, or `dft`. Without `--precision`, UMA and AIMNet2 run in fp32 and ORB and MACE in fp64; AIMNet2 rejects fp64.

Settings apply in the order built-in defaults < `--config` YAML < explicit CLI options. The calculation commands write to `./result_<subcommand>/` by default (for example `./result_all/`, `./result_path_opt/`); change it with `-o/--out-dir`. `extract`, `mm-parm`, and `define-layer` write into `./` and take output file paths instead.

## Cross-cutting pitfalls

- **Wrong charge**: check it before a long job. `extract` prints `Total active site model charge`, and `all --dry-run` runs the preparation and the charge and electron-parity checks in a temporary directory, prints the plan, and skips the calculations. `scan`, `scan2d`, and `scan3d` do not accept `--show-config`. With `--model-indices`, the charge cannot be derived from `-l`; give `-q`.
- **Default backend**: without `-b`, the run uses `uma`. Spell the backend out for production runs.
- **YAML ignored**: explicit CLI values override `--config`; options left at their CLI default do not mask YAML values.
- **Scan literals**: `-s/--scan-lists` takes Python literals. Quote each with single quotes outside and double quotes inside, and watch space- vs backtick-separated atom specs.
- **Hidden options**: an option missing from `--help` may be listed by `--help-advanced`.
- **Out of memory on a Hessian**: narrow the Hessian region with `--hessian-cutoff` (`opt`, `tsopt`, `freq`, `sp`), and keep the default `FiniteDifference` unless a pilot shows `Analytical` is better.
- **`--uma-workers` above 1 with `--hessian-calc-mode Analytical`**: this stops with an error. Use one worker for an analytical Hessian, or `FiniteDifference` with several workers.
- **`--embedcharge`**: the xTB correction runs xTB with and without the MM point charges at every evaluation; keep the ML region to about 200–300 atoms and benchmark first ([backends.md](../mlmm-install-backends/backends.md)). With `-b dft`, it places the MM point charges in the PySCF Hamiltonian instead.

## Where flags and defaults live

- `mlmm <subcommand> --help-advanced` lists every option with its default; the generated reference is `docs/reference/commands/`.
- `--show-config` prints the YAML given with `--config` and its top-level keys, then continues; `all` and `path-search` print the settings after merging defaults, YAML, and CLI, and `sp` prints its merged config and exits.
- `mlmm.core.defaults` holds the built-in defaults, for example `python -c "import mlmm.core.defaults as d; print(d.IRC_KW)"`.

## See also

- [mlmm-overview](../mlmm-overview/SKILL.md) — pick an `all` mode, or run stage by stage.
- [outputs.md](../mlmm-overview/outputs.md) — `summary.json`, `result.json`, and the output tree.
- [ts-strategy.md](../mlmm-overview/ts-strategy.md) — TS candidates, wrong n_imag, and a TS that does not come out.
- [mlmm-model-setup](../mlmm-model-setup/SKILL.md) — what goes into the ML region and the layers.
- [mlmm-structure-io](../mlmm-structure-io/SKILL.md) — formats, residue and atom selectors, charge and multiplicity.
- [mlmm-install-backends](../mlmm-install-backends/SKILL.md) — install, AmberTools, and backends.
- [mlmm-hpc](../mlmm-hpc/SKILL.md) — job scripts for PBS and SLURM.
