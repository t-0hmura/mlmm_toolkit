# Quickstart: `mlmm all --scan-lists`

`mlmm all --scan-lists` (`-s`) builds a reaction path from one structure. It drives the distances you choose to target values under harmonic restraints (a scan) on the full ML/MM system and uses the scan endpoints to search the minimum energy path (MEP). With `--tsopt`, it continues to transition-state (TS) optimization and an intrinsic reaction coordinate (IRC) calculation. As in [Quickstart: `mlmm all`](quickstart-all.md), the same run cuts out the ML region and builds the Amber topology and the layers first. You define the reaction coordinate yourself, so the rest of the system changes less and the intended TS is easier to reach. The starting structure can be a reactant, an intermediate, or a product.

The commands below use `1.R.pdb` from the bundled [`examples/beza/`](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples/beza); run them in that directory.

---

## What it is for

* **No product structure**: drive the forming and breaking bonds to make the product (P) from the reactant (R).
* **Two-step reactions**: drive one step after another, such as a methyl transfer followed by a proton transfer.
* **MEP and TS from the scan**: go on from the scanned path to the MEP and the TS in the same run.

## Choosing a scan command

| Goal | Command | Topology and ML region |
| --- | --- | --- |
| Restrained structures and a scan trajectory only | `mlmm scan` | Give `--parm7` and the ML region (`--model-pdb`, `--model-indices`, or the B-factor layers of the input) |
| Go on from the scan to the MEP, and optionally to TS and IRC | `mlmm all -s ...` | Built by `all` |
| A 2D or 3D energy map over two or three coordinates | `mlmm scan2d` / `scan3d` | Give `--parm7` and the ML region |

## Examples

### 1. Check the input first

`--dry-run` checks the input, the layers, the ML-region extraction, the `-s` syntax, the charge and spin, the atom count and order of the topology, and the AmberTools commands, and then stops without any calculation.

```bash
mlmm all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
  -s '[(4360, 4419, 1.50)]' --dry-run
```

The check has passed when the console shows `[all] Dry-run validation passed: ...` and `[all] Planned stages: extract -> mm_parm -> scan -> path_opt.`, and ends with `[Dry run] --dry-run completed. Input command is valid.` When the check fails, the run stops with an error message instead; look up the message in [Troubleshooting](troubleshooting.md).

### 2. One stage

Drive the distance between SAM CS1 (atom 4360) and GPP C7 (atom 4419, C6 in IUPAC numbering) to 1.50 Å. Give the atoms by number or by name; the two commands are the same scan.

```bash
# Atom numbers (1-based by default)
mlmm all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 -s '[(4360, 4419, 1.50)]' -o ./result_scan

# Atom names
mlmm all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 -s '[("SAM,320,CS1", "GPP,321,C7", 1.50)]' -o ./result_scan
```

The run succeeded when the `====== Pipeline summary ======` block near the end of the console shows `Scientific status: success`; `summary.json` holds the same value in `scientific_status`.

### 3. Two stages

Each bracketed list after `-s` is one stage, called a literal (a list written out as text), and the stages run in order. Stage 1 brings the methyl carbon of SAM (CS1) to C7 of GPP (1.50 Å) and away from SD of SAM (3.30 Å); stage 2 moves H11 of GPP away from C7 (2.90 Å) onto OE2 of Glu186 (1.00 Å).

```bash
mlmm all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -s \
  '[("SAM,320,CS1","GPP,321,C7",1.50),("SAM,320,CS1","SAM,320,SD",3.30)]' \
  '[("GPP,321,C7","GPP,321,H11",2.90),("GLU,186,OE2","GPP,321,H11",1.00)]' \
  -o ./result_scan
```

Each stage holds both the bond that forms and the bond that breaks; for how to split your own reaction into stages, see {ref}`Decide how to split the reaction <mechanism-split>`.

## Writing `--scan-lists`

* **Tuples**: each literal is a list of tuples: distance `(atom1, atom2, target_Å)`, angle `(atom1, atom2, atom3, target_deg)`, or dihedral `(atom1, atom2, atom3, atom4, target_deg)`. Distances are in ångströms, and angles and dihedrals in degrees.
* **Atoms**: give integer atom numbers in the order of the full input structure, counting from 1 (from 0 with `--scan-zero-based`), or atom names in double quotes. The scan runs on the full system, so the numbers are those of the input structure, also with `-c`.
* **Stages**: one literal is one stage, and the tuples in the same literal move together. Several literals run as stages in order, each from the result of the previous one. Give `-s` once and list the literals after it.
* **Quoting**: wrap each literal in single quotes so that the shell leaves the parentheses and spaces alone.
* **Bundled PDBs**: they have an empty chain column, so use the name or number forms such as `"SAM,320,CS1"`, not the forms with a chain ID.

For the full syntax of atom names and quoting, see {ref}`Scan-list spec <scan-list-spec>`.

## Output files

The commands above write:

```text
result_scan/
├── summary.log
├── summary.json                 # Results, with scientific_status
├── mep_trj.pdb                  # MEP over all segments
├── energy_diagram_MEP.png       # MEP energy profile
├── ml_region.pdb                # ML region (reusable with --model-pdb)
├── mm_parm/                     # Amber topology 1.R.parm7 and 1.R.rst7 (reusable with --parm7)
├── layered/                     # 1.R_layered.pdb, with the three layers in the B-factors
└── _work/                       # Intermediate files, including the HEI (TS candidate); kept after the run
    ├── scan/
    │   ├── preopt/              # Optimized starting structure
    │   ├── stage_01/            # Scan stage 1
    │   │   ├── result.{xyz,pdb} # Restrained endpoint (optimized without restraints only with --scan-endopt)
    │   │   ├── scan_trj.xyz     # Scan trajectory
    │   │   └── scan.pdb
    │   ├── stage_02/            # Scan stage 2 (two-stage run)
    │   └── result.json          # Scan results of each stage
    └── path_opt/                # MEP search (path_search/ with --refine-path, the recursive MEP search)
        └── hei_seg_01.{xyz,pdb} # Highest-energy image of segment 1
```

These commands stop after the MEP search and do not create `segments/`. Add `--tsopt` for the R/TS/P structures and IRC of each reactive segment in `segments/seg_NN/`, and `--thermo` for `freq/`.

## Checking the result

1. **Completion**: `scientific_status` is `success` when every requested stage converged; otherwise it is `partial` or `failed`, with the [reasons](json-output.md#execution-and-requested-stage-completion) in `scientific_status_reasons`. With `--tsopt`, two checks are left for you: that the imaginary mode moves the bonds that form or break, and that the endpoints are the intended R and P.
2. **Scan**: open `_work/scan/stage_01/scan_trj.xyz` in a viewer and check that the distances change as intended. At the end of each stage the console prints `[stage 1] Covalent-bond changes (start vs final): Yes` or `No`, and `_work/scan/result.json` records it in `stages[].bond_changes`.
3. **MEP**: open `mep_trj.pdb` and `_work/path_opt/hei_seg_01.pdb`, the highest-energy image (HEI) and TS candidate, and check that `energy_diagram_MEP.png` shows a clear barrier.
4. **TS (with `--tsopt`)**: a successful TS optimization gives one imaginary mode along the reaction coordinate. The console then prints `[microiter] Converged!` and then `[Imaginary modes] n=1 (...)`. Open `segments/seg_01/ts/vib/imag_*_trj.xyz` in a viewer and check that the mode moves the bonds that form or break.
5. **Endpoints (with `--tsopt`)**: open `segments/seg_01/irc/finished_irc_trj.xyz` and the optimized endpoints `segments/seg_01/reactant.pdb` and `product.pdb`, and check that they are the intended R and P. Even if the IRC does not converge, the result is usable when the endpoint optimizations reach the intended R and P.

For how `all` judges each stage, see [Reading the run status](all.md#reading-the-run-status).

## Notes

* **Input**: one full-system PDB or mmCIF structure. Give `-c` to cut out the ML region around the given residues; leave it out to take the ML region from the B-factor layers of the input or from `--model-pdb`.
* **Defaults of `all` and `scan`**: the two commands share the scan engine, but the option names and the defaults for the optimizations before and after the scan differ:

  | Command | Step / restraint | Relaxation limit | Optimization before / after the scan |
  | --- | --- | --- | --- |
  | `mlmm all` | `--scan-max-step-size 0.20` Å; `--scan-restraint-k 300` eV/Å² | `--scan-relax-max-cycles 100000` | before: on (follows `--preopt/--no-preopt`; `--scan-preopt/--no-scan-preopt` sets it separately); after: off (`--scan-endopt/--no-scan-endopt`) |
  | `mlmm scan` | `--max-step-size 0.20` Å; `--restraint-k 300` eV/Å² | `--relax-max-cycles 100000` | before: off (`--preopt/--no-preopt`); after: off (`--endopt/--no-endopt`) |

  The restraint strength can also be set with YAML [`bias.k`](yaml-reference.md#bias).
* **Bond changes and `--refine-path`**: every scan endpoint goes on to the MEP search, whether or not its stage shows a bond change. `--refine-path` alone decides whether the MEP is refined recursively (`path_search/`).
* **Scan endpoints**: a finished scan gives restrained structures. They are not minima or transition states until an optimization without restraints, or a TS optimization and IRC, confirms them.
* **`scan` on its own**: the standalone [`scan`](scan.md) also accepts a YAML or JSON spec file; `all -s` takes inline literals only and stops with an error on a spec file. With `--no-preopt --no-endopt` (the defaults of `scan`), the standalone scan starts from the input structure as given and leaves each stage end restrained.
* **`scan --dry-run`**: it checks the input, the charge and spin, and the parse of `--scan-lists`, prints the planned scan, and stops. The extraction and the topology build of `all` are checked only by `all --dry-run`.

## Next steps

- [Quickstart: TS-only mode](quickstart-tsopt.md): optimize and check a TS candidate, such as the top of a scan
- [Tips for studying reaction mechanisms](mechanism-tips.md): how to split a reaction into stages ("Decide how to split the reaction")
- {ref}`Freeze atoms and restrain distances <freeze-atoms-and-restraints>`: what the restraint strength means
- [`scan`](scan.md): run the scan on its own
- [`all`](all.md): full option reference (also `mlmm all --help-advanced`)
- [Glossary](glossary.md): MEP, TS, IRC, and other terms
- [Troubleshooting](troubleshooting.md): find an error message or symptom and its fix
