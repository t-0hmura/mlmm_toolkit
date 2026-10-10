# Common options and selectors

This page collects the conventions shared by every `mlmm` command: flags, ML/MM options, residue and atom selectors, charge and multiplicity, exit codes, and configuration precedence.

## Boolean options

Turn a stage or behavior on or off with paired flags:

| Form | Example |
|---|---|
| Positive flag | `--tsopt` |
| Negative flag | `--no-tsopt` |

```bash
--tsopt --thermo --no-dft
```

## Progressive help

```bash
mlmm <subcmd> --help               # core options
mlmm <subcmd> --help-advanced      # full option set
```

Every subcommand accepts both.

To check a run before starting it, add `--show-config` to print the YAML file given with `--config` and its top-level keys, and `--dry-run` to check the options and inputs without running the calculation:

```bash
mlmm opt -i input.pdb --parm7 real.parm7 -q -1 --config my_settings.yaml --show-config --dry-run
```

A passing check ends with `[Dry run] --dry-run completed. Input command is valid.`

(verbosity-levels)=

## Verbosity levels

`-v/--verbose LEVEL` is an integer from 0 to 3 (**default 2**) that sets how much each command prints to the console. It is a per-command option, so write it with the subcommand, e.g. `mlmm opt -v 1 ...`. The same four levels apply to every command; command pages describe only what their own command adds.

| Level | What you see |
|---|---|
| `-v 0` | Silent. Confirm success from the exit code and the output files. |
| `-v 1` | Milestones only: version, input summary, key settings, output location, dry-run / final status. No banner, `[command]`, `[mode]`, or configuration printout. |
| `-v 2` | Default. Adds the banner, `[command]`, `[mode]`, stage progress, the main optimizer cycle table, terminal status, the one-line Hessian summary, thermo / DFT summaries, and elapsed time. |
| `-v 3` | Debug: the full configuration in effect, backend DEBUG, raw optimizer and internal-coordinate output, `[HessianTiming]`, and `[HessianVRAM]`. |

The level changes only what is printed, not the exit code; judge a run by its {ref}`exit code <exit-codes>`.

(mlmm-options)=

## ML/MM options

`-q` is the charge of the ML region, not of the whole system. The individual ML/MM commands need the full-system topology (`--parm7`) and the ML region; `all` builds both.

Each command takes the ML region from the first of these that is given: `--model-pdb`, then `--model-indices`, then the ML atoms (B-factor 0) of the input PDB with `--detect-layer`. Even with an explicit ML region, `--detect-layer` still reads the Movable-MM and Frozen-MM layers from the B-factors.

| Option | Meaning | Default |
|---|---|---|
| `--parm7` | Amber parm7 topology of the whole enzyme complex (the full system). | Required (`all` builds it) |
| `--model-pdb` | PDB of the ML region only, without link hydrogens; atom names and order match the full-system PDB and parm7. When given, it defines the ML region. | None |
| `--model-indices` | 1-based atom indices of the ML region, comma-separated; ranges such as `1-5` are allowed. Used when `--model-pdb` is not given. | None |
| `--detect-layer/--no-detect-layer` | Reads B-factors 0 / 10 / 20 as the ML, Movable-MM, and Frozen-MM {ref}`layers <mm-layers>`. With an explicit ML region, only the MM layers are read. | `--detect-layer` |
| `--ref-pdb` | PDB that gives the atom order and residue information for an XYZ input. | None |
| `--movable-cutoff` | Distance (Å) from the ML region: MM atoms within it move, and atoms beyond it are frozen. | None (layers from the B-factors or `--freeze-atoms`) |
| `--mm-backend` | MM engine: `hessian_ff` or `openmm`. MM Hessians use finite differences by default. | `hessian_ff` |
| `--link-atom-method` | Placement of link atoms: `scaled` (g-factor) or `fixed` (1.09 / 1.01 Å). | `scaled` |
| `--cmap/--no-cmap` | Keeps the CMAP terms of the parm7, when present, in both the real-system and model-system MM calculations. | `--cmap` |

`mlmm all` builds the topology when `--parm7` is omitted, and takes the ML region from the extracted model with `-c`.

```bash
mlmm path-search -i R.pdb P.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1
```

The choice of MLIP backend, precision, and `--uma-workers` is described in [MLIP Backends](backends.md).

## Residue selectors

`-c/--center` (on `extract` and `all`) names the residues at the center of the model.

| Form | Example | What it selects |
|---|---|---|
| Chain + name + number (recommended) | `-c 'A:TYR:44'` / `-c 'A:TYR:44,A:SAM:123'` | Exactly one residue per entry. |
| Chain + name | `-c 'A:SAM'` | Every SAM in chain A; a warning is logged when more than one matches. |
| Chain + number | `-c 'A:123'` / `-c 'A:123,B:456'` / `-c 'A:123A'` | Residue 123 of chain A; a trailing letter is the insertion code. |
| Name only | `-c 'SAM,GPP'` / `-c 'LIG'` | Every residue with that name in any chain; a warning is logged when more than one matches. |
| Number only | `-c '123,456'` / `-c '123A'` | The residue with that number in every chain. |
| Structure file | `-c substrate.pdb` / `-c substrate.cif` | The residues whose coordinates match a separate PDB / mmCIF file. |

Long mmCIF chain IDs and residue numbers above 9999 use the same forms. Chain IDs are case-sensitive; residue names are not. A PDB with an empty chain column, such as the bundled example PDBs, takes only the name or number forms (`-c 'SAM,GPP,MG'`, `--selected-resn '44,63,186'`).

```bash
mlmm extract -i complex.cif -c 'LONG_CHAIN:SAM' -o pocket.pdb        # every SAM in chain LONG_CHAIN
mlmm extract -i complex.cif -c 'LONG_CHAIN:SAM:10001' -o pocket.pdb  # one SAM
mlmm extract -i complex.cif -c 'LONG_CHAIN:10001' -o pocket.pdb      # chain + number
```

(selected-resn-takes-ids)=
### `--selected-resn` uses the same residue selectors

`--selected-resn` on `extract` and `all` force-includes residues in the model and accepts the same forms. For example, `A:TYR:44` includes one residue, `A:SAM` every SAM in chain A, and `A:123A` one insertion-code residue. A name without a chain, such as `TYR`, includes every match and warns when more than one is present.

(charge-specification)=

## Charge specification

`-q/--charge` is the charge of the ML region. For PDB/mmCIF inputs, `--ligand-charge/-l` lets you specify charges only for non-standard residues (substrates, cofactors); the ML-region charge is then **automatically derived** by summing standard amino-acid charges, ions, and your ligand charges. Ions in the built-in table take their charge from it (`MG` is +2); listing one in `-l` with the same value is accepted, and a different value is ignored with a warning.

```bash
-l 'SAM:1,GPP:-3'        # per-residue mapping (recommended)
-l 'LIG:-2'              # single mapping
-l -3                    # single integer = total ligand charge
-q 0                     # explicit ML-region charge
```

**Resolution order** (highest priority first):

1. Explicit `-q/--charge`.
2. With `--ligand-charge/-l`: the total of the standard residues, ions, and your ligand charges in the ML region.
3. `calc.model_charge` from `--config`.
4. Otherwise, stop with an error.

The console prints the derived charge as `Total active site model charge`; [Check the model](model-setup.md#check-the-model) lists the lines to read after `extract`.

```{tip}
Always provide `--ligand-charge/-l` for non-standard residues (substrates, cofactors, unusual ligands) to ensure correct charge propagation.
```

## Spin multiplicity

```bash
-m 1    # singlet (default)
-m 2    # doublet
-m 3    # triplet
```

`-m/--multiplicity` is the multiplicity of the ML region; without it, `calc.model_mult` from `--config` is used, and then 1. Use the same value in `all` and per-stage subcommands. `mm-parm` takes the separate `--ligand-mult` for the per-residue multiplicities used in parameterization.

## Atom selectors

Atom selectors name single atoms in `--scan-lists` and in `--distance-restraint` of `opt`; `--freeze-atoms` takes only 1-based atom numbers (see {ref}`Freeze atoms and restrain distances <freeze-atoms-and-restraints>`).

```bash
--scan-lists '[(1, 5, 2.0)]'                                          # 1-based integer indices
--scan-lists '[("SAM,320,CS1", "GPP,321,C7", 1.60)]'                  # residue name, number, atom name
--scan-lists '[("A:SAM:320:CS1", "A:GPP:321:C7", 1.60)]'              # with chain ID
```

A three-field selector gives the residue name, residue number, and atom name in any order, separated by spaces, commas, colons, slashes, backticks, or backslashes (`"SAM,320,CS1"`, `"SAM 320 CS1"`, and `"320,SAM,CS1"` select the same atom). Three fields never include a chain; to name the chain, use the four-field form `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` in this order, with any insertion code after the number (`A:SAM:12B:C1`).

(scan-list-spec)=

### Scan-list spec

`--scan-lists/-s` (on `scan`, `scan2d`, `scan3d`, and `all`) accepts one or more inline Python literals. The standalone `scan` / `scan2d` / `scan3d` commands additionally accept a YAML / JSON spec file path; use a file for complex multi-stage runs, inline literals for short cases.

**YAML / JSON spec file** (root = mapping; key is `stages` for `scan`, `pairs` for `scan2d` / `scan3d`):

```yaml
one_based: true            # optional; defaults to the command's --one-based/--zero-based (1-based)
stages:                    # scan
  - [[1, 5, 1.35]]
  - [[1, 5, 2.20], [2, 8, 1.80]]
```

```yaml
one_based: true
pairs:                     # scan2d (exactly 2 entries) / scan3d (exactly 3 entries)
  - [1, 5, 1.30, 3.10]
  - [2, 8, 1.20, 3.20]
```

In a YAML spec, a `scan` stage may mix distance targets with distance, angle, or dihedral ranges. Each `scan2d` / `scan3d` axis is `(i,j,low,high)`, `(i,j,k,low,high)`, or `(i,j,k,l,low,high)`. Indices may be integers, three-field selectors, or positional `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` selectors. Distances are in Å, angles and dihedrals in degrees.

**Inline literal**: wrap in **single quotes** so the shell does not interpret parens / spaces; use double-quoted PDB selectors inside.

```bash
-s '[(atom1, atom2, target_Å), ...]'             # scan: triples
-s '[(atom1, atom2, low_Å, high_Å), ...]'        # distance range
-s '[(atom1, atom2, atom3, low_deg, high_deg)]'  # angle range
-s '[("SAM,320,CS1","GPP,321,C7",1.60)]'         # quoted selectors
-s "[(\"SAM,320,CS1\",\"GPP,321,C7\",1.60)]"       # avoid: double-quoted outer literal requires escaping inner quotes
```

For `scan`, one literal = one **stage**; multiple stages → multiple literals after a single `--scan-lists` flag. For `scan2d` / `scan3d`, only one literal is accepted (no multi-stage support). A two-stage `scan` that first forms one bond and then moves a proton:

```bash
mlmm scan -i r.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s '[("SAM,320,CS1","GPP,321,C7",1.60)]' '[("GPP,321,H11","GLU,186,OE2",0.90)]'
```

| Command | Accepted scan specification |
| --- | --- |
| `scan` | Inline distance targets `(i,j,target)`, or ranges `(i,j,low,high)`, `(i,j,k,low,high)`, and `(i,j,k,l,low,high)` scanned in both directions from the input geometry ([Bidirectional scan](scan.md#bidirectional-scan-4-tuple)); YAML/JSON is also accepted |
| `all --scan-lists` | Inline targets only: distance `(i,j,target)`, angle `(i,j,k,deg)`, or dihedral `(i,j,k,l,deg)` (no ranges, no YAML/JSON) |
| `scan2d` | One literal/file containing exactly two distance, angle, or dihedral axes |
| `scan3d` | One literal/file containing exactly three distance, angle, or dihedral axes |

A four-element tuple is therefore a distance range in `scan` and an angle target in `all`.

## Input file requirements

- **PDB** — must contain hydrogens (add via `reduce` / `pdb2pqr` / Open Babel / `mlmm mm-parm --add-h`) and element symbols in cols 77–78 (`mlmm add-elem-info` if missing). Multiple PDBs must share identical atoms in the same order.
- **mmCIF** — accepted by the 14 commands listed in {ref}`mmCIF and large structures <mmcif-input>` below.
- **XYZ** — the calculation commands, including `all`, accept XYZ together with `--ref-pdb`, which gives the atom order and residue information.
- **Amber parm7 (`--parm7`)** — force-field parameters of the full system; its atoms must be the same, in the same order, as in every coordinate input.

(mmcif-input)=
### mmCIF and large structures

`all`, `extract`, `define-layer`, `sp`, `opt`, `tsopt`, `freq`, `irc`, `dft`, `scan`, `scan2d`, `scan3d`, `path-opt`, and `path-search` accept `.cif` and `.mmcif`, and so does `--ref-pdb`. The standalone `mm-parm` reads PDB only; `all` converts an mmCIF input before it builds the parameters. Use mmCIF for chain IDs longer than one character, residue numbers beyond four digits, atom serial numbers beyond five digits, or structures with 10,000 or more residues.

`mlmm` reads the first coordinate model and keeps one alternate location (altLoc) per residue, the one with the highest mean occupancy. During the calculation the atoms carry temporary chain IDs and residue numbers; output CIF files restore the original chain IDs, residue numbers, and insertion codes. Large or non-standard PDB files (10,000 or more residues, 99,999 or more atoms, hybrid-36 numbering, or numbers that overflow their columns, for example) are handled the same way.

```bash
mlmm all -i reactant.cif product.cif -c 'enzyme_A:SAM:10001,enzyme_A:GPP:10002' \
    -l 'SAM:1,GPP:-3' --tsopt --thermo -o result
mlmm tsopt -i hei.xyz --ref-pdb full_system.mmcif --parm7 full_system.parm7 -q -2 -o result_tsopt
```

The second command takes the coordinates from the XYZ file and the topology from the mmCIF file. Residue and atom selectors use the original chain IDs and residue numbers:

| Context | Form | Example |
|---|---|---|
| `extract` / `all -c` by number | `CHAIN:RESSEQ[ICODE]` | `enzyme_A:10001B` |
| `extract` / `all -c` by name and number | `CHAIN:RESNAME:RESSEQ[ICODE]` | `enzyme_A:SAM:10001B` |
| Scan atom | `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` | `enzyme_A:SAM:10001B:CS1` |

For the `.cif` files written next to each output, see [Output Directory Layout](output-layout.md).

(trajectory-one-frame)=
### Extract one frame from a trajectory

A `_trj.xyz` file is a plain multi-frame XYZ file (each frame is an atom-count line, a comment line, and the atom lines), so frame k (counted from 1) can be extracted with:

```bash
N=$(head -1 scan_trj.xyz); k=12
sed -n "$(( (k-1)*(N+2)+1 )),$(( k*(N+2) ))p" scan_trj.xyz > frame_12.xyz
```

To continue with the PDB topology, pass the original PDB to `--ref-pdb` of the next calculation command; the coordinates come from the frame.

(exit-codes)=

## Exit codes

| Code | Meaning |
|---|---|
| `0` | Success or usable partial results |
| `1` | Non-convergence, no usable result, runtime exception, or output failure |
| `2` | Invalid input, CLI arguments, or configuration |
| `130` | User interruption (SIGINT) |

Exit codes do not depend on JSON output. Exit code `0` covers both `success` and `partial`; tell them apart by `scientific_status`. `all` and `path-search` write it to `summary.log` (`all` also prints `Scientific status:` on the console), and the other commands record it in `result.json` when run with `--out-json` (see [Execution and requested-stage completion](json-output.md#execution-and-requested-stage-completion)). Without `--out-json`, read the console lines that the command's page lists, for example [Judging the IRC](irc.md#judging-the-irc) for `irc` (intrinsic reaction coordinate).

(opt-mode-semantics)=

## `--opt-mode` (subcommand-dependent)

`--opt-mode` picks the optimizer. L-BFGS (limited-memory BFGS) and RFO (rational function optimization) find minima; Dimer, RS-P-RFO (restricted-step partitioned RFO), RS-I-RFO (restricted-step image RFO), and TRIM (trust-region image minimization) search for a TS (transition state).

| Subcommand | `grad` alias selects | `hess` alias selects | Default |
|---|---|---|---|
| `opt` | L-BFGS (`lbfgs`) | RFO (`rfo`) | `grad` (L-BFGS) |
| `tsopt` | Dimer (`dimer`) | RS-P-RFO (`rsprfo`) | `hess` (RS-P-RFO) |
| `path-opt` (endpoint preopt) | L-BFGS | RFO | `grad` |
| `path-search` (single-structure optimizations) | L-BFGS | RFO | `grad` |
| `scan` / `scan2d` / `scan3d` (per-grid relaxation) | L-BFGS | RFO | `grad` |
| `all` (pre-optimization, scan, and MEP search, `--opt-mode`) | L-BFGS | RFO | `grad` |
| `all` (TS optimization, `--opt-mode-post`) | Dimer | RS-P-RFO | `hess` |
| `all` (endpoint optimization after IRC, `--opt-mode-post`) | L-BFGS | RFO | `hess` |

When `--opt-mode-post` is omitted and `--opt-mode` is given explicitly, `all` uses that value for the TS and endpoint optimizations as well.

The same `--opt-mode` value selects a **different algorithm** on each subcommand, and the defaults differ, so check the table before copying a recipe. Algorithm names are accepted on `opt` (`lbfgs` / `rfo`) and `tsopt` (`dimer` / `rsirfo` / `trim` / `rsprfo`); all other subcommands accept only `grad` / `hess`. On `tsopt`, `--opt-mode grad` is therefore a **Dimer** TS search, not an L-BFGS minimization, and this Dimer periodically computes the Hessian to update its direction. Write `--opt-mode dimer` or `rsirfo` on `tsopt` and `--opt-mode lbfgs` or `rfo` on `opt` to make a recipe unambiguous.

`--microiter/--no-microiter` (on by default) alternates one Hessian-based step of the ML region with an L-BFGS relaxation of the MM atoms: in `opt` with RFO, and in `tsopt` with RS-P-RFO, RS-I-RFO, and TRIM. It is not used with Dimer or with electrostatic embedding (`--embedcharge`).

## CLI ↔ YAML name mismatches

A few CLI flags use slightly different names than their YAML counterparts, and a few are renamed when wrapped in `all`. The main flags and their YAML keys are in {ref}`YAML Reference › Common CLI-to-YAML mapping <common-cli-to-yaml-mapping>`. The most-asked cases:

(pressure-vs-pressure-atm)=
- **`--pressure` (CLI) vs `pressure_atm` (YAML)** — on `freq` the flag is `--pressure FLOAT`; in `all` it is exposed as `--freq-pressure`. YAML key: `thermo.pressure_atm`. Both carry **atm** values (default 1.0).

- **`--step-size` (CLI) vs `step_length` (YAML)** — on `irc` the flag is `--step-size FLOAT` (bohr, default 0.10); in `all` it is `--irc-step-size`. YAML key: `irc.step_length`.

(engine-vs-dft-engine)=
- **`--dft-engine` (CLI) vs `engine` (YAML)** — `--dft-engine` (alias `--engine`) selects `gpu` (GPU4PySCF, the default) or `cpu` (PySCF). YAML key: `dft.engine` for `dft` and the DFT stage of `all`, and `calc.dft.engine` for the other calculation commands run with `--backend dft`.

```bash
mlmm irc -i ts.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 --step-size 0.05
mlmm all -i r.pdb p.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt --irc-step-size 0.05
```

## YAML configuration

```bash
mlmm all -i r.pdb p.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --config my_settings.yaml --out-dir result/
```

(configuration-precedence)=

```
built-in defaults  <  --config (YAML)  <  CLI options
```

`mlmm <subcmd> --help-advanced` and the [Command Reference](reference/commands/index.md) show the built-in default of each option (`[default: …]`). Only *explicitly supplied* CLI values override YAML; options left at their CLI default do not mask YAML values. This order holds for every command that takes `--config`. Full schema: [YAML Reference](yaml-reference.md).

## Output directory

`-o/--out-dir ./my_results/` sets the output directory of a calculation command; each command has its own default, listed in [Output Directory Layout](output-layout.md). The preparation commands take file paths instead: `extract` takes one or more `-o/--output` paths, `mm-parm` an `-o/--out-prefix`, and `define-layer` an `-o/--output` path.

## Notes

* **Put the chain before a residue name and number.** `TYR:44` is read as chain `TYR`, residue 44, and stops with a "not found" error; write `A:TYR:44`.
* **Use one kind of residue selector per list.** The forms in the table fall into three kinds: names only (`SAM`), chain + name with or without a number (`A:SAM`, `A:TYR:44`), and numbers with or without a chain (`A:123`, `123`). Kinds cannot be mixed in one list: `A:TYR:44,A:SAM` works, while `A:SAM,SAM`, `A:44,A:SAM`, and `SAM,TYR:44` stop with an error.
* **Atom selectors on PDB files with an empty chain column** take three fields, such as `'SER:11:HG'` or `'SER 11 HG'`. `_` does not stand for an empty chain, so `'_:SER:11:HG'` matches no atom and stops with an error.
* **`--model-indices` and the derived charge**: the charge cannot be derived from `--ligand-charge` with `--model-indices`; give `-q`, or define the ML region with `--model-pdb` or the B-factor layers (`--detect-layer`).
* **`--movable-cutoff` turns off `--detect-layer`**: with `--movable-cutoff`, the MM layers come from the distance cutoff, and the B-factor layers of the input are not read.
* **B-factor layers need both ML and MM atoms**: the B-factors are read as layers only when at least one atom is ML (0), at least one is MM (10 or 20), and at least 80% of the atoms carry one of these values (±1.0); an all-zero PDB is not a layer assignment.
* **Inline scan literals**: one inline literal, and all literals of one run, hold either distance targets or ranges. To combine them, list them under `stages:` in a YAML/JSON spec.
* **mmCIF and large structures**: up to 619,938 residues can be handled (62 one-character chain IDs × 9,999 residue numbers in the internal PDB used during the calculation). `fix-altloc` and `add-elem-info` read PDB only; for mmCIF, the altLoc is chosen and element symbols are taken from `_atom_site.type_symbol` when the file is read. A row without `_atom_site.type_symbol`, a non-finite coordinate, or a different atom count between the reaction states stops with an error. Coordinates must fit the fixed-width PDB columns; move a structure that lies far from the origin closer to it.

## See Also

- [Installation](installation.md) — setup and dependencies
- [Getting Started](getting-started.md) — the shortest run and which page to read next
- [Output Directory Layout](output-layout.md) — file names and default output directories
- [Troubleshooting](troubleshooting.md) — common errors and fixes
- [YAML Reference](yaml-reference.md) — all configuration options
- [MLIP Backends](backends.md) — backend choice, precision, and workers
- [ML/MM Calculator](mlmm-calc.md) — how the ML/MM energy, forces, and Hessian are computed
