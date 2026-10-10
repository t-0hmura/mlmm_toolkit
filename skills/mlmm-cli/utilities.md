# `mlmm` utility commands

Six small commands around the main workflow:

- [sp](#sp): ML/MM single-point energy and forces, optionally the Hessian.
- [fix-altloc](#fix-altloc): keep one alternate location per residue.
- [add-elem-info](#add-elem-info): fill or fix the PDB element columns.
- [bond-summary](#bond-summary): bonds formed and broken between structures.
- [trj2fig](#trj2fig): energy profile of an XYZ trajectory.
- [energy-diagram](#energy-diagram): state diagram from numbers you give.

## sp

**When to use.** The cheapest stage: the ONIOM energy and forces at one
geometry, and with `--hess` the Hessian, without an optimization. Use it to
spot-check a geometry or extract forces. It does not move atoms; relax with
[opt.md](opt.md) and search a TS with [tsopt.md](tsopt.md).

**Minimal run.**

```bash
mlmm sp -i structure.xyz --ref-pdb structure.pdb --parm7 real.parm7 \
    -q 0 -m 1 -o result_sp
mlmm sp -i my.pdb --parm7 real.parm7 -l 'SAM:1' -b uma -o result_sp
mlmm sp -i my.pdb --parm7 real.parm7 -q -1 -m 1 --hess -o result_sp_hess
```

An XYZ needs `--ref-pdb`, the full-system PDB in the same atom order. `-q`
is the ML-region charge, not the whole-system charge; `-l` derives it from
residue charges. The ML region comes from `--model-pdb`, `--model-indices`,
or the B-factors (`--detect-layer`, on by default). `-b` selects `uma`
(default), `orb`, `mace`, `aimnet2`, or `dft`; `--mm-backend` selects
`hessian_ff` (default) or `openmm`. The default output is `./result_sp/`.

**Judge success.** The console prints
`[sp] energy = … a.u.  |force|_max = … a.u./bohr`, and `forces.npy` holds
the ONIOM forces as an `(N, 3)` array over all atoms (Hartree/bohr; zero on
frozen atoms). With `--hess`, `hessian.npy` holds the Hessian of the moving
atoms, without mass weighting. With `--out-json`, `result.json` reports
`stage`, `execution_status`, `scientific_status`, `mlip_backend`,
`mlip_model`, `mlip_precision`, `mm_backend`, `link_atom_method`,
`use_cmap`, `charge`, `spin`, `energy_au`, `forces_path`, and
`hessian_path` (null without `--hess`). For `--calc-file`, `mlip_backend` is
`custom`, `mlip_model` is `filename:factory`, and `mlip_precision` is null.
`summary.json` has the same content; read `result.json`.

**Pitfalls.**

- A failed run prints a one-line `Error: …`, such as
  `ML region electron count inconsistent`, or
  `Unhandled error during single-point:` with a traceback, and exits nonzero.
  A finite energy that looks wrong points to the ML region or its charge and
  multiplicity:
  [Charge and multiplicity](../mlmm-model-setup/SKILL.md#charge-and-multiplicity).
- `--hessian-calc-mode` is `FiniteDifference` (default) or `Analytical`,
  which uses the analytical Hessian of UMA, ORB, MACE, or AIMNet2 for the ML
  region and cannot run with `--uma-workers` above 1. A backend without it
  stops with an error; `sp` does not switch to `FiniteDifference` by itself.
  The MM part uses finite differences unless YAML `calc.mm_fd: false`.
- `sp -b dft` gives the energy and forces only; for atomic charges, use
  [dft.md](dft.md).
- `--config` YAML sets less common options; `--show-config` prints the
  effective settings.

## fix-altloc

**When to use.** Blank the altLoc column (17) of a PDB and keep one label
per residue, when you need the cleaned file itself: input for `mm-parm`,
which does not resolve altLoc, a file for other programs, or a look at the
choice. `extract`, `define-layer`, and the ML/MM calculation commands apply
the same rule on their own when they read a PDB.

**Minimal run.**

```bash
mlmm fix-altloc -i raw.pdb -o cleaned.pdb
mlmm fix-altloc -i raw_pdbs/ -o cleaned_pdbs/
mlmm fix-altloc -i raw.pdb --inplace --force
```

Without `-o`, the output is `<input>_clean.pdb` (file) or `<input>_clean/`
(directory, same relative paths; `--recursive` adds subdirectories). An
`-o` that does not end in `.pdb` is a directory. `--inplace` overwrites the
input, keeps `<name>.pdb.bak`, and ignores `-o`.

**Judge success.** The console prints `[fix-altloc] Fixed altLoc → <file>`,
or `[fix-altloc] Skipped <file> (no altLoc detected).`; a directory prints
`[fix-altloc] Processed N file(s) → …`. In each residue the label whose
atoms have the highest mean occupancy (columns 55–60) is kept; a label
without readable occupancy ranks last, and a tie goes to the label that
appears first. Blank (shared) atoms and the chosen label's atoms remain
with column 17 blanked; atoms found only in another label are dropped, so a
residue never mixes A and B. Every other column, and ANISOU of the kept
atoms, is unchanged.

**Pitfalls.**

- The occupancy rule is a heuristic. If the kept conformer is not the
  chemistry you want, choose it in a structure editor.
- A file without altLoc is skipped and nothing is written; `--force`
  processes it anyway. An existing output stops the run with
  `Output exists: <path> (use --overwrite to overwrite)`; an existing `.bak`
  is never replaced.
- Serial numbers are not renumbered and CONECT is not updated. MODEL blocks
  are processed one by one.
- On raw RCSB files, repair the element columns with `add-elem-info` first,
  then run `fix-altloc`.

## add-elem-info

**When to use.** Fill blank or fix wrong element symbols (columns 77–78),
for example in PDBs from PyMOL, Maestro, or MD. `all` repairs blank fields
itself; run this before a standalone command such as `extract`, which needs
the element column to classify atoms, or when a symbol is wrong rather than
blank. Resolve altLoc separately with [fix-altloc](#fix-altloc).

**Minimal run.**

```bash
mlmm add-elem-info -i raw.pdb -o cleaned.pdb
```

Without `-o`, the output is `<input>_add_elem.pdb`; `--overwrite` without
`-o` replaces the input, and is required when `-o` names the input. By
default only blank or invalid fields are repaired; `--overwrite-elem`
re-infers valid ones too.

**Judge success.** The console prints `[add-elem-info] Wrote: <file>` with
`total atoms`, `assigned/updated`, `kept existing`, and
`assignment breakdown`; with no `WARNING` line, every atom has an element.
Only columns 77–78 of the repaired records change. The element comes from
the fixed-column atom name and the residue name, in this order:

1. Ion residues: polyatomic ions (NH4, H3O+, …) by atom name (H/D → H,
   N → N, O → O); monatomic metals and halogens by residue name, with a
   CL/BR/I/F atom-name fallback.
2. Protein, nucleic acid, and water: H/D → H, water O/H or EP for virtual
   sites ([4-point water](../mlmm-model-setup/formats.md#amber-parm7-and-rst7)), Se → Se, a first letter of P/N/O/S, then C* → C.
3. Other ligands: the column alignment separates ` NA ` (N) from `NA  `
   (Na); LEaP ` CL1` and ` BR1` are halogens, and `HG11` is H.
4. Otherwise a two-letter, then one-letter match against known elements
   (D → H); still ambiguous means unassigned.

**Pitfalls.**

- `[add-elem-info] WARNING: Could not confidently assign N atoms; left unchanged.`
  lists up to 50 atoms. Type each symbol into columns 77–78 by hand,
  right-aligned; `--overwrite-elem` uses the same rules and will not assign
  them.
- Unusual ligand atom names can be misclassified; spot-check and diff the
  output.
- `--overwrite-elem` decides which fields are re-inferred; `--overwrite`
  decides only where the file is written.

## bond-summary

**When to use.** List the covalent bonds formed and broken between
consecutive structures (R → P, or R → IM1 → IM2 → P), with the criterion
that `irc` and `all` use for their bond changes and that `path-search` uses
to split segments. Use it on R and P, on IRC endpoints, or to see where
`path-search` split a mechanism.

**Minimal run.**

```bash
mlmm bond-summary -i reactant.pdb product.pdb
mlmm bond-summary -i 1.R.pdb -i 3.P.pdb
mlmm bond-summary -i frame_01.xyz frame_05.xyz frame_10.xyz
```

Give two or more XYZ, PDB, or GJF files in order, after one `-i` or with
`-i` repeated. `--bond-factor` (1.20) scales the summed covalent radii;
numbering is `--one-based` (default) or `--zero-based`; `--device` defaults
to `cpu`.

**Judge success.** One block per consecutive pair, written to stdout only:

```
============================================================
  1.R.pdb  →  3.P.pdb
============================================================
Bond formed (2):
  - O14-H106 : 1.502 Å --> 1.011 Å
  - P95-O107 : 3.477 Å --> 1.523 Å
Bond broken (2):
  - P95-O97 : 1.585 Å --> 3.270 Å
  - H106-O107 : 1.034 Å --> 1.673 Å
```

`Bond formed: None` means no bond formed. `--json` prints JSON to stdout
instead, with `scientific_status`, `execution_status`, and `comparisons`
(per pair: `structure_a`, `structure_b`, the counts `bonds_formed` and
`bonds_broken`, and the lists `formed` and `broken`); redirect it to keep it.

**Pitfalls.**

- Every input needs the same atoms in the same order; a pair that differs
  prints `ERROR: Atom types and ordering must be identical.` Neither
  `bond-summary` nor `extract` reorders atoms (`extract` only detects a
  mismatch and stops), so make the inputs from the same topology. A pair
  that cannot be compared gives exit code 1.
- The test is geometric: a pair is bonded at no more than 0.95 × T, where
  T is the summed covalent radii times `--bond-factor`, and a change counts
  only if the distance moves by at least 0.05 × T. It does not tell covalent
  from ionic or hydrogen bonds. Metal–ligand contacts can sit near the
  cutoff; raise `--bond-factor`, for example to 1.30 for metal coordination
  at 2.0–2.4 Å, or to 1.5–1.6 for permissive detection.
- In the `summary.json` of `all` and `path-search`,
  `segments[i]["bond_changes"]` holds the same report as a multi-line
  string.

## trj2fig

**When to use.** Plot the energy along an XYZ trajectory, such as
`optimization_trj.xyz` (`opt --dump`), `scan_trj.xyz`, `mep_trj.xyz`, or
`finished_irc_trj.xyz`, as a figure or a CSV table. For a labeled R/TS/IM/P
diagram from numbers, use [energy-diagram](#energy-diagram).

**Minimal run.**

```bash
mlmm trj2fig -i finished_irc_trj.xyz -o irc_profile.png
mlmm trj2fig -i scan_trj.xyz -o mep.html
mlmm trj2fig -i scan_trj.xyz -o profile.png -o profile.csv
mlmm trj2fig -i trajectory.xyz -q 0 -m 1 -b uma \
    --backend-model uma-s-1p2 --precision fp32 -o profile.png --out-json
```

The default output is `energy.png`. Repeat `-o`, or list more paths after
it; the extension picks the format (`.png`, `.jpg`, `.jpeg`, `.html`,
`.svg`, `.pdf`, or `.csv`). `-r` is `init` (default: the frame at the left
end), `None` (absolute energies), or a 0-based frame index. `--unit` is
`kcal` (default) or `hartree`; `--reverse-x` puts the last frame on the
left.

**Judge success.** The console prints `[trj2fig] Saved figure -> energy.png`.
The CSV has `frame`, `energy_hartree`, and `delta_kcal` or `delta_hartree`
(`energy_kcal` or `energy_hartree` with `-r None`). With `--out-json`,
`result.json` in the directory of the first output has `n_frames`,
`min_energy_hartree`, `max_energy_hartree`, `energy_source`
(`trajectory_comment` or `mlip_recomputed`), and `output_files`; the MLIP
backend, model, and precision are filled only when energies were recomputed.

**Pitfalls.**

- Each comment line needs a lone decimal or scientific value (`-1234.56`,
  `-1.23e3`), an `E=<value>` or `Energy=<value>` field (optional unit `Ha`,
  `Eh`, `hartree`, `eV`, or `kcal/mol`; hartree without one, except
  `energy=` in extended XYZ, which is eV), or the pysisyphus
  `<value> , ...` form. Other text is ambiguous and an integer-only comment
  may be a frame index; both stop the run with an error naming the frame.
- `-q` or `-m` recomputes every frame with the MLIP over all atoms, with no
  ML/MM split, so the result is not the ONIOM energy. With only one given,
  the other is charge 0 or multiplicity 1. For the ONIOM energy of an mlmm
  trajectory, give neither.

## energy-diagram

**When to use.** Draw a state energy diagram from values you already have,
for example R from one run and TS or IM from another, DFT single points, or
the `energy_diagrams` of a `summary.json`. It reads no structure and runs
no calculation.

**Minimal run.**

```bash
mlmm energy-diagram -i 0.0 -i 21.5 -i -0.7 -i 2.2 -i -18.2 -o diagram.png
mlmm energy-diagram -i "[0.0, 21.5, -0.7, 2.2, -18.2]" -o diagram.png
mlmm energy-diagram -i "[0, 12.5, 4.3]" --label-x "['R','TS','P']" \
    --label-y "ΔE (kcal/mol)" -o energy.png
```

Give one `-i` per value, the values after one `-i`, or one list-like
string. The default output is `energy_diagram.png`; `.png`, `.jpg`,
`.jpeg`, `.svg`, and `.pdf` work, and a path without an extension gets
`.png`. Without `--label-x`, states are `S1`, `S2`, … in input order.

**Judge success.** The console prints `[energy-diagram] Saved -> …`. With
`--out-json`, `result.json` next to the image has `n_points` and the image
path, not the values or labels.

**Pitfalls.**

- Fewer than two values stop the run with
  `Provide at least two numeric values with -i/--input.`
- The number of `--label-x` labels must equal the number of values.
- Values are drawn as given, with no unit conversion. Put them all in one
  unit and state it in `--label-y` (default `ΔE (kcal/mol)`).

## Next step

- [opt.md](opt.md), [tsopt.md](tsopt.md), and [freq.md](freq.md): relax,
  find a TS, or run a vibrational analysis after `sp`.
- [dft.md](dft.md): replace the ML level with a DFT single point.
- [extract.md](extract.md): cut the binding pocket from the cleaned PDB.
- [irc.md](irc.md), [path.md](path.md), and [scan.md](scan.md): the
  trajectories and bond changes that `trj2fig` and `bond-summary` read.
- [PDB](../mlmm-model-setup/formats.md#pdb): the altLoc and element
  columns.
- [outputs.md](../mlmm-overview/outputs.md): per-segment bond changes and
  energies in `summary.json`.
