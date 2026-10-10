# `mlmm oniom-export` and `mlmm oniom-import`

Write an ML/MM system as a Gaussian ONIOM or ORCA QM/MM input, and read such
an input back as an XYZ and a layered PDB. Run
`mlmm oniom-export --mode g16 --parm7 real.parm7 -i result_tsopt/final_geometry.pdb --model-pdb ml_region.pdb -o ts_refine.com -q 0 -m 1`.
Success is `[oniom-gaussian] Wrote 'ts_refine.com'` with the QM-atom count
you expect; for `oniom-import`, the two `[oniom-import] wrote:` lines.

## When to use

- `oniom-export`: take a structure from mlmm, such as a TS candidate, into a
  Gaussian ONIOM calculation with a DFT high layer, run the same system with
  the QM/MM module of ORCA, or compare against a reference DFT/MM
  calculation. The ML region becomes the QM region; the MM part comes from
  the `parm7`.
- `oniom-import`: bring a Gaussian or ORCA QM/MM input that you prepared or
  edited back into mlmm, or move a hand-built `oniom(...)` setup to MLIP
  ML/MM. With `--ref-pdb`, an input written by `oniom-export` gets back the
  atom and residue names of the original PDB.

## Minimal run

```bash
mlmm oniom-export --parm7 enzyme.parm7 -i complex_layered.pdb \
    -o complex_oniom.gjf --method 'wB97XD/def2-SVP' -q 0 -m 1 \
    --nproc 16 --mem 32GB
mlmm oniom-export --parm7 enzyme.parm7 -i complex_layered.pdb \
    -o complex_oniom.inp --method 'B3LYP def2-SVP' -q 0 -m 1 \
    --total-charge 0 --total-mult 1
```

The suffix of `-o` picks the program (`.gjf` or `.com` for Gaussian, `.inp`
for ORCA); `--mode g16` or `--mode orca` takes precedence and is needed for
any other suffix. `--method` covers the QM region (default
`wB97XD/def2-TZVPD` for Gaussian, `B3LYP D3BJ def2-SVP` for ORCA).
`-q` is the QM-region charge and `-m` its multiplicity (default 1; set it
for radicals). `--nproc` (8) and `--mem` (16GB, Gaussian) set the resources,
and `--total-charge` and `--total-mult` set the ORCA `Charge_Total` and
`Mult_Total`. The QM region is `--model-pdb` when given, otherwise the atoms
with B-factor 0 in `-i`.

```bash
mlmm oniom-import -i oniom.gjf -o reconstructed
mlmm oniom-import -i oniom.gjf --ref-pdb original.pdb -o reconstructed
mlmm oniom-import -i oniom.inp --mode orca -o reconstructed
```

`-o` is a prefix: the command writes `reconstructed.xyz` and
`reconstructed_layered.pdb` (default prefix: the input stem, in the current
directory). The suffix of `-i` picks the mode unless `--mode` is given.

## Judge success

**oniom-export.** Gaussian prints `[oniom-gaussian] Wrote '<file>'`,
`QM atoms: N, Movable atoms: M`, and `Link boundaries: K`. ORCA prints
`[oniom-orca] Wrote '<file>'`, `QM atoms: N, Active atoms: M`, and
`Link boundaries (auto-capped by ORCA): K`, then
`[oniom-orca] ORCAFF.prms: <path>` when the force-field file is ready. Check
N against your ML region. The Gaussian route is
`#p oniom(<method>:amber=softonly)`, each coordinate row has the movable flag
(`0` movable, `-1` frozen) and the layer (`H` or `L`), and the charge line
holds the whole system (topology total charge and `-m`) and then the QM
region twice. A PDB input also writes the atom-order marker
`MLMM_REF_PDB_ORDER_V1_SHA256=<digest>` (a comment in ORCA).

**oniom-import.** The log prints `[oniom-import] mode=…`, the counts
`atoms=N, qm=…, movable=…, frozen=…`, and the two `wrote:` lines. With
`--ref-pdb` it also prints `ref_order=identity-verified` (the order matched
the export marker), `element-verified`, or `unverified-opt-in`. The XYZ
comment reads `mode=<mode> atoms=N qm=… movable=… q=<charge> m=<multiplicity>`.
The PDB carries B-factors 0, 10, and 20: in Gaussian, `H` atoms are ML, `L`
atoms with `0` are Movable-MM, and `L` atoms with `-1` are Frozen-MM; in
ORCA, `QMAtoms` and `ActiveAtoms` set the same layers.

## Pitfalls and recovery

- A `parm7` with CMAP terms stops the export before any file is written:
  Gaussian ONIOM cannot represent CMAP, and the MM engine of ORCA does not
  apply it. Build a CMAP-free topology for the export with
  `mlmm mm-parm -i input.pdb -l 'LIG:0' --ff-set ff14SB --out-prefix system`
  and confirm `not p.cmaps` in parmed. The ML/MM calculations in mlmm can
  keep CMAP, which they apply in both MM layers.
- `-i` must be a PDB with the same atoms as the `parm7`, in the same order.
  A different atom count stops the export even with `--no-element-check`;
  with the check on, the first different element stops it with
  `Element sequence mismatch at atom index …` (counted from 0). Use
  `--no-element-check` only when you already trust the order.
- The Gaussian whole-system charge, and the default ORCA `Charge_Total`, is
  the sum of the `parm7` partial charges rounded to an integer. A sum more
  than 0.05 from an integer stops the export; in ORCA mode, give
  `--total-charge`.
- In Gaussian, each cut QM–MM bond needs its own MM atom; two QM atoms bonded
  to one MM atom stop the export. `--link-atom-method scaled` (default) places
  the link H with the g-factor, and `fixed` at 1.09 Å (QM carbon) or 1.01 Å
  (QM nitrogen). ORCA builds the caps itself.
- The exported input has no job keyword, so as written it is a single point.
  Add `opt`, `freq`, or the job you want before running it.
- ORCA needs `ORCAFF.prms`: the `--orcaff` file, or
  `<parm7 stem>.ORCAFF.prms` next to `-o`. A missing file is created with
  `orca_mm -convff -AMBER` when `--convert-orcaff` is on (default) and
  `orca_mm` is on `PATH`. Otherwise the `.inp` is still written, the console
  prints `[oniom-orca] NOTE: ORCAFF.prms not found at '<path>'. Run manually: …`,
  and the `.inp` cannot run until you do. The `.inp` names the file by
  absolute path; check it after moving the input to another machine.
- Gaussian and ORCA are not part of mlmm-toolkit; install and license them
  separately.
- Without `--ref-pdb`, every imported atom is named by its element in
  residue `MOL` 1 of chain A. Pass the PDB that went into `oniom-export` to
  recover numbering, chains, residues, and atom names.
- The import checks the order first: elements must agree atom by atom. A
  malformed or repeated marker, or a digest that does not match the
  reference, always stops the run, even with `--allow-unverified-ref-order`;
  use the original PDB or re-export.
- An input without the marker, matched by position, is accepted when no
  element appears twice. With repeated elements, check the order yourself and
  pass `--allow-unverified-ref-order`, which needs `--ref-pdb`; never use it
  for a digest mismatch.
- Gaussian inputs need two-layer ONIOM rows with the movable-flag column, as
  `oniom-export` writes them; a non-ONIOM input without `H` or `L` markers
  stops with an error. ORCA inputs need `QMAtoms {…} end` and
  `ActiveAtoms {…} end`, each on one line inside `%qmmm`, and a
  `* xyz <charge> <multiplicity>` block. If the mode is guessed wrong, pass
  `--mode`.
- `oniom-import` reads input files, not output logs. To continue in mlmm from
  a Gaussian or ORCA optimization, write the final geometry in the topology
  atom order and use it with the original `parm7` and ML region.

## Next step

- [ML region and layers](../mlmm-model-setup/SKILL.md#ml-region-and-layers):
  the B-factor layers that the export reads and the import writes.
- [GJF](../mlmm-model-setup/formats.md#gjf) and
  [Amber parm7 and rst7](../mlmm-model-setup/formats.md#amber-parm7-and-rst7):
  the file formats.
- [mm-parm.md](mm-parm.md): build the topology.
- [define-layer.md](define-layer.md): write the layers into the PDB.
