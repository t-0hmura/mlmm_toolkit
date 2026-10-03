---
name: mlmm-model-setup
description: "Choosing the mlmm-toolkit ML region and layers: what extract puts in the ML region (-c centers, -r radius, waters, backbone, --selected-resn), when to hand-build a link-H-free `model.pdb` and pass it with --model-pdb, --parm7, and -q, how define-layer splits movable and frozen MM by distance, how to cut cost (a smaller ML region, a shorter movable cutoff, an ML-only Hessian), how to enlarge (larger -r, added residues, a longer movable cutoff), and the atom-set rules for R/IM/P and WT/mutant models. TRIGGER on choosing the ML region, a run that is too slow or out of memory, a missing residue or water in the ML region, building `model.pdb`, or changing layers. SKIP for extract or define-layer flag syntax (mlmm-cli), B-factor encoding, parm7, and charge rules (mlmm-structure-io), and TS strategy (mlmm-overview)."
---

# Building the ML region and layers

Let `all -c` choose the ML region (`-r`, 2.6 Å) and the layers (`define-layer`, 8.0 Å movable cutoff), or pass a hand-built `--model-pdb` with `--parm7` and `-q`; trim with a smaller ML region, a shorter `--movable-cutoff`, or an ML-only Hessian, and enlarge with a larger `-r`, `--selected-resn`, or a longer cutoff.

```bash
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3'
```

Before a long run, check: `[all] define-layer [i]: … (ML=…, MovableMM=…, FrozenMM=…)` for each input, `[all] ML structure with link H (N + M; …)`, and `Total active site model charge`. Color by B-factor: the reacting residues should be in the ML region.

## Three layers and what they cost

- ML runs on the MLIP. Movable-MM and Frozen-MM run on the Amber force field; Frozen-MM atoms stay fixed but still count in the MM energy. B-factor encoding: [ML region and layers](../mlmm-structure-io/SKILL.md#ml-region-and-layers).
- `freq` and `tsopt` build the Hessian only over moving atoms (PHVA); frozen atoms are left out.
- Cost falls with a smaller ML region, fewer Movable-MM atoms, and fewer MM atoms in the Hessian.

## Build the ML region

- Give `-c` the substrate, cofactors, metals, and catalytic residues; with chain IDs, write `A:SAM:44`. The bundled examples lack chain IDs and use names.
- A residue joins when any of its atoms lies within `-r` of a `-c` atom; waters join by default. Consecutive amino acids keep their internal main chain; `--exclude-backbone` moves the main chain of non-center amino acids to MM. `--selected-resn` adds residues without a radius.
- `all` writes the first input's ML region to `<out-dir>/ml_region.pdb`; reuse it with `--model-pdb`.
- By hand: cut a link-H-free `model.pdb` from the PDB that `mm-parm` writes (order in [cli/extract.md](../mlmm-cli/extract.md)) and pass `--model-pdb`, `--parm7`, and `-q`. `--model-pdb` overrides `-c` and input B-factors; `--parm7` skips `mm-parm`.
- Automatic extraction derives the charge from residue names, `-l`, and `--modified-residue`; after hand edits to atoms, protonation, or the cut, give `-q`.
- Pitfalls: `--add-linkh` is only for a standalone pocket; the calculator adds link H on `parm7` boundary bonds. A resumed `all` writes `ml_region.pdb` to a temporary directory; keep the first run's copy.

## Check the boundary and the charge

- `model.pdb` selects atoms from the full PDB/`parm7`: keep atom order, names, numbers, and chain IDs; do not renumber or add link H.
- Include every atom in bond or proton transfer, plus covalent partners whose bonding changes.
- End retained backbone fragments at `CA` on both ends; put other cuts on aliphatic C–C single bonds (`CA–CB` or farther). Never cut peptide C–N, polar C–N/C–O, aromatic, disulfide, or metal-coordination bonds; move the boundary instead.
- Check boundary valences and the ML-region charge and multiplicity; `define-layer` cannot fix a bad selection.
- A boundary bond other than C–C, C–N, or N–C stops with `Unsupported ML/MM boundary bond in parm7`; move the cut.

## Set the layers

- `define-layer` puts MM atoms within `--movable-cutoff` (8.0 Å) of the ML region in Movable-MM and the rest in Frozen-MM. It is a freezing threshold: raise it to free more, lower it to lock more.
- With `-c`, `all` runs `define-layer` at 8.0 Å on each input; without `-c`, it keeps the input B-factor layers (`--detect-layer`, on by default).
- For another cutoff, run `define-layer --movable-cutoff` on the `mm-parm` PDB and pass the result to `all` without `-c`, with `--parm7` and `--model-pdb`.
- Pitfalls: `extract` does not assign layers. With `--detect-layer`, `--model-pdb` sets only the ML region; the MM layers come from the PDB B-factors, never from `parm7`. Rerun `define-layer` instead of editing them. `all` without `-c` under `--no-detect-layer` needs `--model-pdb`.

## Trim to lower cost

- Smaller ML region: a smaller `-r`, `--exclude-backbone`, `--no-include-h2o`, `-r 0` with `--selected-resn`, or a trimmed `model.pdb`. Recheck the charge.
- Fewer Movable-MM atoms: a shorter `--movable-cutoff`.
- ML-only Hessian in `freq` and `tsopt`: `--hessian-cutoff 0.0 --active-dof-mode ml-only`. Their analysis covers ML and all Movable-MM by default (`partial`), so a narrower `--hessian-cutoff` alone stops the run. `all` has no `--hessian-cutoff`.
- Microiteration, on by default in `tsopt` and `opt --opt-mode hess`, relaxes MM on the force field alone between ML steps, saving MLIP calls.
- Pitfalls: `--movable-cutoff` replaces the B-factor MM layers; in `opt`, `tsopt`, `freq`, `path-opt`, and `path-search` it also turns off `--detect-layer`, so pass `--model-pdb`.

## Enlarge when the model is too small

- Raise `-r`, add residues with `--radius-het2het`, `--selected-resn`, or `-c`, or add atoms to `model.pdb`. Lengthen `--movable-cutoff` to relax more of the environment.
- The radius is a convergence test: a larger region costs more and is not always better, so compare energies, forces, and barriers over a few sensible regions.

## Same atoms across states and variants

- R/IM/P: every full-system PDB has identical atoms and order, and one `model.pdb` serves all. `all` builds the ML region from the first input and layers every input with it.
- WT/mutant: build and parameterize each system separately, use corresponding ML and movable regions, transfer layer labels only for atoms with a clear match, assign added or deleted atoms explicitly, and set charge and multiplicity per system. Comparing barriers: [Controlled mutant-vs-WT comparison](../mlmm-overview/ts-strategy.md#7-controlled-mutant-vs-wt-comparison).

## Next step

- [cli/extract.md](../mlmm-cli/extract.md) and [cli/define-layer.md](../mlmm-cli/define-layer.md) — the two commands.
- [ML region and layers](../mlmm-structure-io/SKILL.md#ml-region-and-layers) — B-factor encoding.
- [Pick an all mode](../mlmm-overview/SKILL.md#pick-an-all-mode) — which `all` mode.
- [`docs/model-setup.md`](../../docs/model-setup.md) — freezing atoms and distance restraints.
