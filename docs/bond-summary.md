# `bond-summary` (bond changes between structures)

`bond-summary` reports which **covalent bonds form and break** between consecutive structures, such as reactant (R) → product (P), or R → intermediates IM1 → IM2 → P. For *N* input files it prints *N* − 1 comparison blocks (A → B, B → C, …) with each changed bond and its distance before and after, in Å. It writes no file.

---

## What it is for

* **Checking IRC endpoints**: confirm that the two ends of an intrinsic reaction coordinate (IRC) differ by the intended bonds.
* **Screening multistep mechanisms**: list the bonds that change in each step of a chain of intermediates.
* **Checking a workflow by hand**: compare the R, transition state (TS), and P that `all` wrote under [`segments/seg_NN/`](output-layout.md).

---

## Examples

### 1. Compare two structures

Compare the reactant and the product; replace the file names with your own structures.

```bash
mlmm bond-summary -i reactant.xyz product.xyz
```

The report lists the bonds under `Bond formed (k):` and `Bond broken (k):`; `Bond formed: None` means that no bond formed.

```text
============================================================
  reactant.xyz  →  product.xyz
============================================================
Bond formed (2):
  - O14-H106 : 1.502 Å --> 1.011 Å
  - P95-O107 : 3.477 Å --> 1.523 Å
Bond broken (2):
  - P95-O97 : 1.585 Å --> 3.270 Å
  - H106-O107 : 1.034 Å --> 1.673 Å
```

### 2. A multistep chain

Four structures give three blocks: R → IM1, IM1 → IM2, and IM2 → P.

```bash
mlmm bond-summary -i reactant.xyz im1.xyz im2.xyz product.xyz
```

---

## How it works

1. **Reading the structures**:
The files are read in the order given; XYZ, PDB, and GJF are recognized by their extension. At least two files are needed.
2. **Bond criterion**:
A pair is bonded when its distance is at most 0.95 *T*, where *T* is `--bond-factor` (default `1.20`) times the sum of the two covalent radii of [Cordero et al. (2008)](https://doi.org/10.1039/b801115j), with 0.40 Å for H.
3. **Counting a change**:
A pair counts as formed (broken) when it is unbonded (bonded) in the first structure, bonded (unbonded) in the next, and its distance changes by at least 0.05 *T*, so a small move across the threshold is not counted. `irc` and `all` use the same criterion for the bond changes they report.

---

## Output files

`bond-summary` writes no files. It prints one text block per consecutive pair to stdout, as in example 1; atom labels are the element and the atom index (1-based by default). With `--json`, it prints a JSON object to stdout instead, with `scientific_status` (`success`, `partial`, or `failed`), `execution_status`, and `comparisons`: for each pair, `structure_a`, `structure_b`, the counts `bonds_formed` and `bonds_broken`, and the lists `formed` and `broken` with the atom indices, elements, and both distances. Redirect stdout to keep it.

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | Two or more structures in order, listed after one `-i` (`-i` may also be repeated) |
| `--bond-factor` | float | `1.20` | Scale factor for the sum of covalent radii |
| `--json/--no-json` | flag | `False` | Print JSON to stdout instead of the text report |
| `--one-based/--zero-based` | flag | `--one-based` | Atom numbering in the report |

See the [generated CLI reference](reference/commands/bond_summary.md) for every option.

---

## Notes

* **Same atoms in the same order**: every structure must have the same atoms in the same order. A pair that differs prints `ERROR: Atom types and ordering must be identical.` to stderr (see {ref}`Input / extraction <input--extraction>`).
* **Borderline bonds**: to count longer contacts such as metal coordination at 2.0–2.4 Å, raise `--bond-factor` (for example `1.30`).
* **Failed pairs**: if any pair cannot be compared, the run exits with code 1 and the JSON `scientific_status` is `partial` or `failed`.
* **Exit codes**: see {ref}`Exit codes <exit-codes>`.

---

## See also

* [irc](irc.md) — IRC, whose endpoints are checked with the same criterion
* [all](all.md) — the full workflow, which reports bond changes for each step
* [trj2fig](trj2fig.md) — plot the energy profile of a trajectory
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
