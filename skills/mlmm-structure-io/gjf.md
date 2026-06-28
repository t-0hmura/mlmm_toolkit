# Gaussian gjf format (gjf.md)

`gjf` (or `.com`) is Gaussian's input format. In `mlmm-toolkit` it is the
**Gaussian ONIOM exchange format** — not a generic geometry input:

- **written** by `mlmm oniom-export` (builds a g16 ONIOM input from a parm7 + layered PDB), and
- **read** by `mlmm oniom-import` (reconstructs an XYZ + a layer-encoded PDB from a Gaussian/ORCA ONIOM file).

The geometry pipeline does **not** read gjf. `opt` / `tsopt` / `freq` / `irc` /
`scan` / `path-opt` / `path-search` take **PDB / XYZ** via `-i` (other suffixes
are rejected with `Unsupported input format: … Use .pdb or .xyz`), and `dft`
takes a **PDB** (`-i`, with `--ref-pdb` for an XYZ). To feed a Gaussian file into
the pipeline, convert it first with `mlmm oniom-import`.

## Structure

```
%nproc=8
%mem=8GB
%chk=run.chk

# wB97X-D/def2-svp opt freq

  Title (one line, free text)

0 1
  C   0.000   0.000   0.000
  H   0.000   1.090   0.000
  ...

[blank line]
[optional: connectivity table, ECP, basis set, ...]
```

Layout:

| Block | Lines | Content |
|---|---|---|
| `%`-section | 0+ | Link0 commands: `%nproc`, `%mem`, `%chk`, … |
| Route line | 1 | Starts with `#`, declares method, basis, jobtype |
| (blank) | 1 | Required separator |
| Title | 1 | Free text |
| (blank) | 1 | Required separator |
| Charge / spin | 1 | Two integers separated by space: `<charge> <multiplicity>` |
| Coordinates | n | `<element>  <x>  <y>  <z>` (Å), or with frozen flag `<element> -1 <x> <y> <z>` |
| (blank) | 1 | Terminator |
| (optional) | 0+ | Connectivity, ECP block, custom basis, etc. |

## Export: generate a Gaussian ONIOM gjf

```bash
mlmm oniom-export --parm real.parm7 -i pocket.pdb --model-pdb ml.pdb -o out.com
```

`--mode g16` (or a `.gjf`/`.com` output suffix) selects Gaussian; `.inp` /
`--mode orca` selects ORCA. See `mlmm-cli/oniom-export.md`.

## Import: read a Gaussian ONIOM gjf back

```bash
mlmm oniom-import -i system.gjf -o system   # -> system.xyz + a layer-encoded PDB
```

`oniom-import` parses the ONIOM layer assignment / QM-MM regions from the
Gaussian (`.gjf`/`.com`) or ORCA (`.inp`) file and rebuilds the structures the
geometry pipeline consumes. See `mlmm-cli/oniom-import.md`.

## Charge / spin / frozen atoms

These are properties of the geometry pipeline, which reads PDB/XYZ — not gjf.
Charge / spin come from `-q` / `-m` (or the layer-charge summary); declare frozen
atoms with `--freeze-atoms` (CLI, comma-separated 1-based indices) or a YAML
`--config` `freeze_atoms: [<index>, ...]`. A `-1` frozen marker inside a Gaussian
file is **not** a generic-input mechanism; layer/active-atom information is handled
by `oniom-import`. See `charge-multiplicity.md`.

## Notes on output conversion

`--convert-files/--no-convert-files` is a **boolean** toggle (default on) that
writes **PDB companions** for XYZ/TRJ outputs — it takes no format value
(there is no `--convert-files gjf`). gjf output is produced by `oniom-export`,
not by the optimization subcommands.

## See also

- `pdb.md`, `xyz.md` — the geometry-pipeline input formats.
- `mlmm-cli/oniom-export.md`, `mlmm-cli/oniom-import.md` — the gjf producer / consumer.
- `charge-multiplicity.md` — charge / spin resolution.
