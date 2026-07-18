---
name: mlmm-architecture
description: Where the source code lives in `mlmm-toolkit`. 6 physical layer directories (`cli` / `workflows` / `domain` / `backends` / `io` / `core`) + 3 repo-internal forks (`pysisyphus` / `thermoanalysis` / `hessian_ff`). Tells an agent which directory to grep for a given concern (Click option, ONIOM stage runner, MLIP backend, hessian-ff analytical MM Hessian, output writer, chemistry default, link-atom math) before touching code. TRIGGER on questions like "where is X implemented", "which file defines flag Y", "how is the repo organised", "what's safe to refactor". SKIP for usage questions — those belong to `mlmm-cli` / `-overview`.
---

# mlmm-toolkit architecture (one-screen map)

## 6 layers + bundled forks

```
mlmm/                              ← the package body, one folder per layer
├── cli/        # L1 — Click root group, --help-advanced, bool normalisation,
│               #      shared option-decorator factories, subcommand resolver,
│               #      AmberTools preflight, pysisyphus `mlmm` calculator
│               #      registration.
├── workflows/  # L2 — one file per CLI subcommand (`all.py`, `tsopt.py`,
│               #      `freq.py`, `irc.py`, `dft.py`, `extract.py`,
│               #      `define_layer.py`, `mm_parm.py`, `oniom_export.py`,
│               #      `oniom_import.py`, ...).
├── domain/     # L3 — chemistry-aware helpers (bond changes, bond summary,
│               #      element-info repair). No torch / no MLIP dependency.
├── backends/   # L4a — MLIP backend dispatcher + ML/MM ONIOM calculator core
│               #       (`mlmm_calc.py`, monolithic — UMA / Orb / MACE /
│               #       AIMNet2 / xTB + OpenMM + hessian_ff coupling) +
│               #       xTB QM/MM embed-charge correction.
├── io/         # L4b — summary writer, energy diagram, trajectory plot,
│               #       Hessian cache, analytical-Hessian glue, PDB altloc
│               #       fix. (harmonic restraints live in workflows/restraints.py, L2)
└── core/       # L5 — `defaults.py` (single source of truth for every CLI
                #       default), `utils.py` (PDB / XYZ / plot helpers),
                #       `logging.py`, `calc_eval.py`,
                #       `residue_data.py`.

pysisyphus/        ← bundled fork of the optimiser / TS / IRC engine.
                     Slimmed to the subset mlmm actually uses; see its
                     own README for the live divergent-file table. Routine
                     polish is annotation-only; defect/feature logic is gated.

thermoanalysis/    ← bundled fork for ΔG / ZPE / partition functions.
                     QCData.py is the only consumer; same touch restriction
                     as pysisyphus.

hessian_ff/        ← analytical Hessian on the MM force field (AMBER
                     ff14SB-style harmonic + LJ + Coulomb). NO upstream
                     PyPI package — bundling is mandatory. Single consumer
                     is `mlmm/backends/mlmm_calc.py`.
```

Dependency direction: the *design intent* is one-way `L1 → L2 → {L3, L4} → L5`. What is actually enforced is a cycle-free `mlmm` module graph, `core`/`domain` never importing `workflows`, and no bundled fork importing `mlmm` — checked by [`.github/scripts/check_import_graph.py`](../../.github/scripts/check_import_graph.py), which parses the real import edges (`check_engineering_markers.py` does not — it covers chemistry / `DOMAIN_PURE` / MLIP scope only). A few cycle-free back-edges remain today: some `core.utils` helpers reach into `domain`/`io`, and `core.calc_eval` into `backends`. The bundled forks sit outside the layer graph and may be imported from any layer.

## Where to look first

| concern | open this |
|---|---|
| Default for any CLI flag | `mlmm/core/defaults.py` (single source of truth — grep here before any other file) |
| Subcommand body / orchestration | `mlmm/workflows/<subcmd>.py` |
| New MLIP backend | extend `mlmm/backends/mlmm_calc.py` inline (per-backend split is a future refinement) |
| `--help` / option decorator | `mlmm/cli/common_options.py` (shared) or the subcommand file (inline) |
| ONIOM layer assignment (B-factor channels) | `mlmm/workflows/define_layer.py` |
| AMBER parm7 / rst7 generation | `mlmm/workflows/mm_parm.py` |
| ONIOM input deck (Gaussian / ORCA) | `mlmm/workflows/oniom_{export,import}.py` |
| MM analytical Hessian | `hessian_ff/analytical_hessian.py` (consumed by `mlmm_calc.py`) |
| Output schema (summary.json, trajectory, energy diagram) | `mlmm/io/` |
| Chemistry rule (subtractive ONIOM, link-atom Hessian, 5-pass partial Hessian, parm7 indexing) | search `# CHEMISTRY-RULE:` markers (lab-sign-off required to edit) |
| TS / IRC / optimiser internals | `pysisyphus/` (read its live divergence table and validation gate first) |
| MCP server / agent integration | `mlmm/mcp/` — see [`mlmm-mcp`](../mlmm-mcp/SKILL.md) |

## Hidden constraints to remember

1. **`mlmm/cli/app.py:_LAZY_SUBCOMMANDS`** entries MUST use absolute module paths (`"mlmm.workflows.all"`, never `".all"`). Relative dotted paths silently break the resolver if `default_group.py` moves.
2. **VRAM hygiene**: `# DO NOT INLINE` markers around `del calc; gc.collect(); torch.cuda.empty_cache()` between stages are load-bearing — removing them OOMs the next stage on full-protein ONIOM systems.
3. **`pyproject.toml [tool.setuptools.packages.find].include`** and `dependencies` arrays are treated as 0-diff for this release line. Adding a vendor / internal dir or pinning a new runtime dep breaks behaviour-level guarantees and is out of scope.
4. **Bundled-fork edits** use each directory README's live table. Routine polish is annotation-only; logic requires a demonstrated defect or approved numerical feature, focused regression tests, and the relevant HEAVY/GPU validation.

## See also
- Full architecture (~400 lines): [`docs/architecture.md`](../../docs/architecture.md)
- Contributor recipe + per-step gate cycle: [`CONTRIBUTING.md`](../../CONTRIBUTING.md)
- Engineering-marker coverage check (chemistry / `DOMAIN_PURE` / MLIP scope): [`.github/scripts/check_engineering_markers.py`](../../.github/scripts/check_engineering_markers.py)
- Import-graph gate (no cycles / no fork→product / no core·domain→workflows): [`.github/scripts/check_import_graph.py`](../../.github/scripts/check_import_graph.py)
