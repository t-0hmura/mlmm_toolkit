# Installing mlmm itself (core.md)

`mlmm-toolkit` needs no C/C++ build at `pip install` time, but the bundled
`hessian_ff` JIT-compiles its C++ kernels at first run via
`torch.utils.cpp_extension` (Ninja is pulled in as a dependency; a working
C++ compiler must be on `PATH`). The bundled `pysisyphus` (GPU-tensor fork),
`thermoanalysis`, and `hessian_ff` install automatically with
`pip install mlmm-toolkit` as separate top-level packages alongside `mlmm`.

## Prerequisites

- Python ≥ 3.11
- A working PyTorch install matching your CUDA driver — see `env-cuda.md`
- (For DFT) PySCF / GPU4PySCF — see `dft.md`
- (Optional) xtb binary for `--embedcharge` correction — see `xtb.md`
## Install from PyPI (recommended)

```bash
conda activate <YOUR_ENV>
pip install mlmm-toolkit                         # core only (UMA)
pip install 'mlmm-toolkit[orb,aimnet,dft]'        # extras as needed
```

Available extras (canonical list lives in `pyproject.toml`):

| Extra | Pulls in | When you need it |
|---|---|---|
| (none) | UMA via `fairchem-core`, base deps | Default; UMA backend works out of the box |
| `[orb]` | `orb-models` | Using `-b orb` |
| `[aimnet]` | `aimnet>=0.2.0` | Using `-b aimnet2` |
| `[dft]` | `pyscf>=2.13.0`, `gpu4pyscf-cuda12x>=1.7.0` (x86_64), `cupy-cuda12x`, `basis-set-exchange` | `mlmm dft` subcommand |
| `[dev]` | `pytest` family | Contributing |

`[mace]` does **not** exist as an extra because MACE conflicts with
`fairchem-core`'s `e3nn` pin — install `mace-torch` manually in a
**separate environment** (see `mace.md`).

To inspect the live extras list without opening `pyproject.toml`:

```bash
python -c "import importlib.metadata as m; print(m.metadata('mlmm-toolkit').get_all('Provides-Extra'))"
```

## Install from source (development)

```bash
git clone <repo-url-from-the-toolkit-README> mlmm
cd mlmm
pip install -e '.[orb,aimnet,dft]'
```

`pip install -e .` will pick up edits to the source tree without
re-installing. Useful when you are debugging the toolkit itself.

## Verify the install

```bash
mlmm --version
mlmm --help

# Sanity import
python -c "
import mlmm
import mlmm.core.defaults as d
print('version :', mlmm.__version__)
print('defaults:', sorted(n for n in dir(d) if not n.startswith('_'))[:10], '...')
"
```

`mlmm --help` should list 22 subcommands (`all`, `mm-parm`, `extract`,
`path-search`, `path-opt`, `opt`, `sp`, `tsopt`, `freq`, `irc`, `dft`,
`scan`, `scan2d`, `scan3d`, `oniom-export`, `oniom-import`,
`define-layer`, `trj2fig`, `energy-diagram`, `add-elem-info`,
`fix-altloc`, `bond-summary`).

## Where things live after install

```bash
python -c "import mlmm, os; print(os.path.dirname(mlmm.__file__))"
```

Inside that directory:

| File / dir | Purpose |
|---|---|
| `cli/app.py` | Click entry point |
| `core/defaults.py` | All default kwarg dicts (read with `import mlmm.core.defaults`) |
| `backends/` | UMA / Orb / MACE / AIMNet2 calculator factories |
| `workflows/{extract,path_search,tsopt,irc,freq,dft,all}.py` | Subcommand implementations |

`pysisyphus/`, `thermoanalysis/`, and `hessian_ff/` are **not** inside `mlmm/` —
they install as separate top-level packages (siblings of `mlmm/` in
`site-packages/`).

## Upgrading

```bash
pip install --upgrade mlmm-toolkit
mlmm --version                    # confirm new version
```

When upgrading **across minor versions**, also re-check:

- `python -c "import mlmm.core.defaults as d; print(d.RSIRFO_KW)"` — keys may have moved
- `mlmm <subcommand> --help` — flag set may have changed
- `summary.json` schema — see `mlmm-workflows-output/SKILL.md`

## Uninstall / clean rebuild

```bash
pip uninstall mlmm-toolkit
# Or scrap the env entirely:
conda env remove -n <YOUR_ENV>
```

## See also

- `env-cuda.md` — torch / CUDA setup that must come first
- `uma.md`, `orb.md`, `mace.md`, `aimnet2.md` — backend-specific extras
- `dft.md` — `[dft]` install details and aarch64 fallback