# アーキテクチャ: mlmm-toolkit

---

## 1. 概要

`mlmm-toolkit` は、完全なタンパク質環境に対して **ML/MM (ONIOM) 酵素反応経路解析** を実行する Python 製 CLI です。ここでの ML/MM とは、小さな反応コアを機械学習原子間ポテンシャル (ML) で、周囲のタンパク質を分子力学 (MM) 力場で扱い、両者を subtractive ONIOM (Our own N-layered Integrated molecular Orbital and molecular Mechanics) エネルギースキームで結合したハイブリッドモデルを指します。

入力は PDB と基質名です。これらから、本ツールは parm7 トポロジーを自動生成し、ONIOM 領域分割 (ML / Movable MM / Frozen MM) を B-factor チャネルにエンコードします。続いて、マクロ/マイクロ交互スキームによる全系のHessianベース TS (遷移状態) 探索を実行します。

結果として、ステージパイプライン `extract → MM-param → ONIOM model → MEP → tsopt → IRC → freq → dft` による完全な反応経路が生成されます。ここで MEP は最小エネルギー経路、IRC は内在反応座標です。

このパッケージは **6 つの物理レイヤーディレクトリ** (`cli/`、`workflows/`、`domain/`、`backends/`、`io/`、`core/`) として構成されており、それぞれの役割と依存方向は後述の §4 レイヤー表にまとめています。外部コードはレイヤーディレクトリから直接インポートします (`from mlmm.backends.mlmm_calc import MLMMCore`、`from mlmm.core.utils import …`、`import mlmm.io.trj2fig` など)。従来のフラットトップ shim レイヤーは本リリースで廃止されました。

3 つの内蔵フォーク (`pysisyphus/`、`thermoanalysis/`、`hessian_ff/`) は repo-internal モジュールとしてリポジトリのトップに置かれています。これらは意図的に upstream の PyPI 配布物では **ありません** (`hessian_ff/` に至っては upstream がまったく存在しないため、バンドルが必須です)。§6 を参照してください。

---

## 2. レイヤー構造 (6 つの物理ディレクトリ)

### 2.1 レイヤー表

| layer | dir | responsibility | may depend on |
|---|---|---|---|
| **L1 Interface** | `mlmm/cli/` | Click ルートグループ、デコレータファクトリ、`--help-advanced`、bool フラグ正規化、サブコマンドリゾルバ、AmberTools preflight | `workflows/`、`core/` |
| **L2 Application** | `mlmm/workflows/` | サブコマンドごとのオーケストレーション。ステージランナー 1 つにつき 1 ファイル (`all.py`、`path_search.py`、`tsopt.py`、`extract.py`、`oniom_export.py`、`mm_parm.py`、…) | `domain/`、`backends/`、`io/`、`core/` |
| **L3 Domain** | `mlmm/domain/` | 化学を意識したヘルパーロジック (結合変化検出、結合サマリー、元素情報伝播) | `core/` |
| **L4a Infra (MLIP + ONIOM)** | `mlmm/backends/` | MLIP バックエンドディスパッチャ + バックエンドごとのアダプタ + ML/MM ONIOM 計算コア | `core/` |
| **L4b Infra (I/O)** | `mlmm/io/` | 出力レイアウト、サマリー、軌跡、PDB 修正、エネルギー図、Hessianキャッシュ、解析的Hessian glue | `core/` |
| **L5 Foundation** | `mlmm/core/` | defaults (全 CLI デフォルトのソース)、utils (PDB / XYZ / プロットヘルパー)、将来の `errors.py` / `types.py` | (none) |
| (bundle, not a layer) | `<repo>/pysisyphus/`、`<repo>/thermoanalysis/`、`<repo>/hessian_ff/` | repo-internal フォーク (オプティマイザ / 熱化学 / 解析的 MM Hessian) | (sibling, layer-external) |

**依存方向 (一方向)**: `L1 → L2 → {L3, L4} → L5` (§2.1 レイヤー表に準拠)。この方向ルールは CI のマーカーカバレッジ (`.github/scripts/check_engineering_markers.py`) によって強制されます。内蔵フォークはレイヤーグラフの外側に位置し、絶対パッケージパス (`from pysisyphus.X import Y`、`from hessian_ff.analytical_hessian import …`) を通じて任意のレイヤーからインポートできます。

### 2.2 パッケージツリーの ASCII マップ

```
mlmm_toolkit/ [GH: t-0hmura/mlmm_toolkit]
├── pyproject.toml packages.find = ["mlmm*",...] (glob, frozen)
├── README.md / CONTRIBUTING.md / CHANGELOG.md
├── docs/
│ ├── architecture.md ← this file
│ └──... (Sphinx site, unchanged)
├── mlmm/ ← package body, 6-layer physical dir
│ ├── __init__.py PEP 562 lazy: _LAZY_IMPORTS + __getattr__
│ ├── __main__.py `from mlmm.cli.app import cli`
│ ├── _version.py / py.typed
│ │
│ ├── cli/ # === L1 Interface ===
│ │ ├── app.py Click group + _LAZY_SUBCOMMANDS registry (absolute paths)
│ │ ├── common_options.py @add_precision_option / @add_backend_model_option / @add_ml_charge_spin_options et al.
│ │ ├── decorators.py make_is_param_explicit, bool/YAML helpers, render_cli_exception
│ │ ├── help_pages.py --help-advanced pager
│ │ ├── bool_compat.py --flag / --no-flag normalisation
│ │ ├── default_group.py subcommand resolver, lazy module import
│ │ └── preflight.py AmberTools / conda env / GPU preflight
│ │
│ ├── workflows/ # === L2 Application ===
│ │ ├── all.py full pipeline orchestrator (extract → … → DFT)
│ │ ├── path_search.py / path_opt.py MEP search / COS wrapper
│ │ ├── tsopt.py / freq.py / irc.py / dft.py per-stage runners
│ │ ├── opt.py / scan.py / scan2d.py /
│ │ │ scan3d.py / scan_common.py ONIOM geometry opt / scans
│ │ ├── extract.py active-site extraction CLI
│ │ ├── define_layer.py ML / Movable MM / Frozen MM B-factor assignment
│ │ ├── mm_parm.py AmberTools-driven parm7 / rst7 generation
│ │ ├── oniom_export.py ONIOM input writer (Gaussian / ORCA)
│ │ ├── oniom_import.py ONIOM input reader (sanity / atom-name diff)
│ │ └── align_freeze.py Kabsch + frozen-subset rmsd
│ │
│ ├── domain/ # === L3 Domain ===
│ │ ├── bond_changes.py R↔P bond detection
│ │ ├── bond_summary.py post-IRC diagnostic
│ │ └── add_elem_info.py PDB element column normaliser
│ │
│ ├── backends/ # === L4a Infra (MLIP + ONIOM) ===
│ │ ├── __init__.py --precision routing (apply_precision_to_calc_cfg)
│ │ ├── mlmm_calc.py ML/MM ONIOM calculator core (4 MLIP backends UMA / Orb / MACE / AIMNet2
│ │ inline; CHEMISTRY-RULE:1 / 2 / 8 / 9 host)
│ │ │ Future: split into base.py + per-backend uma.py / orb.py
│ │ │ / mace.py / aimnet2.py + ONIOM subdir
│ │ └── xtb_embedcharge_correction.py xTB point-charge embedding correction (--embedcharge)
│ │
│ ├── io/ # === L4b Infra (I/O) ===
│ │ ├── summary.py summary.json / summary.log writer
│ │ ├── energy_diagram.py Plotly diagram
│ │ ├── trj2fig.py trajectory → PNG / HTML / SVG / PDF
│ │ ├── pdb_fix.py altloc resolution
│ │ ├── hessian_cache.py in-memory Hessian cache
│ │ └── hessian_calc.py numerical-Hessian build + frequency / vibrational I/O helpers
│ │
│ ├── core/ # === L5 Foundation ===
│ │ ├── defaults.py C1 single source of truth for every default
│ │ ├── utils.py PDB / XYZ / plot helpers
│ │ ├── logging.py -v / -vv logging wiring
│ │ ├── calc_eval.py per-stage calc evaluation
│ │ └── residue_data.py residue tables
│ │
│ └── mcp/ # non-layer subpackage: MCP server exposing every CLI subcommand
│   ├── server.py / _runner.py
│   └── _tools.py
│
├── tests/ smoke / unit
├── .github/ workflows/ + scripts/ (docs-quality lint helpers; CI-only)
└── (repo-top sibling, layer-external bundled forks)
 pysisyphus/ ~90 file, repo-internal fork (slimmed; CLI driver + QM backends + wavefunction + dead optimisers / IRC / NEB variants removed)
 thermoanalysis/ 5 file, repo-internal fork
 hessian_ff/ 19 file / 4.2k LOC, NO upstream PyPI, mandatory bundling
```

### 2.3 レイヤーごとの責務詳細

**L1 `cli/`**。このレイヤーだけが Click コマンドを構築し argv をパースします。`app.py` はルートの `Click.Group` と `_LAZY_SUBCOMMANDS` レジストリを保持します — すべてのエントリが **絶対モジュールパス** (`mlmm.workflows.all`、`mlmm.io.trj2fig`、…) を使用するため、リゾルバは `default_group.py` 自体の置き場所に依存しません。`mlmm` 固有の `preflight.py` (AmberTools / conda env / GPU preflight) がここにあるのは、CLI 起動時、いかなる L2 ワークフローが呼び出されるよりも前に実行されるためです。

**L2 `workflows/`** (~21 ファイル)。サブコマンドごとに 1 ファイル。各ファイルは `cli` という名前の単一の `@click.command()` とそのプライベートヘルパーを所有します。大きなステージランナー (`all.py` = 4,414 LOC、`path_search.py` = 2,352 LOC、`tsopt.py` = 3,181 LOC、`extract.py` = 2,274 LOC、`oniom_export.py` = 2,027 LOC) は現在のレイアウトでは単一ファイルのまま残されています。将来的にはステージごとのサブディレクトリへ分割する可能性がありますが、これは **opt-in** であり、本リリースラインのスコープ外です。

**L3 `domain/`**。化学を意識したヘルパーロジックで、`torch` / `numpy` / `pysisyphus.constants` (数値バックエンド) はインポートしてよいですが、MLIP ランタイム (`fairchem`、`orb_models`、`mace`、`aimnet`) は **インポートできません**。この deny list は `.github/scripts/check_engineering_markers.py` (`_check_external_library_scope`) によってリポジトリ全体で強制されており、`backends/` 以外のモジュールでこれらのインポートを禁止します。別個の `# DOMAIN_PURE` モジュール docstring マーカーは、これとは異なる CI ゲート (`_check_domain_pure`) です。このマーカーは、MLIP-free を保つ必要があるバックエンド非依存の特定モジュール（`backends/mlmm_calc.py`、`workflows/tsopt.py`、`workflows/freq.py`、および `workflows/sp.py` に存在）を検出します。これ自体は deny-list 機構ではなく、`domain/` のファイルはどれもこのマーカーを持ちません。Domain ヘルパーは任意の L2 ステージランナーから再利用できます。

**L4a `backends/`**。ML/MM ONIOM 計算コア (`mlmm_calc.py` = 2,550 LOC) はバックエンドディスパッチ (`__init__.py`) およびスタンドアロンの xTB 点電荷埋め込み補正 (`xtb_embedcharge_correction.py`、`--embedcharge` で駆動) とともにここに置かれています。現時点では、ML 領域を評価する 4 つの MLIP バックエンド (UMA / Orb / MACE / AIMNet2) と OpenMM / hessian_ff カップリングはすべて `mlmm_calc.py` 内にインラインで配置されています。将来的には、これを MLIP レイヤー用に `backends/{base, uma, orb, mace, aimnet2}.py` へ、ONIOM コア用に `backends/mlmm_calc/` サブディレクトリ (`core.py`、`ase_calc.py`、`embed_charge.py`、`hessianff_calc.py`、`openmm_calc.py`、`facade.py`) へ分割する可能性があります。現在の単一ファイルの `mlmm_calc.py` は化学ルール **#1 (subtractive ONIOM)**、**#2 (link-atom Hessian B-matrix)**、**#8 (3-layer 5-pass partial Hessian)**、**#9 (parm7 atom indexing)** を保持しています — §5.1 を参照してください。

**L4b `io/`** (7 ファイル)。出力側の I/O 関連: ステージごとのサマリーライター、エネルギー図、軌跡レンダリング、PDB altloc 修正、Hessianキャッシュ、数値Hessian構築 + 振動数 / 振動 I/O (`hessian_calc.py`)、調和拘束のセットアップ。`io/` は決して `workflows/` に依存しません。出力フォーマットはここで所有され、ステージランナーから消費されます。

**L5 `core/`**。最下層です。`defaults.py` はすべての CLI デフォルトの **ソース** です — 数値を他のどこかに追加する前にここを grep してください。`utils.py` は PDB / XYZ / プロットヘルパーの ~3,200 LOC の寄せ集めで、将来的には `utils/{pdb,plot,coord,yaml,freeze,input_prep}.py` へ分割する可能性があります。`logging.py` (`-v` / `-vv` の配線)、`calc_eval.py` (ステージごとの calc 評価)、`residue_data.py` (残基テーブル) もここにあります。内部専用モジュール `errors.py`、`types.py` / `_stage.py` は、追加された時点でここに導入されます。

### 2.4 遅延インポート機構 (概念図)

```text
External consumer Package root Layer dir
------------------ ---------------- -----------

from mlmm.core.utils import x ────────────────────────────────────► mlmm/core/utils.py

import mlmm.io.trj2fig ──────────────────────────────────────────► mlmm/io/trj2fig.py

from mlmm.backends.mlmm_calc import ─────────────────────────────► mlmm/backends/mlmm_calc.py
 MLMMCore

from mlmm import MLMMCore ─────► mlmm/__init__.py
 __getattr__("MLMMCore")
 └─► _LAZY_IMPORTS["MLMMCore"]
 = "mlmm.backends.mlmm_calc"
 └─► importlib.import_module(...)
 └─► getattr(module, "MLMMCore")

mlmm myaction ─────────────────► mlmm/cli/app.py
 _LAZY_SUBCOMMANDS["myaction"]
 = ("mlmm.workflows.myaction", "cli", "...")
 └─► importlib.import_module(absolute path)
 └─► getattr(module, "cli") → Click command
```

2 つのインポートサーフェス (フラットトップ shim レイヤーは本リリースで
廃止されました。`from mlmm.<oldmod>` を使っていた下流コードはレイヤー化
パスへ移行する必要があります):

1. **レイヤー化インポートパス**: 外部コードはレイヤーディレクトリから直接インポートします (`from mlmm.backends.mlmm_calc import MLMMCore`、`from mlmm.core.utils import …`、`import mlmm.io.trj2fig` など)。
2. **ルートシンボル属性** (`from mlmm import MLMMCore`) — `mlmm/__init__.py:_LAZY_IMPORTS` + PEP 562 `__getattr__` によって処理されます。再エクスポートされる 5 つのシンボル (`MLMMCore`、`MLMMASECalculator`、`mlmm`、`mlmm_ase`、`mlmm_mm_only`) はすべて `mlmm.backends.mlmm_calc` に解決され、初回アクセス時にロードされるため、`import mlmm` は安価なまま保たれます (eager なのは `__version__` のみ)。ルートのモジュール属性サーフェスは **存在しません** — サブモジュールはトップレベルパッケージの属性としてではなく、フルパス (`import mlmm.io.trj2fig`) で到達します。

CLI サブコマンドリゾルバ (`cli/app.py:_LAZY_SUBCOMMANDS`) は **絶対** モジュールパス (例: `"mlmm.workflows.all"`) を使用するため、`default_group.py` を `cli/` へ移動してもサブコマンド発見が静かに壊れることはありません (レジストリはもはや `__package__` に依存しません)。

---

## 3. 初見者向け 5 ステップナビゲーション (合計 ≈ 40 分)

リポジトリを初めて開くコントリビュータは、上から下へこの経路をたどってください。各ステップは 1 つの関心事を完結させます。

| step | minutes | open | what you learn |
|------|---------|------|-----------------|
| 1 | 3 | [`README.md`](https://github.com/t-0hmura/mlmm_toolkit/blob/main/README.md) | 1 段落のエレベーターピッチ + 単一コマンドの使用法 |
| 2 | 5 | このファイル (`docs/architecture.md`) §2 + §4 | 6 レイヤーのディレクトリツリー、依存方向、各関心事の所在 |
| 3 | 5 | [`mlmm/cli/app.py`](../../mlmm/cli/app.py) | Click ルートグループ、`_LAZY_SUBCOMMANDS` レジストリ (≈ 22 エントリ)、絶対パス解決 |
| 4 | 20 | [`mlmm/workflows/all.py`](../../mlmm/workflows/all.py) (4,414 LOC, skim) | 1 つの完全なサブコマンドを上から下まで。`extract → mm-parm → ONIOM model → MEP → tsopt → IRC → freq → dft` をトレース |
| 5 | 7 | [`CONTRIBUTING.md`](https://github.com/t-0hmura/mlmm_toolkit/blob/main/CONTRIBUTING.md) §3 + §4 | 5 つの add-a-X レシピ + 「触るな」の隠れた制約 |

ステップ 5 のあとは、§4 のファイルインデックスをたどることで他のどのファイルも読めます。このパッケージは意図的に **各レイヤー内でフラット** です — `mlmm/<layer>/` 配下にネストしたパッケージはありません (`backends/mlmm_calc/` をバックエンドごとのモジュールへ分割する将来計画を除く)。したがって 2 ディレクトリより深くナビゲートする必要は決してありません。

---

## 4. ファイルインデックス — 「この関心事はどこにある?」

### 4.1 CLI / エントリ (L1 `cli/`)

| concern | file |
|---|---|
| Click ルートグループ + サブコマンドディスパッチ | `mlmm/cli/app.py` |
| サブコマンドリゾルバ (遅延インポート) | `mlmm/cli/default_group.py` |
| `python -m mlmm` shim | `mlmm/__main__.py` |
| 共有オプションデコレータファクトリ | `mlmm/cli/decorators.py` |
| `--help-advanced` pager | `mlmm/cli/help_pages.py` |
| Bool フラグ互換 (`--flag` / `--no-flag` + 値スタイル) | `mlmm/cli/bool_compat.py` |
| AmberTools / conda env / GPU preflight | `mlmm/cli/preflight.py` |

### 4.2 ワークフローステージランナー (L2 `workflows/`)

以下で用いる略語: MEP = 最小エネルギー経路、GSM = growing-string method、COS = chain-of-states、RSIRFO = restricted-step image-function rational-function optimization (RS-I-RFO とも表記)、Bofill = Bofill Hessian更新式、PHVA = partial Hessian vibrational analysis、IRC = 内在反応座標、Kabsch = Kabsch 剛体アラインメントアルゴリズム。

| concern | file |
|---|---|
| 完全パイプラインオーケストレータ | `mlmm/workflows/all.py` |
| 構造最適化 (ONIOM マクロ/マイクロ pre-opt) | `mlmm/workflows/opt.py` |
| 1D / 2D / 3D スキャン + 共有 | `mlmm/workflows/scan{,2d,3d,_common}.py` |
| MEP 探索 (GSM) | `mlmm/workflows/path_search.py` |
| MEP オプティマイザコア (pysisyphus COS) | `mlmm/workflows/path_opt.py` |
| TS 最適化 (RSIRFO + Bofill + マクロ/マイクロ) | `mlmm/workflows/tsopt.py` |
| 振動解析 (PHVA + UMA active block) | `mlmm/workflows/freq.py` |
| IRC 積分 (マクロ / マイクロ) | `mlmm/workflows/irc.py` |
| 単一点 DFT (gpu4pyscf サブプロセス、ONIOM 埋め込み) | `mlmm/workflows/dft.py` |
| 活性部位抽出 (cluster cap) | `mlmm/workflows/extract.py` |
| ML / Movable MM / Frozen MM 領域割り当て | `mlmm/workflows/define_layer.py` |
| AmberTools 駆動の MM パラメータ生成 | `mlmm/workflows/mm_parm.py` |
| ONIOM 入力ライター (Gaussian / ORCA) | `mlmm/workflows/oniom_export.py` |
| ONIOM 入力リーダー (sanity, atom-name diff) | `mlmm/workflows/oniom_import.py` |
| Kabsch / frozen-subset アラインメント | `mlmm/workflows/align_freeze.py` |

### 4.3 化学ヘルパー (L3 `domain/`)

| concern | file |
|---|---|
| R↔P 結合変化検出 | `mlmm/domain/bond_changes.py` |
| Post-IRC 結合サマリー | `mlmm/domain/bond_summary.py` |
| PDB 元素カラム正規化 | `mlmm/domain/add_elem_info.py` |

### 4.4 MLIP + ONIOM (L4a `backends/`)

| concern | file |
|---|---|
| ML/MM ONIOM 計算コア + 4 つのインライン MLIP バックエンド + ONIOM カップリング | `mlmm/backends/mlmm_calc.py` |
| `--precision` ルーティング (`apply_precision_to_calc_cfg` / `_PRECISION_DISPATCH`) | `mlmm/backends/__init__.py` |
| バックエンドディスパッチ / ファクトリ (`_create_ml_backend`) | `mlmm/backends/mlmm_calc.py` |
| xTB 点電荷埋め込み補正 (`--embedcharge`) | `mlmm/backends/xtb_embedcharge_correction.py` |
| バックエンドごとのアダプタ分割 (計画中、未実装) | `mlmm/backends/{base, uma, orb, mace, aimnet2}.py` |
| ONIOM コアサブディレクトリ (計画中、未実装) | `mlmm/backends/mlmm_calc/{core, ase_calc, embed_charge, hessianff_calc, openmm_calc, facade}.py` |

add-a-backend レシピについては [MLIP Backends](backends.md) を参照してください (現在は計画中のバックエンドごとの分割を前提としています。それが実装されるまで、バックエンドの追加は `mlmm_calc.py` のインライン変更を伴います)。

### 4.5 I/O (L4b `io/`)

| concern | file |
|---|---|
| `summary.json` / `summary.log` ライター | `mlmm/io/summary.py` |
| Plotly エネルギー図 | `mlmm/io/energy_diagram.py` |
| Trajectory → PNG / HTML / SVG / PDF | `mlmm/io/trj2fig.py` |
| PDB altloc 解決 | `mlmm/io/pdb_fix.py` |
| インメモリHessianキャッシュ (run ごとの TTL) | `mlmm/io/hessian_cache.py` |
| 数値Hessian構築 + 振動数 / 振動 I/O | `mlmm/io/hessian_calc.py` |
| 調和拘束のセットアップ | `mlmm/workflows/restraints.py` (L2 ステージヘルパー) |

### 4.6 Foundation (L5 `core/`)

| concern | file |
|---|---|
| **すべての CLI デフォルトのソース** | `mlmm/core/defaults.py` |
| PDB / XYZ / プロットヘルパー | `mlmm/core/utils.py` |
| `-v` / `-vv` ロギングの配線 | `mlmm/core/logging.py` |
| ステージごとの calc 評価 | `mlmm/core/calc_eval.py` |
| 残基テーブル | `mlmm/core/residue_data.py` |
| (将来) 内部 Protocol / TypedDict | `mlmm/core/types.py` |

### 4.7 Repo-internal 内蔵フォーク

| dir | role | divergent files (do NOT replace with upstream) |
|---|---|---|
| `pysisyphus/` | オプティマイザ / TS / IRC エンジン | `irc/IRC.py`、`optimizers/hessian_updates.py`、`tsoptimizers/TSHessianOptimizer.py`、`calculators/*` (合計 4 ファイル) |
| `thermoanalysis/` | 熱化学 (ΔG, ZPE, 分配関数) | `QCData.py` (upstream とのブランディング差分) |
| `hessian_ff/` | MM 力場上の解析的Hessian — **PyPI 404、バンドルは必須** | `analytical_hessian.py` (`mlmm/backends/mlmm_calc.py` が消費する唯一のエントリ) |

各ディレクトリの `README.md` で触るべきでない境界 (touch-restriction boundary) を確認してください。

---

## 5. 隠れた制約 (いかなるパッチの前にもこれを読むこと)

### 5.1 9 つの化学ルール (grep レシピ)

正しさにクリティカルな 9 つのルールが `backends/`、`workflows/`、`core/defaults.py` にまたがって存在します。これらは smoke テストでは検出 **されません** — ここでの静かなドリフトは反応経路の精度を壊します。インラインの `# CHEMISTRY-RULE:N` マーカーと `# DOMAIN_PURE` モジュール docstring マーカーがこれらのルールを識別し、`.github/scripts/check_engineering_markers.py` が CI でマーカーの完全性を強制します。

編集前にすべての化学ルールを見つけるには:

```bash
# List all 9 rule sites in the repo (host file + line)
grep -rnE '# CHEMISTRY-RULE:[0-9]+' mlmm/

# List every # DOMAIN_PURE marker (= chemistry-rule host modules)
grep -rn '# DOMAIN_PURE' mlmm/
```

9 つのルールはすべて `mlmm` に適用されます:

| # | rule | host file |
|---|---|---|
| 1 | Subtractive ONIOM エネルギー式 (`E = mm_real + ml_model − mm_model`) | `mlmm/backends/mlmm_calc.py` |
| 2 | Link-atom Hessian B-matrix 投影 | `mlmm/backends/mlmm_calc.py` |
| 3 | マクロ / マイクロ交互 (RS-I-RFO hess mode microiteration) | `mlmm/workflows/tsopt.py` |
| 4 | gpu4pyscf `rks_lowmem` トリプルガード | `mlmm/workflows/dft.py` |
| 5 | def2 ファミリーの自動 ECP 注入 | `mlmm/workflows/dft.py` |
| 6 | PHVA + UMA active-block partial Hessian | `mlmm/workflows/freq.py` |
| 7 | `bofill_update` advanced-indexing scatter | `mlmm/workflows/tsopt.py` |
| 8 | 3-layer 5-pass partial Hessian アセンブリ | `mlmm/backends/mlmm_calc.py` |
| 9 | parm7 アトムインデックス (1-based / serial gap handling) | `mlmm/backends/mlmm_calc.py` |

これらのいずれかを編集する場合は、`[CHEMISTRY-RULE:N]` コミットプレフィックスと HEAVY 階層の数値ゴールデンゲートの通過が必要です (`CONTRIBUTING.md` §1.1 を参照)。

**推奨される学習順序 (4 つの化学クラスタ)**:

| cluster | rules | shared concern | learn-first file |
|---|---|---|---|
| 5-pass Hessian セット | #1, #2, #8, #9 | subtractive ONIOM + link-atom B-matrix + 3-layer アセンブリ + parm7 インデックス | `mlmm/backends/mlmm_calc.py` (9 ルールのうち 4 つのホスト) |
| TS 最適化セット | #3, #7 | マクロ / マイクロ交互 + Bofill scatter | `mlmm/workflows/tsopt.py` |
| 振動セット | #6 | PHVA + UMA active-block partial Hessian | `mlmm/workflows/freq.py` |
| DFT セット | #4, #5 | gpu4pyscf 低メモリ + def2 ECP 注入 | `mlmm/workflows/dft.py` |

mlmm における実践的なカリキュラムは、まず 5-pass Hessian セット (#1, #2, #8, #9 — すべて `mlmm_calc.py` 内)、次に TS セット (#3, #7)、次に DFT (#4, #5)、最後に振動 (#6) です。

### 5.2 VRAM 管理の不変条件 (`del` チェーンをリファクタしないこと)

IRC / TSopt / Freq ステージは、CUDA メモリを解放するためにステージ間で GPU 常駐オブジェクト (`calc`、`geom`、`hess`) を明示的に `del` します。`all` ワークフローはさらにステージ境界で `gc.collect()` を実行します。**これらの `del` / `gc.collect()` 文をリファクタで取り除かないでください** — 完全なタンパク質環境での長時間 ML/MM `all` ジョブは、これらがないと OOM します。

### 5.3 内蔵フォーク: upstream を併存インストールしないこと

内蔵された `pysisyphus/`、`thermoanalysis/`、`hessian_ff/` パッケージは **フォーク** です (そして `hessian_ff/` の場合は、入手可能な **唯一** の配布物です — PyPI は 404 を返します)。このパッケージの隣に `pip install pysisyphus` や `pip install thermoanalysis` を再インストールすると、次が静かに壊れます:

- `pysisyphus/irc/IRC.py` — 初期変位のメモリ管理
- `pysisyphus/optimizers/hessian_updates.py` — advanced index 上の Bofill scatter、GPU OOM 回避のための CPU 専用 `bofill_update` パス
- `pysisyphus/tsoptimizers/TSHessianOptimizer.py` — RSIRFO kwargs
- `pysisyphus/calculators/...` — GPU を意識したバックエンドフック
- `thermoanalysis/QCData.py` — upstream とのブランディング / I/O 差分
- `hessian_ff/analytical_hessian.py` — `backends/mlmm_calc.py` が消費する唯一のエントリ。**upstream の代替は存在しません**

### 5.4 `pyproject.toml` の配列は 0-diff

`[tool.setuptools.packages.find].include` と `dependencies` は本リリース中、**0-diff 配列** として扱われます。`include` の glob (`mlmm*`) は新しいレイヤーサブパッケージをすでに自動発見します。`vendor/` や `internal/` のコンテナディレクトリを追加したり、新しいランタイム依存をピン留めしたりすると、インストール契約が壊れ、リリーススコープによって禁止されます。リフロー / コメント編集は問題ありません。**配列の内容** は凍結されています。

### 5.5 `_LAZY_SUBCOMMANDS` レジストリは絶対パスを使用すること

`mlmm/cli/app.py:_LAZY_SUBCOMMANDS` は、すべてのサブコマンドを **絶対** モジュールパスで解決します。いずれかのエントリを相対ドット付きインポート (`".all"` など) に戻すと、`default_group.py` が移動した際にサブコマンド発見が静かに壊れます。リゾルバの `__package__` がパッケージルートからドリフトするためです。内部設計ノートを参照してください。

---

## 6. 内蔵フォーク (repo-internal)

`mlmm_toolkit` はリポジトリのトップに **3 つ** の repo-internal モジュールを同梱します:

| dir | upstream PyPI? | purpose | scope of edits allowed |
|---|---|---|---|
| `pysisyphus/` | NO — フォーク、`pip install pysisyphus` を併存させない | オプティマイザ、TS、IRC、COS、calculators | 本リリースでは annotation のみ (docstring + 型ヒント)。ロジック編集は禁止 |
| `thermoanalysis/` | NO — フォーク (ブランディング差分) | ΔG, ZPE, 分配関数, `QCData` | `pysisyphus/` と同じ |
| `hessian_ff/` | **NO — PyPI 404、バンドル必須** | MM 力場上の解析的Hessian | `pysisyphus/` と同じ |

各ディレクトリは、分岐ファイルと触るべきでない境界 (touch-restriction boundary) を列挙する独自の `README.md` を持ちます。レイヤーモデルから見ると、これらのフォークは L1..L5 グラフの **外側** に位置します。任意のレイヤーが絶対パッケージパス (`from pysisyphus.X import Y`、`from hessian_ff.analytical_hessian import …`) を通じてこれらをインポートでき、`L1 → L2 → {L3, L4} → L5` の方向を壊しません。


---

## 7. 推奨される深掘り読書順序 (5〜10 ファイル)

初見者向けツアー (§3) のあとは、この深さ優先の読書順序に従ってください:

1. `mlmm/core/defaults.py` — デフォルト値テーブルを取り込む。下流のすべてがここから読み取る。
2. `mlmm/cli/app.py` — Click ルート + `_LAZY_SUBCOMMANDS` レジストリ。
3. `mlmm/workflows/all.py` — 1 つの完全なパイプラインを上から下まで。
4. `mlmm/workflows/extract.py` + `define_layer.py` — cluster cap + ONIOM レイヤー割り当て。
5. `mlmm/workflows/mm_parm.py` — AmberTools parm7 生成。
6. `mlmm/backends/mlmm_calc.py` — ML/MM の心臓部 (化学ルール #1, #2, #8, #9 がすべてここにある)。
7. `mlmm/workflows/tsopt.py` — RSIRFO + Bofill (CHEMISTRY-RULE:7) + マクロ / マイクロ交互 (CHEMISTRY-RULE:3)。
8. `mlmm/workflows/freq.py` — PHVA + UMA active-block (CHEMISTRY-RULE:6)。
9. `mlmm/workflows/irc.py` — VRAM 管理 + マクロ / マイクロ IRC。
10. `mlmm/core/utils.py` — 共有 PDB / XYZ / プロットヘルパー。

---

## 8. ML/MM (ONIOM) スコープ

`mlmm-toolkit` は ONIOM を介して **完全なタンパク質環境** を扱います:

- **ML 領域**: 基質 + 反応中心残基。4 つの MLIP バックエンド (UMA / Orb / MACE / AIMNet2) のいずれかで評価されます。オプションの xTB 点電荷埋め込み補正 (`--embedcharge`) は MM→ML の環境効果を加えます
- **Movable MM 領域**: ML 領域を取り囲むシェルで、AMBER 力場の下で自由に移動できます
- **Frozen MM 領域**: タンパク質の残りの部分で、剛体として保持されます

この分割は入力 PDB の B-factor チャネルにエンコードされ、`extract → mm-parm → ONIOM model → MEP → tsopt → IRC → freq → dft` を通じて伝播されます。
