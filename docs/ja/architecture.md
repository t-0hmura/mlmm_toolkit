# アーキテクチャ: mlmm-toolkit

## 1. 概要

mlmm-toolkit の開発者向けに、パッケージの層、ファイルの置き場、開発するときに守る制約をまとめたページです。変更したあとは、CONTRIBUTING の [Required validation](https://github.com/t-0hmura/mlmm_toolkit/blob/main/CONTRIBUTING.md#11-required-validation) の検査を走らせてください。計算を実行するだけなら、[はじめに](getting-started.md)から読んでください。

`mlmm-toolkit` は、完全なタンパク質環境に対して **ML/MM (ONIOM) 酵素反応経路解析** を実行する Python 製 CLI です。小さな反応コアを機械学習原子間ポテンシャル (MLIP) で、周囲のタンパク質を分子力学 (MM) 力場で計算し、両者を subtractive ONIOM で合わせます。

`all` workflow は `extract`、`mm-parm`、MEP 探索、TS 最適化、IRC、振動解析、DFT 一点計算をつなぎます。TS、熱化学、DFT は任意の段です。

同梱のフォーク `pysisyphus/`、`thermoanalysis/`、`hessian_ff/` はリポジトリ最上位にあります（§5.3、§6）。

---

## 2. レイヤー構造 (6 つの物理ディレクトリ)

### 2.1 レイヤー表

| 層 | ディレクトリ | 責務 | 依存してよい先 |
|---|---|---|---|
| **L1 Interface** | `mlmm/cli/` | Click ルートグループ、デコレータファクトリ、`--help-advanced`、bool フラグ正規化、サブコマンドリゾルバ、AmberTools preflight | `workflows/`、`core/` |
| **L2 Application** | `mlmm/workflows/` | サブコマンドごとのオーケストレーションと共有ワークフローヘルパー (`_all_helpers.py`、`_opt_freq_common.py`、`_run_session.py`、…) | `domain/`、`backends/`、`io/`、`core/` |
| **L3 Domain** | `mlmm/domain/` | 化学を意識したヘルパーロジック (結合変化検出、結合サマリー、元素情報伝播) | `core/` |
| **L4a Infra (MLIP + ONIOM)** | `mlmm/backends/` | MLIP バックエンドのディスパッチ、インライン実装、および ML/MM ONIOM 計算コア | `core/` |
| **L4b Infra (I/O)** | `mlmm/io/` | 出力レイアウト、サマリー、軌跡、PDB 修正、エネルギー図、Hessian キャッシュ、解析的 Hessian glue | `core/` |
| **L5 Foundation** | `mlmm/core/` | 共有デフォルト、PDB/XYZ/プロットヘルパー、出力・結果確定処理、残基テーブル | `backends/`、`domain/`、`io/`、`cli/` (後述の一部の上向きの import) |
| （レイヤー外の同梱物） | `<repo>/pysisyphus/`、`<repo>/thermoanalysis/`、`<repo>/hessian_ff/` | リポジトリ内フォーク（オプティマイザ / 熱化学 / 解析的 MM Hessian） | （同階層、レイヤー外） |

**依存の向き（設計目標）**: `L1 → L2 → {L3, L4} → L5`。共有の charge/spin 準備とレイヤーヘルパーは `workflows/charge_prep.py` と `workflows/_opt_freq_common.py` にあります。同梱フォークは層の外にあり、どの層からも `from pysisyphus.X import Y` の形で import できます。

CI が検査するのはこの向きの一部だけです。

- `.github/scripts/check_import_graph.py` は、`mlmm` のモジュール間の import の循環、`core` と `domain` から `workflows` への import、同梱フォークから `mlmm` への import を禁止します。
- `.github/scripts/check_engineering_markers.py` は、`# CHEMISTRY-RULE` と `# DOMAIN_PURE` のマーカー（§5.1）と、MLIP ランタイムを `backends/` の下でだけ import していることを検査します。

### 2.2 パッケージツリーの ASCII マップ

```
mlmm_toolkit/ [GH: t-0hmura/mlmm_toolkit]
├── pyproject.toml packages.find = ["mlmm*",...] (パッケージ検出 glob)
├── README.md / CONTRIBUTING.md / CHANGELOG.md
├── docs/
│ ├── architecture.md ← this file
│ └──... (Sphinx ドキュメントサイト)
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
│ │ ├── bool_compat.py --flag / --no-flag normalization
│ │ ├── default_group.py subcommand resolver, lazy module import
│ │ └── preflight.py AmberTools / conda env / GPU preflight
│ │
│ ├── workflows/ # === L2 Application ===
│ │ ├── all.py full pipeline orchestrator (extract → … → DFT)
│ │ ├── path_search.py / path_opt.py MEP search / COS wrapper
│ │ ├── tsopt.py / freq.py / irc.py / dft.py / sp.py per-stage runners
│ │ ├── opt.py / scan.py / scan2d.py /
│ │ │ scan3d.py / scan_common.py ONIOM geometry opt / scans
│ │ ├── extract.py active-site extraction CLI
│ │ ├── define_layer.py ML / Movable-MM / Frozen-MM B-factor assignment
│ │ ├── mm_parm.py AmberTools-driven parm7 / rst7 generation
│ │ ├── oniom_export.py ONIOM input writer (Gaussian / ORCA)
│ │ ├── oniom_import.py ONIOM input reader (sanity / atom-name diff)
│ │ ├── align_freeze.py Kabsch + frozen-subset rmsd
│ │ └── _all_helpers.py / _opt_freq_common.py / _run_session.py /
│ │     restraints.py / charge_prep.py shared workflow helpers
│ │
│ ├── domain/ # === L3 Domain ===
│ │ ├── bond_changes.py R↔P bond detection
│ │ ├── bond_summary.py post-IRC diagnostic
│ │ └── add_elem_info.py PDB element column normalizer
│ │
│ ├── backends/ # === L4a Infra (MLIP + ONIOM) ===
│ │ ├── __init__.py --precision routing (apply_precision_to_calc_cfg)
│ │ ├── mlmm_calc.py ML/MM ONIOM calculator core (4 MLIP backends UMA / ORB / MACE / AIMNet2
│ │ inline; CHEMISTRY-RULE:1 / 2 / 8 host)
│ │ ├── custom.py user ASE calculator loaded from --calc-file (custom backend)
│ │ ├── pyscf_dft.py 任意の PySCF/GPU4PySCF high-level adapter
│ │ └── _determinism.py strict-determinism setup (--deterministic)
│ │
│ ├── io/ # === L4b Infra (I/O) ===
│ │ ├── summary.py summary.json / summary.log writer
│ │ ├── energy_diagram.py Plotly diagram
│ │ ├── trj2fig.py trajectory → PNG / HTML / SVG / PDF
│ │ ├── pdb_fix.py altloc resolution
│ │ ├── pdb_indexing.py parm7 atom indexing (CHEMISTRY-RULE:9)
│ │ ├── hessian_cache.py in-memory Hessian cache
│ │ └── hessian_calc.py numerical-Hessian build + frequency / vibrational I/O helpers
│ │
│ ├── core/ # === L5 Foundation ===
│ │ ├── defaults.py shared workflow/calculator defaults
│ │ ├── dft_settings.py DFT settings (CHEMISTRY-RULE:4)
│ │ ├── utils.py PDB / XYZ / plot helpers
│ │ ├── logging.py -v/--verbose LEVEL（0–3）ロギング配線
│ │ ├── calc_eval.py per-stage calc evaluation
│ │ ├── output.py / result_commit.py output/result commit helpers
│ │ ├── pes_composition.py energy-component composition
│ │ └── residue_data.py residue tables
│ │
│ └── mcp/ # non-layer subpackage: MCP server exposing every CLI subcommand
│   ├── server.py / _runner.py
│   └── _tools.py
│
├── tests/ smoke / unit
├── .github/ workflows/ + scripts/ (CI、リリース、設計、文書チェック)
└── (repo-top sibling, layer-external bundled forks)
 pysisyphus/ リポジトリ内のオプティマイザ、TS、IRC、COS、calculator のフォーク
 thermoanalysis/ リポジトリ内フォーク
 hessian_ff/ リポジトリ内のネイティブ Hessian/MM 支援、上流 PyPI 配布なし、同梱必須
```

### 2.3 レイヤーごとの責務詳細

**L1 `cli/`** は、ルートディスパッチと共通の argv 解析を受け持ちます。各サブコマンドの Click コマンドは、登録先の `workflows/`、`domain/`、`io/` のモジュールで定義されます。`app.py` はルートの `Click.Group` と `_LAZY_SUBCOMMANDS` レジストリを持ち、各エントリは **絶対モジュールパス** を使います (§5.5)。`preflight.py` (AmberTools / conda env / GPU preflight) がここにあるのは、CLI の起動時、L2 のワークフローより前に実行されるためです。

**L2 `workflows/`** にはコマンドモジュールと共有ワークフローヘルパーがあります。`cli/app.py:_LAZY_SUBCOMMANDS` に登録されたモジュールが `cli` という `@click.command()` を所有します。`_all_helpers.py`、`_opt_freq_common.py`、`_run_session.py`、`scan_common.py`、`restraints.py` などは独立したコマンドを持たない共有ヘルパーです。

**L3 `domain/`**。化学を意識したヘルパーロジックで、`torch` / `numpy` / `pysisyphus.constants` (数値バックエンド) はインポートしてよいですが、MLIP ランタイム (`fairchem`、`orb_models`、`mace`、`aimnet`) は **インポートしてはいけません**。Domain ヘルパーは任意の L2 ステージランナーから再利用できます。

**L4a `backends/`**。ML/MM ONIOM 計算コア (`mlmm_calc.py`) とバックエンドディスパッチ (`__init__.py`) はここにあります。ML 領域の UMA / ORB / MACE / AIMNet2 と OpenMM / hessian_ff の連携はこのレイヤーからディスパッチされます。`mlmm_calc.py` は化学ルール #1、#2、#8 を持ちます (§5.1)。

**L4b `io/`**。出力側の I/O には、ステージごとのサマリーライター、エネルギー図、軌跡レンダリング、PDB/altloc 処理、Hessian キャッシュ、数値 Hessian 構築、および振動数・振動 I/O (`hessian_calc.py`) が含まれます。`io/` は `workflows/` に依存しません。出力形式はここで管理され、ステージランナーから使用されます。

**L5 `core/`**。最下層です。`defaults.py` は共有デフォルトの **唯一の出典** です。別の場所に数値を足す前に、まずここを grep し、そのうえで理由があってコマンドごとに置いたデフォルト値を確かめてください。`utils.py` は共有 PDB / XYZ / プロットヘルパーを保持します。

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

インポートの経路は 2 つあります:

1. **レイヤー化インポートパス**: 外部コードはレイヤーディレクトリから直接インポートします。例: `from mlmm.backends.mlmm_calc import MLMMCore`。
2. **ルートシンボル属性** (`from mlmm import MLMMCore`) — `mlmm/__init__.py:_LAZY_IMPORTS` + PEP 562 `__getattr__` によって処理されます。再エクスポートされる 4 つのシンボル `MLMMCore`、`MLMMASECalculator`、`mlmm`、`mlmm_mm_only` はすべて `mlmm.backends.mlmm_calc` に解決され、初回アクセス時にロードされるため、`import mlmm` は軽いままです (eager なのは `__version__` のみ)。サブモジュールはトップレベルパッケージの属性としてではなく、フルパス (`import mlmm.io.trj2fig`) で到達します。

---

## 3. 初見者向け 5 ステップナビゲーション (合計 ≈ 40 分)

リポジトリを初めて開くコントリビュータは、上から下へこの経路をたどってください。各ステップで 1 つずつ要点を押さえます。

| ステップ | 分 | 開くもの | 分かること |
|------|---------|------|-----------------|
| 1 | 3 | [`README.md`](https://github.com/t-0hmura/mlmm_toolkit/blob/main/README.md) | パッケージを 1 段落で説明した概要 + 単一コマンドの使用法 |
| 2 | 5 | このファイル (`docs/architecture.md`) §2 + §4 | 6 レイヤーのディレクトリツリー、依存方向、各関心事の所在 |
| 3 | 5 | [`mlmm/cli/app.py`](../../mlmm/cli/app.py) | Click ルートグループ、`_LAZY_SUBCOMMANDS` レジストリ (≈ 22 エントリ)、絶対パス解決 |
| 4 | 20 | [`mlmm/workflows/all.py`](../../mlmm/workflows/all.py) (skim) | 1 つの完全なサブコマンドを上から下まで。`extract → mm-parm → ONIOM model → MEP → tsopt → IRC → freq → dft` をトレース |
| 5 | 7 | [`CONTRIBUTING.md`](https://github.com/t-0hmura/mlmm_toolkit/blob/main/CONTRIBUTING.md) §3 + §4 | 5 つの add-a-X レシピ + 触ってはいけない制約 |

ステップ 5 のあとは、§4 のファイルインデックスをたどることで他のファイルも読めます。このパッケージは **各レイヤー内でフラット** で、`mlmm/<layer>/` 配下にネストしたパッケージはありません。主要モジュールは `mlmm/` から 2 ディレクトリ以内にあります。

---

## 4. ファイルインデックス — 「この関心事はどこにある?」

### 4.1 CLI / エントリ (L1 `cli/`)

| 関心事 | ファイル |
|---|---|
| Click ルートグループ + サブコマンドディスパッチ | `mlmm/cli/app.py` |
| サブコマンドリゾルバ (遅延インポート) | `mlmm/cli/default_group.py` |
| `python -m mlmm` shim | `mlmm/__main__.py` |
| 共有オプションデコレータファクトリ | `mlmm/cli/common_options.py` |
| Bool/YAML/例外の CLI ヘルパー | `mlmm/cli/decorators.py` |
| `--help-advanced` pager | `mlmm/cli/help_pages.py` |
| Bool フラグの解釈 (`--flag` / `--no-flag` + 値スタイル) | `mlmm/cli/bool_compat.py` |
| AmberTools / conda env / GPU preflight | `mlmm/cli/preflight.py` |

### 4.2 ワークフローステージランナー (L2 `workflows/`)

以下で用いる略語: GSM = growing-string method、COS = chain-of-states、RS-P-RFO = restricted-step partitioned rational-function optimization、RS-I-RFO = restricted-step image-function rational-function optimization、PHVA = partial Hessian vibrational analysis。

| 関心事 | ファイル |
|---|---|
| 完全パイプラインオーケストレータ | `mlmm/workflows/all.py` |
| 構造最適化 (ONIOM マクロ/マイクロ pre-opt) | `mlmm/workflows/opt.py` |
| Scan と 2D/3D energy-landscape grid + 共有 | `mlmm/workflows/scan{,2d,3d,_common}.py` |
| MEP 探索 (GSM / DMF、再帰的) | `mlmm/workflows/path_search.py` |
| MEP オプティマイザコア (pysisyphus COS) | `mlmm/workflows/path_opt.py` |
| TS 最適化 (RS-P-RFO / RS-I-RFO / TRIM + Bofill + マクロ/マイクロ) | `mlmm/workflows/tsopt.py` |
| 振動解析 (PHVA + MLIP active block) | `mlmm/workflows/freq.py` |
| IRC 積分 (マクロ / マイクロ) | `mlmm/workflows/irc.py` |
| 単一点 DFT (ONIOM 埋め込み) | `mlmm/workflows/dft.py` |
| ML/MM の一点エネルギーと力 | `mlmm/workflows/sp.py` |
| 活性部位抽出 (クラスター切り出し + リンク原子キャップ) | `mlmm/workflows/extract.py` |
| ML / 可動 MM / 固定 MM 領域割り当て | `mlmm/workflows/define_layer.py` |
| AmberTools 駆動の MM パラメータ生成 | `mlmm/workflows/mm_parm.py` |
| ONIOM 入力ライター (Gaussian / ORCA) | `mlmm/workflows/oniom_export.py` |
| ONIOM 入力リーダー (sanity, atom-name diff) | `mlmm/workflows/oniom_import.py` |
| Kabsch / frozen-subset アラインメント | `mlmm/workflows/align_freeze.py` |

### 4.3 化学ヘルパー (L3 `domain/`)

| 関心事 | ファイル |
|---|---|
| R↔P 結合変化検出 | `mlmm/domain/bond_changes.py` |
| Post-IRC 結合サマリー | `mlmm/domain/bond_summary.py` |
| PDB 元素列正規化 | `mlmm/domain/add_elem_info.py` |

### 4.4 MLIP + ONIOM (L4a `backends/`)

| 関心事 | ファイル |
|---|---|
| ML/MM ONIOM 計算コア + 4 つのインライン MLIP バックエンド + ONIOM カップリング | `mlmm/backends/mlmm_calc.py` |
| `--precision` ルーティング (`apply_precision_to_calc_cfg` / `_PRECISION_DISPATCH`) | `mlmm/backends/__init__.py` |
| バックエンドディスパッチ / ファクトリ (`_create_ml_backend`) | `mlmm/backends/mlmm_calc.py` |
| 任意の DFT high-level adapter | `mlmm/backends/pyscf_dft.py` |

[MLIP バックエンド](backends.md) ではインストール方法と実行時の挙動を説明します。
バックエンド実装の変更は、現時点では `mlmm_calc.py` とディスパッチャに反映します。

### 4.5 I/O (L4b `io/`)

| 関心事 | ファイル |
|---|---|
| `summary.json` / `summary.log` ライター | `mlmm/io/summary.py` |
| Plotly エネルギー図 | `mlmm/io/energy_diagram.py` |
| Trajectory → PNG / HTML / SVG / PDF | `mlmm/io/trj2fig.py` |
| PDB altloc 解決 | `mlmm/io/pdb_fix.py` |
| parm7 の原子 index（CHEMISTRY-RULE:9） | `mlmm/io/pdb_indexing.py` |
| インメモリ Hessian キャッシュ (run ごとの TTL) | `mlmm/io/hessian_cache.py` |
| 数値 Hessian 構築 + 振動数 / 振動 I/O | `mlmm/io/hessian_calc.py` |
| 調和拘束のセットアップ | `mlmm/workflows/restraints.py` (L2 ステージヘルパー) |

### 4.6 Foundation (L5 `core/`)

| 関心事 | ファイル |
|---|---|
| 共有ワークフロー・calculator デフォルト | `mlmm/core/defaults.py` |
| DFT 設定（CHEMISTRY-RULE:4） | `mlmm/core/dft_settings.py` |
| PDB / XYZ / プロットヘルパー | `mlmm/core/utils.py` |
| `-v/--verbose LEVEL`（0–3）ロギング配線 | `mlmm/core/logging.py` |
| ステージごとの calc 評価 | `mlmm/core/calc_eval.py` |
| 出力・結果確定ヘルパー | `mlmm/core/output.py`、`mlmm/core/result_commit.py` |
| エネルギー成分の合成 | `mlmm/core/pes_composition.py` |
| 残基テーブル | `mlmm/core/residue_data.py` |

### 4.7 Repo-internal 同梱フォーク

| ディレクトリ | 役割 | 上流と異なるファイル（上流の版で置き換えない） |
|---|---|---|
| `pysisyphus/` | オプティマイザ / TS / IRC エンジン | `irc/IRC.py`、`optimizers/hessian_updates.py`、`tsoptimizers/TSHessianOptimizer.py`、`calculators/*` |
| `thermoanalysis/` | 熱化学 (ΔG, ZPE, 分配関数) | `QCData.py` (upstream とのブランディング差分) |
| `hessian_ff/` | MM 力場上の解析的 Hessian — **PyPI 404、バンドルは必須** | `analytical_hessian.py` (`mlmm/backends/mlmm_calc.py` が消費する唯一のエントリ) |

---

## 5. 科学的な不変条件

### 5.1 9 つの化学ルール (grep レシピ)

正しさに関わる 9 つのルールが `backends/`、`workflows/`、`core/`、`io/` にまたがって存在します。インラインの `# CHEMISTRY-RULE:N` マーカーが実装箇所を示し、`.github/scripts/check_engineering_markers.py` がマーカーの完全性を検査します。

編集前にすべての化学ルールを見つけるには:

```bash
# List all 9 rule sites in the repo (host file + line)
grep -rnE '# CHEMISTRY-RULE:[0-9]+' mlmm/

# List every # DOMAIN_PURE marker (modules the CI check requires to carry it)
grep -rn '# DOMAIN_PURE' mlmm/
```

9 つのルールはすべて `mlmm` に適用されます:

| # | ルール | 実装のファイル |
|---|---|---|
| 1 | Subtractive ONIOM エネルギー式 (`E = mm_real + ml_model − mm_model`) | `mlmm/backends/mlmm_calc.py` |
| 2 | Link-atom Hessian B-matrix 投影 | `mlmm/backends/mlmm_calc.py` |
| 3 | Hessian TS オプティマイザのマクロ / マイクロ交互（RS-P-RFO がデフォルト） | `mlmm/workflows/tsopt.py` |
| 4 | gpu4pyscf `rks_lowmem` の closed-shell/GPU/lowmem guard | `mlmm/core/dft_settings.py` |
| 5 | def2 ファミリーの自動 ECP 注入 | `mlmm/workflows/dft.py` |
| 6 | PHVA + MLIP active-block partial Hessian | `mlmm/workflows/freq.py` |
| 7 | `bofill_update` advanced-indexing scatter | `mlmm/workflows/tsopt.py` |
| 8 | 3-layer 5-pass partial Hessian アセンブリ | `mlmm/backends/mlmm_calc.py` |
| 9 | parm7 アトムインデックス (1-based / serial gap handling) | `mlmm/io/pdb_indexing.py` |

これらの実装を変更する場合は、focused regression test と関連する scheduled numerical test を実行してください（`CONTRIBUTING.md` §1.1）。

**推奨される学習順序 (4 つの化学クラスタ)**:

| クラスタ | ルール | 共通の関心事 | 最初に読むファイル |
|---|---|---|---|
| 5-pass Hessian セット | #1, #2, #8, #9 | subtractive ONIOM + link-atom B-matrix + 3-layer アセンブリ + parm7 インデックス | `mlmm/backends/mlmm_calc.py` (9 ルールのうち 3 つのホスト: #1/#2/#8; #9 は `mlmm/io/pdb_indexing.py`) |
| TS 最適化セット | #3, #7 | マクロ / マイクロ交互 + Bofill scatter | `mlmm/workflows/tsopt.py` |
| 振動セット | #6 | PHVA + MLIP active-block partial Hessian | `mlmm/workflows/freq.py` |
| DFT セット | #4, #5 | gpu4pyscf 低メモリ + def2 ECP 注入 | `mlmm/core/dft_settings.py` (#4)、`mlmm/workflows/dft.py` (#5) |

mlmm における実践的なカリキュラムは、まず 5-pass Hessian セット、次に TS セット (#3, #7)、次に DFT (#4, #5)、最後に振動 (#6) です。

### 5.2 VRAM 管理の不変条件 (`del` チェーンをリファクタしないこと)

IRC / TSopt / Freq ステージは、CUDA メモリを解放するためにステージ間で GPU 常駐オブジェクト (`calc`、`geom`、`hess`) を明示的に `del` します。ステージ境界では `gc.collect()` に加え、CUDA allocation がある場合は `torch.cuda.empty_cache()` も実行します。**これらの解放処理をリファクタで取り除かないでください** — 完全なタンパク質環境での長時間 ML/MM `all` ジョブは、これらがないと OOM します。

### 5.3 同梱フォーク: upstream を併存インストールしないこと

同梱された `pysisyphus/`、`thermoanalysis/`、`hessian_ff/` パッケージは **フォーク** です。`hessian_ff/` には PyPI 版がありません。このパッケージの隣に `pip install pysisyphus` や `pip install thermoanalysis` を再インストールすると、次が静かに壊れます:

- `pysisyphus/irc/IRC.py` — 初期変位のメモリ管理
- `pysisyphus/optimizers/hessian_updates.py` — GPU 常駐の in-place rank-two Bofill 更新、オプトインの `PYSIS_BOFILL_CPU_OFFLOAD=1` フォールバック
- `pysisyphus/tsoptimizers/TSHessianOptimizer.py` — Hessian TS オプティマイザ kwargs
- `pysisyphus/calculators/...` — GPU を意識したバックエンドフック
- `thermoanalysis/QCData.py` — upstream とのブランディング / I/O 差分
- `hessian_ff/analytical_hessian.py` — `backends/mlmm_calc.py` が消費する唯一のエントリ。**upstream の代替は存在しません**

### 5.4 パッケージ検出とランタイム依存

`[tool.setuptools.packages.find].include` は `mlmm*` glob でレイヤーサブパッケージを検出します。新しいトップレベルのパッケージ構成は wheel の内容で確認し、インポートするランタイムパッケージはすべて `dependencies` に宣言します。

### 5.5 `_LAZY_SUBCOMMANDS` レジストリは絶対パスを使用すること

`mlmm/cli/app.py:_LAZY_SUBCOMMANDS` は、すべてのサブコマンドを **絶対** モジュールパスで解決します。相対ドット付きインポート (`".all"` など) は、パッケージルートではなくリゾルバモジュールの `__package__` に解決を依存させます。

---

## 6. 同梱フォーク (repo-internal)

`mlmm_toolkit` はリポジトリのトップに **3 つ** の repo-internal モジュールを同梱します:

| ディレクトリ | 上流の PyPI 版か | 用途 | 許される編集の範囲 |
|---|---|---|---|
| `pysisyphus/` | NO — フォーク、`pip install pysisyphus` を併存させない | オプティマイザ、TS、IRC、COS、calculators | 記載された差分を維持し、数値変更は focused test と scheduled numerical test で検証 |
| `thermoanalysis/` | NO — フォーク (ブランディング差分) | ΔG, ZPE, 分配関数, `QCData` | `QCData` の利用側契約を維持。README 参照 |
| `hessian_ff/` | **NO — PyPI 404、バンドル必須** | MM 力場上の解析的 Hessian | 導関数と公開 API の契約を維持。README 参照 |

各ディレクトリの `README.md` は、上流との差分と利用側の契約を列挙します。


---

## 7. 推奨される深掘り読書順序 (5〜10 ファイル)

初見者向け 5 ステップナビゲーション (§3) のあとは、この深さ優先の読書順序に従ってください:

1. `mlmm/core/defaults.py` — デフォルト値のテーブルです。下流のすべてがここから読み取ります。
2. `mlmm/cli/app.py` — Click ルート + `_LAZY_SUBCOMMANDS` レジストリ。
3. `mlmm/workflows/all.py` — 1 つの完全なパイプラインを上から下まで。
4. `mlmm/workflows/extract.py` + `define_layer.py` — クラスター切り出し + リンク原子キャップ + ONIOM の層の割り当て。
5. `mlmm/workflows/mm_parm.py` — AmberTools parm7 生成。
6. `mlmm/backends/mlmm_calc.py` — ML/MM の心臓部 (CHEMISTRY-RULE:1, 2, 8)。
7. `mlmm/workflows/tsopt.py` — Hessian TS オプティマイザ + Bofill (CHEMISTRY-RULE:7) + マクロ / マイクロ交互 (CHEMISTRY-RULE:3)。
8. `mlmm/workflows/freq.py` — PHVA + MLIP active-block (CHEMISTRY-RULE:6)。
9. `mlmm/workflows/irc.py` — VRAM 管理 + マクロ / マイクロ IRC。
10. `mlmm/core/utils.py` — 共有 PDB / XYZ / プロットヘルパー。

---

## 8. ML/MM (ONIOM) スコープ

`mlmm-toolkit` は ONIOM を介して **完全なタンパク質環境** を扱います:

- **ML 領域**: 基質 + 反応中心残基。4 つの MLIP バックエンド (UMA / ORB / MACE / AIMNet2) のいずれかで評価されます
- **可動 MM 領域**: ML 領域を取り囲むシェルで、AMBER 力場の下で自由に移動できます
- **固定 MM 領域**: タンパク質の残りの部分で、剛体として保持されます

この分割は入力 PDB の B-factor チャネルにエンコードされ、`extract → mm-parm → ONIOM model → MEP → tsopt → IRC → freq → dft` を通じて伝播されます。

## 使用上の注意点

- 今は依存の向きを破る `core/` の import がいくつかあります: `core.utils` → `domain.add_elem_info`・`domain.scan_coordinates`・`io.structure_formats`・`cli.completion`（状態の語彙）、`core.calc_eval` → `backends.mlmm_calc`。どれも循環は作りません。
- `_check_domain_pure` ゲートは、`backends/mlmm_calc.py`、`workflows/tsopt.py`、`workflows/freq.py` に `# DOMAIN_PURE` マーカーがあることだけを確かめます。`workflows/sp.py` にも付いていて、`domain/` のファイルには付いていません。
