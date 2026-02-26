# `opt`

## 概要

> **概要:** L-BFGS（`--opt-mode light`、デフォルト）、RFO（`--opt-mode heavy`）、またはハイブリッド（`--opt-mode hybrid`: LBFGS 後にフラットンループ付き RFO リスタート）を使用して単一構造を局所極小に最適化します。`--flatten` で虚数モードのフラットニングを有効化できます。

`mlmm opt` は、ML/MM 計算機（FAIR-Chem UMA + hessian_ff）を使用して単一構造を局所極小に最適化します。L-BFGS（`--opt-mode light`、デフォルト）、RFO（`--opt-mode heavy`）、またはハイブリッド（`--opt-mode hybrid`: LBFGS 後にフラットンループ付き RFO リスタート）が選択できます。入力は `.pdb`、`.xyz`、`_trj.xyz`、または `geom_loader` がサポートする任意の形式が使用可能です。設定の優先順位は **デフォルト < config < 明示CLI < override** です。

入力が PDB の場合、`--convert-files/--no-convert-files`（デフォルト有効）で制御される `.pdb` コンパニオンファイルも書き出されます。PDB 固有の機能:
- `--freeze-links`（デフォルト `True`）により、リンク水素の親原子が検出され `geom.freeze_atoms`（0 始まりインデックス）にマージされます。
- 出力変換で `final_geometry.pdb`（軌跡ダンプ時は `optimization.pdb`）が入力 PDB をトポロジー参照として生成されます。
- B 因子のアノテーション: ML 領域原子 = 100.00、凍結原子 = 50.00、両方に該当する原子 = 150.00。

### 概要
- **用途:** ML/MM で単一酵素構造を局所エネルギー極小に最小化する場合。
- **手法:** L-BFGS（light、デフォルト）、RFO（heavy）、またはハイブリッド（LBFGS 後に RFO フラットンループ）。ML/MM 計算機は ML 領域に FAIR-Chem UMA、MM に hessian_ff を組み合わせます。
- **出力:** `final_geometry.xyz`、`final_geometry.pdb`（PDB 入力時）、任意の軌跡。
- **デフォルト:** `--opt-mode light`、`--thresh gau`、`--max-cycles 10000`、`--freeze-links` 有効。
- **次のステップ:** [freq](freq.md) を実行して構造が真の極小（虚数振動数なし）であることを確認。

## 最小の例

```bash
mlmm opt -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --out-dir ./result_opt
```

## 出力チェックリスト

- `result_opt/final_geometry.xyz`
- `result_opt/final_geometry.pdb`（入力が PDB で変換が有効な場合）
- `result_opt/optimization_trj.xyz`（`--dump` 有効時）

## 使用例

1. 収束を厳しくして軌跡を保存する。

```bash
mlmm opt -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --thresh gau_tight --dump --out-dir ./result_opt_tight
```

2. 最適化中に 1 本の距離拘束をかける。

```bash
mlmm opt -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --dist-freeze "[(12,45,2.20)]" --bias-k 20.0 --out-dir ./result_opt_rest
```

3. heavy モード（RFO）で最適化する。

```bash
mlmm opt -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --opt-mode heavy --out-dir ./result_opt_rfo
```

4. hybrid モードで虚数モードフラットニングを有効にする。

```bash
mlmm opt -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --opt-mode hybrid --flatten --out-dir ./result_opt_hybrid_flat
```

## ワークフロー

1. **入力処理** -- `-i/--input` には PDB ファイル（酵素複合体）が必要です。オプティマイザーは `pysisyphus.helpers.geom_loader` を介してこの PDB から座標を読み取ります。ML/MM レイヤー定義は `--model-pdb`、`--model-indices`、または `--detect-layer`（B 因子エンコーディング: B=0 ML、B=10 Hessian 対象 MM、B=20 凍結 MM）から取得されます。
2. **ML/MM 計算機の構築** -- ML/MM 計算機（FAIR-Chem UMA + hessian_ff）を構築します。`--real-parm7` で Amber MM トポロジーを提供し、`--model-pdb` で ML 領域を定義します。
3. **リンク凍結検出** -- `--freeze-links`（デフォルト）により、リンク水素の親原子が検出され `geom.freeze_atoms`（0 始まり）にマージされます。マージされたセットは ML/MM 計算機（`calc.freeze_atoms`）に伝播されます。
4. **最適化** -- `--opt-mode light` は L-BFGS、`--opt-mode heavy` は RFOptimizer（RFO）、`--opt-mode hybrid` は L-BFGS 後にフラットンループ付き RFO リスタートを実行します。
   - `--flatten` は最適化後の虚数モードフラットニングを有効にします。検出されたすべての虚数モードが各反復でフラットニングされ、なくなるか内部ループ上限に達するまで続きます。
5. **拘束** -- `--dist-freeze` は Python リテラルタプル `(i, j, target_A)` を消費します。`target_A` は目標距離（オングストローム）で、第 3 要素を省略すると開始距離が拘束されます。`--bias-k` はグローバル調和強度（eV/A^2）を設定します。インデックスはデフォルトで 1 始まりですが、`--zero-based` で 0 始まりに変更可能です。
6. **ダンプと変換** -- `--dump` は `optimization_trj.xyz` を書き出します。変換が有効な場合、PDB 入力では軌跡も `.pdb` にミラーリングされます（B 因子アノテーション付き）。`opt.dump_restart` はリスタート YAML スナップショットを出力できます。
7. **終了コード** -- `0` 成功、`2` ゼロステップ（ステップノルム < `min_step_norm`）、`3` オプティマイザーエラー、`130` キーボード割り込み、`1` 予期しないエラー。

## CLIオプション

> **注意:** 表示されるデフォルト値はオプション未指定時に使用されます。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | `geom_loader` が受け付ける入力構造。 | 必須 |
| `--ref-pdb PATH` | 入力が XYZ の場合の参照 PDB トポロジー。 | _None_ |
| `--real-parm7 PATH` | 全酵素の Amber parm7 トポロジー。 | 必須 |
| `--model-pdb PATH` | ML 領域原子を定義する PDB。`--detect-layer` 有効時は省略可。 | _None_ |
| `--model-indices TEXT` | ML 領域のカンマ区切り原子インデックス（範囲指定可、例: `1-5`）。`--model-pdb` の代替。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` のインデックス規約。 | 1 始まり |
| `--detect-layer / --no-detect-layer` | B 因子（0/10/20）から ML/MM レイヤーを自動検出。 | 有効 |
| `-q, --charge INT` | ML 領域の電荷。 | 必須 |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切りインデックス。 | _None_ |
| `--freeze-links/--no-freeze-links` | リンク水素の親原子凍結の切り替え（PDB 入力のみ）。 | `True` |
| `--radius-partial-hessian FLOAT` | ML 領域からの Hessian-MM 原子の距離カットオフ (A)。`--detect-layer` と併用可。 | _None_ |
| `--radius-freeze FLOAT` | ML 領域からの可動 MM 原子の距離カットオフ (A)。これを超える原子は凍結。 | _None_ |
| `--dist-freeze TEXT` | 調和拘束用の Python リテラル `(i, j, target_A)` タプル。 | _None_ |
| `--one-based / --zero-based` | `--dist-freeze` のインデックス規約。 | 1 始まり |
| `--bias-k FLOAT` | 調和バイアス強度 (eV/A^2)。 | `10.0` |
| `--max-cycles INT` | 最適化反復のハードリミット。 | `10000` |
| `--opt-mode {light\|heavy\|hybrid}` | オプティマイザーモード: `light`（LBFGS）、`heavy`（RFO）、または `hybrid`（LBFGS 後にフラットンループ付き RFO リスタート）。 | `light` |
| `--flatten/--no-flatten` | 最適化後の虚数モードフラットンループの有効化/無効化。 | `False` |
| `--micro-step/--no-micro-step` | `--opt-mode heavy` で `--no-micro-step` は `rfo.max_micro_cycles=1` を強制。 | `True` |
| `--dump/--no-dump` | 軌跡ダンプ（`optimization_trj.xyz`）を出力。 | `False` |
| `--convert-files/--no-convert-files` | PDB 入力時の XYZ/TRJ から PDB コンパニオンの有効化/無効化。 | `True` |
| `--out-dir TEXT` | 出力ディレクトリ。 | `./result_opt/` |
| `--thresh TEXT` | 収束プリセットの上書き（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`、`never`）。 | `gau` |
| `--config FILE` | ベース YAML 設定ファイル。 | _None_ |
| `--show-config/--no-show-config` | 実行前に解決済み YAML レイヤー情報を表示。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに設定検証と実行計画表示のみ行う。 | `False` |

### 収束閾値プリセット

力は Hartree/bohr、ステップは bohr 単位。

| プリセット | 用途 | max\|F\| | RMS(F) | max\|step\| | RMS(step) |
| --- | --- | --- | --- | --- | --- |
| `gau_loose` | 粗い事前最適化、ラフな経路探索 | 2.5e-3 | 1.7e-3 | 1.0e-2 | 6.7e-3 |
| `gau` | 標準的な Gaussian 相当の厳密さ | 4.5e-4 | 3.0e-4 | 1.8e-3 | 1.2e-3 |
| `gau_tight` | より厳密; 良好な構造 / freq / TS 精密化向け | 1.5e-5 | 1.0e-5 | 6.0e-5 | 4.0e-5 |
| `gau_vtight` | 非常に厳密; ベンチマーク/高精度最終構造 | 2.0e-6 | 1.0e-6 | 6.0e-6 | 4.0e-6 |
| `baker` | Baker 式規則（`max\|F\| < 3e-4` **かつ** `\|dE\| < 1e-6 または max\|step\| < 3e-4` で収束） | 3.0e-4 | 2.0e-4 | 3.0e-4 | 2.0e-4 |

## 出力

```
out_dir/ (デフォルト: ./result_opt/)
├─ final_geometry.xyz          # 常に書き出し
├─ final_geometry.pdb          # 入力が PDB で変換有効時のみ（B 因子アノテーション付き）
├─ optimization_trj.xyz        # ダンプ有効時のみ
├─ optimization.pdb            # PDB 入力で変換有効時の軌跡 PDB 変換
└─ restart*.yml                # opt.dump_restart 設定時のオプションリスタート
```

コンソールには解決済みの設定ブロック（`geom`、`calc`、`opt`、`lbfgs`）、`print_every` サイクルごとの進捗、最終的な実行時間サマリーが出力されます。

## YAML設定

設定は **デフォルト < config < 明示CLI < override** の順で適用されます。受け付けるセクション:

### `geom`

- `coord_type`（デフォルト `"cart"`）: デカルト座標 vs `"dlc"` 非局在化内部座標。
- `freeze_atoms`（`[]`）: 最適化中に凍結する 0 始まりインデックス。

### `calc` / `mlmm`

- `input_pdb`、`real_parm7`、`model_pdb`: 必須ファイルパス（文字列）。
- `model_charge`（`-q/--charge`、必須）と `model_mult`（`-m/--multiplicity`、デフォルト 1）。
- `link_mlmm`: ML/MM リンクペアを固定する `(ML_atom_id, MM_atom_id)` 文字列のオプションリスト（リンク原子は作成されません）。
- UMA 制御: `uma_model`（デフォルト `"uma-s-1p1"`）、`uma_task_name`（デフォルト `"omol"`）、`ml_hessian_mode`（`"Analytical"` または `"FiniteDifference"`）、`out_hess_torch`（bool）、`H_double`（bool）。
- デバイス選択: `ml_device`（`"auto"`/`"cuda"`/`"cpu"`）、`ml_cuda_idx`、`mm_device`、`mm_cuda_idx`、`mm_threads`。
- MM 有限差分: `mm_fd`（bool）、`mm_fd_dir`（FD 情報の出力ディレクトリ）、`return_partial_hessian`。
- `return_partial_hessian`: `opt` では YAML で明示指定されない限り部分ヘシアンを既定で使用します。完全ヘシアンを強制する場合は `calc.return_partial_hessian: false` を明示してください。
- `freeze_atoms`: `geom.freeze_atoms` から伝播され、ML/MM とオプティマイザーが同じ凍結原子を共有します。

### `opt`

共有オプティマイザー制御:
- `thresh` プリセット（上記の収束テーブルを参照）。
- 共通制御: `max_cycles`（デフォルト 10000）、`print_every`（100）、`min_step_norm`（1e-8）、`assert_min_step` True。
- 収束トグル: `rms_force`、`rms_force_only`、`max_force_only`、`force_only`。
- その他: `converge_to_geom_rms_thresh`、`overachieve_factor`、`check_eigval_structure`。
- ラインサーチ: `line_search`（bool、デフォルト True）。
- 管理項目: `dump`、`dump_restart`、`prefix`、`out_dir`（デフォルト `./result_opt/`）。

### `lbfgs`

L-BFGS 固有の拡張: `keep_last`、`beta`、`gamma_mult`、`max_step`、`control_step`、`double_damp`、`mu_reg`、`max_mu_reg_adaptions`。

### `rfo`

RFOptimizer 固有の拡張: 信頼領域サイジング（`trust_radius`、`trust_min`、`trust_max`、`trust_update`）、`max_energy_incr`、ヘシアン管理（`hessian_update`、`hessian_init`、`hessian_recalc`、`hessian_recalc_adapt`、`small_eigval_thresh`）、マイクロイテレーション制御（`alpha0`、`max_micro_cycles`、`rfo_overlaps`）、DIIS ヘルパー（`gdiis`、`gediis`、閾値、`gdiis_test_direction`）、`adapt_step_func`。

### YAML 例
```yaml
geom:
 coord_type: cart               # 座標タイプ: デカルト vs dlc 内部座標
 freeze_atoms: []               # 0 始まり凍結原子（CLI/リンク検出とマージ）
calc:
 charge: 0                      # 総電荷（CLI 上書き）
 spin: 1                        # スピン多重度 2S+1
mlmm:
 real_parm7: real.parm7         # 全酵素の Amber parm7 トポロジー
 model_pdb: ml_region.pdb       # ML 領域を定義する PDB
 uma_model: uma-s-1p1           # UMA モデルタグ
 uma_task_name: omol             # UMA タスク名
 ml_device: auto                # UMA デバイス選択
 ml_hessian_mode: FiniteDifference  # ヘシアンモード選択
 out_hess_torch: true           # torch 形式ヘシアンを要求
 mm_fd: false                   # MM 有限差分トグル
 return_partial_hessian: true   # 部分ヘシアンを許可（opt のデフォルト）
opt:
 thresh: gau                    # 収束プリセット（Gaussian/Baker 式）
 max_cycles: 10000              # オプティマイザーサイクル上限
 print_every: 100               # ログ出力間隔
 min_step_norm: 1.0e-08         # ステップ受け入れの最小ノルム
 assert_min_step: true          # ステップが閾値以下で停止
 rms_force: null                # 明示的 RMS 力目標
 rms_force_only: false          # RMS 力収束のみに依存
 max_force_only: false          # 最大力収束のみに依存
 force_only: false              # 変位チェックをスキップ
 converge_to_geom_rms_thresh: 0.05  # 参照への収束時の geom RMS 閾値
 overachieve_factor: 0.0        # 閾値を厳しくする係数
 check_eigval_structure: false  # ヘシアン固有値構造の検証
 line_search: true              # ラインサーチを有効化
 dump: false                    # 軌跡/リスタートデータのダンプ
 dump_restart: false            # リスタートチェックポイントのダンプ
 prefix: ""                     # ファイル名プレフィックス
 out_dir: ./result_opt/         # 出力ディレクトリ
lbfgs:
 thresh: gau                    # LBFGS 収束プリセット
 max_cycles: 10000              # 反復上限
 print_every: 100               # ログ出力間隔
 min_step_norm: 1.0e-08         # 受け入れ最小ステップノルム
 assert_min_step: true          # ステップ停滞時にアサート
 rms_force: null                # 明示的 RMS 力目標
 rms_force_only: false          # RMS 力収束のみに依存
 max_force_only: false          # 最大力収束のみに依存
 force_only: false              # 変位チェックをスキップ
 converge_to_geom_rms_thresh: 0.05  # ジオメトリ収束時の RMS 閾値
 overachieve_factor: 0.0        # 閾値を厳しくする
 check_eigval_structure: false  # ヘシアン固有値構造の検証
 line_search: true              # ラインサーチを有効化
 dump: false                    # 軌跡/リスタートデータのダンプ
 dump_restart: false            # リスタートチェックポイントのダンプ
 prefix: ""                     # ファイル名プレフィックス
 out_dir: ./result_opt/         # 出力ディレクトリ
 keep_last: 7                   # LBFGS バッファの履歴サイズ
 beta: 1.0                      # 初期ダンピングベータ
 gamma_mult: false              # 乗算ガンマ更新トグル
 max_step: 0.3                  # 最大ステップ長
 control_step: true             # 適応的ステップ長制御
 double_damp: true              # ダブルダンピングセーフガード
 mu_reg: null                   # 正則化強度
 max_mu_reg_adaptions: 10       # mu 適応の上限
rfo:
 thresh: gau                    # RFOptimizer 収束プリセット
 max_cycles: 10000              # 反復上限
 print_every: 100               # ログ出力間隔
 min_step_norm: 1.0e-08         # 受け入れ最小ステップノルム
 assert_min_step: true          # ステップ停滞時にアサート
 rms_force: null                # 明示的 RMS 力目標
 rms_force_only: false          # RMS 力収束のみに依存
 max_force_only: false          # 最大力収束のみに依存
 force_only: false              # 変位チェックをスキップ
 converge_to_geom_rms_thresh: 0.05  # ジオメトリ収束時の RMS 閾値
 overachieve_factor: 0.0        # 閾値を厳しくする
 check_eigval_structure: false  # ヘシアン固有値構造の検証
 line_search: true              # ラインサーチを有効化
 dump: false                    # 軌跡/リスタートデータのダンプ
 dump_restart: false            # リスタートチェックポイントのダンプ
 prefix: ""                     # ファイル名プレフィックス
 out_dir: ./result_opt/         # 出力ディレクトリ
 trust_radius: 0.1              # 信頼領域半径
 trust_update: true             # 信頼領域更新を有効化
 trust_min: 0.0                 # 最小信頼半径
 trust_max: 0.1                 # 最大信頼半径
 max_energy_incr: null          # ステップごとの許容エネルギー増加
 hessian_update: bfgs           # ヘシアン更新方式
 hessian_init: calc             # ヘシアン初期化ソース
 hessian_recalc: 200            # N ステップごとにヘシアンを再構築
 hessian_recalc_adapt: null     # 適応的ヘシアン再構築上限
 small_eigval_thresh: 1.0e-08   # 安定性のための固有値閾値
 alpha0: 1.0                    # 初期マイクロステップ
 max_micro_cycles: 50           # マイクロイテレーション上限
 rfo_overlaps: false            # RFO オーバーラップを有効化
 gediis: false                  # GEDIIS を有効化
 gdiis: true                    # GDIIS を有効化
 gdiis_thresh: 0.0025           # GDIIS 受け入れ閾値
 gediis_thresh: 0.01            # GEDIIS 受け入れ閾値
 gdiis_test_direction: true     # DIIS 前に降下方向をテスト
 adapt_step_func: true          # 適応的ステップスケーリング
```

## 注意事項

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- 電荷/多重度の運用ルールは [CLI Conventions](cli_conventions.md) に集約しています。
- **デバイス:** `ml_device="auto"` は CUDA が利用可能な場合に選択します。`mm_device` は MM バックエンドのデバイス配置を制御します。
- **ヘシアン:** `calc.out_hess_torch=True` は PyTorch テンソルを返します（`calc.H_double` で任意に倍精度）。
- **凍結原子:** CLI のリンク凍結ロジックは YAML `geom.freeze_atoms` とマージされ、ML/MM 計算機（`calc.freeze_atoms`）に伝播されます。
- **終了コード:** `ZeroStepLength` -> 終了コード **2**、`OptimizationError` -> **3**、`KeyboardInterrupt` -> **130**、その他の未処理例外 -> **1**。
- **優先順位:** 設定は **デフォルト < config < 明示CLI < override** の順序で適用されます。

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細なトラブルシューティングガイド

- [tsopt](tsopt.md) -- 極小ではなく遷移状態（鞍点）を最適化
- [freq](freq.md) -- 最適化が極小に達したことを確認する振動解析
- [all](all.md) -- 端点を事前最適化するエンドツーエンドワークフロー
- [YAML リファレンス](yaml_reference.md) -- `opt`、`lbfgs`、`rfo` の完全な設定オプション
- [用語集](glossary.md) -- L-BFGS、RFO の定義
