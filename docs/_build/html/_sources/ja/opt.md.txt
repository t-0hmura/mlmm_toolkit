# `opt`

## 概要

> **要約:** ML/MM 計算機を使用して L-BFGS（デフォルト）または RFO で単一構造を局所極小に最適化します。入力は PDB ファイルが必須で、B 因子アノテーション付きの XYZ および PDB ジオメトリを出力します。

`mlmm opt` は、ML/MM 計算機（`mlmm_toolkit.mlmm_calc.mlmm`）を使用して PySisyphus LBFGS による単一構造の構造最適化を実行します。計算機は FAIR-Chem UMA（ML 高レイヤー）と OpenMM（MM 低レイヤー）をリンク原子なしで結合します。ML 領域は `--model-pdb` で定義されます。設定は YAML セクション `geom`、`calc`（エイリアス `mlmm`）、`opt`、`lbfgs` を使用します。優先順位は **内部デフォルト < config < 明示CLI < override** です。

## 最小例

```bash
mlmm opt -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --out-dir ./result_opt
```

## 出力の見方

- `result_opt/final_geometry.xyz`
- `result_opt/final_geometry.pdb`（入力が PDB の場合）
- `result_opt/optimization_trj.xyz`（`--dump` 有効時）

## よくある例

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

## 使用法

```bash
mlmm opt -i INPUT.pdb --real-parm7 real.parm7 --model-pdb model.pdb -q CHARGE [-m MULT]
 [--dist-freeze "[(I,J,TARGET_A),...]"] [--one-based|--zero-based] [--bias-k FLOAT]
 [--freeze-atoms "1,3,5"] [--max-cycles N] [--thresh PRESET]
```

### 例

```bash
mlmm opt -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0

mlmm opt -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 \
 --freeze-atoms "1,3,5,7" --thresh gau_tight --dump --out-dir ./result_opt/ \
```

## ワークフロー

- **入力処理:** `-i/--input` には PDB ファイル（酵素複合体）が必要です。オプティマイザーは `pysisyphus.helpers.geom_loader` を介してこの PDB から座標を読み取ります。
- **調和距離拘束** は `--dist-freeze` で利用でき、強度は `--bias-k`（eV/Angstrom^2）で設定します。
- **PDB 変換に関する注記:** XYZ/TRJ から PDB への変換後、B 因子は次のようにアノテーションされます: ML 領域原子 = 100.00、凍結原子 = 50.00、両方に該当する原子 = 150.00。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力 PDB ファイル（酵素複合体）。 | 必須 |
| `--real-parm7 PATH` | 全酵素の Amber parm7 トポロジー。 | 必須 |
| `--model-pdb PATH` | ML 領域原子を定義する PDB。`--detect-layer` 有効時は省略可。 | _None_ |
| `-q, --charge INT` | ML 領域の電荷。 | 必須 |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `--dist-freeze TEXT` | 調和拘束用の Python リテラル `(i, j, target_A)` タプル。 | _None_ |
| `--one-based / --zero-based` | `--dist-freeze` のインデックス規約。 | 1 始まり |
| `--bias-k FLOAT` | 調和バイアス強度 (eV/Angstrom^2)。 | _デフォルト_ |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切りインデックス。 | _None_ |
| `--max-cycles INT` | 最適化反復のハードリミット。 | `10000` |
| `--thresh TEXT` | 収束プリセット（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`）。 | _デフォルト_ |
| `--dump/--no-dump` | 軌跡ダンプ（`optimization_trj.xyz`）を出力。 | `False` |
| `--out-dir PATH` | 出力ディレクトリ。 | `./result_opt/` |
| `--config FILE` | ベースYAML設定ファイル。 | _None_ |
| `--show-config/--no-show-config` | 実行前に解決済みYAMLレイヤ情報を表示。 | `False` |
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

```text
out_dir/ (デフォルト:./result_opt/)
 final_geometry.xyz # 最適化されたジオメトリ（常に書き出し）
 final_geometry.pdb # 入力が PDB の場合に XYZ から変換（B 因子アノテーション付き）
 optimization_trj.xyz # 軌跡（--dump または opt.dump: true の場合）
 optimization.pdb # 入力が PDB でダンプ有効時に TRJ から変換
 restart*.yml # opt.dump_restart サイクルごとのリスタートファイル（有効時）
```

コンソールには解決済みの設定ブロック（`geom`、`calc`、`opt`、`lbfgs`）、`print_every` サイクルごとの進捗、最終的な実行時間サマリーが出力されます。

## YAML 設定

設定は **内部デフォルト < config < 明示CLI < override** の順で適用されます。受け付けるセクション:

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

## 注意事項

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- 電荷/多重度の運用ルールは [CLI Conventions](cli_conventions.md) に集約しています。
- **デバイス:** `ml_device="auto"` は CUDA が利用可能な場合に選択します。`mm_device` は OpenMM の配置を制御します。
- **ヘシアン:** `calc.out_hess_torch=True` は PyTorch テンソルを返します（`calc.H_double` で任意に倍精度）。
- **終了コード:** `ZeroStepLength` -> 終了コード **2**、`OptimizationError` -> **3**、`KeyboardInterrupt` -> **130**、その他の未処理例外 -> **1**。
- **優先順位:** 設定は **内部デフォルト < config < 明示CLI < override** の順序で適用されます。

---

## 関連項目

- [tsopt](tsopt.md) -- 極小ではなく遷移状態（鞍点）を最適化
- [freq](freq.md) -- 最適化が極小に達したことを確認する振動解析
- [all](all.md) -- 端点を事前最適化するエンドツーエンドワークフロー
