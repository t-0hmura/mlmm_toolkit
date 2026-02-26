# `path-search`

## 概要

> **要約:** 再帰的 GSM セグメンテーションにより 2 つ以上の構造から連続した MEP を構築します。結合変化のある領域のみを自動的に精密化し、最高エネルギーイメージ（HEI）を TS 候補としてエクスポートします。

### 概要
- **用途:** R -> ... -> P の構造（2 つ以上の入力）があり、自動精密化付きの単一の連結 MEP が必要な場合に使用。
- **手法:** GSM セグメントを連鎖し、共有結合変化を含むサブ区間のみを再帰的に精密化。
- **出力:** `mep_trj.xyz`（メイン軌跡）、`summary.yaml`（セグメントごとの結果）、有効時にプロット/マージ PDB。
- **デフォルト:** `--opt-mode light`（LBFGS）、`--preopt`、`--align`、`--thresh gau`。
- **次のステップ:** HEI 出力だけでは TS を検証できません。[tsopt](tsopt.md)、[freq](freq.md)、[irc](irc.md) で続行してください。

`mlmm path-search` は GSM を使用して 2 つ以上の構造にわたる連続した最小エネルギー経路（MEP）を構築します。共有結合変化が検出された領域のみを選択的に精密化し、解決されたサブパスを 1 つの軌跡に統合します。

端点が **2 つ**だけで再帰的精密化が不要な場合は、[path-opt](path_opt.md) がより簡単な選択です。

## 最小例

```bash
mlmm path-search -i reactant.pdb product.pdb --real-parm7 real.parm7 \
 --model-pdb ml_region.pdb -q 0 --out-dir ./result_path_search
```

## 出力チェックリスト

- `result_path_search/mep_trj.xyz`
- `result_path_search/summary.yaml`
- `result_path_search/summary.log`
- `result_path_search/mep_plot.png`（プロット生成時）

## よくある例

1. 中間体を含む多段経路を構築する。

```bash
mlmm path-search -i R.pdb IM1.pdb IM2.pdb P.pdb --real-parm7 real.parm7 \
 --model-pdb ml_region.pdb -q -1 --out-dir ./result_path_search_multi
```

2. ポケット軌跡をフルテンプレートへマージする。

```bash
mlmm path-search -i R.pdb IM1.pdb P.pdb --real-parm7 real.parm7 \
 --model-pdb ml_region.pdb -q 0 --ref-pdb holo_template.pdb \
 --out-dir ./result_path_search_merge
```

3. 事前最適化とアライメントを無効にして軽く試す。

```bash
mlmm path-search -i reactant.pdb product.pdb --real-parm7 real.parm7 \
 --model-pdb ml_region.pdb -q 0 --no-preopt --no-align --max-nodes 8 \
 --out-dir ./result_path_search_fast
```

## 使用法

```bash
mlmm path-search -i R.pdb IM1.pdb P.pdb \
 --real-parm7 real.parm7 --model-pdb ml_region.pdb -q CHARGE [-m MULT]
 [--mep-mode gsm|dmf] [--refine-mode peak|minima]
 [--freeze-atoms "1,3,5"] [--max-nodes N] [--max-cycles N] [--climb/--no-climb]
 [--opt-mode light|heavy]
 [--thresh PRESET] [--dump/--no-dump] [--out-dir DIR]
 [--show-config/--no-show-config] [--dry-run/--no-dry-run]
```

### 例

```bash
# ミニマルなポケットのみの 2 状態間 MEP
mlmm path-search -i reactant.pdb product.pdb --real-parm7 real.parm7 \
 --model-pdb ml_region.pdb -q 0

# YAML 上書き、凍結原子、全系マージ出力付きマルチステップ経路
mlmm path-search -i R.pdb IM1.pdb P.pdb --real-parm7 real.parm7 \
 --model-pdb ml_region.pdb -q -1 --freeze-atoms "1,3,5" \
 --ref-pdb holo_template.pdb --out-dir ./run_ps
```

## ワークフロー

1. **初期セグメント（隣接ペア A->B ごと; GSM/DMF）** -- 選択した MEP エンジン（`--mep-mode`）を実行して粗い MEP を取得し、最高エネルギーイメージ（HEI）を特定。
2. **HEI 周辺の局所緩和** -- `--refine-mode`（`peak`: HEI+/-1、`minima`: 最近傍局所極小）で種点を選び、単一構造オプティマイザー（`opt-mode`）で精密化して近傍の極小（`End1`、`End2`）を回復。
3. **キンク vs 精密化の判定**:
 - `End1` と `End2` の間に共有結合変化が検出されない場合、その領域を*キンク*として扱い、`search.kink_max_nodes` 個の線形ノードを挿入して各ノードを個別に最適化。
 - それ以外の場合、`End1` と `End2` の間で**精密化セグメント（GSM）**を起動して障壁を鮮明化。
4. **選択的再帰** -- `(A->End1)` と `(End2->B)` の結合変化を `bond` 閾値で比較。共有結合の更新を含むサブ区間のみに再帰。再帰深度は `search.max_depth` で制限。
5. **統合とブリッジ** -- 解決されたサブパスを連結し、RMSD <= `search.stitch_rmsd_thresh` の重複端点を削除。2 つの統合部分の間の RMSD ギャップが `search.bridge_rmsd_thresh` を超える場合、選択中の `--mep-mode` でブリッジ MEP セグメントを挿入。インターフェース自体に結合変化がある場合、ブリッジの代わりに新たな再帰セグメントを生成。
6. **任意のアライメントとマージ** -- 事前最適化後、`--align` で入力を剛体アラインし凍結を精密化。`--ref-pdb` がある場合、ポケット軌跡を完全テンプレートにマージし、セグメントをプロット/分析用にアノテーション。

結合変化検出は `bond` YAML セクションの閾値を使用する `bond_changes.compare_structures` に依存します。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順の 2 つ以上の完全酵素 PDB。`-i` を繰り返すか、1 つのフラグの後に複数パスを渡す。 | 必須 |
| `--real-parm7 PATH` | 完全酵素複合体の Amber parm7 トポロジー。 | 必須 |
| `--model-pdb PATH` | ML/MM の ML（高レベル）領域原子を定義する PDB。`--detect-layer` または `--model-indices` 利用時は省略可。 | _None_ |
| `-q, --charge INT` | ML 領域の電荷（整数）。 | 必須 |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `--mep-mode [gsm\|dmf]` | セグメント/ブリッジ探索に使う MEP バックエンド。 | `gsm` |
| `--refine-mode [peak\|minima]` | HEI 精密化の種点ルール。 | `gsm` は `peak`、`dmf` は `minima` |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切りインデックス（YAML `geom.freeze_atoms` とマージ）。 | _None_ |
| `--max-nodes INT` | セグメント GSM の内部ノード数。 | `10` |
| `--max-cycles INT` | GSM マクロサイクルの最大数。 | `300` |
| `--climb/--no-climb` | セグメント GSM の TS 精密化を有効化。 | `True` |
| `--opt-mode [light\|heavy]` | 単一構造オプティマイザープリセット（`light` = LBFGS、`heavy` = RFO）。 | `light` |
| `--preopt/--no-preopt` | セグメンテーション前に端点を LBFGS で事前最適化。 | `True` |
| `--align / --no-align` | 事前最適化後に入力を剛体アライメント。 | 有効 |
| `--thresh TEXT` | 収束プリセット（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`）。 | _デフォルト_ |
| `--dump/--no-dump` | オプティマイザーダンプを保存。 | `False` |
| `--out-dir PATH` | 出力ディレクトリ。 | `./result_path_search/` |
| `--ref-pdb PATH...` | 最終マージ用の完全テンプレート PDB。 | _None_ |
| `--config FILE` | 明示 CLI 指定より前に適用されるベース YAML。 | _None_ |
| `--show-config/--no-show-config` | 解決済み設定（YAML レイヤ情報を含む）を表示して実行継続。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画表示のみを行う。 | `False` |

## 出力

```text
out_dir/ (デフォルト:./result_path_search/)
 summary.yaml # MEP レベルの実行サマリー（完全設定ダンプなし）
 summary.log # 人間が読めるサマリー
 mep_trj.xyz # 最終 MEP（常に書き出し）
 mep.pdb # 最終 MEP（参照テンプレート利用可能時は PDB）
 mep_w_ref.pdb # 全系マージ MEP（--ref-pdb 必要）
 mep_w_ref_seg_XX.pdb # セグメントごとのマージ MEP（結合変化セグメント; --ref-pdb 必要）
 mep_seg_XX_trj.xyz / mep_seg_XX.pdb # ポケットのみのセグメント別経路
 hei_seg_XX.xyz / hei_seg_XX.pdb # ポケット HEI と結合変化セグメントごとのオプション PDB
 hei_w_ref_seg_XX.pdb # 結合変化セグメントごとのマージ HEI（--ref-pdb 必要）
 mep_plot.png # イメージインデックスに対する Delta-E プロファイル（trj2fig より）
 energy_diagram_MEP.png # 反応物基準の状態レベルエネルギーダイアグラム（kcal/mol）
 seg_000_*/ # セグメントレベルの GSM と精密化成果物
```

## YAML 設定

マージ順は **defaults < config < 明示指定 CLI < override** です。

YAML ルートはマッピングでなければなりません。受け付けるセクション:

- **`geom`** -- `coord_type`（デフォルト `"cart"`）、`freeze_atoms`（0 始まりインデックス）。
- **`calc` / `mlmm`** -- ML/MM 計算機設定: `input_pdb`、`real_parm7`、`model_pdb`、`model_charge`、`model_mult`、UMA 制御（`uma_model`、`uma_task_name`、`ml_hessian_mode`）、デバイス選択、凍結原子。
- **`gs`** -- Growing String 設定: `max_nodes`、`climb`、`climb_rms`、`climb_fixed`、`reparam_every_full`、`reparam_check`。
- **`opt`** -- StringOptimizer 制御: `max_cycles`、`print_every`、`dump`、`dump_restart`、`out_dir`。
- **`lbfgs`** -- HEI+/-1 精密化用の単一構造オプティマイザー制御: `keep_last`、`beta`、`gamma_mult`、`max_step`、`control_step`、`double_damp`、`mu_reg`、`max_mu_reg_adaptions`。
- **`bond`** -- 結合変化検出: `bond_factor`、`margin_fraction`、`delta_fraction`。
- **`search`** -- 再帰ロジック: `max_depth`、`stitch_rmsd_thresh`、`bridge_rmsd_thresh`、`max_nodes_segment`、`max_nodes_bridge`、`kink_max_nodes`、`max_seq_kink`、`refine_mode`。

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- 入力: `-i/--input` に反応順の完全酵素 PDB を少なくとも 2 つ提供してください。
- preflight チェックで `-i/--input` と `--ref-pdb` のファイル存在を実行前に検証します。
- 電荷/多重度の運用ルールは [CLI Conventions](cli_conventions.md) に集約しています。
- 凍結原子: `--freeze-atoms "1,3,5"` は 0 始まりインデックスとして保存され、YAML `geom.freeze_atoms` とマージされます。
- ノードと再帰: セグメント vs ブリッジのノードは `search.max_nodes_segment` と `search.max_nodes_bridge` で異なります。キンクは `search.kink_max_nodes`（デフォルト 3）の線形ノードを使用します。再帰深度は `search.max_depth`（デフォルト 10）で制限されます。
- オプティマイザー: `--mep-mode gsm` は pysisyphus `GrowingString` + `StringOptimizer`、`--mep-mode dmf` は Direct Max Flux を使用します。単一構造精密化は常に LBFGS です。
- `--align` での最終マージ規則: `--ref-pdb` が提供された場合、最初の参照 PDB がすべてのペアに使用されます。
- コンソールには状態シーケンス（例: `R --> TS1 --> IM1 -->... --> P`）とエネルギーダイアグラム構築に使用されるラベル/エネルギーが出力されます。
- `summary.log` は一部のフィールド欠落時でも生成可能です。内部で次のキーに既定値を補います:
 `root_out_dir`, `path_module_dir`, `pipeline_mode`, `segments`, `energy_diagrams`。

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け

- [path-opt](path_opt.md) -- シングルパス MEP 最適化（再帰的精密化なし）
- [opt](opt.md) -- 単一構造の構造最適化
- [all](all.md) -- 内部で path-search を呼び出すエンドツーエンドワークフロー
- [trj2fig](trj2fig.md) -- MEP 軌跡からエネルギープロファイルをプロット
