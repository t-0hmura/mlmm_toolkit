# `path-search`

## 概要

> **要約:** 再帰的 GSM セグメンテーションにより 2 つ以上の構造から連続した MEP を構築します。結合変化のある領域のみを自動的に精密化し、最高エネルギーイメージ（HEI）を TS 候補としてエクスポートします。

`mlmm path-search` は、反応に沿って順序付けられた 2 つ以上の構造間の連続した最小エネルギー経路（MEP）を構築します。隣接する各ペアは ML/MM 計算機（`mlmm_toolkit.mlmm_calc.mlmm`）を使用した Growing String Method（GSM）で処理されます。ML/MM 計算機は FAIR-Chem UMA と OpenMM をリンク原子なしで結合します。HEI+/-1 イメージは LBFGS で精密化され、共有結合変化が分析され、結合変化を示す領域のみが再帰します。キンク領域は線形補間と単一構造最適化に依存します。マルチ構造入力は、必要に応じて RMSD ブリッジングにより単一の MEP に統合されます。設定の優先順位は **defaults < config < 明示指定 CLI < override**（`--config` → `--override-yaml`）で、`--args-yaml` は legacy alias です。

## 最小例

```bash
mlmm path-search -i reactant.pdb product.pdb --real-parm7 real.parm7 \
  --model-pdb ml_region.pdb -q 0 --out-dir ./result_path_search
```

## 出力の見方

- `result_path_search/summary.md`
- `result_path_search/key_mep.trj` / `result_path_search/key_ts.xyz`（利用可能時）
- `result_path_search/mep.trj`
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
  --model-pdb ml_region.pdb -q 0 --no-pre-opt --no-align --max-nodes 8 \
  --out-dir ./result_path_search_fast
```

## 使用法

```bash
mlmm path-search -i R.pdb IM1.pdb P.pdb \
    --real-parm7 real.parm7 --model-pdb ml_region.pdb -q CHARGE [-m MULT]
    [--freeze-atoms "1,3,5"] [--max-nodes N] [--max-cycles N] [--climb/--no-climb]
    [--thresh PRESET] [--dump/--no-dump] [--out-dir DIR]
    [--config FILE] [--override-yaml FILE] [--args-yaml FILE]
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
    --config base.yaml --override-yaml override.yaml \
    --ref-pdb holo_template.pdb --out-dir ./run_ps
```

## ワークフロー

1. **初期セグメント（隣接ペア A->B ごと）** -- ML/MM で GSM を実行して予備 MEP を取得。
2. **障壁の局所化** -- 最高エネルギーイメージ（HEI）を検出し、LBFGS で HEI+/-1 を最適化して End1 と End2 を取得。
3. **精密化** -- End1-End2 に共有結合変化がない場合（キンク）、`search.kink_max_nodes` 個の線形イメージを挿入して各イメージを最適化。それ以外の場合、End1 と End2 間で精密化 GSM を実行。
4. **選択的再帰** -- (A->End1) と (End2->B) の共有結合変化を評価し、変化のある側のみ再帰。
5. **サブパスの統合** -- RMSD による重複除去でサブ MEP を連結。端点が `search.bridge_rmsd_thresh` を超えて不一致の場合はブリッジ GSM を挿入。共有結合変化のあるインターフェースはブリッジではなく再帰セグメントを生成。
6. **任意のアライメントとマージ** -- 事前最適化後、`--align` で入力を剛体アラインし凍結を精密化。`--ref-pdb` がある場合、ポケット軌跡を完全テンプレートにマージし、セグメントをプロット/分析用にアノテーション。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順の 2 つ以上の完全酵素 PDB。`-i` を繰り返すか、1 つのフラグの後に複数パスを渡す。 | 必須 |
| `--real-parm7 PATH` | 完全酵素複合体の Amber parm7 トポロジー。 | 必須 |
| `--model-pdb PATH` | ML/MM の ML（高レベル）領域原子を定義する PDB。`--detect-layer` または `--model-indices` 利用時は省略可。 | _None_ |
| `-q, --charge INT` | ML 領域の電荷（整数）。省略時は `0`（YAML で上書き可）。 | `0` |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切りインデックス（YAML `geom.freeze_atoms` とマージ）。 | _None_ |
| `--max-nodes INT` | セグメント GSM の内部ノード数。 | `10` |
| `--max-cycles INT` | GSM マクロサイクルの最大数。 | `300` |
| `--climb/--no-climb` | セグメント GSM の TS 精密化を有効化。 | `True` |
| `--pre-opt/--no-pre-opt` | セグメンテーション前に端点を LBFGS で事前最適化。 | `True` |
| `--align / --no-align` | 事前最適化後に入力を剛体アライメント。 | 有効 |
| `--thresh TEXT` | 収束プリセット（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`）。 | _デフォルト_ |
| `--dump/--no-dump` | オプティマイザーダンプを保存。 | `False` |
| `--out-dir PATH` | 出力ディレクトリ。 | `./result_path_search/` |
| `--ref-pdb PATH...` | 最終マージ用の完全テンプレート PDB。 | _None_ |
| `--config FILE` | 明示 CLI 指定より前に適用されるベース YAML。 | _None_ |
| `--override-yaml FILE` | 最後に適用される上書き YAML（YAML 最優先レイヤ）。 | _None_ |
| `--args-yaml FILE` | `--override-yaml` の legacy alias。 | _None_ |
| `--show-config/--no-show-config` | 解決済み設定（YAML レイヤ情報を含む）を表示して実行継続。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画表示のみを行う。 | `False` |

## 出力

```text
out_dir/ (デフォルト: ./result_path_search/)
  summary.yaml                  # MEP レベルの実行サマリー（完全設定ダンプなし）
  summary.log                   # 人間が読めるサマリー
  summary.md                    # 主要成果物へ移動しやすいナビゲーションページ
  key_mep.trj                   # 主要 MEP 軌跡へのショートカット（symlink/copy）
  key_mep.pdb                   # 主要 MEP PDB へのショートカット（symlink/copy）
  key_ts.xyz / key_ts.pdb       # TS 候補スナップショットへのショートカット（利用可能時）
  key_mep_plot.png              # MEP プロファイルへのショートカット（利用可能時）
  key_energy_diagram_MEP.png    # 状態エネルギーダイアグラムへのショートカット（利用可能時）
  mep.trj                       # 最終 MEP（常に書き出し）
  mep.pdb                       # 最終 MEP（参照テンプレート利用可能時は PDB）
  mep_w_ref.pdb                 # 全系マージ MEP（--ref-pdb 必要）
  mep_w_ref_seg_XX.pdb          # セグメントごとのマージ MEP（結合変化セグメント; --ref-pdb 必要）
  mep_seg_XX.trj / mep_seg_XX.pdb  # ポケットのみのセグメント別経路
  hei_seg_XX.xyz / hei_seg_XX.pdb  # ポケット HEI と結合変化セグメントごとのオプション PDB
  hei_w_ref_seg_XX.pdb          # 結合変化セグメントごとのマージ HEI（--ref-pdb 必要）
  mep_plot.png                  # イメージインデックスに対する Delta-E プロファイル（trj2fig より）
  energy_diagram_MEP.png        # 反応物基準の状態レベルエネルギーダイアグラム（kcal/mol）
  seg_000_*/                    # セグメントレベルの GSM と精密化成果物
```

## YAML 設定（`--config`, `--override-yaml`, `--args-yaml`）

マージ順は **defaults < config < 明示指定 CLI < override** です。
`--args-yaml` は `--override-yaml` の legacy alias として維持されています。
YAML ルートはマッピングでなければなりません。受け付けるセクション:

- **`geom`** -- `coord_type`（デフォルト `"cart"`）、`freeze_atoms`（0 始まりインデックス）。
- **`calc` / `mlmm`** -- ML/MM 計算機設定: `input_pdb`、`real_parm7`、`model_pdb`、`model_charge`、`model_mult`、UMA 制御（`uma_model`、`uma_task_name`、`ml_hessian_mode`）、デバイス選択、凍結原子。
- **`gs`** -- Growing String 設定: `max_nodes`、`climb`、`climb_rms`、`climb_fixed`、`reparam_every_full`、`reparam_check`。
- **`opt`** -- StringOptimizer 制御: `max_cycles`、`print_every`、`dump`、`dump_restart`、`out_dir`。
- **`lbfgs`** -- HEI+/-1 精密化用の単一構造オプティマイザー制御: `keep_last`、`beta`、`gamma_mult`、`max_step`、`control_step`、`double_damp`、`mu_reg`、`max_mu_reg_adaptions`。
- **`bond`** -- 結合変化検出: `bond_factor`、`margin_fraction`、`delta_fraction`。
- **`search`** -- 再帰ロジック: `max_depth`、`stitch_rmsd_thresh`、`bridge_rmsd_thresh`、`max_nodes_segment`、`max_nodes_bridge`、`kink_max_nodes`、`max_seq_kink`。

## 注意事項

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- 入力: `-i/--input` に反応順の完全酵素 PDB を少なくとも 2 つ提供してください。
- preflight チェックで `-i/--input` と `--ref-pdb` のファイル存在を実行前に検証します。
- 電荷/多重度の運用ルールは [CLI Conventions](cli-conventions.md) に集約しています。
- 凍結原子: `--freeze-atoms "1,3,5"` は 0 始まりインデックスとして保存され、YAML `geom.freeze_atoms` とマージされます。
- ノードと再帰: セグメント vs ブリッジのノードは `search.max_nodes_segment` と `search.max_nodes_bridge` で異なります。キンクは `search.kink_max_nodes`（デフォルト 3）の線形ノードを使用します。再帰深度は `search.max_depth`（デフォルト 10）で制限されます。
- オプティマイザー: GSM は pysisyphus `GrowingString` + `StringOptimizer` を使用し、単一構造精密化は常に LBFGS を使用します。
- `--align` での最終マージ規則: `--ref-pdb` が提供された場合、最初の参照 PDB がすべてのペアに使用されます。
- コンソールには状態シーケンス（例: `R --> TS1 --> IM1 --> ... --> P`）とエネルギーダイアグラム構築に使用されるラベル/エネルギーが出力されます。
- 実行後に `summary.md` が生成され、主要成果物と `key_*` 直下ショートカットを一覧できます。
- `summary.log` は一部のフィールド欠落時でも生成可能です。内部で次のキーに既定値を補います:
  `root_out_dir`, `path_module_dir`, `pipeline_mode`, `segments`, `energy_diagrams`。

---

## 関連項目

- [path-opt](path_opt.md) -- シングルパス MEP 最適化（再帰的精密化なし）
- [opt](opt.md) -- 単一構造の構造最適化
- [all](all.md) -- 内部で path-search を呼び出すエンドツーエンドワークフロー
- [trj2fig](trj2fig.md) -- MEP 軌跡からエネルギープロファイルをプロット
