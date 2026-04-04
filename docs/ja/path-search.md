# `path-search`

## 概要

> **要約:** 再帰的 GSM セグメンテーションにより 2 つ以上の構造から連続した MEP を構築します。結合変化のある領域のみを自動的に精密化し、最高エネルギーイメージ（HEI）を TS 候補としてエクスポートします。

`mlmm path-search` は GSM を使用して 2 つ以上の構造にわたる連続した最小エネルギー経路（MEP）を構築します。共有結合変化が検出された領域のみを選択的に精密化し、解決されたサブパスを 1 つの軌跡に統合します。

再帰的分解により多段階反応を自動検出し、各素反応ステップの詳細な MEP を構築します。ただし、複雑な多段階反応の検出は困難な場合があり、入力中間体やスキャン仕様、収束閾値の調整など手動での試行錯誤が必要になることがあります。

端点が **2 つ**だけで再帰的精密化が不要な場合は、[path-opt](path-opt.md) がより簡単な選択です。

## 最小例

```bash
mlmm path-search -i reactant.pdb product.pdb --parm real.parm7 \
 --model-pdb ml_region.pdb -q 0 --out-dir ./result_path_search
```

## 出力の見方

- `result_path_search/mep_trj.xyz`
- `result_path_search/summary.json`
- `result_path_search/summary.log`
- `result_path_search/mep_plot.png`（プロット生成時）

## よくある例

1. 中間体を含む多段経路を構築する。

```bash
mlmm path-search -i R.pdb IM1.pdb IM2.pdb P.pdb --parm real.parm7 \
 --model-pdb ml_region.pdb -q -1 --out-dir ./result_path_search_multi
```

2. 事前最適化とアライメントを無効にして軽く試す。

```bash
mlmm path-search -i reactant.pdb product.pdb --parm real.parm7 \
 --model-pdb ml_region.pdb -q 0 --no-preopt --no-align --max-nodes 8 \
 --out-dir ./result_path_search_fast
```

## 使用法

```bash
mlmm path-search -i R.pdb IM1.pdb P.pdb \
 --parm real.parm7 --model-pdb ml_region.pdb -q CHARGE [-m MULT]
 [--mep-mode gsm|dmf] [--refine-mode peak|minima]
 [--freeze-atoms "1,3,5"] [--max-nodes N] [--max-cycles N] [--climb/--no-climb]
 [--opt-mode grad|hess]
 [--thresh PRESET] [--dump/--no-dump] [--out-dir DIR]
 [--show-config/--no-show-config] [--dry-run/--no-dry-run]
```

### 例

```bash
# ミニマルなポケットのみの 2 状態間 MEP
mlmm path-search -i reactant.pdb product.pdb --parm real.parm7 \
 --model-pdb ml_region.pdb -q 0

# YAML 上書き、凍結原子付きマルチステップ経路
mlmm path-search -i R.pdb IM1.pdb P.pdb --parm real.parm7 \
 --model-pdb ml_region.pdb -q -1 --freeze-atoms "1,3,5" \
 --ref-pdb holo_template.pdb --out-dir ./run_ps
```

## ワークフロー

1. **初期セグメント（隣接ペア A->B ごと; GSM/DMF）** -- 選択した MEP エンジン（`--mep-mode`）を実行して粗い MEP を取得し、最高エネルギーイメージ（HEI）を特定。
2. **HEI 周辺の局所緩和** -- `--refine-mode`（`peak`: HEI+/-1、`minima`: 最近傍局所極小）で種点を選び、単一構造オプティマイザー（`opt-mode`）で精密化して近傍の極小（`End1`、`End2`）を回復。
3. **ねじれ vs 精密化の判定**:
 - `End1` と `End2` の間に共有結合変化が検出されない場合、その領域を*ねじれ*として扱い、`search.kink_max_nodes` 個の線形ノードを挿入して各ノードを個別に最適化。
 - それ以外の場合、`End1` と `End2` の間で**精密化セグメント（GSM）**を起動して障壁を鮮明化。
4. **選択的再帰** -- `(A->End1)` と `(End2->B)` の結合変化を `bond` 閾値で比較。共有結合の更新を含むサブ区間のみに再帰。再帰深度は `search.max_depth` で制限。
5. **統合とブリッジ** -- 解決されたサブパスを連結し、RMSD <= `search.stitch_rmsd_thresh` の重複端点を削除。2 つの統合部分の間の RMSD ギャップが `search.bridge_rmsd_thresh` を超える場合、選択中の `--mep-mode` でブリッジ MEP セグメントを挿入。インターフェース自体に結合変化がある場合、ブリッジの代わりに新たな再帰セグメントを生成。
6. **任意のアライメント** -- 事前最適化後、`--align` で入力を剛体アラインし凍結を精密化。セグメントをプロット/分析用にアノテーション。

結合変化検出は `bond` YAML セクションの閾値を使用する `bond_changes.compare_structures` に依存します。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順の 2 つ以上の完全酵素 PDB。`-i` を繰り返すか、1 つのフラグの後に複数パスを渡す。 | 必須 |
| `--parm PATH` | 完全酵素複合体の Amber parm7 トポロジー。 | 必須 |
| `--model-pdb PATH` | ML/MM の ML（高レベル）領域原子を定義する PDB。`--detect-layer` または `--model-indices` 利用時は省略可。 | _None_ |
| `--model-indices TEXT` | ML 領域のカンマ区切り原子インデックス（範囲指定可、例: `1-5`）。`--model-pdb` 省略時に使用。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` を 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer / --no-detect-layer` | 入力 PDB の B 因子（B=0/10/20）から ML/MM レイヤーを検出。無効時は `--model-pdb` または `--model-indices` が必要。 | `True` |
| `-q, --charge INT` | ML 領域の電荷（整数）。`-l` 未指定時は必須。 | _None_ |
| `-l, --ligand-charge TEXT` | 残基ごとの電荷マッピング（例: `SAM:1,PHN:-1`）。`-q` 省略時に合計電荷を導出。PDB 入力または `--ref-pdb` が必要。 | _None_ |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `--mep-mode [gsm\|dmf]` | セグメント/ブリッジ探索に使う MEP バックエンド。 | `gsm` |
| `--refine-mode [peak\|minima]` | HEI 精密化の種点ルール。 | `gsm` は `peak`、`dmf` は `minima` |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切りインデックス（YAML `geom.freeze_atoms` とマージ）。 | _None_ |
| `--hess-cutoff FLOAT` | ML 領域からの Hessian-MM 原子の距離カットオフ (Å)。可動 MM 原子に適用。 | _None_ |
| `--movable-cutoff FLOAT` | ML 領域からの可動 MM 原子の距離カットオフ (Å)。これを超える MM 原子は凍結。指定時は `--detect-layer` が無効化。 | _None_ |
| `--max-nodes INT` | セグメント GSM の内部ノード数。 | `10` |
| `--max-cycles INT` | GSM マクロサイクルの最大数。 | `300` |
| `--climb/--no-climb` | セグメント GSM の TS 精密化を有効化。 | `True` |
| `--opt-mode [grad\|hess]` | 単一構造オプティマイザープリセット（`grad` = LBFGS、`hess` = RFO）。 | `grad` |
| `--preopt/--no-preopt` | セグメンテーション前に端点を LBFGS で事前最適化。 | `False` |
| `--align / --no-align` | 事前最適化後に入力を剛体アライメント。 | 有効 |
| `--thresh TEXT` | 収束プリセット（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`、`never`）。 | _None_（実質: `gau_loose`） |
| `--dump/--no-dump` | オプティマイザーダンプを保存。 | `False` |
| `-o, --out-dir PATH` | 出力ディレクトリ。 | `./result_path_search/` |
| `--ref-pdb PATH...` | XYZ→PDB 変換・トポロジー参照用の完全テンプレート PDB。 | _None_ |
| `--config FILE` | 明示 CLI 指定より前に適用されるベース YAML。 | _None_ |
| `--show-config/--no-show-config` | 解決済み設定（YAML レイヤ情報を含む）を表示して実行継続。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画表示のみを行う。`--help-advanced` に表示。 | `False` |
| `-b, --backend CHOICE` | ML 領域の MLIP バックエンド: `uma`（デフォルト）、`orb`、`mace`、`aimnet2`。 | `uma` |
| `--embedcharge/--no-embedcharge` | xTB 点電荷埋め込み補正の有効化。MM 環境から ML 領域への静電的影響を考慮。 | `False` |
| `--embedcharge-cutoff FLOAT` | xTB 埋め込み用 MM 原子のカットオフ半径（Å）。 | `12.0` |
| `--cmap/--no-cmap` | model parm7 に CMAP（骨格クロスマップ二面角補正）を含めるかどうか。デフォルト: 無効（Gaussian ONIOM と同一）。 | `--no-cmap` |
| `--convert-files/--no-convert-files` | PDB テンプレート利用可能時の XYZ/TRJ から PDB コンパニオン生成の切り替え。 | `True` |

## 出力

```text
out_dir/ (デフォルト:./result_path_search/)
 summary.json # MEP レベルの実行サマリー（完全設定ダンプなし）
 summary.log # 人間が読めるサマリー
 mep_trj.xyz # 最終 MEP（常に書き出し）
 mep.pdb # 最終 MEP（参照テンプレート利用可能時は PDB）
 mep_seg_XX_trj.xyz / mep_seg_XX.pdb # セグメント別経路
 hei_seg_XX.xyz / hei_seg_XX.pdb # 結合変化セグメントごとの HEI
 mep_plot.png # イメージインデックスに対する Delta-E プロファイル（trj2fig より）
 energy_diagram_MEP.png # 反応物基準の状態レベルエネルギーダイアグラム（kcal/mol）
 seg_000_*/ # セグメントレベルの GSM と精密化成果物
```

## YAML 設定

マージ順は **defaults < config < 明示指定 CLI < override** です。

YAML ルートはマッピングでなければなりません。受け付けるセクション:

- **`geom`** -- `coord_type`（デフォルト `"cart"`）、`freeze_atoms`（1 始まりインデックス）。
- **`calc` / `mlmm`** -- ML/MM 計算機設定: `input_pdb`、`real_parm7`、`model_pdb`、`model_charge`、`model_mult`、バックエンド選択（`backend`、`embedcharge`）、UMA 制御（`uma_model`、`uma_task_name`、`hessian_calc_mode`）、デバイス選択、凍結原子。
- **`gs`** -- Growing String 設定: `max_nodes`、`climb`、`climb_rms`、`climb_fixed`、`reparam_every_full`、`reparam_check`。
- **`opt`** -- StringOptimizer 制御: `max_cycles`、`print_every`、`dump`、`dump_restart`、`out_dir`。
- **`lbfgs`** -- HEI+/-1 精密化用の単一構造オプティマイザー制御: `keep_last`、`beta`、`gamma_mult`、`max_step`、`control_step`、`double_damp`、`mu_reg`、`max_mu_reg_adaptions`。
- **`bond`** -- 結合変化検出: `bond_factor`、`margin_fraction`、`delta_fraction`。
- **`search`** -- 再帰ロジック: `max_depth`、`stitch_rmsd_thresh`、`bridge_rmsd_thresh`、`max_nodes_segment`、`max_nodes_bridge`、`kink_max_nodes`、`max_seq_kink`、`refine_mode`。

---

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細なトラブルシューティングガイド

- [path-opt](path-opt.md) -- MEP 最適化（再帰的精密化なし）
- [opt](opt.md) -- 単一構造の構造最適化
- [all](all.md) -- all が内部で path-opt（デフォルト）または path-search（`--refine-path` 指定時）を呼び出す一気通貫ワークフロー
- [trj2fig](trj2fig.md) -- MEP 軌跡からエネルギープロファイルをプロット
