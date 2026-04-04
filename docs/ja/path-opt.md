# `path-opt`

## 概要

> **要約:** ML/MM 計算機で GSM または DMF を使い、**正確に 2 つ**の酵素構造間の MEP を探索します。経路軌跡を書き出し、最高エネルギーイメージ（HEI）を TS 候補としてエクスポートします。

`mlmm path-opt` は、ML/MM 計算機による PySisyphus `GrowingString` を使用して 2 つの酵素状態間の最小エネルギー経路を最適化します。ML/MM 計算機はリンク原子なしで完全な酵素複合体を保持します。ML 領域は `--model-pdb` で定義され、Amber トポロジーは `--parm` から取得され、両端点は全系座標を含む PDB として提供されます。

**2 つ以上**の構造から開始し、反応領域のみを自動精密化するワークフローには、[path-search](path-search.md) を使用してください。

## 最小例

```bash
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --out-dir ./result_path_opt
```

## 出力の見方

- `result_path_opt/final_geometries_trj.xyz`
- `result_path_opt/hei.xyz`
- `result_path_opt/hei.pdb`（PDB 変換が有効な場合）

## よくある例

1. ストリング成長前に両端点を事前最適化する。

```bash
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --preopt --preopt-max-cycles 20000 --out-dir ./result_path_opt_preopt
```

2. まずは高速に確認するため climb を無効化する。

```bash
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --no-climb --max-nodes 8 --out-dir ./result_path_opt_fast
```

3. 凍結原子を指定し、ダンプを保存する。

```bash
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --freeze-atoms "1,3,5,7" --dump --out-dir ./result_path_opt_dump
```

## 使用法

```bash
mlmm path-opt -i REACTANT.pdb PRODUCT.pdb --parm real.parm7 --model-pdb model.pdb \
 -q CHARGE [-m MULT] [--mep-mode gsm|dmf] [--fix-ends/--no-fix-ends] [options]
```

### 例

```bash
# ミニマル呼び出し
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0

# 凍結原子、ノード数増加、YAML 多層設定
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 \
 --freeze-atoms "1,3,5,7" --max-nodes 10 --max-cycles 200 --dump --out-dir ./result_path_opt/

# 端点の事前最適化付き
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0 \
 --preopt --preopt-max-cycles 20000
```

## ワークフロー
1. **端点の読み込み** -- 両方の PDB 構造を読み込み、CLI またはデフォルトから電荷/スピンを解決します。`--parm`、`--model-pdb`、電荷/スピンで ML/MM 計算機を構築します。
2. **事前アライメント** -- 最初の構造以降のすべての端点が最初の構造に Kabsch アライメントされます。`freeze_atoms` が定義されている場合、それらの原子のみが RMSD フィットに参加し、結果の変換がすべての原子に適用されます。
3. **任意の事前最適化** -- `--preopt` の場合、各端点はアライメントとストリング成長の前に LBFGS（同じ ML/MM 計算機を使用）で事前最適化されます。LBFGS サイクル数は `--preopt-max-cycles`（デフォルト: 10000）で制御されます。
4. **経路最適化** -- `--mep-mode gsm` は PySisyphus `GrowingString`（端点込み `(max_nodes + 2)` イメージ）を使用し、`--mep-mode dmf` は Direct Max Flux を使用します。
5. **クライミングイメージ（GSM のみ）** -- `--climb` の場合、ストリングが完全に成長した後にクライミングイメージ精密化が適用され、最高エネルギーイメージ（HEI）が報告されます。
6. **出力** -- 最終経路軌跡と HEI が XYZ および PDB ファイルとして書き出されます。入力が PDB の場合に PDB 変換が実行されます。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH PATH` | 反応物と生成物の PDB 構造。 | 必須 |
| `--parm PATH` | 完全 REAL 系の Amber prmtop。 | 必須 |
| `--model-pdb PATH` | ML 領域を定義する PDB（原子 ID）。`--detect-layer` または `--model-indices` 利用時は省略可。 | _None_ |
| `--model-indices TEXT` | ML 領域のカンマ区切り原子インデックス（範囲指定可、例: `1-5`）。`--model-pdb` 省略時に使用。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` を 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer / --no-detect-layer` | 入力 PDB の B-factor（B=0/10/20）から ML/MM 層を自動検出。無効時は `--model-pdb` または `--model-indices` の指定が必要。 | `True` |
| `-q, --charge INT` | ML 領域の総電荷。 | _None_（`-l` 未指定時は必須） |
| `-l, --ligand-charge TEXT` | 残基ごとの電荷マッピング（例: `SAM:1,PHN:-1`）。`-q` 省略時に合計電荷を導出。PDB 入力または `--ref-pdb` が必要。 | _None_ |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `--mep-mode [gsm\|dmf]` | MEP バックエンド。 | `gsm` |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切り原子インデックス（内部で 0 始まりに変換; YAML `geom.freeze_atoms` とマージ）。 | _None_ |
| `--hess-cutoff FLOAT` | ML 領域からの距離カットオフ (Å)。この範囲内の MM 原子をヘシアン計算に含めます。可動 MM 原子に適用。 | _None_ |
| `--movable-cutoff FLOAT` | ML 領域からの距離カットオフ (Å)。この範囲外の MM 原子を凍結します。`--movable-cutoff` 指定時は `--detect-layer` が無効化されます。 | _None_ |
| `--fix-ends/--no-fix-ends` | GSM 成長中に端点構造を固定（`gs.fix_first/fix_last`）。`defaults.py` の `gs.fix_first` / `gs.fix_last` は共に `True` であるため、実効デフォルトは `True`。 | `False`（CLI）; `True`（GS_KW 実効） |
| `--max-nodes INT` | 内部ストリングノード数（総イメージ = `max_nodes + 2`）。 | `20` |
| `--max-cycles INT` | オプティマイザーマクロ反復上限（成長 + 精密化）。`opt.stop_in_when_full` も設定。 | `300` |
| `--climb/--no-climb` | ストリング完全成長後のクライミングイメージ精密化を有効化。 | `True` |
| `--preopt/--no-preopt` | アライメント/ストリング成長前に各端点を LBFGS で事前最適化。 | `False` |
| `--preopt-max-cycles INT` | 端点事前最適化サイクルの上限。 | `10000` |
| `--thresh TEXT` | 収束プリセット上書き（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`、`never`）。 | _None_（実効: `gau_loose`） |
| `--dump/--no-dump` | `out_dir` 内にオプティマイザー軌跡とリスタートをダンプ。 | `False` |
| `-o, --out-dir TEXT` | 出力ディレクトリ。 | `./result_path_opt/` |
| `--config FILE` | 明示 CLI 指定より前に適用されるベース YAML。 | _None_ |
| `--show-config/--no-show-config` | 解決済み設定（YAML レイヤ情報を含む）を表示して実行継続。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画表示のみを行う。`--help-advanced` に表示。 | `False` |
| `-b, --backend CHOICE` | ML 領域の MLIP バックエンド: `uma`（デフォルト）、`orb`、`mace`、`aimnet2`。 | `uma` |
| `--embedcharge/--no-embedcharge` | xTB 点電荷埋め込み補正の有効化。MM 環境から ML 領域への静電的影響を考慮。 | `False` |
| `--embedcharge-cutoff FLOAT` | xTB 埋め込み用 MM 原子のカットオフ半径（Å）。 | `12.0` |
| `--cmap/--no-cmap` | model parm7 に CMAP（骨格クロスマップ二面角補正）を含めるかどうか。デフォルト: 無効（Gaussian ONIOM と同一）。 | `--no-cmap` |
| `--convert-files/--no-convert-files` | PDB テンプレート利用可能時の XYZ/TRJ から PDB コンパニオン生成の切り替え。 | `True` |

## 出力

```
out_dir/ (デフォルト:./result_path_opt/)
├─ final_geometries_trj.xyz # コメント行にイメージごとのエネルギーを含む XYZ 軌跡
├─ final_geometries.pdb #_trj.xyz と同じだが参照 PDB 順序にマップ
├─ hei.xyz # 最高エネルギーイメージ（XYZ、常に書き出し）
├─ hei.pdb # PDB 形式の HEI（参照 PDB が利用可能な場合）
├─ align_refine/ # 外部アライメント/精密化の成果物
├─ preopt/ # 端点事前最適化出力（--preopt 時）
└─ <optimizer dumps> # --dump または opt.dump_restart > 0 の場合
```

## YAML 設定

マージ順は **defaults < config < 明示指定 CLI < override** です。

### セクション `geom`
- `coord_type`: 座標タイプ（デカルト vs dlc 内部座標）。
- `freeze_atoms`: CLI `--freeze-atoms` とマージされる 1 始まり凍結原子。

### セクション `calc` / `mlmm`
- ML/MM 計算機の設定: `charge`、`spin`、`backend`、`embedcharge`、MLIP モデル設定、`device`、近傍半径、ヘシアンオプション等。

### セクション `gs`
- Growing String 制御: `max_nodes`、`perp_thresh`、再パラメータ化間隔、`max_micro_cycles`、DLC リセット、クライミングトグル/閾値。

### セクション `opt`
- StringOptimizer 設定: `stop_in_when_full`、`scale_step`、`max_cycles`、ダンプフラグ、`reparam_thresh`、`coord_diff_thresh`、`out_dir`、`print_every`。

## 終了コード

| コード | 意味 |
| --- | --- |
| `0` | 成功 |
| `3` | 最適化失敗 |
| `4` | 最終軌跡書き出しエラー |
| `5` | HEI ダンプエラー |
| `130` | キーボード割り込み |
| `1` | 未処理例外 |

---

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [path-search](path-search.md) -- 自動精密化付き再帰的 MEP 探索（2 つ以上の構造用）
- [opt](opt.md) -- 単一構造の構造最適化
- [all](all.md) -- 一気通貫ワークフロー（デフォルトで path-opt を使用、`--refine-path` で再帰的探索）
- [YAML リファレンス](yaml-reference.md) -- `gs`、`opt` の完全な設定オプション
