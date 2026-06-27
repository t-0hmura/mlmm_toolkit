# `path-opt`

`mlmm path-opt` は、PySisyphus `GrowingString`（デフォルト）または DMF（`--mep-mode dmf`）を用いて、**正確に 2 つ**の層付き酵素構造間の最小エネルギー経路（MEP）を最適化します。最適化には ML/MM 計算機を使い、リンク原子なしで完全な酵素複合体を保持します。ML 領域は `--model-pdb` で定義し、Amber トポロジーは `--parm` から取得し、両端点は全系座標を含む PDB として与えます。経路軌跡を書き出し、最高エネルギーイメージ（HEI）を TS 候補としてエクスポートします。2 つの層付き端点が明確で中間体が無いと予想されるときに使います。再帰分割も結合変化駆動の分解も行わない、`path-search` のシンプル版です。**2 つ以上**の構造から開始し、反応領域のみを自動精密化するワークフローには、代わりに [path-search](path-search.md) を使用してください。

## 実行例

```bash
# ミニマル呼び出し
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --out-dir ./result_path_opt
```

```bash
# ストリング成長前に両端点を事前最適化する
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --preopt --preopt-max-cycles 20000 --out-dir ./result_path_opt_preopt
```

```bash
# まずは高速に確認するため climb を無効化する
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --no-climb --max-nodes 8 --out-dir ./result_path_opt_fast
# 凍結原子の指定とダンプの保存: --freeze-atoms "1,3,5,7" --dump
```

コマンド形式:

```bash
mlmm path-opt -i REACTANT.pdb PRODUCT.pdb --parm real.parm7 --model-pdb model.pdb \
 -q CHARGE [-m MULT] [--mep-mode gsm|dmf] [--fix-ends/--no-fix-ends] [options]
```

`mlmm path-opt --help` は主要オプションを、`mlmm path-opt --help-advanced` は全オプション一覧を表示します。

## 処理の流れ
1. **端点の読み込み** -- 両方の PDB 構造を読み込み、CLI またはデフォルトから電荷/スピンを解決します。`--parm`、`--model-pdb`、電荷/スピンで ML/MM 計算機を構築します。
2. **事前アライメント** -- 最初の構造以降のすべての端点が最初の構造に Kabsch アライメントされます。`freeze_atoms` が定義されている場合、それらの原子のみが RMSD フィットに参加し、結果の変換がすべての原子に適用されます。
3. **任意の事前最適化** -- `--preopt` の場合、各端点はアライメントとストリング成長の前に LBFGS（同じ ML/MM 計算機を使用）で事前最適化されます。LBFGS サイクル数は `--preopt-max-cycles`（デフォルト: 10000）で制御されます。
4. **経路最適化** -- `--mep-mode gsm` は PySisyphus `GrowingString`（端点込み `(max_nodes + 2)` イメージ）を使用し、`--mep-mode dmf` は Direct Max Flux を使用します。
5. **クライミングイメージ（GSM のみ）** -- `--climb` の場合、ストリングが完全に成長した後にクライミングイメージ精密化が適用され、最高エネルギーイメージ（HEI）が報告されます。
6. **出力** -- 最終経路軌跡と HEI が XYZ および PDB ファイルとして書き出されます。入力が PDB の場合に PDB 変換が実行されます。

## 出力

結果は通常、以下のファイルを開いて確認します。

- `result_path_opt/final_geometries_trj.xyz`
- `result_path_opt/hei.xyz`
- `result_path_opt/hei.pdb`（PDB 変換が有効な場合）

```text
out_dir/ (デフォルト:./result_path_opt/)
├─ final_geometries_trj.xyz # コメント行にイメージごとのエネルギーを含む XYZ 軌跡
├─ final_geometries.pdb #_trj.xyz と同じだが参照 PDB 順序にマップ
├─ hei.xyz # 最高エネルギーイメージ（XYZ、常に書き出し）
├─ hei.pdb # PDB 形式の HEI（参照 PDB が利用可能な場合）
├─ align_refine/ # 外部アライメント/精密化の成果物
├─ preopt/ # 端点事前最適化出力（--preopt 時）
└─ <optimizer dumps> # --dump または opt.dump_restart > 0 の場合
```

## CLI オプション

全フラグ一覧は生成された[コマンドリファレンス](../reference/commands/index.md)にあります。以下の表は説明を要するオプションを扱います。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH PATH` | 反応物と生成物の PDB 構造（全系座標）。 | 必須 |
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
| `--fix-ends/--no-fix-ends` | 経路成長中に端点構造を固定（`gs.fix_first/fix_last`）。 | `True` |
| `--max-nodes INT` | 内部ストリングノード数（総イメージ = `max_nodes + 2`）。 | `20` |
| `--max-cycles INT` | オプティマイザーマクロ反復上限（成長 + 精密化）。`opt.stop_in_when_full` も設定。 | `300` |
| `--climb/--no-climb` | ストリング完全成長後のクライミングイメージ精密化を有効化。 | `True` |
| `--preopt/--no-preopt` | アライメント/ストリング成長前に各端点を LBFGS で事前最適化。 | `True` |
| `--preopt-max-cycles INT` | 端点事前最適化サイクルの上限。 | `10000` |
| `--thresh TEXT` | 収束プリセット上書き（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`、`never`）。 | _None_（実効: `gau_loose`） |
| `--mm-backend [hessian_ff\|openmm]` | MM バックエンド（解析的ヘシアン vs OpenMM 有限差分）。 | `hessian_ff` |
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

## YAML 設定

マージ順は **defaults < config < 明示指定 CLI < override** です。関連セクションは `geom`（`coord_type`、`freeze_atoms`）、`calc` / `mlmm`（ML/MM 計算機の設定）、`gs`（Growing String 制御）、`opt`（StringOptimizer 設定）です。

完全なスキーマ（全キーとデフォルト）: [YAML リファレンス](yaml-reference.md)。

## 終了コード

| コード | 意味 |
| --- | --- |
| `0` | 成功 |
| `3` | 最適化失敗 |
| `4` | 最終軌跡書き出しエラー |
| `5` | HEI ダンプエラー |
| `130` | キーボード割り込み |
| `1` | 未処理例外 |

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [path-search](path-search.md) -- 自動精密化付き再帰的 MEP 探索（2 つ以上の構造用）
- [opt](opt.md) -- 単一構造の構造最適化
- [all](all.md) -- 一気通貫ワークフロー（既定で単一パス path-opt、`--refine-path` で再帰 path-search）
- [YAML リファレンス](yaml-reference.md) -- `gs`、`opt` の完全な設定オプション
