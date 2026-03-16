# `scan`

## 概要

> **要約:** ML/MM 計算機を使用して、調和拘束による結合距離スキャンで反応座標を駆動します。`-s/--scan-lists` で YAML/JSON スペックファイル（推奨）またはインライン Python リテラルとしてターゲット距離を定義します。

`mlmm scan` は ML/MM 計算機（`mlmm.mlmm_calc.mlmm`）を使用して調和拘束による段階的な結合距離駆動スキャンを実行します。各ステップで一時的なターゲットが更新され、拘束ウェルが適用され、LBFGS で構造が緩和されます。ML/MM 計算機は MLIP バックエンド（デフォルト: UMA、`-b/--backend` で選択）と hessian_ff を結合します。`--embedcharge` で xTB 点電荷埋め込み補正を有効化できます。

## 最小例

```bash
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan.yaml --print-parsed -o ./result_scan
```

## 出力の見方

- `result_scan/stage_01/result.pdb`（または `result.xyz`）
- `result_scan/stage_02/result.pdb`（または `result.xyz`）
- `result_scan/stage_*/scan_trj.xyz` と `scan.pdb`（常に生成）

## よくある例

1. YAML の解釈結果を表示して入力を確認する。

```bash
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan.yaml --print-parsed
```

2. インラインリテラル入力を使う。

```bash
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s "[(12,45,2.20)]"
```

3. ステージごとの軌跡を保存して確認する。

```bash
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan.yaml --dump -o ./result_scan_dump
```

> **注記:** `-s/--scan-lists` の解釈結果を確認したい場合は `--print-parsed` を追加してください。

## 使用法
```bash
mlmm scan -i INPUT.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q CHARGE [-m MULT] \
 [-s scan.yaml | -s "[(I,J,TARGET_ANG)]"] [options]
```

### 例
```bash
# 推奨: YAML/JSON spec
cat > scan.yaml << 'YAML'
one_based: true
stages:
 - [[12, 45, 2.20]]
 - [[10, 55, 1.35], [23, 34, 1.80]]
YAML
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan.yaml --print-parsed

# 代替: インライン Python リテラル
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s "[(12,45,2.20)]"

# ダンプ付き 2 ステージ、凍結原子、YAML 上書き
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q -1 -m 1 --freeze-atoms "1,3,5" -s "[(12,45,2.20)]" \
 "[(10,55,1.35),(23,34,1.80)]" --max-step-size 0.20 --dump
```

## YAML/JSON スペックフォーマット（推奨）

`-s/--scan-lists` は YAML/JSON ファイルを自動検出します。ファイルパスを渡すとスペックモードになります:

```yaml
one_based: true # 任意; デフォルトは CLI の --one-based/--zero-based
stages:
 - [[12, 45, 2.20]]
 - [[10, 55, 1.35], [23, 34, 1.80]]
```

- `stages` は必須です。
- 各ステージは `(i, j, target_A)` の三つ組のリストです。
- インデックスは整数または PDB セレクター（PDB 入力時）が使用可能で、インラインリテラルと同じです。

## インラインリテラルフォーマット

`-s/--scan-lists` がファイルパスでない値を受け取ると、**Python リテラル**文字列として評価されます。シェルクォートに注意してください。

### 基本構造

各リテラルは三つ組 `(atom1, atom2, target_A)` の Python リストです:

```
-s '[(atom1, atom2, target_A),...]'
```

- シェルが括弧やスペースを解釈しないよう、リテラル全体を**シングルクォート**で囲んでください。
- 各三つ組は `atom1`--`atom2` 間の距離を `target_A` に向けて駆動します。
- 1 つのリテラル = 1 つの**ステージ**です。複数ステージの場合、**単一の** `-s/--scan-lists` フラグの後に複数リテラルを渡します（フラグを繰り返さないでください）。

### 原子の指定

原子は**整数インデックス**または **PDB セレクター文字列**で指定できます:

| 方法 | 例 | 備考 |
| --- | --- | --- |
| 整数インデックス | `(1, 5, 2.0)` | デフォルトは 1 始まり（`--one-based`） |
| PDB セレクター | `("TYR,285,CA", "MMT,309,C10", 2.0)` | 残基名、残基番号、原子名 |

PDB セレクターのトークンは、カンマ `,`、スペース、スラッシュ `/`、バッククォート `` ` ``、バックスラッシュ `\` のいずれかで区切れます。トークンの順序は自由です。

```bash
# 以下はすべて同じ原子を指定:
"TYR,285,CA"
"TYR 285 CA"
"TYR/285/CA"
"285,TYR,CA" # 順序は自由
```

### クォート規則

```bash
# 正しい: リスト全体をシングルクォート、内側のセレクター文字列をダブルクォート
-s '[("TYR,285,CA","MMT,309,C10",1.35)]'

# 正しい: 整数インデックスは内側のクォート不要
-s '[(1, 5, 2.0)]'

# 非推奨: 外側をダブルクォートにすると内側のクォートをエスケープする必要あり
-s "[(\"TYR,285,CA\",\"MMT,309,C10\",1.35)]"
```

### 複数ステージ

単一の `-s/--scan-lists` フラグの後に複数リテラルを渡します。各リテラルが 1 ステージになります:

```bash
# ステージ 1: 1 つの結合を 1.35 Å に駆動
# ステージ 2: 2 つの結合を同時に駆動
-s \
 '[("TYR,285,CA","MMT,309,C10",1.35)]' \
 '[("TYR,285,CA","MMT,309,C10",2.20),("TYR,285,CB","MMT,309,C11",1.80)]'
```

ステージは順次実行され、各ステージは前のステージの緩和結果から開始します。**`-s/--scan-lists` フラグを繰り返さないでください** -- 単一のフラグの後にすべてのステージリテラルを供給してください。

## ワークフロー
1. `geom_loader` で構造を読み込み、CLI またはデフォルトから電荷/スピンを解決します。ML/MM 計算機に `--parm`、`--model-pdb`、`-q/--charge`、任意で `-m/--multiplicity` を提供します。
2. 任意でバイアスなし事前最適化（`--preopt`）を実行し、開始点を緩和します。
3. `-s/--scan-lists`（YAML/JSON スペックファイルまたはインラインリテラル）からステージターゲットを解析し、`(i, j)` インデックスを正規化します（デフォルトは 1 始まり）。入力が PDB の場合、各エントリは整数インデックスまたは `'TYR,285,CA'` のような原子セレクター文字列のいずれかで指定可能です。セレクターフィールドはスペース、カンマ、スラッシュ、バッククォート、バックスラッシュで区切ることができ、順序は任意です。
4. 結合ごとの変位を計算してステップに分割します:
 - スキャンタプル `[(i, j, target_A)]` に対し、`delta = target - current_distance_A` を計算。
 - `--max-step-size = h` の場合、ステージは `N = ceil(max(|delta|) / h)` 回のバイアス付き緩和を実行。
 - 各ペアの増分変化は `delta_k = delta_k / N` (Å)。ステップ `s` での一時ターゲットは `r_k(s) = r_k(0) + s * delta_k`。
5. すべてのステップを進み、調和ウェル `E_bias = sum 1/2 * k * (|r_i - r_j| - target_k)^2` を適用して LBFGS で極小化。`k` は `--bias-k`（eV/Å²）から取得され、Hartree/Bohr^2 に一度変換されます。座標は PySisyphus 用に Bohr で保存され、レポート時に内部変換されます。
6. 各ステージの最後のステップ後、任意でバイアスなし緩和（`--endopt`）を実行してから共有結合変化を報告し `result.*` ファイルを書き出します。
7. すべてのステージで繰り返します。任意の軌跡は `--dump` が `True` の場合のみダンプされます。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力 PDB（またはトポロジー用に `--ref-pdb` 付きの XYZ）。 | 必須 |
| `--parm PATH` | 完全 REAL 系の Amber prmtop。 | 必須 |
| `--model-pdb PATH` | ML 領域を定義する PDB（原子 ID）。`--detect-layer` 有効時または `--model-indices` 指定時は省略可能。 | _None_ |
| `--model-indices TEXT` | ML 領域原子インデックス（カンマ区切り、範囲指定可）。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` を 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer / --no-detect-layer` | 入力 PDB の B 因子から ML/MM レイヤーを検出。 | `True` |
| `-q, --charge INT` | ML 領域の総電荷。 | _None_（`-l` 未指定時は必須） |
| `-l, --ligand-charge TEXT` | 残基ごとの電荷マッピング（例: `GPP:-3,SAM:1`）。`-q` 省略時に合計電荷を導出。 | _None_ |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切り原子インデックス（YAML `geom.freeze_atoms` とマージ）。 | _None_ |
| `--hess-cutoff FLOAT` | ML 原子からの MM Hessian 距離カットオフ (Å)。 | _None_ |
| `--movable-cutoff FLOAT` | 可動 MM 距離カットオフ (Å)。指定すると `--detect-layer` を無効化。 | _None_ |
| `-s, --scan-lists TEXT` | スキャンターゲット: YAML/JSON スペックファイルパス（自動検出）または `(i, j, target_A)` タプルを含むインライン Python リテラル。各リテラルが 1 ステージ。単一フラグの後に複数リテラルを供給可能。`i`/`j` は整数インデックスまたは `"TYR,285,CA"` のような PDB 原子セレクターが使用可能。 | 必須 |
| `--one-based/--zero-based` | 原子インデックスを 1 始まり（既定）または 0 始まりとして解釈。 | `True`（1 始まり） |
| `--print-parsed/--no-print-parsed` | `-s/--scan-lists` 解釈後のステージ情報を表示。 | `False` |
| `--max-step-size FLOAT` | ステップごとのスキャン結合の最大変化量 (Å)。積分ステップ数を制御。 | `0.20` |
| `--bias-k FLOAT` | 調和バイアス強度 `k`（eV/Å²）。 | `300` |
| `--opt-mode {grad,hess,lbfgs,rfo,light,heavy}` | `mlmm all` からの転送互換オプション。現状の `scan` 緩和は mode に関わらず LBFGS を使用。 | _None_ |
| `--max-cycles INT` | 各バイアスステップおよび pre/end 最適化ステージの最大 LBFGS サイクル。 | `10000` |
| `--relax-max-cycles INT` | `--max-cycles` の互換エイリアス（指定時は上書き）。 | _None_ |
| `--preopt/--no-preopt` | スキャン前にバイアスなし最適化を実行。 | `True` |
| `--endopt/--no-endopt` | 各ステージ後にバイアスなし最適化を実行。 | `True` |
| `--dump/--no-dump` | ステップごとのオプティマイザー軌跡ファイルをダンプ。注: `scan_trj.xyz`/`scan.pdb` はこのフラグに関係なく常に書き出されます。 | `False` |
| `-o, --out-dir TEXT` | 出力ディレクトリルート。 | `./result_scan/` |
| `--thresh TEXT` | 収束プリセット（`gau_loose\|gau\|gau_tight\|gau_vtight\|baker\|never`）。 | _None_（`gau` を継承） |
| `--config FILE` | ベース YAML 設定ファイル（最初に適用）。 | _None_ |
| `--ref-pdb FILE` | `--input` が XYZ の場合の参照 PDB トポロジー。 | _None_ |
| `-b, --backend CHOICE` | ML 領域の MLIP バックエンド: `uma`（デフォルト）、`orb`、`mace`、`aimnet2`。 | _None_（内部で `uma` を適用） |
| `--embedcharge/--no-embedcharge` | xTB 点電荷埋め込み補正の有効化。MM 環境から ML 領域への静電的影響を考慮。 | `False` |
| `--embedcharge-cutoff FLOAT` | xTB 埋め込み用 MM 原子のカットオフ半径（Å）。 | `12.0` |
| `--dry-run/--no-dry-run` | オプションの検証と実行計画の表示のみ行い、スキャンは実行しない。`--help-advanced` に表示。 | `False` |
| `--convert-files/--no-convert-files` | PDB テンプレート利用可能時の XYZ/TRJ から PDB コンパニオン生成の切り替え。 | `True` |

## 出力
```
out_dir/ (デフォルト:./result_scan/)
├─ scan_trj.xyz              # 全ステージ結合軌跡（常に書き出し）
├─ scan.pdb                  # 結合 PDB コンパニオン（PDB 入力のみ、常に書き出し）
├─ preopt/                   # --preopt が True の場合
│  ├─ result.xyz
│  └─ result.pdb             # PDB 入力のみ
└─ stage_XX/                 # ステージごとに 1 フォルダ（k = 01..K）
   ├─ result.xyz             # 最終（endopt 済みの可能性あり）ジオメトリ
   ├─ result.pdb             # 入力が PDB の場合
   ├─ scan_trj.xyz           # ステージ別バイアスステップフレーム（常に書き出し）
   └─ scan.pdb               # scan_trj.xyz の PDB 版（PDB 入力のみ、常に書き出し）
```

## YAML 設定

- `coord_type`: 座標タイプ（デカルト vs dlc 内部座標）。
- `freeze_atoms`: CLI `--freeze-atoms` とマージされる 0 始まり凍結原子。

### セクション `calc` / `mlmm`
- ML/MM 計算機の設定: `charge`、`spin`、`backend`、`embedcharge`、MLIP モデル設定、`device`、近傍半径、ヘシアンオプション等。

### セクション `opt` / `lbfgs`
- オプティマイザー設定: `thresh`、`max_cycles`、`print_every`、ステップ制御、ラインサーチ、ダンプフラグ。

### セクション `bias`
- `k`（`300`）: 調和強度（eV/Å²）。

### セクション `bond`
- MLIP ベースの結合変化検出:
 - `device`（`"cuda"`）: グラフ分析用の MLIP デバイス。
 - `bond_factor`（`1.20`）: カットオフ用の共有結合半径スケーリング。
 - `margin_fraction`（`0.05`）: 比較用の許容分数。
 - `delta_fraction`（`0.05`）: 結合形成/切断をフラグする最小相対変化。

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [scan2d](scan2d.md) -- 2D 距離グリッドスキャン
- [scan3d](scan3d.md) -- 3D 距離グリッドスキャン
- [opt](opt.md) -- 単一構造の構造最適化
- [all](all.md) -- 単一構造入力の `--scan-lists` 付きエンドツーエンドワークフロー
- [path-search](path_search.md) -- スキャン端点を中間体として使用する MEP 探索
