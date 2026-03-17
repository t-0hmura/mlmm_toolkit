# `scan2d`

## 概要

> **要約:** 調和拘束と ML/MM 緩和による 2 距離（d1, d2）グリッドスキャンを実行します。`-s/--scan-lists` で YAML/JSON スペックファイル（推奨）またはインライン Python リテラルを使用します。

`mlmm scan2d` は `--max-step-size` を使用して 2 つの結合距離の線形グリッドを構築し、適切な拘束を適用して各グリッド点を緩和し、バイアスなしの ML/MM エネルギーを可視化用に記録します。スキャンはまず d1 を反復し d1 拘束のみで構造を緩和し、次に各 d1 値について d2 を反復し両方の拘束を適用します。

各グリッド点のエネルギーはバイアスなしで再評価され、PES グリッドとコンタープロットが作成されます。出力にはグリッド点ごとの XYZ スナップショット、PES をまとめた `surface.csv`、2D コンターマップ（`scan2d_map.png`）、底面投影付き 3D ランドスケープ（`scan2d_landscape.html`）が含まれます。

## 最小例
```bash
mlmm scan2d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan2d.yaml --print-parsed -o ./result_scan2d/
```

## 出力の見方
- `result_scan2d/surface.csv`
- `result_scan2d/grid/point_i000_j000.xyz`
- `result_scan2d/scan2d_map.png` と `result_scan2d/scan2d_landscape.html`

## よくある例
1. YAML spec の解釈結果を先に確認する。
2. インライン `-s` リテラルを使う。
3. `--dump` を有効にして d1 ごとの内側軌跡を保存する。

> **注記:** `-s/--scan-lists` の解釈結果を確認したい場合は `--print-parsed` を追加してください。

## 使用法

```bash
mlmm scan2d -i INPUT.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q CHARGE [-m SPIN] \
 [-s scan2d.yaml | -s "[(I1,J1,LOW1,HIGH1),(I2,J2,LOW2,HIGH2)]"] \
 [--one-based|--zero-based] [--max-step-size FLOAT] [--bias-k FLOAT] \
 [--freeze-atoms "1,3,5"] [--relax-max-cycles INT] [--thresh PRESET] \
 [--dump/--no-dump] [--out-dir DIR] \
 [--preopt/--no-preopt] [--baseline {min|first}] [--zmin FLOAT] [--zmax FLOAT]
```

### 例

```bash
# 推奨: YAML/JSON spec
cat > scan2d.yaml << 'YAML'
one_based: true
pairs:
 - [12, 45, 1.30, 3.10]
 - [10, 55, 1.20, 3.20]
YAML
mlmm scan2d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan2d.yaml --print-parsed

# 代替: インライン Python リテラル
mlmm scan2d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s "[(12,45,1.30,3.10),(10,55,1.20,3.20)]"

# TRJ ダンプ付き LBFGS スキャン、コンタープロットの固定カラースケール
mlmm scan2d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s "[(12,45,1.30,3.10),(10,55,1.20,3.20)]" \
 --max-step-size 0.20 --dump -o ./result_scan2d/ --preopt --baseline min \
 --zmin 0.0 --zmax 40.0
```

## YAML/JSON スペックフォーマット（推奨）

`-s/--scan-lists` は YAML/JSON ファイルを自動検出します。ファイルパスを渡すとスペックモードになります:

```yaml
one_based: true # 任意; デフォルトは CLI の --one-based/--zero-based
pairs:
 - [12, 45, 1.30, 3.10]
 - [10, 55, 1.20, 3.20]
```

- `pairs` は必須で、正確に 2 つの四つ組を含む必要があります。
- 各四つ組は `(i, j, low_A, high_A)` です。
- インデックスは整数または PDB セレクター（インラインリテラルと同じ）が使用可能です。

## インラインリテラルフォーマット

`-s/--scan-lists` がファイルパスでない値を受け取ると、**単一の Python リテラル**文字列として評価されます。シェルクォートに注意してください。

### 基本構造

リテラルは正確に **2 つ**の四つ組 `(atom1, atom2, low_A, high_A)` の Python リストです:

```
-s '[(atom1, atom2, low_A, high_A), (atom3, atom4, low_A, high_A)]'
```

- シェルが括弧やスペースを解釈しないよう、リテラル全体を**シングルクォート**で囲んでください。
- 各四つ組は 1 つのスキャン軸を定義します: `atom1`--`atom2` 間の距離を `low_A` から `high_A` までスキャンします。
- `scan` と異なり、**1 つのリテラル**のみ受け付けます（マルチステージ非対応）。

### 原子の指定

原子は**整数インデックス**または **PDB セレクター文字列**で指定できます:

| 方法 | 例 | 備考 |
| --- | --- | --- |
| 整数インデックス | `(1, 5, 1.30, 3.10)` | デフォルトは 1 始まり（`--one-based`） |
| PDB セレクター | `("TYR,285,CA", "MMT,309,C10", 1.30, 3.10)` | 残基名、残基番号、原子名 |

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
-s '[("TYR,285,CA","MMT,309,C10",1.30,3.10),("TYR,285,CB","MMT,309,C11",1.20,3.20)]'

# 正しい: 整数インデックスは内側のクォート不要
-s '[(1, 5, 1.30, 3.10), (2, 8, 1.20, 3.20)]'

# 非推奨: 外側をダブルクォートにすると内側のクォートをエスケープする必要あり
-s "[(\"TYR,285,CA\",\"MMT,309,C10\",1.30,3.10),...]"
```

## ワークフロー
1. **入力と事前最適化** -- 酵素 PDB を読み込み、電荷/スピンを解決し、ML/MM 計算機（MLIP バックエンド + hessian_ff）を構築し、`--preopt` の場合は任意でバイアスなし事前最適化を実行。`-b/--backend` で ML バックエンドを選択（デフォルト: `uma`）、`--embedcharge` で xTB 点電荷埋め込み補正を有効化可能。
2. **グリッド構築** -- `-s/--scan-lists`（YAML/JSON スペックファイルまたはインラインリテラル）からターゲットを 2 つの四つ組に解析し、インデックスを正規化（デフォルト 1 始まりまたは `"TYR,285,CA"` のような PDB 原子セレクター）。`ceil(|high - low| / h) + 1` 点の線形グリッドを構築（`h = --max-step-size`）。
3. **外側ループ（d1）** -- 各 d1 値について、**d1 拘束のみ**で系を緩和。
4. **内側ループ（d2）** -- 現在の d1 での各 d2 値について、最も近い収束済み構造から開始し**両方の拘束**で緩和。
5. **エネルギー評価** -- 各 (i, j) ペアで ML/MM エネルギーをバイアスなしで評価し `surface.csv` に記録。
6. **可視化** -- `scan2d_map.png`（2D コンター）と `scan2d_landscape.html`（3D サーフェス）を書き出し。`--zmin/--zmax` でカラースケールをクランプ。ベースライン: `--baseline min` は最小エネルギーをゼロに; `--baseline first` は (i=0, j=0) グリッド点をゼロに。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力酵素複合体 PDB（必須）。 | 必須 |
| `--parm PATH` | 酵素の Amber parm7 トポロジー（必須）。 | 必須 |
| `--model-pdb PATH` | ML 領域を定義する PDB。`--detect-layer` 有効時はオプション。 | _None_ |
| `--model-indices TEXT` | ML 領域のカンマ区切り原子インデックス（範囲指定可）。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` を 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer / --no-detect-layer` | 入力 PDB の B 因子から ML/MM レイヤーを検出。 | `True` |
| `-q, --charge INT` | ML 領域の総電荷。 | _None_（`-l` 未指定時は必須） |
| `-l, --ligand-charge TEXT` | 残基ごとの電荷マッピング（例: `GPP:-3,SAM:1`）。`-q` 省略時に合計電荷を導出。 | _None_ |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切りインデックス。 | _None_ |
| `--hess-cutoff FLOAT` | ML 領域からの距離カットオフ (Å) — ヘシアン計算に含める MM 原子を指定。`--detect-layer` と併用可能。 | _None_ |
| `--movable-cutoff FLOAT` | ML 領域からの可動 MM 原子の距離カットオフ (Å)。指定すると `--detect-layer` が無効化されます。 | _None_ |
| `-s, --scan-lists TEXT` | スキャンターゲット: YAML/JSON スペックファイルパス（自動検出、`pairs` に 2 四つ組）またはインライン Python リテラル `"[(i1,j1,low1,high1),(i2,j2,low2,high2)]"`。インデックスは整数または PDB 原子セレクター。 | 必須 |
| `--one-based / --zero-based` | `-s/--scan-lists` の `(i,j)` インデックスを 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--print-parsed/--no-print-parsed` | `-s/--scan-lists` 解釈後のペア情報を表示。 | `False` |
| `--max-step-size FLOAT` | ステップごとの最大距離増分 (Å)。グリッド密度を決定。 | `0.20` |
| `--bias-k FLOAT` | 調和ウェル強度 k (eV/Å²)。 | `300.0` |
| `--relax-max-cycles INT` | バイアス緩和ごとの最大 LBFGS サイクル。 | `10000` |
| `--dump/--no-dump` | d1 スライスごとの内側 d2 スキャン TRJ を書き出し。 | `False` |
| `-o, --out-dir TEXT` | 基本出力ディレクトリ。 | `./result_scan2d/` |
| `--thresh TEXT` | 収束プリセット（`gau_loose\|gau\|gau_tight\|gau_vtight\|baker\|never`）。 | `baker` |
| `--config FILE` | ベース YAML 設定ファイル（最初に適用）。 | _None_ |
| `--ref-pdb FILE` | `--input` が XYZ の場合の参照 PDB トポロジー。 | _None_ |
| `--preopt/--no-preopt` | スキャン前にバイアスなし事前最適化を実行。 | `False` |
| `--baseline {min,first}` | 相対エネルギーの基準（kcal/mol）。 | `min` |
| `--zmin FLOAT` | コンターカラースケールの下限（kcal/mol）。 | 自動スケール |
| `--zmax FLOAT` | コンターカラースケールの上限（kcal/mol）。 | 自動スケール |
| `-b, --backend CHOICE` | ML 領域の MLIP バックエンド: `uma`（デフォルト）、`orb`、`mace`、`aimnet2`。 | _None_（内部で `uma` を適用） |
| `--embedcharge/--no-embedcharge` | xTB 点電荷埋め込み補正の有効化。MM 環境から ML 領域への静電的影響を考慮。 | `False` |
| `--embedcharge-cutoff FLOAT` | xTB 埋め込み用 MM 原子のカットオフ半径（Å）。 | `12.0` |
| `--convert-files/--no-convert-files` | PDB テンプレート利用可能時の XYZ/TRJ から PDB コンパニオン生成の切り替え。 | `True` |

## 出力

```
out_dir/ (デフォルト:./result_scan2d/)
├── surface.csv # PES グリッド: i, j, d1_A, d2_A, energy_hartree, energy_kcal, bias_converged
├── scan2d_map.png # 2D コンターマップ
├── scan2d_landscape.html # 3D サーフェス可視化（Plotly）
├── grid/
│ ├── point_i###_j###.xyz # 各 (i, j) ペアの緩和ジオメトリ
│ ├── point_i###_j###.pdb # PDB コンパニオン（入力が PDB の場合）
│ ├── preopt_i###_j###.xyz # 事前最適化構造（--preopt 時）
│ └── inner_path_d1_###_trj.xyz # d1 スライスごとの内側 d2 軌跡（--dump 時）
└── (stdout) # 進捗とエネルギーサマリー
```

## YAML 設定

```yaml
geom:
 coord_type: cart
 freeze_atoms: []
calc:
 charge: 0
 spin: 1
mlmm:
 real_parm7: real.parm7
 model_pdb: ml_region.pdb
opt:
 thresh: baker
 max_cycles: 10000
 dump: false
 out_dir: ./result_scan2d/
lbfgs:
 max_step: 0.3
 out_dir: ./result_scan2d/
bias:
 k: 300.0
```

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細なトラブルシューティングガイド

- [scan](scan.md) -- 1D 結合距離駆動スキャン
- [scan3d](scan3d.md) -- 3D 距離グリッドスキャン
- [opt](opt.md) -- 単一構造の構造最適化
- [all](all.md) -- end-to-endワークフロー
