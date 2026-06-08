# `scan3d`

調和拘束と ML/MM 緩和による 3 距離（d1, d2, d3）グリッドスキャンを実行します。`mlmm scan3d` は d1、d2、d3 のネストループを実行し、ML/MM 計算機（`mlmm.backends.mlmm_calc.mlmm`）を使用して適切な拘束で各点を緩和します。ML 領域は `--model-pdb` から、Amber パラメータは `--parm` から読み取られ、MLIP バックエンドは `-b/--backend` で選択（デフォルト: `uma`）、オプティマイザーは PySisyphus LBFGS です。`-s/--scan-lists` で YAML/JSON スペックファイル（推奨）またはインライン Python リテラルを使用します。

## 使いどころ

- 3 つの結合距離を跨ぐ 3D PES のマッピング。
- `--csv` で事前計算した surface を読み込めば、スキャンを再実行せずに再描画のみ可能。

## 実行例

```bash
# 最小: YAML spec から 3D スキャンを実行
mlmm scan3d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan3d.yaml -o ./result_scan3d/
```
（`--print-parsed` を追加すると、解釈されたスキャンスペックを検証し、GPU 計算を実行せずに終了します。）

```bash
# 推奨: YAML/JSON spec
cat > scan3d.yaml << 'YAML'
one_based: true
pairs:
 - [12, 45, 1.30, 3.10]
 - [10, 55, 1.20, 3.20]
 - [15, 60, 1.10, 3.00]
YAML
mlmm scan3d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan3d.yaml --print-parsed
```

```bash
# インライン Python リテラル、事前最適化・--dump・カスタム出力ディレクトリ付き
mlmm scan3d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s "[(12,45,1.30,3.10),(10,55,1.20,3.20),(15,60,1.10,3.00)]" \
 --max-step-size 0.20 --dump -o ./result_scan3d/ \
 --preopt --baseline min
```

## 入力

コマンド形式:

```bash
mlmm scan3d -i INPUT.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q CHARGE [-m MULT] \
 [--csv precomputed_surface.csv] \
 [-s scan3d.yaml | -s "[(I1,J1,LOW1,HIGH1),(I2,J2,LOW2,HIGH2),(I3,J3,LOW3,HIGH3)]"] \
 [--one-based|--zero-based] [--max-step-size FLOAT] [--bias-k FLOAT] \
 [--freeze-atoms "1,3,5"] [--relax-max-cycles INT] [--thresh PRESET] \
 [--dump/--no-dump] [--out-dir DIR] \
 [--preopt/--no-preopt] [--baseline {min|first}] [--zmin FLOAT] [--zmax FLOAT]
```

| 入力 | 必須 | 備考 |
| --- | --- | --- |
| `-i, --input` | はい（`--csv` 指定時を除く） | 完全酵素 PDB（リンク原子なし）。 |
| `--parm` | はい（`--csv` 指定時を除く） | 完全酵素の Amber parm7 トポロジー。 |
| `--model-pdb` | 推奨 | ML 領域を定義する PDB（または `--model-indices` を使用）。 |
| `-s, --scan-lists` | はい | YAML/JSON スペックファイルパス（`pairs` に 3 四つ組）またはインライン Python リテラル（3 つの四つ組 `(i,j,low,high)`）。 |
| `-q, --charge` | はい（`-l` または `--csv` 指定時を除く） | ML 領域の総電荷。 |
| `--csv` | 任意 | 事前計算済み `surface.csv` を読み込み、スキャンを実行せずにプロットを生成。 |

### 入力構文 — YAML/JSON スペックフォーマット（推奨）

`-s/--scan-lists` は YAML/JSON ファイルを自動検出します。ファイルパスを渡すとスペックモードになります:

```yaml
one_based: true # 任意; デフォルトは CLI の --one-based/--zero-based
pairs:
 - [12, 45, 1.30, 3.10]
 - [10, 55, 1.20, 3.20]
 - [15, 60, 1.10, 3.00]
```

- `pairs` は必須で、正確に 3 つの四つ組を含む必要があります。
- 各四つ組は `(i, j, low_A, high_A)` です。
- インデックスは整数または PDB セレクター（`--scan-lists` と同じ）が使用可能です。

### 入力構文 — インラインリテラルフォーマット

`-s/--scan-lists` がファイルパスでない値を受け取ると、**単一の Python リテラル**文字列として評価されます。シェルクォートに注意してください。

リテラルは正確に **3 つ**の四つ組 `(atom1, atom2, low_A, high_A)` の Python リストです:

```
-s '[(atom1, atom2, low_A, high_A), (atom3, atom4, low_A, high_A), (atom5, atom6, low_A, high_A)]'
```

- シェルが括弧やスペースを解釈しないよう、リテラル全体を**シングルクォート**で囲んでください。
- 各四つ組は 1 つのスキャン軸を定義します: `atom1`--`atom2` 間の距離を `low_A` から `high_A` までスキャンします。
- `scan` と異なり、**1 つのリテラル**のみ受け付けます（マルチステージ非対応）。

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

クォート規則:

```bash
# 正しい: リスト全体をシングルクォート、内側のセレクター文字列をダブルクォート
-s '[("TYR,285,CA","MMT,309,C10",1.30,3.10),("TYR,285,CB","MMT,309,C11",1.20,3.20),("TYR,285,CG","MMT,309,C12",1.10,3.00)]'

# 正しい: 整数インデックスは内側のクォート不要
-s '[(1, 5, 1.30, 3.10), (2, 8, 1.20, 3.20), (3, 12, 1.10, 3.00)]'

# 非推奨: 外側をダブルクォートにすると内側のクォートをエスケープする必要あり
-s "[(\"TYR,285,CA\",\"MMT,309,C10\",1.30,3.10),...]"
```

## 処理の流れ
1. `geom_loader` で構造を読み込み、CLI から電荷/スピンを解決し、`--preopt` の場合は任意でバイアスなし事前最適化を実行。
2. `-s/--scan-lists`（YAML/JSON スペックファイルまたはインラインリテラル）からターゲットを解析して 3 つの四つ組にします（デフォルト 1 始まりインデックス、`--zero-based` 指定時は 0 始まり）。PDB 入力の場合、各原子エントリは整数インデックスまたは `"TYR,285,CA"` のようなセレクター文字列が使用可能。区切り文字はスペース、カンマ、スラッシュ、バッククォート、バックスラッシュ。
3. 外側ループ `d1[i]`: d1 拘束のみで緩和。d1 値が最も近い以前のスキャン済みジオメトリから開始。
4. 中間ループ `d2[j]`: d1 と d2 の拘束で緩和。最も近い (d1, d2) ジオメトリから開始。
5. 内側ループ `d3[k]`: 3 つの拘束すべてで緩和。バイアスなしエネルギーを測定（評価時にバイアス除去）し、拘束ジオメトリと収束フラグを書き出し。
6. スキャン完了後、`surface.csv` を組み立て、kcal/mol ベースラインシフト（`--baseline {min|first}`）を適用し、3D RBF 補間アイソサーフェスプロット（`scan3d_density.html`）を生成（`--zmin/--zmax` を尊重）。

## 出力

```text
out_dir/ (デフォルト:./result_scan3d/)
 surface.csv # グリッドメタデータ（d1, d2, d3, energy, convergence）
 scan3d_density.html # 3D エネルギーアイソサーフェス可視化
 grid/point_i###_j###_k###.xyz # 各グリッド点の緩和ジオメトリ
 grid/point_i###_j###_k###.pdb # PDB コンパニオン（B 因子: ML=0, Movable-MM=10, Frozen=20）
 grid/inner_path_d1_###_d2_###_trj.xyz # --dump が True の場合のみ
```

ファイル名タグ `i###_j###_k###` は Å の 100 分の 1 の整数（d1×100, d2×100, d3×100）であり、ステップ番号ではありません。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 完全酵素 PDB（リンク原子なし）。 | `--csv` 指定時を除き必須 |
| `--parm PATH` | 完全酵素の Amber parm7 トポロジー。 | `--csv` 指定時を除き必須 |
| `--model-pdb PATH` | ML 領域を定義する PDB。 | _None_ |
| `--model-indices TEXT` | 明示的な ML 領域原子インデックス（`--model-pdb` の代替）。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` のインデックス規約。 | `True`（1 始まり） |
| `--detect-layer / --no-detect-layer` | B 因子から ML/MM レイヤーを自動検出。 | `True` |
| `-q, --charge INT` | ML 領域の総電荷。 | _None_（`-l` または `--csv` 未指定時は必須） |
| `-l, --ligand-charge TEXT` | 残基ごとの電荷マッピング（例: `GPP:-3,SAM:1`）。`-q` 省略時に合計電荷を導出。 | _None_ |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `--freeze-atoms TEXT` | 1 始まりカンマ区切りの凍結原子インデックス。 | _None_ |
| `--hess-cutoff FLOAT` | ML 領域からの距離カットオフ (Å) — ヘシアン計算に含める MM 原子を指定。`--detect-layer` と併用可能。 | _None_ |
| `--movable-cutoff FLOAT` | ML 領域からの可動 MM 原子の距離カットオフ (Å)。指定すると `--detect-layer` が無効化されます。 | _None_ |
| `-s, --scan-lists TEXT` | スキャンターゲット: YAML/JSON スペックファイルパス（自動検出、`pairs` に 3 四つ組）またはインライン Python リテラル。`i`/`j` は整数インデックスまたは PDB 原子セレクター。 | 必須 |
| `--csv FILE` | 事前計算済み `surface.csv` を読み込みスキャンなしでプロット生成。 | _None_ |
| `--one-based / --zero-based` | `(i, j)` インデックスを 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--print-parsed/--no-print-parsed` | `-s/--scan-lists` 解釈後のペア情報を表示。 | `False` |
| `--max-step-size FLOAT` | ステップごとの最大距離増分 (Å)。グリッド密度を制御。 | `0.20` |
| `--bias-k FLOAT` | 調和ウェル強度 k (eV/Å²)。 | `300.0` |
| `--relax-max-cycles INT` | バイアス緩和ごとの最大オプティマイザーサイクル。 | `10000` |
| `--dump/--no-dump` | (d1, d2) スライスごとの内側 d3 スキャン TRJ を書き出し。 | `False` |
| `-o, --out-dir TEXT` | グリッドとプロットの出力ディレクトリルート。 | `./result_scan3d/` |
| `--thresh TEXT` | 収束プリセット上書き（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`、`never`）。 | `baker` |
| `--config FILE` | ベース YAML 設定ファイル（最初に適用）。 | _None_ |
| `--ref-pdb FILE` | 非 PDB 入力用の参照 PDB トポロジー。 | _None_ |
| `--preopt/--no-preopt` | スキャン前にバイアスなし最適化を実行。 | `False` |
| `--baseline {min,first}` | kcal/mol エネルギーをグローバル最小値または `(i,j,k)=(0,0,0)` がゼロになるようシフト。 | `min` |
| `--zmin FLOAT` | アイソサーフェスカラーバンドの手動下限（kcal/mol）。 | 自動スケール |
| `--zmax FLOAT` | アイソサーフェスカラーバンドの手動上限（kcal/mol）。 | 自動スケール |
| `-b, --backend CHOICE` | ML 領域の MLIP バックエンド: `uma`、`orb`、`mace`、`aimnet2`。 | `uma` |
| `--embedcharge/--no-embedcharge` | xTB 点電荷埋め込み補正の有効化。MM 環境から ML 領域への静電的影響を考慮。 | `False` |
| `--embedcharge-cutoff FLOAT` | xTB 埋め込み用 MM 原子のカットオフ半径（Å）。 | `12.0` |
| `--cmap/--no-cmap` | model parm7 に CMAP（骨格クロスマップ二面角補正）を含めるかどうか。デフォルト: 無効（Gaussian ONIOM と同一）。 | `--no-cmap` |
| `--convert-files/--no-convert-files` | PDB テンプレート利用可能時の XYZ/TRJ から PDB コンパニオン生成の切り替え。 | `True` |

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
 out_dir: ./result_scan3d/
lbfgs:
 max_step: 0.3
 out_dir: ./result_scan3d/
bias:
 k: 300.0
```

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [scan](scan.md) -- 1D 結合距離駆動スキャン
- [scan2d](scan2d.md) -- 2D 距離グリッドスキャン
- [opt](opt.md) -- 構造最適化（スキャン前に実行する場合が多い）
- [all](all.md) -- 一気通貫ワークフロー
- [YAML リファレンス](yaml-reference.md) -- スキャンの完全な設定オプション
