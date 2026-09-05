# `scan`

`mlmm scan` は、レイヤー分けした単一の酵素構造から、調和拘束によって1つ以上の原子間距離を目標値へ駆動し、各ステップで L-BFGS により構造を緩和します。この ML/MM スキャンで粗い反応軌跡と、下流の MEP 精密化用の中間体・生成物候補を生成します。入力には PDB/mmCIF、または `--ref-pdb` を伴う XYZ を使用できます。`-s/--scan-lists` で、YAML/JSON スペックファイル（推奨）またはインライン Python リテラルとして目標距離を定義します。

## スキャン座標のステージ構成

3 要素タプルの入力では、1 リテラルまたは YAML の 1 つの `stages` 要素が 1 ステージを定義します。
同一ステージ内の複数距離tupleは協奏的に駆動し、複数のリテラル/要素は
多段階scanとして前ステージの端点から順次実行されます。`scan2d` / `scan3d`は独立な
距離軸を用いて energy landscape を探索し、PESを描画します。

## 実行例

以下の例では、`pocket.pdb` が `real.parm7` に対応する全系構造で、`ml_region.pdb` が ML 領域（リンク水素なし）です。

コマンド形式:

```bash
mlmm scan -i INPUT.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q CHARGE [-m MULT] \
 (-s scan.yaml | -s "[(I,J,TARGET_ANG)]") [options]
```

スペックファイルによるスキャン（`--print-parsed` を追加すると、解釈したスキャンスペックを検証し GPU 計算を実行せずに終了します）:

```bash
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan.yaml -o ./result_scan
```

インライン Python リテラル:

```bash
# インライン Python リテラル
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s "[(12,45,2.20)]"
```

ステージごとの軌跡を保存して確認する:

```bash
# ステージごとの軌跡を保存して確認する
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan.yaml --dump -o ./result_scan_dump
```

## 処理の流れ

1. `geom_loader` で構造を読み込み、CLI またはデフォルトから電荷/スピンを解決します。ML/MM calculatorに `--parm`、`--model-pdb`、`-q/--charge`、任意で `-m/--multiplicity` を提供します。
2. 任意でバイアスなし事前最適化（`--preopt`）を実行し、開始点を緩和します。
3. `-s/--scan-lists`（YAML/JSON スペックファイルまたはインラインリテラル）からステージターゲットを解析し、`(i, j)` インデックスを正規化します（デフォルトは 1 始まり）。PDB メタデータを利用できる場合、各エントリは整数インデックスまたは `'TYR,285,CA'` のような原子セレクター文字列のいずれかで指定可能です。セレクターフィールドはスペース、カンマ、スラッシュ、バッククォート、バックスラッシュで区切ることができ、順序は任意です。
4. 結合ごとの変位を計算してステップに分割します:
 - スキャンタプル `[(i, j, target_A)]` に対し、`delta = target - current_distance_A` を計算。
 - `--max-step-size = h` の場合、ステージは `N = ceil(max(|delta|) / h)` 回のバイアス付き緩和を実行。
 - 各ペアの増分変化は `step_k = delta / N` (Å)。ステップ `s` での一時ターゲットは `r_k(s) = r_k(0) + s * step_k`。
5. すべてのステップを進み、調和拘束ポテンシャル `E_bias = sum 1/2 * k * (|r_i - r_j| - target_k)^2` を適用して L-BFGS で極小化。`k` は `--bias-k`（eV/Å²）から取得され、Hartree/Bohr^2 に一度だけ変換されます。座標は PySisyphus 用に Bohr で保存され、レポート時に内部変換されます。
6. 各ステージの最後のステップ後、任意でバイアスなし緩和（`--endopt`）を実行してから共有結合変化を報告し `result.*` ファイルを書き出します。
7. すべてのステージで繰り返します。

## 出力

各ステージは最終ジオメトリとバイアスステップ軌跡を `stage_XX/` 配下に書き出し、ルートに連結軌跡を生成します。最初に確認するファイルはステージ別の `result.xyz` と、常に生成される `scan_trj.xyz` です。PDB companion は `--convert-files` が有効で PDB テンプレートを利用できる場合に生成されます。

```
out_dir/ (デフォルト:./result_scan/)
├─ scan_trj.xyz              # 全ステージ連結軌跡（常に書き出し）
├─ scan.pdb                  # --convert-files と PDB テンプレートがある場合
├─ preopt/                   # --preopt が True の場合
│  ├─ result.xyz
│  └─ result.pdb             # --convert-files と PDB テンプレートがある場合
└─ stage_XX/                 # ステージごとに 1 フォルダ（k = 01..K）
   ├─ result.xyz             # 最終（endopt 済みの可能性あり）ジオメトリ
   ├─ result.pdb             # --convert-files と PDB テンプレートがある場合
   ├─ scan_trj.xyz           # ステージ別バイアスステップフレーム（常に書き出し）
   └─ scan.pdb               # --convert-files と PDB テンプレートがある場合
```

## CLI オプション

完全なフラグ一覧は生成された [コマンドリファレンス](../reference/commands/index.md) にあります。下表は説明が必要なオプションを扱います。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力 PDB/mmCIF、またはトポロジー用に `--ref-pdb` を伴う XYZ。 | 必須 |
| `--parm PATH` | 完全 REAL 系の Amber prmtop。 | 必須 |
| `--model-pdb PATH` | ML 領域を定義する PDB（原子 ID）。`--detect-layer` 有効時または `--model-indices` 指定時は省略可能。 | _None_ |
| `--model-indices TEXT` | ML 領域原子インデックス（カンマ区切り、範囲指定可）。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` を 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer / --no-detect-layer` | 入力 PDB の B 因子から ML/MM レイヤーを自動検出。 | 有効 |
| `-q, --charge INT` | ML 領域の総電荷。 | _None_（`-l` 未指定時は必須） |
| `-l, --ligand-charge TEXT` | 残基ごとの電荷マッピング（例: `GPP:-3,SAM:1`）。`-q` 省略時に合計電荷を導出。 | _None_ |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切り原子インデックス（YAML `geom.freeze_atoms` とマージ）。 | _None_ |
| `--movable-cutoff FLOAT` | 可動 MM 距離カットオフ (Å)。指定すると `--detect-layer` を無効化。 | _None_ |
| `-s, --scan-lists TEXT` | スキャンターゲット: YAML/JSON スペックファイルパス（自動検出）または `(i, j, target_A)` 3 要素タプルもしくは `(i, j, start, end)` 4 要素タプル（双方向スキャン）を含むインライン Python リテラル。単一フラグの後に複数リテラルを供給可能。`i`/`j` は整数インデックスまたは `"TYR,285,CA"` のような PDB 原子セレクターが使用可能。 | 必須 |
| `--one-based/--zero-based` | 原子インデックスを 1 始まり（デフォルト）または 0 始まりとして解釈。 | `True`（1 始まり） |
| `--print-parsed/--no-print-parsed` | 解釈したスキャン対象を表示し、計算せず終了。 | `False` |
| `--max-step-size FLOAT` | ステップごとのスキャン結合の最大変化量 (Å)。積分ステップ数を制御。 | `0.20` |
| `--bias-k FLOAT` | 調和バイアス強度 `k`（eV/Å²）。 | `300` |
| `--max-cycles INT` | 各バイアスステップおよび pre/end 最適化ステージの L-BFGS サイクル上限。 | `100000` |
| `--relax-max-cycles INT` | `--max-cycles` の互換エイリアス（指定時は上書き）。 | `--max-cycles`を継承 |
| `--preopt/--no-preopt` | スキャン前にバイアスなし最適化を実行。 | `False` |
| `--endopt/--no-endopt` | 各ステージ後にバイアスなし最適化を実行。 | `False` |
| `--dump/--no-dump` | ステップごとのオプティマイザ軌跡ファイルをダンプ。`scan_trj.xyz` は常に書き出され、PDB/CIF companion には `--convert-files` と参照トポロジーが必要です。 | `False` |
| `-o, --out-dir TEXT` | 出力ディレクトリルート。 | `./result_scan/` |
| `--thresh TEXT` | 収束プリセット（`gau_loose\|gau\|gau_tight\|gau_vtight\|baker\|never`）。 | _None_（`gau` を継承） |
| `--config FILE` | ベース YAML 設定ファイル（最初に適用）。 | _None_ |
| `--ref-pdb FILE` | `--input` が XYZ の場合の参照 PDB トポロジー。 | _None_ |
| `-b, --backend CHOICE` | ML 領域の MLIP バックエンド: `uma`、`orb`、`mace`、`aimnet2`。 | `uma` |
| `--cmap/--no-cmap` | REAL と MODEL の両 MM 層で CMAP を保持します。 | `--cmap` |
| `--mm-backend [hessian_ff\|openmm]` | MM バックエンド。Hessian 構築法は `calc.mm_fd` が別に制御します（デフォルト `true`: 有限差分）。 | `hessian_ff` |
| `--link-atom-method [scaled\|fixed]` | リンク原子の配置: scaled（$g$ 因子）または固定 1.09/1.01 Å。 | `scaled` |
| `--out-json/--no-out-json` | `result.json` を `out_dir` に書き出す。 | `False` |
| `--dry-run/--no-dry-run` | オプションの検証と実行計画の表示のみ行い、スキャンは実行しない。`--help-advanced` に表示。 | `False` |
| `--convert-files/--no-convert-files` | PDB テンプレートが利用可能な場合に、XYZ/TRJ から対応する PDB を生成するかどうかを切り替え。 | `True` |

## スキャン対象の構文

**YAML/JSON スペックフォーマット（推奨）**

`-s/--scan-lists` は YAML/JSON ファイルを自動検出します。ファイルパスを渡すとスペックモードになります:

```yaml
one_based: true # 任意; デフォルトは CLI の --one-based/--zero-based
stages:
 - [[12, 45, 2.20]]
 - [[10, 55, 1.35], [23, 34, 1.80]]
```

- `stages` は必須です。
- 各ステージは `(i, j, target_A)` の 3 要素タプルのリストです。
- インデックスは整数、または PDB メタデータを利用できる場合は PDB セレクターが使用可能で、インラインリテラルと同じです。

**インラインリテラルフォーマット**

`-s/--scan-lists` がファイルパスでない値を受け取ると、**Python リテラル**文字列として評価されます。シェルクォートに注意してください。

各リテラルは 3 要素タプル `(atom1, atom2, target_A)` の Python リストです:

```
-s '[(atom1, atom2, target_A),...]'
```

- シェルが括弧やスペースを解釈しないよう、リテラル全体を**シングルクォート**で囲んでください。
- 各 3 要素タプルは `atom1`--`atom2` 間の距離を `target_A` に向けて駆動します。
- 1 つのリテラル = 1 つの**ステージ**です。複数ステージの場合、**単一の** `-s/--scan-lists` フラグの後に複数リテラルを渡します（フラグを繰り返さないでください）。

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

クォート規則:

```bash
# 正しい: リスト全体をシングルクォート、内側のセレクター文字列をダブルクォート
-s '[("TYR,285,CA","MMT,309,C10",1.35)]'

# 正しい: 整数インデックスは内側のクォート不要
-s '[(1, 5, 2.0)]'

# 非推奨: 外側をダブルクォートにすると内側のクォートをエスケープする必要あり
-s "[(\"TYR,285,CA\",\"MMT,309,C10\",1.35)]"
```

単一の `-s/--scan-lists` フラグの後に複数リテラルを渡します。各リテラルが 1 ステージになります:

```bash
# ステージ 1: 1 つの結合を 1.35 Å に駆動
# ステージ 2: 2 つの結合を同時に駆動
-s \
 '[("TYR,285,CA","MMT,309,C10",1.35)]' \
 '[("TYR,285,CA","MMT,309,C10",2.20),("TYR,285,CB","MMT,309,C11",1.80)]'
```

ステージは順次実行され、各ステージは前のステージの緩和結果から開始します。

**同期ステージと順次ステージの例**

```bash
# 協奏的: 1 ステージ、2 つの距離を同時に駆動
mlmm scan -i r.pdb --parm enzyme.parm7 -l 'LIG:Q' \
    -s '[(1,5,1.40),(7,9,1.60)]' -o result_concerted

# 段階的: 2 つの順次ステージ
mlmm scan -i r.pdb --parm enzyme.parm7 -l 'LIG:Q' \
    -s '[(1,5,1.40)]' \
       '[(7,9,0.95)]' -o result_staged
```

協奏scanの後に[`path-search`](path-search.md)を用いると、多段階候補の
segmentを構築できます。4-tupleは別の双方向構文であり、2ステージへ
展開されます。

**双方向スキャン（4-tuple）**

3-tuple `(i, j, target)` の代わりに **4-tuple** `(i, j, start, end)` を指定すると、現在の構造から両方向にスキャンします。CLI は各 4-tuple を自動的に 2 ステージに展開します:

1. **パス 1:** `i`--`j` の距離を現在の値から `start` に向けて駆動。
2. **パス 2:** 初期構造を復元し、`i`--`j` の距離を `end` に向けて駆動。

連結軌跡は `start → 初期構造 → end` の順に並び、出発構造を通る連続的な経路が得られます。

```bash
# 双方向スキャン: 結合 12--45 を現在の構造から
# 1.35 Å（パス 1）と 2.50 Å（パス 2）に向けて駆動
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s '[(12, 45, 1.35, 2.50)]'
```

これは 2 つの手動ステージの間にジオメトリリセットを行うのと同等ですが、スクリプトを書く必要がありません。同じリテラル内で 3-tuple と 4-tuple を混在させることもできます。

## バリアの方向を読む

`scan` は各点のエネルギーを保存します。反応障壁は、TS・IRC で検証した構造のエネルギーから求めます。

| 量 | 式 |
| --- | --- |
| 順方向バリア | `E(TS) − E(reactant)` |
| 逆方向の反応障壁 | `E(TS) − E(product)` |

反応物・生成物の割り当ては、IRC 後に最適化した端点構造で確認してください。

## YAML 設定

スキャンは共有の `geom`（`coord_type`、`freeze_atoms`）、`calc` / `mlmm`（ML/MM calculator設定）、`opt` / `lbfgs`（オプティマイザ）の各セクションに加え、`bias`（`k`、調和強度（eV/Å²））と MLIP ベースの結合変化検出用 `bond` セクションを読み込みます。

- `coord_type`: 共有キーですが、拘束付き scan はマージ後に `cart` へ正規化するため DLC は有効になりません。
- `freeze_atoms`: CLI `--freeze-atoms` とマージされる 1 始まり凍結原子。

### セクション `calc` / `mlmm`
- ML/MM calculatorの設定: `model_charge`、`model_mult`、`backend`、MLIP モデル設定、`device`、近傍半径、Hessian オプション等。

### セクション `opt` / `lbfgs`
- オプティマイザ設定: `thresh`、`max_cycles`、`print_every`、ステップ制御、ラインサーチ、ダンプフラグ。

### セクション `bias`
- `k`（`300`）: 調和強度（eV/Å²）。

### セクション `bond`
- MLIP ベースの結合変化検出:
 - `device`（`"auto"`）: グラフ分析用の MLIP デバイス。
 - `bond_factor`（`1.20`）: カットオフ用の共有結合半径スケーリング。
 - `margin_fraction`（`0.05`）: 比較用の許容割合。
 - `delta_fraction`（`0.05`）: 結合形成/切断と判定する最小相対変化。

全スキーマ（すべてのキーとデフォルト）: [YAML リファレンス](yaml-reference.md)。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [scan2d](scan2d.md) -- 2D 距離グリッドスキャン
- [scan3d](scan3d.md) -- 3D 距離グリッドスキャン
- [opt](opt.md) -- 単一構造の構造最適化
- [all](all.md) -- 単一構造入力の `--scan-lists` 付き一気通貫ワークフロー
- [path-search](path-search.md) -- スキャン端点を中間体として使用する MEP 探索
