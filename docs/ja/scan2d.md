# `scan2d`（2 次元の拘束付きグリッドスキャン）

`scan2d` サブコマンドは、層付き酵素構造の 2 つの座標の格子の各点を調和拘束で保って ML/MM 計算機で緩和し、拘束を外したエネルギーを記録して、反応の 2D エネルギーマップを作ります。各軸は範囲で指定し、距離 `(i,j,low,high)`（Å）、角度 `(i,j,k,low,high)`、二面角 `(i,j,k,l,low,high)`（度）を使えます。

---

## 主な用途

* **TS 領域の見当づけ**: 経路探索や TS 最適化の前に、反応物側と生成物側の谷を結ぶ鞍点がどのあたりにあるかを調べる
* **MEP の前の地形の確認**: 最小エネルギー経路（MEP）を精密化する前に、結合形成とプロトン移動のような 2 つの変化が同時に起きるか、順に起きるかを確かめる
* **結果の判定**: [`surface.csv`](#surfacecsv-の読み方と判定) に各点の収束の可否が入り、`result.json` には `scientific_status` が入る

ML 領域の計算バックエンドにはデフォルトの **UMA**（Meta）のほか、`-b/--backend` オプションで **ORB**、**MACE**、**AIMNet2**、DFT（`dft`）も選択可能です。1 つ以上の座標を動かして 1 本の経路を作るには [`scan`](scan.md) を、3 つの座標の格子には [`scan3d`](scan3d.md) を使います。

---

## 基本的な実行例

例の `pocket.pdb` は `real.parm7` に対応する全系の構造で、`ml_region.pdb` はそのうちの ML 領域（リンク水素なし）を選びます。

原子は同梱の酵素の例（`examples/beza/1.R.pdb`）のもので、この PDB は chain の欄が空です。そのため、原子は chain を省いた 3 項目（残基名・残基番号・原子名）を任意の順序で、カンマか空白で区切って書きます（`"SAM,320,CS1"`）。

### 1. YAML スペックファイルからの実行

2 つの範囲を `pairs:` に書きます。デフォルトの刻み幅 0.2 Å では、このファイルから 9 × 9 の格子ができます。

```yaml
# scan2d.yaml
pairs:
  - ["SAM,320,CS1", "GPP,321,C7", 1.50, 3.00]
  - ["GPP,321,H11", "GLU,186,OE2", 0.90, 2.50]
```

構造・電荷と、このファイルを渡します。`--out-json` を付けると `result.json` も出力します。

```bash
mlmm scan2d -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s scan2d.yaml --out-json -o ./result_scan2d/
```

等高線図は `result_scan2d/scan2d_map.png`、3D 曲面は `scan2d_landscape.html` で確認できます。`result.json` には `scientific_status` と使える点の数が入ります。

### 2. インラインリテラルでの指定

同じ 2 つの範囲を、1 つのリテラルとしてコマンドラインに書けます。

```bash
mlmm scan2d -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s '[("SAM,320,CS1","GPP,321,C7",1.50,3.00),("GPP,321,H11","GLU,186,OE2",0.90,2.50)]'
```

### 3. 事前最適化・軌跡の保存・基準の指定

スキャンの前に入力構造を最適化し、内側ループの軌跡を保存して、相対エネルギーを使える点の最小値から測ります。`--max-step-size 0.20`・`--opt-mode grad`・`--baseline min` はデフォルトで、ここでは明示しています。

```bash
mlmm scan2d -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s '[("SAM,320,CS1","GPP,321,C7",1.50,3.00),("GPP,321,H11","GLU,186,OE2",0.90,2.50)]' \
    --max-step-size 0.20 --dump -o ./result_scan2d/ --opt-mode grad \
    --preopt --baseline min
```

---

## 処理の仕組みと計算仕様

1. **開始構造と格子**:
ML 領域の {ref}`電荷 <ja-charge-specification>` は `-q` または `-l` から決まります。`--preopt` を付けると、まず拘束なしで入力構造を最適化します。収束しなかった場合は入力構造を使います。各軸には両端を含めて ceil(|high − low| / h) + 1 個の等間隔の値ができます。h は距離では `--max-step-size`（Å）、角度と二面角では `--max-angle-step-size` と `--max-dihedral-step-size`（度）です。値は開始構造に近いものから順に計算します。
2. **外側と内側のループ**:
d₁ の各値で、d₁ の拘束だけをかけて構造を緩和します。続く内側ループで、両方の拘束をかけて d₂ を走査します。各点は、すでに収束した点のうち最も近いものから始めます。
3. **各点の緩和**:
調和拘束 E = ½ k (q − q_target)² が各座標 q を目標値に保ち（k は `--restraint-k`）、残りの構造を ML/MM 計算機で L-BFGS（`--opt-mode grad`、デフォルト）または RFO（有理関数最適化、`--opt-mode hess`）により緩和します。凍結 MM 層の原子は動きません。そのあと拘束を外してエネルギーを計算し、構造を `grid/` に書き出します。
4. **表と図**:
最後の点のあと、全点を `surface.csv` にまとめます。使える点を 50 × 50 の格子上で動径基底関数（RBF）で補間し、等高線図（`scan2d_map.png`）と 3D 曲面（`scan2d_landscape.html`）を描きます。

---

## surface.csv の読み方と判定

`surface.csv` には格子点ごとの行と、基準の行が 1 つ入ります。

| 列 | 内容 |
| --- | --- |
| `i`, `j` | 格子の番号。開始構造に最も近い値が 0 なので、値の昇順ではなく計算した順の番号 |
| `d1_A`, `d2_A`（`q1`, `q2` も同じ値） | 緩和の後に測った座標の値。どの軸でも列名は `_A` のままで、角度の軸には度が入る。単位は `q1_unit`, `q2_unit`（`angstrom` か `degree`） |
| `target_d1_A`, `target_d2_A`（`target_q1`, `target_q2` も同じ値） | その点の拘束の目標値 |
| `energy_hartree` | 拘束を外したエネルギー（Hartree） |
| `energy_kcal` | 基準からの相対エネルギー（kcal/mol） |
| `bias_converged` | 拘束付きの緩和が収束したか |
| `is_preopt` | 基準の行だけ `true` |
| `d1_label`, `d2_label` | 図に使う軸の名前 |

* **基準の行**: `i = j = -1`、`is_preopt = true` の行は開始構造です。表には残りますが、格子点・エネルギーの基準・図の点には使いません。
* **使える点**: 緩和が収束し、エネルギーが有限で、構造ファイルが書けた点を「使える点」とします。エネルギーの基準と図には、使える点だけを使います。
* **判定**: `result.json`（`--out-json`）の `scientific_status` は、すべての格子点が使える点なら `success`、一部だけなら `partial`（{ref}`終了コード <ja-exit-codes>` 0）、1 つも無ければ `failed`（終了コード 1）です。点の数は `n_points_attempted` と `n_points_usable` に入ります。
* **次の段階**: 図は補間なので、鞍点に近い計算点の構造 `grid/point_*.pdb` を [`tsopt`](tsopt.md) に渡します。2 つの谷の点は [`path-search`](path-search.md) の入力にできます。

---

## 主な出力ファイル

`--out-dir` に次のファイルを書きます。

```text
result_scan2d/
├─ surface.csv                   # 基準の行を含む格子の表
├─ scan2d_map.png                # 2D 等高線図
├─ scan2d_landscape.html         # 3D 曲面（ブラウザで開く）
├─ grid/
│  ├─ point_i150_j090.xyz        # 各格子点の緩和後の構造
│  ├─ preopt_iDDD_jDDD.xyz       # 開始構造（基準の行）
│  └─ inner_path_d1_000_trj.xyz  # d₁ の値ごとの内側ループの軌跡（--dump 指定時）
└─ result.json                   # 結果の要約（--out-json 指定時）。summary.json も同じ内容
```

まず `surface.csv` と 2 つの図を確認し、各点の構造は `grid/` を見てください。`result.json` の `grid_points[]` には、各格子点の番号・値・目標値・エネルギー・収束の可否・構造ファイルが入ります。ファイル名から値を読み取らずに、この対応を使ってください。

* **ファイル名**: `i`・`j` の後の数字（タグ `DDD`）は目標値の 100 倍（Å、角度では度）を 3 桁以上に 0 で埋めた数で、`surface.csv` の格子の番号ではありません。`d1 = 1.50 Å, d2 = 0.90 Å` なら `point_i150_j090.xyz`、角度 120° なら `12000` です。丸めたタグが別の点と重なると、後のファイル名には 0 始まりの番号 `_grid_III_JJJ` が付きます。
* **ほかの形式**: 各構造は `.pdb` でも書き、B-factor の欄に各原子の層を入れます：ML 領域 0、可動 MM 10、凍結 MM 20。`--no-convert-files` で止められます。{ref}`mmCIF の入力 <ja-mmcif-input>` と、PDB の欄に入りきらない大きな PDB の入力では、元の識別子を保った `.cif` も書きます。

---

## 主な CLI オプション

ML/MM の計算コマンドに共通のオプションは {ref}`ML/MM の共通オプション <ja-mlmm-options>` に 1 か所でまとめてあります。下の表は `scan2d` に固有のものだけです。

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | 全系の構造ファイル（`.pdb`, `.cif`, `.mmcif`、または `--ref-pdb` を付けた `.xyz`） |
| `-q, --charge` | 整数 | `None` | ML 領域の電荷。`-l` を使う場合のほかは必須 |
| `-m, --multiplicity` | 整数 | `1` | ML 領域のスピン多重度（2S+1） |
| `-l, --ligand-charge` | 文字列 | `None` | 未知のリガンド残基の総電荷（例: `-1`）または残基名ごとの電荷（例: `'GPP:-3,SAM:1'`）。`-q` を省いたときに ML 領域の電荷を求めるのに使用（PDB/mmCIF 入力または `--ref-pdb`） |
| `-s, --scan-lists` | 文字列 | （必須） | YAML/JSON スペックファイルまたは 1 つのインラインリテラルで 2 つの範囲を指定。距離 `(i,j,low,high)`、角度 `(i,j,k,low,high)`、二面角 `(i,j,k,l,low,high)` |
| `-o, --out-dir` | パス | `./result_scan2d/` | 出力先ディレクトリ |
| `--max-step-size` | 浮動小数点数 | `0.2` | 距離の軸の格子間隔の上限（Å） |
| `--max-angle-step-size` | 浮動小数点数 | `5.0` | 角度の軸の格子間隔の上限（度） |
| `--max-dihedral-step-size` | 浮動小数点数 | `10.0` | 二面角の軸の格子間隔の上限（度） |
| `--restraint-k` | 浮動小数点数 | `300.0` | 拘束の強さ k（距離は eV/Å²、角度は eV/rad²）。別名 `--bias-k` |
| `--preopt/--no-preopt` | フラグ | `False` | スキャンの前に入力構造を拘束なしで最適化 |
| `--opt-mode` | `grad` / `hess` | `grad` | 各点の緩和の方法：L-BFGS / RFO |
| `--dump/--no-dump` | フラグ | `False` | d₁ の値ごとの内側ループの軌跡を `grid/` に出力 |
| `--baseline` | `min` / `first` | `min` | `energy_kcal` の 0 点：使える点の最小値、または点 `(0, 0)` |
| `--zmin`, `--zmax` | 浮動小数点数 | 曲面の最小値 / 最大値 | カラースケールの下限と上限（kcal/mol） |
| `--thresh` | 文字列 | `baker` | 各緩和の収束プリセット（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） |
| `--out-json/--no-out-json` | フラグ | `False` | 結果の要約を `result.json` に出力（[JSON 出力の一覧](json-output.md)） |

全オプションは [自動生成のオプションの一覧（英語のみ）](../reference/commands/scan2d.md) を参照してください。

> **補足:** YAML（`--config`）では、`--opt-mode grad` のとき `lbfgs` か `opt.lbfgs` の節が、`--opt-mode hess` のとき `rfo` か `opt.rfo` の節が効きます。軌跡を書くかどうかは `--dump` だけで決まり、YAML の `opt.dump` は使いません。

---

## 使用上の注意点

* **範囲は 1 つのリテラルに 2 つ**: `-s` には、1 つのインラインリテラル、または YAML/JSON ファイルの `pairs:` で、ちょうど 2 つの範囲を渡します。複数ステージのスキャンには [`scan`](scan.md) を使ってください。
* **chain のある PDB**: {ref}`位置固定の 4 項目の形 <ja-scan-list-spec>` `A:SAM:320:CS1` を使うと原子を一意に指定できます。
* **拘束の単位**: `--restraint-k` の単位は、距離では eV/Å²、角度と二面角では eV/rad² です。オプションを省くと YAML の `bias.k` が使われます。
* **拘束を外したエネルギー**: 記録するエネルギーはすべて拘束を外して計算した値なので、`surface.csv` をそのままフィッティングや作図のスクリプトに使えます。
* **エネルギーの基準**: `--baseline min`（デフォルト）は使える点の最小値を 0 にします。`--baseline first` は点 `(i, j) = (0, 0)` を 0 にし、`(0, 0)` が使える点でなければ使える点の最小値を使います。
* **凍結原子**: `--freeze-atoms` か YAML の `geom.freeze_atoms` で指定した原子と、凍結 MM 層の原子は、どの緩和でも固定されます。スキャンする座標の原子がすべて {ref}`凍結原子 <ja-freeze-atoms-and-restraints>` だとエラーになります。
* **計算せずに指定を確かめる**: `--dry-run` は入力・電荷とスピン・`-s` を読み、計画を表示して、最適化をせずに終わります。
* **サイクル数の上限**: `--relax-max-cycles`（デフォルト `100000`）が各緩和のサイクル数を制限します。指定すると YAML の `opt.max_cycles` より優先されます。
* **使える点が少ないとき**: 使える点が 3 つ未満か、すべて 1 本の直線上にあるときは、図だけを省きます。`[plot] NOTE: Plots skipped: …` を表示し、`surface.csv`（`--out-json` 指定時は `result.json` も）を書いて、終了コード 0 で終わります。使える点が 1 つも無いときは `[plot] No finite data for plotting.` を表示し、終了コード 1 で終わります。
* **PNG の書き出し**: PNG は Plotly と Kaleido で書き出します。書き出しに失敗したときは `[plot] NOTE: PNG export skipped: …` を表示して HTML の曲面は書き、`result.json` には PNG を載せません。

---

## 関連ドキュメント

* [scan](scan.md) — 1 つの構造からの、1 つ以上の座標の段階的スキャン
* [scan3d](scan3d.md) — 3 つの座標のエネルギー格子
* [opt](opt.md) — スキャンの前後の単一構造の最適化
* [tsopt](tsopt.md) — 鞍点に近い構造からの TS 最適化
* [path-search](path-search.md) — マップから取った構造を通る MEP 探索
* [all](all.md) — 一貫ワークフロー
* [トラブルシューティング](troubleshooting.md) — 異常終了時の原因切り分けと対処法
