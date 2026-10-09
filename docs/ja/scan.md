# `scan`（拘束付き座標スキャン）

`scan` サブコマンドは、層付き酵素構造の中の距離・角度・二面角を調和拘束で少しずつ動かし、各点でそれ以外の自由度を ML/MM 計算機で緩和して、1 つの構造から反応経路の候補を作ります。1 つのリテラル（または YAML の 1 つのステージ）に書いた座標は 1 つの**ステージ**として一緒に動きます。リテラルを複数並べるとステージが順に実行され、各ステージは前のステージの緩和後の構造から始まります。

---

## 主な用途

* **1 つの構造からの経路づくり**: 反応物の反応する結合を動かして、中間体や生成物に近い構造を作り、[`path-search`](path-search.md) に渡す
* **反応の順序の検討**: 結合形成とプロトン移動を 1 つのステージで動かす場合と、別のステージに分ける場合とで、エネルギーの変化を比べる
* **`all` のスキャン段の単独実行**: [`all`](all.md) が `-s` で行うスキャンを、刻み幅や拘束を変えて単独で実行し直す
* **{ref}`結果の判定 <ja-scan-checking-result>`**: 各ステージで共有結合ができたか切れたかが出力され、`result.json` には `scientific_status` が入る

ML 領域の計算バックエンドにはデフォルトの **UMA**（Meta）のほか、`-b/--backend` オプションで **ORB**、**MACE**、**AIMNet2**、DFT（`dft`）も選択可能です。独立した 2 つまたは 3 つの座標でエネルギーの格子を作るには、[`scan2d`](scan2d.md) または [`scan3d`](scan3d.md) を使います。

---

## 基本的な実行例

例の `pocket.pdb` は `real.parm7` に対応する全系の構造で、`ml_region.pdb` はそのうちの ML 領域（リンク水素なし）を選びます。

原子は同梱の酵素の例（`examples/beza/1.R.pdb`）のもので、この PDB は chain の欄が空です。そのため、原子は chain を省いた 3 項目（残基名・残基番号・原子名）を任意の順序で、カンマか空白で区切って書きます。

### 1. YAML スペックファイルからの実行

ステージをファイルに書き、`--out-json` を付けて結果の要約も出力します。

```yaml
# scan.yaml
stages:
  - [["SAM,320,CS1", "GPP,321,C7", 1.60]]
  - [["GPP,321,H11", "GLU,186,OE2", 0.90]]
```

```bash
mlmm scan -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s scan.yaml --out-json -o ./result_scan
```

端末には各ステージで `[stage k] Covalent-bond changes (start vs final): Yes`（または `No`）が出て、最後に全ステージの `Summary` と `====== Scan finished ======` が出ます。`result_scan/result.json` の `scientific_status` には、すべてのステージが収束すると `success` が入ります。

### 2. インラインリテラルでの指定

単純な 1 ステージのスキャンは、コマンドラインに直接書けます。

```bash
mlmm scan -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s '[("SAM,320,CS1","GPP,321,C7",1.60)]'
```

### 3. 2 つの座標を 1 つのステージで動かす

同じリテラルの中の座標は一緒に（協奏的に）動きます。

```bash
mlmm scan -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s '[("CS1 SAM 320","GPP 321 C7",1.60),("GPP 321 H11","GLU 186 OE2",0.90)]' -o ./result_concerted
```

### 4. 2 つのステージを順に実行する

1 つの `-s` の後にリテラルを複数並べます。ステージ 2 はステージ 1 の緩和後の構造から始まります。

```bash
mlmm scan -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s '[("SAM,320,CS1","GPP,321,C7",1.60)]' '[("GPP,321,H11","GLU,186,OE2",0.90)]' -o ./result_staged
```

### 5. 双方向スキャン

[4-tuple](#双方向スキャン4-tuple) を使うと、1 つの距離を入力構造から両方向にスキャンします。

```bash
mlmm scan -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s '[(12, 45, 1.35, 2.50)]'
```

### 6. 軌跡の保存

`--dump` を付けると、各ステップの最適化の軌跡も保存します。

```bash
mlmm scan -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s scan.yaml --dump -o ./result_scan_dump
```

---

## 処理の仕組みと計算仕様

1. **構造の読み込み**:
ML 領域の {ref}`電荷 <ja-charge-specification>` は `-q` または `-l` から決まります。`--preopt` を付けると、まず拘束なしで構造を最適化します。収束しなかった場合は入力構造を使います。
2. **ステージのステップ分割**:
座標ごとに変化量 Δ = 目標値 − 現在値 を求め、ステージを N = ceil(max(|Δ| / h)) ステップに分けます。h は距離では `--max-step-size`（Å）、角度では `--max-angle-step-size`、二面角では `--max-dihedral-step-size`（度）です。各座標は 1 ステップに Δ / N ずつ動くので、ステージ内のすべての座標が同時に目標値に着きます。
3. **拘束付きの緩和**:
各ステップで、調和拘束 E = ½ k (q − q_target)² がスキャンする座標 q をそのステップの目標値に保ち（k は `--restraint-k`）、残りの構造を ML/MM 計算機で L-BFGS（`--opt-mode grad`、デフォルト）または RFO（有理関数最適化、`--opt-mode hess`）により緩和します。凍結 MM 層の原子は動きません。各ステップのエネルギーは、拘束を外して計算した ML/MM のエネルギーを記録します。ML/MM のスキャンでは Cartesian 座標（`geom.coord_type: cart`）がデフォルトで、推奨です。YAML で `dlc` を選ぶこともできますが、収束までにずっと長くかかることがあります。
4. **ステージの終わり**:
`--endopt` を付けると、ステージの最後の構造を拘束なしでもう一度最適化します。そのあと、ステージの最初と最後の構造を比べて共有結合の変化を調べ、ステージの結果を書き出します。
5. **次のステージ**:
次のステージはこの結果から始まります。最後のステージが終わると、全ステージの軌跡を 1 つのファイルにつなぎます。

### 双方向スキャン（4-tuple）

目標値 `(i, j, target)` の代わりに範囲 `(i, j, low, high)` を指定すると、入力構造から両方向にスキャンします。範囲は 2 つのステージに展開されます。

1. **パス 1**: `i`–`j` の距離を現在の値から `low` に向けて動かす。
2. **パス 2**: 入力構造に戻し、`i`–`j` の距離を `high` に向けて動かす。

つないだ軌跡は `low → 入力構造 → high` の順になり、出発構造を通る連続した経路になります。角度の範囲 `(i, j, k, low, high)` と二面角の範囲 `(i, j, k, l, low, high)` も同じようにスキャンします。

(ja-section-bond)=
### 結合変化の検出

両原子の共有結合半径の和に `bond_factor`（デフォルト `1.20`）を掛けた値を T とします。2 原子の距離が T − `margin_fraction` × T 以下なら、結合しているとみなします。結合の形成・切断として報告するのは、距離が `delta_fraction` × T 以上変わった組だけです。`margin_fraction` と `delta_fraction` のデフォルトはどちらも `0.05` です。`path-search` も同じ基準を使います。キーは YAML の [`bond`](yaml-reference.md#bond) の節にあります。

---

(ja-scan-direction-barrier-sign)=
## スキャン方向とバリアの符号

(ja-scan-checking-result)=
### 結果の判定

| 確認する場所 | 見るもの |
| --- | --- |
| 端末（各ステージ） | `[stage k] Covalent-bond changes (start vs final): Yes` とできた結合・切れた結合の一覧、または `No` と `(no covalent changes detected)` |
| 端末（実行の最後） | `Summary`：各ステージの目標値・初期値・座標ごとの刻み・ステップ数・結合変化。続いて `====== Scan finished ======` |
| `result.json`（`--out-json`） | `scientific_status`：すべてのステージの全ステップが収束し（`--preopt` と `--endopt` を付けたときはそれらの最適化も収束し）、エネルギーが有限なら `success`、一部だけなら `partial`、1 つも無ければ `failed` |
| `result.json`（`--out-json`） | `stages[].converged`、`stages[].bond_changes.changed`、`stages[].final_energy_hartree`、各ステップのエネルギー `stages[].energies_hartree` |

`partial` の {ref}`終了コード <ja-exit-codes>` は 0、`failed` は 1 です。収束して狙った結合変化が起きたスキャンは経路の候補になり、エネルギーが最も高いステップは [`tsopt`](tsopt.md) に渡す TS 候補になります。このステップは `scan_trj.xyz` から {ref}`取り出せます <ja-trajectory-one-frame>`。

### バリアの向き

`scan` はエネルギーを記録しますが、バリアは出力しません。**生成物側**から始めたスキャン（またはそこから作った経路や TS 候補）からバリアを読む場合、開始構造との差は**逆方向**のバリア `E(TS) − E(product)` です。順方向のバリアは反応物から計算します。

| 実行内容 | 順方向バリア |
| --- | --- |
| 反応物から始めたスキャン | `E(TS) − E(reactant)`。開始構造との差がそのまま順方向のバリア |
| 生成物から始めたスキャン | `E(TS) − E(reactant)`。開始構造との差では**ない**。E(reactant) は最適化した反応物のエネルギー（例: [`opt`](opt.md) で最適化した IRC の端点） |

これを切り替えるオプションはありません。バリアを引用する前に、スキャンがどちらの端点から始まったかを確認してください。結晶構造の生成物複合体から始めた場合は特に注意してください。

---

## 主な出力ファイル

`--out-dir` に次のファイルを書きます。

```text
result_scan/
├─ preopt/
│  └─ result.xyz                    # 事前最適化した構造（--preopt 指定時）
├─ stage_01/                        # ステージごとのディレクトリ（stage_NN）
│  ├─ result.xyz                    # ステージの final geometry
│  ├─ scan_trj.xyz                  # ステージ内の各ステップの構造とエネルギー
│  └─ scan_s0001_optimization_trj.xyz  # 各ステップの最適化の軌跡（--dump 指定時）
├─ scan_trj.xyz                     # 全ステージをつないだ軌跡
└─ result.json                      # 結果の要約（--out-json 指定時）。summary.json も同じ内容
```

構造と軌跡は同じ名前の PDB（`result.pdb`、`scan.pdb`）でも書きます。`--no-convert-files` で止められます。{ref}`mmCIF の入力 <ja-mmcif-input>` と、PDB の欄に入りきらない大きな PDB の入力では、元の識別子を保った `.cif` も書きます。

* **ステージの結果**: `stage_NN/result.*` はステージ NN の終わりの構造です。[`path-search`](path-search.md) には、開始構造に続けて `stage_NN/result.*` をステージの順に渡します。
* **エネルギーの変化**: `scan_trj.xyz` の各フレームのコメント行には、拘束を外したエネルギー（Hartree）が入っています。[`trj2fig`](trj2fig.md) で図にできます。

---

## 主な CLI オプション

ML/MM の計算コマンドに共通のオプションは {ref}`ML/MM の共通オプション <ja-mlmm-options>` に 1 か所でまとめてあります。下の表は `scan` に固有のものだけです。

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | 全系の構造ファイル（`.pdb`, `.cif`, `.mmcif`、または `--ref-pdb` を付けた `.xyz`） |
| `-q, --charge` | 整数 | `None` | ML 領域の電荷。`-l` を使う場合のほかは必須 |
| `-m, --multiplicity` | 整数 | `1` | ML 領域のスピン多重度（2S+1） |
| `-l, --ligand-charge` | 文字列 | `None` | 未知のリガンド残基の総電荷（例: `-1`）または残基名ごとの電荷（例: `'GPP:-3,SAM:1'`）。`-q` を省いたときに ML 領域の電荷を求めるのに使用（PDB/mmCIF 入力または `--ref-pdb`） |
| `-s, --scan-lists` | 文字列 | （必須） | YAML/JSON スペックファイル、または 1 つ以上のインラインリテラル（1 つが 1 ステージ）。距離の目標値 `(i,j,target)`、または距離 `(i,j,low,high)`・角度 `(i,j,k,low,high)`・二面角 `(i,j,k,l,low,high)` の範囲 |
| `-o, --out-dir` | パス | `./result_scan/` | 出力先ディレクトリ |
| `--one-based/--zero-based` | フラグ | `--one-based` | `-s` の原子インデックスを 1 始まり / 0 始まりとして読む |
| `--max-step-size` | 浮動小数点数 | `0.2` | 1 ステップあたりの距離の最大変化量（Å） |
| `--max-angle-step-size` | 浮動小数点数 | `5.0` | 1 ステップあたりの角度の最大変化量（度） |
| `--max-dihedral-step-size` | 浮動小数点数 | `10.0` | 1 ステップあたりの二面角の最大変化量（度） |
| `--restraint-k` | 浮動小数点数 | `300.0` | 拘束の強さ k（距離は eV/Å²、角度は eV/rad²）。別名 `--bias-k` |
| `--preopt/--no-preopt` | フラグ | `False` | スキャンの前に入力構造を拘束なしで最適化 |
| `--endopt/--no-endopt` | フラグ | `False` | 各ステージの結果を拘束なしで最適化 |
| `--dump/--no-dump` | フラグ | `False` | 各ステップの最適化の軌跡を出力 |
| `--opt-mode` | `grad` / `hess` | `grad` | 緩和の方法：L-BFGS / RFO（`tsopt` では同じ語が別の最適化法を指す。{ref}`コマンドごとの --opt-mode <ja-opt-mode-semantics>` を参照） |
| `--freeze-atoms` | 文字列 | `None` | 凍結する原子の 1 始まりのインデックス（カンマ区切り）。YAML の `geom.freeze_atoms` と凍結 MM 層に加えられる |
| `--out-json/--no-out-json` | フラグ | `False` | 結果の要約を `result.json` に出力（[JSON 出力の一覧](json-output.md)） |

全オプションは [自動生成のオプションの一覧（英語のみ）](../reference/commands/scan.md) を参照してください。

> **補足:** YAML（`--config`）では、`--restraint-k` を指定しないときの拘束の強さを [`bias.k`](yaml-reference.md#bias) で、結合変化の閾値 `bond_factor`・`margin_fraction`・`delta_fraction` を [`bond`](yaml-reference.md#bond) の節で設定できます。

---

## 使用上の注意点

* **`--preopt` は呼び出し方で変わる**: `scan` を単独で実行したときは、`--preopt` を付けない限り事前最適化をしません。`all` の中では、`all --preopt`（デフォルトで有効）に従って事前最適化し、[`all --scan-preopt/--no-scan-preopt`](../reference/commands/all.md) で上書きできます。
* **インラインでは目標値と範囲を混ぜない**: 1 つのインラインリテラルの中でも、1 回の実行のリテラルどうしでも、目標値 `(i,j,target)` と範囲のどちらか一方だけを使います。両方を組み合わせるときは、YAML/JSON スペックの `stages:` に並べてください。
* **範囲を使うときのステージ番号**: 範囲 1 つは `low` 向きと `high` 向きの 2 つのステージになります（4-tuple 1 つなら `stage_01/` と `stage_02/`）。インラインでは、1 つのリテラルの範囲がすべてこの 2 つのステージで一緒に動きます。YAML の `stages:` では、範囲を含むステージの項目がそれぞれ別のステージになり、目標値は 1 つ、範囲は 2 つのステージになります。
* **目標の距離は正の値**にしてください。また、1 つのステージに同じ座標を 2 回書くことはできません。
* **計算せずに指定を確かめる**: `--dry-run` は入力・電荷とスピン・`-s` を読み、ステージの数を表示して、最適化をせずに終了します。
* **凍結原子**: `--freeze-atoms` か YAML の `geom.freeze_atoms` で指定した原子と、凍結 MM 層の原子は、どの緩和でも固定されます。スキャンする座標の原子がすべて {ref}`凍結原子 <ja-freeze-atoms-and-restraints>` だとエラーになります。
* **サイクル数の上限**: `--relax-max-cycles`（デフォルト `100000`）が各緩和のサイクル数を制限します。指定すると YAML の `opt.max_cycles` より優先されます。

---

## 関連ドキュメント

* {ref}`スキャンリスト仕様 <ja-scan-list-spec>` — YAML/JSON スペックファイル、インラインリテラル、原子の指定
* [scan2d](scan2d.md) — 2 つの座標のエネルギーマップ
* [scan3d](scan3d.md) — 3 つの座標のエネルギー格子
* [path-search](path-search.md) — スキャンの結果からの最小エネルギー経路（MEP）探索
* [all](all.md) — 1 つの構造と `-s` からのスキャンを含む一貫ワークフロー
* [トラブルシューティング](troubleshooting.md) — 異常終了時の原因切り分けと対処法
