# `define-layer`（ML・MM の層の割り当て）

## 概要

`define-layer` は、全系を ML 領域の周りの 3 つの層に分け、各原子の層を PDB の B-factor の欄に書き込みます。計算コマンドは、この B-factor から層を読み戻します。

| 層 | B-factor | 原子 | 計算での扱い |
| --- | --- | --- | --- |
| ML | 0.0 | ML 領域 | MLIP のエネルギー・力・Hessian |
| Movable-MM | 10.0 | ML 領域から `--movable-cutoff`（既定 8.0 Å）以内の MM 原子 | MM。動ける |
| Frozen-MM | 20.0 | それより遠い MM 原子 | MM。座標は固定し、MM のエネルギーには入る |

### 主な用途

* **計算コマンドの入力を作る**: `opt`・`tsopt`・`freq` などが `--parm7` と一緒に受け取る、層付きの全系の PDB を書き出します。
* **動く殻を変える**: `--movable-cutoff` で Movable-MM の層を広げたり狭めたりします。
* **層の大きさを確かめる**: 計算の前に、各層の原子数を見ます。

`all` に `-c` を付けると、`all` が `define-layer` を実行します。

---

## 基本的な実行例

### 1. model PDB で ML 領域を渡す

全系と、`extract` か `all` が書いた ML 領域の PDB を渡します。

```bash
mlmm define-layer -i system.pdb --model-pdb ml_region.pdb -o labeled.pdb
```

端末の `Layer Summary` の下に各層の原子数が出ます。

### 2. 原子番号で ML 領域を渡す

ML 領域の原子を番号で並べます。

```bash
mlmm define-layer -i system.pdb --model-indices "0,1,2,3,4" --zero-based -o labeled.pdb
```

### 3. 動く殻を広げる

ML 領域から 10.0 Å 以内の MM 残基をすべて動けるようにします。

```bash
mlmm define-layer -i system.pdb --model-pdb ml_region.pdb \
    --movable-cutoff 10.0 -o labeled.pdb
```

---

## 処理の仕組みと計算仕様

1. **ML 領域**: `--model-pdb` の原子を、chain・残基番号・挿入コード・残基名・原子名で入力の原子と照合します。chain が空の原子はほかの欄で照合し、複数の chain の原子に当てはまるとエラーで止まります。`--model-pdb` が無ければ、`--model-indices` に並べた原子が ML 領域になります。
2. **距離**: ML 領域の外の各原子について、いちばん近い ML 原子までの距離を求めます。
3. **割り当て**: ML 原子を持たない残基は、いちばん近い原子の距離で、残基ごと 1 つの層に入ります。`--movable-cutoff` 以内なら Movable-MM、それより遠ければ Frozen-MM です。ML 原子を含む残基では、ML でない原子を 1 つずつ距離で割り当てます。
4. **出力**: 入力の B-factor の欄だけを 0・10・20 に書き換えて書き出し、Layer Summary を表示します。

---

## 主な出力ファイル

```text
./
├─ <input>_layered.pdb   # -o が無いとき、入力と同じディレクトリ
└─ <input>_layered.cif   # mmCIF の入力か、PDB の桁に収まらない PDB の入力のとき
```

PDB は入力のレコードをすべて保ち、B-factor だけを変えます。`-o` の名前が `.cif` か `.mmcif` で終わっても、拡張子を `.pdb` にした名前の PDB で書きます。端末には `Layer Summary` の下に、`Layer 1 (ML, B=0):`・`Layer 2 (Movable MM, B=10):`・`Layer 3 (Frozen MM, B=20):` の行で各層の原子数が出て、続いて `Total atoms:` が出ます。ビューアで出力を B-factor で色分けすると、3 つの層が見えます。

---

## 主な CLI オプション

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | 全系の PDB か mmCIF |
| `--model-pdb` | パス | `None` | ML 領域の原子の PDB か mmCIF |
| `--model-indices` | 文字列 | `None` | ML 原子の番号（`'1,2,3,4'`、`'1-10,15,20-25'` など）。既定は 1 始まりで、`--zero-based` で 0 始まり。`--model-pdb` が無いときに使う |
| `--movable-cutoff` | 浮動小数点数 | `8.0` | ML 領域からこの距離（Å）以内の MM 原子を Movable-MM にする。遠い原子は Frozen-MM |
| `-o, --output` | パス | `<input>_layered.pdb` | 出力 PDB |
| `--one-based/--zero-based` | フラグ | `--one-based` | `--model-indices` の読み方 |

全オプションの一覧は [自動生成 CLI リファレンス](../reference/commands/define_layer.md) を参照してください。

---

## 使用上の注意点

* **ML 領域の指定は必須**: `--model-pdb` も `--model-indices` も無いと、`ERROR: Either --model-pdb or --model-indices must be provided.` を出して終了コード 2 で止まります。
* **`--model-pdb` が優先**: 両方を与えると、計算コマンドと同じく `--model-pdb` が使われます。
* **トポロジーと対応した PDB を使う**: `-i` には `mm-parm` が書く PDB を渡し、層付きの PDB の原子が `parm7` と同じ順に並ぶようにしてください（[mm-parm の例 4](mm-parm.md#基本的な実行例)）。入力に無い原子が `--model-pdb` にあると、エラーで止まります。
* **マルチ MODEL の入力**は、最初の MODEL だけを使い、警告を出します。
* **しきい値の選び方**: `--movable-cutoff` を小さくすると計算は軽くなり、大きくすると周りの環境がより緩和できます。計算コマンドに `--movable-cutoff` を渡すと、B-factor の層の代わりにその距離で層を決めます。詳しくは [ML 領域と層の組み方](model-setup.md) を見てください。

---

## 関連ドキュメント

* [ML 領域と層の組み方](model-setup.md) — ML 領域、動く殻、Hessian の範囲を決める
* [extract](extract.md) — `--model-pdb` に渡す ML 領域を切り出す
* [mm-parm](mm-parm.md) — トポロジーと、層を付ける対応した PDB を作る
* [all](all.md) — 一括のワークフロー。`-c` で `define-layer` を実行する
* [opt](opt.md) — 層付きの系を最適化する
* [トラブルシューティング](troubleshooting.md) — 層と原子の順のエラー
