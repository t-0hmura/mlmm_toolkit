# `oniom-import`（ONIOM 入力から XYZ と層付き PDB へ）

## 概要

`oniom-import` サブコマンドは、**Gaussian ONIOM または ORCA QM/MM の入力ファイルから、XYZ ファイルと層付き PDB を作り直します**。

### 主な用途

* **mlmm の外で直した入力**: 自分で作った、または直した Gaussian・ORCA の QM/MM 入力を mlmm に戻す
* **名前つきの往復**: [`oniom-export`](oniom-export.md) で書き出した入力を、`--ref-pdb` で元の PDB の原子名・残基名つきで戻す

---

## 基本的な実行例

### 1. ORCA の入力

ORCA QM/MM の入力から構造を作り直します。

```bash
mlmm oniom-import -i ts_guess.inp -o ts_guess_imported
```

端末に原子数と層ごとの件数が出て、続いて `ts_guess_imported.xyz` と `ts_guess_imported_layered.pdb` の 2 行の `[oniom-import] wrote:` が出ます。

### 2. Gaussian の入力（拡張子からモードを決める）

拡張子が `.gjf` と `.com` なら Gaussian モード、`.inp` なら ORCA モードです。

```bash
mlmm oniom-import -i model.gjf -o model_imported
```

### 3. モードを明示する

`--mode` は拡張子より優先され、ファイル名が `.gjf`・`.com`・`.inp` のどれでもないときに必要です。

```bash
mlmm oniom-import -i model.inp --mode orca -o model_imported
```

### 4. 参照 PDB の名前を残す

`oniom-export` に渡した PDB の原子名・残基名を、読み込んだ座標に書き写します。

```bash
mlmm oniom-import -i model.inp --ref-pdb complex_layered.pdb -o model_imported
```

ログの `[oniom-import] ref_order=identity-verified` は、原子の並びが `oniom-export` の書いたマーカーと一致したことを示します。

---

## 処理の仕組みと計算仕様

1. **モード**:
`--mode`、無ければ `-i` の拡張子で決めます。
2. **座標と層**:
Gaussian モードでは、座標の各行を `<atom> <0|-1> x y z H|L` として読みます。`H` の原子は QM、`L` で `0` の原子は可動 MM、`L` で `-1` の原子は固定 MM です。ORCA モードでは、`%qmmm` ブロックの `QMAtoms` と `ActiveAtoms`、`* xyz` ブロックの座標を読みます。QM 原子は常に可動として数えます。
3. **QM の電荷と多重度**:
Gaussian では、6 つの整数が並ぶ電荷と多重度の行の 3 番目と 4 番目（高レベルの QM 領域）、ORCA では `* xyz <charge> <multiplicity>` の行から読みます。どちらも XYZ のコメント行に `q=<charge> m=<multiplicity>` として書きます。
4. **XYZ**:
入力のすべての原子を `<out_prefix>.xyz` に書きます。
5. **層付き PDB**:
B-factor を 0・10・20 にして `<out_prefix>_layered.pdb` を書きます。`--ref-pdb` が無いときは、各原子に元素名をつけ、chain A の残基 `MOL` 1 にまとめます。`--ref-pdb` があるときは参照 PDB の行を残し、座標と B-factor だけを書き換えます。その前に原子の並びを確かめます。元素が 1 原子ずつ一致することに加えて、`oniom-export` のマーカーがある入力はハッシュの一致（`identity-verified`）、マーカーの無い入力はどの元素も 1 回しか現れないこと（`element-verified`）が条件で、どちらでもなければ `--allow-unverified-ref-order` が必要です（`unverified-opt-in`）。

---

## 主な出力ファイル

* **`<out_prefix>.xyz`**: すべての原子の座標。コメント行は `mode=<mode> atoms=N qm=… movable=… q=<charge> m=<multiplicity>` です。
* **`<out_prefix>_layered.pdb`**: 同じ座標に、層を B-factor で持たせた PDB です。
* **ログ**: `[oniom-import] mode=…`、件数の `atoms=N, qm=…, movable=…, frozen=…`、2 行の `wrote:` が出ます。`--ref-pdb` 指定時は `ref_order=identity-verified` / `element-verified` / `unverified-opt-in` のどれかも出ます。

---

## 主な CLI オプション

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | ONIOM の入力。`.gjf` / `.com`（g16）または `.inp`（ORCA） |
| `--mode` | `g16` / `orca` | `-i` の拡張子から | 入力の形式 |
| `-o, --out-prefix` | パス | カレントディレクトリに入力の stem | 出力ファイルの接頭辞 |
| `--ref-pdb` | パス | `None` | 原子名・残基名を出力に書き写す PDB。同じ原子を同じ順に並べたもの |
| `--allow-unverified-ref-order/--no-allow-unverified-ref-order` | フラグ | `False` | マーカーでも元素の一意性でも並びを確かめられないときに、`--ref-pdb` を位置で対応づける |

全オプションの一覧は [自動生成 CLI リファレンス](../reference/commands/oniom_import.md) を参照してください。

---

## 使用上の注意点

* **読むのは入力ファイル**: `oniom-import` が読むのは入力ファイルで、Gaussian・ORCA の出力から最適化構造を取り出す機能ではありません。mlmm で計算を続けるには、外部プログラムで得た final geometry をトポロジーと同じ原子順で書き出し、元の parm7 と ML 領域と組み合わせて使ってください。
* **読める書式**: Gaussian の入力は、`oniom-export` が書く形の、可動フラグの列を持つ 2 層 ONIOM の座標行が必要です。ORCA の入力は、`%qmmm` の中にそれぞれ 1 行で書いた `QMAtoms {…} end` と `ActiveAtoms {…} end`、および `* xyz <charge> <multiplicity>` ブロックが必要です。ほかの書式はエラーで止まります。
* **マーカーの照合**: `oniom-export` が書いたマーカーは厳しく照合します。壊れたマーカー、重複したマーカー、参照 PDB と一致しないハッシュは、`--allow-unverified-ref-order` を付けても止まります。
* **`--allow-unverified-ref-order`**: マーカーが無く、同じ元素が複数あって並びを確かめられない入力で、並びを自分で確かめたときだけ使ってください。`--ref-pdb` と一緒に指定する必要があります。
* **終了コード**: {ref}`終了コード <ja-exit-codes>`を参照してください。

---

## 関連ドキュメント

* [oniom-export](oniom-export.md) — Gaussian ONIOM・ORCA QM/MM の入力の書き出し
* [define-layer](define-layer.md) — 層の B-factor
* [トラブルシューティング](troubleshooting.md) — 実行に失敗したときの対処
