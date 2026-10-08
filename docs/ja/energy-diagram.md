# `energy-diagram`（状態エネルギー図）

`energy-diagram` サブコマンドは、与えた数値から**状態エネルギー図を描きます**。構造ファイルを読まず、ML/MM の計算も行いません。`all` や `path-search` が書き出す `summary.json` の `energy_diagrams`（[JSON 出力リファレンス](json-output.md)）や論文の表など、すでに手元にあるエネルギーの作図に向いています。

---

## 主な用途

* **既知のエネルギーの作図**: ワークフローや DFT で得た反応物（R）・遷移状態（TS）・中間体（IM）・生成物（P）のエネルギー
* **論文やスライド用の図**: SVG・PDF のベクター形式での出力
* **手早い確認**: 作図のコードを書かずに、少数の値を図にする

---

## 基本的な実行例

### 1. 値を 1 つのリストで渡す

すべての値を、引用符で囲んだ 1 つのリストとして渡します。

```bash
mlmm energy-diagram -i "[0, 12.5, 4.3]" -o energy.png --out-json
```

端末に `[energy-diagram] Saved -> …/energy.png` が出て、画像と同じ場所の `result.json` に `n_points: 3` があれば成功です。

### 2. 値ごとに `-i` を付ける

値の数だけ `-i` を繰り返します。

```bash
mlmm energy-diagram -i 0 -i 12.5 -i 4.3 -o energy.png
```

### 3. 状態と軸のラベル

x 軸の状態に名前を付け、y 軸のラベルを指定します。

```bash
mlmm energy-diagram -i "[0, 12.5, 4.3]" \
  --label-x "['R','TS','P']" --label-y "ΔE (kcal/mol)" -o energy.png
```

---

## 処理の仕組みと計算仕様

1. **値の読み込み**:
`-i` から値を読みます。値ごとに `-i` を繰り返すか、リスト形式の文字列 1 つ（`"[0, 12.5, 4.3]"` または `"0, 12.5, 4.3"`）で渡します。数値はそのまま描くので、相対エネルギーでも絶対エネルギーでもかまいません。
2. **ラベル**:
`--label-x` で状態ごとのラベルを、繰り返しかリスト形式の文字列 1 つで与えます。省略すると `S1`、`S2`、… になります。
3. **作図**:
各状態をそのエネルギーの高さの短い横棒で描き、隣り合う横棒を点線で結び、最初の状態のエネルギーの高さに薄い灰色の点線を引きます。
4. **保存**:
`-o` の拡張子で形式が決まります。拡張子の無いパスには `.png` を付け、親ディレクトリが無ければ作ります。

---

## 主な出力ファイル

```text
energy_diagram.png   # 状態エネルギー図（デフォルト名。-o で指定）
result.json          # execution_status、scientific_status、n_points、files（--out-json 指定時）
summary.json         # result.json の写し。result.json を読む（--out-json 指定時）
```

`result.json` と `summary.json` は画像と同じディレクトリに書き出します。記録するのは点の数と画像のパスで、値とラベルは含みません。

---

## 主な CLI オプション

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | 文字列 | （必須） | エネルギーの値。値ごとに `-i` を繰り返すか、リスト形式の文字列 1 つで指定 |
| `-o, --output` | パス | `energy_diagram.png` | 出力画像（`.png`, `.jpg`, `.jpeg`, `.svg`, `.pdf`） |
| `--label-x` | 文字列 | `S1, S2, …` | x 軸の状態ラベル。状態ごとに繰り返すか、リスト形式の文字列 1 つで指定 |
| `--label-y` | 文字列 | `ΔE (kcal/mol)` | y 軸のラベル |
| `--out-json/--no-out-json` | フラグ | `False` | 画像の隣に `result.json` と `summary.json` を出力 |

全オプションの一覧は [自動生成 CLI リファレンス](../reference/commands/energy_diagram.md) を参照してください。

---

## 使用上の注意点

* **値は 2 つ以上**: 1 つ以下では `Provide at least two numeric values with -i/--input.` で止まります。
* **ラベルの数**: `--label-x` のラベルの数は値の数と同じにしてください。
* **順序**: 入力の順がそのまま x 軸の順になります。
* **単位**: 値はそのまま描くので、単位は `--label-y` に書いてください。
* **終了コード**: {ref}`終了コード <ja-exit-codes>`を参照してください。

---

## 関連ドキュメント

* [trj2fig](trj2fig.md) — 軌跡の各フレームからエネルギープロファイルを作図
* [all](all.md) — エネルギー図も自動で描く全工程のワークフロー
* [JSON 出力リファレンス](json-output.md#energy-diagram) — `result.json` の欄
* [トラブルシューティング](troubleshooting.md) — 実行に失敗したときの対処。画像の書き出しの失敗は {ref}`図のエクスポートが失敗する <ja-plot-export-fails-chrome-missing>`
