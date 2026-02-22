# `energy-diagram`

## 概要

> **要約:** 数値エネルギーだけを入力として状態エネルギーダイアグラムを描画します（構造ファイル不要、ML/MM 計算なし）。

### クイックリファレンス
- **入力:** `-i/--input` で与える数値列（複数トークンまたはリスト形式文字列）。
- **出力:** 画像ファイル 1 つ（`.png` / `.jpg` / `.jpeg` / `.svg` / `.pdf`）。
- **デフォルト出力:** `energy_diagram.png`。
- **状態ラベル:** `--label-x` で任意指定。省略時は `S1`, `S2`, ...。
- **用途:** 既にエネルギー値があり、図だけ作成したいとき。

`mlmm energy-diagram` は与えた数値を可視化するだけで、PDB/XYZ を読み込まず、`--thermo` / `--dft` のような計算も実行しません。

## 使用法
```bash
mlmm energy-diagram -i VALUES... [-o OUTPUT] [--label-x ...] [--label-y ...]
```

## 例
```bash
# 数値を複数引数で指定
mlmm energy-diagram -i 0 12.5 4.3 -o energy.png

# リスト文字列で指定
mlmm energy-diagram -i "[-205.1, -190.4, -198.7]" -o energy.png

# X/Yラベルを指定
mlmm energy-diagram -i 0 12.5 4.3 --label-x R TS P --label-y "ΔE (kcal/mol)" -o energy.png
```

## ワークフロー
1. `-i/--input` から値を収集します（繰り返し指定、1 フラグ後の複数値、リスト文字列に対応）。
2. 全値を float として解釈し、2 点未満なら早期にエラーを返します。
3. 任意の `--label-x` を解釈します。未指定時は `S1`, `S2`, ... を自動生成します。
4. `--label-x` の個数と値の個数の一致を検証し、図を描画します。
5. `-o/--output` に画像を保存し、保存先パスを表示します。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input TEXT ...` | 数値入力（複数引数またはリスト形式文字列） | 必須 |
| `-o, --output PATH` | 出力画像パス（`.png/.jpg/.jpeg/.svg/.pdf`） | `energy_diagram.png` |
| `--label-x TEXT ...` | X軸状態ラベル（入力値と同じ個数が必要） | `S1, S2, ...` |
| `--label-y TEXT` | Y軸ラベル | `ΔE (kcal/mol)` |

## 出力
```
OUTPUT.(png|jpg|jpeg|svg|pdf)
```
- `-o/--output` を省略した場合、カレントディレクトリに `energy_diagram.png` を出力します。
- 出力拡張子がない場合は `.png` が自動で補完されます。
- 必要なら親ディレクトリを自動作成します。

## 注意事項
- 入力順がそのまま描画順になります。
- 入力値は最低2点必要です。
- 構造ファイルの読み込みやエネルギー計算は行いません。

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [trj2fig](trj2fig.md) -- 軌跡エネルギーからプロファイルを描画
- [all](all.md) -- エネルギーダイアグラム出力を含むエンドツーエンド実行
