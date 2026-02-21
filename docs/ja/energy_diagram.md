# `energy-diagram`

## 概要
`energy-diagram` は、数値だけを入力して状態エネルギーダイアグラムを描画します。

- 構造ファイルは不要
- エネルギー計算は実行しない
- `--thermo` / `--dft` は使用しない

`-i/--input` は次のいずれかで指定できます。
- 数値を複数引数で渡す
- リスト形式の数値文字列を1引数で渡す

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

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input TEXT ...` | 数値入力（複数引数またはリスト形式文字列） | 必須 |
| `-o, --output PATH` | 出力画像パス（`.png/.jpg/.jpeg/.svg/.pdf`） | `energy_diagram.png` |
| `--label-x TEXT ...` | X軸状態ラベル（入力値と同じ個数が必要） | `S1, S2, ...` |
| `--label-y TEXT` | Y軸ラベル | `ΔE (kcal/mol)` |

## 注意事項
- 入力順がそのまま描画順になります。
- 入力値は最低2点必要です。

---

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [trj2fig](trj2fig.md) -- 軌跡エネルギーからプロファイルを描画
- [all](all.md) -- エネルギーダイアグラム出力を含むエンドツーエンド実行
