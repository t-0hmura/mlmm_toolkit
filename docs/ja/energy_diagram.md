# `energy-diagram`

## 概要

> **要約:** 数値エネルギーだけを入力として状態エネルギーダイアグラムを描画します（構造ファイル不要、ML/MM 計算なし）。

`mlmm energy-diagram` は指定された数値のみを可視化します。PDB/XYZ 構造ファイルの読み込みは行わず、熱化学計算（`--thermo`）や DFT 計算（`--dft`）のステップも実行しません。

## 最小例

```bash
mlmm energy-diagram -i 0 12.5 4.3 -o energy.png
```

## 出力の見方

- `OUTPUT.(png|jpg|jpeg|svg|pdf)` -- 描画されたエネルギーダイアグラム画像

## よくある例

1. リスト文字列で入力。

```bash
mlmm energy-diagram -i "[-205.1, -190.4, -198.7]" -o energy.png
```

2. X/Y ラベルを指定。

```bash
mlmm energy-diagram -i 0 12.5 4.3 --label-x R TS P --label-y "ΔE (kcal/mol)" -o energy.png
```

## ワークフロー
1. `-i/--input` から値を収集します（繰り返し指定、1 フラグ後の複数値、リスト文字列に対応）。
2. 全値を float として解釈し、2 点未満なら早期にエラーを返します。
3. 任意の `--label-x` を解釈します。未指定時は `S1`, `S2`,... を自動生成します。
4. `--label-x` の個数と値の個数の一致を検証し、図を描画します。
5. `-o/--output` に画像を保存し、保存先パスを表示します。

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input TEXT...` | 数値入力（複数引数またはリスト形式文字列） | 必須 |
| `-o, --output PATH` | 出力画像パス（`.png/.jpg/.jpeg/.svg/.pdf`） | `energy_diagram.png` |
| `--label-x TEXT...` | X 軸状態ラベル（入力値と同じ個数が必要） | `S1, S2,...` |
| `--label-y TEXT` | Y 軸ラベル | `ΔE (kcal/mol)` |

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [trj2fig](trj2fig.md) -- 軌跡エネルギーからプロファイルを描画
- [all](all.md) -- エネルギーダイアグラム出力を含むエンドツーエンド実行
