# `energy-diagram`

数値エネルギーだけを入力として状態エネルギーダイアグラムを描画します（構造ファイル不要、ML/MM 計算なし）。既にエネルギー値を持っていて構造ベースの ML/MM 計算が不要なときに、数値から直接ダイアグラムを描きます。`mlmm energy-diagram` は与えられた数値を可視化するだけで、PDB/XYZ 構造を読み込まず、熱化学（`--thermo`）や DFT（`--dft`）の処理も実行しません。

## 実行例

相対エネルギーをリスト形式文字列で指定:

```bash
mlmm energy-diagram -i "[0, 12.5, 4.3]" -o energy.png
```

絶対エネルギー（任意の数値）:

```bash
mlmm energy-diagram -i "[-205.1, -190.4, -198.7]" -o energy.png
```

`-i` フラグの繰り返しで値を指定:

```bash
mlmm energy-diagram -i 0 -i 12.5 -i 4.3 -o energy.png
```

X 軸の状態ラベルと Y 軸ラベルを指定:

```bash
mlmm energy-diagram -i "[0, 12.5, 4.3]" --label-x R TS P --label-y "ΔE (kcal/mol)" -o energy.png
```

## 処理の流れ
1. `-i/--input` から値を収集します（繰り返し指定、1 フラグ後の複数値、リスト文字列に対応）。
2. 全値を float として解釈し、2 点未満なら早期にエラーを返します。
3. 任意の `--label-x` を解釈します。未指定時は `S1`, `S2`,... を自動生成します。
4. `--label-x` の個数と値の個数の一致を検証し、図を描画します。
5. `-o/--output` に画像を保存し、保存先パスを表示します。

## 出力

- `OUTPUT.(png|jpg|jpeg|svg|pdf)` -- 描画されたエネルギーダイアグラム画像

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input TEXT...` | 数値入力（複数引数またはリスト形式文字列）。最低 2 点が必要。 | 必須 |
| `-o, --output PATH` | 出力画像パス（`.png/.jpg/.jpeg/.svg/.pdf`） | `energy_diagram.png` |
| `--label-x TEXT...` | X 軸状態ラベル（入力値と同じ個数が必要） | `S1, S2,...` |
| `--label-y TEXT` | Y 軸ラベル | `ΔE (kcal/mol)` |

すべてのフラグの一覧は生成された[コマンドリファレンス](../reference/commands/index.md)を参照してください。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [trj2fig](trj2fig.md) -- 軌跡エネルギーからプロファイルを描画
- [all](all.md) -- エネルギーダイアグラム出力を含む一気通貫実行
