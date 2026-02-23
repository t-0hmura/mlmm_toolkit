# `trj2fig`

## 概要

> **要約:** XYZ 軌跡のコメント行からエネルギーを抽出し、相対または絶対エネルギープロファイルを計算し、Plotly 図と CSV テーブルとしてエクスポートします。

### 概要
- **入力:** 2 行目（コメント）に Hartree エネルギーが格納された XYZ 軌跡。
- **基準モード:** 最初のフレーム（`init`）、基準なし（`None`）、または明示的な 0 始まりフレームインデックス。
- **出力形式:** PNG（デフォルト）、JPEG、HTML、SVG、PDF、CSV。
- **単位:** kcal/mol（デフォルト）または Hartree。
- **X 軸反転:** `--reverse-x` で軸を反転し、最後のフレームが左側に表示されます。

`mlmm trj2fig` は XYZ 軌跡の各フレームのコメント行にエンコードされた Hartree エネルギーを読み取り、kcal/mol または Hartree に変換し、任意で選択したフレームを基準にすべての値を参照し、結果の系列を静的/インタラクティブ図と CSV テーブルとしてエクスポートします。図は太字の目盛り、一貫したフォント、マーカー、平滑化スプライン曲線を使用します（タイトルなし）。

## 使用法

```bash
mlmm trj2fig -i TRAJECTORY.xyz [-o OUTPUTS...] [-r REFERENCE] [--unit {kcal|hartree}] [--reverse-x]
```

### 例

```bash
# デフォルト PNG、最初のフレームを基準とした相対エネルギー
mlmm trj2fig -i traj.xyz

# 基準フレーム #5 の CSV + SVG、Hartree で報告
mlmm trj2fig -i traj.xyz -o energy.csv energy.svg -r 5 --unit hartree

# X 軸反転付きの複数出力を一度に生成
mlmm trj2fig -i traj.xyz -o energy.png energy.html energy.pdf --reverse-x
```

## ワークフロー

1. XYZ 軌跡を解析し、各フレームのコメント行に見つかる最初の浮動小数点数から Hartree エネルギーを抽出します。エネルギーが見つからないか、フレームコメントに解析可能なエネルギーがない場合はエラーが発生します。
2. 基準仕様を正規化します:
 - `init` -- フレーム `0`（`--reverse-x` が有効な場合は最後のフレーム）。
 - `None`/`none`/`null` -- 絶対エネルギー（基準なし）。
 - 整数リテラル -- 対応する 0 始まりフレームインデックス。
3. エネルギーを kcal/mol（デフォルト）または Hartree に変換し、基準が有効な場合は基準値を減算して delta-E を生成します。
4. Plotly 図（太字目盛り、スプライン補間、マーカー、タイトルなし）を構築し、要求されたすべての拡張子にエクスポートします。
5. 任意で `frame`、`energy_hartree`、および要求された単位の適切な delta-E または絶対 E カラムを含む CSV テーブルを出力します。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 2 行目にエネルギーが格納された XYZ 軌跡。 | 必須 |
| `-o, --out PATH` | 繰り返し可能な出力ファイル名。`.png`、`.jpg`/`.jpeg`、`.html`、`.svg`、`.pdf`、`.csv` をサポート。 | `energy.png` |
| _追加引数_ | オプション後に列挙された位置ファイル名。`-o` リストとマージ。 | _None_ |
| `--unit {kcal,hartree}` | プロット/エクスポートされる値のターゲット単位。 | `kcal` |
| `-r, --reference TEXT` | 基準仕様（`init`、`None`、または 0 始まり整数）。 | `init` |
| `--reverse-x` | X 軸を反転し、最後のフレームが左側に表示されます（`init` は最後のフレームになります）。 | `False` |

## 出力

```
<output>.[png|jpg|jpeg|html|svg|pdf] # 要求された拡張子ごとの Plotly エクスポート（デフォルトは energy.png）
<output>.csv # CSV 要求時のオプションエネルギーテーブル
```

- `-o` も位置出力も提供されない場合、カレントディレクトリに `energy.png` が 1 つ書き出されます。
- CSV エクスポートには `frame`、`energy_hartree`、および delta-E カラム（`delta_kcal`/`delta_hartree`）または絶対カラム（基準適用なし時の `energy_kcal`/`energy_hartree`）が含まれます。
- PNG は高解像度のため `scale=2` で Plotly の PNG エクスポートを使用します。

## 注意事項

- エネルギーは各コメントの最初の十進数から取得されます。不正なコメントではエラーが発生します。
- `-o` でサポートされていないファイル拡張子はエラーの原因となります。
- `--reverse-x` は軸方向と `-r init` の動作の両方を反転し、可視化された経路が逆方向に読まれるようにします。
- レガシーの `--output-peak` オプションは削除されました。

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [path_search](path_search.md) -- 再帰的 MEP 探索（trj2fig に適した XYZ 軌跡を生成）
- [irc](irc.md) -- TS からの IRC（エネルギープロファイリング用軌跡を生成）
- [all](all.md) -- エンドツーエンドワークフロー
