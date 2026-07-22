# `trj2fig`

XYZ 軌跡の各フレームのコメント行に格納された Hartree エネルギーを抽出し、kcal/mol または Hartree に変換します。任意のフレームを基準にして、静的・インタラクティブ図と CSV テーブルを出力できます。`-q/--charge` または `-m/--multiplicity` を指定すると、代わりに選択した MLIP バックエンドで全フレームを再計算します。この再計算は各フレームの MLIP による直接再評価であり、ML/MM ONIOM エネルギーではありません。図は太字目盛り、統一した書体、マーカー、スプライン平滑曲線を使用します（タイトルなし）。

## 実行例

デフォルト PNG、最初のフレームを基準とした相対エネルギー:

```bash
# デフォルト PNG、最初のフレームを基準とした相対エネルギー
mlmm trj2fig -i traj.xyz
```

基準フレーム #5 の CSV + SVG、Hartree で報告:

```bash
# 基準フレーム #5 の CSV + SVG、Hartree で報告
mlmm trj2fig -i traj.xyz -o energy.csv energy.svg -r 5 --unit hartree
```

X 軸反転付きの複数出力を一度に生成:

```bash
# X 軸反転付きの複数出力を一度に生成
mlmm trj2fig -i traj.xyz -o energy.png energy.html energy.pdf --reverse-x
```

バックエンド設定を明示して全フレームを再計算し、来歴を JSON に保存:

```bash
mlmm trj2fig -i traj.xyz -q 0 -m 1 -b uma --backend-model uma-s-1p2 \
    --precision fp32 -o energy.png energy.csv --out-json
```

## 処理の流れ

1. XYZ 軌跡を解析します。`-q/--charge` と `-m/--multiplicity` のどちらも指定しない場合は、各フレームのコメント行から Hartree エネルギーを抽出します。いずれかを指定すると選択したバックエンドで全フレームを再計算し、省略した値は電荷 0、スピン多重度 1 として解決されます。
2. 基準仕様を正規化します:
 - `init` -- フレーム `0`（`--reverse-x` が有効な場合は最後のフレーム）。
 - `None`/`none`/`null` -- 絶対エネルギー（基準なし）。
 - 整数リテラル -- 対応する 0 始まりフレームインデックス。
3. エネルギーを kcal/mol（デフォルト）または Hartree に変換し、基準が有効な場合は基準値を減算して delta-E を生成します。
4. Plotly 図（太字目盛り、スプライン補間、マーカー、タイトルなし）を構築し、要求されたすべての拡張子にエクスポートします。
5. 任意で `frame`、`energy_hartree`、および要求された単位の適切な delta-E または絶対 E カラムを含む CSV テーブルを出力します（カラム構成は「出力」を参照）。

## 出力

```
<output>.[png|jpg|jpeg|html|svg|pdf] # 要求された拡張子ごとの Plotly エクスポート（デフォルトは energy.png）
<output>.csv # CSV 要求時のオプションエネルギーテーブル
result.json # --out-json 指定時の機械可読な結果とエネルギーの来歴
summary.json # --out-json 指定時の同一内容の機械可読ミラー
```

- `-o` も位置出力も提供されない場合、カレントディレクトリに `energy.png` が 1 つ書き出されます。
- CSV エクスポートには `frame`、`energy_hartree`、および delta-E カラム（`delta_kcal`/`delta_hartree`）または絶対カラム（基準適用なし時の `energy_kcal`/`energy_hartree`）が含まれます。
- PNG は高解像度のため `scale=2` で Plotly の PNG エクスポートを使用します。
- コメントモードの JSON は `energy_source: trajectory_comment` を記録し、`mlip_backend`、`mlip_model`、`mlip_precision`、`charge`、`multiplicity` は null です。再計算時は `energy_source: mlip_recomputed` と解決済みの値を記録します。
- `--out-json` 指定時は `result.json` と `summary.json` に同一のペイロードが書き込まれます。
- JSON では、順序を保持する `output_files` を正規の出力一覧として使用してください。後方互換用の `files` はベース名をキーにするため、別ディレクトリにある同名出力を 2 件とも表現できません。

## CLI オプション

全フラグの一覧は生成された [コマンドリファレンス](../reference/commands/index.md) にあります。以下の表は説明が必要なオプションを扱います。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 2 行目にエネルギーが格納された XYZ 軌跡。 | 必須 |
| `-o, --out PATH` | 繰り返し可能な出力ファイル名。`.png`、`.jpg`/`.jpeg`、`.html`、`.svg`、`.pdf`、`.csv` をサポート。 | `energy.png` |
| _追加引数_ | オプション後に列挙された位置ファイル名。`-o` リストとマージ。 | _None_ |
| `--unit {kcal,hartree}` | プロット/エクスポートされる値のターゲット単位。 | `kcal` |
| `-r, --reference TEXT` | 基準仕様（`init`、`None`、または 0 始まり整数）。 | `init` |
| `-q, --charge INT` | MLIP 再計算に使う総電荷。指定時に再計算を実行。 | _None_ |
| `-m, --multiplicity INT` | MLIP 再計算に使うスピン多重度 (2S+1)。指定時に再計算を実行。 | _None_ |
| `-b, --backend {uma,orb,mace,aimnet2}` | 再計算で使用する MLIP バックエンド。 | `uma` |
| `--backend-model TEXT` | 選択したバックエンドのモデル。 | バックエンドのデフォルト |
| `--precision {fp32,fp64}` | 大文字・小文字を区別しない、バックエンド共通の再計算精度。 | バックエンドのデフォルト |
| `--out-json/--no-out-json` | 最初の出力と同じディレクトリに `result.json` を書き出す。 | `False` |
| `--reverse-x/--no-reverse-x` | X 軸を反転し、最後のフレームを左側に表示します（`init` は最後のフレームになります）。 | `False` |

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [path-search](path-search.md) -- 再帰的 MEP 探索（trj2fig に適した XYZ 軌跡を生成）
- [irc](irc.md) -- TS からの IRC（エネルギープロファイリング用軌跡を生成）
- [all](all.md) -- 一気通貫ワークフロー
