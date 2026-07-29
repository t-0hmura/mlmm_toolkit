# `add-elem-info`

PDB の元素記号（列 77-78）を追加または修復します。固定列の原子名と残基情報から元素を推定し、ATOM/HETATM の元素列だけを置換します。それ以外の列とレコードは入力のまま保持されます。既存の元素フィールドも再推定する場合は `--overwrite` を指定します。

## 実行例

コマンド形式:

```bash
mlmm add-elem-info -i INPUT [-o OUTPUT] [--inplace] [--overwrite]
```

元素列を追加または修復し、安全なデフォルト出力へ書き出す:

```bash
mlmm add-elem-info -i 1abc.pdb
```

このコマンドは `1abc_add_elem.pdb` を書き出します。入力ファイルを明示的に
置換する場合は、次のように指定します:

```bash
mlmm add-elem-info -i 1abc.pdb --inplace
```

結果を別の出力ファイルに書き出す:

```bash
mlmm add-elem-info -i 1abc.pdb -o 1abc_fixed.pdb
```

既存の元素フィールドを再推定して上書きする:

```bash
mlmm add-elem-info -i 1abc.pdb --overwrite
```

## 処理の流れ

1. PDB の生レコードを読み、`extract.py` と同じ残基定義（`AMINO_ACIDS`、`WATER_RES`、`ION`）で分類します。
2. 各原子について、原子名、残基名、レコードが HETATM かどうかを組み合わせて元素を推定します:
 - **イオン残基:** 残基由来の元素を優先します。多原子イオン（例: NH4、H3O+）は原子ごとに割り当て（H/N/O）。
 - **タンパク質、核酸、水:** H/D は H に、水素・酸素原子は H/O に、P/N/O/S は先頭文字から、Se は専用判定で、炭素ラベル（CA/CB/CG/...）は C に割り当てます。
 - **リガンド/補因子:** 原子名接頭辞（C*/P*、CL を除く）と 2 文字/1 文字の正規化を使用。ハロゲン（Cl/Br/I/F）を認識。
3. ATOM/HETATM の列 77–78 だけを置換し、その他の列とレコードを変更せずに書き出します:
 - `-o/--out` 未指定: `<input>_add_elem.pdb` へ書き出し。
 - `--inplace` を指定し、`-o/--out` を省略: 入力ファイルを置換。
 - `-o/--out` 指定: 指定パスに書き出し。
4. 処理結果のサマリーを出力します: 総原子数、新規割り当て数、既存保持数、上書き数（`--overwrite` 時）、元素ごとのカウント、未解決原子（最大 50 件、モデル/鎖/残基/原子/シリアル番号）。

## 出力

- 元素列（77-78）が正しく設定された PDB ファイル
- コンソールに処理/割り当て済み原子の合計、元素ごとのカウント、未解決原子（最大 50 件）を報告

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力 PDB ファイル。 | 必須 |
| `-o, --out PATH` | 出力 PDB パス。`--inplace` より優先。 | _None_ → `<input>_add_elem.pdb` |
| `--inplace/--no-inplace` | `-o/--out` の省略時に入力ファイルを置換。 | `False` |
| `--overwrite/--no-overwrite` | 既存の元素フィールドがあっても再推定して上書き（デフォルトでは既存値を保持）。 | `False` |

修復対象の ATOM/HETATM レコードの列 77–78 を除き、各入力行はそのまま
保持されます。HEADER、REMARK、CONECT、ANISOU と従来形式の電荷列も
保持されます。

すべてのフラグの一覧は生成された [コマンドリファレンス](../reference/commands/index.md) を参照してください。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [mm-parm](mm-parm.md) -- AMBER トポロジー構築（正しい元素列が必要）
- [extract](extract.md) -- タンパク質-リガンド複合体から活性部位ポケットを抽出
