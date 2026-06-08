# `fix-altloc`

PDB ファイルから代替位置（altLoc）指示子を削除し、各原子について占有率に基づいて最良のコンフォーマーを選択して重複を除去します。altLoc 文字を含む構造を下流の ML/MM 準備の前にクリーンアップしたいとき（通常は [add-elem-info](add-elem-info.md) で元素カラムを修復した後）に使用します。altLoc カラム（列 17、1-based）はスペース 1 文字に置換され（1 文字の置換のみで、フォーマットのシフトや再構成は行いません）、同じ原子が複数の altLoc 状態で出現する場合は占有率が最も高いコピーを保持します（同点の場合、または占有率が欠損している場合はファイル内で最初に出現したもの）。`ATOM` / `HETATM` レコードは altLoc の選択とブランク化の対象となり、`ANISOU` レコードは対応する ATOM/HETATM 行（同じシリアル番号）が保持される場合のみ保持されます。

## 実行例

コマンド形式:

```bash
mlmm fix-altloc -i INPUT [-o OUTPUT] [options]
```

単一ファイルの altLoc を解決する（`<input>_clean.pdb` を出力）:

```bash
mlmm fix-altloc -i 1abc.pdb
```

出力名を明示して単一ファイルの altLoc を解決する:

```bash
mlmm fix-altloc -i 1abc.pdb -o 1abc_fixed.pdb
```

ディレクトリを再帰的に処理して新しい出力ディレクトリへ書き出す:

```bash
mlmm fix-altloc -i ./structures -o ./cleaned --recursive
```

ディレクトリを再帰的に処理してファイルをその場で上書きする:

```bash
mlmm fix-altloc -i ./structures --inplace --recursive
```

altLoc が検出されなくても強制的に処理するには `--force` を使用します。

```bash
mlmm fix-altloc -i 1abc.pdb -o 1abc_fixed.pdb --force
```

## 処理の流れ

1. 入力ファイルに非空白の altLoc 文字（列 17）が含まれているかチェック。
 - altLoc が見つからず `--force` が設定されていない場合、ファイルをスキップ。
2. 各 ATOM/HETATM レコードについて、altLoc フィールドを無視した同一性キーを構築:
 - レコード名、原子名、残基名、チェーン ID、残基番号、挿入コード、segID
3. 同じ同一性キーを持つ原子の中から、以下の基準で選択:
 - 最高の占有率（列 55-60）
 - 同点の場合、ファイル内で最初に出現したもの
4. 出力を書き込み:
 - 選択された原子のみを保持
 - altLoc カラム（17）を空白（スペース 1 文字）に置換
 - ANISOU レコードは保持された原子に一致するもののみフィルタリング

### altLoc 状態間で原子数が異なる場合の処理

異なる altLoc 状態で異なる原子が含まれている場合（例：altLoc A には N, CA, CB, CG、
altLoc B には N, CA, CB, CD がある場合）、`fix-altloc` は以下のように正しく処理します：

- **重複原子**（複数の altLoc で同じ残基＋原子名、例：N, CA, CB）：
  占有率に基づいて最良のものを選択（最高値を優先、同点の場合はファイル内で最初のもの）
- **ユニーク原子**（1 つの altLoc にのみ存在、例：A の CG、B の CD）：
  **すべてのユニーク原子が出力に保持されます**

これにより、出力構造にはすべての altLoc 状態の原子が含まれ、
真の重複のみが単一のコンフォーマーに解決されます。

**例:**
```
入力:
 ATOM 1 N AALA A 1... 0.50 # altLoc A
 ATOM 2 CA AALA A 1... 0.50 # altLoc A
 ATOM 3 CG AALA A 1... 0.50 # altLoc A のみ
 ATOM 4 N BALA A 1... 0.40 # altLoc B
 ATOM 5 CA BALA A 1... 0.40 # altLoc B
 ATOM 6 CD BALA A 1... 0.40 # altLoc B のみ

出力:
 ATOM 1 N ALA A 1... 0.50 # A から（高 occupancy）
 ATOM 2 CA ALA A 1... 0.50 # A から（高 occupancy）
 ATOM 3 CG ALA A 1... 0.50 # 保持（A のみ）
 ATOM 6 CD ALA A 1... 0.40 # 保持（B のみ）
```

## 出力

- 代替位置が削除された PDB ファイル:
 - 入力がファイル: デフォルトは `<input>_clean.pdb`（`-o/--out` が省略された場合）
 - 入力がディレクトリ: デフォルトは `<input>_clean/`（サブパスを保持）
 - `-o/--out` 指定時: `OUTPUT.pdb`
 - `--inplace` 設定時: 元のファイルを上書き（バックアップは `<input>.pdb.bak` として保存）

元のファイルは変更されません（`--inplace` が設定されていない限り）。

## Python API

プログラム的な使用のために、モジュールは以下をエクスポートします:
```python
from pathlib import Path
from mlmm.io.pdb_fix import has_altloc, clean_pdb_file

# ファイルに altLoc があるかチェック
if has_altloc(Path("input.pdb")):
    # altLoc を解決した PDB を書き出す (出力は常に上書き)
    clean_pdb_file(Path("input.pdb"), Path("output.pdb"))
```

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力 PDB ファイルまたはディレクトリ | 必須 |
| `-o, --out PATH` | 出力ファイル（入力がファイルの場合）またはディレクトリ（入力がディレクトリの場合） | 入力がファイル: `<input>_clean.pdb`、入力がディレクトリ: `<input>_clean/` |
| `--recursive/--no-recursive` | 入力がディレクトリの場合、`*.pdb` ファイルを再帰的に処理 | `False` |
| `--inplace/--no-inplace` | 入力ファイルをその場で上書き（`.bak` バックアップを作成） | `False` |
| `--overwrite/--no-overwrite` | 既存の出力ファイルの上書きを許可 | `False` |
| `--force/--no-force` | altLoc が検出されなくてもファイルを処理 | `False` |

全フラグの一覧は生成された[コマンドリファレンス](../reference/commands/index.md)を参照してください。

## 注記

- altLoc 文字を含まないファイルは `--force` を設定しない限りスキップされます。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [add-elem-info](add-elem-info.md) -- altLoc 修正前に PDB 元素カラムを修復
- [extract](extract.md) -- altLoc 解決後に活性部位ポケットを抽出
- [all](all.md) -- ML/MM 一気通貫ワークフロー（入力に altLoc がある場合は事前に `fix-altloc` を実行）
