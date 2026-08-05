# `fix-altloc`

PDB ファイルの代替位置（altLoc）を、原子単位ではなく残基単位で
解決します。各残基について、ラベル付き原子の平均 occupancy が最大の
非空白ラベルを一つ選び、同値の場合は最初に現れるラベルを選びます。
空白（共通）原子は残し、非選択ラベルの原子は削除します。これにより、
実在しない A/B 混合残基の生成を防ぎます。

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
2. ラベル付き ATOM/HETATM レコードを site（chain ID、残基番号、
   insertion code、segID）ごとにまとめる。
3. 各残基で、解析可能なoccupancy（列55–60）の平均が最大のラベルを選ぶ。
   occupancyを1件も解析できないラベルは、解析可能な平均を持つラベルより下位になる。
   scoreが同じ場合（全ラベルでoccupancy欠損の場合を含む）は最初の出現順で決める。
4. 空白（共通）原子と選択ラベルの原子を残し、残る不正な重複は occupancy
   と出現順で解決する。
5. 出力を書き込み:
 - 空白（共通）原子と選択した残基 conformer のみを保持
 - altLoc 列（17）を空白（スペース 1 文字）に置換
 - ANISOU レコードは保持された原子に一致するもののみフィルタリング

### altLoc 状態間で原子数が異なる場合の処理

異なる altLoc 状態で異なる原子が含まれている場合（例：altLoc A には N, CA, CB, CG、
altLoc B には N, CA, CB, CD がある場合）、`fix-altloc` は以下のように処理します：

選択した残基ラベルに属する原子だけを残します。非選択ラベルにしかない
原子は削除します。

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
 ATOM 1 N ALA A 1... 0.50 # A から（占有率が高い）
 ATOM 2 CA ALA A 1... 0.50 # A から（占有率が高い）
 ATOM 3 CG ALA A 1... 0.50 # 保持（A のみ）
```

## 出力

- 代替位置が削除された PDB ファイル:
 - 入力がファイル: デフォルトは `<input>_clean.pdb`（`-o/--out` が省略された場合）
 - 入力がディレクトリ: デフォルトは `<input>_clean/`（サブパスを保持）
 - `-o/--out` 指定時: `OUTPUT.pdb`
 - `--inplace` 設定時: 元のファイルを上書き（バックアップは `<input>.pdb.bak` として保存）

元のファイルは変更されません（`--inplace` が設定されていない限り）。

## Python API

プログラムから利用する場合、モジュールは以下をエクスポートします:
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

- [add-elem-info](add-elem-info.md) -- altLoc 修正前に PDB 元素列を修復
- [extract](extract.md) -- altLoc 解決後に活性部位ポケットを抽出
- [all](all.md) -- ML/MM 一気通貫ワークフロー（入力に altLoc がある場合は事前に `fix-altloc` を実行）
