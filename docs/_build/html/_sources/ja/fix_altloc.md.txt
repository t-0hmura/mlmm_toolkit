# `fix-altloc`

## 概要
`fix-altloc` は、PDB ファイルから代替位置（altLoc）指示子を削除します。
各原子について占有率（occupancy）に基づいて最良のコンフォーマーを選択し、
重複を除去します。

### 機能
1. PDBのaltLocカラム（列17、1-based）を空白（スペース1文字）に置換。
 - 1文字の置換のみで、フォーマットのシフトや再構成は行いません。
2. 代替位置（A/B/... やカスタムラベル H/L など）により同じ原子が複数回
 出現する場合、「最良」のものを保持:
 - 占有率（occupancy）が最も高いものを優先
 - 同点の場合（または占有率が欠損している場合）、ファイル内で最初に出現したものを保持

### 処理対象レコード
- `ATOM` / `HETATM`: altLocの選択とブランク化
- `ANISOU`: 対応するATOM/HETATMレコード（同じシリアル番号）が保持される場合のみ保持

### 自動スキップ動作
デフォルトでは、ファイルに**altLoc文字が含まれていない**場合（列17がすべて空白）、
ファイルは**スキップ**され、出力は書き込まれません。altLocの有無に関わらず
処理を行うには `--force` を使用してください。

## 使用法
```bash
mlmm fix-altloc -i INPUT.pdb [-o OUTPUT.pdb] [OPTIONS]
```

## 例
```bash
# 単一ファイルを処理（出力: INPUT_clean.pdb）
mlmm fix-altloc -i 1abc.pdb

# 出力ファイルを指定
mlmm fix-altloc -i 1abc.pdb -o 1abc_fixed.pdb

# ディレクトリを再帰的に処理
mlmm fix-altloc -i ./structures -o ./cleaned --recursive

# 入力ファイルをその場で上書き（.bakバックアップを作成）
mlmm fix-altloc -i ./structures --inplace --recursive

# altLocが検出されなくても強制的に処理
mlmm fix-altloc -i 1abc.pdb -o 1abc_fixed.pdb --force
```

## ワークフロー
1. 入力ファイルに非空白のaltLoc文字（列17）が含まれているかチェック。
 - altLocが見つからず `--force` が設定されていない場合、ファイルをスキップ。
2. 各ATOM/HETATMレコードについて、altLocフィールドを無視した同一性キーを構築:
 - レコード名、原子名、残基名、チェーンID、残基番号、挿入コード、segID
3. 同じ同一性キーを持つ原子の中から、以下の基準で選択:
 - 最高の占有率（列55–60）
 - 同点の場合、ファイル内で最初に出現したもの
4. 出力を書き込み:
 - 選択された原子のみを保持
 - altLocカラム（17）を空白（スペース1文字）に置換
 - ANISOUレコードは保持された原子に一致するもののみフィルタリング

## CLI オプション
| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力 PDB ファイルまたはディレクトリ | 必須 |
| `-o, --out PATH` | 出力ファイル（入力がファイルの場合）またはディレクトリ（入力がディレクトリの場合） | 入力がファイル: `<input>_clean.pdb`、入力がディレクトリ: `<input>_clean/` |
| `--recursive/--no-recursive` | 入力がディレクトリの場合、`*.pdb` ファイルを再帰的に処理 | `False` |
| `--inplace/--no-inplace` | 入力ファイルをその場で上書き（`.bak` バックアップを作成） | `False` |
| `--overwrite/--no-overwrite` | 既存の出力ファイルの上書きを許可 | `False` |
| `--force/--no-force` | altLocが検出されなくてもファイルを処理 | `False` |

## 出力
- 代替位置が削除されたPDB ファイル:
 - 入力がファイル: デフォルトは `<input>_clean.pdb`（`-o/--out` が省略された場合）
 - 入力がディレクトリ: デフォルトは `<input>_clean/`（サブパスを保持）
 - `-o/--out` 指定時: `OUTPUT.pdb`
 - `--inplace` 設定時: 元のファイルを上書き（バックアップは `<input>.pdb.bak` として保存）

## `all` ワークフローとの統合
`mlmm all` ワークフローを実行する際、`fix-altloc` は **`add-elem-info`の後**
（元素フィールドが欠損していた場合）、**ポケット抽出の前**に自動的に実行されます。
これにより:
1. まず元素記号が補完される
2. 代替位置が単一のコンフォーマーに解決される
3. 抽出ステップがクリーンで曖昧さのない座標を受け取る

altLoc文字が検出された場合のみファイルが処理され、そうでない場合は元のファイルが
そのまま渡されます。

## altLoc状態間で原子数が異なる場合の処理

異なるaltLoc状態で異なる原子が含まれている場合（例：altLoc Aには N, CA, CB, CG、
altLoc Bには N, CA, CB, CD がある場合）、`fix-altloc` は以下のように正しく処理します：

- **重複原子**（複数のaltLocで同じ残基＋原子名、例：N, CA, CB）：
 占有率に基づいて最良のものを選択（最高値を優先、同点の場合はファイル内で最初のもの）
- **ユニーク原子**（1つのaltLocにのみ存在、例：AのCG、BのCD）：
 **すべてのユニーク原子が出力に保持されます**

これにより、出力構造にはすべてのaltLoc状態の原子が含まれ、
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
 ATOM 1 N ALA A 1... 0.50 # Aから（高occupancy）
 ATOM 2 CA ALA A 1... 0.50 # Aから（高occupancy）
 ATOM 3 CG ALA A 1... 0.50 # 保持（Aのみ）
 ATOM 6 CD ALA A 1... 0.40 # 保持（Bのみ）
```

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- 原子シリアル番号は**再番号付けされません**（重複削除後にギャップが残る場合があります）。
- `CONECT` およびその他の接続/注釈レコードは**更新されません**。
- 変更されるのは列17（altLoc）のみ。座標、占有率、B因子、電荷、挿入コード、
 レコード順序は保持されます（重複削除を除く）。
- MODEL/ENDMDLブロックは独立して処理されます。

## API使用法
プログラム的な使用のために、モジュールは以下をエクスポートします:
```python
from mlmm.fix_altloc import has_altloc, fix_altloc_file

# ファイルにaltLocがあるかチェック
if has_altloc(Path("input.pdb")):
 # altLocを修正
 was_processed = fix_altloc_file("input.pdb", "output.pdb", overwrite=True)
```
