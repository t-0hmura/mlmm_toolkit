# `fix-altloc`（PDB の代替位置の解決）

`fix-altloc` サブコマンドは、PDB ファイルから**代替位置（altLoc）を取り除きます**。残基ごとに平均占有率が最も高い altLoc ラベルを 1 つ選ぶので、各残基は実際に登録された 1 つのコンフォマーになります。`extract`・`define-layer` と ML/MM の計算コマンドも、PDB を読み込むときに同じ規則を自動で適用します。整理した PDB ファイルそのものが必要なときに `fix-altloc` を使ってください。

---

## 主な用途

* **保存用の整理した PDB**: 残基ごとに 1 つのコンフォマーにしたファイルを、ほかのプログラムや記録に使う
* **`mm-parm` の入力**: `mm-parm` は altLoc を解決しないので、先に整理する
* **多数のファイルの一括処理**: ディレクトリ内のすべての `.pdb`（必要ならサブディレクトリも）
* **選択の確認**: 計算の前に、どのコンフォマーが残るかを確かめる

---

## 基本的な実行例

### 1. 1 つのファイル

1 つのファイルを整理し、`1abc_clean.pdb` に書き出します。

```bash
mlmm fix-altloc -i 1abc.pdb
```

端末に `[fix-altloc] Fixed altLoc → 1abc_clean.pdb` が出れば成功です。altLoc の無いファイルでは `[fix-altloc] Skipped 1abc.pdb (no altLoc detected).` が出ます。

### 2. 出力ファイルを指定する

```bash
mlmm fix-altloc -i 1abc.pdb -o 1abc_fixed.pdb
```

### 3. ディレクトリを再帰的に処理する

`./structures` 以下のすべての `.pdb` を整理し、同じサブディレクトリ構成で `./cleaned` に書き出します。

```bash
mlmm fix-altloc -i ./structures -o ./cleaned --recursive
```

### 4. バックアップを残して上書きする

入力ファイルを上書きし、元のファイルを `<name>.pdb.bak` として残します。

```bash
mlmm fix-altloc -i ./structures --inplace --recursive
```

---

## 処理の仕組みと計算仕様

1. **altLoc の検出**:
各ファイルで、空白でない altLoc の文字（17 列目）を探します。
2. **残基ごとのまとめ**:
ラベルの付いた ATOM・HETATM レコードを、chain ID・残基番号・挿入コード・segID で残基ごとにまとめます。残基名はキーに含めません。
3. **残基ごとに 1 つのラベルを選ぶ**:
原子の平均占有率（55–60 列）が最も高いラベルを選びます。同点ならファイル内で先に出たラベルを選びます。
4. **書き出し**:
空白（共通）の原子と、選んだラベルの原子を残し、17 列目を空白にします。選んだラベルに同じ原子があるときは、空白の原子のほうを除きます。

ANISOU レコードは、残った原子（同じシリアル番号）の分だけを残し、そのほかのレコードはそのまま書き出します。

### altLoc の間で原子数が違う場合

altLoc の状態ごとに原子が違う場合も、選んだラベルの原子だけを残し、選ばなかったラベルにしか無い原子は削除します。1 つの残基の中で A と B の原子が混ざることはありません。

```text
入力:
 ATOM 1 N ALYS A 1... 0.50 # altLoc A
 ATOM 2 CA ALYS A 1... 0.50 # altLoc A
 ATOM 3 CB ALYS A 1... 0.50 # altLoc A
 ATOM 4 CG ALYS A 1... 0.50 # altLoc A
 ATOM 5 N BLYS A 1... 0.40 # altLoc B
 ATOM 6 CA BLYS A 1... 0.40 # altLoc B
 ATOM 7 CB BLYS A 1... 0.40 # altLoc B
 ATOM 8 CG BLYS A 1... 0.40 # altLoc B
 ATOM 9 CD BLYS A 1... 0.40 # altLoc B のみ

出力:
 ATOM 1 N LYS A 1... 0.50 # A から（占有率が高い）
 ATOM 2 CA LYS A 1... 0.50 # A から
 ATOM 3 CB LYS A 1... 0.50 # A から
 ATOM 4 CG LYS A 1... 0.50 # A から
 （altLoc B にしか無い CD は削除）
```

---

## 主な出力ファイル

* **ファイル入力**: デフォルトは `<input>_clean.pdb`。`-o` を指定するとそのパスです。`-o` が `.pdb` で終わらないときはディレクトリとして扱い、その中に入力と同じ名前で書き出します。
* **ディレクトリ入力**: デフォルトは `<input>_clean/`。`-o` を指定するとそのディレクトリで、入力と同じ相対パスに書き出します。端末に `[fix-altloc] Processed N file(s) → …` と、altLoc の無いファイルについて `Skipped N file(s)` が出ます。
* **`--inplace`**: 入力ファイルを上書きし、元のファイルを `<name>.pdb.bak` として保存します。

---

## 主な CLI オプション

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | 入力 PDB ファイルまたはディレクトリ |
| `-o, --output` | パス | `None` | 出力ファイル（ファイル入力）またはディレクトリ（ディレクトリ入力）。省略時は `<input>_clean.pdb` または `<input>_clean/` |
| `--recursive/--no-recursive` | フラグ | `False` | ディレクトリ入力で、サブディレクトリの `.pdb` も処理 |
| `--inplace/--no-inplace` | フラグ | `False` | 入力ファイルを上書き（`.bak` のバックアップを作成） |
| `--overwrite/--no-overwrite` | フラグ | `False` | 既存の出力ファイルの上書きを許可。無いときに出力がすでにあると `Output exists: <path> (use --overwrite to overwrite)` で止まる |
| `--force/--no-force` | フラグ | `False` | altLoc が見つからないファイルも処理 |

全オプションは [自動生成のオプションの一覧（英語のみ）](../reference/commands/fix_altloc.md) を参照してください。

---

## 使用上の注意点

* **altLoc の無いファイル**: 17 列目がすべて空白のファイルはスキップし、何も書き出しません。`--force` を付けると処理します。
* **`--inplace` と `-o`**: `--inplace` のときは `-o` を無視します。`.bak` がすでにあれば置き換えないので、最初に上書きする前のファイルが残ります。
* **シリアル番号**は振り直さないので、原子を削除した箇所に欠番が残ることがあります。`CONECT` などの結合・注釈のレコードも更新しません。
* **残したレコード**は 17 列目のほかは書き換えないので、座標・占有率・B-factor・電荷・挿入コード・並び順はそのままです。
* **MODEL/ENDMDL ブロック**はブロックごとに別々に処理します。
* **占有率の規則は経験則です**: 化学的な接触や登録されたアンサンブルの解釈で活性部位のコンフォマーを選ぶ必要があるときは、構造エディタで自分で選び、目で確かめてください。
* **終了コード**: {ref}`終了コード <ja-exit-codes>`を参照してください。

---

## 関連ドキュメント

* [extract](extract.md) — PDB の読み込み時に同じ altLoc の規則を適用する活性部位モデルの抽出
* [mm-parm](mm-parm.md) — 整理した PDB から Amber のトポロジーを作る
* [add-elem-info](add-elem-info.md) — PDB の元素列（77–78 列）を埋める
* [all](all.md) — 全工程のワークフロー
* [トラブルシューティング](troubleshooting.md) — 実行に失敗したときの対処
