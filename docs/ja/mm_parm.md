# `mm-parm`

## 概要

> **要約:** PDB から Amber prmtop/rst7 トポロジーファイルを構築します。GAFF2 リガンドの自動パラメータ化、ジスルフィド結合検出、PDBFixer による任意の水素付加に対応。

### 概要
- **用途:** ML/MM 計算に必要な AMBER トポロジーおよび座標ファイルが必要な場合に使用します。
- **入力:** 1 つの PDB ファイル（`--add-h True` でない限りそのまま使用）。
- **デフォルト力場:** ff19SB（タンパク質）+ OPC3（水）+ GAFF2（有機分子）+ Lipid21 + GLYCAM_06j-1 + OL3/OL21。
- **出力:** `<prefix>.parm7`、`<prefix>.rst7`、および任意で `<prefix>.pdb`（LEaP 出力）。
- **要件:** AmberTools（tleap、antechamber、parmchk2）が PATH 上にあること。PDBFixer + OpenMM は `--add-h True` の場合のみ必要。

`mlmm mm-parm` は PDB から Amber トポロジー/座標ファイルを生成します。入力 PDB はデフォルトでは構造修正なしにそのまま使用されます。不明な残基は antechamber（GAFF2、AM1-BCC 電荷）と parmchk2 で自動的にパラメータ化されます。ジスルフィド結合は SG-SG（または S-S）距離が 2.3 A 以内であることから幾何学的に推定されます。非標準アミノ酸（選択された力場で認識されない N/CA/C を含む残基）は自動処理**されません** -- 手動でパラメータ化するようメッセージを表示してビルドを中断します。

`--add-h True` の場合、指定された `--pH` で PDBFixer を使用して tleap 処理前に水素が付加されます。他の構造修正は行われません。`--ff-set ff14SB` を使用すると、力場は ff14SB（タンパク質）+ TIP3P（水）（+ phosaa14SB）に切り替わります。

## 使用法

```bash
mlmm mm-parm -i INPUT.pdb [--out-prefix PREFIX] \
             [--ligand-charge "RES=-3,ABC:+1"] \
             [--ligand-mult "RES=2,ABC:1"] \
             [--keep-temp] \
             [--add-ter True|False] \
             [--add-h True|False] [--pH 7.0] \
             [--ff-set ff19SB|ff14SB]
```

### 例

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
    --ligand-charge "GPP=-3,MMT=-1" --ligand-mult "GPP=1,MMT=1" \
    --add-ter True --ff-set ff19SB --add-h True --pH 7.0
```

## 説明

### ジスルフィド結合検出

CYS/CYM/CYX のジスルフィド結合は PDB 座標から純粋に幾何学的に推定されます: SG-SG（または S-S）距離が 2.3 A 以下の残基ペアが結合として検出されます。

### 力場

- **ff19SB**（デフォルト）: ff19SB（タンパク質）+ OPC3（水）+ GAFF2（一般有機分子）+ Lipid21 + GLYCAM_06j-1 + OL3/OL21（+ phosaa19SB）。
- **ff14SB**: ff14SB（タンパク質）+ TIP3P（水）（+ phosaa14SB）。

### 不明な残基

LEaP が不明な残基を報告した場合、antechamber（GAFF2、AM1-BCC 電荷）と parmchk2 で自動的にパラメータ化されます。形式電荷とスピン多重度は `--ligand-charge` と `--ligand-mult` で制御できます。

### 自動 TER 挿入

`--add-ter True`（デフォルト）の場合、`--ligand-charge` に名前が含まれる残基、水残基、またはイオンの連続ブロックの**前後**に TER レコードが挿入されます。連続する同種の残基間には TER レコードは挿入されません。

### 失敗時の動作

ビルドが失敗した場合でも、`--add-h True` で水素付加が成功していれば、水素付加済み PDB は LEaP PDB と同じ命名規則（`--out-prefix` 省略時は `<input_stem>_parm.pdb`、指定時は `<out_prefix>.pdb`）でディスクに書き出されます。

### 出力の命名規則

- `<prefix>.parm7` -- prmtop トポロジー
- `<prefix>.rst7` -- ASCII inpcrd（LEaP が生成した `complex.inpcrd` のコピー）
- `<prefix>.pdb` -- LEaP `savepdb` 出力:
  - `--out-prefix` 指定時: `<out_prefix>.pdb`
  - `--out-prefix` 省略かつ `--add-h True`: `<input_stem>_parm.pdb`
  - `--out-prefix` 省略かつ `--add-h False`: PDB は書き出されません

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --pdb PATH` | 入力 PDB（`--add-h True` でない限りそのまま使用）。 | 必須 |
| `--out-prefix TEXT` | parm7/rst7/pdb ファイルの出力接頭辞。 | 入力 PDB のファイル名幹 |
| `--ligand-charge TEXT` | 残基名と形式電荷のマッピング（例: `"GPP=-3,MMT=-1"`）。 | _None_ |
| `--ligand-mult TEXT` | 残基名とスピン多重度（`-m`）のマッピング（例: `"HEM=1,NO:2"`）。 | `1` |
| `--keep-temp` | 作業ディレクトリの中間ファイル/ログを保持（デバッグ用）。 | `False` |
| `--add-ter {True\|False}` | リガンド/水/イオンブロックの前後に TER を挿入。 | `True` |
| `--add-h {True\|False}` | PDBFixer で `--pH` に基づいて水素を付加。 | `False` |
| `--pH FLOAT` | PDBFixer の水素付加用 pH（`--add-h True` の場合のみ使用）。 | `7.0` |
| `--ff-set {ff19SB\|ff14SB}` | 力場セット: ff19SB（デフォルト）または ff14SB。 | `ff19SB` |

## 出力

```
<prefix>.parm7         # Amber prmtop トポロジー
<prefix>.rst7          # Amber ASCII inpcrd 座標
<prefix>.pdb           # LEaP savepdb 出力（上記の命名規則を参照）
```

## 要件

- **AmberTools**（`tleap`、`antechamber`、`parmchk2`）が PATH 上に利用可能であること。
- **PDBFixer** + **OpenMM** は `--add-h True` の場合**のみ**必要。

---

## 関連項目

- [all](all.md) -- エンドツーエンドワークフロー（内部で mm-parm を呼び出し）
- [extract](extract.md) -- パラメータ化前に活性部位ポケットを抽出
- [define-layer](define_layer.md) -- トポロジー構築後に ML/MM レイヤーを定義
