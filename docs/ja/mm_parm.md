# `mm-parm`

## 概要

> **要約:** PDB から Amber prmtop/rst7 トポロジーファイルを構築します。GAFF2 リガンドの自動パラメータ化、ジスルフィド結合検出、PDBFixer による任意の水素付加に対応。

### 概要
- **用途:** ML/MM 計算に必要な AMBER トポロジーおよび座標ファイルが必要な場合に使用します。
- **手法:** tleap + antechamber (GAFF2, AM1-BCC) + parmchk2、任意で PDBFixer による水素付加。
- **出力:** `<prefix>.parm7`、`<prefix>.rst7`、および任意で `<prefix>.pdb`（LEaP 出力）。
- **デフォルト:** ff19SB（タンパク質）+ OPC3（水）+ GAFF2（有機分子）+ Lipid21 + GLYCAM_06j-1 + OL3/OL21。
- **次のステップ:** [define-layer](define_layer.md) で ML/MM レイヤーを割り当て。

`mlmm mm-parm` は PDB から Amber トポロジー/座標ファイルを生成します。入力 PDB はデフォルトでは構造修正なしにそのまま使用されます。不明な残基は antechamber（GAFF2、AM1-BCC 電荷）と parmchk2 で自動的にパラメータ化されます。ジスルフィド結合は SG-SG（または S-S）距離が 2.3 Å 以内であることから幾何学的に推定されます。非標準アミノ酸（選択された力場で認識されない N/CA/C を含む残基）は自動処理**されません**。手動でパラメータ化するようメッセージを表示してビルドを中断します。

`--add-h` の場合、指定された `--ph` で PDBFixer を使用して tleap 処理前に水素が付加されます。他の構造修正は行われません。`--ff-set ff14SB` を使用すると、力場は ff14SB（タンパク質）+ TIP3P（水）（+ phosaa14SB）に切り替わります。

## 最小例

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 --ligand-charge "GPP=-3,MMT=-1" --ligand-mult "GPP=1,MMT=1"
```

## 出力の見方

- `<prefix>.parm7` -- Amber prmtop トポロジー
- `<prefix>.rst7` -- Amber ASCII inpcrd 座標
- `<prefix>.pdb` -- LEaP savepdb 出力（下記の命名規則を参照）

## よくある例

1. 水素付加付きの基本的なトポロジー構築。

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 --ligand-charge "GPP=-3,MMT=-1" --ligand-mult "GPP=1,MMT=1" \
 --add-ter --ff-set ff19SB --add-h --ph 7.0
```

2. 水素付加なしのビルド。

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 --ligand-charge "GPP=-3" --no-add-h
```

## ワークフロー

1. **入力準備** -- 入力 PDB はそのまま読み込まれます（構造修正なし）。`--add-h` が設定されている場合、PDBFixer で指定 `--ph` にて水素が付加されます。
2. **TER 挿入** -- `--add-ter`（デフォルト）の場合、リガンド/水/イオン残基の連続ブロックの前後に TER レコードが挿入されます。
3. **不明残基のパラメータ化** -- 力場で認識されない残基は antechamber（GAFF2、AM1-BCC）と parmchk2 で自動パラメータ化されます。形式電荷とスピン多重度は `--ligand-charge` と `--ligand-mult` で制御されます。
4. **ジスルフィド検出** -- CYS/CYM/CYX ペアで SG-SG（または S-S）距離が 2.3 Å 以下のものが自動的に結合されます。
5. **トポロジー構築** -- tleap が選択された力場セットを使用して parm7/rst7/pdb ファイルを生成します。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力 PDB（`--add-h` でない限りそのまま使用）。 | 必須 |
| `--out-prefix TEXT` | parm7/rst7/pdb ファイルの出力接頭辞。 | 入力 PDB のファイル名幹 |
| `--ligand-charge TEXT` | 残基名と形式電荷のマッピング（例: `"GPP=-3,MMT=-1"`）。 | _None_ |
| `--ligand-mult TEXT` | 残基名とスピン多重度（`-m`）のマッピング（例: `"HEM=1,NO=2"`）。 | `1` |
| `--allow-nonstandard-aa` | 非標準アミノ酸様残基（N/CA/C を含む）の antechamber パラメータ化を許可。 | `False` |
| `--keep-temp` | 作業ディレクトリの中間ファイル/ログを保持（デバッグ用）。 | `False` |
| `--add-ter/--no-add-ter` | リガンド/水/イオンブロックの前後に TER を挿入。 | `True` |
| `--add-h/--no-add-h` | PDBFixer で `--ph` に基づいて水素を付加。 | `False` |
| `--ph FLOAT` | PDBFixer の水素付加用 pH（`--add-h` の場合のみ使用）。 | `7.0` |
| `--ff-set {ff19SB\|ff14SB}` | 力場セット: ff19SB（デフォルト）または ff14SB。 | `ff19SB` |

## 注意事項
- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- **AmberTools**（`tleap`、`antechamber`、`parmchk2`）が PATH 上に利用可能であること。
- **PDBFixer** + **OpenMM** は `--add-h` の場合**のみ**必要。
- ビルドが失敗した場合でも、`--add-h` で水素付加が成功していれば、水素付加済み PDB はディスクに書き出されます。

### 力場

- **ff19SB**（デフォルト）: ff19SB（タンパク質）+ OPC3（水）+ GAFF2（一般有機分子）+ Lipid21 + GLYCAM_06j-1 + OL3/OL21（+ phosaa19SB）。
- **ff14SB**: ff14SB（タンパク質）+ TIP3P（水）（+ phosaa14SB）。

### 出力の命名規則

- `<prefix>.parm7` -- prmtop トポロジー
- `<prefix>.rst7` -- ASCII inpcrd（LEaP が生成した `complex.inpcrd` のコピー）
- `<prefix>.pdb` -- LEaP `savepdb` 出力:
 - `--out-prefix` 指定時: `<out_prefix>.pdb`
 - `--out-prefix` 省略かつ `--add-h`: `<input_stem>_parm.pdb`
 - `--out-prefix` 省略かつ `--no-add-h`: PDB は書き出されません

### 失敗時の最短復旧レシピ

1. AmberTools コマンドが見つからない。

```bash
which tleap antechamber parmchk2
mlmm mm-parm -i input.pdb --keep-temp
```

2. 非標準残基でビルドが停止する。

```bash
mlmm mm-parm -i input.pdb --keep-temp
# keep-temp で残した tleap ログを確認し、残基パラメータを追加して再実行
```

3. `--add-h` が環境依存で失敗する。

```bash
python -c "import pdbfixer, openmm; print('ok')"
mlmm mm-parm -i input.pdb --no-add-h
```

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [all](all.md) -- エンドツーエンドワークフロー（内部で mm-parm を呼び出し）
- [extract](extract.md) -- パラメータ化前に活性部位ポケットを抽出
- [define-layer](define_layer.md) -- トポロジー構築後に ML/MM レイヤーを定義
