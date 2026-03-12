# `mm-parm`

## 概要

> **要約:** PDB から Amber prmtop/rst7 トポロジーファイルを構築します。GAFF2 リガンドの自動パラメータ化、ジスルフィド結合検出、PDBFixer による任意の水素付加に対応。

`mlmm mm-parm` は PDB から Amber トポロジー/座標ファイルを生成します。入力 PDB はデフォルトでは構造修正なしにそのまま使用されます。不明な残基は antechamber（GAFF2、AM1-BCC 電荷）と parmchk2 で自動的にパラメータ化されます。`--ligand-charge` に明示した残基は ligand/cofactor 定義として扱われ、GAFF2 パラメータ化が最優先されます。ジスルフィド結合は SG-SG（または S-S）距離が 2.3 Å 以内であることから幾何学的に推定されます。`AMINO_ACIDS` に載っているアミノ酸系残基で、選択された力場に認識されないものは自動処理**されません**。手動でパラメータ化するようメッセージを表示してビルドを中断します。

`--add-h` の場合、指定された `--ph` で PDBFixer を使用して tleap 処理前に水素が付加されます。他の構造修正は行われません。`--ff-set ff14SB` を使用すると、力場は ff14SB（タンパク質）+ TIP3P（水）（+ phosaa14SB）に切り替わります。

## 最小例

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 -l "GPP=-3,MMT=-1" --ligand-mult "GPP=1,MMT=1"
```

## 出力の見方

- `<prefix>.parm7` -- Amber prmtop トポロジー
- `<prefix>.rst7` -- Amber ASCII inpcrd 座標
- `<prefix>.pdb` -- LEaP savepdb 出力（下記の命名規則を参照）

## よくある例

1. 水素付加付きの基本的なトポロジー構築。

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 -l "GPP=-3,MMT=-1" --ligand-mult "GPP=1,MMT=1" \
 --add-ter --ff-set ff19SB --add-h --ph 7.0
```

2. 水素付加なしのビルド。

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 -l "GPP=-3" --no-add-h
```

## ワークフロー

1. **入力準備** -- 入力 PDB はそのまま読み込まれます（構造修正なし）。`--add-h` が設定されている場合、PDBFixer で指定 `--ph` にて水素が付加されます。
2. **TER 挿入** -- `--add-ter`（デフォルト）の場合、リガンド/水/イオン残基の連続ブロックの前後に TER レコードが挿入されます。
3. **不明残基のパラメータ化** -- 力場で認識されない残基は antechamber（GAFF2、AM1-BCC）と parmchk2 で自動パラメータ化されます。`--ligand-charge` に名前を列挙した残基はこの経路が最優先されます。形式電荷とスピン多重度は `--ligand-charge` と `--ligand-mult` で制御されます。
4. **ジスルフィド検出** -- CYS/CYM/CYX ペアで SG-SG（または S-S）距離が 2.3 Å 以下のものが自動的に結合されます。
5. **トポロジー構築** -- tleap が選択された力場セットを使用して parm7/rst7/pdb ファイルを生成します。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力 PDB（`--add-h` でない限りそのまま使用）。 | 必須 |
| `--out-prefix TEXT` | parm7/rst7/pdb ファイルの出力接頭辞。 | 入力 PDB のファイル名幹 |
| `-l, --ligand-charge TEXT` | 残基名と形式電荷のマッピング（例: `"GPP=-3,MMT=-1"`）。 | _None_ |
| `--ligand-mult TEXT` | 残基名とスピン多重度（`-m`）のマッピング（例: `"HEM=1,NO=2"`）。 | `1` |
| `--keep-temp` | 作業ディレクトリの中間ファイル/ログを保持（デバッグ用）。 | `False` |
| `--add-ter/--no-add-ter` | リガンド/水/イオンブロックの前後に TER を挿入。 | `True` |
| `--add-h/--no-add-h` | PDBFixer で `--ph` に基づいて水素を付加。 | `False` |
| `--ph FLOAT` | PDBFixer の水素付加用 pH（`--add-h` の場合のみ使用）。 | `7.0` |
| `--ff-set {ff19SB\|ff14SB}` | 力場セット: ff19SB（デフォルト）または ff14SB。 | `ff19SB` |

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [all](all.md) -- エンドツーエンドワークフロー（内部で mm-parm を呼び出し）
- [extract](extract.md) -- パラメータ化前に活性部位ポケットを抽出
- [define-layer](define_layer.md) -- トポロジー構築後に ML/MM レイヤーを定義
