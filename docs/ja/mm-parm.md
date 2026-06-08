# `mm-parm`

`mlmm mm-parm` は PDB から Amber トポロジー/座標ファイル（parm7/rst7/pdb）を生成します。不明な残基は GAFF2（AM1-BCC 電荷）で自動的にパラメータ化され、ジスルフィド結合検出や PDBFixer による任意の水素付加にも対応します。パイプライン全体は「処理の流れ」、力場や水素付加のフラグは「CLI オプション」を参照してください。金属酵素、糖鎖、非標準アミノ酸・翻訳後修飾、MD スナップショット入力には不向きで、これらは外部で用意したトポロジーを `--parm` から供給します（[注記](#注記)を参照）。

## 実行例

基本的なトポロジー構築（リガンドの電荷・多重度を指定）。

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 -l "GPP=-3,MMT=-1" --ligand-mult "GPP=1,MMT=1"
```

TER レコード追加、ff19SB、pH 7 での水素付加。

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 -l "GPP=-3,MMT=-1" --ligand-mult "GPP=1,MMT=1" \
 --add-ter --ff-set ff19SB --add-h --ph 7.0
```

水素付加をスキップ（入力がすでにプロトン化済み）。

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 -l "GPP=-3" --no-add-h
```

## 処理の流れ

1. **入力準備** -- 入力 PDB はそのまま読み込まれます（構造修正なし）。`--add-h` が設定されている場合、PDBFixer で指定 `--ph` にて水素が付加されます。
2. **TER 挿入** -- `--add-ter`（デフォルト）の場合、リガンド/水/イオン残基の連続ブロックの前後に TER レコードが挿入されます。
3. **不明残基のパラメータ化** -- 力場で認識されない残基は antechamber（GAFF2、AM1-BCC）と parmchk2 で自動パラメータ化されます。`--ligand-charge` に名前を列挙した残基はこの経路が最優先されます。形式電荷とスピン多重度は `--ligand-charge` と `--ligand-mult` で制御されます。
4. **ジスルフィド検出** -- CYS/CYM/CYX ペアで SG-SG（または S-S）距離が 2.3 Å 以下のものが自動的に結合されます。
5. **トポロジー構築** -- tleap が選択された力場セットを使用して parm7/rst7/pdb ファイルを生成します。

## 出力

- `<prefix>.parm7` -- Amber prmtop トポロジー
- `<prefix>.rst7` -- Amber ASCII inpcrd 座標
- `<prefix>.pdb` -- LEaP savepdb 出力

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力 PDB（`--add-h` でない限りそのまま使用）。 | 必須 |
| `--out-prefix TEXT` | parm7/rst7/pdb ファイルの出力接頭辞。 | 入力 PDB のファイル名幹 |
| `-l, --ligand-charge TEXT` | 残基名と形式電荷のマッピング（例: `"GPP=-3,MMT=-1"`）。 | _None_ |
| `--ligand-mult TEXT` | 残基名とスピン多重度のマッピング（例: `"HEM=1,NO=2"`）。未指定の残基はデフォルトで一重項（1）。 | _None_ |
| `--keep-temp/--no-keep-temp` | 作業ディレクトリの中間ファイル/ログを保持（デバッグ用）。 | `False` |
| `--add-ter/--no-add-ter` | リガンド/水/イオンブロックの前後に TER を挿入。 | `True` |
| `--add-h/--no-add-h` | PDBFixer で `--ph` に基づいて水素を付加。 | `False` |
| `--ph FLOAT` | PDBFixer の水素付加用 pH（`--add-h` の場合のみ使用）。 | `7.0` |
| `--ff-set {ff19SB\|ff14SB}` | 力場セット: ff19SB（デフォルト）または ff14SB。 | `ff19SB` |

全フラグの一覧は生成された[コマンドリファレンス](../reference/commands/index.md)にあります。

## 注記

`mm-parm` は AmberTools の tleap と GAFF2 自動パラメータ化に依存しており、基質が**典型的な有機分子**である場合にうまく機能します。以下のケースでは、外部でトポロジーを自作し（例: tleap, MCPB.py, glycam.org ツール）、各サブコマンドの `--parm` フラグから入力することを強く推奨します。

- **金属酵素** -- 金属中心には専用の結合/非結合パラメータが必要です（例: MCPB.py, bonded model, ZAFF）。GAFF2 の自動パラメータ化では金属-配位子の配位を扱えません。
- **糖鎖（Glycan）を含む系** -- 糖鎖結合には GLYCAM 力場パラメータが必要であり、標準の GAFF2/ff19SB セットアップには含まれていません。
- **非標準アミノ酸・翻訳後修飾** -- リン酸化、メチル化などの修飾残基にはカスタム `frcmod`/`lib` ファイルが必要な場合があります。
- **MD スナップショットからの初期構造** -- MD トラジェクトリのスナップショットから出発する場合、MD シミュレーションで使用した同じ `.parm7` ファイルを再利用するのが最も妥当です。これにより、ML/MM で使用する MM エネルギー面と前段の MD で使用したものとの間の整合性が保たれ、再パラメータ化による人工的な差異（例: 部分電荷や原子タイプの割り当ての違い）を回避できます。
- `AMINO_ACIDS` に載っているアミノ酸系残基で、選択された力場に認識されないものは自動処理されません。手動でパラメータ化するようメッセージを表示してビルドを中断します。
- `--ff-set ff14SB` を使用すると、力場は ff14SB（タンパク質）+ TIP3P（水）（+ phosaa14SB）に切り替わります。それ以外の場合はデフォルトの `ff19SB` セットが使用されます。

```bash
# 例: MD で構築済みのトポロジーを供給
mlmm opt -i snapshot_layered.pdb --parm md_system.parm7 -q -1 -m 1 \
  --opt-mode grad --out-dir result
```

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [all](all.md) -- 一気通貫ワークフロー（内部で mm-parm を呼び出し）
- [extract](extract.md) -- パラメータ化前に活性部位ポケットを抽出
- [define-layer](define-layer.md) -- トポロジー構築後に ML/MM レイヤーを定義
