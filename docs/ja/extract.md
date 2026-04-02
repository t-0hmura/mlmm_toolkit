# `extract`

## 概要

> **要約:** タンパク質-リガンド PDB から活性部位ポケットを抽出し、ML 領域と周辺 MM 環境を定義します。`-c` で残基名、残基 ID、または PDB パスにより基質を指定します。`--add-linkh` 有効時は切断結合にリンク水素を付加します。非標準残基の電荷には `--ligand-charge` を使用します。

`mlmm extract` は、タンパク質-リガンド PDB から活性部位ポケットを抽出し、ML 領域を定義します。基質周辺の残基を選択し、主鎖/側鎖規則に従ってモデルを切断し、任意で切断結合にリンク水素を付加し、単一構造またはアンサンブルを処理できます。

これは通常、mlmm-toolkit ワークフローの**最初のステップ**であり、完全なタンパク質-リガンド複合体からより小さな計算可能なモデルを生成します。

誤分類が発生した場合（例: 非標準的な残基/原子命名）、以下の付録の命名要件と内部参照リストを参照してください。

## 最小例

```bash
mlmm extract -i complex.pdb -c A:123 -o pocket.pdb -l -3
```

## 出力の見方

- `pocket.pdb`（または `-o` によるカスタムパス）
- INFO で記録される電荷サマリー（アミノ酸、イオン、基質、合計）

## よくある例

1. ID ベースの基質と明示的リガンド電荷によるミニマル実行。

```bash
mlmm extract -i complex.pdb -c A:123 -o pocket.pdb -l -3
```

2. PDB で基質を指定; 残基名ごとの電荷マッピング。

```bash
mlmm extract -i complex.pdb -c substrate.pdb -o pocket.pdb \
 -l "GPP:-3,MMT:-1"
```

3. マルチ構造から単一マルチモデル出力、ヘテロ-ヘテロ近接有効。

```bash
mlmm extract -i complex1.pdb complex2.pdb -c A:123 \
 -o pocket_multi.pdb --radius-het2het 2.6 -l -3 --verbose
```

## 使用法
```bash
mlmm extract -i COMPLEX.pdb [COMPLEX2.pdb...]
 -c SUBSTRATE_SPEC
 [-o POCKET.pdb [POCKET2.pdb...]]
 [--radius Å] [--radius-het2het Å]
 [--include-h2o/--no-include-h2o]
 [--exclude-backbone/--no-exclude-backbone]
 [--add-linkh/--no-add-linkh]
 [--selected-resn LIST]
 [-l, --ligand-charge MAP_OR_NUMBER]
 [--verbose/--no-verbose]
```

### 例
```bash
# ID ベースの基質と明示的リガンド電荷によるミニマル実行
mlmm extract -i complex.pdb -c '123' -o pocket.pdb -l -3

# PDB で基質を指定; 残基名ごとの電荷マッピング（その他は 0）
mlmm extract -i complex.pdb -c substrate.pdb -o pocket.pdb -l 'GPP:-3,SAM:1'

# 名前ベースの基質選択（すべてのマッチを含む; WARNING がログに出力）
mlmm extract -i complex.pdb -c 'GPP,SAM' -o pocket.pdb -l 'GPP:-3,SAM:1'

# マルチ構造から単一マルチモデル出力、ヘテロ-ヘテロ近接有効
mlmm extract -i complex1.pdb complex2.pdb -c 'GPP,SAM' -o pocket_multi.pdb --radius-het2het 2.6 -l 'GPP:-3,SAM:1'

# マルチ構造から複数出力、ヘテロ-ヘテロ近接有効
mlmm extract -i complex1.pdb complex2.pdb -c 'GPP,SAM' -o pocket1.pdb pocket2.pdb --radius-het2het 2.6 -l 'GPP:-3,SAM:1'
```

## ワークフロー

### 残基包含
- `-c/--center` で指定された基質残基は常に含まれます。
- **標準カットオフ（`--radius`、デフォルト 2.6 Å）:**
 - `--no-exclude-backbone` の場合、カットオフ内の任意の原子で残基を包含。
 - `--exclude-backbone` の場合、アミノ酸残基は**非主鎖**原子（N/H*/CA/HA*/C/O 以外）で基質に接触する必要あり。非アミノ酸は任意の原子で可。
- **独立したヘテロ-ヘテロカットオフ（`--radius-het2het`）:** 基質ヘテロ原子（C/H 以外）が指定 Å 以内のタンパク質ヘテロ原子に近接する残基を追加。主鎖除外有効時はタンパク質原子が非主鎖である必要あり。
- **水分子処理:** HOH/WAT/H2O/DOD/TIP/TIP3/SOL はデフォルトで含まれます（`--include-h2o`）。
- **強制包含:** `--selected-resn` は鎖/インサーションコード付き ID を受け付けます（例: `A:123A`）。
- **隣接セーフガード:**
 - 主鎖除外オフで残基が主鎖原子で基質に接触する場合、ペプチド隣接 N/C 残基を自動包含（C-N <= 1.9 Å）。末端はキャップ（N/H* または C/O/OXT）を保持。
 - ジスルフィド結合（SG-SG <= 2.5 Å）は両方の Cys を包含。
 - 非末端 PRO 残基は常に N 側のアミノ酸を包含; 主鎖原子が除去されても CA は保持。`--exclude-backbone` の場合、隣接残基の C/O/OXT はペプチド結合維持のため残留。

### 切断/キャッピング
- 孤立残基は側鎖原子のみ保持; アミノ酸主鎖原子（N, CA, C, O, OXT および N/CA 水素）は PRO/HYP セーフガードを除いて除去。
- 連続ペプチド区間は内部主鎖を保持; 末端キャップ（N/H* または C/O/OXT）のみ除去。TER 認識により鎖切断をまたぐキャッピングを防止。
- `--exclude-backbone` の場合、すべての**非基質**アミノ酸の主鎖原子を削除（PRO/HYP セーフガードと PRO 隣接残基保持に従う）。
- 非アミノ酸残基は主鎖的名称（N/CA/HA/H/H1/H2/H3）の原子を失わない。

### リンク水素（`--add-linkh`）
- 切断結合ベクトル（CB-CA, CA-N, CA-C; PRO/HYP は CA-C のみ）に沿って 1.09 Å のカーボンオンリーリンク水素を追加。
- `TER` 後に連続した `HETATM` レコードとして挿入。残基名 `LKH`（鎖 `L`）の原子名 `HL`。シリアル番号はメインブロックから続行。
- マルチ構造モードでは、同じ結合がすべてのモデルでキャッピングされます; 座標はモデルごとに異なります。

### 電荷サマリー（`--ligand-charge`）
- アミノ酸と一般的なイオンは内部辞書の電荷を使用; 水分子はゼロ。
- 未知残基は `--ligand-charge` が総電荷（未知基質残基に分配、または未知基質がなければ全未知に分配）または残基名別マッピング（`GPP:-3,SAM:1`）を指定しない限りデフォルト 0。
- verbose モード有効時、最初の入力について（タンパク質/リガンド/イオン/合計の）サマリーがログに出力されます。

### 基質指定（`-c/--center`）
- PDB パス: 座標が先頭入力と正確に一致する必要あり（許容 1e-3 Å）; 残基 ID を他の構造に伝播。
- 残基 ID: `'123,456'`、`'A:123,B:456'`、`'123A'`、`'A:123A'`（インサーションコード対応）。
- 残基名: カンマ区切りリスト（大文字小文字を区別しない）。同名残基が複数ある場合、**すべて**含まれ WARNING がログに出力。

### マルチ構造アンサンブル
- 複数入力 PDB を受け付けます（同一原子順序を各ファイルの先頭/末尾で検証）。各構造は独立に処理され、選択残基の**和集合**がすべてのモデルに適用されるため、出力は一貫性を保ちます。
- 出力規則:
 - `-o` なし、複数入力 -> ファイルごとの `pocket_<original_basename>.pdb`。
 - `-o` 1 つ -> 単一マルチ MODEL PDB。
 - N 個の出力（N == 入力数）-> N 個の個別 PDB。
- 診断では、モデルごとの元/保持原子数と残基 ID をエコー。

## CLI オプション

> **注記:** 表示されるデフォルト値は、オプション未指定時に使用されます。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 1 つ以上のタンパク質-リガンド PDB ファイル（同一原子順序が必要）。 | 必須 |
| `-c, --center SPEC` | 基質指定（PDB パス、残基 ID、または残基名）。 | 必須 |
| `-o, --output PATH...` | ポケット PDB 出力。1 パス => マルチ MODEL、N パス => 入力ごと。 | 自動（`pocket.pdb` または `pocket_<input>.pdb`） |
| `-r, --radius FLOAT` | 包含の原子間距離カットオフ (Å)。 | `2.6` |
| `--radius-het2het FLOAT` | 独立したヘテロ-ヘテロカットオフ (Å, C/H 以外)。 | `0.0`（内部でゼロの場合 0.001 Å） |
| `--include-h2o/--no-include-h2o` | HOH/WAT/H2O/DOD/TIP/TIP3/SOL 水分子を含める。 | `True` |
| `--exclude-backbone/--no-exclude-backbone` | 非基質アミノ酸の主鎖原子を除去（PRO/HYP セーフガード）。 | `False` |
| `--add-linkh/--no-add-linkh` | 切断結合に 1.09 Å のカーボンオンリーリンク水素を付加。 | `False` |
| `--selected-resn TEXT` | 強制包含する残基（鎖/インサーションコード付き ID）。 | `""` |
| `-l, --ligand-charge TEXT` | 総電荷または残基名別マッピング（例: `GPP:-3,SAM:1`）。 | _None_ |
| `-v, --verbose/--no-verbose` | INFO レベルのログ出力（`True`）または WARNING のみ（`False`）。 | `True` |

## 出力
```text
<output>.pdb # TER レコード後にリンク水素を含む可能性のあるポケット PDB
 # 単一入力 -> デフォルトで pocket.pdb
 # 複数入力で -o なし -> 構造ごとの pocket_<original_basename>.pdb
 # 複数入力で -o 1 つ -> 単一マルチ MODEL PDB
 # 出力ディレクトリは自動作成されません; 事前に存在を確認してください
```
- verbose モード有効時、モデル #1 の電荷サマリー（タンパク質/リガンド/イオン/合計）がログに出力されます。
- プログラム利用（`extract_api`）では `{"outputs": [...], "counts": [...], "charge_summary": {...}}` を返します。

## MCPB 等で生成された非標準残基を含む系

Amber の `MCPB.py`（Metal Center Parameter Builder）等で金属配位残基のパラメータを生成した場合、金属配位アミノ酸に非標準の残基名（`HD1`, `HE1`, `CM1`, `AP1` 等）が割り当てられます。これらは `extract` の内部辞書 `AMINO_ACIDS` に含まれないため、**主鎖原子の切断・リンク水素の付加が正しく行われません**。

このような残基が検出された場合、`extract` は以下の警告を表示します:

```
[extract] WARNING: Residue HD1 83 may be an amino acid (has N, CA, C, O)
but is not recognized as a standard residue name.
Backbone truncation was not applied.
Consider preparing the pocket model manually.
```

```{important}
非標準残基を含む系では、**ポケットモデルを手動で構築する**ことを推奨します。
手順:

1. 活性部位周辺の残基を選定し、切断箇所を決定する
2. 切断された共有結合の親原子（残る側の原子）に、リンク水素を付加する
3. リンク水素は残基名 `LKH`（チェーン `L`）、原子名 `HL` で記述する
4. 結合方向に沿って **1.09 Å** の位置に配置する
```

## 付録: PDB 命名要件と参照リスト

この付録は主に、**非標準的な残基/原子命名**により `extract` が残基を誤分類する場合のデバッグ用です。標準 PDB 規約に従った入力であれば通常スキップ可能です。

```{important}
`extract` が正しく動作するには、**入力 PDB の残基名と原子名が標準 PDB 命名規約に準拠**している必要があります。ツールは内部辞書を使用してアミノ酸、イオン、水分子、主鎖原子を認識します。非標準命名は残基の誤分類や電荷の不正な割り当てを引き起こします。
```

以下の内部定数が認識される名前を定義します:

### `AMINO_ACIDS`

残基名を公称整数電荷にマッピングする辞書。この辞書への所属が、主鎖処理、切断、電荷計算においてアミノ酸として扱われるかを決定します。

**標準 20 アミノ酸**（電荷は生理的 pH を反映）:
- 中性: `ALA`, `ASN`, `CYS`, `GLN`, `GLY`, `HIS`, `ILE`, `LEU`, `MET`, `PHE`, `PRO`, `SER`, `THR`, `TRP`, `TYR`, `VAL`
- 正電荷 (+1): `ARG`, `LYS`
- 負電荷 (-1): `ASP`, `GLU`

**カノニカル追加:**
- `SEC`（セレノシステイン, 0）、`PYL`（ピロリシン, +1）

**プロトン化/互変異性体バリアント**（Amber/CHARMM 形式）:
- `HIP`（+1, 完全プロトン化 His）、`HID`（0, Nd プロトン化 His）、`HIE`（0, Ne プロトン化 His）
- `ASH`（0, 中性 Asp）、`GLH`（0, 中性 Glu）、`LYN`（0, 中性 Lys）、`ARN`（0, 中性 Arg）
- `TYM`（-1, 脱プロトン化 Tyr フェノレート）

**リン酸化残基:**
- 二価アニオン (-2): `SEP`, `TPO`, `PTR`
- 一価アニオン (-1): `S1P`, `T1P`, `Y1P`
- リン酸化 His (phosaa19SB): `H1D` (0), `H2D` (-1), `H1E` (0), `H2E` (-1)

**システインバリアント:**
- `CYX` (0, ジスルフィド), `CSO` (0, スルフェン酸), `CSD` (-1, スルフィン酸), `CSX` (0, 一般誘導体)
- `OCS` (-1, システイン酸), `CYM` (-1, 脱プロトン化 Cys)

**リシンバリアント / カルボキシル化:**
- `MLY` (+1), `LLP` (+1), `KCX` (-1, Nz カルボキシル酸), `DLY` (+1)

**D-アミノ酸** (19 残基):
- `DAL`, `DAR`, `DSG`, `DAS`, `DCY`, `DGN`, `DGL`, `DHI`, `DIL`, `DLE`, `DLY`, `MED`, `DPN`, `DPR`, `DSN`, `DTH`, `DTR`, `DTY`, `DVA`

**その他の修飾残基:**
- `CGU` (-2, ガンマカルボキシグルタミン酸), `CGA` (-1), `PCA` (0, ピログルタミン酸), `MSE` (0, セレノメチオニン), `OMT` (0, メチオニンスルホン), `HYP` (0, ヒドロキシプロリン)
- その他: `ASA`, `CIR`, `FOR`, `MVA`, `IIL`, `AIB`, `HTN`, `SAR`, `NMC`, `PFF`, `NFA`, `ALY`, `AZF`, `CNX`, `CYF`

**N 末端バリアント**（接頭辞 `N`）: `NALA` (+1), `NARG` (+2), `NASP` (0), `NGLU` (0), `NLYS` (+2) 等、及び `ACE` (0), `NTER` (+1, 汎用)

**C 末端バリアント**（接頭辞 `C`）: `CALA` (-1), `CARG` (0), `CASP` (-2), `CGLU` (-2), `CLYS` (0) 等、及び `NHE` (0), `NME` (0), `CTER` (-1, 汎用)

### `BACKBONE_ATOMS`

アミノ酸の主鎖原子とみなされる原子名のセット。`--exclude-backbone` 時に非基質残基からどの原子を除去するかの判定に使用されます:

```
N, C, O, CA, OXT, H, H1, H2, H3, HN, HA, HA2, HA3
```

### `ION`

イオン残基名を形式電荷にマッピングする辞書。認識されたイオンは電荷サマリーで正しい電荷が自動的に割り当てられます。

| 電荷 | 残基名 |
|------|--------|
| +1 | `LI`, `NA`, `K`, `RB`, `CS`, `TL`, `AG`, `CU1`, `Ag`, `K+`, `Na+`, `NH4`, `H3O+`, `HE+`, `HZ+`, `Tl` |
| +2 | `MG`, `CA`, `SR`, `BA`, `MN`, `FE2`, `CO`, `NI`, `CU`, `ZN`, `CD`, `HG`, `PB`, `Be`, `PD`, `PT`, `Sn`, `Ra`, `YB2`, `V2+` |
| +3 | `FE`, `AU3`, `AL`, `GA`, `IN`, `CE`, `Ce`, `CR`, `Cr`, `Dy`, `EU`, `EU3`, `Er`, `GD3`, `LA`, `LU`, `Nd`, `PR`, `SM`, `Sm`, `TB`, `Tm`, `Y`, `Pu` |
| +4 | `U4+`, `Th`, `Hf`, `Zr` |
| -1 | `F`, `CL`, `BR`, `I`, `Cl-`, `IOD` |

### `WATER_RES`

水分子として認識される残基名のセット。水分子はデフォルトで含まれ（`--include-h2o`）、電荷はゼロが割り当てられます:

```
HOH, WAT, H2O, DOD, TIP, TIP3, SOL
```

---

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [はじめに](getting-started.md) -- インストールと初回実行
- [概念とワークフロー](concepts.md) -- 全系 vs. ML 領域
- [CLI 規約](cli-conventions.md) -- 残基セレクタと電荷指定
- [mm-parm](mm-parm.md) -- 抽出したポケットから Amber トポロジーを生成
- [define-layer](define-layer.md) -- 3 層 ML/MM 分割の付与
