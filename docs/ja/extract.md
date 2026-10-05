# `extract`（ML 領域の切り出し）

## 概要

`extract` は、タンパク質–リガンドの PDB/mmCIF から基質の周りの残基を切り出し、決まった規則で主鎖を切って、切り出した領域の電荷を数えます。mlmm-toolkit では、この領域が ML 領域になり、残りのタンパク質は MM として計算に残ります。

### 主な用途

* **ML 領域を決める**: `define-layer` と計算コマンドが `--model-pdb` で受け取る原子の選択を書き出します。
* **計算の前に ML 領域を確かめる**: どの残基が入るか、原子数、電荷を見ます。
* **複数の状態を同じ境界で切る**: 原子の並びが同じ反応物と生成物を 1 回で渡すと、どの出力も同じ残基になります。
* **非標準残基を扱う**: MCPB.py などが付けた残基名を、`--modified-residue` でアミノ酸として登録します。

`all` に `-c` を付けると、`all` が `extract` を実行します。ML 領域の大きさの決め方は [ML 領域と層の組み方](model-setup.md) を見てください。

---

## 基本的な実行例

### 1. 残基 ID と総電荷で選ぶ

基質を chain:残基名:番号 で、その総電荷を 1 つの数で渡します。

```bash
mlmm extract -i complex.pdb -c 'A:GPP:301' -o pocket.pdb -l -3 --out-json
```

端末に `[extract] Atoms after truncation: N` が出ます。領域の原子数は N です。電荷は `[extract] Total active site model charge` の行に出ます。`result.json` の `n_atoms_extracted`・`total_charge` にも、同じ数が入ります。`pocket.pdb` をビューアで開き、反応に関わる残基が入っているかを確かめてください。

### 2. 基質を PDB ファイルで渡す

基質の PDB ファイルを中心にし、残基名ごとに電荷を渡します。

```bash
mlmm extract -i complex.pdb -c substrate.pdb -o pocket.pdb -l 'GPP:-3,SAM:1'
```

基質のファイルの座標は、複合体の座標と 0.001 Å 以内で一致している必要があります。

### 3. 残基名で選ぶ

残基名を並べると、その名前の残基がすべて中心になります。

```bash
mlmm extract -i complex.pdb -c 'GPP,SAM' -o pocket.pdb -l 'GPP:-3,SAM:1'
```

### 4. 複数の構造を 1 回で切る

反応物と生成物を 1 つの `-i` の後に並べると、両方が同じ残基になり、1 つのマルチ MODEL の PDB に書き出されます。

```bash
mlmm extract -i complex_R.pdb complex_P.pdb -c 'A:GPP:301,A:SAM:302' \
    -o pocket_multi.pdb -l 'GPP:-3,SAM:1'
```

(ja-extract-modified-residue)=
### 5. 非標準残基（`--modified-residue`）

Amber の MCPB.py などは、金属に配位する残基に非標準の名前（`HD1`、`HE1`、`CM1`、`AP1`）を付けます。`extract` はこの名前を知らないので、主鎖を切らず、次の警告を出します。

```text
[extract] WARNING: Residue HD1 83 may be an amino acid (has N, CA, C, O) but is not recognized as a standard residue name. Backbone truncation was not applied. Consider preparing the active site model manually.
```

この名前を `NAME:charge` の形で、電荷と一緒にアミノ酸として登録します。

```bash
mlmm extract -i complex.pdb -c 'A:SUB:301' -o pocket.pdb \
    --modified-residue 'HD1:0,HE1:0,CM1:0,AP1:0'
```

---

## 処理の仕組みと計算仕様

1. **中心**: `-c` に基質・補因子・金属を並べます。各項目は [残基の指定](cli-conventions.md#残基セレクタ) で、推奨の形は `A:TYR:44`（chain:残基名:番号）です。`A:SAM`、`SAM` のような名前、番号、基質の PDB/mmCIF ファイルも使えます。`--selected-resn` は同じ形で残基を足し、そこからは距離での探索を始めません。
2. **隣の残基**: 中心の原子から `-r`（既定 2.6 Å）以内、または両側とも C・H 以外の原子どうしで `--radius-het2het` 以内に原子がある残基を入れます（どちらか一方で足ります）。水は `--no-include-h2o` を付けない限り数え、`--exclude-backbone` ではアミノ酸の主鎖の原子による接触を数えません。さらに、選んだシステインの S–S 結合の相手（S–S ≤ 2.5 Å）、選んだプロリンの N 側の隣、`--exclude-backbone` でないときは主鎖の原子が中心に触れたアミノ酸のペプチド結合の両隣を足します。
3. **主鎖の切断**: つながったアミノ酸の並びは内側の主鎖を残し、両端が CA で終わるように切ります。1 残基だけで入った残基は側鎖だけになります。`-c` のアミノ酸も同じ規則で切りますが、ペプチド結合でつながった隣の残基が `-r` 以内に入って領域に加わるので、全部の原子が残ります。プロリンは環を残します。`--exclude-backbone` では、ペプチド結合でつながった `-c` のアミノ酸どうしの間を除き、アミノ酸の主鎖の原子をすべて除きます。水とアミノ酸でない残基は切りません。
4. **電荷**: アミノ酸とイオンは組み込みの表から、水は 0、ほかの残基は `-l` で渡さない限り 0 とします（次の小節）。
5. **キャップ水素（`--add-linkh` のときだけ）**: 切断で CA か CB の結合相手が無くなった所（CB–CA、CA–N、CA–C。プロリンは CA–C だけ）に、その炭素から元の結合の向きに 1.09 Å の位置へ水素を置きます。キャップ水素は `TER` の後に、残基 `LKH`・chain `L` の `HETATM` 原子 `HL` として書かれます。`--model-pdb` のファイルには入れないでください。ML/MM 計算機が、ML/MM の境界をまたぐ `parm7` の結合に自分でリンク水素を付けるので、キャップ付きのポケットを単独で使うとき以外は `--add-linkh` を付けません。

境界を自分で決めて確かめるときは、{ref}`信頼できる model.pdb の作り方 <ja-model-pdb-selection>` を見てください。

### 電荷の内訳

`-l` には `'GPP:-3,SAM:1'` のような残基名ごとの電荷か、1 つの数を渡します。未知の残基とは、付録でアミノ酸・イオン・水のどれにも挙がっていない残基です。数を渡すと、`-c` の中の未知の残基に均等に割り、`-c` に未知の残基が無ければ、すべての未知の残基に割ります。`-l -3` を 2 残基に割ると −1.5 ずつになり、合計は −3 のままです。残基名で渡したときは、書かなかった未知の残基は 0 です。既定の詳細度で、端末にはタンパク質・リガンド・イオンの電荷に続いて `Total active site model charge` が出ます。入力が複数のときは、最初の入力の内訳です。

### 複数の構造

入力が複数のときは、構造ごとに残基を選び、その和集合をすべての構造に当てるので、どの出力も同じ原子になります。座標はモデルごとのものです。端末には、モデルごとに `[extract:multi] Atoms after truncation (model k): N` が出ます。

---

## 主な出力ファイル

```text
./
├─ pocket.pdb    # 切り出した領域（キャップ水素は --add-linkh のときだけ、TER の後）
├─ pocket.cif    # mmCIF の入力か、PDB の桁に収まらない PDB の入力のとき
├─ result.json   # --out-json のとき。最初の出力ファイルと同じディレクトリ
└─ summary.json  # result.json の写し。result.json を読む（--out-json のとき）
```

| 入力 | `-o` | 出力 |
| --- | --- | --- |
| 1 つ | なし | `pocket.pdb` |
| 複数 | なし | 入力ごとに `pocket_<入力の名前>.pdb` |
| 複数 | 1 つ | マルチ MODEL の PDB 1 つ |
| 複数 | 入力と同じ数 | 入力ごとに PDB 1 つ |

`-o` がこれ以外の数のときと、出力先が入力のファイルそのもののときは、エラーで止まります。出力先の親ディレクトリは自動で作られます。`result.json` には、原子数（`n_atoms_raw`、`n_atoms_extracted`、`n_link_hydrogens`）、電荷（`total_charge`、`protein_charge`、`ligand_total_charge`、`ion_total_charge`）と使った設定が入ります（[JSON 出力リファレンス](json-output.md)）。mmCIF の入力と、PDB の桁に収まらない大きな PDB の入力では、元の ID のままの `.cif` も出ます（{ref}`mmCIF の入力 <ja-mmcif-input>`）。

---

## 主な CLI オプション

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス（複数可） | （必須） | タンパク質–リガンドの PDB/mmCIF。1 つの `-i` の後に並べても、`-i` を繰り返してもよい。原子が同じ順に並んでいること |
| `-c, --center` | 文字列 | （必須） | 中心の残基か、基質の PDB/mmCIF ファイル（例: `'A:TYR:44,A:SAM:301'`） |
| `-o, --output` | パス（複数可） | 上の表 | 出力する PDB のパス |
| `-r, --radius` | 浮動小数点数 | `2.6` | 中心の原子からの距離のしきい値（Å）。`0` では距離で隣の残基を足さない（[使用上の注意点](#使用上の注意点)） |
| `--radius-het2het` | 浮動小数点数 | `0`（無効） | C・H 以外の原子どうしの 2 つ目のしきい値（Å） |
| `--selected-resn` | 文字列 | `""` | 距離で探さずに足す残基。`-c` と同じ形 |
| `--include-h2o/--no-include-h2o` | フラグ | `True` | 水（HOH、WAT、H2O、DOD、TIP、TIP3、SOL）を入れる |
| `--exclude-backbone/--no-exclude-backbone` | フラグ | `False` | アミノ酸から主鎖の原子を除く。ペプチド結合でつながった `-c` のアミノ酸どうしの間は残す |
| `--add-linkh/--no-add-linkh` | フラグ | `False` | 切断で CA か CB の相手が無くなった所にキャップ水素を付ける。`--model-pdb` に使うファイルでは付けない |
| `--modified-residue` | 文字列 | `""` | アミノ酸として扱う残基名。`NAME:charge`（電荷なしの `NAME` は組み込みの表にある名前だけ） |
| `-l, --ligand-charge` | 文字列 | `None` | 総電荷か、残基名ごとの電荷（例: `'GPP:-3,SAM:1'`） |
| `--out-json/--no-out-json` | フラグ | `False` | `result.json` と `summary.json` を書き出す |

全オプションの一覧は [自動生成 CLI リファレンス](../reference/commands/extract.md) を参照してください。

---

## 使用上の注意点

* **`-r 0`** では距離で隣の残基を足さず（内部では 0.001 Å）、`-c` と `--selected-resn` の残基に、手順 2 で足す S–S 結合の相手とプロリンの N 側の隣を加えた領域になります。`--radius-het2het 0` も同じです。
* **領域の大きさ**: ML 領域を広げても結果が変わらないことを系ごとに確かめてください。`-r` を大きくすると計算は重くなり、精度が上がるとは限りません（{ref}`モデルを広げる <ja-model-setup-larger>`）。
* **使い回す `--model-pdb` を手で作るとき**は、`mm-parm` が書く PDB から切り出して、原子を `parm7` とそろえてください（[mm-parm の例 4](mm-parm.md#基本的な実行例)）。
* **名前はすべての chain に当たる**: `TYR` のような名前は、どの chain の TYR もすべて選び、複数あれば警告を出します。
* **`TYR:44` は chain TYR と読まれる**: 2 つの欄では最初の欄が必ず chain で、2 つ目は番号（`TYR:44`）か名前（`A:SAM`）なので、`A:TYR:44` と書いてください。同梱例のように chain の欄が空の PDB では、名前か番号だけを使います。
* **1 つの list に 1 つの形**: `'SAM,44'` のように名前と番号を混ぜた list はエラーで止まります。
* **境界の警告**: 手順 5 の CA と CB の切断のほかで、非金属の原子どうしの結合を切ると、`extract` は警告を出します。計算の前に、境界・電荷・多重度を確かめてください。境界の選び方は [ML 領域と層の組み方](model-setup.md) にあります。
* **どの入力も同じ原子**: 原子の数や並びが違う入力は `[multi] Atom count mismatch` か `[multi] Atom order mismatch` で止まります（[トラブルシューティング](troubleshooting.md)）。
* **マルチ MODEL の入力**は、最初の MODEL だけを使い、警告を出します。
* **altLoc（別位置の配座）**: `extract` は残基ごとに、平均占有率がいちばん高い altLoc を 1 つ残します。整えたファイルそのものが要るときは [fix-altloc](fix-altloc.md) を使ってください。
* **`--modified-residue`**: 電荷を書かない `NAME` を使えるのは、組み込みの表（付録）にある名前だけで、表の電荷のままになります（`SEP` なら −2）。表に無い名前を電荷なしで渡すと、`NAME:charge` で書くよう求めるエラーで止まります。`NAME:charge` は組み込みの電荷もこの実行に限って上書きします（中性の Lys なら `LYS:0`）。`--modified-residue` で足りないときは、ML 領域の原子を手で選びます。
* **組み込みの残基名**は Amber/CHARMM の命名です。PDB の残基が別の化合物と同じ名前を持つときは、`--modified-residue NAME:charge` で意図する電荷を渡してください。

---

## 関連ドキュメント

* [ML 領域と層の組み方](model-setup.md) — ML 領域を削る・広げる、MM の層を決める、原子を固定する
* [all](all.md) — 一括のワークフロー。`-c` で `extract` を実行する
* [mm-parm](mm-parm.md) — Amber の topology と、切り出し元にする対応した PDB を作る
* [define-layer](define-layer.md) — 切り出した領域から ML・Movable-MM・Frozen-MM の層を付ける
* [fix-altloc](fix-altloc.md) — 残基ごとに別位置の配座を 1 つにした PDB を書く
* [add-elem-info](add-elem-info.md) — 切り出しの前に元素の欄を埋める
* [共通オプションと残基・原子の指定](cli-conventions.md) — 残基の指定と電荷
* [トラブルシューティング](troubleshooting.md) — 切り出しのエラー
* [用語集](glossary.md) — ML 領域、リンク原子

## 付録: PDB 命名規則と参照リスト

PDB の残基名や原子名が標準でないために、`extract` が残基を取り違えたり電荷を誤ったりするときに、この付録を使ってください。標準の PDB の名前なら読み飛ばしてかまいません。

```{important}
`extract` は、PDB の残基名と原子名でアミノ酸・イオン・水・主鎖の原子を見分けます。入力は標準の PDB の化合物名に従っている必要があり、ほかの名前では残基の取り違えや電荷の誤りが起きます。
```

### アミノ酸

アミノ酸として扱う残基名と、その名目の電荷です。主鎖の切断とアミノ酸の電荷は、この残基にだけ適用します。

**標準 20 アミノ酸**（電荷は生理的 pH）:

- 中性: `ALA`, `ASN`, `CYS`, `GLN`, `GLY`, `HIS`, `ILE`, `LEU`, `MET`, `PHE`, `PRO`, `SER`, `THR`, `TRP`, `TYR`, `VAL`
- 正電荷 (+1): `ARG`, `LYS`
- 負電荷 (−1): `ASP`, `GLU`

**追加の標準アミノ酸:** `SEC`（セレノシステイン, 0）、`PYL`（ピロリシン, 0）

**プロトン化状態・互変異性体**（Amber/CHARMM の命名）: `HIP`（+1, 完全にプロトン化した His）、`HID`（0, Nδ がプロトン化した His）、`HIE`（0, Nε がプロトン化した His）、`ASH`（0, 中性の Asp）、`GLH`（0, 中性の Glu）、`LYN`（0, 中性の Lys）、`ARN`（0, 中性の Arg）、`TYM`（−1, 脱プロトン化した Tyr のフェノラート）

**リン酸化残基:** 2 価のアニオン (−2) `SEP`, `TPO`, `PTR`。1 価のアニオン (−1) `S1P`, `T1P`, `Y1P`。リン酸化 His（phosaa19SB）`H1D` (0), `H2D` (−1), `H1E` (0), `H2E` (−1)

**システインの変異体:** `CYX`（0, ジスルフィド）、`CSO`（0, スルフェン酸）、`CSD`（−1, スルフィン酸）、`CSX`（0, 一般の誘導体）、`OCS`（−1, システイン酸）、`CYM`（−1, 脱プロトン化した Cys）

**リシンの変異体・カルボキシル化:** `MLY` (+1), `LLP` (0), `KCX`（−1, Nz のカルボン酸）

**D-アミノ酸**（19 残基）: `DAL`, `DAR`, `DSG`, `DAS`, `DCY`, `DGN`, `DGL`, `DHI`, `DIL`, `DLE`, `DLY`, `MED`, `DPN`, `DPR`, `DSN`, `DTH`, `DTR`, `DTY`, `DVA`

**その他の修飾残基:** `CGU`（−2, γ-カルボキシグルタミン酸）、`CGA` (−1)、`PCA`（0, ピログルタミン酸）、`MSE`（0, セレノメチオニン）、`OMT`（0, メチオニンスルホン）、`HYP`（0, ヒドロキシプロリン）。ほかに `ASA`, `CIR`, `FOR`, `MVA`, `IIL`, `AIB`, `HTN`, `SAR`, `NMC`, `PFF`, `NFA`, `ALY`, `AZF`, `CNX`, `CYF`（いずれも 0）

**N 末端の変異体**（接頭辞 `N`）: `NALA` (+1), `NARG` (+2), `NASP` (0), `NGLU` (0), `NLYS` (+2) など。ほかに `ACE` (0), `NTER`（+1, 汎用）

**C 末端の変異体**（接頭辞 `C`）: `CALA` (−1), `CARG` (0), `CASP` (−2), `CGLU` (−2), `CLYS` (0) など。ほかに `NHE` (0), `NME` (0), `CTER`（−1, 汎用）

接頭辞 `N`・`C` の Amber の名前は標準の残基として読み（`NALA` → `ALA`）、末端の電荷は、モデルに N 末端の H1〜H3（プロリンは H2 と H3）か OXT が残るときだけ数えます。

### 主鎖の原子

アミノ酸の主鎖として扱う原子名です。`--exclude-backbone` では、ペプチド結合でつながった `-c` のアミノ酸どうしの間を除き、この原子を除きます。

```
N, C, O, CA, OXT, H, H1, H2, H3, HN, HA, HA2, HA3
```

### イオン

イオンとして扱う残基名と、その形式電荷です。

| 電荷 | 残基名 |
|------|--------|
| +1 | `LI`, `NA`, `K`, `RB`, `CS`, `TL`, `AG`, `CU1`, `K+`, `NA+`, `NH4`, `H3O+`, `H3O`, `HE+`, `HZ+` |
| +2 | `MG`, `CA`, `SR`, `BA`, `MN`, `FE2`, `CO`, `NI`, `CU`, `ZN`, `CD`, `HG`, `PB`, `BE`, `PD`, `PT`, `SN`, `RA`, `YB2`, `V2+` |
| +3 | `FE`, `AU3`, `AL`, `GA`, `IN`, `CE`, `CR`, `DY`, `EU`, `EU3`, `ER`, `GD3`, `LA`, `LU`, `ND`, `PR`, `SM`, `TB`, `TM`, `Y`, `PU` |
| +4 | `U4+`, `TH`, `HF`, `ZR` |
| −1 | `F`, `CL`, `BR`, `I`, `CL-`, `IOD` |

### 水

水として扱う残基名です（既定の `--include-h2o` で入り、電荷は 0）。

```
HOH, WAT, H2O, DOD, TIP, TIP3, SOL
```
