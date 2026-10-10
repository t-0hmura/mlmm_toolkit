# `mm-parm`（Amber トポロジーの作成）

`mm-parm` は、酵素–基質複合体の全系の PDB から、AmberTools の tleap で Amber のトポロジー（`parm7`）・座標（`rst7`）と、それに対応する PDB を作ります。基質や補因子のように力場が知らない残基には、GAFF2 のパラメータと AM1-BCC 電荷を付けます。ML/MM の計算コマンドは、どれもこの `parm7` を `--parm7` で読みます。

---

## 主な用途

* **MM のトポロジーを作る**: 全系の `parm7`・`rst7`・PDB を書き出します。
* **リガンドのパラメータを作る**: 未知の残基に、`-l` と `--ligand-mult` の形式電荷と多重度で GAFF2 のパラメータを付けます。
* **モデルを手で組む**: `mm-parm` が書く PDB を `extract` と `define-layer` に渡します。この PDB の原子は `parm7` と同じ順に並んでいます。

---

## 基本的な実行例

### 1. リガンドの電荷と多重度を渡す

リガンドごとの形式電荷とスピン多重度を残基名で渡します。

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
    -l 'GPP:-3,MMT:-1' --ligand-mult 'GPP:1,MMT:1'
```

端末に `complex.pdb`・`complex.parm7`・`complex.rst7` の `[mm-parm] Wrote:` が出れば成功です。

### 2. pH 7.0 で水素を付ける

作成の前に PDBFixer で水素を付けます。

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
    -l 'GPP:-3,MMT:-1' --ligand-mult 'GPP:1,MMT:1' \
    --add-ter --ff-set ff19SB --add-h --ph 7.0
```

### 3. 水素が付いている入力

入力をそのまま使います。

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
    -l 'GPP:-3' --no-add-h
```

### 4. モデルを手で組む

トポロジーを作り、`mm-parm` が書いた PDB から ML 領域を切り出し、同じ PDB に層を付けます。

```bash
mlmm mm-parm -i input.pdb -l 'LIG:0' --out-prefix system
mlmm extract -i system.pdb -c LIG -l 'LIG:0' -o model.pdb
mlmm define-layer -i system.pdb --model-pdb model.pdb -o system_layered.pdb
```

計算コマンドには `system_layered.pdb` と `--parm7 system.parm7` を渡します。

### 5. mmCIF の入力に使うトポロジー

`mm-parm` が読むのは PDB だけです。mmCIF の入力では、全系を原子の並びと元素を変えずに PDB に書き出し、そこからトポロジーを作って、その `parm7` を mmCIF の構造と一緒に `all` に渡します。

```bash
mlmm mm-parm -i reactant_topology.pdb -l 'SAM:1,GPP:-3' \
    --out-prefix full_system
mlmm all -i reactant.cif product.cif --parm7 full_system.parm7 \
    -c 'enzyme_A:SAM:10001,enzyme_A:GPP:10002' \
    -l 'SAM:1,GPP:-3' --tsopt --thermo -o result
```

---

## 処理の仕組みと計算仕様

1. **入力**: PDB をそのまま使います。`--add-h` では、PDBFixer が `--ph` で水素を付けます。欠けた重原子や残基は足しません。
2. **TER レコード**: `--add-ter`（デフォルト）では、`-l` に書いた残基・水・イオンのまとまりの前後に `TER` を入れます。こうした残基が続く所は分けません。ペプチドの C–N 結合でつながっていない隣のアミノ酸の間（別の chain か、C–N > 1.9 Å）にも `TER` を入れます。
3. **ジスルフィド結合**: SG 原子どうしが 2.5 Å 以内の CYS/CYX の組を結合し、結合した CYS の名前を CYX に変えて、tleap が HG を外すようにします。`--no-auto-disulfide` では、もとから CYX という名前の残基だけを結合します。
4. **未知の残基**: まず力場だけで tleap を実行します。tleap が未知と報告した残基名ごとに、ファイルの中でその名前の最初の残基から、antechamber（GAFF2、AM1-BCC）と parmchk2 でパラメータを作ります。電荷は `-l`（無ければ 0）、多重度は `--ligand-mult`（無ければ 1）です。antechamber の前に、残基の電子数がこの電荷・多重度と合うかを確かめます。
5. **トポロジー**: 新しいパラメータを読んで tleap をもう一度実行し、トポロジー・座標・PDB を書きます。`mm-parm` は PDB の空の元素の欄を `parm7` から埋め、ほかのレコードと原子の順は変えません。

---

## 主な出力ファイル

```text
./
├─ <prefix>.parm7   # Amber のトポロジー
├─ <prefix>.rst7    # Amber の座標（ASCII）
└─ <prefix>.pdb     # 元素の欄を埋めた tleap の PDB。原子と順序は parm7 と同じ
```

`<prefix>` のデフォルト値は入力のファイル名から拡張子を除いたもので、現在のディレクトリに書かれます。PDB は `--out-prefix` を指定したときに書かれ、`--out-prefix` なしで `--add-h` を指定したときは `<入力の名前>_parm.pdb` になります。どちらでもなければ `parm7` と `rst7` だけです。`<prefix>.pdb` が入力を置き換えないよう、接頭辞には入力と違う名前を指定してください。`--add-h` の後で作成が失敗したときは、その PDB のパスにまだファイルが無ければ、水素を付けた構造をそこに書きます。`--keep-temp` では、tleap のログを含む作業ディレクトリ `parm7build_*` が現在のディレクトリに残ります。

---

## 主な CLI オプション

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | 入力 PDB。`--add-h` が無ければそのまま使う |
| `-o, --out-prefix` | 文字列 | 入力の名前（拡張子なし） | 出力ファイルの接頭辞 |
| `-l, --ligand-charge` | 文字列 | `None` | 残基名ごとの形式電荷（例: `'GPP:-3,MMT:-1'`） |
| `--ligand-mult` | 文字列 | `1` | 残基名ごとのスピン多重度（例: `'HEM:1,NO:2'`） |
| `--keep-temp/--no-keep-temp` | フラグ | `False` | 作業ディレクトリと tleap のログを残す |
| `--add-ter/--no-add-ter` | フラグ | `True` | リガンド・水・イオンのまとまりの前後と、ペプチド結合でつながらないアミノ酸の間に `TER` を入れる |
| `--auto-disulfide/--no-auto-disulfide` | フラグ | `True` | SG–SG ≤ 2.5 Å の CYS/CYX の組を結合し、結合した CYS の名前を CYX に変える。オフでは、もとから CYX の残基だけを結合する |
| `--add-h/--no-add-h` | フラグ | `False` | PDBFixer で、`--ph` の pH に合わせて水素を付ける |
| `--ph` | 浮動小数点数 | `7.0` | `--add-h` の pH |
| `--ff-set` | `ff19SB` か `ff14SB` | `ff19SB` | [力場の組](#使用上の注意点) |

全オプションは [自動生成のオプションの一覧（英語のみ）](../reference/commands/mm_parm.md) を参照してください。

---

## `oniom-export` 用の CMAP を含まないトポロジー

ff14SB で作成し、CMAP 項がないことを確認します。

```bash
mlmm mm-parm -i input.pdb -l 'LIG:0' --ff-set ff14SB --out-prefix system
python -c "import parmed as pmd; p=pmd.load_file('system.parm7'); assert not p.cmaps"
```

---

## 使用上の注意点

* **トポロジーを自分で作る系**: `mm-parm` は、基質が典型的な有機分子のときに向いています。次の系では、トポロジーを自分で作り、`--parm7` で渡してください。
  * **金属酵素**: 金属中心には専用の結合・非結合パラメータ（MCPB.py、bonded model、ZAFF）が要ります。GAFF2 では金属–配位子の配位を表せません。
  * **糖鎖**: 力場の組は GLYCAM_06j-1 を読み込みますが、`mm-parm` はジスルフィド以外の結合を作りません。グリコシド結合など残基の間の共有結合には、tleap の `bond` コマンドを自分で書く必要があります。
  * **非標準アミノ酸・翻訳後修飾**: 修飾残基には専用の `frcmod`/`lib` ファイルが要ることがあります。
  * **MD から取った構造**: MD の `parm7` を使い回してください。そうすれば、ML/MM の計算は MD と同じ MM のエネルギー面を使い、部分電荷や原子タイプも変わりません。

  ```bash
  # MD で作ったトポロジーを使う
  mlmm opt -i snapshot_layered.pdb --parm7 md_system.parm7 -q -1 -m 1 \
    --opt-mode grad --out-dir result
  ```
* **`parm7` は原子の順で対応づける**: どの座標の入力も、全系の原子を `parm7` と同じ順に持っている必要があります。ML/MM 計算機は計算の前に、原子数と、原子ごとの元素・原子名（`1HB` と `HB1` は同じ）・残基名・残基の順番を比べ、最初に食い違った所で止まります。反応物・中間体・生成物の構造を作るときは、原子名・残基名・原子の順を保ってください。`--model-pdb` のファイルは、全系から変えずに取り出した部分集合で、ML 領域の原子を選ぶためだけに使います。
* **出力 PDB の残基番号**: tleap は残基を出てくる順に 1, 2, … と数えるので、`mm-parm` が書く PDB の残基番号は入力と違うことがあります。`examples/beza/1.R.pdb` では、ARG 38 が 1 番目、SAM 320 が 283 番目の残基です。この PDB を手で切るときは、残基を名前で選ぶか、番号をファイルで確かめてください。`all` は元の入力を切るので、元の番号のまま指定できます。
* **力場が知らないアミノ酸**: `extract` がアミノ酸として扱う残基（`extract` の付録）を tleap が知らないと、作成は止まります。メッセージは 3 つの対処を示します。`-l` に書いて GAFF2 のパラメータを付ける、入力の残基を変える、tleap でトポロジーを自分で作る、です。
* **リガンドの電荷と水素**: 残基の水素の数が電荷・多重度と合わないと、電子数の確認で止まります（SAM では水素 22 個が電荷 0、23 個が +1）。tleap がもう知っている残基への `-l`・`--ligand-mult` は使われず、警告が出ます。
* **力場の組**: `ff19SB` は ff19SB と phosaa19SB・ff19SB_modAA、OPC3 の水とそのイオンのパラメータを読み込みます。`ff14SB` は ff14SB と phosaa14SB・ff14SB_modAA、TIP3P の水とそのイオンのパラメータを読み込みます。どちらも lipid21・RNA.OL3・DNA.OL21・GLYCAM_06j-1・GAFF2 を読み込みます。
* **仮想サイトを持つ水**: デフォルトの MM バックエンド `hessian_ff` は、質量の無い仮想サイトを持つ水（OPC・TIP4P/-Ew・TIP5P）のトポロジーを拒み、その数と原子番号を表示します。3 点の水を使うか、計算を `--mm-backend openmm` で実行してください。
* **必要なもの**: [AmberTools](installation.md) の tleap・antechamber・parmchk2 が `PATH` にあることが必要で、`--add-h` には PDBFixer も要ります。

---

## 関連ドキュメント

* [ML 領域と層の組み方](model-setup.md) — `mm-parm` が書く PDB で ML 領域と MM の層を決める
* [all](all.md) — 一気通貫のワークフロー。`--parm7` が無いと `mm-parm` を実行する
* [extract](extract.md) — トポロジーと対応した PDB から ML 領域を切り出す
* [define-layer](define-layer.md) — トポロジーと対応した PDB に ML・可動 MM・固定 MM の層を付ける
* [oniom-export](oniom-export.md) — Gaussian ONIOM・ORCA QM/MM の入力を書き出す。上の CMAP を含まないトポロジーが要る
* [トラブルシューティング](troubleshooting.md) — トポロジーと原子の順のエラー
