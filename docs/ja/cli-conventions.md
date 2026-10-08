# 共通オプションと残基・原子の指定

全コマンドに共通するフラグ、ML/MM の共通オプション、残基と原子の指定、電荷と多重度、終了コード、設定の優先順位をまとめたページです。

## ブール値オプション

段階や動作のオン・オフは、対になったフラグで指定します。

| 形 | 例 |
|---|---|
| 有効にする | `--tsopt` |
| 無効にする | `--no-tsopt` |

```bash
--tsopt --thermo --no-dft
```

## 段階的ヘルプ

```bash
mlmm <subcmd> --help               # 主要オプションのみ
mlmm <subcmd> --help-advanced      # 全オプション
```

どのサブコマンドも両方を受け付けます。

実行の前に設定を確かめるときは、`--show-config` で `--config` に渡した YAML とその最上位のキーを表示し、`--dry-run` で計算をせずにオプションと入力を検査できます。

```bash
mlmm opt -i input.pdb --parm7 real.parm7 -q -1 --config my_settings.yaml --show-config --dry-run
```

確認が通ると、最後に `[Dry run] --dry-run completed. Input command is valid.` と出ます。

(ja-verbosity-levels)=

## ログ詳細度 (verbosity)

`-v/--verbose LEVEL` は 0〜3 の整数 (**デフォルト 2**) で、各コマンドのコンソール出力量を決めます。コマンドごとのオプションなので、サブコマンドと一緒に指定します (例: `mlmm opt -v 1 ...`)。4 段階は全コマンド共通で、各コマンドページはそのコマンド固有の出力だけを説明します。

| レベル | 表示内容 |
|---|---|
| `-v 0` | 無出力。成功は終了コードと出力ファイルで確認します。 |
| `-v 1` | マイルストーンのみ: バージョン、入力要約、主要設定、出力先、dry-run / 最終ステータス。バナー・`[command]`・`[mode]`・設定の一覧は出ません。 |
| `-v 2` | デフォルト。バナー、`[command]`、`[mode]`、ステージ進捗、主要なオプティマイザのサイクル表、終了ステータス、Hessian 1 行要約、熱化学 / DFT 要約、経過時間を追加します。 |
| `-v 3` | デバッグ: 実際に使う設定の全体、バックエンドの DEBUG、オプティマイザ・内部座標の詳細、`[HessianTiming]`、`[HessianVRAM]`。 |

レベルで変わるのは表示だけで、終了コードは変わりません。成否は {ref}`終了コード <ja-exit-codes>` で判断してください。

(ja-mlmm-options)=

## ML/MM の共通オプション

`-q` は系全体ではなく ML 領域の電荷です。ML/MM の個別のコマンドには、全系のトポロジー（`--parm7`）と ML 領域が要ります。`all` はどちらも自分で作ります。

各コマンドは、次のうち最初に与えられたものから ML 領域を決めます：`--model-pdb`、次に `--model-indices`、次に `--detect-layer` での入力 PDB の ML 原子（B-factor 0）。ML 領域を明示したときも、`--detect-layer` は B-factor から Movable-MM と Frozen-MM の層を読みます。

| オプション | 意味 | デフォルト |
|---|---|---|
| `--parm7` | 酵素複合体全体（全系）の Amber parm7 トポロジー。 | 必須（`all` は自分で作る） |
| `--model-pdb` | ML 領域だけの PDB（リンク水素を含まない）。原子の名前と順序は全系の PDB・parm7 と同じにします。与えると、これで ML 領域が決まります。 | なし |
| `--model-indices` | ML 領域の原子番号（1 始まり、カンマ区切り。`1-5` のような範囲も可）。`--model-pdb` が無いときに使います。 | なし |
| `--detect-layer/--no-detect-layer` | B-factor 0 / 10 / 20 を ML・Movable-MM・Frozen-MM の層として読みます（{ref}`MM の層 <ja-mm-layers>`）。ML 領域を明示したときは MM の層だけを読みます。 | `--detect-layer` |
| `--ref-pdb` | XYZ 入力のときに、原子の順序と残基の情報を与える PDB。 | なし |
| `--movable-cutoff` | ML 領域からの距離（Å）。この内側の MM 原子は動き、外側は固定されます。 | なし（層は B-factor か `--freeze-atoms` から） |
| `--mm-backend` | MM の計算エンジン：`hessian_ff` か `openmm`。MM の Hessian はデフォルトで有限差分です。 | `hessian_ff` |
| `--link-atom-method` | リンク原子の置き方：`scaled`（g-factor）か `fixed`（1.09 / 1.01 Å）。 | `scaled` |
| `--cmap/--no-cmap` | parm7 に CMAP 項があれば、real system と model system の両方の MM 計算に残します。 | `--cmap` |

`mlmm all` は、`--parm7` を省くとトポロジーを作り、`-c` で抽出したモデルから ML 領域を決めます。

```bash
mlmm path-search -i R.pdb P.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1
```

MLIP バックエンドの選び方、精度、`--uma-workers` は [MLIP バックエンド](backends.md) にあります。

## 残基セレクタ

`extract` と `all` の `-c/--center` は、モデルの中心にする残基を指定します。下の表は、範囲の狭い形から順に並べています。

| 形 | 例 | 選ばれる残基 |
|---|---|---|
| chain＋残基名＋番号（推奨） | `-c 'A:TYR:44'` / `-c 'A:TYR:44,A:SAM:123'` | 1 項目につき 1 残基だけを確実に選べます。 |
| chain＋残基名 | `-c 'A:SAM'` | chain A の SAM をすべて選びます。複数あるときは警告をログに出します。 |
| chain＋番号 | `-c 'A:123'` / `-c 'A:123,B:456'` / `-c 'A:123A'` | chain A の残基 123 を選びます。末尾の英字は挿入コードです。 |
| 残基名だけ | `-c 'SAM,GPP'` / `-c 'LIG'` | どの chain でも、同じ名前の残基をすべて選びます。複数あるときは警告をログに出します。 |
| 番号だけ | `-c '123,456'` / `-c '123A'` | すべての chain から同じ番号の残基を選びます。 |
| 構造ファイル | `-c substrate.pdb` / `-c substrate.cif` | 別の PDB / mmCIF ファイルの座標と一致する残基を選びます。 |

mmCIF の長い chain ID や 9999 を超える残基番号も、同じ形で指定できます。chain ID は大文字と小文字を区別します。残基名は区別しません。同梱の例の PDB のように chain 欄が空の PDB では、残基名か番号の形だけを使えます（`-c 'SAM,GPP,MG'`、`--selected-resn '44,63,186'`）。

```bash
mlmm extract -i complex.cif -c 'LONG_CHAIN:SAM' -o pocket.pdb        # chain LONG_CHAIN の SAM をすべて
mlmm extract -i complex.cif -c 'LONG_CHAIN:SAM:10001' -o pocket.pdb  # SAM を 1 つだけ
mlmm extract -i complex.cif -c 'LONG_CHAIN:10001' -o pocket.pdb      # chain＋番号
```

(ja-selected-resn-takes-ids)=
### `--selected-resn` も同じ残基セレクタを使う

`extract` と `all` の `--selected-resn` は、指定した残基をモデルに必ず含めます。指定の形は上の表と同じです。たとえば `A:TYR:44` は 1 残基、`A:SAM` は chain A の SAM すべて、`A:123A` は挿入コード付きの 1 残基を含めます。chain を付けない `TYR` は一致する残基をすべて含め、複数あるときは警告を出します。

(ja-charge-specification)=
## 電荷の指定

`-q/--charge` は ML 領域の電荷です。PDB/mmCIF 入力では、`--ligand-charge/-l` を使うと**非標準残基（基質・補因子など）の電荷だけ**を指定すれば、標準アミノ酸やイオンの電荷と合算して ML 領域の電荷が自動計算されます。組み込みの表にあるイオンは表の電荷を使います（`MG` は +2）。`-l` に同じ値で書いても構わず、違う値は警告を出して無視します。

```bash
-l 'SAM:1,GPP:-3'        # 残基ごとの指定（推奨）
-l 'LIG:-2'              # 1 残基の指定
-l -3                    # 整数 1 つ = リガンドの総電荷
-q 0                     # ML 領域の電荷を明示
```

**電荷の決まり方**（上ほど優先）:

1. 明示的な `-q/--charge`
2. `--ligand-charge/-l` があるとき：ML 領域の標準残基・イオン・指定したリガンドの電荷の合計
3. `--config` 内の `calc.model_charge`
4. どれでも決まらなければ実行を中断

導出した電荷は、端末の `Total active site model charge` の行に出ます。`extract` の後に読む行は [モデルを確かめる](model-setup.md#モデルを確かめる) にあります。

```{tip}
非標準の残基（基質、補因子、特殊なリガンド）には必ず `--ligand-charge/-l` を指定し、電荷が正しく伝播するようにしてください。
```

## スピン多重度

```bash
-m 1    # 一重項 (singlet)（デフォルト）
-m 2    # 二重項 (doublet)
-m 3    # 三重項 (triplet)
```

`-m/--multiplicity` は ML 領域の多重度です。指定しないときは `--config` 内の `calc.model_mult`、それも無ければ 1 を使います。`all` と各サブコマンドで同じ値を指定してください。`mm-parm` は、パラメータ化に使う残基ごとの多重度を別の `--ligand-mult` で受け取ります。

## 原子セレクタ

原子セレクタは `--scan-lists` と `opt` の `--distance-restraint` で 1 つの原子を指します。`--freeze-atoms` は 1 始まりの原子番号だけを受け付けます（{ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>`）。

```bash
--scan-lists '[(1, 5, 2.0)]'                                          # 1 始まりの整数インデックス
--scan-lists '[("SAM,320,CS1", "GPP,321,C7", 1.60)]'                  # 残基名、残基番号、原子名
--scan-lists '[("A:SAM:320:CS1", "A:GPP:321:C7", 1.60)]'              # chain ID 付き
```

3 項目のセレクタは、残基名・残基番号・原子名を任意の順序で並べ、空白・カンマ・コロン・スラッシュ・バッククォート・バックスラッシュのどれでも区切れます（`"SAM,320,CS1"`、`"SAM 320 CS1"`、`"320,SAM,CS1"` は同じ原子）。3 項目に chain は入りません。chain を指定するときは、4 項目の形 `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` をこの順で書き、挿入コードは残基番号に続けます（`A:SAM:12B:C1`）。

(ja-scan-list-spec)=
### スキャンリスト仕様

`-s/--scan-lists`（`scan` / `scan2d` / `scan3d` / `all` で使用）は 1 個以上のインライン Python リテラルを受け付けます。スタンドアロンの `scan` / `scan2d` / `scan3d` はこれに加えて YAML/JSON スペックファイルパスも受け付けます。複数ステージの実行にはファイルが、短い指定にはインラインリテラルが適しています。

**YAML/JSON スペックファイル**（ルートはマッピング。キーは `scan` では `stages`、`scan2d` / `scan3d` では `pairs`）:

```yaml
one_based: true            # 任意。未指定時はコマンドの --one-based/--zero-based（デフォルトは 1 始まり）に従う
stages:                    # scan 用
  - [[1, 5, 1.35]]
  - [[1, 5, 2.20], [2, 8, 1.80]]
```

```yaml
one_based: true
pairs:                     # scan2d（要素は 2 つちょうど） / scan3d（要素は 3 つちょうど）
  - [1, 5, 1.30, 3.10]
  - [2, 8, 1.20, 3.20]
```

YAML スペックでは、`scan` の 1 ステージに距離の目標値と、距離・角度・二面角の範囲を混ぜられます。`scan2d` / `scan3d` の各軸は `(i,j,low,high)`、`(i,j,k,low,high)`、`(i,j,k,l,low,high)` のいずれかです。インデックスは整数、3 項目のセレクタ、または位置固定の `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` で指定できます。距離の単位は Å、角度と二面角は度です。

**インラインリテラル**: シェルがカッコや空白を解釈しないように、リスト全体を**シングルクォート**で囲み、中の PDB セレクタはダブルクォートで書いてください。

```bash
-s '[(atom1, atom2, target_Å), ...]'             # scan: 3 要素タプル
-s '[(atom1, atom2, low_Å, high_Å), ...]'        # 距離の range
-s '[(atom1, atom2, atom3, low_deg, high_deg)]'  # 角度の range
-s '[("SAM,320,CS1","GPP,321,C7",1.60)]'         # クォートしたセレクタ
-s "[(\"SAM,320,CS1\",\"GPP,321,C7\",1.60)]"       # 非推奨: 外側をダブルクォートで囲むと内側のエスケープが必要
```

`scan` では **1 リテラル = 1 ステージ**です。複数ステージを実行するには、**1 つの `--scan-lists` フラグの後に複数リテラル**を並べます。`scan2d` / `scan3d` ではリテラルは **1 つだけ** を受け付けます（複数ステージは非対応）。結合を 1 つ作ってからプロトンを動かす 2 ステージの `scan` は次のとおりです。

```bash
mlmm scan -i r.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s '[("SAM,320,CS1","GPP,321,C7",1.60)]' '[("GPP,321,H11","GLU,186,OE2",0.90)]'
```

| コマンド | 受け付けるスキャン仕様 |
| --- | --- |
| `scan` | インラインの距離の目標値 `(i,j,target)`、または入力構造から両方向にスキャンする範囲 `(i,j,low,high)`・`(i,j,k,low,high)`・`(i,j,k,l,low,high)`（[双方向スキャン](scan.md#双方向スキャン4-tuple)）。YAML/JSON にも対応 |
| `all --scan-lists` | インラインの目標値のみ：距離 `(i,j,target)`、角度 `(i,j,k,deg)`、二面角 `(i,j,k,l,deg)`（範囲と YAML/JSON は不可） |
| `scan2d` | 距離・角度・二面角の軸を 2 本含む 1 リテラル/ファイル |
| `scan3d` | 距離・角度・二面角の軸を 3 本含む 1 リテラル/ファイル |

したがって 4 要素のタプルは、`scan` では距離の範囲、`all` では角度の目標値です。

## 入力ファイル要件

- **PDB** — 水素原子を含み（`reduce`、`pdb2pqr`、Open Babel、`mlmm mm-parm --add-h` で追加）、列 77–78 に元素記号が必要です（欠けている場合は `mlmm add-elem-info`）。複数の PDB は同じ原子を同じ順序で持つ必要があります。
- **mmCIF** — 下の {ref}`mmCIF と大きな構造 <ja-mmcif-input>` に挙げた 14 のコマンドが受け付けます。
- **XYZ** — 計算コマンド（`all` を含む）は、XYZ を `--ref-pdb` と一緒に受け付けます。`--ref-pdb` が原子の順序と残基の情報を与えます。
- **Amber parm7（`--parm7`）** — 全系の力場パラメータです。原子は、どの座標入力とも同じものが同じ順序で並んでいる必要があります。

(ja-mmcif-input)=
### mmCIF と大きな構造

`all`、`extract`、`define-layer`、`sp`、`opt`、`tsopt`、`freq`、`irc`、`dft`、`scan`、`scan2d`、`scan3d`、`path-opt`、`path-search` は `.cif` と `.mmcif` を受け付け、`--ref-pdb` も同じです。単独の `mm-parm` は PDB だけを読みます。`all` は mmCIF の入力を変換してからパラメータを作ります。2 文字以上の chain ID、4 桁を超える残基番号、5 桁を超える原子番号を持つ構造や、残基が 10,000 以上の構造には mmCIF を使ってください。

`mlmm` は最初の座標モデルを読み、altLoc（別位置の配座）は残基ごとに平均占有率の最も高いものを 1 つ残します。計算の間は一時的な chain ID と残基番号を使い、出力の CIF で元の chain ID・残基番号・挿入コードに戻します。大きな PDB や標準の欄に収まらない PDB（残基が 10,000 以上、原子が 99,999 以上、hybrid-36 の番号、桁あふれした数字の欄など）も同じように扱います。

```bash
mlmm all -i reactant.cif product.cif -c 'enzyme_A:SAM:10001,enzyme_A:GPP:10002' \
    -l 'SAM:1,GPP:-3' --tsopt --thermo -o result
mlmm tsopt -i hei.xyz --ref-pdb full_system.mmcif --parm7 full_system.parm7 -q -2 -o result_tsopt
```

2 つ目のコマンドは、座標を XYZ ファイルから、トポロジーを mmCIF ファイルから取ります。残基と原子のセレクタには、元の chain ID と残基番号を使います。

| 場面 | 形 | 例 |
|---|---|---|
| `extract` / `all -c` の番号指定 | `CHAIN:RESSEQ[ICODE]` | `enzyme_A:10001B` |
| `extract` / `all -c` の残基名＋番号指定 | `CHAIN:RESNAME:RESSEQ[ICODE]` | `enzyme_A:SAM:10001B` |
| スキャンの原子 | `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` | `enzyme_A:SAM:10001B:CS1` |

出力の横に書かれる `.cif` については [出力ディレクトリのレイアウト](output-layout.md) を参照してください。

(ja-trajectory-one-frame)=
### 軌跡から 1 フレームを取り出す

`_trj.xyz` はふつうの複数フレームの XYZ（各フレームは原子数の行・コメントの行・原子の行）なので、k 番目（1 から数える）のフレームは次のように取り出せます。

```bash
N=$(head -1 scan_trj.xyz); k=12
sed -n "$(( (k-1)*(N+2)+1 )),$(( k*(N+2) ))p" scan_trj.xyz > frame_12.xyz
```

PDB のトポロジーを使って続けるときは、元の PDB を次の計算コマンドの `--ref-pdb` に渡してください。座標はフレームのものを使います。

(ja-exit-codes)=
## 終了コード

| コード | 意味 |
|---|---|
| `0` | 成功、または利用できる一部の結果あり |
| `1` | 未収束、利用できる結果なし、実行中の例外、出力失敗 |
| `2` | 入力・CLI 引数・設定の指定ミス |
| `130` | ユーザー中断（SIGINT） |

JSON の有無で終了コードは変わりません。終了コード `0` には `success` と `partial` の両方が入るので、`scientific_status` で見分けます。`all` と `path-search` はこれを `summary.log` に書き（`all` は端末にも `Scientific status:` として出します）、ほかのコマンドは `--out-json` を付けたときに `result.json` に記録します（[実行と要求段階の完了状況](json-output.md#実行と要求段階の完了状況)）。`--out-json` を付けないときは、そのコマンドのページにある端末の行で見分けてください（例：`irc`（固有反応座標）では [IRC の成否の判定](irc.md#irc-の成否の判定)）。

(ja-opt-mode-semantics)=

## `--opt-mode`（サブコマンド依存）

`--opt-mode` は最適化法を選びます。L-BFGS（記憶制限 BFGS）と RFO（有理関数最適化）は極小を探し、Dimer、RS-P-RFO（制限ステップ分割 RFO）、RS-I-RFO（制限ステップイメージ RFO）、TRIM（信頼領域イメージ最小化）は TS（遷移状態）を探します。

| サブコマンド | `grad` エイリアス | `hess` エイリアス | デフォルト |
|------------|------------------|------------------|-----------|
| `opt` | L-BFGS (`lbfgs`) | RFO (`rfo`) | `grad` (L-BFGS) |
| `tsopt` | Dimer (`dimer`) | RS-P-RFO (`rsprfo`) | `hess` (RS-P-RFO) |
| `path-opt`（端点の事前最適化） | L-BFGS | RFO | `grad` |
| `path-search`（1 構造ずつの最適化） | L-BFGS | RFO | `grad` |
| `scan` / `scan2d` / `scan3d`（格子の各点の緩和） | L-BFGS | RFO | `grad` |
| `all`（前処理の最適化・スキャン・MEP 探索、`--opt-mode`） | L-BFGS | RFO | `grad` |
| `all`（TS 最適化、`--opt-mode-post`） | Dimer | RS-P-RFO | `hess` |
| `all`（IRC 後の端点の最適化、`--opt-mode-post`） | L-BFGS | RFO | `hess` |

`--opt-mode-post` を省いて `--opt-mode` を明示すると、`all` はその値を TS 最適化と端点の最適化にも使います。

同じ `--opt-mode` の値でも、サブコマンドによって**選ばれる最適化法が異なり**、デフォルトも異なります。レシピをコピーする前に表を確認してください。アルゴリズム名を受け付けるのは `opt`（`lbfgs` / `rfo`）と `tsopt`（`dimer` / `rsirfo` / `trim` / `rsprfo`）だけで、ほかのサブコマンドは `grad` / `hess` だけを受け付けます。したがって `tsopt` の `--opt-mode grad` は L-BFGS 最小化ではなく **Dimer TS 探索**で、この Dimer は Hessian を周期的に計算してダイマーの方向を更新します。曖昧さを避けるには、`tsopt` では `--opt-mode dimer` か `rsirfo`、`opt` では `--opt-mode lbfgs` か `rfo` と書いてください。

`--microiter/--no-microiter`（デフォルトで有効）は、ML 領域の Hessian を使う 1 ステップと、MM 原子の L-BFGS 緩和を交互に行います。`opt` では RFO、`tsopt` では RS-P-RFO・RS-I-RFO・TRIM のときに使います。Dimer と静電埋め込み（`--embedcharge`）では使いません。

## CLI ↔ YAML 名称の不一致

一部の CLI フラグは YAML の対応キーと微妙に名前が異なり、`all` でラップされたときにリネームされるものもあります。主なフラグと YAML キーの対応は {ref}`YAML 設定の一覧の主要な CLI→YAML マッピング <ja-common-cli-to-yaml-mapping>` にあります。特によく聞かれるケースを以下に示します:

(ja-pressure-vs-pressure-atm)=
- **`--pressure` (CLI) と `pressure_atm` (YAML)** — `freq` のフラグは `--pressure FLOAT`、`all` では `--freq-pressure` です。YAML キーは `thermo.pressure_atm` です。どちらも値は **atm** 単位です（デフォルト 1.0）。

- **`--step-size` (CLI) と `step_length` (YAML)** — `irc` のフラグは `--step-size FLOAT`（bohr、デフォルト 0.10）、`all` では `--irc-step-size` です。YAML キーは `irc.step_length` です。

(ja-engine-vs-dft-engine)=
- **`--dft-engine` (CLI) と `engine` (YAML)** — `--dft-engine`（別名 `--engine`）は `gpu`（GPU4PySCF、デフォルト）か `cpu`（PySCF）を選びます。YAML キーは、`dft` と `all` の DFT 段では `dft.engine`、`--backend dft` で実行するほかの計算コマンドでは `calc.dft.engine` です。

```bash
mlmm irc -i ts.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 --step-size 0.05
mlmm all -i r.pdb p.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt --irc-step-size 0.05
```

## YAML 設定

```bash
mlmm all -i r.pdb p.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --config my_settings.yaml --out-dir result/
```

(ja-configuration-precedence)=

```
組み込みデフォルト  <  --config (YAML)  <  CLI オプション
```

各オプションの組み込みデフォルトは、`mlmm <subcmd> --help-advanced` と [コマンドの一覧（英語のみ）](../reference/commands/index.md) の `[default: …]` で確かめられます。YAML を上書きするのは*明示的に指定した* CLI の値だけで、CLI のデフォルトのままのオプションは YAML の値を隠しません。この順序は `--config` を受け付けるすべてのコマンドに共通です。全設定は [YAML 設定の一覧](yaml-reference.md) を参照してください。

## 出力ディレクトリ

`-o/--out-dir ./my_results/` で計算コマンドの出力先を指定します。コマンドごとのデフォルトは [出力ディレクトリのレイアウト](output-layout.md) にあります。準備のコマンドはファイルパスを取ります：`extract` は `-o/--output` に 1 つ以上のパス、`mm-parm` は `-o/--out-prefix`、`define-layer` は `-o/--output` に 1 つのパスです。

## 使用上の注意点

* **残基名と番号には chain を付けてください**: `TYR:44` は「chain `TYR` の残基 44」と読まれ、「見つからない」というエラーで止まります。`A:TYR:44` と書いてください。
* **1 つのリストには 1 種類の残基セレクタだけを使ってください**: 表の形は、残基名だけ（`SAM`）、chain＋残基名（番号の有無を問わない。`A:SAM`・`A:TYR:44`）、番号（chain の有無を問わない。`A:123`・`123`）の 3 種類に分かれ、種類の違う形は 1 つのリストに混ぜられません。`A:TYR:44,A:SAM` は通りますが、`A:SAM,SAM`・`A:44,A:SAM`・`SAM,TYR:44` はエラーで止まります。
* **chain 欄が空の PDB での原子セレクタ**: `'SER:11:HG'` や `'SER 11 HG'` のように 3 項目で書きます。`_` は空の chain の意味にならないので、`'_:SER:11:HG'` はどの原子にも一致せず、エラーで止まります。
* **`--model-indices` と電荷の導出**: `--model-indices` では `--ligand-charge` から電荷を導出できません。`-q` を指定するか、`--model-pdb` か B-factor の層（`--detect-layer`）で ML 領域を決めてください。
* **`--movable-cutoff` は `--detect-layer` を無効にします**: `--movable-cutoff` を指定すると、MM の層は距離で決まり、入力の B-factor の層は読みません。
* **B-factor の層には ML 原子と MM 原子の両方が要ります**: B-factor を層として読むのは、ML（0）の原子と MM（10 か 20）の原子がそれぞれ 1 個以上あり、原子の 80% 以上がこれらの値（±1.0）を持つときだけです。全原子が 0 の PDB は層の指定になりません。
* **インラインのスキャンリテラル**: 1 つのインラインリテラルと、1 回の実行のすべてのリテラルは、距離の目標値か範囲のどちらか一方だけを持ちます。両方を使うときは、YAML/JSON スペックの `stages:` に並べてください。
* **mmCIF と大きな構造の制限**: 扱える残基は 619,938 個までです（計算の間に使う内部の PDB の、1 文字の chain ID 62 種 × 残基番号 9,999）。`fix-altloc` と `add-elem-info` は PDB だけを読みます。mmCIF では、読み込みのときに altLoc を選び、元素記号を `_atom_site.type_symbol` から取ります。`_atom_site.type_symbol` の無い行、有限でない座標、反応の状態の間で原子数が違う入力はエラーで止まります。座標は PDB の固定幅の欄に収まる必要があるので、原点から遠い構造は原点の近くに移してください。

## 関連ドキュメント

- [インストール](installation.md) — セットアップと依存関係
- [はじめに](getting-started.md) — 最短の実行と次に読むページ
- [出力ディレクトリのレイアウト](output-layout.md) — ファイル名とデフォルトの出力ディレクトリ
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法
- [YAML 設定の一覧](yaml-reference.md) — 全設定オプション
- [MLIP バックエンド](backends.md) — バックエンドの選び方、精度、workers
- [ML/MM 計算機](mlmm-calc.md) — ML/MM のエネルギー・力・Hessian の計算方法
