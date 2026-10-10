# `all`（一気通貫ワークフロー）

`all` サブコマンドは、活性部位のまわりの ML 領域の選択、全系の Amber トポロジーと 3 層の作成、最小エネルギー経路（MEP）の探索を 1 回の実行で行います。指定すれば、各反応段の遷移状態（TS）の最適化と、固有反応座標（IRC）・振動数・DFT の計算まで行います。

`--tsopt` を付けない場合は、TS 候補（各 MEP セグメントで最もエネルギーの高い点、HEI）までで終わります。ML 領域の計算バックエンドにはデフォルトの **UMA**（Meta が公開した学習済みの[機械学習原子間ポテンシャル（MLIP）](backends.md)）のほか、`-b/--backend` オプションで **ORB**、**MACE**、**AIMNet2**、[DFT](dft-backend.md)（`dft`）も選択可能です。酵素の残りの部分は Amber 力場で計算し、両者を ONIOM で合わせます。

---

## 主な用途

与える入力でモードが決まります。

* **R と P から経路とエネルギー図を作る（Endpoint モード）**: 反応順に並べた全系の 2 構造以上（反応物、中間体、生成物）を与えると、隣り合う構造の間の MEP を求め、エネルギー図を描きます。
* **反応物 1 つから経路を作る（Scan-list モード）**: 1 構造と、作る結合・切れる結合を `-s` で与えると、段階的スキャンで中間体を作り、それらを通る MEP を求めます。
* **TS 候補 1 つを確かめる（TS-only モード）**: 1 構造に `--tsopt` を付け、`-s` を付けずに与えると、TS を最適化し、そこから IRC をたどります。n_imag = 1 で、IRC が狙った R と P に着けば TS と確かめられます。

---

## 基本的な実行例

例は GPP C6-メチル基転移酵素 BezA（[Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)）の系で、スクリプト一式は [`examples/beza/`](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples/beza) にあります。`1.R.pdb`（反応物）・`2.IM.pdb`（中間体）・`3.P.pdb`（生成物）は、水素原子をすべて含む酵素全体の構造です。自分の構造にも[水素原子](getting-started.md)が要ります。例 1〜3 の流れと結果の確かめ方は [クイックスタート: `all`](quickstart-all.md)、[クイックスタート: `--scan-lists`](quickstart-scan.md)、[クイックスタート: TS-only モード](quickstart-tsopt.md) にあります。

### 1. MEP を求め、TS 最適化・熱化学・DFT まで計算する

`-c` で ML 領域の中心にする残基を、`-l` で非標準残基の電荷を指定します。`--refine-path` は結合が変わる所で経路を分けるので、化学的な段がそれぞれ 1 つのセグメントになります。

```bash
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --refine-path --tsopt --thermo --dft --out-dir ./result_mep
```

端末に各 TS の `[Imaginary modes] n=1 (...)` が出て、最後の `====== Pipeline summary ======` の下に `Scientific status: success` が出れば、指定した段はすべて終わっています。続けて [実行結果の判定](#実行結果の判定) のとおり端点を確かめてください。最適化した構造は `result_mep/segments/seg_NN/` にあります。

### 2. 反応物から段階的スキャンで経路を作る

段 1 で SAM のメチル炭素（CS1）を GPP の C7 に近づけ（1.50 Å）、SAM の SD から離します（3.30 Å）。段 2 で GPP の H11 を C7 から離し（2.90 Å）、Glu186 の OE2 に移します（1.00 Å）。

```bash
mlmm all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.50),("CS1 SAM 320","SD SAM 320",3.30)]' \
       '[("C7 GPP 321","H11 GPP 321",2.90),("OE2 GLU 186","H11 GPP 321",1.00)]' \
    --tsopt --thermo --out-dir ./result_scan
```

1 つのリテラルの中の目標は、同じ段で一緒に動きます。リテラルを並べると順に別の段として実行し、各段は前の段の終わりの構造から始まります。各段の終わりの構造が MEP 探索の入力になります。`-s` は 1 回だけ書き、その後にすべてのリテラルを並べてください。反応の分け方は [反応機構を調べるコツ](mechanism-tips.md) を参照してください。chain が空の PDB では、原子を残基名・残基番号・原子名の 3 つで順不同に指定します（`"CS1 SAM 320"`）。chain があるときは `A:SAM:320:CS1` と書きます。指定できる形はすべて [共通オプションと残基・原子の指定](cli-conventions.md) にあります。

### 3. TS 候補を確かめる（TS-only モード）

入力を 1 つにして `--tsopt` を付け、`-s` を付けないと、MEP 探索を省きます。

```bash
mlmm all -i TS_candidate.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft
```

最適化した R・TS・P は `result_all/segments/seg_01/` に書き出されます。

### 4. 途中のセグメントから後処理をやり直す

セグメント N から後処理をやり直すときは、元のコマンドと同じ入力・抽出・層・経路・計算バックエンドのオプションと同じ `--out-dir` を指定し、`--resume-segment N` を足してください。`--tsopt-max-cycles` などの後処理のオプションは変えられます。

```bash
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --refine-path --tsopt --thermo --dft \
    --resume-segment 1 --out-dir ./result_mep
```

N より前のセグメントはそのまま残し、セグメント N 以降の後処理と、要約・エネルギー図を書き直します。

---

## 処理の仕組みと計算仕様

```text
全系の構造 (PDB / mmCIF、または --ref-pdb 付きの XYZ)
  ├─ (-c のとき) ML 領域の選択: extract
  │   └─ ml_region.pdb
  ├─ 全系の Amber トポロジー: mm-parm（--parm7 のときは省く）
  │   └─ mm_parm/<入力名>.parm7
  ├─ B-factor による 3 層の割り当て: define-layer
  │   └─ layered/<入力名>_layered.pdb
  ├─ (1 構造 + -s のとき) 段階的スキャン: scan
  │   └─ 各段の終わりの構造 = 中間体
  ├─ MEP 探索: path-opt（デフォルト）または path-search（--refine-path）
  │   └─ mep_trj.xyz と energy_diagram_MEP.png
  └─ (--tsopt のとき) TS 最適化と IRC: tsopt → irc
      ├─ (--thermo のとき) 振動数と熱化学: freq
      └─ (--dft のとき) DFT 一点計算: dft
```

1. **入力と ML 領域の準備**: 元素の欄が空の PDB では、元素記号を補います。`-c` を指定すると、指定した残基のまわりの活性部位モデルを切り出します。最初の入力のモデルが ML 領域になり、`ml_region.pdb` に書き出されます。
2. **トポロジーと層の作成**: `mm-parm` が、最初の入力から AmberTools で全系の Amber トポロジーを作ります。力場は ff19SB、非標準残基は GAFF2 です。`--parm7` を指定するとこの段を省きます。続いて `define-layer` が、3 つの層を B-factor に書き込みます。ML（0）、可動 MM（10）、固定 MM（20）で、可動 MM は ML 領域から 8 Å 以内の残基です。以降の計算はすべて全系を ML/MM で行い、ML/MM の境界で切れる結合はリンク水素で埋めます。
3. **経路の作成**: まず入力構造を最適化します（`--preopt`）。`-s` を指定すると、段階的スキャンで中間体を作ります。続いて `path-opt` が、隣り合う構造の間の MEP を GSM（growing string method、デフォルト）か DMF（direct max flux）で求めます。`--refine-path` では、再帰的な `path-search` が経路を詰め、結合が変わる所で段に分けます。各段の HEI がその段の TS 候補です。
4. **TS の最適化と IRC の追跡**（`--tsopt`）: 各 HEI をデフォルトでは RS-P-RFO（restricted-step partitioned rational function optimization）で最適化し、ML と可動 MM の原子の最後の Hessian（PHVA、部分 Hessian 振動解析）から n_imag を求めます。TS から EulerPC（Euler 予測子–修正子法による積分）で IRC を両方向へたどり、両端を極小まで最適化します。これがそのセグメントの R と P になります。
5. **熱化学と DFT**: `--thermo` では R・TS・P で `freq` を実行して ML/MM のギブズエネルギーを求め、`--dft` では同じ構造の ML 領域を DFT で計算し、MM のエネルギーと合わせます。それぞれのエネルギー図も描きます。

---

## 実行結果の判定

TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます（n_imag = 1）。`all` が TS から IRC へ進むのは、TS 最適化が収束し、最後の Hessian を計算でき、n_imag ≥ 1 のときだけです。

| TS の結果 | `all` の次の動作 |
| --- | --- |
| 収束、n_imag = 1 | IRC を実行し、IRC の両端を最適化します。 |
| 収束、n_imag ≥ 2 | 警告を出し、MEP の方向に最もよく合う虚振動（合うものが無ければ最も低い虚振動）に沿って IRC を実行します。結果は `partial` です。 |
| 収束、n_imag = 0 | IRC の前で止まります。 |
| サイクル上限で停止 | 最後の Hessian を計算せず、IRC の前で止まります。 |
| `--stop-plateau` で停止（stalled） | 最後の Hessian を計算して n_imag を出し、IRC の前で止まります。 |
| `--skip-final-freq`、Hessian の失敗 | IRC の前で止まります。 |

TS 最適化の終わり方の一覧は [`tsopt` の「TS の判定」](tsopt.md#ts-の判定) にあります。

IRC が収束しなくても、端点の最適化で狙った R と P に着けば、その結果は使えます。

結果は次の 3 か所で確かめます。

* **端末**: 虚振動が 1 つの TS では、`[Imaginary modes] n=1 (...)` に虚振動数が出ます。`====== Pipeline summary ======` の下に `Execution status:` と `Scientific status:` が出ます。結果が `success` でないときは、`RESULT WARNING:` の行に理由が出ます。
* **`summary.log`**: ヘッダーに `Pipeline mode`（`MEP`、`Scan`、`TS-only`）と 2 つのステータスが出ます。[1] は MEP の概要、[2] は各セグメントの MEP 上の障壁 ΔE‡・反応エネルギー ΔE・結合の変化、[3] は各セグメントの後処理で、`TS imaginary freq:` の下に n_imag が出ます。[4] はエネルギー図の表、[5] は出力のツリーです。
* **`summary.json`**: `scientific_status` には、指定した段がすべて収束すると `success`、そうでなければ `partial` か `failed` が入り、理由は `scientific_status_reasons` に出ます。各 TS の n_imag は `post_segments[].tsopt.n_imaginary_modes` です。ステータスの欄の説明は [実行と要求段階の完了状況](json-output.md#実行の完了と指定した段の完了) にあります。
  * **障壁**: `--tsopt` のとき、各セグメントの障壁は、最適化した TS と R の ML/MM のエネルギーの差で、`post_segments[].mlip.barrier_kcal` に入ります。`--thermo` と `--dft` のときは、ほかの手法の障壁が同じ形で `gibbs_mlip`・`dft`・`gibbs_dft_mlip` に入ります。`segments[].barrier_kcal` は TS 最適化の前の MEP 上の障壁です。TS-only モードでは TS − R です。

端点が狙った R と P かは自分で確かめてください。`summary.log` の [2] の結合の変化と、`segments/seg_NN/reactant.*`・`product.*` の構造を、狙った R と P と比べます。n_imag が 1 でないときや、IRC の端点が狙いと違うときは {ref}`TS が取れないとき <ja-ts-search-fails>` を参照してください。

---

## 主な出力ファイル

`all` は `--out-dir` に次のファイルを書き出します。

```text
result_all/
├─ summary.log                  # 結果の要約（テキスト）
├─ summary.json                 # 機械可読な結果（常に出力。all に --out-json はありません）
├─ mep_trj.xyz                  # 全セグメントの MEP 軌跡
├─ mep_trj.pdb                  # 同じ軌跡の PDB
├─ mep_plot.png                 # MEP 軌跡に沿った ML/MM のエネルギー
├─ energy_diagram_MEP.png       # 全セグメントの MEP のエネルギー
├─ energy_diagram_*_all.png     # 全セグメントの R → TS → P の図（--tsopt、--thermo、--dft）
├─ irc_plot_all.png             # 全セグメントの IRC のエネルギー（--tsopt）
├─ ml_region.pdb                # ML 領域（--model-pdb で再利用できる）
├─ ml_region_without_linkH.xyz  # リンク水素を付ける前と後の ML 領域
├─ ml_region_with_linkH.xyz     #   （PDB 入力のときは .pdb も書く）
├─ mm_parm/                     # Amber トポロジー <入力名>.parm7 と .rst7（--parm7 で再利用できる）
├─ layered/                     # 3 層を B-factor に書き込んだ全系の構造
├─ segments/
│  └─ seg_NN/                   # 反応の 1 段: seg_01, seg_02, ...
│     ├─ reactant.*             # 最適化した R・TS・P（入力と同じ形式、--tsopt）
│     ├─ ts.*
│     ├─ product.*
│     ├─ energy_diagram_*.png   # この段の R → TS → P の図
│     ├─ ts/                    # TS 最適化。vib/imag_*_trj.xyz は虚振動のアニメーション
│     ├─ irc/                   # IRC の軌跡と irc_plot.png
│     ├─ endpoint_opt/          # 端点の最適化（--dump のとき、または端点が収束しなかったときに残る）
│     ├─ freq/{R,TS,P}/         # 振動数と熱化学（--thermo）
│     └─ dft/{R,TS,P}/          # DFT 一点計算（--dft）
└─ _work/                       # 途中のファイル（TS 候補の HEI を含む）
   ├─ pockets/                  # 抽出したモデル pocket_<入力名>.pdb（-c のとき）
   ├─ scan/                     # 段階的スキャン（-s のとき）
   └─ path_opt/                 # MEP 探索と hei_seg_NN.*（--refine-path のときは path_search/）
```

* **報告に使う構造**: `segments/seg_NN/reactant.*`・`ts.*`・`product.*` を使ってください。`seg_NN/` の下の各ディレクトリには、各段の計算のファイルが入っています。
* **`.cif` ファイル**: {ref}`mmCIF の入力 <ja-mmcif-input>` と、PDB の欄に入りきらない大きな PDB の入力では、元の識別子を保った `.cif` も書きます。
* **TS-only モード**: MEP 探索が無いので、MEP のファイルと `_work/path_opt/` はありません。R・TS・P は `segments/seg_01/` に入ります。

エネルギー図のファイル名は手法を表します。

| ファイル名 | 生成されるとき | 内容 |
| --- | --- | --- |
| `energy_diagram_MEP.png` | MEP 探索の完了時 | 全セグメントの MEP のエネルギー |
| `energy_diagram_MLIP.png` | `--tsopt` | R → TS → P、ML/MM のエネルギー |
| `energy_diagram_G_MLIP.png` | `--thermo` | R → TS → P、ML/MM のギブズエネルギー |
| `energy_diagram_DFT.png` | `--dft` | R → TS → P、ML/MM の構造での ML 領域の DFT のエネルギー |
| `energy_diagram_G_DFT_plus_MLIP.png` | `--dft` と `--thermo` | R → TS → P、ML(DFT)/MM のエネルギーに ML/MM の熱補正を足した値 |
| `energy_diagram_*_all.png` | `_all` の無い図と同じ | 全セグメントをまとめた図（出力ディレクトリの直下） |
| `irc_plot.png`（`seg_NN/irc/` の中）、`irc_plot_all.png` | `--tsopt` | 1 つのセグメントと全セグメントの IRC のエネルギー |

図のエネルギーは、最初の状態（反応物）を基準にした kcal/mol です。

---

## 主な CLI オプション

`all` は Amber のトポロジーと層を自分で作るので、{ref}`ML/MM の共通オプション <ja-mlmm-options>` は前の実行のトポロジーと `ml_region.pdb` を使い回すときだけ渡します。下の表は `all` に固有のものだけです。

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス（複数可） | （必須） | 反応順に並べた全系の 2 構造以上、または `-s` か `--tsopt` を付けた 1 構造（`.pdb`、`.cif`、`.mmcif`、`--ref-pdb` 付きの `.xyz`）。1 つの `-i` の後に並べるか、`-i` を繰り返す |
| `-c, --center` | 文字列 | `None` | ML 領域の中心にする残基（通常は基質と触媒残基）。残基名（`'SAM,GPP'`）、残基 ID（`'123,124'`、`'A:123,B:456'`）、PDB ファイル。省略すると入力の B-factor か `--model-pdb` から ML 領域を決める |
| `-l, --ligand-charge` | 文字列 | `None` | 非標準残基の電荷（例: `'SAM:1,GPP:-3'`）、またはその合計の数値。ML 領域の電荷とトポロジーの両方に使う |
| `-q, --charge` | 整数 | `None` | ML 領域の電荷。ML 領域から求める。明示すると求めた値より優先し、警告を出す |
| `-m, --multiplicity` | 整数 | `1` | ML 領域のスピン多重度（2S+1） |
| `-b, --backend` | 文字列 | `uma` | ML 領域の計算バックエンド（`uma`, `orb`, `mace`, `aimnet2`, `dft`） |
| `-r, --radius` | 浮動小数点数 | `2.6` | 中心原子からの抽出の半径（Å）。`0` では `-c` と `--selected-resn` の残基だけを残す |
| `--selected-resn` | 文字列 | `""` | 半径による拡張なしで入れる残基。ID（`'123'`、`'A:123A'`）、名前（`'SAM'`）、chain 付きの名前（`'A:SAM'`、`'A:SAM:123'`） |
| `--auto-mm-ff-set` | `ff19sb` / `ff14sb` | `ff19sb` | トポロジーの力場: ff19SB と OPC3 水、または ff14SB と TIP3P 水 |
| `--auto-mm-add-ter/--no-auto-mm-add-ter` | フラグ | `True` | トポロジーを作る前に、リガンド・水・イオンのブロックの前後と、つながっていないペプチドの間に TER を入れる |
| `--auto-mm-disulfide/--no-auto-mm-disulfide` | フラグ | `True` | SG–SG の距離で見つけたシステインを結合し、名前を CYX にする。無効のときは、名前がすでに CYX の残基だけを結合する |
| `--auto-mm-ligand-mult` | 文字列 | `None`（すべてのリガンドで 1） | トポロジーに使うリガンドのスピン多重度（例: `'GPP:2,SAM:1'`） |
| `--auto-mm-keep-temp` | フラグ | `False` | トポロジー作成の一時ディレクトリを残す |
| `-s, --scan-lists` | 文字列 | `None` | 1 構造の段階的スキャンの目標。1 つのリテラルが 1 段（例: `'[("A:SAM:320:CS1","A:GPP:321:C7",1.50)]'`。書き方は {ref}`スキャンリスト仕様 <ja-scan-list-spec>`） |
| `--tsopt/--no-tsopt` | フラグ | `False` | 各セグメントの TS を最適化し、IRC を実行 |
| `--thermo/--no-thermo` | フラグ | `False` | R・TS・P の振動数と ML/MM の熱化学（`--tsopt` が必要） |
| `--dft/--no-dft` | フラグ | `False` | R・TS・P の ML 領域の DFT 一点計算（`--tsopt` が必要） |
| `--refine-path/--no-refine-path` | フラグ | `False` | 隣り合う組ごとの `path-opt` の代わりに、再帰的な `path-search` を実行 |
| `--mep-mode` | `gsm` / `dmf` | `gsm` | MEP の手法: GSM または DMF |
| `--opt-mode` | `grad` / `hess` | `grad` | 単一構造の最適化とスキャンのオプティマイザ: `grad` = L-BFGS、`hess` = RFO。明示し、`--opt-mode-post` を指定しないときは、TS と端点にも使う |
| `--opt-mode-post` | `grad` / `hess` | `hess`（`--opt-mode` を明示したときはその値） | TS と IRC 後の端点のオプティマイザ: `grad` = TS は Dimer・端点は L-BFGS、`hess` = TS は RS-P-RFO・端点は RFO |
| `--preopt/--no-preopt` | フラグ | `True` | スキャンと MEP 探索の前に入力構造を最適化 |
| `--flatten/--no-flatten` | フラグ | `False` | TS 最適化の後に残った余分な虚振動を消す |
| `--stop-plateau/--no-stop-plateau` | フラグ | `False` | 収束の前にエネルギーが変わらなくなったら最適化を止める。収束ではなく stalled として報告する。MM のマイクロイテレーションはこれでは止めない |
| `--tsopt-max-cycles` | 整数 | `100000` | TS 最適化のサイクル上限 |
| `--resume-segment` | 整数 | `None` | `--out-dir` の MEP を使い、セグメント N から後処理をやり直す（例 4） |
| `--dry-run/--no-dry-run` | フラグ | `False` | 一時ディレクトリで入力の準備と事前の検査を行い、計画を表示して、計算の段は実行しない |
| `-o, --out-dir` | パス | `./result_all/` | 出力先ディレクトリ |

全オプションは `mlmm all --help-advanced` か [自動生成のオプションの一覧（英語のみ）](../reference/commands/all.md) を参照してください。

> **補足:** YAML（`--config`）では、上の表のオプションに無い設定もできます。コマンドラインで指定したオプションはファイルの値より優先されます。節とキーは [YAML 設定の一覧](yaml-reference.md) にあります。

---

## 使用上の注意点

* **`--dft` と `-b dft`**: 一緒には使えず、実行の始めにエラーで止まります。`-b dft` の計算の後に DFT の一点計算を足すときは、別のジョブで `mlmm sp -b dft` を実行してください。
* **`--dft` の費用**: 必要なメモリは ML 領域の大きさ・基底・汎関数・精度・ソフトウェアの構成で変わります。対象の計算ノードで代表的な構造を試し、最大メモリ使用量を見てください。ML 領域が大きいときは、MLIP の計算を先に終え、DFT の一点計算を別のジョブで実行してください。
* **TS-only モードの R と P**: IRC のエネルギーの高いほうの端を反応物と呼びます。R・P の名前、ファイル名、障壁、反応エネルギーはこのエネルギーの順に従うもので、化学的に分かった反応の向きではありません。P からの障壁は `barrier_kcal − delta_kcal` です。`summary.json` の `endpoint_assignment` にこの規則が記録され、`chemical_direction_known: false` となります。
* **TS-only モードの `summary.log`**: [1] は TS と IRC の概要で、[2] は最適化した TS と端点から求めます。
* **IRC の前で止まったとき**: TS のファイルは `segments/seg_NN/ts/` に残り、後のセグメントの後処理は行いません。
* **端点の最適化**: 片方の端点の最適化が収束しない場合、結果は `partial` になり、`segments/seg_NN/endpoint_opt/` を確認用に残します。端点の最適化がエラーで失敗した場合は、エラーを `segments/seg_NN/endpoint_opt/failure.json` に記録し、そのセグメントは振動数と DFT の段の前で止まります。TS と IRC の構造は残ります。
* **熱化学のファイル**: `all` は熱化学量を `thermoanalysis.yaml` から読むので、`--thermo` では `--no-dump` を指定してもこのファイルを残します。
* **抽出の半径**: `-r 0` では半径による拡張を無効にし、`-c` と `--selected-resn` で選んだ残基からモデルを組みます。構造上の安全策として、ジスルフィド結合の相手や隣の残基の主鎖が加わることはあります。半径 0 は内部で 0.001 Å として扱います。
* **`-c` を省いたとき**: 抽出を行わず、入力構造の全体を使います。ML 領域は入力の B-factor（`--detect-layer`、デフォルト有効）か `--model-pdb` から決めます。`--no-detect-layer` で `--model-pdb` も無いときは、エラーで止まります。1 構造のときは、このときも `-s` か `--tsopt` が必要です。
* **AmberTools**: `--parm7` を指定しないと、AmberTools が見つからないときにエラーで止まります。
* **入力の形式**: `all` は PDB と mmCIF を読みます。XYZ の入力には、同じ原子を持つ PDB を `--ref-pdb` で指定してください。1 回の実行のすべての構造は、同じ原子を同じ順に持つ必要があります。
* **電荷と多重度**: `-q` と `-m` は酵素全体ではなく ML 領域の値です。`-c` のときの ML 領域の電荷は、最初の入力から抽出したモデルの合計で、アミノ酸・イオン・水は組み込みの値、そのほかの残基は `-l` の値、`-l` に無い残基は 0 として数えます。`-c` が無いときは、B-factor か `--model-pdb` で決めた ML 領域で同じように合計します。`-q` は求めた値よりも優先し、警告を出します。値を求められないときは、YAML の `calc.model_charge` を使います。多重度は `-m`、無ければ YAML の `calc.model_mult`、それも無ければ 1 です。詳しくは [共通オプションと残基・原子の指定](cli-conventions.md) を参照してください。
* **固定原子と剛体運動**: 振動解析では、固定原子を動かさない並進と回転だけを射影で除きます。固定 MM の層があれば、除かれる運動はふつうありません。詳しくは [freq の「固定境界での剛体モード」](freq.md#固定境界での剛体モード) を参照してください。
* **別々に用意した構造**: 入力構造を別々に用意すると、反応座標の外の構造の違いも障壁に入ります。障壁を読む前に構造を比べてください。組成が同じ 2 つの機構を比べるときは、両方の経路で共通の原子の集合と順序を使います。
* **`--resume-segment`**: `--tsopt`・`--thermo`・`--dft` のどれかが必要で、`--dry-run` とは一緒に使えません。保存した入力・ML 領域・トポロジー・層を割り当てた構造・MEP がコマンドと合わないときは、エラーで止まります。

### 変異体と野生型の比較

1 つの経路の中では、すべての構造が同じ原子を同じ順に持ちます。変異体と野生型（WT）では残基が違い、原子数も変わることが多いので、両者の全エネルギーをそのまま差し引くことはできません。代わりに、それぞれの系の中で求めた障壁を比べます。

`ΔΔG‡ = (G_TS − G_R)_mutant − (G_TS − G_R)_WT`

* 2 つの系で、ML 領域に入れる残基と層の決め方をそろえ、狙った変異だけが違うようにします。半径で別々に選ぶと、境界の残基が片方の ML 領域にだけ入ることがあるので、2 つの `ml_region.pdb` を比べてください。
* プロトン化の決め方、電荷の決め方、力場、バックエンドとモデル、精度、拘束、熱化学の条件をそろえます。変異で ML 領域のプロトン化の状態や形式電荷が変わる場合は ML 領域の電荷も違うので、両方に同じ `-q` を当てはめないでください。

2 つの実行は、入力と出力先のほかは同じオプションにします。R が化学的な反応物になるよう、それぞれの系の R と P を与えます（Endpoint モード）。`G_TS − G_R` は `post_segments[].gibbs_mlip.barrier_kcal` です。

```bash
mlmm all -i wt_R.pdb wt_P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo --out-dir ./result_wt
mlmm all -i mutant_R.pdb mutant_P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo --out-dir ./result_mutant
```

---

## 関連ドキュメント

* [extract](extract.md) — ML 領域の選択
* [mm-parm](mm-parm.md) — 全系の Amber トポロジー
* [define-layer](define-layer.md) — B-factor による 3 層の割り当て
* [scan](scan.md) — 距離・角度・二面角の段階的スキャン
* [path-opt](path-opt.md) — 2 構造の間の MEP（GSM / DMF）
* [path-search](path-search.md) — 経路を段に分ける再帰的な MEP 探索
* [tsopt](tsopt.md) — 遷移状態（TS）の構造最適化
* [irc](irc.md) — TS からの IRC
* [freq](freq.md) — 振動解析と熱化学
* [dft](dft.md) — ML 領域の DFT 一点計算
* [MLIP の TS を DFT で確かめる](dft-backend.md) — `-b dft` と `--dft`
* [反応機構を調べるコツ](mechanism-tips.md) — 反応の分け方、TS の確かめ方、TS が取れないときの次の手
* [トラブルシューティング](troubleshooting.md) — 計算が失敗したとき
* [はじめに](getting-started.md) — 最短の実行と次に読むページ
