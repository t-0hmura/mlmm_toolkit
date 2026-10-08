# クイックスタート: `mlmm all`

`mlmm all` は、反応物（R）と生成物（P）の 2 つの構造から、1 回の実行で反応経路を作ります。基質のまわりから ML 領域を切り出し、全系の Amber トポロジー（`mm-parm`）と 3 つの層（`define-layer`）を組んで、R と P の間の最小エネルギー経路（MEP）を ML/MM で探索します。`--tsopt --thermo --dft` を付けると、同じ実行のまま遷移状態（TS）の最適化・固有反応座標（IRC）の計算・振動数・DFT 一点計算まで進みます。反応座標の定義が手間なときに向いています。全系の中間体や生成物の構造を PyMOL や GaussView で自分で作って入力できます。反応座標を自分で決めずに MEP を探索するので、新しい反応機構の候補が見つかることもあります。

以下のコマンドは、ゲラニル二リン酸（GPP）の C6 位をメチル化する酵素 BezA の同梱例（[`examples/beza/`](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples/beza)）を使います。`1.R.pdb` が反応物、`2.IM.pdb` が中間体、`3.P.pdb` が生成物です。どれもすべての水素を含む酵素全体の構造です。同梱例は `git clone https://github.com/t-0hmura/mlmm_toolkit && cd mlmm_toolkit/examples/beza` で取得できます。自分の反応では、全系の構造に置き換えてください。

---

## 主な用途

* **全工程を初めて通す**: 同梱例で、すべての段を 1 回実行
* **R と P の間の MEP を作る**: 経路と、その最高エネルギーのイメージ（HEI、TS の候補）を取得
* **同じ実行で TS・IRC・振動数・DFT まで進める**: `--tsopt --thermo --dft` を付けて TS の候補を確認

## 最小コマンド

R と P を反応の順に渡し、ML 領域の中心にする残基（`-c`）とリガンドの電荷（`-l`）を指定します。

```bash
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --out-dir ./result_all
```

端末の最後のほうの `====== Pipeline summary ======` の下に `Scientific status: success` と出れば成功で、`summary.json` の `scientific_status` にも同じ値が入ります。

### （オプション）同一実行で後処理まで行う

`--tsopt` で[反応セグメント](glossary.md)（ここでは `seg_01`）ごとの TS 最適化と IRC を、`--thermo` で振動数と熱化学を、`--dft` で R・TS・P の ML 領域の DFT 一点計算を追加します。

```bash
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft --out-dir ./result_all
```

## 実行の前に

構造にはすべての水素原子が要り、R と P は同じ原子を同じ順に並べている必要があります。詳しくは [入力構造に関する重要事項](getting-started.md#入力構造に関する重要事項) を参照してください。トポロジーを組むために、さらに次の 2 点が必要です。

* **電荷と水素の数**: `-l` には、ファイルの中の水素の数に合う電荷を書いてください。同梱例の SAM は水素が 23 個なので `SAM:1` で、22 個なら `SAM:0` です。合わないと、`mm-parm` は `antechamber` を実行する前に電子数のエラーで止まります。
* **AmberTools**: `all` は AmberTools（`tleap`、`antechamber`、`parmchk2`）でトポロジーを組み、見つからないとエラーで止まります。既存のトポロジーを `--parm7` で渡すと、この段を省けます。

## 主な出力ファイル

最小コマンドは次のファイルを書き出します。

```text
result_all/
├── summary.log                  # 実行の要約
├── summary.json                 # 結果（scientific_status を含む）
├── mep_trj.pdb                  # 全セグメントの MEP
├── energy_diagram_MEP.png       # 全セグメントの MEP のエネルギープロファイル
├── ml_region.pdb                # ML 領域（--model-pdb で再利用できる）
├── mm_parm/                     # Amber トポロジー 1.R.parm7 と 1.R.rst7（--parm7 で再利用できる）
├── layered/                     # 3 つの層を B-factor に書いた 1.R_layered.pdb と 3.P_layered.pdb
└── _work/                       # 途中のファイル（TS 候補の HEI を含む。実行後も残る）
    └── path_opt/                # MEP 探索（MEP を再帰的に詰める --refine-path のときは path_search/）
        ├── hei_seg_01.{xyz,pdb} # セグメント 1 の最高エネルギーのイメージ
        └── summary.json         # MEP 探索の結果
```

最小コマンドは MEP 探索で終わるため、`segments/` は作られません。`--tsopt` を付けると、反応セグメントごとに `segments/seg_NN/` ができ、R/TS/P の構造（`reactant.pdb`・`ts.pdb`・`product.pdb`）、`ts/`、`irc/` が入ります。`--thermo` を付けると `freq/`、`--dft` を付けると `dft/` も加わります。mmCIF を入力したときは、`mep_trj.cif` などの `.cif` ファイルも書き出されます。

## 結果の確認

1. **完了状況**: `scientific_status` には、求めた段がすべて収束すると `success`、そうでなければ `partial` か `failed` が入り、[理由](json-output.md#実行と要求段階の完了状況)は `scientific_status_reasons` に出ます。`--tsopt` のとき、虚振動のモードができる結合と切れる結合を動かすかと、端点が狙った R と P かの 2 つは自分で確かめてください。
2. **TS の候補**: 最初のセグメントの HEI `_work/path_opt/hei_seg_01.pdb` を開きます。層を B-factor に持つ全系の PDB です。`--tsopt` のときは、最適化した TS の `segments/seg_01/ts.pdb` も開きます。
3. **エネルギープロファイル**: `energy_diagram_MEP.png` で、R と P の間にはっきりした障壁があるかを確かめます。
4. **TS（`--tsopt` のとき）**: TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます。このとき端末に `[microiter] Converged!` が出て、続いて `[Imaginary modes] n=1 (...)` の括弧の中に虚振動の波数（cm⁻¹）が出ます。`segments/seg_01/ts/vib/imag_*_trj.xyz` をビューアで開き、できる結合と切れる結合に沿って原子が動くかを確認してください。
5. **端点（`--tsopt` のとき）**: `segments/seg_01/irc/finished_irc_trj.xyz` と、最適化した端点の `segments/seg_01/reactant.pdb`・`product.pdb` を開き、狙った R と P かを確かめます。IRC が収束しなくても、端点の最適化で狙った R と P に着けば、その結果は使えます。

`all` が各段をどう判定するかは [実行結果の判定](all.md#実行結果の判定) を参照してください。

## 使用上の注意点

* **準備の結果の再利用**: `mm_parm/1.R.parm7` を `--parm7` で、`ml_region.pdb` を `--model-pdb` で次の実行や個別のコマンドに渡すと、トポロジーを組み直さずに同じ系を計算できます。
* **DFT と GPU のメモリ**: `--dft` には DFT 用の追加パッケージ（{ref}`詳細なインストール手順 <ja-step-by-step-installation>` の手順 7）が要ります。GPU メモリについては [MLIP の TS を DFT で確かめる](dft-backend.md#使用上の注意点) の使用上の注意点を参照してください。
* **`summary.json` の障壁**: `segments[].barrier_kcal` は TS 最適化の前の、MEP の上の障壁です。`--tsopt` を付けると、最適化した TS と端点から求めた ML/MM の障壁が `post_segments[].mlip.barrier_kcal` に入り、`--thermo` で `post_segments[].gibbs_mlip.barrier_kcal`、`--dft` で `post_segments[].dft.barrier_kcal` が加わります。`rate_limiting_step.barrier_kcal` は、すべてのセグメントにそろっている最も高いレベル（`DFT//MLIP/MM_Gibbs` > `DFT` > `MLIP_Gibbs` > `MLIP` > `MEP`）で比べた、最も高い障壁です。使ったレベルは `rate_limiting_step.method` に入ります。
* **実行時間**: 系の大きさ、ML 領域の大きさ、GPU、求めた段によって変わります。

## 次のステップ

- [クイックスタート: scan](quickstart-scan.md): 生成物の構造が無いとき、1 つの構造から始める
- [クイックスタート: TS-only モード](quickstart-tsopt.md): 手元の TS 候補を最適化して確かめる
- [ML 領域と層の組み方](model-setup.md): ML 領域を小さくする、残基が足りないときに広げる
- [反応機構を調べるコツ](mechanism-tips.md): 計算の計画と、TS が取れないときに試すこと
- [MLIP の TS を DFT で確かめる](dft-backend.md): TS を DFT/MM で詰めて確かめる
- [`all`](all.md): オプションの説明。`mlmm all --help-advanced` でも見られます
- [JSON 出力の一覧](json-output.md): `summary.json` の欄
- [トラブルシューティング](troubleshooting.md): エラーメッセージや症状から対処を探す
