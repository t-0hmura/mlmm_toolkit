# MLIP の TS を DFT で確かめる

MLIP/MM で妥当な経路が見つかったら、その TS をそのまま DFT/MM での TS 構造最適化にもっていくことにも `mlmm-toolkit` は対応しています。TS 最適化 → IRC → 端点の最適化 → 振動数計算のワークフローを、GPU4PySCF を用いることで GPU で高速化された DFT 計算により実行可能です。DFT で計算するのは ML 領域だけで、周りのタンパク質は MM のままです。

主役は MLIP/MM による経路探索で、DFT/MM は MLIP/MM で見つけた TS 候補を確かめるための追加の機能です。

---

## 主な用途

- **TS を DFT/MM で詰める**：ML 領域を DFT にして、TS 最適化 → IRC → 端点の最適化 → 振動数計算を実行します（`-b dft`）。
- **MLIP/MM の計算に DFT のエネルギーを足す**：MLIP/MM で得た R・TS・P に DFT の一点計算を行います（`--dft`）。
- **1 構造の DFT 一点計算**：DFT/MM のエネルギーとポピュレーション解析（`mlmm dft`）、またはエネルギーと力（`mlmm sp -b dft`）を求めます。

## 処理の流れ

1. **MLIP/MM で探す**：経路を作り、条件を変えて試し、いちばん有望な TS 候補を選びます。
2. **DFT/MM で詰める**：その TS を入力にして、`-b dft` 付きの [TS-only モード](#基本的な実行例)を実行します。MM のトポロジーと ML 領域は最初の計算のものを使い回します。
3. **確かめる**：見る点は MLIP/MM のときと同じです。TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます（log には `[Imaginary modes] n=1`）。求めた段がすべて収束すると `====== Pipeline summary ======` の下に `Scientific status: success` と出ます。IRC が収束しなくても、端点の最適化で狙った R と P に着けば、その結果は使えます。

## 基本的な実行例

### 1. 小さな ML 領域で MLIP/MM の探索を行う

DFT で扱える大きさの ML 領域で MLIP/MM の探索を行います。出力ディレクトリ `result_all/` に、例 2 で使い回す TS・MM のトポロジー・ML 領域が入ります。

```bash
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -r 0 --selected-resn '44,63,186' --tsopt
```

`-r 0` で切り出しの半径を 0 Å にすると、距離で近くの残基を足すのを止め、`-c` と `--selected-resn` の残基から ML 領域を組みます。[`examples/beza/`](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples/beza) の同梱例の 44・63・186 番は、SAM のメチル炭素（CS1）にいちばん近い 3 残基です。同梱の PDB は chain の欄が空です。chain が空の PDB では、残基を名前か番号で指定してください。自分の系では、反応に関わる残基を選んでください。

### 2. TS を DFT/MM で詰める

例 1 の TS を 1 つだけ入力にすると、TS-only モードになります。`ts.pdb` は全系の構造で、`--parm7` と `--model-pdb` で例 1 のトポロジーと ML 領域を使い回すので、2 つの計算は同じ系を扱います。DFT 用の追加パッケージが要ります（[使用上の注意点](#使用上の注意点)）。

```bash
mlmm all -i result_all/segments/seg_01/ts.pdb \
    --parm7 result_all/mm_parm/1.R.parm7 --model-pdb result_all/ml_region.pdb \
    -l 'SAM:1,GPP:-3' --tsopt --thermo -b dft -o ./result_dft
```

`seg_01` は最初の反応セグメントです。セグメントが複数あるときは、詰めたい段のものを選んでください。

## ML 領域は 300 原子くらいまでに

`-b dft` では ML 領域全体を DFT で計算するので、リンク水素を含めて 300 原子くらいまでに収めます。

小さくする方法は [ML 領域を小さくする](model-setup.md#ml-領域を小さくする) を見てください。

## `-b dft` と `--dft` の違い

| オプション | DFT で計算するもの | 使う場面 |
|---|---|---|
| `-b dft` | 計算のすべて（MEP 探索、TS 最適化、IRC、端点の最適化、振動数）の ML 領域 | TS 候補を DFT/MM で詰めて確かめる |
| `--dft` | MLIP/MM で得た R・TS・P の一点計算（`all` だけ） | MLIP/MM の構造で DFT のエネルギーを得る |

`-b dft` は 11 のコマンドで使えます：`all`、`opt`、`tsopt`、`irc`、`freq`、`scan`、`scan2d`、`scan3d`、`path-opt`、`path-search`、`sp`。

## 主な出力ファイル

`-b dft` の出力の並びは、同じモードの MLIP/MM の計算と同じです（TS-only モードは[期待される出力](quickstart-tsopt.md#期待される出力)）。`--dft` を付けると、次のファイルが増えます。

| ファイル | 内容 |
|---|---|
| `segments/seg_NN/dft/{R,TS,P}/` | 各状態の DFT 一点計算の結果 |
| `segments/seg_NN/energy_diagram_DFT.png` | MLIP/MM の構造での DFT のエネルギー図 |
| `segments/seg_NN/energy_diagram_G_DFT_plus_MLIP.png` | DFT のエネルギーに MLIP/MM の熱補正を足した図（`--thermo` のとき） |
| `energy_diagram_DFT_all.png`、`energy_diagram_G_DFT_plus_MLIP_all.png` | 全セグメントをまとめた同じ図（出力ディレクトリの直下） |

## 主な CLI オプション

| オプション | 説明 | 既定値 |
|---|---|---|
| `-b, --backend dft` | ML 領域を DFT で計算します（GPU4PySCF。`--dft-engine cpu` で CPU の PySCF）。 | `uma` |
| `--parm7 FILE`、`--model-pdb FILE` | 前の計算の MM のトポロジーと ML 領域を使い回します。省くと、`all` が入力から作り直します。 | — |
| `--func-basis TEXT` | 汎関数と基底を `FUNCTIONAL/BASIS` の形で指定します。`-b dft` と `--dft` の両方に効きます。 | `wb97m-v/def2-svp` |
| `--embedcharge/--no-embedcharge`、`--embedcharge-cutoff FLOAT` | `-b dft` のとき、ML 領域から cutoff（Å）以内の MM の点電荷を DFT のハミルトニアンに入れます。 | `--no-embedcharge`、`12.0` |
| `--dft/--no-dft` | R・TS・P に DFT の一点計算を足します（`all` だけ）。 | `--no-dft` |

ほかの DFT のオプション（`--dft-engine`、`--dft-low-memory/--no-dft-low-memory`、`--dft-nprocs`、`--dft-memory`、SCF のチェックポイント）は [`all` のオプションの一覧（英語のみ）](../reference/commands/all.md)にあります。

> **補足:** YAML では同じ設定を `calc.dft` に書きます。`calc.dft.pyscf` では PySCF のオブジェクトの名前ごとに属性を渡せます（例：収束しにくい SCF に `mf: {level_shift: 0.2}`）。キーは [YAML 設定の一覧](yaml-reference.md)にあります。

## 使用上の注意点

- **DFT 用の追加パッケージ**：PyTorch の wheel が `cu130` か `cu132` なら `pip install "mlmm-toolkit[dft]"`、`cu126` なら `pip install "mlmm-toolkit[dft-cuda12]"` で入れてください（{ref}`詳細なインストール手順 <ja-step-by-step-installation>`の手順 7）。GPU が無いときは `--dft-engine cpu` を付けてください。
- **`--parm7` と `--model-pdb` を外さない**：外すと、DFT/MM の計算は TS の構造からトポロジーを作り直し、ML 領域をその B-factor の層から取ります。parm7 の名前は例 1 の最初の入力から付きます（ここでは `mm_parm/1.R.parm7`）。名前は `result_all/mm_parm/` で確かめてください。
- **電荷**：ML 領域を小さくすると、ふつう ML 領域の電荷も変わります。DFT/MM を流す前に、例 1 の端末の出力の `Total active site model charge` を確かめてください。
- **300 原子は目安**：code の上限ではありません。出力ディレクトリの `ml_region_with_linkH.xyz` の 1 行目が、リンク水素を含む ML 領域の原子数です。
- **組み合わせ**：`-b dft` と `--dft` は一緒に使えず、実行の始めにエラーで止まります。`-b dft` の計算の後に DFT の一点計算を足すときは、別のジョブで `mlmm sp -b dft` か `mlmm dft` を実行してください。`--dft` と `--thermo` には `--tsopt` が必要です。
- **図のファイル名とキーの名前**：`-b dft` でも、図のファイル名は `energy_diagram_MLIP.png` と `energy_diagram_G_MLIP.png`（`--thermo` のとき）、`summary.json` のブロックの名前は `mlip` と `gibbs_mlip` のままです。中身は DFT/MM の値で、図の題には DFT/MM と出ます。
- **メモリとスレッド**：`--dft-memory` の値は PySCF が使うホスト RAM で、GPU の VRAM ではありません。GPU のメモリが足りないときは ML 領域を小さくしてください。`--dft` でメモリが足りないときは、`--dft` を外して `mlmm dft` を別に実行してください。

## 関連ドキュメント

- [`dft`](dft.md)：DFT 一点計算とポピュレーション解析
- [`sp`](sp.md)：任意のバックエンドでの一点のエネルギーと力
- [クイックスタート: TS-only モード](quickstart-tsopt.md)：`all --tsopt` で TS 候補を確かめる
- [ML 領域と層の組み方](model-setup.md)：ML 領域を組む・小さくする・広げる
- {ref}`インストール <ja-step-by-step-installation>`：手順 7 で DFT 用の追加パッケージを入れる
- [MLIP バックエンド](backends.md)：バックエンドの選び方
- [トラブルシューティング](troubleshooting.md)：計算が失敗したとき
