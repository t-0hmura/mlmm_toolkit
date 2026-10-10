# 反応機構を調べるコツ

mlmm-toolkit で得られる TS と経路は、仮説として立てた機構の候補です。このページでは、仮説を立てて計算を組み、TS を確かめ、候補の機構を比べて障壁を読むまでの手順をまとめます。

## 早見表

| 目的・症状 | 次の手 | 節 |
| --- | --- | --- |
| 仮説の機構から計算を始める | 作る結合・切れる結合・移る H を書き出し、入力モードを選ぶ | {ref}`仮説から始める <ja-mechanism-hypothesis>` |
| 協奏か段階かで計算を組む | `-s` のリテラルをまとめる・分ける | {ref}`反応の分け方を決める <ja-mechanism-split>` |
| TS が取れたかを見る | n_imag と IRC の端点を確かめる | {ref}`TS を確かめる <ja-mechanism-check-ts>` |
| n_imag が 2 以上で、余分なモードは反応部位の外 | `--flatten` を付ける。余分なモードが ML 領域と可動 MM 原子にまたがるなら `tsopt --no-microiter --flatten`、Hessian の範囲を絞っていたら広げ直す | {ref}`TS が取れないとき <ja-ts-search-fails>` |
| n_imag が 2 以上で、2 つのモードがどちらも反応の結合を動かす | 段に分けて試す | {ref}`TS が取れないとき <ja-ts-search-fails>`、{ref}`反応の分け方を決める <ja-mechanism-split>` |
| n_imag が 0、または候補が R・P の側へ滑り落ちた | `--refine-path` を付ける、別の初期構造から始める | {ref}`TS が取れないとき <ja-ts-search-fails>` |
| 候補の機構を比べる、または IRC の端点が狙った R・P と違う | 段の順を入れ替える、段階と協奏を比べる、ML 領域の範囲を見直す | {ref}`候補の機構を比べる <ja-mechanism-compare>` |
| 障壁を読む | その段の直前の極小から数える | {ref}`障壁を読む <ja-mechanism-barrier>` |
| 最適化が max cycles で止まる | オプティマイザを切り替える（`tsopt --opt-mode` / `all --opt-mode-post`）、ステップを小さくする、別の初期構造から始める | {ref}`トラブルシューティング：TS 最適化 <ja-troubleshooting-ts>` |

(ja-mechanism-hypothesis)=
## 仮説から始める

計算の前に、仮説の機構で作る結合・切れる結合・移る H を書き出してください。これがそのまま `-s` の座標になり、IRC の端点で確かめる点になります。書き出した原子は、すべて ML 領域に入れてください。

[入力モード](getting-started.md#3-反応経路を探す)は手元の構造で選びます。

- **R と P（あれば中間体も）がある**：[`all`](quickstart-all.md) に並べます。
- **R だけがある**：R から [scan](quickstart-scan.md) で経路を作ります。
- **TS 候補だけがある**：[TS-only モード](quickstart-tsopt.md)を使います。

(ja-mechanism-split)=
## 反応の分け方を決める

`-s` のリテラル（`-s` の後の角括弧のリスト）1 つが 1 つの段です。1 つのリテラルの中の座標は、同じ段で一緒に動きます（協奏）。リテラルを並べると、1 つずつ順に動かし（段階）、次の段に進むと前の段の拘束は外れます。各段の刻みの数は、座標ごとに変化量（目標値と初期値の差）を刻み幅の上限で割り、いちばん多くの刻みが要る座標で決まります。上限は、距離が `--scan-max-step-size`（0.20 Å。`scan` では `--max-step-size`）、角度が 5°、二面角が 10° です。

次の 2 つのコマンドは、同梱例の同じ 4 つの座標を動かします。作る C–C 結合、切れる C–S 結合、GPP から Glu186 へ移る H で、H の移動は 2 つの距離（C7–H11 と OE2–H11）で書きます。

4 つの座標を 1 つの段で協奏的に動かします。

```bash
mlmm all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -s '[("SAM,320,CS1","GPP,321,C7",1.50),("SAM,320,CS1","SAM,320,SD",3.30),("GPP,321,C7","GPP,321,H11",2.90),("GLU,186,OE2","GPP,321,H11",1.00)]' \
    --tsopt --thermo -o ./result_concerted
```

C–C 結合の形成と C–S 結合の切断を先に動かし、H の移動を 2 つ目の段にします。

```bash
mlmm all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -s '[("SAM,320,CS1","GPP,321,C7",1.50),("SAM,320,CS1","SAM,320,SD",3.30)]' \
       '[("GPP,321,C7","GPP,321,H11",2.90),("GLU,186,OE2","GPP,321,H11",1.00)]' \
    --tsopt --thermo -o ./result_stepwise
```

同梱の `1.R.pdb` は chain の欄が空なので、残基名・残基番号・原子名で指定しています（`"SAM,320,CS1"`）。chain がある PDB では `"A:SAM:320:CS1"` と書きます。指定できる形はすべて {ref}`スキャンリスト仕様 <ja-scan-list-spec>` にあります。

- **`-s` は 1 回だけ書く**：1 つの `-s` の後にリテラルを並べてください。この形は `all` でも `scan` でも通ります。
- **動く座標は同じ段に全部入れる**：作る結合だけでなく、切れる結合と移る H もその段に入れてください。2 つの座標だけを動かして、残りがついてくることを期待しないでください。
- **scan を使わない方法**：R と P（あれば中間体も）を用意できるなら、それらを `-i` に並べた MEP 探索にすることもできます。隣り合う 2 構造ごとに 1 つのセグメントになり、`--refine-path` を付けると、結合が変わる所で経路が段に分かれます。

(ja-mechanism-check-ts)=
## TS を確かめる

TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます。得られるのは TS の候補で、IRC の両端が狙った R と P に着くことで確かめます。IRC が収束しなくても、端点の最適化で狙った R と P に着けば、その結果は使えます。

- **終わり方**：収束した `tsopt` は `[microiter] Converged!`（マイクロイテレーションを使わないときは `[tsopt] Numerical optimization converged.`）を出し、続けて `[Imaginary modes] n=1 (...)` のように虚振動を出します。n_imag が 1 でないときは、`[tsopt] WARNING: Higher-order stationary point (n_imag=N). Try --flatten or all --refine-path.` か `[tsopt] No imaginary mode detected. Try all --refine-path.` も出ます。`--out-json` を付けると `result.json` も書きます。`all --tsopt` は TS ごとに同じ行を出し、指定した段がすべて収束すると、最後の `====== Pipeline summary ======` の下に `Scientific status: success` と出ます。`summary.log` の [3] には、セグメントごとの n_imag と IRC の出力がまとまります。判定の読み方は [tsopt の TS の判定](tsopt.md#ts-の判定) と [all の実行結果の判定](all.md#実行結果の判定) にあります。
- **max cycles とプラトー停止**：TS 最適化が max cycles に達して未収束のときは Hessian を計算しないので、n_imag は出ません。エネルギーが変わらなくなって止まったとき（`--stop-plateau`。{ref}`プラトー停止 <ja-optimizer-stalls-with-flat-energy--forces-just-above-threshold-mlip-force-noise-floor>`）は、必ず Hessian を計算して n_imag を出します。
- **端点**：終了コードが 0 であるだけでは、狙った TS が取れた証拠になりません。端点の最適化の後、IRC の両端で共有結合と移る H の付き先が狙った R と P と同じかを確かめてください。
- **診断用の IRC**：n_imag が 2 以上でも、最適化が数値的に収束し（プラトー停止は除く）、最後の PHVA（部分 Hessian 振動解析）が終わっていれば、`all` は IRC を続けます。反応の方向に合う虚振動のモードを使い、選べないときはいちばん低い虚振動のモードを使います。ログには `[all] WARNING: continuing diagnostic IRC from a numerically converged higher-order saddle …` と出ます。
- **モード**：`vib/imag_*_trj.xyz` をビューアで開き、各モードでどの原子が動くかを見てください。

(ja-ts-search-fails)=
## TS が取れないとき

### 虚振動が 2 つ以上残るとき

- **余分なモードが反応部位の外にある**（側鎖や水の回転など）：`--flatten` を付けてやり直してください。`tsopt`・`opt`・`all` で使えます。最適化の後に n_imag が 2 以上なら、余分なモードに沿って構造を少しずらして最適化をやり直し、デフォルトでは最大 {ref}`50 回 <ja-flatten-precedence-caveat>` まで繰り返します。
- **`--flatten` の後**：IRC の両端をもう一度確かめてください。n_imag が 1 になっても、別の反応の TS 候補に移ることがあります。
- **余分なモードが ML 領域と可動 MM 原子にまたがる**（基質と一緒に動く水など）：`tsopt` はデフォルトで {ref}`マイクロイテレーション <ja-microiteration>` を使うので、TS のステップで動くのは ML 原子とそれに結合した MM 原子だけです。`tsopt` に `--no-microiter --flatten` を付けてやり直し、TS のオプティマイザが動く原子のすべてを動かすようにしてください。`all` には `--microiter` が無いので、`all` の実行の `--parm7` と `--model-pdb` を `tsopt` に渡します。
- **Hessian の範囲を絞っていた**：`tsopt` に小さい [`--hessian-cutoff`](model-setup.md#hessian-の範囲を絞る) を付けていたら、広げるか外して、デフォルト（可動 MM 原子のすべて）に戻してください。
- **`--flatten` を使わない方法**：余分なモードのファイル `vib/imag_*_trj.xyz` と、B-factor に層を保った `vib/imag_*.pdb` には 20 のフレームがあり、6 番目と 16 番目が両方向にいちばん大きくずらした構造です。ビューアでこの 2 つをそれぞれ別の PDB に保存し、同じ電荷とスピンで、それぞれを TS-only モードか `tsopt` の初期構造にします。
  ```bash
  mlmm all -i <frame>.pdb --parm7 <run>/mm_parm/<name>.parm7 --model-pdb <run>/ml_region.pdb -l ... --tsopt
  ```
- **余分なモードも反応の結合を動かしている**：2 つの段が 1 つの候補に重なっているかもしれません。段に分けて試してください。

ほかの手（精度、座標の種類）は {ref}`tsopt：最適化後に n_imag が 1 でない場合 <ja-wrong-imaginary-mode-count>` にあります。

### 虚振動が無い・候補が滑り落ちたとき

- **HEI（最高エネルギーのイメージ）を見る**：HEI で作る結合と切れる結合がどちらも長い、または HEI が、結合が切れた後の構造になっていると、TS 最適化で反応のモードを失いやすくなります。HEI の結合長を確かめてください。
- **`--refine-path`**：`all` に付けると、1 回の `path-opt` の代わりに再帰的な `path-search` で経路を詰め、HEI を選び直します。デフォルトで無効です。先に粗い MEP を見てから使ってください。
- **経路の設定**：MEP の点の数（`--max-nodes`、デフォルト 20）、MEP の方法（`--mep-mode dmf`）、GSM のノードの置き方を変えると、HEI の位置が変わります。`--gsm-param energy` は、エネルギーの高い所にノードを多く置きます。
- **別の初期構造**：分割前の MEP の HEI（`path-opt` の `hei.xyz`）、段ごとの HEI（`path-search` の `hei_seg_NN.xyz`）、scan の最高点の近くのフレームから TS-only モードをやり直してください。`all` では、HEI のファイルは `_work/path_opt/`（`--refine-path` では `path_search/`）の下の `hei_seg_NN.{xyz,pdb}` です。

(ja-mechanism-compare)=
## 候補の機構を比べる

- **段の順**：段の順が、MEP 探索に渡す構造の並び（R → 1 段目の終わり → 2 段目の終わり …）になります。順を入れ替えると別の経路になるので、考えられる順をそれぞれ回し、エネルギー図で比べてください。
- **段階か協奏か**：段階で回し、拘束を外した最適化（`--scan-endopt`）の後も途中の結合の状態が残るかを見ます。残らずに元へ戻るなら、協奏のほうが合います。協奏で回しても、`--refine-path` で 2 つの段に分かれて同じ中間体を経るなら、段階のほうが合います。
- **ML 領域**：IRC の端点が狙いと違う、または反応に関わる残基・水・補因子が ML 領域の外にあるときは、{ref}`ML 領域 <ja-model-setup-larger>` を広げてください。周りの環境をもっと緩和させたいときは、ML 領域ではなく可動 MM を広げます。
- **DFT**：候補の TS を DFT/MM で確かめるときは、[MLIP の TS を DFT で確かめる](dft-backend.md) を参照してください。

(ja-mechanism-barrier)=
## 障壁を読む

- scan の最高点は TS の候補です。障壁は TS 最適化の後の値で読んでください。
- 障壁は、その段の直前の極小から数えます。`all` は各セグメントの障壁を、そのセグメントの E(TS) − E(R) として出します。TS-only モードでは、IRC の両端のうちエネルギーの高いほうが R です。候補の機構を比べるときは、どの候補でも基準の R をそろえてください。
- `--tsopt` の後、最適化した TS からの各セグメントの障壁は、`summary.log` の [3]（`Per-segment post-processing (TSOPT / Thermo / DFT)`）と `summary.json` の `post_segments[].mlip.barrier_kcal`（ML/MM のエネルギー）にあります。MEP 探索では、`summary.log` の [2]（`Segment-level MEP summary (ML/MM path)`）が TS 最適化の前の MEP の障壁です。[4]（`Energy diagrams (overview)`）はエネルギー図の表です。

## 使用上の注意点

- 診断用の IRC は、モードがどこへ向かうかの手がかりで、TS の確認ではありません。
- `--flatten` は余分な虚振動を消すだけで、無い反応のモードは作れません。
- 粗い経路では、`--refine-path` が反応を不要な段に分け、段ごとの MEP・TS 最適化・IRC・振動計算が増えることがあります。
- n_imag を数える閾値を変えても、構造は TS に近づきません。n_imag の数え方は [tsopt](tsopt.md) にあります。
- 座標はデフォルトのデカルト座標（`--coord-type cart`）のままにしてください。`tsopt` は `dlc` などの内部座標も、`all` は `cart` と `dlc` を受けますが、ML/MM のモデルでは内部座標を組むのに時間がかかります。

## 関連ドキュメント

- [はじめに](getting-started.md) — 入力モードの選び方
- [tsopt](tsopt.md) — TS の結果の読み方
- [all](all.md) — ワークフロー全体と出力
- [scan](scan.md) — 拘束付きスキャンと段のリテラル
- [path-opt](path-opt.md) — 2 構造間の MEP と HEI
- [path-search](path-search.md) — セグメントに分ける再帰的な MEP 探索
- [irc](irc.md) — 反応モードを R と P までたどる
- [freq](freq.md) — 虚振動の数え方
- [ML 領域と層の組み方](model-setup.md) — ML 領域と層を確かめる・広げる
- [define-layer](define-layer.md) — ML と MM の層の割り当て
- [ML/MM 計算機](mlmm-calc.md) — マイクロイテレーションと Hessian の範囲
- [MLIP の TS を DFT で確かめる](dft-backend.md) — 候補の TS を DFT/MM で確かめる
- [トラブルシューティング](troubleshooting.md) — エラーと収束の問題
- [共通オプションと残基・原子の指定](cli-conventions.md) — `-s` の原子の指定
