# クイックスタート: `mlmm all --tsopt`（TS-only モード）

TS-only モードは、手元の遷移状態（TS）の候補 1 つを、最小エネルギー経路（MEP）の探索を省いて確かめます。`mlmm all --tsopt` は ML/MM の全系の上で TS を最適化し、固有反応座標（IRC）を両方向に追跡して、両端の構造（反応物 R と生成物 P）を最適化します。`--thermo` で振動解析と熱化学補正を、`--dft` で R・TS・P の ML 領域の DFT 一点計算を追加できます。

---

## 主な用途

* **スキャンや MEP の候補を TS に詰める**: スキャンの最高点や、MEP の最高エネルギーのイメージ（HEI）を TS まで最適化
* **別の方法で作った候補を確かめる**: 他のプログラムで得た構造や手で組んだ構造が、狙った R と P を結ぶ TS（n_imag = 1）かを確認
* **DFT に渡す前に MLIP の TS を確かめる**: 機械学習原子間ポテンシャル（MLIP）で得た TS を、[DFT/MM](dft-backend.md) で詰める前に確認

## 最小コマンド

全系の TS 候補を 1 つ渡し、`--tsopt` を付けます。同梱例には TS 候補が無いので、次のコマンドは [`all` のクイックスタート](quickstart-all.md) の実行で得た HEI を使います。自分の反応では、自分で用意した候補を渡してください。

```bash
mlmm all -i result_all/_work/path_opt/hei_seg_01.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo -o ./result_ts_only
```

端末の最後のほうの `====== Pipeline summary ======` の下に `Scientific status: success` と出れば成功で、虚振動が 1 つの TS では `[Imaginary modes] n=1 (...)` と出ます。状態は `summary.json` の `scientific_status` にも入ります。

前の実行とまったく同じ系で TS を計算するときは、そのトポロジーを `--parm7` で、ML 領域を `--model-pdb` で渡し、`-c` を省きます。

```bash
mlmm all -i result_all/_work/path_opt/hei_seg_01.pdb \
    --parm7 result_all/mm_parm/1.R.parm7 --model-pdb result_all/ml_region.pdb \
    -l 'SAM:1,GPP:-3' --tsopt --thermo -o ./result_ts_only
```

### （任意）DFT 一点計算を追加

`--dft` で R・TS・P の ML 領域の DFT 一点計算を追加し、`--func-basis` で汎関数と基底を指定します。

```bash
mlmm all -i result_all/_work/path_opt/hei_seg_01.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft --func-basis 'wb97m-v/def2-tzvpd' \
    -o ./result_ts_only
```

TS そのものを DFT/MM（`-b dft`）で最適化する方法と、DFT 用の追加パッケージのインストールや GPU メモリについては [MLIP の TS を DFT で確かめる](dft-backend.md) を参照してください。

## 実行の前に

* **入力**: 全系の TS 候補 1 つで、PDB か mmCIF、または同じ原子の PDB を `--ref-pdb` で添えた XYZ です。`-c` を付けると、指定した残基のまわりから ML 領域を切り出します。`-c` を省くと、`--model-pdb` か入力の B-factor の層から ML 領域を取ります。`all` が書き出した HEI の PDB は、この層を持っています。
* **電荷と多重度**: `-q` は ML 領域の電荷で、全系の電荷ではありません。`-q` を省くと、`all` は ML 領域の残基と `-l` から電荷を求めます（{ref}`電荷の指定 <ja-charge-specification>` を参照）。多重度は `-m` → YAML の `calc.model_mult` → 1 の順で決まります。
* **TS-only モードになる条件**: 入力が 1 つで、`--tsopt` を付け、`--scan-lists` を付けないとき。入力が 2 つ以上なら MEP 探索、入力 1 つに `--scan-lists` を付けるとスキャンになります。

## 期待される出力

成功すると、次のファイルが書き出されます。

```text
result_ts_only/
├── summary.log                     # 実行の要約
├── summary.json                    # 結果（scientific_status を含む）
├── ml_region.pdb                   # ML 領域（--model-pdb で再利用できる）
├── mm_parm/                        # 候補から組んだ Amber トポロジー（--parm7 のときは無い）
├── layered/                        # hei_seg_01_layered.pdb：3 つの層を B-factor に書いた候補（-c のとき）
└── segments/
    └── seg_01/
        ├── reactant.pdb            # R/TS/P の構造（入力と同じ形式）
        ├── ts.pdb
        ├── product.pdb
        ├── energy_diagram_MLIP.png # R–TS–P の ML/MM エネルギー図（--thermo で energy_diagram_G_MLIP.png も）
        ├── ts/
        │   ├── final_geometry.{xyz,pdb}
        │   └── vib/imag_*_trj.xyz  # 虚振動モードごとのアニメーション
        ├── irc/
        │   └── {forward,backward,finished}_irc_trj.xyz
        ├── freq/{R,TS,P}/          # --thermo のとき
        │   ├── frequencies_cm-1.txt
        │   └── thermoanalysis.yaml
        └── dft/{R,TS,P}/           # --dft のとき
            └── result.yaml
```

## 結果の確認

1. **完了状況**: `scientific_status` には、求めた段がすべて収束すると `success`、そうでなければ `partial` か `failed` が入り、[理由](json-output.md#実行と要求段階の完了状況)は `scientific_status_reasons` に出ます。虚振動のモードができる結合と切れる結合を動かすかと、端点が狙った R と P かの 2 つは自分で確かめてください。
2. **TS のモード**: TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます。このとき端末に `[microiter] Converged!` が出て、続いて `[Imaginary modes] n=1 ([-593.1])` のように虚振動の波数（cm⁻¹）が出ます。`segments/seg_01/ts/vib/imag_*_trj.xyz` をビューアで開き、できる結合と切れる結合に沿って原子が動くかを確認してください。
3. **端点**: `segments/seg_01/irc/finished_irc_trj.xyz` と R/TS/P の構造（`reactant.pdb`・`ts.pdb`・`product.pdb`）を開き、`segments[0].bond_changes` を読みます。端点は狙った R と P のはずです。IRC が収束しなくても、端点の最適化で狙った R と P に着けば、その結果は使えます。
4. **端点の振動数**: `--thermo` のとき、`segments/seg_01/freq/{R,TS,P}/frequencies_cm-1.txt` に符号つきの全振動数が出ます。R と P には虚振動（−5.00 cm⁻¹ より小さい値）が無いはずです。
5. **エネルギー**: `post_segments[0].mlip.barrier_kcal` が ΔE‡（TS − R）、`.delta_kcal` が ΔE（P − R）で、単位は kcal/mol です。どちらも最適化した TS と端点の ML/MM エネルギーから求めます。MEP が無いので、`segments[0].barrier_kcal` と `.delta_kcal` にも同じ値が入ります。`--thermo` のときは `post_segments[0].gibbs_mlip.barrier_kcal` と `.delta_kcal` が ΔG‡ と ΔG、`--dft` のときは `post_segments[0].dft.barrier_kcal` と `.delta_kcal` が DFT の値です。

`all` が各段をどう判定するかは [実行結果の判定](all.md#実行結果の判定) を参照してください。

| 結果 | 次に試すこと |
|---|---|
| n_imag = 0 | MEP の HEI やスキャンの最高点など、よりよい候補から始めます。TS-only モードには手がかりにする経路がありません。 |
| n_imag ≥ 2 | 各虚振動のモードを見ます。`--flatten` を付けて最適化し直すか、`all --thresh-post gau_tight`（既定の [`baker`](tsopt.md#処理の仕組みと計算仕様) より厳しい）または `tsopt --thresh gau_tight` で収束を厳しくします。 |
| `bond_changes` が空、または端点が狙いと違う | TS のモードと IRC を確認します。経路が別の極小点どうしを結んでいる可能性があります。 |
| R や P に虚振動が残る | 端点の構造とモードを確認し、`--thresh-post gau_tight` で端点の最適化を厳しくします。 |

n_imag が 1 でないときや、IRC の端点が狙いと違うときは、[反応機構を調べるコツ](mechanism-tips.md) の {ref}`TS を確かめる <ja-mechanism-check-ts>` と {ref}`TS が取れないとき <ja-ts-search-fails>` を参照してください。

## 使用上の注意点

* **IRC に進む条件**: `all` が IRC に進むのは、TS 最適化が収束し、最後の Hessian の計算が終わり、n_imag が 1 以上のときだけです。n_imag が 2 以上のときの IRC は、虚振動の 1 つに沿った診断用の計算で、構造が一次の鞍点になるわけではありません。
* **余分な虚振動**: `--flatten` は、余分な虚振動のモードに沿って構造をずらして最適化し直すことを、最大 50 回くり返します。{ref}`--flatten を使うとき <ja-flatten-precedence-caveat>` を参照してください。
* **Hessian の計算法**: 既定の `--hessian-calc-mode FiniteDifference` のまま使ってください。`--hessian-calc-mode Analytical` は、自分の系の代表的な構造で既定と速度・メモリ・結果を比べてから指定します。
* **R と P の名付け**: MEP が無いと反応の向きが分からないため、TS-only モードはエネルギーが高いほうの IRC 端点を R、低いほうを P と呼び、この決まりを `summary.json` の `endpoint_assignment` に記録します。この名前は化学的な反応の向きではありません。P からの障壁は `barrier_kcal − delta_kcal` です。
* **`tsopt` と `freq` を単独で使うとき**: `--opt-mode`、`--max-cycles`、`--no-microiter`、`--hessian-cutoff` などの Hessian のオプションを変えるときは、この実行の `--parm7` と `--model-pdb` を渡して [`tsopt`](tsopt.md) を単独で実行してください。final geometry は `final_geometry.{xyz,pdb}` に、虚振動のアニメーションは `vib/` に出ます。その構造の全振動数と熱化学は、`final_geometry.pdb` に同じ `--parm7`・`--model-pdb`・`-q`・`-m` を付けて [`freq`](freq.md) を実行すると得られます。`all` の全オプションは `mlmm all --help-advanced` で確認できます。

## 次のステップ

- [MLIP の TS を DFT で確かめる](dft-backend.md): TS を DFT/MM で詰めて確かめる
- [反応機構を調べるコツ](mechanism-tips.md): TS の確かめ方と、TS が取れないときに試すこと
- [`tsopt`](tsopt.md)・[`irc`](irc.md)・[`freq`](freq.md): 各段を単独で実行する
- [クイックスタート: `mlmm all`](quickstart-all.md): R と P から MEP を作る
- [クイックスタート: scan](quickstart-scan.md): 1 つの構造から経路を作る
- [`all`](all.md)・[`dft`](dft.md): 全オプションのリファレンス
- [トラブルシューティング](troubleshooting.md): エラーメッセージや症状から対処を探す
