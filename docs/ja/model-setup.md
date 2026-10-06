# ML 領域と層の組み方

mlmm-toolkit は酵素の全体を計算します。基質の周りの ML 領域を MLIP で、残りのタンパク質を Amber の力場（MM）で扱います。計算の重さは 3 つで決まります。ML 領域の原子の数、動く MM 原子（Movable-MM）の数、Hessian に入る原子の範囲です。このページでは、ML 領域と層の組み方、計算を軽くするための削り方、残基が足りないときの広げ方、原子の固定と距離の拘束をまとめます。

## 早見表

| 目的 | すること | 節 |
| --- | --- | --- |
| 既定のモデルを組む | `all` の `-c` に基質・補因子・金属を並べる | [既定のモデルを組む](#既定のモデルを組む) |
| 重い計算の前に原子数と電荷を見る | `define-layer` の Layer Summary か `all` のログを読む | [モデルを確かめる](#モデルを確かめる) |
| 計算を軽くする | ML 領域を小さくする、動く MM の殻を薄くする、ML だけの Hessian にする | [モデルを削る](#モデルを削る) |
| 反応に関わる残基が入っていない | `-r` を広げる、`-c` か `--selected-resn` に足す | {ref}`モデルを広げる <ja-model-setup-larger>` |
| 自分で組んだモデルを使う | `--parm7` と `--model-pdb` を渡す | [自分で組んだモデルを使う](#自分で組んだモデルを使う) |
| `model.pdb` を手で作る | チェックリストに沿って確かめる | {ref}`信頼できる model.pdb の作り方 <ja-model-pdb-selection>` |
| 原子を固定する・距離を拘束する | Frozen-MM の層、`--freeze-atoms`、`--distance-restraint` | {ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>` |

## 既定のモデルを組む

### ML 領域

`all` の `-c` に基質を書き、補因子と金属も同じリストに入れてください。

```bash
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3'
```

これは [`examples/beza/`](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples/beza) の同梱例です。`all` は `extract` で ML 領域を切り出します。`-c` の原子から `-r`（既定 2.6 Å）以内に原子が 1 つでもある残基を入れ、水も入れます。規則の全体は [extract](extract.md#処理の仕組みと計算仕様) にあります。`all` は 1 つ目の入力の ML 領域を出力ディレクトリの `ml_region.pdb` に書くので、次の実行ではこれを `--model-pdb` に渡すと同じ ML 領域を使えます。

(ja-mm-layers)=
### MM の層

ML 領域の外の原子は、2 つの MM の層に分かれます。各原子の層は PDB の B-factor の欄に書かれます。

| 層 | B-factor | 役割 |
| --- | --- | --- |
| **ML** | 0 | 反応する領域。エネルギー・力・Hessian を MLIP で計算する |
| **Movable-MM** | 10 | 最適化で動く MM 原子 |
| **Frozen-MM** | 20 | 座標を固定する MM 原子。MM のエネルギーには加わる |

[`define-layer`](define-layer.md) は、ML 原子から `--movable-cutoff`（既定 8.0 Å）以内の MM 原子を Movable-MM に、それより外を Frozen-MM にします。`-c` があると、`all` は入力ごとに 8.0 Å で `define-layer` を実行し、層を付けた PDB を出力ディレクトリの `layered/` に書きます。

## モデルを確かめる

重い計算の前に、原子数と電荷を確かめてください。

- **層**：`define-layer` は Layer Summary に `Layer 1 (ML, B=0):`、`Layer 2 (Movable MM, B=10):`、`Layer 3 (Frozen MM, B=20):`、`Total atoms:` の行を出します。`all` は入力ごとに `[all] define-layer [i]: … (ML=…, MovableMM=…, FrozenMM=…)` を出します。
- **リンク水素**：`all` は `[all] ML structure with link H (N + M; …)` を出します。N が ML 領域の原子数、M がリンク水素の数です。
- **電荷**：`extract` と、`-c` を付けた `all` は、ML 領域の電荷を `Total active site model charge` の行に出します。多重度も確かめてください。
- **目で見る**：層を付けた PDB をビューアで開いて B-factor で色分けし、反応に関わる残基が ML 領域に入っているかを確かめてください。

## モデルを削る

3 つの重さには、それぞれ別の設定があります。

| 決めるもの | 重くなる理由 | `all` で | 個別のコマンドで | 既定 |
| --- | --- | --- | --- | --- |
| ML 領域 | 毎ステップ MLIP（`-b dft` では DFT）で計算する | `-r`、`--exclude-backbone`、`--no-include-h2o`、`--selected-resn` | `--model-pdb` | `-r 2.6` Å |
| Movable-MM | 最適化で動かす原子になる | オプションなし（`-c` があると 8.0 Å） | `define-layer --movable-cutoff`、または `opt`・`tsopt`・`freq`・`scan`・`scan2d`・`scan3d`・`path-opt`・`path-search`・`sp` の `--movable-cutoff` | 8.0 Å |
| Hessian の範囲 | ML 領域と Movable-MM にまたがる密な行列になる | オプションなし | `opt`・`tsopt`・`freq`・`sp` の `--hessian-cutoff` | ML と Movable-MM の全部 |

### ML 領域を小さくする

同梱例（`examples/beza/1.R.pdb`、`-c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3'`）では、次のオプションで下の大きさの ML 領域になります。原子数は `extract --add-linkh` で数えた値で、キャップ水素を含みます。

| ML 領域の作り方 | オプション | 原子数 |
| --- | --- | --- |
| 既定 | — | 632 |
| 主鎖の原子を MM 側に回す | `--exclude-backbone` | 397 |
| さらに水も MM 側に回す | `--exclude-backbone --no-include-h2o` | 367 |
| 残基を自分で選ぶ | `-r 0 --selected-resn '44,63,186'` | 129 |

- **`--exclude-backbone`** と **`--no-include-h2o`**：主鎖の原子と水を ML 領域から外します。外した原子は MM として系に残ります。
- **`-r 0 --selected-resn`**：距離で残基を足さず、`-c` の残基と、指定した残基だけで ML 領域を作ります。同梱例の 44・63・186 番は、SAM のメチル炭素（CS1）にいちばん近い 3 残基です。自分の系では、反応に関わる残基を選んでください。

`all` でも同じオプションを使えます。例は [MLIP の TS を DFT で確かめる](dft-backend.md#基本的な実行例) の例 1 に、DFT/MM で計算するときの大きさは [ML 領域は 300 原子くらいまでに](dft-backend.md#ml-領域は-300-原子くらいまでに) にあります。ML 領域を小さくすると電荷が変わることが多いので、毎回確かめてください。

### 動く MM の殻を薄くする

`--movable-cutoff` を小さくすると、Movable-MM が減って Frozen-MM が増え、最適化で動かす原子が減ります。`define-layer --movable-cutoff` で層を付け直した PDB を個別のコマンドに渡すか、上の表のコマンドに `--movable-cutoff` を付けてください。

### Hessian の範囲を絞る

`--hessian-cutoff` を付けると、ML 領域からその距離以内の Movable-MM 原子だけを Hessian に入れます。`freq` の振動解析と、`tsopt` の最後の振動解析では、解析する原子（`--active-dof-mode`。既定の `partial` は ML 領域と Movable-MM の全部）を Hessian が覆っている必要があります。Hessian のほうが狭いと計算の前にエラーで止まるので、解析する原子も一緒に絞ってください。ML だけの Hessian にするときは `--hessian-cutoff 0.0 --active-dof-mode ml-only` を渡します。`all` には `--hessian-cutoff` がありません。Hessian の MM の部分の計算方法は [MM Hessian](mlmm-calc.md#mm-hessian) にあります。

動く MM 原子が多い系では、{ref}`マイクロイテレーション <ja-microiteration>`（`tsopt` と `opt --opt-mode hess` で既定で有効）が ML のステップの合間に MM 原子を力場だけで緩和するので、MLIP を呼ぶ回数が減ります。

(ja-model-setup-larger)=
## モデルを広げる

反応に関わる残基・水・補因子が ML 領域の外にあるときは、ML 領域を広げてください。タンパク質の残りはもう MM に入っているので、ML 領域に足すのは、結合・プロトン化・電荷が変わる原子と、その共有結合の相手です。

- **`-r` を広げる**（既定 2.6 Å）。
- **`--selected-resn`**：残基を足します。足した残基からは距離での探索を始めません。
- **`--radius-het2het`**（既定 0 で無効）：中心の側も相手の側も C・H 以外の原子だけで測る、2 つ目の距離を足します。全体の半径を広げずに、近くの N・O の相手を拾えます。
- **残基を `-c` に足す**：残基を丸ごと残したいときに使います。`-c` のアミノ酸からも距離での探索が始まり、`-r` が 0 より大きければペプチド結合でつながった隣の残基が加わるので、`--exclude-backbone` なしでは全部の原子が残ります。

周りの環境を動けるようにしたいときは、`define-layer --movable-cutoff 10.0` などで、ML 領域ではなく Movable-MM を広げてください。Hessian を ML 領域だけに絞っていたら、`--hessian-cutoff` を外して既定に戻します（[反応機構を調べるコツ](mechanism-tips.md)）。

半径は、ML 領域の大きさに対して結果が収束したかを確かめるためのパラメータです。ML 領域を広げると計算は重くなり、精度が上がるとは限らないので、化学的に妥当な何通りかの ML 領域で、エネルギー・力・障壁を比べてください。反応物・中間体・生成物では同じ ML 領域を使い、`ml_region.pdb` を `--model-pdb` で渡します。

## 自分で組んだモデルを使う

個別のコマンドには `--parm7` と ML 領域の指定が要り、`all` はこの 2 つを自分で作ります。モデルを手で組むときは、3 つのコマンドを実行します。`mm-parm` がトポロジーと、同じ原子を同じ順に並べた PDB を作り、`extract` がその PDB から ML 領域を切り出し、`define-layer` が同じ PDB に層を書きます。

```bash
mlmm mm-parm -i input.pdb -l 'LIG:0' --out-prefix system
mlmm extract -i system.pdb -c LIG -l 'LIG:0' -o model.pdb
mlmm define-layer -i system.pdb --model-pdb model.pdb -o system_layered.pdb
```

個別のコマンドか、`-c` を省いた `all` に、`system_layered.pdb` と `--parm7 system.parm7`、`--model-pdb model.pdb` を渡してください。どのコマンドも、ML 領域を `--model-pdb`、`--model-indices`、B-factor の層の順に探します（{ref}`ML/MM の共通オプション <ja-mlmm-options>`）。電荷は、ML 領域の電荷を `-q` で渡すか、残基名ごとに `-l` で渡します（PDB/mmCIF の入力）。領域を手で切ったときや、プロトン化を手で変えたときは、残基名から電荷を求められないので `-q` で渡してください。

(ja-model-pdb-selection)=
### 信頼できる `model.pdb` の作り方

`model.pdb` は、独立に組み直したクラスターではなく、全系の PDB / `parm7` から原子を選んだ**原子選択ファイル**です。原子名、残基名と番号、chain ID、全系での原子の順を変えないでください。原子番号の振り直し・並べ替え、リンク水素の手での追加、別に水素を付けたモデルの書き出しはしません。

- 反応中心、共有結合した補因子や相手、経路の途中でプロトン化や結合が変わる原子は、すべて ML 領域に入れます。
- タンパク質の主鎖の断片を残すときは、主鎖の両端がどちらも Cα（`CA`）で終わる範囲を選びます。境界の原子価は ML/MM のリンク原子の処理に任せます。
- 側鎖・リガンド・補因子の境界は、なるべく脂肪族の **C–C 単結合**（`CA–CB` か、反応中心からさらに外側）に置きます。ペプチドの C–N、極性の C–N・C–O、芳香環や共役の結合、S–S 結合、金属への配位結合は切りません。結合の相手を入れるか、境界を動かします。
- 反応物・中間体・生成物では、全系の原子とその順を同じにし、同じ `model.pdb` の選択を使います。状態ごとに別々にモデルを作ると、原子の対応と障壁の比較が成り立ちません。
- 本番の計算の前に、すべての境界を目で見て、ML 領域の電荷と多重度を確かめます。`define-layer` は層を付けるだけで、化学的にまずい境界は直しません。

PyMOL で `model.pdb` を保存するときは、書き出しの画面で **Original atom order** に印を付けてください。

(ja-freeze-atoms-and-restraints)=
## 原子の固定と距離の拘束

固定は、原子の力を 0 にしてその場から動かさないことです。拘束は、2 原子の距離を調和ポテンシャルで目標の値へ引くことです。ML/MM では Frozen-MM の層がタンパク質の外側をもう固定しているので、`--freeze-atoms` には、ほかに止めたい原子だけを書いてください。

### 原子を固定する 3 つの方法

- **Frozen-MM の層**（B-factor 20）：自動で固定します（{ref}`MM の層 <ja-mm-layers>`）。
- **`--freeze-atoms 'i,j,k'`**：全系の 1 始まりの原子番号です。`all`・`opt`・`tsopt`・`irc`・`freq`・`scan`・`scan2d`・`scan3d`・`path-opt`・`path-search`・`sp` で使えます。
- **YAML の `geom.freeze_atoms`**（`--config` で渡す）：リストが長いときや、ほかの設定と一緒に残したいときに使います。

```yaml
geom:
  freeze_atoms: [12, 15, 28, 29, 42]   # 1 始まり
```

実行では 3 つの和集合を固定し、どれかが他を置き換えることはありません。

### 固定の効果

- 固定した原子の力を 0 にするので、その原子は動きません。
- 固定した原子を Hessian から外すので、`freq` は残りの原子で部分 Hessian 振動解析（PHVA）を行います。取り除く剛体運動は [凍結境界での剛体モード](freq.md#凍結境界での剛体モード) にあります。
- `path-opt` と `path-search` の `--mep-mode dmf`（Direct Max Flux）は、力を 0 にする代わりに調和拘束（k = 300 eV/Å²、YAML の `dmf.k_fix`）で固定した原子を留めるので、座標が少し動くことがあります（[`path-opt` の注意点](path-opt.md#使用上の注意点)）。

### 距離を拘束する

`opt` の `--distance-restraint` は、2 原子のあいだに調和の拘束を足します。`(i, j, 目標の距離)`（Å）を渡すか、`(i, j)` で最初の距離を保ちます。

```bash
mlmm opt -i system_layered.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --distance-restraint '[(12,45,2.20)]' --restraint-k 20.0 \
    --out-dir ./result_opt_rest
```

- 原子番号は 1 始まりで、`--zero-based` を付けると 0 始まりになります。`--restraint-k` は力の定数です（既定 300 eV/Å²）。
- `scan`・`scan2d`・`scan3d` では、`--restraint-k` が各段の力の定数です（距離は eV/Å²、角度は eV/rad²）。YAML の `bias.k` より CLI の値が優先されます。各段はその段の座標だけを拘束します。`all` では `--scan-restraint-k` で渡します。スキャンリストの書き方は {ref}`スキャンリスト仕様 <ja-scan-list-spec>` にあります。
- 出力のエネルギーは、拘束のエネルギーを除いた値です。

## 使用上の注意点

- **`--movable-cutoff` は B-factor の MM の層を置き換えます**：このオプションを持つどのコマンドでも同じです。`opt`・`tsopt`・`freq`・`path-opt`・`path-search` では `--detect-layer` も無効になるので、ML 領域を `--model-pdb` か `--model-indices` で渡してください。
- **`--model-pdb` が決めるのは ML 領域だけ**：`--detect-layer` が有効なら、Movable-MM と Frozen-MM は入力の B-factor から決まります。
- **層は `define-layer` で変える**：B-factor を手で書き換えず、`define-layer` をやり直してください。
- **同梱の PDB は chain の欄が空です**。chain が空の PDB では、残基を名前か番号で指定してください。chain のある PDB では、`-c 'A:TYR:44'` のように chain:残基名:番号 で書きます（[残基セレクタ](cli-conventions.md#残基セレクタ)）。
- **リンク水素は自動で付きます**。計算機は、`parm7` の結合のうち片方の端だけが ML 領域にあるものに 1 つずつ水素を置くので、`model.pdb` には入れません。詳しくは [リンク原子の再分配](mlmm-calc.md#リンク原子の再分配) にあります。
- **境界の結合**：リンク水素を置けるのは C–C、C–N、N–C の結合だけです。ほかの `parm7` の結合が境界をまたぐと、`Unsupported ML/MM boundary bond in parm7` で止まります。境界を、置ける結合へ動かしてください。
- **金属・糖鎖・MD の snapshot**：`parm7` を自分で作り、`--parm7` で渡してください（[`mm-parm` の注意点](mm-parm.md#使用上の注意点)）。

## 関連ドキュメント

- [`extract`](extract.md) — 切り出しのオプション、`-c` の残基の指定、非標準の残基名
- [`define-layer`](define-layer.md) — ML・Movable-MM・Frozen-MM の層を付ける
- [`mm-parm`](mm-parm.md) — Amber のトポロジーと、対応する PDB を作る
- [`all`](all.md) — 一括のワークフロー。`-c` で ML 領域と層を作る
- [`opt`](opt.md) — 距離の拘束を付けた構造最適化
- [`scan`](scan.md) — 拘束を付けた段階的なスキャン
- [`freq`](freq.md) — 固定した原子があるときの PHVA と剛体モード
- [反応機構を調べるコツ](mechanism-tips.md) — モデルを広げるとき
- [MLIP の TS を DFT で確かめる](dft-backend.md) — DFT で扱える大きさの ML 領域
- [ML/MM 計算機](mlmm-calc.md) — リンク原子、マイクロイテレーション、MM の Hessian
- {ref}`ML/MM の共通オプション <ja-mlmm-options>` — `--parm7`、`--model-pdb`、`--detect-layer`、`--movable-cutoff`
- [デバイス設定 & HPC セットアップ](device-hpc.md) — GPU のメモリとモデルの大きさ
- [トラブルシューティング](troubleshooting.md) — 切り出しと層のエラー
