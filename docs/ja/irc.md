# `irc`（固有反応座標）

`irc` サブコマンドは、ML/MM の系で最適化した遷移状態（TS）から、EulerPC（Euler 予測子–修正子法）で固有反応座標（IRC）を両方向へたどります。各分岐の軌跡と、2 つの端点の候補を書き出します。この端点を [`opt`](opt.md) で最適化すると、TS がどの反応物（R）と生成物（P）をつなぐかが分かります。

---

## 主な用途

* **TS の確認**: [`tsopt`](tsopt.md) と [`freq`](freq.md)（n_imag = 1）の後に、TS が意図した R と P をつなぐかを確かめる
* **R と P の取得**: 端点を [`opt`](opt.md) で最適化し、この TS の R と P の構造を得る
* **`all` の IRC 段のやり直し**: [`all`](all.md) の IRC を、設定を変えて単独でたどり直す

ML 領域の計算バックエンドにはデフォルトの **UMA**（Meta）のほか、`-b/--backend` で **ORB**、**MACE**、**AIMNet2**、**DFT** も選べます。MM 原子には `--parm7` の Amber 力場を使います。

---

## 基本的な実行例

### 1. 両方向の IRC

TS の構造 `ts.pdb` から両方向へたどり、`--out-json` で結果の要約も書き出します。

```bash
mlmm irc -i ts.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --out-json --out-dir ./result_irc
```

### 2. 順方向だけ

順方向の分岐だけをたどります。

```bash
mlmm irc -i ts.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --no-backward --out-dir ./result_irc_forward
```

### 3. 解析 Hessian

最初の Hessian を、ML バックエンドに有限差分ではなく解析的に計算させます。

```bash
mlmm irc -i ts.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --hessian-calc-mode Analytical --out-dir ./result_irc_analytical
```

### 4. 小さいステップでの再試行

分岐が数フレームで止まるときは、最大ステップを 0.05 bohr にして再試行します。

```bash
mlmm irc -i ts.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --step-size 0.05 --out-dir ./result_irc_small_step
```

### 5. サイクルの上限までたどる

`--never-stop` を付けると、勾配とエネルギーによる停止の条件を無視し、各分岐を `--max-cycles` までたどります。

```bash
mlmm irc -i ts.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --step-size 0.05 --never-stop --max-cycles 250 \
    --out-dir ./result_irc_continue
```

---

## 処理の仕組みと計算仕様

1. **ML/MM の系の組み立て**: `-i` から TS の構造を、`--parm7` から Amber のトポロジーを、`--model-pdb` から ML 領域を読みます（{ref}`ML/MM の共通オプション <ja-mlmm-options>` を参照）。`-q` と `-m` は ML 領域の電荷とスピン多重度です。
2. **出発の方向**: TS で Hessian を計算するか `--read-hess` のファイルから読み、剛体運動を [`freq`](freq.md#凍結境界での剛体モード) と同じように除いてから、`--root` 番目（デフォルト `0`）の固有ベクトルを反応モードとします。そのモードが虚振動でなければ、エラーで止まります。
3. **EulerPC による積分**: 各分岐（順方向、次に逆方向）は TS から始まります。各ステップでは、質量加重の最急降下方向に沿って Euler 予測子で進み、続いて DWI（距離加重補間）面の上で修正 Bulirsch–Stoer 修正子をかけます。予測子の勾配は、Bofill 式で更新する現在の Hessian を使った 2 次の Taylor 展開で見積もります。分岐は、TS の近くを出た後に RMS 勾配が 1 × 10⁻³ hartree/bohr を下回ったとき、エネルギーが上がったとき、1 ステップのエネルギー変化が 1 × 10⁻⁶ hartree 以下になったとき、または `--max-cycles`（デフォルト 125）に達したときに止まります。
4. **経路の書き出し**: 各分岐、TS を通る経路全体、端の構造を書き出します。PDB/mmCIF の入力か `--ref-pdb` があるときは、軌跡と 2 つの端点の候補を PDB にも変換します。

---

## IRC の成否の判定

IRC が収束しなくても、端点の最適化で狙った R と P に着けば、その結果は使えます。

| 確かめること | 見る場所 |
| --- | --- |
| 出発点が TS か | 端末の `Transition vector is mode 0 with wavenumber … cm⁻¹.` の行の波数が負 |
| 各分岐の止まり方 | `result.json` の `forward_integration_converged` / `backward_integration_converged`。RMS 勾配が閾値を下回ったときは `true`、エネルギーで止まったときやサイクルの上限では `false`。理由は `forward_integration_stop_reason` / `backward_integration_stop_reason` に出る |
| 経路に沿って変わる結合 | `result.json` の `bond_changes`（`finished_first` から `finished_last` への `formed` と `broken`） |
| どちらの端が R でどちらが P か | `forward_first.xyz` と `backward_last.xyz` を [`opt`](opt.md) で最適化し、意図した R と P と比べる。順方向 / 逆方向の別では決まらない |

`irc` は端点を判定しないので、端点が狙った R と P かは自分で確かめてください。

端点は `.xyz` なので、原子の順と層を与える TS の PDB を `--ref-pdb` で渡して、両方の端点を `opt` で最適化してください。

```bash
mlmm opt -i result_irc/forward_first.xyz --ref-pdb ts.pdb --parm7 real.parm7 \
    --model-pdb ml_region.pdb -q 0 -m 1 --out-dir ./result_opt_forward
mlmm opt -i result_irc/backward_last.xyz --ref-pdb ts.pdb --parm7 real.parm7 \
    --model-pdb ml_region.pdb -q 0 -m 1 --out-dir ./result_opt_backward
```

端点が意図した R と P でないときは、{ref}`TS が取れないとき <ja-ts-search-fails>` を参照してください。

---

## 主な出力ファイル

実行が終わると、`--out-dir`（デフォルト: `./result_irc/`）に次のファイルができます。

```text
result_irc/
├─ finished_irc_trj.xyz    # TS を通る IRC 経路全体
├─ finished_irc.pdb        # 同じ経路の PDB
├─ finished_first.xyz      # 経路全体の最初のフレーム（順方向を実行したときは forward_first.xyz と同じ構造）
├─ finished_last.xyz       # 経路全体の最後のフレーム（逆方向を実行したときは backward_last.xyz と同じ構造）
├─ forward_irc_trj.xyz     # TS から順方向の分岐（実行したとき）
├─ forward_irc.pdb         # 同じ分岐の PDB
├─ forward_first.xyz       # 順方向の分岐の端（端点の候補）
├─ forward_first.pdb       # 同じ構造の PDB
├─ backward_irc_trj.xyz    # TS から逆方向の分岐（実行したとき）
├─ backward_irc.pdb        # 同じ分岐の PDB
├─ backward_last.xyz       # 逆方向の分岐の端（端点の候補）
├─ backward_last.pdb       # 同じ構造の PDB
└─ result.json             # 結果の要約（--out-json）
```

`.pdb` は、PDB/mmCIF の入力か `--ref-pdb` があるときに書きます。mmCIF の入力と、PDB の欄に入りきらない大きな PDB の入力では、元の識別子を保った `.cif` も書きます（{ref}`mmCIF の入力 <ja-mmcif-input>` を参照）。

* **端点の候補**: `forward_first.xyz` と `backward_last.xyz` を [`opt`](opt.md) で最適化します。各分岐は TS 側のもう一方の端（`forward_last.xyz`、`backward_first.xyz`）も書きます。
* **経路**: `finished_irc_trj.xyz` か `finished_irc.pdb` を PyMOL や VMD で開くと、反応の動きを見られます。
* **要約**: `--out-json` を付けると、`result.json` に各分岐のフレーム数（`n_frames_forward`、`n_frames_backward`）、各分岐の止まり方、`bond_changes`、両端と TS のエネルギー（`energy_first_hartree`、`energy_ts_hartree`、`energy_last_hartree`）、`rigid_projection` に除いた剛体運動と最初の Hessian の情報が記録されます（[JSON 出力の一覧](json-output.md) を参照）。
* **端末**: 各分岐のステップの表と実行時間が出ます。

> **補足:** YAML で `irc.prefix: trial` とすると、`result.json` 以外のファイルの名前が `trial_finished_irc_trj.xyz` のように `trial_` で始まり、`result.json` の `files` にも接頭辞つきの名前が記録されます。YAML の `irc.dump_every` に正の整数を指定すると、実行中に HDF5 のチェックポイント `irc_data.h5` も書きます（デフォルトは書きません）。

---

## 主な CLI オプション

ML/MM の計算コマンドに共通のオプションは {ref}`ML/MM の共通オプション <ja-mlmm-options>` に 1 か所でまとめてあります。下の表は `irc` に固有のものだけです。

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | TS の構造（`.pdb`, `.cif`, `.mmcif`、または `--ref-pdb` と組み合わせた `.xyz`） |
| `-q, --charge` | 整数 | `None` | ML 領域の電荷。`-l` を使う場合のほかは必須 |
| `-l, --ligand-charge` | 文字列 | `None` | 未知のリガンド残基の総電荷（例: `-1`）または残基名ごとの電荷（例: `'GPP:-3,SAM:1'`）。`-q` を省いたときに ML 領域の電荷を求めるのに使用（PDB/mmCIF 入力または `--ref-pdb`） |
| `-m, --multiplicity` | 整数 | `1` | ML 領域のスピン多重度（2S+1） |
| `--max-cycles` | 整数 | `125` | 分岐ごとの IRC ステップの上限 |
| `--step-size` | 実数 | `0.10` | 最大ステップ長（bohr、質量加重しない Cartesian 座標） |
| `--root` | 整数 | `0` | 反応モードとする Hessian の固有ベクトル。固有値の昇順に 0 から数える |
| `--forward/--no-forward` | フラグ | `True` | 順方向の分岐を実行 |
| `--backward/--no-backward` | フラグ | `True` | 逆方向の分岐を実行 |
| `--never-stop/--no-never-stop` | フラグ | `False` | 勾配とエネルギーによる停止の条件を無視し、`--max-cycles` までたどる |
| `-o, --out-dir` | パス | `./result_irc/` | 出力先ディレクトリ |
| `--hessian-calc-mode` | `FiniteDifference` / `Analytical` | `FiniteDifference` | ML バックエンドが最初の Hessian を計算する方法 |
| `-b, --backend` | 文字列 | `uma` | ML 領域のバックエンド（`uma`, `orb`, `mace`, `aimnet2`, `dft`） |
| `--read-hess` | パス | `None` | Hessian を計算せず、`.npy` ファイル（`freq` や `tsopt --dump-hess` で書いたものなど）から読んで始める |
| `--out-json/--no-out-json` | フラグ | `False` | 結果の要約を `result.json` に出力（[JSON 出力の一覧](json-output.md)） |

全オプションは [自動生成のオプションの一覧（英語のみ）](../reference/commands/irc.md) を参照してください。

> **補足:** YAML（`--config`）の `irc` ブロックのキーは、YAML 設定の一覧の {ref}`irc <ja-irc-section>` にすべて載っています。

---

## 使用上の注意点

* **すぐ止まる分岐**: 分岐がサイクルの上限より前に 3 フレーム以下で終わると、端末に `[irc] IRC stopped after only a few frames in …` の警告が出ます。ステップが大きすぎると EulerPC が不安定になることがあるので、ほかの設定を変える前に、小さい `--step-size`（例: `0.05`）で再試行してください。
* **`--never-stop` はデフォルト無効**: 有効にすると、物理的な端点を過ぎてもサイクルの上限まで進みます。数値的な失敗や外部からの中断では止まります。軌跡を確かめて端点を最適化し、先の経路が役に立つときだけ `--max-cycles` を増やしてください。
* **`--root` は 0 から数える**: TS 最適化が成功すると、反応モードの虚振動が 1 つ出るので、n_imag = 1 の TS では `--root 0`（ただ 1 つの負の固有値）のままにしてください。`1`、`2` などは、反応モードより固有値の小さい（より負の）疑似モードがあると分かっているときだけ使います。
* **Cartesian 座標**: YAML の `geom.coord_type` にかかわらず、`irc` は Cartesian 座標を使います。
* **`--read-hess` のファイル**: [`freq`](freq.md) と同じ `.npy` ファイルで、単位は Hartree/bohr²、全原子か Hessian の計算に入る原子だけの分を持ちます。同じ構造・電荷・多重度・計算機で計算した Hessian を渡してください。`irc.hessian_init: calc`（デフォルト）が必要です。ファイルを使ったときは、`result.json["rigid_projection"]["hessian_source"]` が `"file"` になります。
* **解析 Hessian と `--uma-workers`**: UMA では、`--hessian-calc-mode Analytical` は 1 より大きい `--uma-workers` と併用できず、エラーで止まります。解析 Hessian には `--uma-workers 1` を使ってください（[バックエンド](backends.md) を参照）。速度とメモリ量はバックエンドと系によって変わるので、先に対象の系で両方を比べてください。
* **凍結原子**: 凍結 MM 層のほかに、`--freeze-atoms` でほかの原子（1 始まり）も凍結できます。選び方は {ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>` を参照してください。
* **大きな系**: `--hess-device cpu` を付けると、最初の Hessian と IRC の Hessian の演算を CPU で行い、GPU のメモリに収めます。
* **分岐は少なくとも 1 つ**: `--no-forward` と `--no-backward` を両方付けると、エラーで止まります。
* **1 回に 1 構造**: `-i` には 1 つの構造を指定します。軌跡からは、使うフレームを先に `.xyz` に切り出し、`--ref-pdb` と一緒に渡してください。
* **設定の優先順位**: デフォルト < YAML < コマンドライン（[CLI 規約](cli-conventions.md) を参照）。

---

## 関連ドキュメント

* [tsopt](tsopt.md) — IRC の前に TS を最適化する
* [freq](freq.md) — TS の虚振動が 1 つ（n_imag = 1）であることを確かめる
* [opt](opt.md) — IRC の端点を R と P へ最適化する
* [all](all.md) — `tsopt` の後に IRC を実行し、端点まで最適化する一連のワークフロー
* [トラブルシューティング](troubleshooting.md) — 実行が失敗したときの切り分け
* [YAML 設定の一覧](yaml-reference.md) — `irc` のすべての設定
* [用語集](glossary.md) — IRC などの用語
* [終了コード](cli-conventions.md#終了コード) — 終了ステータスの意味
