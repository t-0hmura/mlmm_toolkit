# `tsopt`（遷移状態の構造最適化）

`tsopt` サブコマンドは、層を定義した ML/MM の酵素モデルで遷移状態（TS）の候補構造を 1 次の鞍点へ最適化し、final geometry で Hessian を計算して虚振動数の本数（n_imag）を数えます。TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます。

---

## 主な用途

* **TS 候補の仕上げ**: [`path-opt`](path-opt.md) / [`path-search`](path-search.md) の最高エネルギーのイメージ（HEI）や [`scan`](scan.md) の頂点を、最適化した TS に仕上げる
* **自作の構造の検証**: 手で作った候補が TS か（n_imag = 1）を確かめ、反応モードをアニメーションで確認する
* **`all` の TS 段のやり直し**: [`all`](all.md) で得た TS を、設定を変えて単独で最適化し直す

ML 領域の計算バックエンドのデフォルトは、Meta が公開した学習済みの[機械学習原子間ポテンシャル（MLIP）](backends.md) の **UMA** です。`-b/--backend` で **ORB**、**MACE**、**AIMNet2**、**DFT** も選べます。MM 領域には `--parm7` の Amber パラメータを使います。

候補がまだ無い場合は、先に次のコマンドで作ってください。

| 手元にあるもの | 候補を作るコマンド |
| --- | --- |
| 反応物**と**生成物 | [`path-opt`](path-opt.md)（2 構造。`hei.xyz`）または [`path-search`](path-search.md)（2 構造以上。結合が変わるセグメントごとに `hei_seg_NN.xyz`） |
| 反応物だけ、または動かしたい結合がある | [`scan`](scan.md) で反応する距離を少しずつ動かし、ほかの自由度を緩和する |

候補を `tsopt` で最適化し、得られた TS から [`irc`](irc.md) を実行してください。[`all`](all.md) を使えば、これらを一度に実行できます。

---

## 基本的な実行例

以下の例では、`ts_guess.pdb` が `real.parm7` に対応する全系の候補構造、`ml_region.pdb` が ML 領域の定義です（[ML 領域と層の組み方](model-setup.md) を参照）。

### 1. 標準の実行（RS-P-RFO）

ML 領域の電荷とスピン多重度を明示して実行します。

```bash
mlmm tsopt -i ts_guess.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --out-dir ./result_tsopt
```

### 2. Dimer 法

完全な Hessian を繰り返し計算するのが重い場合や、難しい候補で別の方法を試したい場合に使います。

```bash
mlmm tsopt -i ts_guess.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --opt-mode dimer --out-dir ./result_tsopt_dimer
```

### 3. 余分な虚振動の除去

候補に虚振動が 2 つ以上ある場合は `--flatten` を付けます。

```bash
mlmm tsopt -i ts_guess.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --flatten --out-dir ./result_tsopt_flatten
```

### 4. 保存した Hessian から開始

同じ構造で `freq` などが `--dump-hess` で保存した Hessian を読み込み、計算し直さずに始めます。

```bash
mlmm tsopt -i ts_guess.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --read-hess ts_guess_hess.npy --out-dir ./result_tsopt
```

---

## 処理の仕組みと計算仕様

1. **モデルの読み込みと境界の凍結**: ML 領域を `--model-pdb` から、可動 MM 層と凍結 MM 層を入力 PDB の B-factor から読み込みます。ML 領域の電荷は `-q` または `-l` から決まります（{ref}`電荷の指定 <ja-charge-specification>` を参照）。凍結 MM 層と `--freeze-atoms` で指定した原子は固定したままで、最後の Hessian は `--active-dof-mode` で選んだ原子だけで扱います（PHVA: 部分 Hessian 振動解析）。`.xyz` の候補には、原子の順序と層を与える全系の PDB を `--ref-pdb` で渡してください。
2. **最適化法の選択**（`--opt-mode`）: `hess`（デフォルト）は完全な Hessian を使う **RS-P-RFO**（制限ステップ分割有理関数最適化）を実行し、`rsirfo` と `trim` はそれぞれ RS-I-RFO（restricted-step image RFO）と TRIM（trust-region image minimization）を選びます。`dimer`（または `grad`）は **Hessian-guided Dimer** 法で、勾配を使って最低固有モードを追い、ときどき厳密な Hessian で方向を更新します。
3. **反応モードに沿った探索**: 反応モードの方向にはエネルギーを上り、それ以外の方向には下りながら、収束条件（`--thresh`）を満たすまで構造を動かします。デフォルトの `baker` は、力の最大値 3 × 10⁻⁴ 未満、力の RMS 2 × 10⁻⁴ 未満、ステップの最大値 3 × 10⁻⁴ 未満、ステップの RMS 2 × 10⁻⁴ 未満（原子単位）、エネルギー変化 10⁻⁶ hartree 未満の 5 つをすべて同時に求め、どれも Gaussian の既定（`gau`）より厳しい条件です。RS-P-RFO は Bofill 式で Hessian を更新し、1 ステップを信頼半径 0.1 bohr（`rsirfo.trust_max`）以内に収めます。マイクロイテレーション（`--microiter`）はデフォルトで有効で、各 macro ステップで ML 領域とそれに結合した境界の MM 原子を動かし、続いて可動 MM 原子を L-BFGS で緩和します。マイクロイテレーションでは最適化の行が `[microiter]` で始まり、Dimer 法と `--no-microiter` では `[tsopt]` の形で出ます。
4. **最後の確認**: 収束すると、final geometry で Hessian を計算して n_imag を数え、各虚振動モードをアニメーションとして書き出します。凍結原子の扱いは [`freq`](freq.md#凍結境界での剛体モード) と同じです。
5. **`--flatten` による余分な虚振動の除去**: 虚振動が 2 つ以上残る場合は、余分なモードに沿って構造をずらして最適化し直し、1 つになるか回数の上限に達するまで繰り返します。Dimer 法では、各回で厳密な Hessian を計算し直してダイマー方向を更新し、短い Dimer + L-BFGS の区間を実行します。

---

## TS の判定

結果は、実行がどう終わったかで決まります。

| 終わり方 | 端末の判定の行 | `tsopt` が残すもの | `all` の次の動作 |
| --- | --- | --- | --- |
| 収束 | `[microiter] Converged!` または `[tsopt] Numerical optimization converged.` の後に `[tsopt] Wrote N final imaginary mode(s).`。n_imag ≥ 2 では `[tsopt] WARNING: Higher-order stationary point (n_imag=N). …` も出て、n_imag = 0 では `[tsopt] No imaginary mode detected. …` | final geometry、n_imag、虚振動モード | n_imag ≥ 1 なら IRC へ進み、n_imag = 0 なら IRC の前で止まる |
| エネルギーが変わらなくなって停止 | `[microiter] Stalled (energy plateau; not converged)` または `[tsopt] Stalled (energy plateau; not converged)` | final geometry と n_imag | IRC の前で止まる |
| `--max-cycles` に達して未収束 | `[microiter] Reached max macro iterations (M).` または `[tsopt] Reached max cycles (N/M).` | final geometry（Hessian なし） | IRC の前で止まる |
| `--skip-final-freq` を付けて収束 | `[tsopt] WARNING: TS saddle-point order is not verified (--skip-final-freq).` | final geometry（Hessian なし） | 反応モードを確かめられないため、IRC の前で止まる |
| 最後の Hessian の計算に失敗 | `[tsopt] ERROR: Terminal PHVA failed.` | final geometry と、`hessian_status: failed` とその理由 | IRC の前で止まる |

n_imag は次のように読みます。

| n_imag | 意味 |
| --- | --- |
| 1 | 1 次の鞍点。モードが狙った原子を動かしているかを確かめてから、[`irc`](irc.md) を実行してください |
| 0 | 虚振動なし。構造が極小点の側へ緩和しています |
| 2 以上 | 高次の鞍点。余分な虚振動が残っています。`all` はそれでも、最適化が追っていた虚振動のモード（それが使えないときは最も低い虚振動のモード）に沿って IRC を流すので、IRC の端点でそのモードがどこへつながるかを確かめられます |

ν < −5.00 cm⁻¹ のモードを虚振動として数えます。−5.00 以上 0 cm⁻¹ 未満の値は数値誤差として扱います。

(ja-wrong-imaginary-mode-count)=
### 最適化後に虚振動数の本数が誤っている場合

n_imag が 1 でない場合や、モードが狙った反応の原子を動かしていない場合は、次を試してください。これらは組み合わせて使えます。

| 結果 | 試すこと |
| --- | --- |
| n_imag = 0 | 候補が鞍点から遠い状態です。経路探索や scan でよりよい候補を作ってください。`all` では `--refine-path` で再帰的な `path-search` を実行し、HEI を細かく求め直せます。増えた素過程のそれぞれに TS 最適化と IRC がかかるため、計算量は増えます |
| n_imag ≥ 2 | 各モードの動きを確かめてください。`--flatten` を付けて最適化し直すか、この候補で `--precision fp32` / `fp64` を比べてください |
| モードは 1 つだが動きが違う | どの原子が動くかを確かめ、狙った反応に近い候補から始めてください |

例として、fp64 で、flatten を有効にしてやり直す場合は次のようにします。

```bash
mlmm tsopt -i ts_guess.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 \
    --precision fp64 --flatten -o result_tsopt
```

ほかの手は {ref}`TS が取れないとき <ja-ts-search-fails>` に、そのほかの失敗は [トラブルシューティング](troubleshooting.md) にあります。

---

## 主な出力ファイル

実行が終わると、`--out-dir`（デフォルト: `./result_tsopt/`）に次のファイルができます。

```text
result_tsopt/
├─ final_geometry.xyz               # final geometry（常に出力）
├─ final_geometry.pdb               # 同じ構造の PDB（B-factor が層を表す）
├─ vib/
│  ├─ imag_01_-385.20cm-1_trj.xyz   # 虚振動モードごとのアニメーション
│  └─ imag_01_-385.20cm-1.pdb       # 同じアニメーションの PDB
├─ optimization_all_trj.xyz         # 最適化の軌跡（--dump）
├─ optimization_all.pdb             # 同じ軌跡の PDB（--dump）
├─ .dimer_mode.dat                  # Dimer の今の方向（Dimer 法のとき）
└─ result.json                      # 結果の要約（--out-json）
```

mmCIF の入力と、PDB の欄に入りきらない大きな PDB の入力では、元の識別子を保った `.cif` も書きます（{ref}`mmCIF の入力 <ja-mmcif-input>` を参照）。

* **final geometry**: `final_geometry.*` を、[`irc`](irc.md) に渡す TS として使います。`final_geometry.pdb` の B-factor は、ML 領域が 0、可動 MM 原子が 10、凍結 MM 原子が 20 です。
* **反応モード**: `vib/imag_*_trj.xyz` を PyMOL や VMD で開き、生成・切断される結合に沿って原子が動いているかを確かめてください。
* **要約**: `--out-json` を付けると、`result.json` に終わり方（`optimization_status`: `converged`、`stalled`、`not_converged`）、`hessian_status`、n_imag が記録されます（[JSON 出力リファレンス](json-output.md#tsopt) を参照）。

---

## 主な CLI オプション

ML/MM の計算コマンドに共通のオプションは {ref}`ML/MM の共通オプション <ja-mlmm-options>` に 1 か所でまとめてあります。下の表は `tsopt` に固有のものだけです。

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | 全系の構造 1 つ（`.pdb`, `.cif`, `.mmcif`、または `--ref-pdb` と組み合わせた `.xyz`）。軌跡は 1 フレームを `.xyz` に切り出してから指定（{ref}`軌跡から 1 フレームを取り出す <ja-trajectory-one-frame>` を参照） |
| `-q, --charge` | 整数 | `None` | ML 領域の電荷。`-l` を使う場合のほかは必須 |
| `-m, --multiplicity` | 整数 | `1` | ML 領域のスピン多重度（2S+1） |
| `-l, --ligand-charge` | 文字列 | `None` | 未知のリガンドの総電荷、または残基名ごとの電荷（例: `'GPP:-3,SAM:1'`）。`-q` を省いたときに ML 領域の電荷を求めるのに使用（PDB 入力または `--ref-pdb`） |
| `-o, --out-dir` | パス | `./result_tsopt/` | 出力先ディレクトリ |
| `-b, --backend` | 文字列 | `uma` | ML 領域のバックエンド（`uma`, `orb`, `mace`, `aimnet2`, `dft`） |
| `--opt-mode` | `hess` / `dimer` / `rsirfo` / `trim` | `hess` | 最適化法: RS-P-RFO / Dimer / RS-I-RFO / TRIM（`rsprfo` = `hess`、`grad` = `dimer`）。`opt` では `grad` は L-BFGS を指す（{ref}`コマンドごとの --opt-mode <ja-opt-mode-semantics>` を参照） |
| `--ref-mode` | パス | `None` | 反応モードの参照方向（`.npz`, `.npy`, テキスト）。`all` が MEP から渡すもので、通常は指定しない。Dimer 法では使わない |
| `--microiter/--no-microiter` | フラグ | `True` | 各 macro TS ステップと、可動 MM 原子の L-BFGS 緩和を交互に実行。Dimer 法では使わず、`--embedcharge` では通常の最適化に切り替える |
| `--hessian-calc-mode` | `FiniteDifference` / `Analytical` | `FiniteDifference` | ML 領域の Hessian の計算法（有限差分 / 解析的） |
| `--flatten/--no-flatten` | フラグ | `False` | 余分な虚振動を除く |
| `--freeze-atoms` | 文字列 | `None` | 追加で凍結する原子（1 始まり、カンマ区切り: 例 `'1,3,5'`） |
| `--active-dof-mode` | `all` / `ml-only` / `partial` / `unfrozen` | `partial` | 最後の振動解析に入れる原子: 全原子 / ML 領域だけ / ML 領域と可動 MM 原子 / 凍結層以外のすべての原子 |
| `--thresh` | プリセット | `baker` | 収束条件（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） |
| `--max-cycles` | 整数 | `100000` | 最適化サイクルの上限 |
| `--stop-plateau/--no-stop-plateau` | フラグ | `False` | エネルギーが変わらなくなったら（直近 50 サイクルの幅が 1e-4 hartree 未満）止め、Hessian を計算 |
| `--skip-final-freq/--no-skip-final-freq` | フラグ | `False` | 収束後の最後の Hessian を省く |
| `--read-hess` | パス | `None` | Hessian を計算せず `.npy` ファイルから読んで開始（Cartesian、Hartree/bohr²、全原子または可動原子だけ） |
| `--dump-hess` | パス | `None` | final geometry の Hessian を `.npy` ファイルに保存（`freq`・`tsopt`・`irc` の `--read-hess` 用）。最後の Hessian を計算したときだけ書く |
| `--precision` | `fp32` / `fp64` | バックエンドごと（`uma`: `fp32`、`orb`・`mace`: `fp64`） | バックエンドの精度。`aimnet2` は `fp64` を受け付けない（{ref}`MLIP バックエンド: 精度 <ja-precision>` を参照） |
| `--coord-type` | `cart` / `redund` / `dlc` / `tric` | `cart` | 最適化に使う座標系：デカルト座標 / 冗長内部座標 / 非局在化内部座標（DLC）/ 並進・回転を含む内部座標（TRIC） |
| `--config` | パス | `None` | コマンドラインのオプションより前に適用する YAML ファイル |
| `--dump/--no-dump` | フラグ | `False` | 最適化の軌跡 `optimization_all_trj.xyz` を書き出す |
| `--out-json/--no-out-json` | フラグ | `False` | 結果の要約を `result.json` に出力（[JSON 出力リファレンス](json-output.md)） |

全オプションは `mlmm tsopt --help-advanced` または [自動生成 CLI リファレンス](../reference/commands/tsopt.md) を参照してください。

> **補足:** YAML では、Dimer 法は `hessian_dimer:` ブロックを読み、RS-P-RFO・RS-I-RFO・TRIM は `rsirfo:` ブロックを共用し、マイクロイテレーションは `microiter:` ブロックを読みます。キーの一覧は YAML リファレンスの [`rsirfo`](yaml-reference.md#rsirfo)、[`hessian_dimer`](yaml-reference.md#hessian_dimer)、[`microiter`](yaml-reference.md#microiter) にあります。

> **補足:** 最適化の途中で反応モードが別の Hessian 固有ベクトル（root）に入れ替わる場合は、`rsirfo.track_mode_by_overlap: true` を設定してください。

> **補足:** 収束が遅い場合は、`rsirfo.hessian_recalc`（デフォルト `500`）を 50〜200 に下げてください。厳密な Hessian を計算し直す間隔が短くなり、計算は増えますが収束しやすくなります。

---

## 使用上の注意点

(ja-flatten-precedence-caveat)=
### `--flatten` を使うとき

`--flatten` はデフォルトで無効です。flatten の回数は、Dimer 法でも RS-P-RFO・RS-I-RFO・TRIM でも、YAML の 1 つのキー `hessian_dimer.flatten_max_iter` で決まります。

| コマンドライン | flatten の回数 |
| --- | --- |
| `--flatten` も `--no-flatten` も付けない | `0`（無効）。YAML で `hessian_dimer.flatten_max_iter` を指定した場合はその値 |
| `--flatten` | YAML の値（正の値の場合）、無ければ `50` |
| `--no-flatten` | YAML に値があっても `0` |

`--flatten` は余分な虚振動を除きますが、欠けている反応モードを作ることはできません。n_imag = 0 の場合は、よりよい候補を作ってください。

### そのほかの注意

* **上り方向のステップは常に許可**: 鞍点探索では反応モードに沿ってエネルギーを上る必要があるため、YAML で指定しても `tsopt` は `reject_uphill: false` を保ちます。`--reject-uphill/--no-reject-uphill` は、`opt` と `all` の端点の最適化で使うフラグです。
* **生成物側から scan した障壁**: この候補を作った scan が生成物から始まった場合、障壁の読み方は {ref}`scan: スキャン方向とバリアの符号 <ja-scan-direction-barrier-sign>` を参照してください。
* **追う固有ベクトル（root）は 1 つ**: 最適化は 1 つの固有ベクトルに沿って上ります（`0` が最も低い固有値）。`rsirfo.roots: [0]` のように 1 要素のリストで指定します。Dimer 法では `hessian_dimer.root` を使います。`tsopt` に `--root` フラグはありません。
* **そのほかの RS-P-RFO の設定**: `trust_norm: max_atom` はステップ全体ではなく原子ごとの変位を制限し（Cartesian 座標だけ）、`hessian_update: ts_bfgs` は Bofill の代わりに TS-BFGS で Hessian を更新します。どちらも信頼半径は変えません。
* **追加の探索は指定したときだけ**: 収束後は、n_imag が 1 でなくても、自動では追加の探索をしません。`--flatten` を使うか、`rsirfo.saddle_recovery_max_cycles` を `0` より大きくしてください（デフォルト `0`）。後者では、厳密な Hessian に虚振動が無いとき、RS-P-RFO・RS-I-RFO・TRIM がエネルギーを上る向きにステップを進めます。
* **併用できない組み合わせ**: `--skip-final-freq` と `--dump-hess`、2 以上の `--uma-workers`（MLIP の並列ワーカー数）と `--hessian-calc-mode Analytical`（[ワーカーと Hessian の計算方式](backends.md#ワーカーと-hessian-の計算方式) を参照）。
* **`--skip-final-freq` と `--flatten`**: RS-P-RFO・RS-I-RFO・TRIM では、`--skip-final-freq` を付けると、最後の Hessian を使う `--flatten` も省かれます。
* **`--read-hess` を RS-P-RFO・RS-I-RFO・TRIM で使う場合**: ファイルの Hessian が最初の厳密な Hessian の代わりになるので、`rsirfo.hessian_init` はデフォルトの `calc` のままにしてください。ほかの値ではエラーで止まります。
* **マイクロイテレーションの MM 緩和の収束条件**: `microiter.micro_thresh` で MM 緩和の収束条件を指定します。指定しないときは macro ステップと同じ条件を使います。
* **エネルギーの停滞とマイクロイテレーション**: `--stop-plateau` が見るのは macro ステップで、MM 緩和を止めることはありません。
* **`--ref-mode` と凍結原子**: `--ref-mode` は MEP から反応の方向を与えるだけで、凍結境界の扱いは変えません。
* **全原子凍結の禁止**: すべての原子を凍結すると、`tsopt` はエラーで停止します。
* **エラーでの停止**: 入力や構造が不正な場合や、最適化を続けられないエラーが起きた場合は、エラーメッセージを出して止まります。それまでに書いたファイルだけが残ります。
* **設定の優先順位**: デフォルト < YAML < コマンドライン（{ref}`設定の優先順位 <ja-configuration-precedence>` を参照）。

---

## 関連ドキュメント

* [irc](irc.md) — 最適化した TS からの反応経路の追跡
* [freq](freq.md) — 完全な振動解析と熱化学補正
* [path-opt](path-opt.md) / [path-search](path-search.md) / [scan](scan.md) — TS 候補の作成
* [all](all.md) — モデル作成・MEP・TS 最適化・IRC・振動解析を一度に実行するワークフロー
* [反応機構を調べるコツ](mechanism-tips.md) — TS が取れないときに試すこと
* [トラブルシューティング](troubleshooting.md) — 実行が失敗したときの切り分け
* [YAML リファレンス](yaml-reference.md) — `rsirfo`・`hessian_dimer`・`microiter` のすべての設定
* [用語集](glossary.md) — TS、Dimer、Hessian などの用語
* {ref}`終了コード <ja-exit-codes>` — 終了ステータスの意味
