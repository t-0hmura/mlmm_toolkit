# `opt`（構造最適化）

`opt` サブコマンドは、層を定義した ML/MM の酵素モデルの構造 1 つを、局所極小点へ最適化します。最適化法は L-BFGS（`--opt-mode grad`、デフォルト）と RFO（`--opt-mode hess`）から選べます。

---

## 主な用途

* **R・P・中間体の準備**: 経路探索や振動解析の前に、反応物・生成物・中間体の構造を緩和し、[`freq`](freq.md) で極小点（n_imag = 0）であることを確かめる
* **距離を保った緩和**: 選んだ原子の組の距離を保ったまま、ほかの自由度を緩和する
* **IRC の端点から R と P へ**: [`irc`](irc.md) の端点を、それぞれがつながる極小点まで最適化する
* **MM による事前緩和**: ML/MM の最適化の前に、MM 力場だけで全系を緩和する（`--mm-only`）

ML 領域の計算バックエンドにはデフォルトの **UMA**（Meta）のほか、`-b/--backend` で **ORB**、**MACE**、**AIMNet2**、**DFT** も選べます。MM 原子には `--parm7` の Amber 力場を使います。

---

## 基本的な実行例

### 1. 標準の最小化

全系 `system_layered.pdb` を、Amber のトポロジー `real.parm7` と ML 領域 `ml_region.pdb` で最適化し、`--out-json` で結果の要約も書き出します。

```bash
mlmm opt -i system_layered.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --out-json --out-dir ./result_opt
```

端末に `[opt] Converged!` が出て、`result_opt/result.json` の `"optimization_status"` が `"converged"` であれば収束しています。

### 2. 厳しい収束条件と軌跡の保存

収束条件を `gau_tight` にし、最適化の軌跡を残します。

```bash
mlmm opt -i system_layered.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --thresh gau_tight --dump --out-dir ./result_opt_tight
```

### 3. 距離拘束

弱い調和拘束（20 eV·Å⁻²）で、原子 12 と 45 の距離を 2.20 Å へ近づけます。

```bash
mlmm opt -i system_layered.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --distance-restraint '[(12,45,2.20)]' --restraint-k 20.0 \
    --out-dir ./result_opt_rest
```

### 4. RFO とマイクロイテレーション

`--opt-mode hess` で、厳密な Hessian から始める RFO に切り替えます。マイクロイテレーションはデフォルトで有効です。

```bash
mlmm opt -i system_layered.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --opt-mode hess --out-dir ./result_opt_rfo
```

---

## 処理の仕組みと計算仕様

1. **ML/MM の系の組み立て**: `-i` から全系の構造を、`--parm7` から Amber のトポロジーを、`--model-pdb` から ML 領域を読みます。残りの原子は {ref}`可動 MM 原子か凍結 MM 原子 <ja-mlmm-options>` になります。`-q` と `-m` は ML 領域の電荷とスピン多重度です。`--freeze-atoms` でほかの原子も凍結できます。
2. **最適化法の選択**（`--opt-mode`）: `grad`（別名 `lbfgs`）は勾配だけを使う **L-BFGS** を実行します。`hess`（別名 `rfo`）は **RFO** を実行し、厳密な Hessian から始めて TS-BFGS 式で更新し（YAML の [`rfo.hessian_update`](yaml-reference.md#rfo) のデフォルト）、500 サイクルごとに計算し直します。`hess` のマイクロイテレーションでは、ML 原子とリンク原子の MM 側の親原子を動かす RFO の 1 ステップと、ほかの可動 MM 原子を MM の力だけで動かす L-BFGS の緩和とを交互に行います。Gaussian のマイクロイテレーションと同じ方式です。
3. **距離拘束の追加**（`--distance-restraint`）: `(i, j, target)` のそれぞれが、力の定数 `--restraint-k`（eV·Å⁻²）の調和項を加え、原子 i と j の距離を `target`（Å）へ引き寄せます。`(i, j)` は最初の距離を保ちます。番号は 1 始まりで、`--zero-based` を付けると 0 始まりになります。
4. **最小化**: 収束条件を満たすか `--max-cycles` に達するまで構造を動かします。デフォルトの `--thresh gau` は、力の最大値が 4.5 × 10⁻⁴、RMS が 3.0 × 10⁻⁴ hartree/bohr 未満、ステップの最大値が 1.8 × 10⁻³、RMS が 1.2 × 10⁻³ bohr 未満を求め、Gaussian の既定と同じ条件です。
5. **`--flatten` による虚振動の除去**: 最適化の後に Hessian を計算し、すべての虚振動モード（ν < −5.00 cm⁻¹）に沿って構造を 0.10 Å ずらして最適化し直します。虚振動が無くなるか 50 回に達するまで繰り返します。`--flatten` では、各回の後に端末の `[Imaginary modes] n=…` の行に n_imag が出て、最後の回の後にも虚振動が残ると `[flatten] WARNING: Remaining imaginary modes after the flatten loop: N` が出ます。

---

## 収束の判定

実行の終わり方は、端末と `result.json`（`--out-json`）に出ます。

| 終わり方 | `optimization_status` | 端末の行 | `scientific_status` / 終了コード |
| --- | --- | --- | --- |
| 収束 | `converged` | `[opt] Converged!` | `success` / 0 |
| `--max-cycles` に達して未収束 | `not_converged` | `[opt] Reached max cycles (N/M).` | `failed` / 1 |
| エネルギーが変わらなくなって停止（`--stop-plateau`） | `stalled` | `[opt] Stalled (energy plateau; not converged)` | `failed` / 1 |

どの行の後にも `[opt] Total cycles: N` が出ます。`stalled` は収束ではありません。力の収束条件を満たさないまま、エネルギーが変わらなくなった状態です。

収束して得られるのは停留点で、極小点とは限りません。`opt` は `--flatten` のとき以外は最後の Hessian を計算しないので、final geometry に [`freq`](freq.md) を実行し、n_imag = 0 を確かめてください。

---

## 主な出力ファイル

実行が終わると、`--out-dir`（デフォルト: `./result_opt/`）に次のファイルができます。

```text
result_opt/
├─ final_geometry.xyz        # final geometry（常に出力）
├─ final_geometry.pdb        # 同じ構造の PDB（PDB/mmCIF 入力または --ref-pdb）
├─ optimization_trj.xyz      # 最適化の軌跡（--dump）
├─ optimization.pdb          # 同じ軌跡の PDB（--dump）
├─ optimization_all_trj.xyz  # 最適化の全ステップをつないだ軌跡（--dump）
├─ optimization_all.pdb      # 同じ軌跡の PDB（--dump）
├─ restart_NNN.yaml          # オプティマイザの状態（--dump と YAML の opt.dump_restart）
└─ result.json               # 結果の要約（--out-json）
```

{ref}`mmCIF の入力 <ja-mmcif-input>` と、PDB の欄に入りきらない大きな PDB の入力では、元の識別子を保った `.cif` も書きます。

* **final geometry**: `final_geometry.*` が最適化した構造です。[`freq`](freq.md) や経路探索に渡してください。
* **要約**: `--out-json` を付けると、[`result.json`](json-output.md) に `optimization_status`、最後のエネルギー `energy_hartree`（拘束のエネルギーを除いた値）、サイクル数 `n_opt_cycles` が記録されます。マイクロイテレーションでは、MM の緩和のサイクル数 `n_micro_cycles` も記録されます。
* **端末**: サイクルごとの表と実行時間が出ます。

---

## 主な CLI オプション

ML/MM の計算コマンドに共通のオプションは {ref}`ML/MM の共通オプション <ja-mlmm-options>` に 1 か所でまとめてあります。下の表は `opt` に固有のものだけです。

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | 全系の入力構造ファイル（`.pdb`, `.cif`, `.mmcif`、または `--ref-pdb` と組み合わせた `.xyz`） |
| `-q, --charge` | 整数 | `None` | ML 領域の電荷。`-l` を使う場合のほかは必須 |
| `-l, --ligand-charge` | 文字列 | `None` | 未知のリガンド残基の総電荷（例: `-1`）または残基名ごとの電荷（例: `'GPP:-3,SAM:1'`）。`-q` を省いたときに ML 領域の電荷を求めるのに使用（PDB/mmCIF 入力または `--ref-pdb`） |
| `-m, --multiplicity` | 整数 | `1` | ML 領域のスピン多重度（2S+1） |
| `-b, --backend` | 文字列 | `uma` | ML 領域のバックエンド（`uma`, `orb`, `mace`, `aimnet2`, `dft`） |
| `--opt-mode` | `grad` / `hess` | `grad` | 最適化法: L-BFGS / RFO（`lbfgs` と `rfo` は別名） |
| `--microiter/--no-microiter` | フラグ | `True` | `hess` で、ML 領域の RFO のステップと可動 MM 原子の L-BFGS の緩和とを交互に行う |
| `--mm-only/--no-mm-only` | フラグ | `False` | MM 力場だけで全系を最小化する（`grad` のときだけ） |
| `--thresh` | プリセット | `gau` | 収束条件（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`。下の表を参照） |
| `--max-cycles` | 整数 | `100000` | 最適化サイクルの上限。`--flatten` の各回と共有 |
| `--coord-type` | `cart` / `redund` / `dlc` / `tric` | `cart` | 最適化の座標系。ML/MM では `cart` のままにする |
| `--dump/--no-dump` | フラグ | `False` | 軌跡 `optimization_trj.xyz` と `optimization_all_trj.xyz` を書き出す |
| `--distance-restraint` | 文字列 | `None` | 調和の距離拘束。直接書く（`'[(i,j,target_Å),...]'`）か、YAML/JSON ファイルで指定。`(i,j)` は最初の距離を保つ |
| `--restraint-k` | 実数 | `300` | 距離拘束の力の定数（eV·Å⁻²） |
| `--one-based/--zero-based` | フラグ | `--one-based` | `--distance-restraint` の番号を 1 から数えるか 0 から数えるか |
| `--freeze-atoms` | 文字列 | `None` | 凍結する原子（1 始まり、カンマ区切り: 例 `'1,3,5'`） |
| `--hessian-cutoff` | 実数 | `None` | ML 領域からこの距離（Å）以内の可動 MM 原子だけを Hessian に入れる。デフォルトでは可動 MM 原子すべて |
| `--flatten/--no-flatten` | フラグ | `False` | 最適化の後に虚振動を除く |
| `--reject-uphill/--no-reject-uphill` | フラグ | `False` | `hess` で、エネルギーが 1e-4 hartree を超えて上がる RFO のステップを捨て、信頼半径を縮める |
| `--stop-plateau/--no-stop-plateau` | フラグ | `False` | エネルギーが変わらなくなったら（直近 50 サイクルの幅が 1e-4 hartree 未満）止め、`stalled` と報告 |
| `-o, --out-dir` | パス | `./result_opt/` | 出力先ディレクトリ |

全オプションは [自動生成のオプションの一覧（英語のみ）](../reference/commands/opt.md) を参照してください。

`--thresh` のプリセットは次の上限を決めます（力は hartree/bohr、ステップは bohr）。

| プリセット | 力の最大値 | 力の RMS | ステップの最大値 | ステップの RMS |
| --- | --- | --- | --- | --- |
| `gau_loose` | 2.5e-3 | 1.7e-3 | 1.0e-2 | 6.7e-3 |
| `gau` | 4.5e-4 | 3.0e-4 | 1.8e-3 | 1.2e-3 |
| `gau_tight` | 1.5e-5 | 1.0e-5 | 6.0e-5 | 4.0e-5 |
| `gau_vtight` | 2.0e-6 | 1.0e-6 | 6.0e-6 | 4.0e-6 |
| `baker` | 3.0e-4 | 2.0e-4 | 3.0e-4 | 2.0e-4 |

`baker` では、サイクル間のエネルギー変化が 1e-6 hartree 未満であることも求めます。`never` は収束を報告しないので、`--max-cycles` まで続きます。

> **補足:** YAML（`--config`）のキーは、YAML 設定の一覧の [`geom`](yaml-reference.md#geom)、[`opt`](yaml-reference.md#opt)、[`lbfgs`](yaml-reference.md#lbfgs)、[`rfo`](yaml-reference.md#rfo)、[`microiter`](yaml-reference.md#microiter) にあります。

---

## 使用上の注意点

* **マイクロイテレーション**: `--distance-restraint` があると通常の RFO を、`--embedcharge` では標準の最適化を使います。MM だけのステップには電荷埋め込みの力が入らないためです。MM の緩和は `--thresh` と同じプリセットで収束を判定します。別のプリセットは YAML の `microiter.micro_thresh` で指定できます。
* **`--mm-only` は `grad` だけ**: `--opt-mode hess` と組み合わせるとエラーで止まります（終了コード 2）。可動 MM 層と凍結 MM 層の区別はそのまま使います。
* **プラトーでの停止**: `--stop-plateau` は、力のノイズで力の収束条件に届かないときにサイクルを節約できますが、エネルギーが平坦であることは停留点の証拠になりません。実質的な上限は `--max-cycles` です。マイクロイテレーションの MM の緩和はこの判定では止めません。エネルギーの幅とサイクル数は `--stop-plateau-thresh` と `--stop-plateau-window` で指定できます。
* **拘束の強さ**: デフォルトの力の定数 300 eV·Å⁻² は距離を強く保ちます。実行例の 20 eV·Å⁻² は、目標の距離へゆるやかに導きます。
* **`--reject-uphill` は `hess` だけで有効**: `grad`（L-BFGS）では無視されます。
* **1 回に 1 構造**: `-i` には 1 つの構造を指定します。`.xyz` の入力には、原子の順と層を与える `--ref-pdb` が要ります。軌跡からは、使うフレームを先に `.xyz` に切り出してください。
* **`--flatten` はほぼ収束した構造に使う**: 虚振動が 25 本を超えると、`opt` は虚振動の除去を飛ばして警告を出します。先に構造を最適化してから、`--flatten` を付けて実行し直してください。
* **凍結原子があるときの剛体運動**: `--flatten` は、剛体運動を [`freq`](freq.md#凍結境界での剛体モード) と同じように扱い、`result.json` の `rigid_projection` に記録します。
* **凍結原子と拘束の全体**: 凍結する原子や拘束の選び方は、{ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>` を参照してください。
* **オプティマイザの状態の書き出し**: `--dump` を付け、YAML の `opt.dump_restart` に正の整数 N を指定すると、N サイクルごとに `restart_NNN.yaml` を書きます。mlmm-toolkit はこのファイルを読み戻さないので、止まった計算は final geometry から `opt` をやり直してください。
* **モデルと精度**: `--backend-model` でバックエンドのモデルを、`--precision` で精度を選べます。詳しくは自動生成のオプションの一覧（英語のみ）を参照してください。

---

## 関連ドキュメント

* [freq](freq.md) — 最適化した構造が極小点（n_imag = 0）かの確認
* [tsopt](tsopt.md) — 極小点ではなく TS（鞍点）の最適化
* [irc](irc.md) — TS から反応経路をたどり、最適化する端点を得る
* [define-layer](define-layer.md) — 最適化の前に ML 層と MM 層を B-factor に書き込む
* [all](all.md) — IRC の端点の最適化まで含む一連のワークフロー
* [トラブルシューティング](troubleshooting.md) — 実行が失敗したときの切り分け
* [YAML 設定の一覧](yaml-reference.md) — `opt`、`lbfgs`、`rfo`、`microiter` のすべての設定
* [用語集](glossary.md) — L-BFGS、RFO などの用語
* {ref}`終了コード <ja-exit-codes>` — 終了ステータスの意味
