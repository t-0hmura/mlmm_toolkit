# `dft`（DFT 一点計算）

`dft` サブコマンドは、1 つの ML/MM 構造の **ML 領域に対して GPU4PySCF（GPU）または PySCF（CPU）で DFT（密度汎関数理論）一点計算**を行い、MM のエネルギーと組み合わせて **ML(DFT)/MM の総エネルギー** `E_total = E_REAL_low + E_ML(DFT) - E_MODEL_low` を求めます。ML 領域の**原子電荷**も出力します。求めるのはエネルギーだけで、力は計算しません。

---

## 主な用途

* **ML/MM 構造での DFT エネルギー**: ML 領域を MLIP で計算して最適化した反応物（R）・遷移状態（TS）・生成物（P）の一点計算
* **電荷分布の把握**: ML 領域の原子ごとの電荷と、開殻系のスピン密度
* **タンパク質の静電場の効果**: `--embedcharge` で MM の点電荷を DFT のハミルトニアンに入れる

---

## 基本的な実行例

### 1. GPU での一点計算

中性の一重項の ML 領域について、GPU でエネルギーと電荷を計算します。`enzyme.pdb` は全系、`real.parm7` はその Amber トポロジー、`ml_region.pdb` は DFT で計算する原子です。`-q` と `-m` は ML 領域の電荷と多重度です。

```bash
mlmm dft -i enzyme.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 \
    --out-dir ./result_dft
```

端末に `E_DFT (Hartree): …` と `E_total ML(dft)/MM (Hartree): …` が出て、`result_dft/result.yaml` に `energy.converged: true` があれば成功です。

### 2. SCF を厳しくし、基底を大きくする

SCF（自己無撞着場）の収束を厳しくし、基底を大きくします。

```bash
mlmm dft -i enzyme.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 \
    --func-basis 'wb97m-v/def2-tzvpd' --scf-tol 1e-10 --scf-max-cycles 200 \
    --out-dir ./result_dft_tight
```

### 3. CPU だけで計算する

GPU の無いマシンでは、CPU の PySCF で計算できます。

```bash
mlmm dft -i enzyme.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 \
    --dft-engine cpu --out-dir ./result_dft_cpu
```

### 4. リガンドの電荷から ML 領域の電荷を求める

`-q` を省略して `-l` でリガンドの形式電荷を与えると、`dft` は ML 領域にあるアミノ酸残基とイオンの電荷を足して ML 領域の電荷を求め、その内訳を端末に表示します。

```bash
mlmm dft -i enzyme.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -l 'SAM:1,GPP:-3' -m 1 --out-dir ./result_dft_ligand
```

---

## 処理の仕組みと計算仕様

1. **ML 領域を組む**:
`-i` から全系を、`--parm7` から Amber トポロジーを、`--model-pdb`・`--model-indices`・入力の B-factor のどれかから ML 領域を読みます。XYZ 入力では、PDB/mmCIF のトポロジーを `--ref-pdb` で与えます。ML/MM の境界で切れる `--parm7` の結合をリンク水素でふさぎ、ML 領域をリンク水素なしとありの 2 通りで保存します。
2. **SCF**:
`--func-basis` で汎関数と基底を選びます。名前が `def2` で始まる基底には、対応する def2 の有効内殻ポテンシャル（ECP）を付けます。`--dft-engine` で GPU4PySCF（`gpu`、デフォルト）か PySCF（`cpu`）を選びます。閉殻は RKS、開殻は UKS で計算します。デフォルトで有効な低メモリモードでは、密度フィッティング（density fitting）を使わずに J と K を直接組み立てます。このとき、GPU で閉殻系を計算する場合は GPU4PySCF の低メモリ版 RKS を使います。`--no-dft-low-memory` では密度フィッティングを使います。`--embedcharge` を付けると、ML 領域から `--embedcharge-cutoff` 以内にある `--parm7` の MM 点電荷を DFT のハミルトニアンに入れます。
3. **ML(DFT)/MM のエネルギー**:
リンク水素を含む ML 領域の DFT エネルギーを、ONIOM の和の ML のエネルギーの代わりに使います。`E_REAL_low` と `E_MODEL_low` は、全系と ML 領域の MM エネルギーです。
4. **電荷と結果ファイル**:
SCF の後に ML 領域の Mulliken・meta-Löwdin・IAO（内在的原子軌道）の電荷とスピン密度を求め、エネルギー（Hartree と kcal/mol）とともに `result.yaml` に書き出します。失敗した解析の列は `null` になります。

---

## 主な出力ファイル

`--out-dir` に以下のファイルを書き出します。

```text
result_dft/
├─ ml_region_without_linkH.xyz   # 選んだままの ML 領域（リンク水素なし）
├─ ml_region_with_linkH.xyz      # リンク水素を付けた ML 領域（PySCF に渡した構造）
├─ ml_region_without_linkH.pdb   # 同じ構造の PDB（PDB 入力で --convert-files のとき）
├─ ml_region_with_linkH.pdb      # 同じ構造の PDB（PDB 入力で --convert-files のとき）
├─ result.yaml                   # エネルギー、収束、エンジン、原子ごとの電荷とスピン密度
├─ result.json                   # 機械可読な要約（--out-json 指定時）
└─ summary.json                  # result.json と同じ内容（--out-json 指定時）
```

* **`energy`**（`result.yaml`）: ML 領域の DFT エネルギーの `hartree`・`kcal_per_mol`、`converged`、使ったエンジンの `engine`（`gpu4pyscf(rks_lowmem)`・`gpu4pyscf`・`pyscf(cpu)`）・`used_gpu`・`used_lowmem`。
* **`mlmm_energy`**（`result.yaml`）: MM のエネルギー `E_real_low_hartree` と `E_model_low_hartree`、総エネルギー `E_total_ml_dft_mm_hartree`。kcal/mol の値もあります。
* **`charges [index, element, mulliken, lowdin, iao]`**: リンク水素を含む ML 領域の 1 原子 1 行の表で、`index` は 0 始まりです。端末にも同じ表が出ます。
* **`spin_densities [index, element, mulliken, lowdin, iao]`**: 同じ形の表です。`result.yaml` には常に書き出し、端末には開殻のときだけ表示します。
* **`result.json`**: エネルギー、`mulliken`・`lowdin`・`iao` の配列での電荷とスピン密度、電荷・多重度・汎関数・基底・SCF の設定を持ちます。[JSON 出力の一覧](json-output.md#dft) を参照してください。

---

## 主な CLI オプション

ML/MM の計算コマンドに共通のオプションは {ref}`ML/MM の共通オプション <ja-mlmm-options>` に 1 か所でまとめてあります。下の表は `dft` に固有のものだけです。

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | 全系の入力構造ファイル（`.pdb`, `.cif`、または `--ref-pdb` と組み合わせた `.xyz`） |
| `-q, --charge` | 整数 | `None` | ML 領域の電荷。`-l` か YAML の `calc.model_charge` が無ければ必須 |
| `-m, --multiplicity` | 整数 | `1` | ML 領域のスピン多重度（2S+1） |
| `-l, --ligand-charge` | 文字列 | `None` | 残基ごとの形式電荷（例: `'SAM:1,GPP:-3'`）またはリガンドの総電荷。`-q` を省いたときに ML 領域の電荷を求めるのに使用（PDB/mmCIF 入力または `--ref-pdb`） |
| `--func-basis` | 文字列 | `wb97m-v/def2-svp` | 汎関数と基底（`汎関数/基底` の形） |
| `--scf-tol` | 浮動小数点数 | `1e-9` | SCF の収束閾値（Hartree） |
| `--scf-max-cycles` | 整数 | `100` | SCF の最大反復回数 |
| `--dft-grid-level` | 整数 | `3` | 数値積分グリッドのレベル（PySCF の `grids.level`） |
| `--dft-engine` | `gpu` / `cpu` | `gpu` | GPU4PySCF か CPU の PySCF |
| `--dft-low-memory/--no-dft-low-memory` | フラグ | `True` | J と K を直接組み立てる。`--no-dft-low-memory` で密度フィッティングを使用 |
| `--scf-stepwise-grid/--no-scf-stepwise-grid` | フラグ | `True` | SCF をまず[粗いグリッド](dft-backend.md#使用上の注意点)で収束させ、その密度から最終のグリッドで収束させる |
| `--dft-nprocs` | 整数 | auto | PySCF の CPU スレッド数（スケジューラとホストから自動検出） |
| `--dft-memory` | 文字列 | auto | PySCF のホスト RAM の上限（例: `64GB`）。GPU メモリではない |
| `--embedcharge/--no-embedcharge` | フラグ | `False` | MM の点電荷を DFT のハミルトニアンに入れる |
| `--embedcharge-cutoff` | 浮動小数点数 | `12.0` | 点電荷を入れる MM 原子の、ML 領域からの距離（Å） |
| `--convert-files/--no-convert-files` | フラグ | `True` | ML 領域を PDB でも書き出す（PDB 入力のときだけ） |
| `-o, --out-dir` | パス | `./result_dft/` | 出力先ディレクトリ |

全オプションは [自動生成のオプションの一覧（英語のみ）](../reference/commands/dft.md) を参照してください。

> **補足:** YAML（`--config`）では、{ref}`dft <ja-dft-section>` の節で同じ設定を指定できます。`dft.pyscf` は PySCF のオブジェクトに属性を名前で渡し、SCF が収束しにくいときは `pyscf: {mf: {level_shift: 0.2}}` のように使えます。ML 領域の電荷と多重度は `calc.model_charge`・`calc.model_mult` に書きます。優先されるのは `-q`・`-l`・`-m`、YAML の順です。

---

## 使用上の注意点

* **必要なもの**: DFT 用の追加パッケージが必要で、PyTorch の wheel が `cu130`・`cu132` なら `pip install "mlmm-toolkit[dft]"`、`cu126` なら `pip install "mlmm-toolkit[dft-cuda12]"` で入れます。
* **基底のコスト**: `def2-tzvpd` は `def2-svp` よりはるかに重い計算です。原子数や GPU メモリの決まった上限は無く、コストは基底関数の数・元素・汎関数・グリッド（`--dft-grid-level`）・GPU で決まります。まず代表構造を 1 つ計算し、メモリの最大使用量を確かめてください。足りないときは、基底を小さくするか、メモリの大きい GPU を使ってください。
* **GPU**: GPU4PySCF が動かないとき、`dft` は CPU のエンジンを勧めるエラーで止まり、自動では CPU に切り替えません。Blackwell（RTX 50xx）のような新しい世代の GPU では、メモリ不足や未対応カーネルのエラーがメモリ量ではなく GPU4PySCF と CuPy の版から来ることがあるので、まず版とトレースバックを確かめてください。
* **CPU**: `--dft-engine cpu` は GPU を必要としません。実用になる ML 領域の大きさは手法とマシンで変わるため、代表構造の一点計算で時間を測ってください。
* **x86 以外のマシン**: GPU4PySCF のビルド済みホイールが対応しないことがあります。その場合はソースからビルドしてください（https://github.com/pyscf/gpu4pyscf）。
* **ECP**: 名前が `def2` で始まる基底には、元素によらず def2 の ECP を付け、端末に `[dft] Using ECP: …` と表示します。
* **IAO 解析**は難しい系で失敗することがあり、そのとき `result.yaml` の該当列は `null` になります。
* **SCF が収束しないとき**: `dft` は `WARNING: SCF did not converge.` を表示し、`converged: false` として `result.yaml`（`--out-json` なら `result.json` も）を書いた上で、終了コード 1 で終わります。低メモリモードでは、メモリに余裕があれば密度フィッティング `--no-dft-low-memory`（別名 `--no-lowmem`）で再実行するよう提案します。
* **多重度**: 1 未満は受け付けません。
* **前回の結果**: 実行の最初に、出力ディレクトリに残っている `result.yaml`・`result.json`・`summary.json` と 4 つの `ml_region_*` ファイルを削除します。
* **終了コード**: {ref}`終了コード <ja-exit-codes>`を参照してください。

---

## 関連ドキュメント

* [MLIP の TS を DFT で確かめる](dft-backend.md) — ワークフローでの `-b dft` と `--dft`、DFT の設定、GPU メモリ
* [sp](sp.md) — `-b dft` を含む任意のバックエンドでの ML/MM の一点エネルギーと力
* [all](all.md) — 全工程のワークフロー。`--dft` で R・TS・P に DFT 一点計算を追加
* [MLIP バックエンド](backends.md) — バックエンドの選び方
* [トラブルシューティング](troubleshooting.md) — 実行に失敗したときの対処
