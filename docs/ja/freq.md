# `freq`（振動解析・熱化学計算）

`freq` サブコマンドは、層を定義した ML/MM の酵素モデルについて、**調和振動数解析**と、ゼロ点振動エネルギー（ZPE）・エンタルピー・ギブズ自由エネルギーなどの**熱化学補正量**を計算します。

---

## 主な用途

* **停留点の検証**: 最適化された構造が極小点（n_imag = 0）か、遷移状態（TS、n_imag = 1）かを確認
* **熱化学量の算出**: QRRHO（準剛体ローター・調和振動子）モデルに基づき、自由エネルギーなどの熱力学諸量を算出
* **振動モードの可視化**: 虚振動（反応座標方向の運動）や特定モードの原子変位アニメーションを、酵素全体の `.xyz` / `.pdb` / `.cif` として出力

ML 領域の計算バックエンドのデフォルトは、Meta が公開した学習済みの[機械学習原子間ポテンシャル（MLIP）](backends.md) の **UMA** です。`-b/--backend` で **ORB**、**MACE**、**AIMNet2**、DFT も選択可能です。MM 領域には `--parm7` の Amber パラメータを使います。

---

## 基本的な実行例

以下の例では、`pocket.pdb` が `real.parm7` に対応する全系の構造、`ml_region.pdb` が [ML 領域](model-setup.md)の定義です。

### 1. 最小構成での実行（電荷と多重度を明示）

```bash
mlmm freq -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --out-dir ./result_freq
```

端末の熱化学の要約の `Number of Imaginary Freq = N` の行に n_imag が出ます。

### 2. 凍結原子の追加と熱化学の詳細ファイルの出力

`--freeze-atoms` で凍結 MM 層に加えて原子を凍結し、`--dump` で熱化学解析の詳細ファイル `thermoanalysis.yaml` も出力します。

```bash
mlmm freq -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --freeze-atoms "1,3,5,7" --dump --out-dir ./result_freq_phva
```

### 3. 解析的 Hessian（Analytical）の指定

有限差分の変位幅による誤差を避けたい場合に指定します。GPU メモリを多く使うので、先に対象の系で試してください。

```bash
mlmm freq -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --hessian-calc-mode Analytical --out-dir ./result_freq_analytical
```

---

## 処理の仕組みと計算仕様

1. **層の読み込みと凍結（PHVA）**:
ML 領域を `--model-pdb` から、可動 MM 層と凍結 MM 層を入力 PDB の B-factor から読み込みます。振動解析に入れる原子は `--active-dof-mode` で選び、デフォルトの `partial` は ML 領域と可動 MM 原子です。凍結 MM 層と `--freeze-atoms` で指定した原子は固定したままです。
2. **Hessian の計算**:
`--hessian-calc-mode` で、ML 領域の Hessian を `FiniteDifference`（有限差分、デフォルト）または `Analytical`（解析的）で計算します。`--hess-device` は、計算した Hessian を置いて対角化するデバイスを選びます。
3. **熱化学ポリシー（QRRHO）**:
振動数に、ローター閾値 100 cm⁻¹ の QRRHO 法を適用し、ギブズ自由エネルギー補正（`G_corr`）を算出します。端末の熱化学の要約に、`G_corr` とギブズエネルギー `E + G_corr = G`（E は電子エネルギー）が Hartree 単位で出ます。`--dump` では `thermoanalysis.yaml` に、`--out-json` では `result.json` の `thermochemistry` にも書き出します。凍結原子があるときは、振動の項に PHVA の振動数を使います。G には、正の振動数による振動の項のほか、構造全体の並進と回転の項を常に含めます。点群と回転対称数は構造から自動で判定し、YAML の `thermo.symmetry_number` で上書きできます。
4. **振動モードの書き出し**:
低振動数側または虚振動側から順に、原子振動アニメーション（軌跡ファイル）を出力します。出力本数は `--max-write`、並び順は `--sort` で指定します。

### 凍結境界での剛体モード

凍結原子が無いときは、振動ではない剛体運動（並進 3 つと回転 3 つ）の 6 つを取り除いてから振動数を出します。凍結原子があるときは、すべての凍結原子をその場に残す剛体運動だけを取り除きます。凍結 MM 層を持つふつうの ML/MM モデルのように、一直線に並ばない凍結原子が 3 つ以上あれば、取り除く剛体運動は 0 で、可動原子の振動モードはすべて残ります。凍結原子が 1 つなら 3 つ（その原子のまわりの回転）、2 つなら 1 つ（2 原子を結ぶ軸まわりの回転）を取り除きます。

`irc`、`tsopt` の TS の振動数の確認と Dimer の向きの計算、`opt`・`tsopt` の `--flatten`（余分な虚振動を消す処理）も、剛体運動を同じように扱います。`--out-json` を付けると、取り除いた剛体運動の数と使った Hessian が `result.json` の `rigid_projection` に記録されます（[JSON 出力の一覧](json-output.md#剛体モードの射影の記録)）。

---

## 振動数の出力と判定基準

`frequencies_cm-1.txt` と JSON の記録での扱いは次のとおりです。

| 項目 / パラメータ | 基準値・挙動 | 説明 |
| --- | --- | --- |
| **虚振動の表示** | 負の値（ν < 0 cm⁻¹） | 虚振動モードは負の振動数として表記されます。 |
| **虚振動の判定閾値** | ν < −5.00 cm⁻¹ | このモードを虚振動（n_imag。JSON の欄は `n_imaginary`）として数えます。閾値は YAML の `freq.zero_cutoff_cm`（デフォルト `5.0`）です。 |
| **微小な負のモード** | −5.00 ≤ ν < 0 cm⁻¹ | 数値誤差による微小な負モードは n_imag に数えません。`n_negative_modes` は、これも含めたすべての負の振動数の数です。 |
| **熱化学計算での扱い** | 反転・フロア処理なし | 虚振動数の反転や微小正振動数の底上げ（floor）は行いません。QRRHO は正の振動数のモードだけを使うので、虚振動モードは ZPE と G に入りません。 |

---

## 主な出力ファイル

実行が終わると、`--out-dir` に次のファイルができます。

```text
result_freq/
├─ frequencies_cm-1.txt          # 全振動数のリスト（cm⁻¹）
├─ mode_0001_-385.20cm-1_trj.xyz # 振動モードごとの変位アニメーション（XYZ）
├─ mode_0001_-385.20cm-1.pdb     # PDB 形式の変位アニメーション
├─ mode_0001_-385.20cm-1.cif     # mmCIF 形式の変位アニメーション（mmCIF 入力または大きな PDB 入力のとき）
├─ thermoanalysis.yaml           # 詳細熱化学ログ（--dump 指定時）
└─ result.json                   # 結果の要約（--out-json 指定時）
```

* **極小か TS か**: 端末の熱化学の要約の `Number of Imaginary Freq = N` の行に n_imag が出ます。`freq` はこの値を判定しないので、`result.json` の `scientific_status` は n_imag によらず `success` です。極小は n_imag = 0 です。TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます。`frequencies_cm-1.txt` の先頭に有意な虚振動（反応座標に沿った負の値）が**1 つだけ**あり、2 番目以降が正の値（または許容誤差内）であれば TS で、次は [`irc`](irc.md) に進みます。極小のはずの構造に虚振動が出たら [`opt`](opt.md) の `--flatten` で最適化し直し、TS に虚振動が無いか 2 つ以上あるときは {ref}`TS が取れないとき <ja-ts-search-fails>` を見てください。
* **可視化**: 生成された `mode_*_trj.xyz` や `.pdb` を PyMOL や VMD などの分子ビューアで開くと、原子が振動するアニメーションを確認できます。

---

## 主な CLI オプション

ML/MM の計算コマンドに共通のオプションは {ref}`ML/MM の共通オプション <ja-mlmm-options>` に 1 か所でまとめてあります。下の表は `freq` に固有のものだけです。

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | `--parm7` に対応する全系の構造（`.pdb`, `.cif`, `.mmcif`、または `--ref-pdb` と組み合わせた `.xyz`） |
| `-q, --charge` | 整数 | `None` | ML 領域の電荷。`-l` を使う場合のほかは必須 |
| `-m, --multiplicity` | 整数 | `1` | ML 領域のスピン多重度（2S+1） |
| `-l, --ligand-charge` | 文字列 | `None` | 未知のリガンドの総電荷、または残基名ごとの電荷（例: `'GPP:-3,SAM:1'`）。`-q` を省いたときに ML 領域の電荷を求めるのに使用（PDB 入力または `--ref-pdb`） |
| `-o, --out-dir` | パス | `./result_freq/` | 出力先ディレクトリ |
| `-b, --backend` | 文字列 | `uma` | ML 領域の計算バックエンド（`uma`, `orb`, `mace`, `aimnet2`, `dft`） |
| `--hessian-calc-mode` | `FiniteDifference` / `Analytical` | `FiniteDifference` | ML 領域の Hessian の計算法（有限差分 / 解析的） |
| `--hess-device` | `auto` / `cuda` / `cpu` | `auto` | 計算した Hessian を置いて対角化するデバイス。`cpu` では対角化の前に GPU から移す |
| `--active-dof-mode` | `all` / `ml-only` / `partial` / `unfrozen` | `partial` | 解析に入れる原子: 全原子 / ML 領域だけ / ML 領域と可動 MM 原子 / 凍結層以外のすべての原子 |
| `--freeze-atoms` | 文字列 | `None` | 追加で凍結する原子インデックス（1 始まり、カンマ区切り: 例 `'1,3,5'`） |
| `--read-hess` | パス | `None` | Hessian を計算せず、`.npy` ファイル（`freq`・`tsopt` の `--dump-hess` で保存したものなど）から読む |
| `--dump-hess` | パス | `None` | Hessian を `.npy` ファイルに保存する（`freq`・`tsopt`・`irc` の `--read-hess` で使える） |
| `--max-write` | 整数 | `10` | アニメーション出力する振動モードの最大数 |
| `--sort` | `value` / `abs` | `value` | モードの並び順（値順 / 絶対値順） |
| `--temperature` | 浮動小数点数 | `298.15` | 熱化学計算の温度（K） |
| `--pressure` | 浮動小数点数 | `1.0` | 熱化学計算の圧力（atm） |
| `--dump/--no-dump` | フラグ | `False` | 詳細熱化学ファイル（`thermoanalysis.yaml`）を出力 |
| `--out-json/--no-out-json` | フラグ | `False` | 結果要約を `result.json` に出力 |

全オプションは [自動生成のオプションの一覧（英語のみ）](../reference/commands/freq.md) を参照してください。

> **補足:** YAML（`--config`）では、{ref}`freq <ja-freq-section>` の節で虚振動の判定閾値 `zero_cutoff_cm` や書き出すモードの本数・振幅を、[`thermo`](yaml-reference.md#thermo) の節で温度と圧力を設定できます。

---

## 使用上の注意点

* **`tsopt` との使い分け**: `tsopt` コマンドには内部で虚振動数チェックが含まれています。別途 `freq` を単独実行するのは、主に詳細な熱化学量（ZPE、ギブズエネルギー）の取得や、変位ベクトルの可視化ファイルを出力したい場合です。
* **全原子凍結の禁止**: すべての原子を凍結指定すると、可動な振動自由度（DOF）が存在しなくなるためエラーで停止します。
* **Hessian モードの優先度**: `--hessian-calc-mode` の設定は **「デフォルト値 < 設定 YAML < コマンドライン引数」** の順で優先されます。
* **解析的 Hessian と `--uma-workers`**: UMA で `--uma-workers`（MLIP の並列ワーカー数）を 2 以上にすると、`--hessian-calc-mode Analytical` は使えずエラーで停止します。[解析的 Hessian](backends.md#ワーカーと-hessian-の計算方式)には `--uma-workers 1` を指定してください。
* **`all --thermo` の熱化学ファイル**: `all` は `thermoanalysis.yaml` からギブズエネルギーの図を作るので、`--thermo` で実行すると、その `freq` の段は `--no-dump` を指定してもこのファイルを書き出します。
* **`--read-hess`・`--dump-hess` のファイル**: `numpy.save` で書いた配列 1 つで、中身は質量重み付けなしの Cartesian の Hessian（Hartree/bohr²）です。原子は入力の順で、全原子の 3N×3N か、`--active-dof-mode` で選んだ原子の分だけを持ちます。`--read-hess` は、正方・有限・対称で、この 2 つの大きさのどちらかであることしか確かめないので、同じ構造・電荷・多重度・層・計算設定で求めた Hessian を渡してください。

---

## 関連ドキュメント

* [opt](opt.md) — 極小点への構造最適化
* [tsopt](tsopt.md) — 遷移状態（TS）の構造最適化
* [irc](irc.md) — 遷移状態からの固有反応座標（IRC）追跡
* [dft](dft.md) — 構造に対する高精度な DFT 一点エネルギー計算
* [all](all.md) — モデル作成・経路探索・TS 最適化・振動解析を一貫実行するワークフロー
* [YAML 設定の一覧](yaml-reference.md) — 設定ファイルの記法
* [トラブルシューティング](troubleshooting.md) — 異常終了時の原因切り分けと対処法
* {ref}`終了コード <ja-exit-codes>` — 終了コードの意味
