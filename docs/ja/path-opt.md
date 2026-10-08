# `path-opt`（2 構造間の MEP 探索）

`path-opt` は、反応物と生成物の 2 つの層付き酵素構造の間の最小エネルギー経路（MEP）を、系全体に ML/MM 計算機を使って GSM（Growing String Method、デフォルト）または DMF（Direct Max Flux）で 1 回だけ求めます。最もエネルギーの高いイメージ（HEI）を TS 候補として書き出します。

---

## 主な用途

* **R と P からの MEP の初案**: 2 つの端点構造から、再帰的な精密化なしで経路とエネルギープロファイルを得る
* **`tsopt` に渡す TS 候補**: `hei.pdb`（または `hei.xyz`）を [`tsopt`](tsopt.md) の初期構造にする
* **GSM と DMF の比較**: 同じ 2 構造を `--mep-mode gsm` と `--mep-mode dmf` で計算し、経路を見比べる

ML 領域の計算にはデフォルトの **UMA**（Meta）を使います。残りの部分は `--parm7` の Amber 力場で計算します。2 つ以上の構造から反応領域だけを自動で精密化したい場合は、[`path-search`](path-search.md) を使ってください。

---

## 基本的な実行例

### 1. 2 つの端点からの実行

反応物と生成物を 1 つの `-i` の後に並べ、ML 領域の電荷とスピン多重度を明示します。`reactant.pdb` と `product.pdb` は `real.parm7`（Amber のトポロジー）に対応する全系の構造で、`ml_region.pdb` はそのうちの ML 領域を指定します。

```bash
mlmm path-opt -i reactant.pdb product.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --out-json --out-dir ./result_path_opt
```

端末に `[write] Wrote '…/hei.xyz'` で始まる行が出れば、TS 候補が書き出されています。`result.json` の `scientific_status` には、求めた段（端点の事前最適化と MEP）がすべて収束すると `success`、そうでなければ `partial` か `failed` が入ります。`barrier_kcal` は最初のイメージを基準にした HEI のエネルギー、`hei_index` は経路上の HEI の位置です（下の「HEI の判定」を参照）。

### 2. 端点の事前最適化の上限を変える

両端点はデフォルトで事前最適化され、`--preopt-max-cycles` で各回のサイクル数の上限を変えられます。端点が最適化済みなら `--no-preopt` で省けます。

```bash
mlmm path-opt -i reactant.pdb product.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --preopt-max-cycles 20000 --out-dir ./result_path_opt_preopt
```

### 3. GSM の代わりに DMF を使う

DMF には `cyipopt` と `pydmf` が必要です（「使用上の注意点」を参照）。この例では可動なイメージの数も減らしています。

```bash
mlmm path-opt -i reactant.pdb product.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --mep-mode dmf --max-nodes 12 --out-dir ./result_path_opt_dmf
```

### 4. 短時間の確認（クライミングなし、少ないノード）

クライミングイメージ探索を省き、可動なイメージの数も減らして、経路の形をすばやく確かめます。

```bash
mlmm path-opt -i reactant.pdb product.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --no-climb --max-nodes 8 --out-dir ./result_path_opt_quick
```

---

## 処理の仕組みと計算仕様

1. **端点の準備**:
各端点を事前最適化します（デフォルトは L-BFGS）。続いて凍結原子（無ければ全原子）の Kabsch フィットで生成物を反応物に剛体で重ね合わせ、凍結原子を少しずつ反応物側の位置へ動かしながら残りの原子を緩和します。凍結原子は Frozen-MM 層と `--freeze-atoms` で指定した原子を合わせたものです。
2. **経路の成長と精密化**:
GSM は 2 つの端点の間に `--max-nodes` 個の可動なイメージのストリングを成長させ、`--thresh-gsm` まで最適化します。`--climb`（デフォルト有効）では、続くクライミングイメージ探索で最も高いイメージを鞍点へ押し上げます。DMF は補間で作った経路を IPOPT（内点法の最適化ソルバー）で `--dmf-tol` まで最適化します。
3. **HEI の書き出し**:
最終経路でエネルギーが最も高いイメージを HEI とし、コメント行にエネルギーを付けて `hei.xyz` に書き出します。

---

## HEI の判定

| HEI の位置 | 意味 | 次の操作 |
| --- | --- | --- |
| 経路の内側（`hei_index` が `1` 〜 `n_images − 2`） | TS 候補 | [`tsopt`](tsopt.md) で最適化し、[`irc`](irc.md) を実行する |
| 端点（`hei_index` が `0` または `n_images − 1`） | TS 候補ではない（端点の間に、高いほうの端点より上にあるイメージがない） | 端点を見直すか、別の方法で候補を作る（{ref}`TS が取れないとき <ja-ts-search-fails>` を参照） |

HEI は近似的な経路の頂点であり、TS そのものではありません。TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます。`tsopt` で n_imag = 1 になり、その TS からの IRC が狙った R と P に着けば、HEI から TS が得られたことになります。

---

## 主な出力ファイル

実行完了後、`--out-dir`（デフォルト: `./result_path_opt/`）内に以下のファイル群が生成されます。

```text
result_path_opt/
├─ final_geometries_trj.xyz   # 最終経路の全イメージ（コメント行にエネルギー）
├─ final_geometries.pdb       # 同じ経路の PDB
├─ hei.xyz                    # HEI（TS 候補）。コメント行にエネルギー
├─ hei.pdb                    # 同じ HEI の PDB
├─ preopt/endNN/              # 端点の事前最適化のファイル（--preopt）
├─ align_refine/              # 端点の重ね合わせと緩和のファイル
├─ dmf_initial_trj.xyz        # 補間で作った初期経路（DMF のときのみ）
├─ result.json                # 結果の要約（--out-json）
└─ summary.json               # result.json と同じ内容（--out-json）
```

経路は `final_geometries_trj.xyz` を開いて確かめてください。PDB ファイルの B-factor の列には層が入っている（ML 0、Movable-MM 10、Frozen-MM 20）ので、`tsopt` には `hei.pdb` を同じ `--parm7` と `--model-pdb` とともに渡してください。`hei.xyz` を渡すときは `--ref-pdb` を付けてください。mmCIF 入力と、PDB の列に収まらない大きな PDB 入力では、元の識別子を保った `.cif` も書き出します（{ref}`mmCIF 入力 <ja-mmcif-input>` を参照）。`--no-convert-files` では `.xyz` だけを書き出します。DMF では IPOPT のログ `dmf_fbenm_ipopt.out` と `dmf_ipopt.out` も書き出します。`--dump` を付けると、オプティマイザの軌跡も残します。

端末には MEP の進行状況がサイクルごとに、所要時間とともに出ます。

---

## 主な CLI オプション

ML/MM の計算コマンドに共通のオプションは {ref}`ML/MM の共通オプション <ja-mlmm-options>` に 1 か所でまとめてあります。下の表は `path-opt` に固有のものだけです。

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス 2 つ | （必須） | 反応物と生成物をこの順に 1 つの `-i` の後に並べる（`.pdb`, `.cif`, `.mmcif`、または `--ref-pdb` を付けた `.xyz`） |
| `-q, --charge` | 整数 | `None` | ML 領域の電荷。`-l` を使う場合のほかは必須 |
| `-m, --multiplicity` | 整数 | `1` | スピン多重度（2S+1） |
| `-l, --ligand-charge` | 文字列 | `None` | リガンドの総電荷（例: `-1`）または残基名ごとの電荷（例: `'GPP:-3,SAM:1'`）。`-q` を省いたときに ML 領域の電荷を求めるのに使う（PDB/mmCIF 入力か `--ref-pdb` のとき） |
| `-b, --backend` | 文字列 | `uma` | ML 領域の計算バックエンド（`uma`, `orb`, `mace`, `aimnet2`, `dft`） |
| `-o, --out-dir` | パス | `./result_path_opt/` | 出力先ディレクトリ |
| `--mep-mode` | `gsm` / `dmf` | `gsm` | 経路の手法（Growing String Method / Direct Max Flux） |
| `--dmf-backend` | `gpu` / `cpu` | `gpu` | DMF の計算バックエンド（`--mep-mode dmf` のときのみ）。CUDA 上の PyTorch / NumPy |
| `--max-nodes` | 整数 | `20` | 端点の間の可動なイメージの数。経路のイメージは全部で `max_nodes + 2` 個 |
| `--preopt/--no-preopt` | フラグ | `True` | 重ね合わせの前に各端点を事前最適化 |
| `--preopt-max-cycles` | 整数 | `100000` | 端点の事前最適化 1 回あたりの最大サイクル数 |
| `--opt-mode` | `grad` / `hess` | `grad` | 端点の事前最適化の最適化法: L-BFGS / RFO |
| `--thresh-gsm` | プリセット | `gau_loose` | GSM のストリングの収束条件（`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`） |
| `--dmf-tol` | 文字列 | `tight` | DMF の経路の IPOPT の許容値: `tight`（0.04）、`middle`（0.10）、`loose`（0.20）、または正の数。別名 `--thresh-dmf` |
| `--fix-ends/--no-fix-ends` | フラグ | `True` | GSM のストリングの最適化の間、端点を固定する（DMF では使わない） |
| `--climb/--no-climb` | フラグ | `True` | 経路の成長後に GSM のクライミングイメージ探索を行う（DMF では使わない） |
| `--freeze-atoms` | 文字列 | `None` | すべてのイメージで凍結する原子の番号（1 始まり、カンマ区切り）。YAML の `geom.freeze_atoms` と Frozen-MM 層に加わる（{ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>` を参照） |
| `--out-json/--no-out-json` | フラグ | `False` | 結果の要約を `result.json` に出力（[JSON 出力リファレンス](json-output.md)） |

全オプションの一覧は [自動生成 CLI リファレンス](../reference/commands/path_opt.md) を参照してください。

> **補足:** YAML（`--config`）では、[`gs`](yaml-reference.md#gs) の節で GSM のストリングを、[`dmf`](yaml-reference.md#dmf) の節で DMF の経路を、[`stopt`](yaml-reference.md#stopt) の節でストリングのオプティマイザを設定できます。`stopt.lbfgs`・`stopt.rfo` でも、`opt.lbfgs`・`opt.rfo` と同じように端点のオプティマイザを設定できます。

---

## 使用上の注意点

* **DMF では凍結原子が少し動く**: DMF は凍結原子を固定せず、調和拘束（k = 300 eV/Å²、YAML の `dmf.k_fix`）で保持するので、参照位置からわずかにずれることがあります。GSM では固定されたままです。{ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>` を参照してください。
* **DMF には `cyipopt` と `pydmf` が必要**: どちらも `mlmm-toolkit` と一緒にはインストールされません。`--mep-mode dmf` を使う前に、`conda install -c conda-forge cyipopt -y` を実行し、続けて GPU バックエンドなら `pip install 'pydmf[torch]>=1.2'`、CPU バックエンドなら `pip install 'pydmf>=1.2'` でインストールしてください。デフォルトの `--dmf-backend gpu` は CUDA が使えないとエラーで停止します。その場合と GPU のメモリ不足のときは `--dmf-backend cpu` を指定してください。
* **DMF で使われないオプション**: `--climb`・`--fix-ends` は、指定しても DMF では使われません。
* **XYZ の端点のテンプレートは 1 つ**: `--ref-pdb` には全系の PDB を 1 つだけ指定し、両方の `.xyz` 端点に使います。
* **YAML のオプティマイザ設定の矛盾**: 同じ YAML ファイルの中で、同じキーを `opt:` と実際に動くオプティマイザの節（`lbfgs:`・`opt.lbfgs:`・`stopt.lbfgs:`、または `rfo` の同等の節）とで別の値にすると、エラーで停止します。
* **オプションの優先度**: **「デフォルト値 < 設定 YAML < コマンドライン引数」** の順です（{ref}`設定の優先順位 <ja-configuration-precedence>` を参照）。
* **終了コード**: {ref}`終了コード <ja-exit-codes>` を参照してください。

---

## 関連ドキュメント

* [path-search](path-search.md) — 2 つ以上の構造を通り、結合が変わる区間を精密化する MEP 探索
* [tsopt](tsopt.md) — HEI から TS を最適化
* [irc](irc.md) — TS が狙った R と P につながるかを確認
* [all](all.md) — 一貫実行のワークフロー。MEP の段は `path-opt` を使い、`--refine-path` で `path-search` に切り替えられます（[all の「主な CLI オプション」](all.md#主な-cli-オプション)）
* [YAML リファレンス](yaml-reference.md) — `gs`・`dmf`・`stopt` の全設定
* [用語集](glossary.md) — MEP、GSM、DMF、HEI などの用語
* [トラブルシューティング](troubleshooting.md) — 異常終了時の原因切り分けと対処法
* {ref}`終了コード <ja-exit-codes>` — 終了コードの意味
