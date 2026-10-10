# `path-search`（2 つ以上の構造を通る再帰的な MEP 探索）

`path-search` は、反応の順（R → … → P）に並べた **2 つ以上**の層付き酵素構造を通る、1 本につながった最小エネルギー経路（MEP）を、系全体に ML/MM 計算機を使って作ります。共有結合が変わる区間（セグメント）だけを再帰的に細かくし、各区間の経路は GSM（Growing String Method、デフォルト）または DMF（Direct Max Flux）で求めます。

---

## 主な用途

* **R → P を反応の区間に分ける**: 反応が 1 段か多段かが分からないときに、結合が変わる区間を見つける
* **中間体を通る多段階の経路**: R と P の間に既知の中間体を入れ、1 本につないだ経路を得る
* **区間ごとの TS 候補**: 結合が変わる区間ごとに HEI（最もエネルギーの高いイメージ）を `hei_seg_NN.xyz` に書き出し、[`tsopt`](tsopt.md) に渡す

ML 領域の計算にはデフォルトの **UMA**（Meta）を使います。残りの部分は `--parm7` の Amber 力場で計算します。2 つの端点だけで再帰的な精密化が要らない場合は、[`path-opt`](path-opt.md) のほうが簡単です。

---

## 基本的な実行例

### 1. 2 つの端点からの実行

反応物と生成物を 1 つの `-i` の後に並べ、ML 領域の電荷とスピン多重度を明示します。`reactant.pdb` と `product.pdb` は `real.parm7`（Amber のトポロジー）に対応する全系の構造で、`ml_region.pdb` はそのうちの ML 領域を指定します。

```bash
mlmm path-search -i reactant.pdb product.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --out-dir ./result_path_search
```

実行が終わったら `summary.log` の `[2] Segment-level MEP summary` の節を開くか、`summary.json` を読みます。`scientific_status` には、事前最適化とすべての経路の計算が収束すると `success`、そうでなければ `partial` か `failed` が入ります。`segments` には区間ごとの `index`・`tag`・`kind`・`bond_changes`・`converged`・`barrier_kcal` が並びます。結合が変わる区間には、TS 候補の `hei_seg_NN.xyz` も書き出されます。

### 2. 中間体を入れて多段階の経路を作る

構造を反応の順に 1 つの `-i` の後に並べます。隣り合う組ごとに探索し、得られた経路を 1 本につなぎます。

```bash
mlmm path-search -i R.pdb IM1.pdb IM2.pdb P.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q -1 -m 1 --out-dir ./result_path_search_multi
```

### 3. 事前最適化と重ね合わせを省いて軽く試す

入力の事前最適化と重ね合わせを省き、可動なイメージの数も減らします。入力がすでに最適化され、重ね合わせてある場合に使えます。

```bash
mlmm path-search -i reactant.pdb product.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --no-preopt --no-align --max-nodes 8 --out-dir ./result_path_search_fast
```

---

## 処理の仕組みと計算仕様

探索の前に、デフォルトでは各入力を事前最適化し（`--preopt`）、1 つ前の構造に重ね合わせます（`--align`）。固定原子は少しずつ位置を合わせ、そのあいだ残りの原子を緩和します。

1. **隣り合う組ごとの粗い MEP**:
隣り合う入力の組（A → B）ごとに、GSM または DMF で粗い MEP を作り、その HEI を求めます。
2. **HEI のまわりの緩和**:
`--refine-mode peak` では HEI の両隣のイメージ（HEI ± 1）を、`minima` では HEI から外側へたどった両側の最も近い極小を最適化し、近くの 2 つの極小 End1 と End2 を得ます。`--refine-mode` を省くと、GSM では `peak`、DMF では `minima` になります。
3. **キンク（kink）か反応区間か**:
End1 と End2 の間で共有結合が変わらなければ、その区間は *キンク* とみなし、線形補間のノードを数個入れて 1 つずつ最適化します。結合が変われば *反応区間* とみなし、End1 と End2 の間に新しく GSM または DMF の経路を作って障壁をはっきりさせます。
4. **結合の変化が残る区間だけを再帰**:
A → End1 と End2 → B の部分で結合の変化を調べ、変化が残る部分だけを、`--max-depth` の階層まで探索し直します。
5. **区間をつなぐ**:
得られた部分を 1 本の経路につなぎます。重なる端点は除き、隣り合う部分の端どうしで結合がまだ違えば、その間を新しい区間として探索します。それ以外のすき間は短い接続経路で埋めます。

結合が変わったかどうかは YAML の `bond` 節の閾値で決まり、判定の仕方は {ref}`scan <ja-section-bond>` と同じです。

---

## セグメントの判定

| 見えるもの | 意味 | 次の操作 |
| --- | --- | --- |
| 結合が変わる区間と、その `hei_seg_NN.xyz` | その段の TS 候補 | [`tsopt`](tsopt.md) で最適化して虚振動が 1 つかを確かめ、[`irc`](irc.md) を実行する |
| `tag` が `seg_NNN_maxdepth` または `seg_NNN_kinklimit` の区間 | 階層の上限に達した（`_maxdepth`）か、キンクの区間が続いた（`_kinklimit`）ため、そこより先は分けていない | 複数の段を含むことがある。上と同じように確かめるか、`--max-depth` を上げるか、中間体を入れる |
| `tag` が `_kink` で終わる区間しかない、または `HEI is at an endpoint` の警告 | 結合の変化が見つからないか、経路の端と端の間に頂点がない | 入力を見直すか、中間体を入れる（例 2） |

区間の分け方は結合距離による結合の判定にもとづく目安です。1 つの区間が 1 つの素過程であることも、TS をちょうど 1 つ含むことも保証しません。TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます。各 HEI を `tsopt`（n_imag = 1）と IRC で確かめてから、機構の 1 段として扱ってください。

---

## 主な出力ファイル

`--out-dir` に次のファイルを書き出します。

```text
result_path_search/
├─ mep_trj.xyz               # つないだ MEP 全体（コメント行にエネルギー）
├─ mep_trj.pdb               # 同じ経路の PDB
├─ mep_plot.png              # 経路に沿った ΔE（kcal/mol、反応物基準）
├─ energy_diagram_MEP.png    # MEP の状態エネルギー図（反応物基準）
├─ summary.json              # 区間ごとの障壁と分類の要約
├─ summary.log               # 同じ要約のテキスト
├─ mep_seg_NN_trj.xyz        # 反応区間 NN の経路（PDB は mep_seg_NN.pdb）
├─ hei_seg_NN.xyz            # 反応区間 NN の HEI（TS 候補。PDB は hei_seg_NN.pdb）
├─ hei_mode_seg_NN.*         # その HEI での反応方向の推定。all が tsopt に渡す
├─ align_refine/             # 入力の重ね合わせと緩和のファイル（--align）
├─ initNN_*_opt/             # 各入力の事前最適化（--preopt）
└─ seg_NNN_*/                # GSM・DMF の実行と HEI の両側の最適化ごとの作業ファイル
```

`summary.json` は常に書き出され、ほかのコマンドの `result.json` とは別の構造です。{ref}`path-search と all の summary.json <ja-summary-json-path-search-all>` を参照してください。`mep_seg_NN_*` と `hei_seg_NN.*` は、結合が変わる区間にだけ書き出します。NN は `summary.json` の区間の `index`（最終経路で 01 から数える）で、`seg_NNN` のタグやディレクトリの NNN は GSM・DMF の実行を 000 から数えた番号なので、両者は一致しません。{ref}`mmCIF 入力 <ja-mmcif-input>` と、PDB の列に収まらない大きな PDB 入力では、元の識別子を保った `.cif` も書き出します。`--no-convert-files` では `.xyz` だけを書き出します。

---

## 主な CLI オプション

ML/MM の計算コマンドに共通のオプションは {ref}`ML/MM の共通オプション <ja-mlmm-options>` に 1 か所でまとめてあります。下の表は `path-search` に固有のものだけです。

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | 反応の順に並べた 2 つ以上の構造（`.pdb`, `.cif`, `.mmcif`、または `--ref-pdb` を付けた `.xyz`）を 1 つの `-i` の後に並べる（ファイルごとに `-i` を繰り返してもよい） |
| `-q, --charge` | 整数 | `None` | ML 領域の電荷。`-l` を使う場合のほかは必須 |
| `-m, --multiplicity` | 整数 | `1` | スピン多重度（2S+1） |
| `-l, --ligand-charge` | 文字列 | `None` | リガンドの総電荷（例: `-1`）または残基名ごとの電荷（例: `'GPP:-3,SAM:1'`）。`-q` を省いたときに ML 領域の電荷を求めるのに使う（PDB/mmCIF 入力か `--ref-pdb` のとき） |
| `-b, --backend` | 文字列 | `uma` | ML 領域の計算バックエンド（`uma`, `orb`, `mace`, `aimnet2`, `dft`） |
| `-o, --out-dir` | パス | `./result_path_search/` | 出力先ディレクトリ |
| `--mep-mode` | `gsm` / `dmf` | `gsm` | 経路の手法（Growing String Method / Direct Max Flux） |
| `--dmf-backend` | `gpu` / `cpu` | `gpu` | DMF の計算バックエンド（`--mep-mode dmf` のときのみ）。CUDA 上の PyTorch / NumPy |
| `--refine-mode` | `peak` / `minima` | GSM では `peak`、DMF では `minima` | HEI のまわりの緩和のしかた（HEI ± 1 / 最も近い極小） |
| `--max-depth` | 整数 | `10` | 再帰的な分割の最大階層数。`0` で分割しない |
| `--max-nodes` | 整数 | `20` | 区間ごとの可動なイメージの数。区間のイメージは全部で `max_nodes + 2` 個 |
| `--preopt/--no-preopt` | フラグ | `True` | 探索の前に各入力を事前最適化 |
| `--align/--no-align` | フラグ | `True` | 探索の前に各入力を 1 つ前の構造に重ね合わせる |
| `--freeze-atoms` | 文字列 | `None` | {ref}`固定する原子 <ja-freeze-atoms-and-restraints>` の番号（1 始まり、カンマ区切り）。YAML の `geom.freeze_atoms` と固定 MM 層に加わる |
| `--climb/--no-climb` | フラグ | `True` | 反応区間で GSM のクライミングイメージ探索を行う。接続経路では常に行わない |

全オプションは [自動生成のオプションの一覧（英語のみ）](../reference/commands/path_search.md) を参照してください。

> **補足:** YAML（`--config`）では、`--max-depth` を指定しないときに `search.max_depth` が階層の上限になり、`search.kink_max_nodes`（デフォルト `3`）がキンクに入れるノードの数を、`bond.bond_factor`（デフォルト `1.20`）が結合の変化の判定に使う共有結合半径の倍率を決めます。すべてのキーは YAML 設定の一覧の [`search`](yaml-reference.md#search) と [`bond`](yaml-reference.md#bond) にあります。[`stopt`](yaml-reference.md#stopt) の `stopt.lbfgs`・`stopt.rfo` でも、`opt.lbfgs`・`opt.rfo` と同じように単一構造のオプティマイザを設定できます。

---

## 使用上の注意点

* **入力**: 2 つ以上の構造を、`--parm7` と同じ原子を同じ順に並べて与えてください。2 つ未満ではエラーで停止します。
* **XYZ 入力のテンプレート**: `--ref-pdb` には全系の PDB を入力ごとに 1 つ、`-i` と同じ順に並べます。テンプレートの無い `.xyz` 入力はエラーで停止します。
* **区間をつなぐ経路ではクライミングしない**: `--climb` は反応区間に効き、隣り合う部分をつなぐ短い経路は常にクライミングなしで計算します。
* **入力ファイルの保護**: 決まった名前の出力（`mep_trj.*`, `mep_plot.png`, `energy_diagram_MEP.png`, `summary.json`, `summary.log`）が入力ファイルを置き換えそうなときは、何も書き出す前に停止します。
* **YAML のオプティマイザ設定の矛盾**: 同じ YAML ファイルの中で、同じキーを `opt:` と実際に動くオプティマイザの節（`lbfgs:`・`opt.lbfgs:`・`stopt.lbfgs:`、または `rfo` の同等の節）とで別の値にすると、エラーで停止します。
* **DMF では固定原子が少し動く**: DMF は固定原子を調和拘束（k = 300 eV/Å²、YAML の `dmf.k_fix`）で保持するので、わずかにずれることがあります。[path-opt の「使用上の注意点」](path-opt.md#使用上の注意点) と {ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>` を参照してください。
* **DMF には `cyipopt` と `pydmf` が必要**: どちらも `mlmm-toolkit` と一緒にはインストールされません。`--mep-mode dmf` を使う前に、[path-opt の「使用上の注意点」](path-opt.md#使用上の注意点) の手順でインストールしてください。
* **複雑な機構**では、中間体、スキャンの設定、収束の閾値を調整する必要があることがあります。

---

## 関連ドキュメント

* [path-opt](path-opt.md) — 2 構造間の MEP を 1 回だけ求める
* [scan](scan.md) — 結合を段階的に動かして経路や TS 候補を作る
* [tsopt](tsopt.md) — 区間ごとの HEI から TS を最適化
* [ML 領域と層の組み方](model-setup.md) — 入力に使う全系の PDB、`real.parm7`、ML 領域を作る
* [all](all.md) — 一貫実行のワークフロー。`all --refine-path` で MEP の段に `path-search` を使います
* [YAML 設定の一覧](yaml-reference.md) — `search`・`bond`・`gs`・`dmf` の全設定
* [用語集](glossary.md) — MEP、GSM、DMF、HEI、キンクなどの用語
* [トラブルシューティング](troubleshooting.md) — 異常終了時の原因切り分けと対処法
* {ref}`終了コード <ja-exit-codes>` — 終了コードの意味
