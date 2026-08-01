# `tsopt`

`mlmm tsopt` はレイヤー分けした酵素 PDB の遷移状態*候補*を一次サドル点まで精密化します。単独の TS（遷移状態）推測構造でも、[`path-search`](path-search.md) が抽出する最高エネルギー像（HEI）でも実行できます。

オプティマイザは 2 種類あり、`--opt-mode` で選択します。

- **RS-I-RFO**（`--opt-mode hess`）はデフォルトで、Hessian 計算のコストを許容できる場合の保守的な選択肢です。マイクロイテレーション（`--microiter`、デフォルト有効）が ML 1 ステップ RS-I-RFO と MM L-BFGS 緩和を交互に実行します。
- **Hessian-Guided Dimer**（`--opt-mode grad`）はより軽量な代替で、低コストな探索や複数の TS 推測構造からの素早い反復に向きます。`--ml-only-hessian-dimer` を付けると ML 領域のみの Hessian を Dimer 方向決定に使用できます（高速）。

`tsopt` は、YAML 上書き後も鞍点探索を担う RFO 系および Dimer
optimizer の `reject_uphill` を常に `false` に固定します。TS 探索では
反応モードに沿った物理エネルギー上昇を許す必要があるためです。
`--reject-uphill/--no-reject-uphill` は最小値最適化（`opt` と `all` の
IRC 後エンドポイント再最適化）だけに適用されます。マイクロイテレーション
内部の MM-only 緩和は最小化の部分問題なので、この区別を維持します。

収束後は `--flatten` の余剰虚モード除去ループが質量重み付け変位で余分な負のモードを整理します。検証済み TS は**正確に 1 つ**の虚振動数を示すべきで、必ず [`freq`](freq.md) / [`irc`](irc.md) でモードと結合性を確認してください。

## まず TS 候補を構築する

`tsopt` は*候補*を精密化するもので、ゼロから探索するわけではありません。手元にある情報に合った経路を選び、その結果を `tsopt → irc → freq`（または `mlmm all --tsopt`）に渡します。

| 経路 | サブコマンド | 内容 | 使う場面 |
| --- | --- | --- | --- |
| (a) MEP / 経路探索 | [`path-search`](path-search.md)（1 セグメントなら [`path-opt`](path-opt.md)） | 再帰的 GSM/DMF 最小エネルギー経路探索。端点間で TS を挟み込み、セグメント間のギャップを橋渡しし、セグメントごとに TS を 1 つ出力。 | 反応物（任意で生成物や中間体）があり、経路を*探索*させたいとき。 |
| (b) 距離拘束による構築 | [`scan`](scan.md) | 反応ペアごとに調和拘束 `E = ½·k·(r_ij − target)²` を加え、残りを L-BFGS で緩和しながら反応距離をバリアへ向けて駆動。 | 使える第 2 端点も TS 推測構造もないとき — 反応結合を直接駆動する。 |

```bash
# 経路 (a): 経路を探索し、その最高エネルギー像を精密化
mlmm path-search -i r.pdb p.pdb --parm enzyme.parm7 -l 'LIG:Q' -o result_mep

# 経路 (b): 反応距離を駆動して TS 候補を構築
mlmm scan -i r.pdb --parm enzyme.parm7 -l 'LIG:Q' \
    --scan-lists '[(1,5,1.40)]' -o result_scan
```

```{note}
`opt --restraint` フラグは存在しません。素の [`opt`](opt.md) は*拘束なし*のオプティマイザです。TS 候補の拘束付き構築は [`scan`](scan.md)（距離を駆動）または [`path-search`](path-search.md)（経路 a）で行います。
```

## 虚振動数の数が正しくないとき

クリーンな一次サドルは反応座標に沿った**支配的な虚モードを正確に 1 つ**だけ持ちます。よくある失敗は、余計な 2 つ目の小さな虚モードが出る、あるいは支配的な反応モードがまったく出ないことです。

| 症状 | 対処 |
| --- | --- |
| `n_imag = 0`（極小へ崩壊） | 失敗として扱い、TS 初期構造または MEP を改善する。`--flatten` は余分な負モードを除くだけで、失われた反応方向を作れない。 |
| `n_imag > 1` | バックエンドの本計算精度で再計算し、`--coord-type dlc` と残留モード向け `--flatten` を検討する。 |
| 1 モードだが原子運動が意図と異なる | 経路/初期構造を改善し IRC で接続性を確認する。モード数だけでは目的反応を同定できない。 |

`--flatten` は余分な虚モードの平坦化ループを実行します（`grad`:
Dimer loop、`hess`: RS-I-RFO 後処理）。`--no-flatten` は
`flatten_max_iter=0` に強制します。追加 Hessian 評価を伴うため opt-in
です。経路自体が粗い場合は TS 最適化の前に `all --refine-path`
（または `path-search`）で MEP を精密化します。再帰精密化は悪い経路を
複数セグメントに分割して計算量を大きく増やす場合があるため、これも
デフォルトでは無効です。

```bash
mlmm tsopt -i ts_guess.pdb --parm enzyme.parm7 -l 'LIG:Q' -b uma \
    --precision fp64 --coord-type dlc -o result_ts
```

`--coord-type` は最適化の座標系（`cart` | `redund` | `dlc` | `tric`、デフォルト `cart`）を選びます。`dlc`（非局在化内部座標）は低速ですが、ねじれの多い系やクリーンな一次サドルへの収束でより堅牢です。

```{warning}
`--coord-type dlc` は**Hessianベース**のオプティマイザが必要です。デフォルトの L-BFGS（`--opt-mode grad`）の [`opt`](opt.md) では警告を出して `cart` に戻ります。`tsopt`（RFO / RS-I-RFO）または `opt --opt-mode hess` で使ってください。`path-opt` / `path-search` は `cart` と `dlc` のみ受け付けます。`DLC + リンク原子` と `DLC + 3 層凍結 MM` は数値的に未検証なため、`cart` がデフォルトです。
```

同じ症状からの切り分けについては [典型エラー別レシピ — レシピ 4](recipes-common-errors.md#レシピ-4-収束後処理で止まる) を参照してください。

### 高度な MEP 参照モード

`--ref-mode` は MEP tangent を指定する非ゼロ Cartesian 3N vector
（`.npy` または空白区切り text）を読み込みます。これは一気通貫 workflow
向けの内部的な高度オプションで、`mlmm all` が MEP から生成して自動的に
渡します。通常の単独 `mlmm tsopt` では省略してください。同じ原子順の
vector を明示的に構築済みの場合だけ指定します。

## 対照を揃えた変異体 vs 野生型（あるいは機構 vs 機構）の比較

```{important}
変異体 vs 野生型（あるいは機構 vs 機構）のバリア比較では、**比較するすべてのモデルが同じ原子集合 — 原子数と残基が同一 — を使わなければなりません。** さもないとエネルギー基準が異なり、対照実験として成立しません。変異体上で幾何学的に再導出した ML/可動/凍結分割は、見かけの軟モードも生みます（`tsopt.n_imaginary ≥ 2`、いずれも微小 → IRC が中断）。
```

mlmm では、**野生型の B 因子レイヤーエンコーディングを変異構造に移植**し `--detect-layer` で実行することで、野生型の ML/MM レイヤリングを保持します。

| 手順 | 操作 |
| --- | --- |
| 1 | 変異体を構築。残基集合は野生型と同じに保つ（変異残基の種類だけが異なる）。 |
| 2 | 野生型の原子ごとの B 因子レイヤーコードを `(resid, atom-name)` で変異体にコピー: **ML = 0.0、MovableMM = 10.0、FrozenMM = 20.0**。 |
| 3 | `--detect-layer`（デフォルト有効）で実行し、*同一の*レイヤー割り当てを再利用する。 |

```bash
mlmm all -i mutant_layered.pdb -l 'LIG:Q' \
    --tsopt --thermo -o result_mutant
```

| フラグ | 操作 | 理由 |
| --- | --- | --- |
| `--detect-layer` | 維持（デフォルト `True`） | 移植した B 因子を読み、ML / 可動 / 凍結レイヤーを野生型と完全に同一にする |
| `-c/--center`、`-r/--radius` | **省略** | 幾何学的抽出は変異体上で*異なる*ポケットを再導出してしまう。`-c` を省くと抽出をスキップし、構造をそのまま使う |
| `-l/--ligand-charge` | 維持 | 非標準リガンドの電荷は標準アミノ酸表にない。`-l 'RES:Q'` が総電荷を自動導出（ハードコードの `-q` より推奨） |
| `--movable-cutoff` | 渡さない | `--detect-layer` を無効化する |

同じ原子集合の原則は、同一酵素上で 2 つの機構を比較するときも同様に適用されます。両モデルで原子集合を同一に保ち、反応座標だけを変えます。

## 実行例

コマンド形式は `mlmm tsopt -i TS_GUESS --parm PARM7 --model-pdb ML_REGION -q CHARGE -m MULT [options]` です。`mlmm tsopt --help` で主要オプション、`mlmm tsopt --help-advanced` で全オプション一覧を表示します。

デフォルト実行:

```bash
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_tsopt
```

light モード（Dimer）+ 解析的 Hessian:

```bash
# VRAM に余裕がある場合に light モード + 解析的Hessianで実行する
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --opt-mode grad --hessian-calc-mode Analytical --out-dir ./result_tsopt_grad
```

heavy モード（RS-I-RFO）+ YAML 上書き:

```bash
# heavy モードを YAML 上書きと併用する
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --opt-mode hess --config tsopt.yaml --out-dir ./result_tsopt_hess
# --dump で最適化軌跡を保存、--backend mace で MACE バックエンドを使用
```

## 処理の流れ

1. **入力処理** — 酵素 PDB、Amber トポロジー、ML 領域定義を読み込みます。電荷/スピンを解決します。CLI と YAML の凍結原子がマージされます。
2. **ML/MM calculatorの構築** — ML/MM calculator（MLIP バックエンド + hessian_ff）を構築します。`-b/--backend` で ML バックエンドを選択し（デフォルト: `uma`）、`--hessian-calc-mode` は MLIP が Hessian を解析的に評価するか有限差分で評価するかを制御します。v0.3.3 は機械的埋め込みを使用し、電子埋め込みの要求は構築前に拒否します。
3. **Light モード（Dimer）:**
   - Hessian Dimer ステージはアクティブ部分空間の部分 Hessian を評価して Dimer 方向を定期的に更新します。固定の constrained 処理は凍結 anchor と両立する全系剛体運動だけを除去します。保存・回転・試行する全方向で凍結Cartesian成分をゼロに保ち、中心外のforce評価でも凍結座標を中心imageと厳密に一致させます。
   - 平坦化ループが有効な場合（`--flatten`）、保存されたアクティブ Hessian は変位と勾配差分を使用した Bofill 更新により更新されます。各ループで虚振動数モードを推定し、1 回平坦化し、Dimer 方向を更新し、Dimer + L-BFGS マイクロセグメントを実行します。
4. **Heavy モード（RS-I-RFO）:**
   - RS-I-RFO オプティマイザを、`rsirfo` YAML セクションで定義されたオプションの Hessian 参照ファイルとマイクロサイクル制御とともに実行します。
   - `--flatten` が有効で収束後に 2 つ以上の虚振動数モードが残る場合、余分なモードを平坦化し、1 つだけ残るか平坦化反復上限に達するまで RS-I-RFO を再実行します。
5. **モードエクスポートと変換** — 収束した虚振動数モードは常に `vib/imag_*_trj.xyz` に書き出され、入力が PDB で変換が有効な場合は `.pdb` にもミラーリングされます。最適化軌跡と最終ジオメトリも `--dump` 時に入力テンプレート経由で PDB に変換されます。

## 出力

最終振動解析を行う場合、`result.json` の `status` が `converged` になるのは、
optimizer が収束し、かつ最終 Hessian の虚モードが正確に 1 つの場合だけです。
0 または複数なら `not_converged`、`--skip-final-freq` で鞍点次数を検証しない
場合は `unverified` です。

最適化が成功すると 3 種類の成果物が `result_tsopt/` に出力されます。

- `result_tsopt/final_geometry.pdb`（または `final_geometry.xyz`）
- `result_tsopt/vib/imag_*_trj.xyz`
- `result_tsopt/vib/imag_*.pdb`

```
out_dir/ (デフォルト: ./result_tsopt/)
├── result.json                    # --out-json 時。rigid_projection provenance を含む
├── final_geometry.xyz             # 常に書き出し
├── final_geometry.pdb             # 入力が PDB の場合
├── optimization_all_trj.xyz       # 連結 Dimer セグメント（--dump 時）
├── optimization_all.pdb           # 対応する PDB（--dump かつ入力が PDB の場合）
├── vib/
│   ├── imag_NN_±XXXX.XXcm-1_trj.xyz  # 虚振動数モード軌跡
│   └── imag_NN_±XXXX.XXcm-1.pdb      # 虚振動数モードに対応する PDB
└── .dimer_mode.dat                # Dimer 方向シード（light モード）
```

## CLI オプション

完全なフラグ一覧は生成された[コマンドリファレンス](../reference/commands/index.md)を参照してください。以下の表は説明が必要なオプションのみを扱い、網羅一覧の手動複製は行いません。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| **入力と電荷** | | |
| `-i, --input PATH` | 開始ジオメトリ（PDB または XYZ）。XYZ の場合はトポロジーに `--ref-pdb` を使用。 | 必須 |
| `--ref-pdb FILE` | 入力が XYZ の場合の参照 PDB トポロジー。 | _None_ |
| `--parm PATH` | 全酵素の Amber parm7 トポロジー。 | 必須 |
| `--model-pdb PATH` | ML 領域原子を含む PDB。`--detect-layer` 有効時はオプション。 | _None_ |
| `--model-indices TEXT` | ML 領域のカンマ区切り原子インデックス（範囲指定可）。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` を 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer / --no-detect-layer` | 入力 PDB の B 因子から ML/MM レイヤーを検出。 | `True` |
| `-q, --charge INT` | ML 領域の総電荷。 | _None_（`-l` 未指定時は必須） |
| `-l, --ligand-charge TEXT` | 残基ごとの電荷マッピング（例: `GPP:-3,SAM:1`）。`-q` 省略時に合計電荷を導出。PDB 入力または `--ref-pdb` が必要。 | _None_ |
| `-m, --multiplicity INT` | ML 領域のスピン多重度 (2S+1)。 | _None_（デフォルト 1） |
| **アクティブ領域の凍結** | | |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切りインデックス（YAML `geom.freeze_atoms` とマージ）。 | _None_ |
| `--hess-cutoff FLOAT` | ML 領域からの Hessian-MM 原子の距離カットオフ (Å)。未指定時は最終解析に必要な可動 MM 原子をすべて含めます。`0.0` で ML のみを評価する場合は、最終周波数解析も `--active-dof-mode ml-only` にします。エイリアス: `--radius-hessian`。 | _None_ |
| `--movable-cutoff FLOAT` | 可動 MM 原子の距離カットオフ (Å)。 | _None_ |
| **TS 探索とオプティマイザモード** | | |
| `--hessian-calc-mode CHOICE` | MLIP Hessian モード: `Analytical` または `FiniteDifference`。 | `FiniteDifference` |
| `--ref-mode PATH` | 高度な Cartesian 3N 経路方向ヒント。`all` が MEP から自動供給し、通常の単独 `tsopt` では省略。 | _None_ |
| `--max-cycles INT` | 最大総オプティマイザサイクル。 | `10000` |
| `--opt-mode CHOICE` | TS オプティマイザモード（Choice: `grad` / `hess` / `light` / `heavy` / `dimer` / `rsirfo` / `trim` / `rsprfo`）。`grad`/`light`/`dimer` → Hessian-Guided Dimer; `hess`/`heavy`/`rsirfo` → RS-I-RFO（デフォルト）; `trim` → TRIM（Helgaker）; `rsprfo` → RS-P-RFO（Banerjee）。Hessian TS オプティマイザ3種（`rsirfo`/`rsprfo`/`trim`）はいずれも microiter 対応。 | `hess` |
| `--microiter/--no-microiter` | マイクロイテレーション: 1 ステップの macro TS 移動（RS-I-RFO / RS-P-RFO / TRIM）+ MM 緩和（L-BFGS）を交互に実行。任意の Hessian モード（`hess`/`rsirfo`/`rsprfo`/`trim`）で有効。 | `True` |
| `--ml-only-hessian-dimer/--no-ml-only-hessian-dimer` | `grad` モードで Dimer 方向決定に ML 領域のみの Hessian を使用。高速だが精度は低下。 | `False` |
| **収束と平坦化** | | |
| `--thresh TEXT` | 収束プリセット（`gau_loose\|gau\|gau_tight\|gau_vtight\|baker\|never`）。 | _None_ |
| `--flatten/--no-flatten` | 余分な虚振動数モード平坦化ループの有効化/無効化。`--flatten` はデフォルト反復回数（50）を使用、`--no-flatten` は 0 に強制。light と heavy の両モードに適用。 | _None_（CLI デフォルトは無効 = `flatten_max_iter` 0; `--flatten` または YAML/config で初めて有効化され、その場合 50 回） |
| `--partial-hessian-flatten / --full-hessian-flatten` | 平坦化ループでの虚振動数モード検出に部分 Hessian（ML のみ）を使用。 | `True`（部分） |
| `--active-dof-mode CHOICE` | 最終振動解析のアクティブ自由度: `all`、`ml-only`、`partial`、`unfrozen`。 | `partial` |
| `--skip-final-freq/--no-skip-final-freq` | 収束後の振動解析と虚振動数モード平坦化をスキップ。大規模非凍結系で Hessian 対角化が高コストな場合に有用。TS の鞍点次数は未検証のままになる。 | `False` |
| **バックエンドと計算** | | |
| `-b, --backend CHOICE` | ML 領域の MLIP バックエンド: `uma`（デフォルト）、`orb`、`mace`、`aimnet2`。 | `uma` |
| `--precision [fp32\|fp64]` | MLIP バックエンド精度。省略時は UMA/AIMNet2 fp32、ORB/MACE fp64。AIMNet2 は fp64 を拒否。 | バックエンド依存 |
| `--workers INT` | UMA predictor worker 数。2 以上は `fairchem-core[extras]` が必要で、`Analytical` と併用不可。 | `1` |
| `--workers-per-node INT` | UMA 並列 predictor のノード当たり worker 数。 | _None_ |
| `--allow-charge-mult-mismatch` | 警告を出した上で ML 領域の電荷・多重度の電子パリティ検証を省略。意図した不一致の場合のみ使用。 | off |
| `--embedcharge/--no-embedcharge` | v0.3.3 では使用不可。旧コマンドを明示的に拒否するためにのみ残されています。 | `False` |
| `--embedcharge-cutoff FLOAT` | 廃止した電子埋め込み経路とともに使用不可。 | — |
| `--cmap/--no-cmap` | REAL と MODEL の両 MM 層で CMAP を保持します。 | `--cmap` |
| `--mm-backend [hessian_ff\|openmm]` | MM backend。 | `hessian_ff` |
| `--link-atom-method [scaled\|fixed]` | link atom 配置方式。 | `scaled` |
| **出力と設定** | | |
| `--dump/--no-dump` | 連結軌跡 `optimization_all_trj.xyz` を書き出し。 | `False` |
| `--convert-files/--no-convert-files` | PDB 入力時の XYZ/TRJ から対応する PDB の生成を切り替え。 | `True` |
| `-o, --out-dir TEXT` | 出力ディレクトリ。 | `./result_tsopt/` |
| `--config FILE` | 明示 CLI オプションより前に適用するベース YAML 設定ファイル。 | _None_ |
| `--show-config/--no-show-config` | 解決後の設定レイヤーを表示して実行を継続。 | `False` |
| `--out-json/--no-out-json` | machine-readable `result.json` を出力。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに入力/設定を検証し、実行計画を表示。`--help-advanced` に表示。 | `False` |

## YAML 設定

設定は **デフォルト < config < 明示CLI < override** の順で適用されます。
共有セクションは [YAML リファレンス](yaml-reference.md) を再利用します。ワークフローに合致している場合は以下のブロック全体をそのまま保持し、変更が必要な値のみ調整してください。

```yaml
geom:
 coord_type: cart                  # 座標タイプ: デカルト vs dlc 内部座標
 freeze_atoms: []                  # 1 始まり凍結原子（CLI/リンク検出とマージ）
 tr_projection: constrained        # 固定の内部 PHVA 処理
calc:
 charge: 0                         # 総電荷（CLI 上書き）
 spin: 1                           # スピン多重度 2S+1
mlmm:
 real_parm7: real.parm7            # Amber parm7 トポロジー
 model_pdb: ml_region.pdb          # ML 領域定義
 backend: uma                      # ML バックエンド (uma/orb/mace/aimnet2)
 embedcharge: false                # 互換性用。true は拒否される
 uma_model: uma-s-1p2              # uma-s-1p2 | uma-m-1p1
 uma_task_name: omol                # UMA タスク名 (backend=uma 時)
 ml_device: auto                   # ML デバイス選択
 hessian_calc_mode: Analytical        # Hessianモード（デフォルトは FiniteDifference; Analytical は opt-in）
opt:
 thresh: baker                     # 収束プリセット（Gaussian/Baker 式）
 max_cycles: 10000                 # オプティマイザサイクル上限
 print_every: 100                  # ログ出力間隔
 min_step_norm: 1.0e-08            # ステップ受け入れの最小ノルム
 assert_min_step: true             # ステップが閾値以下で停止
 rms_force: null                   # 明示的 RMS 力目標
 rms_force_only: false             # RMS 力収束のみに依存
 max_force_only: false             # 最大力収束のみに依存
 force_only: false                 # 変位チェックをスキップ
 converge_to_geom_rms_thresh: 0.05  # 参照への収束時の geom RMS 閾値
 overachieve_factor: 0.0           # 閾値を厳しくする係数
 check_eigval_structure: false     # Hessian固有値構造の検証
 line_search: true                 # ラインサーチを有効化
 dump: false                       # 軌跡/リスタートデータのダンプ
 dump_restart: false               # リスタートチェックポイントのダンプ
 prefix: ""                        # ファイル名プレフィックス
 out_dir: ./result_tsopt/          # 出力ディレクトリ
hessian_dimer:
 thresh_loose: gau_loose           # ゆるい収束プリセット
 thresh: baker                     # メイン収束プリセット
 update_interval_hessian: 500      # Hessian再構築間隔
 neg_freq_thresh_cm: 5.0           # 負の振動数閾値 (cm^-1)
 flatten_amp_ang: 0.1              # フラットニング振幅 (Å)
 flatten_max_iter: 50              # フラットニング反復上限（--no-flatten 時は無効）
 flatten_sep_cutoff: 0.0           # 代表原子間の最小距離 (Å)
 flatten_k: 10                     # モードごとにサンプルされる代表原子数
 flatten_loop_bofill: false        # フラットニング変位に対する Bofill 更新
 mem: 100000                       # ソルバーのメモリ上限
 device: auto                      # 固有値ソルバーのデバイス選択
 root: 0                           # 対象 TS ルートインデックス
 dimer:
  length: 0.0189                   # Dimer 間隔 (Bohr)
  rotation_max_cycles: 15          # 最大回転反復数
  rotation_method: fourier         # 回転オプティマイザ手法
  rotation_thresh: 0.0001          # 回転収束閾値
  rotation_tol: 1                  # 回転許容係数
  rotation_max_element: 0.001      # 回転行列の最大要素
  rotation_interpolate: true       # 回転ステップの補間
  rotation_disable: false          # 回転を完全に無効化
  rotation_disable_pos_curv: true  # 正の曲率検出時に無効化
  rotation_remove_trans: true      # 選択した剛体null成分を除去
  trans_force_f_perp: true         # 並進に垂直な力を射影
  bonds: null                      # 拘束用の結合リスト
  N_hessian: null                  # Hessianサイズの上書き
  bias_rotation: false             # 回転探索にバイアス
  bias_translation: false          # 並進探索にバイアス
  bias_gaussian_dot: 0.1           # ガウスバイアスの内積
  seed: null                       # 回転用の RNG シード
  write_orientations: true         # 回転方向を書き出し
  forward_hessian: true            # Hessianを前方伝播
 lbfgs:
  thresh: baker                    # L-BFGS 収束プリセット
  max_cycles: 10000                # 反復上限
  print_every: 100                 # ログ出力間隔
  min_step_norm: 1.0e-08           # 受け入れ最小ステップノルム
  assert_min_step: true            # ステップ停滞時にアサート
  max_step: 0.3                    # 最大ステップ長
  control_step: true               # 適応的ステップ長制御
  double_damp: true                # ダブルダンピングセーフガード
  keep_last: 7                     # L-BFGS バッファの履歴サイズ
  beta: 1.0                        # 初期ダンピングベータ
  mu_reg: null                     # 正則化強度
  max_mu_reg_adaptions: 10         # mu 適応の上限
rsirfo:
 thresh: baker                     # RS-IRFO 収束プリセット
 trust_radius: 0.10                # 初期信頼半径（ONIOM 向けに小さめ）
 trust_update: true                # 適応的信頼半径更新
 trust_min: 1.0e-04                # 最小信頼半径
 trust_max: 0.10                   # 最大信頼半径（ML/MM 安定性のため調整）
 max_energy_incr: null             # ステップごとの最大許容エネルギー増加
 hessian_update: bofill            # Hessian更新方式の上書き
 hessian_init: calc                # 初期Hessianソース
 hessian_recalc: 500               # N ステップごとにHessian再計算
 hessian_recalc_adapt: null        # 適応的Hessian再計算
 small_eigval_thresh: 1.0e-08      # 小固有値の閾値
 alpha0: 1.0                       # 初期シフトパラメータ
 max_micro_cycles: 50              # マクロサイクルごとのマイクロイテレーション
 track_mode_by_overlap: false      # ステップ間の固有ベクトル重なりで TS モードを追跡
```

```{tip}
最適化中に TS モードが別のルートに切り替わる場合（例: 複数の虚振動数が存在する場合）は `rsirfo.track_mode_by_overlap: true` を設定してください。
```

```{tip}
TS 収束が遅い場合や最適化中に TS モードが失われる場合は、`rsirfo` セクションの `hessian_recalc` を小さくしてみてください（例: 50--200）。正確なHessian再計算の頻度を上げることで、追加のHessian評価コストと引き換えに堅牢性が向上します。
```

## 注記

凍結境界 PHVA と質量加重 TR 処理は `freq.py` と共通です。
`constrained`（デフォルト）は、凍結 anchor をすべて動かさない全系剛体運動だけを
除去します。一般的な有効 rank は anchor が 0/1/2/非共線の 3 個以上のとき
6/3/1/0 で、実用的な ML/MM 境界では通常 0 です。全原子凍結は明示的なエラーになります。
Dimer は中心imageが変わるたびにこの基底を再構築して方向と回転forceに適用し、
凍結境界に対する有限曲率運動であるactive fragmentの並進は差し引きません。

固定の constrained 剛体モード処理と `--ref-mode` は別の機能です。後者は鞍点回復用の
高度な 3N MEP 接線を与えます。`geom.tr_projection` の古い非constrained値は明示的に拒否されます。
`--out-json` 時は
`result.json.rigid_projection` に treatment、有効 rank、Hessian source、Hessian shape を記録します。

```{note}
`rsirfo.trust_max` のデフォルトは 0.10 bohr です。TS 近傍での ML/MM 安定性が改善します。

共有 `opt` ブロックには **エネルギープラトー・フォールバック**（`energy_plateau: true`、`energy_plateau_thresh: 1.0e-4` au を `energy_plateau_window: 50` ステップにわたって適用、いずれもデフォルト）も備わっています。MLIP の力ノイズフロアが勾配ベースの `thresh` プリセットに到達できない場合、残りのサイクルを費やす代わりにプラトーで探索を停止し、`status: "stalled"` を報告します（`converged` とは区別される非収束の結果で、決して `converged` にはなりません）。終端の厳密 Hessian は実行されるため鞍点診断は報告されますが、flatten/retry ループは実行されません（stalled な root は検証済みの TS モードではないため）。詳細は [yaml-reference](yaml-reference.md#opt) を参照してください。
```

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) — 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) — 詳細なトラブルシューティングガイド

- [opt](opt.md) — 単一構造の構造最適化
- [freq](freq.md) — 検証済み TS の単一虚振動数を確認
- [irc](irc.md) — 最適化された TS からの反応経路追跡
- [all](all.md) — ML/MM モデル構築 -> MEP 探索 -> tsopt -> IRC -> freq を連鎖させる一気通貫ワークフロー
- [YAML リファレンス](yaml-reference.md) — `hessian_dimer`（Hessian ガイド付き Dimer）と `rsirfo` の完全な設定オプション
- [用語集](glossary.md) — TS、Dimer、RS-I-RFO、Hessian の定義
