# YAML 設定の一覧

YAML 設定ファイル（`--config`）に書けるキーとデフォルト値を、セクションごとに引くページです。セクションの一覧、優先順位、CLI フラグと YAML キーの対応も最初にまとめています。

| セクション | 説明 | 使用されるコマンド |
|---------|-------------|---------|
| [`geom`](#geom) | ジオメトリと座標設定 | all, opt, scan, scan2d, scan3d, tsopt, freq, irc, sp, dft, path-opt, path-search |
| [`calc`](#calc) | ML/MM 計算機の設定（別名: `mlmm:`） | all, opt, scan, scan2d, scan3d, tsopt, freq, irc, sp, dft, path-opt, path-search |
| [`sp`](#sp-セクション) | 一点計算の出力設定 | sp |
| [`opt`](#opt) | 最適化の共通設定 | all, opt, scan, scan2d, scan3d, tsopt, path-opt, path-search |
| [`lbfgs`](#lbfgs) | L-BFGS の設定 | all, opt, scan, scan2d, scan3d, tsopt（マイクロイテレーションの MM 緩和）, path-opt, path-search |
| [`rfo`](#rfo) | RFO の設定 | all, opt, scan, scan2d, scan3d, path-opt, path-search |
| [`gs`](#gs) | GSM（Growing String Method）設定 | all, path-opt, path-search |
| [`dmf`](#dmf) | DMF（Direct Max Flux）設定 | all, path-opt, path-search |
| [`irc`](#ja-irc-section) | IRC 積分設定 | all, irc |
| [`freq`](#ja-freq-section) | 振動解析設定 | all, freq |
| [`thermo`](#thermo) | 熱化学設定 | all, freq |
| [`dft`](#ja-dft-section) | DFT 計算設定 | all, dft |
| [`bias`](#bias) | 調和バイアス設定 | all, scan, scan2d, scan3d |
| [`bond`](#bond) | 結合変化検出設定 | all, scan, path-search |
| [`search`](#search) | 再帰的経路探索設定 | all, path-search |
| [`hessian_dimer`](#hessian_dimer) | Hessian Dimer による TS 最適化 | all, tsopt |
| [`rsirfo`](#rsirfo) | Hessian TS 最適化設定 | all, tsopt |
| [`stopt`](#stopt) | ストリング最適化の設定 | all, path-opt, path-search |
| [`microiter`](#microiter) | マイクロイテレーション（MM 緩和）の設定 | all, opt, tsopt |

(ja-yaml-configuration-precedence)=
## 設定の優先順位

設定は以下の順序で適用されます（後のものが前のものを上書き）:

```
組み込みデフォルト  <  --config (YAML)  <  CLI フラグ
```

1. **組み込みデフォルト** — `mlmm <subcmd> --help-advanced` と [コマンドの一覧（英語のみ）](../reference/commands/index.md) の `[default: …]` に出る値。
2. **`--config`** — デフォルトを上書きする YAML ファイル（例: `--config my_settings.yaml`）。
3. **CLI フラグ** — コマンドラインで明示的に指定したオプション（例: `-q -1`, `--thresh gau_loose`）。*明示的に指定された*値のみが YAML を上書きし、CLI デフォルトのままのオプションは YAML の値を上書きしません。

例: YAML で `calc.model_charge: 0` を設定し、CLI で `-q -1` を渡した場合、ML 領域の電荷は `-1` になります。

この優先順位は、`--config` を持つすべてのコマンドに共通です。

実行で使われる値は {ref}`-v 3 <ja-verbosity-levels>` で確かめられます。各セクションを、名前、`-` の下線、実際に使う値の順に表示します:

```text
opt
---
thresh: gau
max_cycles: 100000
…
```

セクション名を書き間違えると `[config] WARNING: YAML section(s) … are not recognized and were ignored.` が出て、そのセクションを使わずに実行を続けます。

(ja-common-cli-to-yaml-mapping)=
## 主要な CLI→YAML マッピング

| CLI フラグ | YAML キー | セクション |
|----------|----------|---------|
| `-q` / `--charge` | `model_charge` | `calc` |
| `-m` / `--multiplicity` | `model_mult` | `calc` |
| `-b` / `--backend` | `backend` | `calc` |
| `--backend-model` | `uma_model`・`orb_model`・`mace_model`・`aimnet2_model`（選んだバックエンドのキー） | `calc` |
| `--precision` | `uma_precision`・`orb_precision`・`mace_dtype`（選んだバックエンドのキー） | `calc` |
| `--workers` | `workers` | `calc` |
| `--link-atom-method` | `link_atom_method` | `calc` |
| `--mm-backend` | `mm_backend` | `calc` |
| `--cmap/--no-cmap` | `use_cmap` | `calc` |
| `--detect-layer/--no-detect-layer` | `use_bfactor_layers` | `calc` |
| `--hess-cutoff` | `hess_cutoff` | `calc` |
| `--movable-cutoff` | `movable_cutoff`（`use_bfactor_layers: false` にもする） | `calc` |
| `--embedcharge/--no-embedcharge` | `embedcharge` | `calc` |
| `--embedcharge-cutoff` | `embedcharge_cutoff` | `calc` |
| `--thresh` | `thresh` | `opt`・`tsopt`・スキャンのコマンドは `opt`、`path-opt`・`path-search` は `lbfgs` と `rfo` |
| `--thresh-gsm` | `thresh` | `stopt` |
| `--dmf-tol` / `--thresh-dmf` | `tol` | `dmf` |
| `--max-cycles` | `max_cycles` | コマンド別: `opt`/`tsopt`/`scan` は `opt`、`irc` は `irc` |
| `--max-cycles-gsm` | `max_cycles` | `stopt`（`stopt.stop_in_when_full` も設定） |
| `--dmf-max-iterations` / `--max-cycles-dmf` | `max_cycles` | `dmf` |
| `--gsm-param` | `param` | `gs` |
| `--max-nodes` | `max_nodes` | `gs` |
| `--preopt-max-cycles` | `max_cycles` | `lbfgs` と `rfo` |
| `--dump` | `dump` | `opt`（opt、tsopt、scan）、`stopt`（path-opt、path-search）、`thermo`（freq） |
| `--step-size`（irc） | `step_length` | `irc` |
| `--freeze-atoms` | `freeze_atoms`（YAML のリストと合わせる） | `geom` |
| `--coord-type` | `coord_type` | `geom` |
| `--temperature`（freq、`all --freq-temperature`） | `temperature` | `thermo` |
| `--pressure`（freq、`all --freq-pressure`） | `pressure_atm` | `thermo` |
| `--dft-engine` / `--engine` | `engine` | `dft` コマンドは `dft`、`--backend dft` は `calc.dft` |

### サブコマンド別の `--thresh` デフォルト

`--thresh` のデフォルトはサブコマンドごとに異なります。

| サブコマンド | デフォルト `--thresh` |
|------------|---------------------|
| `opt` | `gau` |
| `tsopt`（Hessian Dimer、RS-P-RFO、RS-I-RFO、TRIM） | `baker` |
| `scan` | `gau` |
| `scan2d`, `scan3d` | `baker` |
| `path-opt`、`path-search`（端点の事前最適化と、整列の後の緩和） | `gau` |
| `path-opt`、`path-search`（GSM のストリング: `--thresh-gsm`、`stopt.thresh`） | `gau_loose` |
| `path-opt`、`path-search`（DMF の経路: `--dmf-tol`、`dmf.tol`） | `tight`（0.04。Gaussian のプリセットではない） |
| `all`（`--thresh`: 単一構造の最適化とスキャンの緩和） | `gau` |
| `all`（`--thresh-post`: TS と IRC の後の端点の最適化） | `baker` |

受け付ける値: `gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`。実行ごとに `--thresh <preset>` または YAML の `opt.thresh` で上書きできます。マイクロイテレーションの MM 緩和は `microiter.micro_thresh` を使い、設定しないときはマクロステップのプリセットに従います。

```{note}
**`--thresh` を持たないサブコマンド。** `irc`、`freq`、`dft`、`sp` には `--thresh` が**ありません**:

- [`irc`](#ja-irc-section) — 収束は `irc.rms_grad_thresh`、`irc.energy_thresh`、`irc.max_cycles` で制御されます。IRC は予測子・修正子で経路を積分するので、最適化用の収束プリセットは使いません。
- `freq` と `sp` — 最適化ステップが無いため `--thresh` は存在しません。
- `dft` — SCF 収束は `dft.conv_tol`（デフォルト `1e-9` Hartree）と `dft.max_cycle` で制御されます。`gau`/`baker` のプリセットは使用しません。[`dft` セクション](#ja-dft-section) を参照してください。
```

## 共通セクション

### `geom`

ジオメトリ読み込みと座標系の設定。

```yaml
geom:
 coord_type: cart # opt と tsopt では "cart"（デカルト）・"redund"（冗長内部座標）・"dlc"（非局在化内部座標）・"tric"（並進・回転を含む内部座標）、all・path-opt・path-search では cart か dlc のみ
 freeze_atoms: [] # 1 始まりの固定原子インデックス
```

**注記:**
- YAML の `freeze_atoms` は、`--freeze-atoms` で指定した原子と合わせて使われます。
- 固定 MM 層の原子は力がゼロに設定され、Hessian の対応する列もゼロになります。
- デカルト座標の PHVA（部分 Hessian 振動解析）では、固定原子を動かさない全系の剛体運動だけを除きます。詳細は [freq](freq.md#固定境界での剛体モード) を参照してください。
- `irc` では、YAML や CLI の指定によらず `geom.coord_type` は `cart` です。

---

### `calc`

ML/MM 計算機の設定：入力ファイル、ML 領域とその電荷、MLIP バックエンド、MM バックエンド、層、Hessian。

```yaml
calc:
 # --- 入力ファイル ---
 input_pdb: null # 入力 PDB ファイルパス (CLI --input から設定)
 real_parm7: null # 全系の Amber parm7 トポロジー (CLI --parm7)
 model_pdb: null # ML 領域を定義する PDB (CLI --model-pdb)
 model_indices: null # model_pdb を省いたときの ML 原子のインデックス
 model_indices_base: 1 # YAML の model_indices だけに適用する 1 または 0
 model_charge: 0 # ML 領域の電荷。デフォルトなし。-q・-l・このキーのどれも無いと止まる
 model_mult: 1 # ML 領域のスピン多重度 (CLI -m で上書き)
 link_mlmm: null # null: parm7 結合から自動決定; list: 明示上書き
 link_atom_method: scaled    # リンク原子配置: "scaled" (g-factor) または "fixed" (1.09/1.01 Å)

 # --- MLIP バックエンドの選択 ---
 backend: uma # 高レベルバックエンド: uma, orb, mace, aimnet2, dft

 # --- UMA バックエンドの設定 ---
 uma_model: uma-s-1p2 # UMA モデル名: uma-s-1p2, uma-m-1p1
 uma_task_name: omol # UMA バッチに記録されるタスクタグ (backend=uma 時)
 uma_precision: fp32 # fp32 | fp64 (UMA バックエンドの数値精度)

 # --- ORB バックエンドの設定 ---
 orb_model: orb_v3_conservative_omol  # ORB モデル名 (backend=orb 時)
 orb_precision: float64  # ORB 浮動小数点精度 (backend=orb 時; "float32-high" は TF32 matmul で --precision fp32 でも選択可、"float32" も受け付ける)

 # --- MACE バックエンドの設定 ---
 mace_model: MACE-OMOL-0 # MACE モデル名 (backend=mace 時)
 mace_dtype: float64      # MACE 浮動小数点精度 (backend=mace 時)

 # --- AIMNet2 バックエンドの設定 ---
 aimnet2_model: aimnet2   # AIMNet2 モデル名 (backend=aimnet2 時)

 # --- PySCF/GPU4PySCF の高レベルバックエンド ---
 dft:
  func_basis: wb97m-v/def2-svp
  engine: gpu # gpu | cpu
  lowmem: true # DF テンソルを保持しない direct JK
  scf_stepwise_grid: true # 最初の SCF を粗いグリッドから（--scf-stepwise-grid）
  density_fit: false # --no-dft-low-memory のときはデフォルトで有効
  nprocs: auto # スケジューラ/affinity から PySCF のスレッド数を決定
  memory: auto # ホストの RAM 上限（例 64GB、GPU の VRAM ではない）
  save_scf_checkpoint: false
  checkpoint_path: null # 有効時のデフォルト（all 以外）: <out-dir>/_work/dft_scf/state.chk
  pyscf:
   mol: {}
   mf: {}
   grids: {}
   density_fit: {}
   with_df: {}

 # --- ML のデバイスと Hessian ---
 ml_device: auto # ML デバイス: "cuda", "cpu", "auto"
 ml_cuda_idx: 0 # CUDA デバイスインデックス
 hessian_calc_mode: FiniteDifference # ML Hessianモード: "FiniteDifference" または "Analytical"

 # --- 静電埋め込み（任意） ---
 embedcharge: false # MLIP は高コストな xTB 補正、dft は PySCF への直接静電埋込み
 embedcharge_cutoff: 12.0 # ML 領域からの MM 点電荷カットオフ (Å)
 embedcharge_step: 0.001 # MLIP/xTB 補正用の数値 Hessian ステップ（dft では未使用）
 xtb_cmd: xtb # xTB 実行コマンド
 xtb_acc: 0.2 # xTB 精度パラメータ
 xtb_workdir: tmp # xTB 作業ディレクトリ
 xtb_keep_files: false # xTB 一時ファイルを保持
 xtb_ncores: 4 # xTB プロセス数

 # --- MM バックエンドの設定 ---
 mm_backend: hessian_ff # MM バックエンド: "hessian_ff" | "openmm"。Hessian の方法は下の mm_fd で選ぶ
 use_cmap: true         # parm7 の CMAP を REAL と MODEL の両 MM 層で保持
 mm_device: cpu # MM デバイス (hessian_ff は CPU のみ、OpenMM は CUDA/CPU 対応)
 mm_cuda_idx: 0 # MM CUDA インデックス (OpenMM のみ)
 mm_threads: 16 # MM 計算のスレッド数
 workers: 1 # ローカル ML worker process 数（対応 backend のみ）
 workers_per_node: 1 # workers > 1 のときのノードあたりの worker 数（UMA の並列 predictor）
 mm_fd: true # MM Hessianに有限差分を使用
 mm_hessian_mode: null # 明示指定は finite_difference/analytical。null は mm_fd に従う
 mm_fd_dir: null # MM Hessianログの出力ディレクトリ
 mm_fd_delta: 0.001 # MM Hessian の有限差分の変位（Å）

 # --- Hessian の出力設定 ---
 out_hess_torch: true # Hessianを torch.Tensor で返す
 H_double: true # Hessianを float64 で組み立て・返却
 symmetrize_hessian: true # 最終Hessianを 0.5*(H+H^T) で対称化
 return_partial_hessian: true # アクティブブロック部分Hessian（CLI ラッパー側で true デフォルトを適用）

 # --- 層の設定 ---
 freeze_atoms: [] # geom.freeze_atoms から継承
 hess_cutoff: null # Å: null = 可動 MM をすべて Hessian 対象に含める (デフォルト)、>0.0 で ML 周辺の指定距離内 MM のみに限定
 movable_cutoff: null # Å: ML からこの距離以内の MM を可動にする（null は freeze_atoms に従う）
 use_bfactor_layers: true # 入力 PDB の B-factor から層を読み取り
 hess_mm_atoms: null # 明示的 Hessian 対象 MM 原子インデックス（1 始まり、カットオフより優先）
 movable_mm_atoms: null # 明示的 可動 MM 原子インデックス（1 始まり、カットオフより優先）
 frozen_mm_atoms: null # 明示的 固定 MM 原子インデックス（1 始まり、カットオフより優先）

 # --- 診断 ---
 print_timing: true # ML/MM Hessianのタイミング内訳を表示
 print_vram: true # CUDA VRAM 使用量を表示
```

**注記:**
- CLI の `--model-indices` は常に 1 始まりです。YAML の `calc.model_indices` にはリストか範囲の文字列を書けます。0 始まりの YAML データのときだけ `calc.model_indices_base: 0` を指定します。
- `backend` は高レベルのバックエンドを選びます。`uma`（デフォルト）、`orb`、`mace`、`aimnet2`、PySCF/GPU4PySCF の `dft` から選べます。
- `backend: dft` で `embedcharge: true` のときは、MM の点電荷を PySCF に直接入れ、MM 原子にかかる力も計算します。Hessian は ML/MM 全体の力の有限差分なので、ML と MM の間の応答も含みます。
- バックエンド固有のモデルキーは、対応するバックエンドが選択されている場合にのみ有効です:
  - `uma_model`、`uma_task_name` — UMA バックエンドのみ
  - `orb_model`、`orb_precision` — ORB バックエンドのみ
  - `mace_model`、`mace_dtype` — MACE バックエンドのみ
  - `aimnet2_model` — AIMNet2 バックエンドのみ
- `hessian_calc_mode: Analytical` はバックエンドの解析 Hessian を明示的に要求します。UMA、ORB、MACE、AIMNet2 がこの経路を実装しており、インストール済みバックエンドが非対応なら計算法を暗黙に変更せずエラーになります。`workers > 1` との併用もエラーです。
- `opt`/`tsopt`/`irc`/`freq` は、YAML で `calc.return_partial_hessian` を明示しない場合に部分 Hessian をデフォルトで使用します。
- これらのコマンドで完全 Hessian を強制するには `calc.return_partial_hessian: false` を明示してください。
- `mm_fd: true` は有限差分 MM Hessian、`false` は `hessian_ff` の解析 MM Hessian を使います。`mm_hessian_mode` は同じ選択を名前（`finite_difference` か `analytical`）で指定し、`null` の場合は `mm_fd` で決まります。
- `use_cmap: true`（デフォルト）は parm7 に含まれる CMAP を REAL と MODEL の両 MM 層で保持します。明示的な改変力場計算だけ `false` を指定してください。この場合は両層から CMAP を除去します。
- standalone ML/MM 計算には `real_parm7` が必須です。ML 領域は `model_pdb`、明示的な model index、または有効な B-factor layer から指定できます。
- `-q` と `-m` を明示すると、YAML の `model_charge` と `model_mult` より優先されます。省いた値は YAML から取ります。

---

### `opt`

L-BFGS/RFO で共通の最適化設定。ここに書いた全キーが、`--microiter` の有無によらず `opt` コマンドのオプティマイザに届きます。`tsopt` のマクロオプティマイザにも届き、その上に [`rsirfo`](#rsirfo) / [`hessian_dimer`](#hessian_dimer) が重なります。転送されるのは**実際に変更した値だけ**なので、触っていないキーについてはオプティマイザ固有セクションが優先されます。同じ YAML ファイルで同じ設定に違う値を書くとエラーで止まります。共通の `opt` のキーと、選んだオプティマイザのセクションの同じキーの組も対象です。

```yaml
opt:
 thresh: gau # 収束プリセット: gau_loose, gau, gau_tight, gau_vtight, baker, never
 max_cycles: 100000 # オプティマイザサイクル上限
 print_every: 100 # ログ出力間隔
 min_step_norm: 1.0e-08 # 最小ステップノルム
 assert_min_step: true # ステップが閾値以下で停止
 rms_force: null # 明示的 RMS 力ターゲット
 rms_force_only: false # RMS 力のみで収束判定
 max_force_only: false # 最大力のみで収束判定
 force_only: false # 変位チェックをスキップ
 converge_to_geom_rms_thresh: 0.05 # 参照ジオメトリへの収束 RMS 閾値
 overachieve_factor: 0.0 # 0.0 で無効。正の値では力が閾値/係数を下回ると step 基準なしで収束（baker では不使用）
 check_eigval_structure: false # Hessian固有値構造の検証
 energy_plateau: false # opt-in（--stop-plateau）: エネルギーが停滞したら stalled として停止 (収束扱いにはしない)
 energy_plateau_thresh: 1.0e-4 # エネルギー変動許容幅 au（約 0.06 kcal/mol）
 energy_plateau_window: 50 # プラトー判定に用いる直近ステップ数
 line_search: true # ラインサーチを有効化
 dump: false # 軌跡/リスタートデータの出力
 dump_restart: false # リスタートチェックポイントの出力
 prefix: "" # ファイル名の接頭辞
 out_dir: ./result_opt/ # 出力ディレクトリ
```

**収束プリセット**（デカルト座標で、力は Hartree/Bohr、ステップは Bohr）:

| プリセット | Max Force | RMS Force | Max Step | RMS Step |
|-----------|-----------|-----------|----------|----------|
| `gau_loose` | 2.5e-3 | 1.7e-3 | 1.0e-2 | 6.7e-3 |
| `gau` | 4.5e-4 | 3.0e-4 | 1.8e-3 | 1.2e-3 |
| `gau_tight` | 1.5e-5 | 1.0e-5 | 6.0e-5 | 4.0e-5 |
| `gau_vtight` | 2.0e-6 | 1.0e-6 | 6.0e-6 | 4.0e-6 |
| `baker` | 3.0e-4 | 2.0e-4 | 3.0e-4 | 2.0e-4 |

`baker` は表の 4 列すべてに加えて、直前のサイクルとの `|delta E| < 1e-6` Hartree を要求します。これは Bakken と Helgaker（*J. Chem. Phys.* **117**, 9160 (2002)）が示した Baker 基準（`max(|force|) <= 3e-4` **かつ**（`|delta E| < 1e-6` **または** `max(|step|) <= 3e-4`））より厳しい条件です。文献の形では RMS の力が残った構造も収束とみなされ、機械学習ポテンシャルの面では高次の鞍点で止まることがあるため、厳しい形を使います。`min_step_norm` 以下のステップは、エネルギーの条件を満たすとみなします。

ほかのプリセットでは、`overachieve_factor` を 0 より大きくすると、`max(force)` と `rms(force)` がともに `閾値 / overachieve_factor` を下回った時点で、ステップの条件を満たしていなくても収束とします。デフォルトは `0.0`（無効）で、`baker` では使いません。

**エネルギープラトー停止（opt-in、デフォルト無効）:**

`energy_plateau` のデフォルトは `false` です。`opt` / `tsopt` / `all` の `--stop-plateau` で有効化し、`--stop-plateau-thresh` / `--stop-plateau-window` が上記の 2 つの値を設定します。有効時、直近 `energy_plateau_window` ステップのエネルギー範囲 `max(E) - min(E)` が `energy_plateau_thresh` を下回ると、オプティマイザを `optimization_status: "stalled"` で停止します。

MLIP の力のノイズで力が収束閾値を下回らないとき、サイクルを節約できます。ただしエネルギーの平坦化は停留点の証拠ではないため、この停止は明示的に指定したときだけ働き、実行の実質的な上限は常に `max_cycles` です。

GSM・DMF などの chain-of-states 最適化と、`--microiter` の MM 緩和では、プラトー判定を行いません。MM 緩和で止めると周辺環境が緩和されないまま終わるためです。`microiter.micro_max_cycles` は MM 緩和の回数上限で、通常は収束した時点で終了します。

---

### `lbfgs`

L-BFGS の設定（`opt` を拡張）。

```yaml
lbfgs:
 keep_last: 7 # L-BFGS バッファの履歴サイズ
 beta: 1.0 # 初期ダンピング beta
 gamma_mult: false # 乗法的 gamma 更新
 max_step: 0.3 # 最大ステップ長
 control_step: true # 適応的ステップ長制御
 double_damp: true # 二重ダンピング安全装置
 mu_reg: null # 正則化強度
 max_mu_reg_adaptions: 10 # mu 適応の上限
 reject_uphill: false # 許容値を超えるエネルギー上昇の拒否を明示的に有効化
 uphill_tolerance: 0.0001 # エネルギー上昇の許容値（Hartree）
 rejection_step_floor: 1.0e-07 # 再試行ステップの下限
 max_rejections_at_floor: 3 # 下限での連続拒否後に停止
```

---

### `rfo`

RFO（Rational Function Optimizer）の設定（`opt` を拡張）。

```yaml
rfo:
 trust_radius: 0.10 # 信頼領域半径
 trust_update: true # 信頼領域更新を有効化
 trust_min: 0.0001 # 最小信頼半径
 trust_max: 0.10 # 最大信頼半径（ML/MM 安定性のため調整）
 max_energy_incr: null # ステップあたりの許容エネルギー増加
 reject_uphill: false # 許容値を超えるエネルギー上昇の拒否を明示的に有効化
 uphill_tolerance: 0.0001 # エネルギー上昇の許容値（Hartree）
 rejection_trust_floor: 1.0e-07 # 再試行の信頼半径の下限
 max_rejections_at_floor: 3 # 下限での連続拒否後に停止
 hessian_update: ts_bfgs # Hessian更新スキーム: ts_bfgs, bfgs, bofill 等
 hessian_init: calc # Hessian初期化: calc, unit 等
 hessian_recalc: 500 # N ステップごとにHessianを再構築
 hessian_recalc_adapt: null # 適応的Hessian再構築係数
 small_eigval_thresh: 1.0e-08 # 安定性のための固有値閾値
 alpha0: 1.0 # 初期マイクロステップ
 max_micro_cycles: 50 # 1 step 内の RS 反復の上限（ML/MM のマイクロイテレーションとは別）
 rfo_overlaps: false # RFO オーバーラップを有効化
 gediis: false # GEDIIS を有効化
 gdiis: true # GDIIS を有効化
 gdiis_thresh: 0.0025 # GDIIS 受容閾値
 gediis_thresh: 0.01 # GEDIIS 受容閾値
 gdiis_test_direction: true # DIIS 前に降下方向をテスト
 adapt_step_func: true # 適応的ステップスケーリング
```

---

### `microiter`

ML/MM 最適化のマイクロイテレーションの設定です。`--microiter` を有効にすると、ML 領域のマクロステップの合間に、ML 原子を固定したまま MM 領域を L-BFGS で緩和します。

```yaml
microiter:
 micro_thresh: null       # MM緩和の収束プリセット（L-BFGS）; null → マクロステップと同じ
 micro_max_cycles: 100000 # マイクロイテレーションサイクル上限
```

**注記:**
- CLI の `--microiter` / `--no-microiter` で切り替えます（デフォルトは有効）。
- `opt --opt-mode hess` と、すべての Hessian TS mode（`hess`, `rsirfo`, `rsprfo`, `trim`）で使用可能
- ML 原子を固定したまま、L-BFGS で MM 領域の力を最小化します
- `micro_thresh` には `opt.thresh` と同じプリセットを書けます。`null` か省略のときは、マクロステップと同じ閾値を使います。

---

## 経路最適化セクション

### `gs`

Growing String Method（GSM）の設定。

```yaml
gs:
 fix_first: true # 最初の端点を固定
 fix_last: true # 最後の端点を固定
 max_nodes: 20 # 最大ストリングノード数
 perp_thresh: 0.005 # 垂直変位閾値
 reparam_check: rms # 再パラメータ化チェック指標
 reparam_every: 1 # 再パラメータ化間隔
 reparam_every_full: 1 # 完全再パラメータ化間隔
 param: equi # パラメータ化スキーム
 max_micro_cycles: 10 # 1 step 内の RS 反復の上限（ML/MM のマイクロイテレーションとは別）
 reset_dlc: true # 各ステップで非局在化座標を再構築
 climb: true # クライミングイメージを有効化
 climb_rms: 0.0005 # クライミング RMS 閾値
 climb_lanczos: true # クライミングの Lanczos 精密化
 climb_lanczos_rms: 0.0005 # Lanczos RMS 閾値
 climb_fixed: false # クライミングイメージを固定
 scheduler: null # オプションのスケジューラバックエンド
```

`gs.param` には `equi` か `energy` を書けます。`energy` は GSM のストリングが伸びきった後にだけ効き、エネルギーの高い領域にノードを寄せます。CLI では `--gsm-param` で指定します。

---

### `dmf`

Direct Max Flux（DMF）による MEP 最適化。

```yaml
dmf:
 max_cycles: 3000 # DMF/IPOPT反復上限
 tol: tight # IPOPT dual_inf_tol: tight(0.04) | middle(0.10) | loose(0.20) または正の float（--dmf-tol で上書き）
 correlated: true # 相関 DMF 伝搬
 sequential: true # 逐次 DMF 実行
 fbenm_only_endpoints: false # 端点を超えて FB-ENM を実行
 fbenm_options:
   delta_scale: 0.2 # FB-ENM 変位スケーリング
   bond_scale: 1.25 # 結合カットオフスケーリング
   fix_planes: true # 平面拘束の強制
 cfbenm_options:
   bond_scale: 1.25 # CFB-ENM 結合カットオフスケーリング
   corr0_scale: 1.1 # corr0 の相関スケール
   corr1_scale: 1.5 # corr1 の相関スケール
   corr2_scale: 1.6 # corr2 の相関スケール
   eps: 0.05 # 相関イプシロン
   pivotal: true # ピボット残基の処理
   single: true # 単一原子ピボット
   remove_fourmembered: true # 四員環の除去
 dmf_options:
   remove_rotation_and_translation: false # 剛体運動を保持
   mass_weighted: false # 質量重み付けの切替
   parallel: false # 並列 DMF を有効化
   eps_vel: 0.01 # 速度許容値
   eps_rot: 0.01 # 回転許容値
   beta: 10.0 # DMF の beta パラメータ
   update_teval: false # 遷移評価の更新
 k_fix: 300.0 # 拘束の調和定数
```

`dmf.tol` は DMF ソルブが最後に適用する許容値なので、同じファイル内の `ipopt_options.dual_inf_tol` より優先されます。生の IPOPT オプションを固定したい場合は `dmf.tol` を書かず `ipopt_options.dual_inf_tol` のみを指定してください。`gau_tight` などの Gaussian プリセットはここでは拒否され、`--thresh` / `--thresh-gsm` の担当です。

---

### `search`

再帰的経路探索（path-search のみ）。

```yaml
search:
 max_depth: 10 # 許可する再帰分割の階層数（0 = 分割しない）
 stitch_rmsd_thresh: 0.0001 # セグメント縫合の RMSD 閾値
 bridge_rmsd_thresh: 0.0001 # ブリッジノードの RMSD 閾値
 max_nodes_segment: 20 # セグメントあたりの最大ノード数
 max_nodes_bridge: 5 # ブリッジあたりの最大ノード数
 kink_max_nodes: 3 # キンク（kink）最適化の最大ノード数
 max_seq_kink: 2 # 連続キンクの上限
 refine_mode: null # 精密化戦略: peak, minima, null (自動)
```

---

### `stopt`

ストリング最適化（GS/DMF）の設定。path-opt と path-search で使用。

```yaml
stopt:
 type: string           # 最適化タイプのラベル
 thresh: gau_loose      # ストリング最適化の収束プリセット（--thresh-gsm で上書き）
 stop_in_when_full: 300 # ストリングが満杯になったときの早期停止閾値
 align: false           # アライメントトグル
 scale_step: global     # ステップスケーリングモード
 max_cycles: 300         # ストリング最適化サイクル上限
 dump: false            # 軌跡/リスタートデータ出力
 dump_restart: false    # リスタートチェックポイントの出力
 reparam_thresh: 0.0    # 再パラメータ化閾値
 coord_diff_thresh: 0.0 # 座標差分閾値
 out_dir: ./result_path_opt/  # 出力ディレクトリ
 print_every: 10        # ログ出力間隔
 lbfgs:
   # 単一構造最適化用（端点の事前最適化、HEI±1、キンクノード）
   thresh: gau
   # max_cycles: 100000 # 任意の上書き
   # ...（詳細は lbfgs セクション参照）
```

**注記:**
- `stopt.lbfgs` / `stopt.rfo` は端点の事前最適化、HEI±1 精密化、キンクノードに使う単一構造オプティマイザを設定します。
- 外側の `stopt` キーはストリング最適化を制御します。
- path 系のコマンドは、実行ごとの出力先と接頭辞を自分で決めます。

---

## TS 最適化セクション

### `hessian_dimer`

`tsopt --opt-mode grad` の Hessian Dimer による TS 最適化の設定です。`opt.thresh` と `hessian_dimer.thresh` を両方書くときは同じ値にします。片方だけならその値を、どちらも無ければ Dimer のデフォルト値を使います。

```yaml
hessian_dimer:
 thresh_loose: gau_loose # 緩い収束プリセット
 thresh: baker # メイン収束プリセット
 update_interval_hessian: 500 # Hessian再構築間隔
 flatten_amp_ang: 0.1 # flattening 振幅 (Å)
 flatten_max_iter: 50 # flattening を有効にしたとき（--flatten かこのキー）の反復上限。tsopt はデフォルトでは flattening を行わない
 flatten_sep_cutoff: 0.0 # 代表原子間の最小距離
 flatten_k: 10 # モードあたりのサンプル代表原子数
 flatten_loop_bofill: false # flattening 変位に Bofill 更新
 mem: 100000 # ソルバーのメモリ上限
 device: auto # 固有値ソルバーのデバイス選択
 root: 0 # ターゲット TS ルートインデックス
 partial_hessian_flatten: true # 部分Hessianを虚モード検出に使用
 ml_only_hessian_dimer: false # Dimer 方向決定に ML 領域のみのHessianを使用
 dimer:
   length: 0.0189 # Dimer 間隔 (Bohr)
   rotation_max_cycles: 15 # 最大回転反復数
   rotation_method: fourier # 回転最適化手法
   rotation_thresh: 0.0001 # 回転収束閾値
   rotation_tol: 1 # 回転許容係数
   rotation_max_element: 0.001 # 回転行列の最大要素
   rotation_interpolate: true # 回転ステップの補間
   rotation_disable: false # 回転を完全に無効化
   rotation_disable_pos_curv: true # 正曲率検出時に回転を無効化
   rotation_remove_trans: true # 選択した剛体null成分を除去
   trans_force_f_perp: true # 並進に垂直な力の投影
   bonds: null # 拘束用の結合リスト
   N_hessian: null # Hessianサイズの上書き
   bias_rotation: false # 回転探索のバイアス
   bias_translation: false # 並進探索のバイアス
   bias_gaussian_dot: 0.1 # ガウスバイアスの内積
   seed: null # 回転の乱数シード
   write_orientations: false # 回転方向の書き出し（明示的な true も可）
   forward_hessian: true # Hessianの前方伝搬
 lbfgs:
   # lbfgs セクションと同じキー
   thresh: baker
   line_search: false # 必須: Dimer の有効力は物理エネルギーと共役でない
```

**注記:**
- `--flatten` を付けない `tsopt` は、YAML で有効にしない限り flattening を行いません。有効にしたときの反復の上限は `flatten_max_iter`（デフォルト 50）です。
- `tsopt` と `all` の `--flatten` はデフォルトの `flatten_max_iter` で flattening を有効にし、`--no-flatten` は `flatten_max_iter` を 0 にします。`--flatten` と同時に YAML で `flatten_max_iter` を明示した場合は、YAML の値が優先されます。
- 内側の L-BFGS 固有設定は、最上位の `lbfgs` ではなく `hessian_dimer.lbfgs` に置きます。共通の `print_every` と `energy_plateau*` は上記の競合規則に従います。`line_search` は `false` 固定で、Dimer の有効力は表示する物理エネルギーの勾配ではないため `true` は拒否されます。`max_cycles` は設定できず、各 segment には `opt.max_cycles` の残り cycle 数が渡されます。

---

### `rsirfo`

Hessian TS 最適化の共通設定です。デフォルトの RS-P-RFO（`tsopt --opt-mode hess` / `rsprfo`）と、明示的な `rsirfo` / `trim` に適用されます。

```yaml
rsirfo:
 thresh: baker # Hessian TS 収束プリセット
 max_cycles: 100000 # opt.max_cycles と共有するサイクル上限
 print_every: 100 # ログ出力間隔
 min_step_norm: 1.0e-08 # 最小ステップノルム
 assert_min_step: true # ステップ停滞時にアサート
 roots: [0] # 追跡する root は 1 個のみ（空のリスト・複数の root は拒否）
 hessian_ref: null # 参照Hessian
 rx_modes: null # 反応モード定義
 prim_coord: null # 監視する主座標
 rx_coords: null # 監視する反応座標
 hessian_update: bofill # Hessian更新スキーム
 hessian_recalc_reset: true # 正確なHessian後に再計算カウンタをリセット
 hessian_init: calc # Hessian初期化
 hessian_recalc: 500 # Hessian再構築間隔
 max_micro_cycles: 50 # 1 step 内の RS 反復の上限（ML/MM のマイクロイテレーションとは別）
 augment_bonds: false # 結合解析に基づく反応経路の拡張
 min_line_search: false # 常に false: RS-P-RFO は line search を使わない
 max_line_search: false # 常に false: RS-P-RFO は line search を使わない
 assert_neg_eigval: false # 収束時に負の固有値を要求
 track_mode_by_overlap: false # mlmm 固有: オーバーラップでターゲットモードを追跡
 trust_radius: 0.10 # 信頼領域半径
 trust_update: true # 信頼領域更新
 trust_min: 0.0001 # 最小信頼半径
 trust_max: 0.10 # 最大信頼半径（ML/MM 安定性のため調整）
 small_eigval_thresh: 1.0e-08 # 安定性のための固有値閾値
 out_dir: ./result_tsopt/ # 出力ディレクトリ
```

RS-P-RFO は line search を使いません。`min_line_search` または `max_line_search` に `true` を書くと警告を出して `false` に戻します。

`opt` と `rsirfo` に同じ設定を明示する場合は値を一致させてください。片方だけならその値、無指定なら `rsirfo` のデフォルトを使います。

---

## IRC セクション

(ja-irc-section)=
### `irc` セクション

IRC 積分設定。

```yaml
irc:
 step_length: 0.1 # 積分ステップ長
 never_stop: false # 物理的端点判定を無視してmax_cyclesまで追跡
 max_cycles: 125 # IRCステップ上限
 forward: true # 順方向に伝搬
 backward: true # 逆方向に伝搬
 root: 0 # 基準振動モードのルートインデックス
 hessian_init: calc # Hessian初期化ソース
 hessian_update: bofill # Hessian更新スキーム
 hessian_recalc: null # Hessian再構築間隔
 displ: energy # 変位構築方法
 displ_energy: 0.001 # エネルギーベースの変位スケーリング
 displ_length: 0.1 # 長さベースの変位フォールバック
 rms_grad_thresh: 0.001 # RMS 勾配の収束閾値
 hard_rms_grad_thresh: null # ハード RMS 勾配停止閾値
 energy_thresh: 0.000001 # エネルギー変化閾値
 energy_increase_thresh: 0.0   # 通常のモードでは、1 ステップでもエネルギーが上がると停止
 imag_below: 0.0 # 虚振動数カットオフ
 force_inflection: true # 変曲点検出の強制
 check_bonds: false # 伝搬中の結合チェック
 out_dir: ./result_irc/ # 出力ディレクトリ
 prefix: "" # ファイル名の接頭辞
 dump_fn: irc_data.h5 # IRC データファイル名
 dump_every: null # デフォルトでは無効。有効化する場合のみ正の間隔を指定
 max_pred_steps: 500 # 予測子-修正子の最大ステップ数
 loose_cycles: 3 # 引き締め前の緩いサイクル数
 corr_func: mbs # EulerPC の修正子関数
```

---

## 振動解析セクション

(ja-freq-section)=
### `freq` セクション

振動解析設定。

```yaml
freq:
 active_dof_mode: partial # アクティブ原子の選択: "all" | "ml-only" | "partial" | "unfrozen"
 zero_cutoff_cm: 5.0 # ν < −zero_cutoff_cm のモードを虚振動とみなす
 amplitude_ang: 0.8 # モード変位振幅 (Å)
 n_frames: 20 # モードの軌跡のフレーム数
 max_write: 10 # 書き出すモードの最大数
 sort: value # ソート順: "value" または "abs"
 out_dir: ./result_freq/ # 出力ディレクトリ
```

`freq.zero_cutoff_cm` は、単体の `freq`、`opt` の flatten、Dimer、Hessian を使う TS 最適化で共通です。`hessian_dimer.neg_freq_thresh_cm` と `rsirfo.saddle_imaginary_threshold_cm` は同じ閾値の別名で、異なる値を書くとエラーで止まります。

デフォルトでは ν < −5.00 cm⁻¹ を虚振動と分類します。n_imag は鞍点の次数を表し、`n_negative_modes` はすべての負の振動数を数えます。どちらも最適化の収束を変えません。符号付き物理モードと熱化学に使う正のモードはすべて保持します。

**注記:**
- `active_dof_mode`: 振動解析に参加させる原子集合を選択します。`all` は全原子、`ml-only` は ML 領域のみ、`partial`（デフォルト）は ML + 可動 MM、`unfrozen` は固定されていない全原子を使用します。CLI フラグ `--active-dof-mode` が明示された場合は YAML 値より優先されます。

---

### `thermo`

熱化学設定。

```yaml
thermo:
 temperature: 298.15 # 熱化学温度 (K)
 pressure_atm: 1.0 # 熱化学圧力 (atm)
 symmetry_number: null # 自動判定。正整数は高度な上書き指定
 dump: false # thermoanalysis.yaml の書き出し
```

---

## 一点計算セクション

### `sp` セクション

一点計算の設定。`mlmm sp` だけが読み込みます。

```yaml
sp:
 hess: false # active-coordinate ONIOM Hessian block も計算
 hessian_calc_mode: FiniteDifference # "FiniteDifference" | "Analytical"
 out_dir: ./result_sp/
```

対応する CLI の `--hess`、`--hessian-calc-mode`、`-o/--out-dir` を明示した場合は CLI が上書きします。

---

## DFT セクション

(ja-dft-section)=
### `dft` セクション

DFT 計算設定。

```yaml
dft:
 func_basis: wb97m-v/def2-svp # 汎関数/基底関数の組み合わせ文字列
 conv_tol: 1.0e-09 # SCF 収束許容値 (Hartree)
 max_cycle: 100 # SCF反復上限
 grid_level: 3 # PySCF グリッドレベル
 engine: gpu # 計算エンジン: "gpu"（gpu4pyscf）または "cpu"（pyscf）。CLI --dft-engine が優先
 ecp: null # ECP 基底名。null の場合は def2-* 基底から自動導出
 lowmem: true # 低メモリの direct JK。false で density fitting
 scf_stepwise_grid: true # 最初のSCFを粗いグリッドから（--scf-stepwise-grid）
 nprocs: auto # scheduler/affinityからPySCF thread数を決定
 memory: auto # host RAM上限（例64GB、GPU VRAMではない）
 verbose: 0 # PySCF 出力詳細レベル; CLI -v 2/3 では実行時 PySCF verbosity が >=4
 out_dir: ./result_dft/ # 出力ディレクトリ
```

**注記:**
- `engine`: `gpu` は gpu4pyscf で実行し、`lowmem: true` の閉殻計算は `rks_lowmem.RKS` を使います。`cpu` は標準の PySCF RKS/UKS を使います。
- `ecp`: 基底名が `def2-` で始まり `ecp` が `null` の場合、ECP として同名の基底が自動的に使用されます。明示的に上書きするには値を設定してください。

---

## スキャン関連セクション

### `bias`

調和バイアス設定。

```yaml
bias:
 k: 300.0 # 調和バイアス強度 (eV/Å^2)
```

---

### `bond`

MLIP ベースの結合変化検出。

```yaml
bond:
 device: auto # MLIP デバイス
 bond_factor: 1.2 # 共有結合半径スケーリング
 margin_fraction: 0.05 # 比較の分率許容値
 delta_fraction: 0.05 # 結合形成/切断を検出する最小相対変化
```

---

## 例: 設定ファイルの全体例

```yaml
# mlmm configuration example

geom:
 coord_type: cart
 freeze_atoms: []

calc:
 model_charge: 0
 model_mult: 1
 backend: uma                  # 高レベルbackend: uma | orb | mace | aimnet2 | dft
 uma_model: uma-s-1p2          # uma-s-1p2 | uma-m-1p1
 ml_device: auto
 hessian_calc_mode: Analytical   # 試しの計算で FiniteDifference と比べる
 mm_device: cpu
 mm_fd: true
 use_bfactor_layers: true # 入力 PDB の B-factor から層を読み取り

gs:
 max_nodes: 20
 climb: true
 climb_lanczos: true

opt:
 thresh: gau
 max_cycles: 100000 # オプティマイザサイクル上限
 dump: false
 out_dir: ./result_all/

stopt:
 thresh: gau_loose
 max_cycles: 300 # ストリング最適化サイクル上限
 lbfgs:
   thresh: gau
   # max_cycles: 100000 # 任意の上書き

bond:
 bond_factor: 1.2
 delta_fraction: 0.05

search:
 max_depth: 10
 max_nodes_segment: 20

freq:
 max_write: 10
 amplitude_ang: 0.8

thermo:
 temperature: 298.15
 pressure_atm: 1.0
 symmetry_number: null

dft:
 func_basis: wb97m-v/def2-svp
 grid_level: 3
```

## 使用上の注意点

- 計算機のセクション名には `calc:` と `mlmm:` のどちらも使えます。両方あるときはキーを合わせて使い、同じキーに異なる値があるとエラーで止まります。
- `opt.lbfgs` / `opt.rfo` は `lbfgs` / `rfo`、`freq.thermo` は `thermo` の別の書き方です。熱化学設定は `all --thermo` でも使用します。
- `mlmm all` は、選択して有効化した stage の section だけを使用します。
- `--show-config`（`scan`・`scan2d`・`scan3d` には無い）は、読み込んだ YAML ファイルと最上位のキーを表示してから、そのまま実行を続けます。

## 関連ドキュメント

- [all](all.md) - 一気通貫ワークフロー
- [opt](opt.md) - 単一構造最適化
- [tsopt](tsopt.md) - 遷移状態最適化
- [path-search](path-search.md) - 再帰的 MEP 探索
- [freq](freq.md) - 振動解析
- [dft](dft.md) - DFT 計算
- [ML/MM 計算機](mlmm-calc.md) - ML/MM の層、リンク原子、ONIOM のエネルギー
- [MLIP バックエンド](backends.md) - MLIP バックエンドの詳細
- [トラブルシューティング](troubleshooting.md) - 実行に失敗したときの対処
