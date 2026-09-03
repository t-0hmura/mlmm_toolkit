# YAML 設定リファレンス

## 概要

`mlmm all` は、選択して有効化した stage の section だけを使用します。

| セクション | 説明 | 使用されるコマンド |
|---------|-------------|---------|
| [`geom`](#geom) | ジオメトリと座標設定 | all, opt, scan, scan2d, scan3d, tsopt, freq, irc, path-opt, path-search |
| [`calc`](#calc) | ML/MM calculatorの設定（別名: `mlmm:`） | all, opt, scan, scan2d, scan3d, tsopt, freq, irc, path-opt, path-search |
| [`opt`](#opt) | 最適化の共通設定 | all, opt, scan, scan2d, scan3d, tsopt, path-opt, path-search |
| [`lbfgs`](#lbfgs) | L-BFGSの設定 | all, opt, scan, scan2d, scan3d, tsopt（マイクロイテレーションの MM 緩和）, path-opt, path-search |
| [`rfo`](#rfo) | RFOの設定 | all, opt |
| [`gs`](#gs) | GSM（Growing String Method）設定 | all, path-opt, path-search |
| [`dmf`](#dmf) | DMF（Direct Max Flux）設定 | all, path-opt, path-search |
| [`irc`](#ja-irc-section) | IRC 積分設定 | all, irc |
| [`freq`](#ja-freq-section) | 振動解析設定 | all, freq |
| [`thermo`](#thermo) | 熱化学設定 | all, freq |
| [`dft`](#ja-dft-section) | DFT 計算設定 | all, dft |
| [`bias`](#bias) | 調和バイアス設定 | all, scan, scan2d, scan3d |
| [`bond`](#bond) | 結合変化検出設定 | all, scan, path-search |
| [`search`](#search) | 再帰的経路探索設定 | all, path-search |
| [`hessian_dimer`](#hessian_dimer) | Hessian・ダイマーTS 最適化 | all, tsopt |
| [`rsirfo`](#rsirfo) | Hessian TS 最適化設定 | all, tsopt |
| [`stopt`](#stopt) | ストリング最適化（StringOptimizer）設定 | all, path-opt, path-search |
| [`microiter`](#microiter) | マイクロイテレーション（MM緩和）設定 | all, opt, tsopt |

---

## 共通セクション

### `geom`

ジオメトリ読み込みと座標系の設定。

```yaml
geom:
 coord_type: cart # 座標タイプ: "cart" (デカルト) または "dlc" (非局在化内部座標)
 freeze_atoms: [] # 1 始まりの凍結原子インデックス
 tr_projection: constrained # 固定の内部 Cartesian PHVA 処理
```

**注記:**
- Frozen 層の原子は力がゼロに設定され、Hessian の対応する列もゼロになります。
- 固定の `tr_projection: constrained` 処理は、凍結 anchor をすべて動かさない
  全系剛体運動だけを除去します。一般的な有効 rank は anchor が
  0/1/2/非共線の 3 個以上のとき 6/3/1/0 で、実用的な ML/MM 境界では
  通常 0 です。全原子凍結は明示的なエラーになります。
- 古い非constrained値は明示的に拒否されます。
- `tr_projection` は `freq`、`irc`、`tsopt`、`opt --flatten` が使う内部の
  凍結境界 PHVA fieldであり、user-selectableな処理ではありません。
  `tsopt --ref-mode` の MEP 接線とは無関係です。
- `irc` では `geom.coord_type` が YAML/CLI マージ後に `cart` へ強制されます。

---

### `calc`

ML/MM calculator（MLIP バックエンド + hessian_ff）の設定。

```yaml
calc:
 input_pdb: null # 入力 PDB ファイルパス (CLI --input から設定)
 real_parm7: null # 全系の Amber parm7 トポロジー (CLI --parm)
 model_pdb: null # ML 領域を定義する PDB (CLI --model-pdb)
 model_charge: 0 # ML 領域の電荷 (CLI -q で上書き)
 model_mult: 1 # ML 領域のスピン多重度 (CLI -m で上書き)
 link_mlmm: null # null: parm7 結合から自動決定; list: 明示上書き
 link_atom_method: scaled    # リンク原子配置: "scaled" (g-factor) または "fixed" (1.09/1.01 Å)
 backend: uma # ML バックエンド: "uma" (デフォルト), "orb", "mace", "aimnet2"
 embedcharge: false # 実験的: MLIP は高コストな xTB 補正、dft は PySCF への直接静電埋込み
 embedcharge_step: 0.001 # MLIP/xTB 補正用の数値 Hessian ステップ（dft では未使用）
 embedcharge_cutoff: 12.0 # ML 領域からの MM 点電荷カットオフ (Å)
 xtb_cmd: xtb # xTB 実行コマンド
 xtb_acc: 0.2 # xTB 精度パラメータ
 xtb_workdir: tmp # xTB 作業ディレクトリ
 xtb_keep_files: false # xTB 一時ファイルを保持
 xtb_ncores: 4 # xTB プロセス数
 uma_model: uma-s-1p2 # UMA モデル名: uma-s-1p2, uma-m-1p1
 uma_task_name: omol # UMA バッチに記録されるタスクタグ (backend=uma 時)
 uma_precision: fp32 # fp32 | fp64 (UMA バックエンドの数値精度)
 orb_model: orb_v3_conservative_omol  # ORB モデル名 (backend=orb 時)
 orb_precision: float64  # ORB 浮動小数点精度のデフォルト (backend=orb 時; "float32-high" は TF32 matmul で --precision fp32 でも選択可、レガシー "float32" alias は受理)
 mace_model: MACE-OMOL-0 # MACE モデル名 (backend=mace 時)
 mace_dtype: float64      # MACE 浮動小数点精度 (backend=mace 時)
 aimnet2_model: aimnet2   # AIMNet2 モデル名 (backend=aimnet2 時)
 hessian_calc_mode: FiniteDifference # ML Hessianモード: "FiniteDifference" または "Analytical"
 out_hess_torch: true # Hessianを torch.Tensor で返す
 H_double: true # Hessianを float64 で組み立て・返却
 ml_device: auto # ML デバイス: "cuda", "cpu", "auto"
 ml_cuda_idx: 0 # CUDA デバイスインデックス
 mm_backend: hessian_ff # MM バックエンド: "hessian_ff" (解析的) | "openmm" (FD Hessian)
 use_cmap: true         # parm7 の CMAP を REAL と MODEL の両 MM 層で保持
 mm_device: cpu # MM デバイス (hessian_ff は CPU のみ、OpenMM は CUDA/CPU 対応)
 mm_cuda_idx: 0 # MM CUDA インデックス (OpenMM のみ)
 mm_threads: 16 # MM 計算のスレッド数
 workers: 1 # ローカル ML worker process 数（対応 backend のみ）
 workers_per_node: null # 未設定。UMA parallel predictor 使用時の実効値は 1
 mm_fd: true # MM Hessianに有限差分を使用
 mm_hessian_mode: null # 明示指定は finite_difference/analytical。null は mm_fd に従う
 mm_fd_dir: null # MM Hessianログの出力ディレクトリ
 mm_fd_delta: 0.001 # 有限差分ステップ（保持）
 symmetrize_hessian: true # 最終Hessianを 0.5*(H+H^T) で対称化
 print_timing: true # ML/MM Hessianのタイミング内訳を表示
 print_vram: true # CUDA VRAM 使用量を表示
 return_partial_hessian: true # アクティブブロック部分Hessian（CLI ラッパー側で true デフォルトを適用）
 freeze_atoms: [] # geom.freeze_atoms から継承
 # 層設定:
 hess_cutoff: null # Å: null = 可動 MM をすべて Hessian 対象に含める (デフォルト)、>0.0 で ML 周辺の指定距離内 MM のみに限定
 movable_cutoff: null # Å: movable MM の距離カットオフ
 use_bfactor_layers: true # 入力 PDB の B-factor から層を読み取り
 hess_mm_atoms: null # 明示的 Hessian 対象 MM 原子インデックス (1始まり)
 movable_mm_atoms: null # 明示的 movable MM 原子インデックス (1始まり)
 frozen_mm_atoms: null # 明示的 frozen MM 原子インデックス (1始まり)
```

**注記:**
- セクション名は `calc:` が正式名で、`mlmm:` は互換用の別名として受け付けます（`opt`、`sp`、`tsopt`、`freq`、`irc`、`dft`、`path-opt`、`path-search`、`scan`、`scan2d`、`scan3d` で認識）。両方が存在する場合は `calc:` が優先されます。
- `backend`: ML バックエンドを選択します。`uma`（デフォルト）、`orb`、`mace`、`aimnet2` から選択可能です。UMA 以外のバックエンドを使用するには、対応するオプション依存パッケージのインストールが必要です（例: `pip install "mlmm-toolkit[orb]"`）。
- バックエンド固有のモデルキーは、対応するバックエンドが選択されている場合にのみ有効です:
  - `uma_model`、`uma_task_name` — UMA バックエンドのみ
  - `orb_model`、`orb_precision` — ORB バックエンドのみ
  - `mace_model`、`mace_dtype` — MACE バックエンドのみ
  - `aimnet2_model` — AIMNet2 バックエンドのみ
- `hessian_calc_mode: Analytical` はバックエンドの解析 Hessian を明示的に要求します。UMA、ORB、MACE、AIMNet2 がこの経路を実装しており、インストール済みバックエンドが非対応なら計算法を暗黙に変更せずエラーになります。`workers > 1` との併用もエラーです。
- `hess_cutoff` のデフォルト `null` は可動 MM 原子をすべて Hessian 対象に含めることを意味します（freq/irc/opt はすべての可動原子を解析します）。値（>0.0）を指定すると、その距離以内の MM 原子のみに Hessian 対象を限定します。`movable_cutoff` を指定しない場合は `freeze_atoms` の指定に従います。
- `use_bfactor_layers: true` を設定すると、`define-layer` で書き込んだ B-factor から層割り当てを読み取ります。
- 明示的インデックス（`hess_mm_atoms` 等）が設定された場合、カットオフや B-factor よりも優先されます。
- `opt`/`tsopt`/`irc`/`freq` は、YAML で `calc.return_partial_hessian` を明示しない場合に部分 Hessian をデフォルトで使用します。
- これらのコマンドで完全 Hessian を強制するには `calc.return_partial_hessian: false` を明示してください。
- `mm_fd: true` は有限差分 MM Hessian、`false` は `hessian_ff` の解析
  MM Hessian を使います。`mm_hessian_mode` は
  `finite_difference`/`analytical` の明示形で、`null` の場合は互換用の
  `mm_fd` に従います。
- `use_cmap: true`（デフォルト）は parm7 に含まれる CMAP を REAL と MODEL の両 MM 層で保持します。明示的な改変力場計算だけ `false` を指定してください。この場合は両層から CMAP を除去します。
- standalone ML/MM 計算には `real_parm7` が必須です。ML 領域は `model_pdb`、明示的な model index、または有効な B-factor layer から指定できます。
- `irc` は YAML の設定にかかわらず `geom.coord_type = cart` を強制します。

---

### `opt`

L-BFGS/RFO で共通の最適化設定。ここに書いた全キーが `opt` コマンドの optimizer に
届きます（`--microiter` の有無に関わらず。microiteration の macro step がその
optimizer です）。`tsopt` の macro optimizer にも届き、その上に
[`rsirfo`](#rsirfo) / [`hessian_dimer`](#hessian_dimer) が重なります。
転送されるのは**実際に変更した値だけ**なので、触っていないキーについては
optimizer 固有セクションが優先されます。

```yaml
opt:
 thresh: gau # 収束プリセット: gau_loose, gau, gau_tight, gau_vtight, baker, never
 align: false # StringOptimizer 専用: alignment の有効/無効
 max_cycles: 100000 # オプティマイザサイクル上限
 print_every: 100 # ログ出力間隔
 min_step_norm: 1.0e-08 # 最小ステップノルム
 assert_min_step: true # ステップが閾値以下で停止
 rms_force: null # 明示的 RMS 力ターゲット
 rms_force_only: false # RMS 力のみで収束判定
 max_force_only: false # 最大力のみで収束判定
 force_only: false # 変位チェックをスキップ
 converge_to_geom_rms_thresh: 0.05 # 参照ジオメトリへの収束 RMS 閾値
 overachieve_factor: 0.0 # 閾値の引き締め係数
 check_eigval_structure: false # Hessian固有値構造の検証
 energy_plateau: false # opt-in（--stop-plateau）: エネルギーが停滞したら stalled として停止 (収束扱いにはしない)
 energy_plateau_thresh: 1.0e-4 # エネルギー変動許容幅 au（約 0.06 kcal/mol）
 energy_plateau_window: 50 # プラトー判定に用いる直近ステップ数
 line_search: true # ラインサーチを有効化
 dump: false # 軌跡/リスタートデータの出力
 dump_restart: false # リスタートチェックポイントの出力
 reparam_thresh: 0.0 # StringOptimizer 専用: 再パラメータ化閾値
 coord_diff_thresh: 0.0 # StringOptimizer 専用: 座標差分閾値
 prefix: "" # ファイル名プレフィックス
 out_dir: ./result_opt/ # 出力ディレクトリ
```

**収束プリセット:**

| プリセット | Max Force | RMS Force | Max Step | RMS Step |
|-----------|-----------|-----------|----------|----------|
| `gau_loose` | 2.5e-3 | 1.7e-3 | 1.0e-2 | 6.7e-3 |
| `gau` | 4.5e-4 | 3.0e-4 | 1.8e-3 | 1.2e-3 |
| `gau_tight` | 1.5e-5 | 1.0e-5 | 6.0e-5 | 4.0e-5 |
| `gau_vtight` | 2.0e-6 | 1.0e-6 | 6.0e-6 | 4.0e-6 |
| `baker` | 3.0e-4 | 2.0e-4 | 3.0e-4 | 2.0e-4 |

`baker` は5基準すべてを要求します。`|delta E| < 1e-6` に加え、max/RMS force と max/RMS step がすべて閾値を満たす必要があります。

**エネルギープラトー停止（opt-in、デフォルト無効）:**

`energy_plateau` のデフォルトは `false` です。`opt` / `tsopt` / `all` の
`--stop-plateau` で有効化し、`--stop-plateau-thresh` / `--stop-plateau-window` が
上記の 2 つの値を設定します。有効時、直近 `energy_plateau_window` ステップ
（デフォルト 50）のエネルギー範囲 `max(E) - min(E)` が `energy_plateau_thresh`
（デフォルト `1.0e-4` au、約 0.06 kcal/mol）を下回ると、オプティマイザを `status: "stalled"` で停止します（`converged` とは
区別される非収束の結果で、決して `converged` にはなりません）。

これは ML/MM 最適化で cycle を節約するための機構です。MLIP の力には数値精度に起因する
ノイズフロアがあり、これが `gau`/`baker` などの勾配ベース収束閾値を
上回ると、ジオメトリが実質的に停止していても力が閾値を下回らないことがあります。
エネルギー自体が MLIP の数値精度内で平坦化した段階では、追加ステップを回しても
残差力はノイズフロア以下にならない場合があります。ただしエネルギーの平坦化は
停留点の証拠ではないため、この停止は明示的に指定したときだけ働き、
実行の実質的な上限は常に `max_cycles` です。

Chain-of-states（COS）最適化（GS/DMF ストリング最適化等）では、
プラトー判定は自動的にスキップされます。`--microiter` の **MM micro 反復** でも
常にスキップされます（力が閾値を超えたまま MM エネルギーが平坦なのは MM 平衡ではなく
停滞した micro 緩和であり、そこで止めると周辺環境が緩和されないまま macro/micro
交互計算が終了してしまうためです）。micro 側の上限は `microiter.micro_max_cycles` で、これは予定回数ではなく backstop です（micro は収束で抜けます。リリーススモークの ML/MM TS レーンでは中央値 56 cycle、791 回中3回の過渡だけが 10^4 台を要しました）。

---

### `lbfgs`

L-BFGSの設定（`opt` を拡張）。

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
 reject_uphill: false # 許容値を超えるenergy上昇の拒否を明示的に有効化
 uphill_tolerance: 0.0001 # energy上昇許容値（Hartree）
 rejection_step_floor: 1.0e-07 # retry stepの下限
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
 reject_uphill: false # 許容値を超えるenergy上昇の拒否を明示的に有効化
 uphill_tolerance: 0.0001 # energy上昇許容値（Hartree）
 rejection_trust_floor: 1.0e-07 # retry trust radiusの下限
 max_rejections_at_floor: 3 # 下限での連続拒否後に停止
 hessian_update: bfgs # Hessian更新スキーム: bfgs, bofill 等
 hessian_init: calc # Hessian初期化: calc, unit 等
 hessian_recalc: 500 # N ステップごとにHessianを再構築
 hessian_recalc_adapt: null # 適応的Hessian再構築係数
 small_eigval_thresh: 1.0e-08 # 安定性のための固有値閾値
 alpha0: 1.0 # 初期マイクロステップ
 max_micro_cycles: 50 # マイクロイテレーションの上限
 rfo_overlaps: false # RFO オーバーラップを有効化
 gediis: false # GEDIIS を有効化
 gdiis: true # GDIIS を有効化
 gdiis_thresh: 0.0025 # GDIIS 受容閾値
 gediis_thresh: 0.01 # GEDIIS 受容閾値
 gdiis_test_direction: true # DIIS 前に降下方向をテスト
 adapt_step_func: true # 適応的ステップスケーリング
```

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
 max_micro_cycles: 10 # マイクロイテレーションの上限
 reset_dlc: true # 各ステップで非局在化座標を再構築
 climb: true # クライミングイメージを有効化
 climb_rms: 0.0005 # クライミング RMS 閾値
 climb_lanczos: true # クライミングの Lanczos 精密化
 climb_lanczos_rms: 0.0005 # Lanczos RMS 閾値
 climb_fixed: false # クライミングイメージを固定
 scheduler: null # オプションのスケジューラバックエンド
```

`gs.param` は `equi` または `energy` を受け付けます。energy weighting はGSMストリングの完全成長後にのみ適用され、高エネルギー領域へノード密度を寄せます。対応するCLIオプションは `--gsm-param` です。

---

### `dmf`

Direct Max Flux（DMF）による MEP 最適化。

```yaml
dmf:
 max_cycles: 3000 # DMF/IPOPT反復上限
 tol: tight # IPOPT dual_inf_tol: tight(0.04) | middle(0.10) | loose(0.20) または正の float（--thresh-dmf で上書き）
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
 kink_max_nodes: 3 # ねじれ最適化の最大ノード数
 max_seq_kink: 2 # 連続ねじれの上限
 refine_mode: null # 精密化戦略: peak, minima, null (自動)
```

---

## TS 最適化セクション

### `hessian_dimer`

Hessian・ダイマー TS 最適化（`tsopt --opt-mode grad`）。`opt.thresh` と `hessian_dimer.thresh` を両方明示する場合は同じ値にしてください。片方だけならその値、無指定なら Dimer のデフォルトを使います。

```yaml
hessian_dimer:
 thresh_loose: gau_loose # 緩い収束プリセット
 thresh: baker # メイン収束プリセット
 update_interval_hessian: 500 # Hessian再構築間隔
 flatten_amp_ang: 0.1 # フラット化振幅 (Å)
 flatten_max_iter: 50 # フラット化反復上限（デフォルト 50、--no-flatten で 0 に設定）
 flatten_sep_cutoff: 0.0 # 代表原子間の最小距離
 flatten_k: 10 # モードあたりのサンプル代表原子数
 flatten_loop_bofill: false # フラット化変位に Bofill 更新
 mem: 100000 # ソルバーのメモリ上限
 device: auto # 固有値ソルバーのデバイス選択
 root: 0 # ターゲット TS ルートインデックス
 partial_hessian_flatten: true # 部分Hessianを虚モード検出に使用
 ml_only_hessian_dimer: false # ダイマー方向決定に ML 領域のみのHessianを使用
 dimer:
 length: 0.0189 # ダイマー間隔 (Bohr)
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
- 通常の TSOPT 省略時は flattening が無効で実効反復数は 0 です。有効化した場合の上限を `flatten_max_iter` が制御し、そのデフォルトは 50 です。
- CLI フラグ `--flatten` / `--no-flatten`（`tsopt` および `all`）はこの設定と連動します。`--flatten` はデフォルトの `flatten_max_iter`（50）でflatteningループを有効化し、`--no-flatten` は `flatten_max_iter` を 0 に強制してループを無効化します。`--flatten` と同時に YAML で `flatten_max_iter` を明示指定した場合は、YAML の値が優先されます。
- 内側の L-BFGS 固有設定は、最上位の `lbfgs` ではなく `hessian_dimer.lbfgs` に置きます。共通の `print_every` と `energy_plateau*` は上記の競合規則に従います。`line_search` は `false` 固定で、Dimer の有効力は表示する物理エネルギーの勾配ではないため `true` は拒否されます。`max_cycles` は設定できず、各 segment には `opt.max_cycles` の残り cycle 数が渡されます。

---

### `rsirfo`

Hessian TS 最適化の共通設定です。デフォルトの RS-P-RFO
（`tsopt --opt-mode hess` / `rsprfo`）と、明示的な `rsirfo` / `trim` に適用されます。

```yaml
rsirfo:
 thresh: baker # Hessian TS 収束プリセット
 max_cycles: 100000 # opt.max_cycles と共有するサイクル上限
 print_every: 100 # ログ出力間隔
 min_step_norm: 1.0e-08 # 最小ステップノルム
 assert_min_step: true # ステップ停滞時にアサート
 roots: [0] # 追跡するrootは1個のみ（空list・複数rootは拒否）
 hessian_ref: null # 参照Hessian
 rx_modes: null # 反応モード定義
 prim_coord: null # 監視する主座標
 rx_coords: null # 監視する反応座標
 hessian_update: bofill # Hessian更新スキーム
 hessian_recalc_reset: true # 正確なHessian後に再計算カウンタをリセット
 hessian_init: calc # Hessian初期化
 hessian_recalc: 500 # Hessian再構築間隔
 max_micro_cycles: 50 # マクロサイクルあたりのマイクロイテレーション数
 augment_bonds: false # 結合解析に基づく反応経路の拡張
 min_line_search: false # RS-P-RFO のみ: 最小化部分空間で補間
 max_line_search: false # RS-P-RFO のみ: 最大化部分空間で補間
 assert_neg_eigval: false # 収束時に負の固有値を要求
 track_mode_by_overlap: false # mlmm 固有: オーバーラップでターゲットモードを追跡
 trust_radius: 0.10 # 信頼領域半径
 trust_update: true # 信頼領域更新
 trust_min: 0.0001 # 最小信頼半径
 trust_max: 0.10 # 最大信頼半径（ML/MM 安定性のため調整）
 small_eigval_thresh: 1.0e-08 # 安定性のための固有値閾値
 out_dir: ./result_tsopt/ # 出力ディレクトリ
```

`min_line_search` と `max_line_search` を使用するのは
`--opt-mode rsprfo` だけで、YAML の明示値を反映します。

`opt` と `rsirfo` に同じ設定を明示する場合は値を一致させてください。片方だけならその値、無指定なら `rsirfo` のデフォルトを使います。

---

### `stopt`

ストリング最適化（GS/DMF）の設定。path-opt と path-search で使用。

```yaml
stopt:
 type: string           # 最適化タイプラベル（StringOptimizer用）
 thresh: gau_loose      # ストリング最適化の収束プリセット（--thresh-gsm で上書き）
 stop_in_when_full: 300 # ストリングが満杯時の早期停止閾値
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
   # 単一構造最適化用（HEI±1、ねじれノード）
   thresh: gau
   # max_cycles: 100000 # 任意の上書き
   # ...（詳細は lbfgs セクション参照）
```

**注意:**
- `stopt.lbfgs` は HEI±1 端点最適化およびねじれノード最適化に使用される単一構造最適化（L-BFGS）の設定です。この入れ子レベルでは L-BFGS のみが参照されるため、`stopt.rfo:` ブロックは無視されます。
- 外側の `stopt` キーはストリング最適化（GS または DMF ラッパー）を制御します。

---

## IRC セクション

(ja-irc-section)=
### `irc` (section)

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
 energy_increase_thresh: 0.0   # 通常modeでは1 stepでもenergyが上昇すれば停止
 imag_below: 0.0 # 虚振動数カットオフ
 force_inflection: true # 変曲点検出の強制
 check_bonds: false # 伝搬中の結合チェック
 out_dir: ./result_irc/ # 出力ディレクトリ
 prefix: "" # ファイル名プレフィックス
 dump_fn: irc_data.h5 # IRC データファイル名
 dump_every: null # デフォルトでは無効。有効化する場合のみ正の間隔を指定
 max_pred_steps: 500 # 予測子-修正子の最大ステップ数
 loose_cycles: 3 # 引き締め前の緩いサイクル数
 corr_func: mbs # 相関関数の選択
```

---

## 振動解析セクション

(ja-freq-section)=
### `freq` (section)

振動解析設定。

```yaml
freq:
 active_dof_mode: partial # アクティブ原子の選択: "all" | "ml-only" | "partial" | "unfrozen"
 zero_cutoff_cm: 5.0 # |振動数| がこの値以下のモードを除外（cm^-1）
 amplitude_ang: 0.8 # モード変位振幅 (Å)
 n_frames: 20 # モードtrajectoryのフレーム数
 max_write: 10 # 書き出すモードの最大数
 sort: value # ソート順: "value" または "abs"
 out_dir: ./result_freq/ # 出力ディレクトリ
```

`freq.zero_cutoff_cm` は standalone `freq`、`opt` flatten、Dimer、
Hessian系TS最適化が共有します。旧`hessian_dimer.neg_freq_thresh_cm` と
`rsirfo.saddle_imaginary_threshold_cm` は互換aliasですが、競合する値は
エラーになります。

**注記:**
- `active_dof_mode`: 振動解析に参加させる原子集合を選択します。`all` は全原子、`ml-only` は ML 領域のみ、`partial`（デフォルト）は ML + Movable-MM、`unfrozen` は凍結されていない全原子を使用します。CLI フラグ `--active-dof-mode` が明示された場合は YAML 値より優先されます。

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

### `sp` (section)

single-point 設定。`mlmm sp` だけが読み込みます。

```yaml
sp:
 hess: false # active-coordinate ONIOM Hessian block も計算
 hessian_calc_mode: FiniteDifference # "FiniteDifference" | "Analytical"
 out_dir: ./result_sp/
```

対応する CLI の `--hess`、`--hessian-calc-mode`、`-o/--out-dir` を
明示した場合は CLI が上書きします。

---

### `microiter`

ML/MM最適化用のマイクロイテレーション設定。`--microiter` 有効時、MLリージョンの
マクロステップ間でMMリージョンをL-BFGSで緩和（ML原子は凍結）します。

```yaml
microiter:
 micro_thresh: null       # MM緩和の収束プリセット（L-BFGS）; null → マクロステップと同じ
 micro_max_cycles: 100000 # マイクロイテレーションサイクル上限
```

**注意:**
- CLIフラグ `--microiter` / `--no-microiter` で有効化（デフォルト: 有効）
- `opt --opt-mode hess` と、すべての Hessian TS mode（`hess`, `rsirfo`, `rsprfo`, `trim`）で使用可能
- `micro_thresh` は `opt.thresh` と同じプリセット（gau_loose, gau, gau_tight等）を受け付けます。`null` または省略時はマクロステップの閾値と同じになります

---

## DFT セクション

(ja-dft-section)=
### `dft` (section)

DFT 計算設定。

```yaml
dft:
 func_basis: wb97m-v/def2-tzvpd # 汎関数/基底関数の組み合わせ文字列
 conv_tol: 1.0e-09 # SCF 収束許容値 (Hartree)
 max_cycle: 100 # SCF反復上限
 grid_level: 3 # PySCF グリッドレベル
 engine: gpu # 計算エンジン: "gpu"（gpu4pyscf）または "cpu"（pyscf）。CLI --engine が優先
 ecp: null # ECP 基底名。null の場合は def2-* 基底から自動導出
 lowmem: true # closed-shell GPU で gpu4pyscf rks_lowmem.RKS を使用
 verbose: 0 # PySCF 出力詳細レベル; CLI -v 2/3 では実行時 PySCF verbosity が >=4
 out_dir: ./result_dft/ # 出力ディレクトリ
```

**注記:**
- `engine`: `gpu` は gpu4pyscf 経由で実行（closed-shell かつ `lowmem: true` のとき `rks_lowmem.RKS` を使用）。`cpu` は標準 PySCF の RKS/UKS にフォールバックします。CLI フラグ `--engine` が明示された場合は YAML 値より優先されます。
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

## 例: 複数セクションを含む設定ファイル

```yaml
# mlmm configuration example

geom:
 coord_type: cart
 freeze_atoms: []
 tr_projection: constrained

calc:
 model_charge: 0
 model_mult: 1
 backend: uma                  # ML バックエンド: uma | orb | mace | aimnet2
 uma_model: uma-s-1p2          # uma-s-1p2 | uma-m-1p1
 ml_device: auto
 hessian_calc_mode: Analytical   # 代表的な pilot で FiniteDifference と比較
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
 func_basis: wb97m-v/def2-tzvpd
 grid_level: 3
```

---

## 参照

- [all](all.md) - メインワークフロー
- [opt](opt.md) - 単一構造最適化
- [tsopt](tsopt.md) - 遷移状態最適化
- [path-search](path-search.md) - 再帰的 MEP 探索
- [freq](freq.md) - 振動解析
- [dft](dft.md) - DFT 計算
- [概念とワークフロー](concepts.md) - ML/MM 3層システムと ONIOM エネルギー分解
- [ML/MM calculator](mlmm-calc.md) - ML/MM calculatorの詳細
