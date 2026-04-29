# YAML 設定リファレンス

## 概要

| セクション | 説明 | 使用されるコマンド |
|---------|-------------|---------|
| [`geom`](#geom) | ジオメトリと座標設定 | all, opt, scan, scan2d, scan3d, tsopt, freq, irc, path-opt, path-search |
| [`calc`](#calc) | ML/MM 計算機の設定 | all, opt, scan, scan2d, scan3d, tsopt, freq, irc, path-opt, path-search |
| [`opt`](#opt) | 最適化の共通設定 | opt, scan, scan2d, scan3d, tsopt, path-opt, path-search |
| [`lbfgs`](#lbfgs) | L-BFGSの設定 | opt, scan, scan2d, scan3d, path-search |
| [`rfo`](#rfo) | RFOの設定 | opt, scan, scan2d, scan3d, path-search |
| [`gs`](#gs) | GSM（Growing String Method）設定 | path-opt, path-search |
| [`dmf`](#dmf) | DMF（Direct Max Flux）設定 | path-opt, path-search |
| [`irc`](#ja-irc-section) | IRC積分設定 | irc |
| [`freq`](#ja-freq-section) | 振動解析設定 | freq |
| [`thermo`](#thermo) | 熱化学設定 | freq |
| [`dft`](#ja-dft-section) | DFT計算設定 | dft |
| [`bias`](#bias) | 調和バイアス設定 | scan, scan2d, scan3d |
| [`bond`](#bond) | 結合変化検出設定 | scan, path-search |
| [`search`](#search) | 再帰的経路探索設定 | path-search |
| [`hessian_dimer`](#hessian_dimer) | ヘシアン・ダイマーTS 最適化 | tsopt |
| [`rsirfo`](#rsirfo) | RS-I-RFO TS 最適化 | tsopt |
| [`stopt`](#stopt) | ストリング最適化（StringOptimizer）設定 | path-opt, path-search |
| [`microiter`](#microiter) | マイクロイテレーション（MM緩和）設定 | opt, tsopt |

---

## 共通セクション

### `geom`

ジオメトリ読み込みと座標系の設定。

```yaml
geom:
 coord_type: cart # 座標タイプ: "cart" (デカルト) または "dlc" (非局在化内部座標)
```

**注記:**
- ML/MM モードでは、Frozen 層の原子が自動的に `freeze_atoms` に追加されます。
- `irc` では `geom.coord_type` が YAML/CLI マージ後に `cart` へ強制されます。

---

### `calc`

ML/MM 計算機（MLIP バックエンド + hessian_ff）の設定。

```yaml
calc:
 input_pdb: null # 入力 PDB ファイルパス (CLI --input から設定)
 real_parm7: null # 全系の Amber parm7 トポロジー (CLI --parm)
 model_pdb: null # ML 領域を定義する PDB (CLI --model-pdb)
 model_charge: 0 # ML 領域の電荷 (CLI -q で上書き)
 model_mult: 1 # ML 領域のスピン多重度 (CLI -m で上書き)
 link_mlmm: null # リンク原子ペアの明示指定 (null で自動検出)
 link_atom_method: scaled    # リンク原子配置: "scaled" (g-factor) または "fixed" (1.09/1.01 Å)
 backend: uma # ML バックエンド: "uma" (デフォルト), "orb", "mace", "aimnet2"
 embedcharge: false # xTB 点電荷埋め込み補正 (CLI --embedcharge で有効化)
 embedcharge_step: 0.001 # 埋め込み補正の数値ヘシアンステップ (Å)
 embedcharge_cutoff: 12.0 # xTB 埋め込み用 MM 点電荷のカットオフ距離 (Å)
 xtb_cmd: xtb # xTB 実行コマンド
 xtb_acc: 0.2 # xTB 精度パラメータ
 xtb_workdir: tmp # xTB 作業ディレクトリ
 xtb_keep_files: false # xTB 一時ファイルを保持
 xtb_ncores: 4 # xTB のコア数
 uma_model: uma-s-1p1 # UMA モデル名: uma-s-1p1, uma-m-1p1
 uma_task_name: omol # UMA バッチに記録されるタスクタグ (backend=uma 時)
 orb_model: orb_v3_conservative_omol  # ORB モデル名 (backend=orb 時)
 orb_precision: float32  # ORB 浮動小数点精度 (backend=orb 時)
 mace_model: MACE-OMOL-0 # MACE モデル名 (backend=mace 時)
 mace_dtype: float64      # MACE 浮動小数点精度 (backend=mace 時)
 aimnet2_model: aimnet2   # AIMNet2 モデル名 (backend=aimnet2 時)
 hessian_calc_mode: FiniteDifference # ML ヘシアンモード: "FiniteDifference" または "Analytical"
 out_hess_torch: true # ヘシアンを torch.Tensor で返す
 H_double: true # ヘシアンを float64 で組み立て・返却
 ml_device: auto # ML デバイス: "cuda", "cpu", "auto"
 ml_cuda_idx: 0 # CUDA デバイスインデックス
 mm_backend: hessian_ff # MM バックエンド: "hessian_ff" (解析的) | "openmm" (FD ヘシアン)
 use_cmap: false        # true にすると model parm7 に CMAP 項を含める。false (デフォルト) は Gaussian ONIOM と同一
 mm_device: cpu # MM デバイス (hessian_ff は CPU のみ、OpenMM は CUDA/CPU 対応)
 mm_cuda_idx: 0 # MM CUDA インデックス (OpenMM のみ)
 mm_threads: 16 # MM 計算のスレッド数
 mm_fd: true # MM ヘシアンに有限差分を使用
 mm_fd_dir: null # MM ヘシアンログの出力ディレクトリ
 mm_fd_delta: 0.001 # 有限差分ステップ（保持）
 symmetrize_hessian: true # 最終ヘシアンを 0.5*(H+H^T) で対称化
 print_timing: true # ML/MM ヘシアンのタイミング内訳を表示
 print_vram: true # CUDA VRAM 使用量を表示
 return_partial_hessian: true # アクティブブロック部分ヘシアン（CLI ラッパー側で true 既定を適用）
 freeze_atoms: [] # geom.freeze_atoms から継承
 # 層設定:
 hess_cutoff: null # Å: Hessian 対象 MM の距離カットオフ
 movable_cutoff: null # Å: movable MM の距離カットオフ
 use_bfactor_layers: true # 入力 PDB の B-factor から層を読み取り
 hess_mm_atoms: null # 明示的 Hessian 対象 MM 原子インデックス (1始まり)
 movable_mm_atoms: null # 明示的 movable MM 原子インデックス (1始まり)
 frozen_mm_atoms: null # 明示的 frozen MM 原子インデックス (1始まり)
```

**注記:**
- `backend`: ML バックエンドを選択します。`uma`（デフォルト）、`orb`、`mace`、`aimnet2` から選択可能です。UMA 以外のバックエンドを使用するには、対応するオプション依存パッケージのインストールが必要です（例: `pip install "mlmm-toolkit[orb]"`）。
- バックエンド固有のモデルキーは、対応するバックエンドが選択されている場合にのみ有効です:
  - `uma_model`、`uma_task_name` — UMA バックエンドのみ
  - `orb_model`、`orb_precision` — ORB バックエンドのみ
  - `mace_model`、`mace_dtype` — MACE バックエンドのみ
  - `aimnet2_model` — AIMNet2 バックエンドのみ
- `embedcharge`: `true` に設定すると、xTB 点電荷埋め込み補正が有効化されます。MM 領域の部分電荷を点電荷として ML 計算に埋め込み、MM 環境から ML 領域への静電的影響（分極効果）を考慮します。デフォルトは `false` です。`$PATH` 上に `xtb` 実行ファイルが必要です。
- `xtb_cmd`、`xtb_acc`、`xtb_ncores`、`xtb_workdir`、`xtb_keep_files` は `embedcharge` が有効な場合に xTB サブプロセスを設定します。
- `hessian_calc_mode: Analytical` が推奨です（VRAM に余裕がある場合、ML 原子 300 以上では 24 GB 以上推奨）。UMA バックエンドでのみ利用可能で、他のバックエンドでは自動的に `FiniteDifference` が使用されます。- `hess_cutoff`/`movable_cutoff` を指定しない場合、ML 以外の全原子が Hessian 対象 MM に分類されます。
- `use_bfactor_layers: true` を設定すると、`define-layer` で書き込んだ B-factor から層割り当てを読み取ります。
- 明示的インデックス（`hess_mm_atoms` 等）が設定された場合、カットオフや B-factor よりも優先されます。
- `opt`/`tsopt`/`irc`/`freq` は、YAML で `calc.return_partial_hessian` を明示しない場合に部分ヘシアンを既定で使用します。
- これらのコマンドで完全ヘシアンを強制するには `calc.return_partial_hessian: false` を明示してください。
- `mm_fd: true` は MM ヘシアンに有限差分を使用します。解析的 MM ヘシアン（hessian_ff）を使用するには `false` に設定してください。
- `use_cmap: false`（デフォルト）は model parm7 から CMAP 項（骨格クロスマップ二面角補正）を除外します。これは Gaussian ONIOM の挙動（CMAP を model MM に含めない）と一致します。`true` に設定すると model parm7 に CMAP が含まれ、ONIOM 差し引きで骨格 CMAP がキャンセルされます。ML 領域に骨格原子を含まない典型的な活性部位モデルでは、どちらの設定も実質的に同じ結果になります。
- `real_parm7` と `model_pdb` は ML/MM 計算に必須です。
- `irc` は YAML の設定にかかわらず `geom.coord_type = cart` を強制します。

---

### `opt`

L-BFGS/RFO で共通の最適化設定。

```yaml
opt:
 type: string # StringOptimizer 専用: optimizer type label
 thresh: gau # 収束プリセット: gau_loose, gau, gau_tight, gau_vtight, baker, never
 stop_in_when_full: 300 # StringOptimizer 専用: string 完了時の早期停止閾値
 align: false # StringOptimizer 専用: alignment の有効/無効
 scale_step: global # StringOptimizer 専用: step scaling モード
 max_cycles: 10000 # 最大反復回数
 print_every: 100 # ログ出力間隔
 min_step_norm: 1.0e-08 # 最小ステップノルム
 assert_min_step: true # ステップが閾値以下で停止
 rms_force: null # 明示的 RMS 力ターゲット
 rms_force_only: false # RMS 力のみで収束判定
 max_force_only: false # 最大力のみで収束判定
 force_only: false # 変位チェックをスキップ
 converge_to_geom_rms_thresh: 0.05 # 参照ジオメトリへの収束 RMS 閾値
 overachieve_factor: 0.0 # 閾値の引き締め係数
 check_eigval_structure: false # ヘシアン固有値構造の検証
 energy_plateau: true # フォールバック: エネルギーが停滞したら収束と判定
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

**エネルギープラトー・フォールバック:**

`energy_plateau: true`（デフォルト）の場合、直近 `energy_plateau_window` ステップ
（デフォルト 50）のエネルギー範囲 `max(E) - min(E)` が `energy_plateau_thresh`
（デフォルト `1.0e-4` au、約 0.06 kcal/mol）を下回ると、オプティマイザは収束したと
判定します。

これは ML/MM 最適化特有の安全網です。MLIP の力には数値精度に起因する
ノイズフロアがあり、これが `gau`/`baker` などの勾配ベース収束閾値を
上回ると、ジオメトリが実質的に停止していても力が閾値を下回らず、
最適化が延々と回り続けることがあります。エネルギー自体が MLIP の
数値精度内で平坦化した段階では、追加ステップを回しても残差力は
ノイズフロア以下にならないため、プラトー判定によってクリーンに
終了させます。

フォールバックを無効化するには `energy_plateau: false` を設定します
（この場合は `thresh` プリセットのみで収束判定されます）。
Chain-of-states（COS）最適化（GS/DMF ストリング最適化等）では、
プラトー判定は自動的にスキップされます。

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
```

---

### `rfo`

RFO（Rational Function Optimizer）の設定（`opt` を拡張）。

```yaml
rfo:
 trust_radius: 0.10 # 信頼領域半径
 trust_update: true # 信頼領域更新を有効化
 trust_min: 0.0001 # 最小信頼半径
 trust_max: 0.10 # 最大信頼半径（v0.2.8 で ML/MM 安定性のため厳格化）
 max_energy_incr: null # ステップあたりの許容エネルギー増加
 hessian_update: bfgs # ヘシアン更新スキーム: bfgs, bofill 等
 hessian_init: calc # ヘシアン初期化: calc, unit 等
 hessian_recalc: 500 # N ステップごとにヘシアンを再構築
 hessian_recalc_adapt: null # 適応的ヘシアン再構築係数
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
 max_micro_cycles: 10 # マイクロ反復の上限
 reset_dlc: true # 各ステップで非局在化座標を再構築
 climb: true # クライミングイメージを有効化
 climb_rms: 0.0005 # クライミング RMS 閾値
 climb_lanczos: true # クライミングの Lanczos 精密化
 climb_lanczos_rms: 0.0005 # Lanczos RMS 閾値
 climb_fixed: false # クライミングイメージを固定
 scheduler: null # オプションのスケジューラバックエンド
```

---

### `dmf`

Direct Max Flux（DMF）による MEP 最適化。

```yaml
dmf:
 max_cycles: 300 # DMF/IPOPT の最大反復数（--max-cycles で上書き）
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

---

### `search`

再帰的経路探索（path-search のみ）。

```yaml
search:
 max_depth: 10 # 再帰深度の上限
 stitch_rmsd_thresh: 0.0001 # セグメント縫合の RMSD 閾値
 bridge_rmsd_thresh: 0.0001 # ブリッジノードの RMSD 閾値
 max_nodes_segment: 10 # セグメントあたりの最大ノード数
 max_nodes_bridge: 5 # ブリッジあたりの最大ノード数
 kink_max_nodes: 3 # ねじれ最適化の最大ノード数
 max_seq_kink: 2 # 連続ねじれの上限
 refine_mode: null # 精密化戦略: peak, minima, null (自動)
```

---

## TS 最適化セクション

### `hessian_dimer`

ヘシアン・ダイマー TS 最適化（`tsopt --opt-mode grad`）。

```yaml
hessian_dimer:
 thresh_loose: gau_loose # 緩い収束プリセット
 thresh: baker # メイン収束プリセット
 update_interval_hessian: 500 # ヘシアン再構築間隔
 neg_freq_thresh_cm: 5.0 # 負振動数閾値 (cm^-1)
 flatten_amp_ang: 0.1 # フラット化振幅 (Å)
 flatten_max_iter: 50 # フラット化反復上限（デフォルト 50、--no-flatten で 0 に設定）
 flatten_sep_cutoff: 0.0 # 代表原子間の最小距離
 flatten_k: 10 # モードあたりのサンプル代表原子数
 flatten_loop_bofill: false # フラット化変位に Bofill 更新
 mem: 100000 # ソルバーのメモリ上限
 device: auto # 固有値ソルバーのデバイス選択
 root: 0 # ターゲット TS ルートインデックス
 partial_hessian_flatten: true # 部分ヘシアンを虚モード検出に使用
 ml_only_hessian_dimer: false # ダイマー方向決定に ML 領域のみのヘシアンを使用
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
 rotation_remove_trans: true # 並進成分の除去
 trans_force_f_perp: true # 並進に垂直な力の投影
 bonds: null # 拘束用の結合リスト
 N_hessian: null # ヘシアンサイズの上書き
 bias_rotation: false # 回転探索のバイアス
 bias_translation: false # 並進探索のバイアス
 bias_gaussian_dot: 0.1 # ガウスバイアスの内積
 seed: null # 回転の乱数シード
 write_orientations: true # 回転方向の書き出し
 forward_hessian: true # ヘシアンの前方伝搬
 lbfgs:
 # lbfgs セクションと同じキー
 thresh: baker
 max_cycles: 10000
```

**注記:**
- `flatten_max_iter` は虚振動数モードフラットニングの最大反復回数を制御します。デフォルト値は 50 です。
- CLI フラグ `--flatten` / `--no-flatten`（`tsopt` および `all`）はこの設定と連動します。`--flatten` はデフォルトの `flatten_max_iter`（50）でフラットニングループを有効化し、`--no-flatten` は `flatten_max_iter` を 0 に強制してループを無効化します。`--flatten` と同時に YAML で `flatten_max_iter` を明示指定した場合は、YAML の値が優先されます。

---

### `rsirfo`

RS-I-RFO TS 最適化（`tsopt --opt-mode hess`）。

```yaml
rsirfo:
 thresh: baker # RS-I-RFO 収束プリセット
 max_cycles: 10000 # 反復上限
 print_every: 100 # ログ出力間隔
 min_step_norm: 1.0e-08 # 最小ステップノルム
 assert_min_step: true # ステップ停滞時にアサート
 roots: [0] # ターゲットルートインデックス（pysisyphus デフォルト; mlmm では未設定）
 hessian_ref: null # 参照ヘシアン
 rx_modes: null # 反応モード定義
 prim_coord: null # 監視する主座標
 rx_coords: null # 監視する反応座標
 hessian_update: bofill # ヘシアン更新スキーム
 hessian_recalc_reset: true # 正確なヘシアン後に再計算カウンタをリセット
 hessian_init: calc # ヘシアン初期化
 hessian_recalc: 500 # ヘシアン再構築間隔
 max_micro_cycles: 50 # マクロサイクルあたりのマイクロイテレーション数
 augment_bonds: false # 結合解析に基づく反応経路の拡張
 min_line_search: false # 虚モードに沿ったラインサーチ（pysisyphus デフォルト）
 max_line_search: false # 最小化部分空間でのラインサーチ（pysisyphus デフォルト）
 assert_neg_eigval: false # 収束時に負の固有値を要求
 trust_radius: 0.10 # 信頼領域半径
 trust_update: true # 信頼領域更新
 trust_min: 0.0001 # 最小信頼半径
 trust_max: 0.10 # 最大信頼半径（v0.2.8 で ML/MM 安定性のため厳格化）
 small_eigval_thresh: 1.0e-08 # 安定性のための固有値閾値
 out_dir: ./result_tsopt/ # 出力ディレクトリ
```

---

### `stopt`

ストリング最適化（GS/DMF）の設定。path-opt と path-search で使用。

```yaml
stopt:
 type: string           # 最適化タイプラベル（StringOptimizer用）
 thresh: gau_loose      # ストリング最適化の収束プリセット
 stop_in_when_full: 300 # ストリングが満杯時の早期停止閾値
 align: false           # アライメントトグル
 scale_step: global     # ステップスケーリングモード
 max_cycles: 300        # ストリング最適化の最大反復数
 dump: false            # 軌跡/リスタートデータ出力
 dump_restart: false    # リスタートチェックポイントの出力
 reparam_thresh: 0.0    # 再パラメータ化閾値
 coord_diff_thresh: 0.0 # 座標差分閾値
 out_dir: ./result_path_opt/  # 出力ディレクトリ
 print_every: 10        # ログ出力間隔
 lbfgs:
   # 単一構造最適化用（HEI±1、ねじれノード）
   thresh: gau
   max_cycles: 10000
   # ...（詳細は lbfgs セクション参照）
 rfo:
   # 単一構造最適化用
   thresh: gau
   max_cycles: 10000
   # ...（詳細は rfo セクション参照）
```

**注意:**
- `stopt.lbfgs` / `stopt.rfo` は HEI±1 端点最適化およびねじれノード最適化に使用される単一構造最適化の設定
- 外側の `stopt` キーはストリング最適化（GS または DMF ラッパー）を制御

---

## IRC セクション

(ja-irc-section)=
### `irc` (section)

IRC 積分設定。

```yaml
irc:
 step_length: 0.1 # 積分ステップ長
 max_cycles: 125 # IRC の最大ステップ数
 downhill: false # 下り方向のみ追跡
 forward: true # 順方向に伝搬
 backward: true # 逆方向に伝搬
 root: 0 # 基準振動モードのルートインデックス
 hessian_init: calc # ヘシアン初期化ソース
 hessian_update: bofill # ヘシアン更新スキーム
 hessian_recalc: null # ヘシアン再構築間隔
 displ: energy # 変位構築方法
 displ_energy: 0.001 # エネルギーベースの変位スケーリング
 displ_length: 0.1 # 長さベースの変位フォールバック
 rms_grad_thresh: 0.001 # RMS 勾配の収束閾値
 hard_rms_grad_thresh: null # ハード RMS 勾配停止閾値
 energy_thresh: 0.000001 # エネルギー変化閾値
 imag_below: 0.0 # 虚振動数カットオフ
 force_inflection: true # 変曲点検出の強制
 check_bonds: false # 伝搬中の結合チェック
 out_dir: ./result_irc/ # 出力ディレクトリ
 prefix: "" # ファイル名プレフィックス
 dump_fn: irc_data.h5 # IRC データファイル名
 dump_every: 5 # ダンプ間隔
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
 amplitude_ang: 0.8 # モード変位振幅 (Å)
 n_frames: 20 # モードアニメーションのフレーム数
 max_write: 10 # 書き出すモードの最大数
 sort: value # ソート順: "value" または "abs"
 out_dir: ./result_freq/ # 出力ディレクトリ
```

---

### `thermo`

熱化学設定。

```yaml
thermo:
 temperature: 298.15 # 熱化学温度 (K)
 pressure_atm: 1.0 # 熱化学圧力 (atm)
 dump: false # thermoanalysis.yaml の書き出し
```

---

### `microiter`

ML/MM最適化用のマイクロイテレーション設定。`--microiter` 有効時、MLリージョンの
マクロステップ間でMMリージョンをL-BFGSで緩和（ML原子は凍結）します。

```yaml
microiter:
 micro_thresh: null       # MM緩和の収束プリセット（L-BFGS）; null → マクロステップと同じ
 micro_max_cycles: 10000  # マイクロイテレーションあたりの最大L-BFGS反復数
```

**注意:**
- CLIフラグ `--microiter` / `--no-microiter` で有効化（デフォルト: 有効）
- `opt`（`--opt-mode hess`）および `tsopt`（`--opt-mode hess`）で使用可能
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
 max_cycle: 100 # 最大 SCF 反復数
 grid_level: 3 # PySCF グリッドレベル
 lowmem: true # closed-shell GPU で gpu4pyscf rks_lowmem.RKS を使用
 verbose: 4 # PySCF 出力詳細レベル
 out_dir: ./result_dft/ # 出力ディレクトリ
```

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

calc:
 model_charge: 0
 model_mult: 1
 backend: uma                  # ML バックエンド: uma | orb | mace | aimnet2
 embedcharge: false            # xTB 点電荷埋め込み補正
 uma_model: uma-s-1p1          # uma-s-1p1 | uma-m-1p1
 ml_device: auto
 hessian_calc_mode: Analytical   # VRAM に余裕がある場合に推奨
 mm_device: cpu
 mm_fd: true
 use_bfactor_layers: true # 入力 PDB の B-factor から層を読み取り

gs:
 max_nodes: 12
 climb: true
 climb_lanczos: true

opt:
 thresh: gau
 max_cycles: 300
 dump: false
 out_dir: ./result_all/

stopt:
 thresh: gau_loose
 max_cycles: 300
 lbfgs:
   thresh: gau
   max_cycles: 10000
 rfo:
   thresh: gau
   max_cycles: 10000

bond:
 bond_factor: 1.2
 delta_fraction: 0.05

search:
 max_depth: 10
 max_nodes_segment: 10

freq:
 max_write: 10
 amplitude_ang: 0.8

thermo:
 temperature: 298.15
 pressure_atm: 1.0

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
- [ML/MM 計算機](mlmm-calc.md) - ML/MM 計算機の詳細
