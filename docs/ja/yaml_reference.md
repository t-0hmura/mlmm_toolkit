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
| [`freq`](#ja-freq-section) | 振動数解析設定 | freq |
| [`thermo`](#thermo) | 熱化学設定 | freq |
| [`dft`](#ja-dft-section) | DFT計算設定 | dft |
| [`bias`](#bias) | 調和バイアス設定 | scan, scan2d, scan3d |
| [`bond`](#bond) | 結合変化検出設定 | scan, path-search |
| [`search`](#search) | 再帰的経路探索設定 | path-search |
| [`hessian_dimer`](#hessian_dimer) | ヘシアン・ダイマーTS 最適化 | tsopt |
| [`rsirfo`](#rsirfo) | RS-I-RFO TS 最適化 | tsopt |
| [`sopt`](#sopt) | path-search用単一構造最適化 | path-search |

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

ML/MM 計算機（UMA + hessian_ff）の設定。

```yaml
calc:
 input_pdb: null # 入力 PDB ファイルパス (CLI --input から設定)
 real_parm7: null # 全系の Amber parm7 トポロジ (CLI --real-parm7)
 model_pdb: null # ML 領域を定義する PDB (CLI --model-pdb)
 model_charge: 0 # ML 領域の電荷 (CLI -q で上書き)
 model_mult: 1 # ML 領域のスピン多重度 (CLI -m で上書き)
 link_mlmm: null # リンク原子ペアの明示指定 (null で自動検出)
 uma_model: uma-s-1p1 # UMA 事前学習モデル名
 uma_task_name: omol # UMA バッチに記録されるタスクタグ
 ml_hessian_mode: Analytical # ML ヘシアンモード: "Analytical" または "FiniteDifference"
 hessian_calc_mode: null # ml_hessian_mode のエイリアス
 out_hess_torch: true # ヘシアンを torch.Tensor で返す
 H_double: false # ヘシアンを float64 で組み立て・返却
 ml_device: auto # ML デバイス: "cuda", "cpu", "auto"
 ml_cuda_idx: 0 # CUDA デバイスインデックス
 mm_backend: hessian_ff # MM バックエンド: "hessian_ff" (解析的) | "openmm" (FD ヘシアン)
 mm_device: cpu # MM デバイス (hessian_ff は CPU のみ、OpenMM は CUDA/CPU 対応)
 mm_cuda_idx: 0 # MM CUDA インデックス (OpenMM のみ)
 mm_threads: 16 # MM 計算のスレッド数
 mm_fd: true # MM ヘシアンを計算するか
 mm_fd_dir: null # MM ヘシアンログの出力ディレクトリ
 mm_fd_delta: 0.001 # 有限差分ステップ（保持）
 symmetrize_hessian: true # 最終ヘシアンを 0.5*(H+H^T) で対称化
 print_timing: true # ML/MM ヘシアンのタイミング内訳を表示
 print_vram: true # CUDA VRAM 使用量を表示
 return_partial_hessian: false # 計算機の基底既定値（CLI ラッパー側で true 既定を適用する場合あり）
 freeze_atoms: [] # geom.freeze_atoms から継承
 # 層設定:
 hess_cutoff: null # Å: Hessian 対象 MM の距離カットオフ
 movable_cutoff: null # Å: movable MM の距離カットオフ
 use_bfactor_layers: false # 入力 PDB の B-factor から層を読み取り
 hess_mm_atoms: null # 明示的 Hessian 対象 MM 原子インデックス (0始まり)
 movable_mm_atoms: null # 明示的 movable MM 原子インデックス (0始まり)
 frozen_mm_atoms: null # 明示的 frozen MM 原子インデックス (0始まり)
```

**注記:**
- `ml_hessian_mode: Analytical` が推奨です（VRAM に余裕がある場合）。
- `hess_cutoff`/`movable_cutoff` を指定しない場合、ML 以外の全原子が Hessian 対象 MM に分類されます。
- `use_bfactor_layers: true` を設定すると、`define-layer` で書き込んだ B-factor から層割り当てを読み取ります。
- 明示的インデックス（`hess_mm_atoms` 等）が設定された場合、カットオフや B-factor よりも優先されます。
- `opt`/`tsopt`/`irc`/`freq` は、YAML で `mlmm.return_partial_hessian` を明示しない場合に部分ヘシアンを既定で使用します。
- これらのコマンドで完全ヘシアンを強制するには `mlmm.return_partial_hessian: false` を明示してください。

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
 line_search: true # ラインサーチを有効化
 dump: false # 軌跡/リスタートデータの出力
 dump_restart: false # リスタートチェックポイントの出力
 prefix: "" # ファイル名プレフィックス
 out_dir:./result_opt/ # 出力ディレクトリ
```

**収束プリセット:**

| プリセット | Max Force | RMS Force | Max Step | RMS Step |
|-----------|-----------|-----------|----------|----------|
| `gau_loose` | 2.5e-3 | 1.7e-3 | 1.0e-2 | 6.7e-3 |
| `gau` | 4.5e-4 | 3.0e-4 | 1.8e-3 | 1.2e-3 |
| `gau_tight` | 1.5e-5 | 1.0e-5 | 6.0e-5 | 4.0e-5 |
| `gau_vtight` | 2.0e-6 | 1.0e-6 | 6.0e-6 | 4.0e-6 |
| `baker` | 3.0e-4 | 2.0e-4 | 3.0e-4 | 2.0e-4 |

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
 trust_radius: 0.1 # 信頼領域半径
 trust_update: true # 信頼領域更新を有効化
 trust_min: 0.0 # 最小信頼半径
 trust_max: 0.1 # 最大信頼半径
 max_energy_incr: null # ステップあたりの許容エネルギー増加
 hessian_update: bfgs # ヘシアン更新スキーム: bfgs, bofill 等
 hessian_init: calc # ヘシアン初期化: calc, unit 等
 hessian_recalc: 200 # N ステップごとにヘシアンを再構築
 hessian_recalc_adapt: null # 適応的ヘシアン再構築係数
 small_eigval_thresh: 1.0e-08 # 安定性のための固有値閾値
```

---

## 経路最適化セクション

### `gs`

Growing String Method（GSM）の設定。

```yaml
gs:
 fix_first: true # 最初の端点を固定
 fix_last: true # 最後の端点を固定
 max_nodes: 10 # 最大ストリングノード数
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
 max_cycles: 300 # DMF/IPOPT の最大反復数
 correlated: true # 相関 DMF 伝搬
 sequential: true # 逐次 DMF 実行
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
 kink_max_nodes: 3 # キンク最適化の最大ノード数
 max_seq_kink: 2 # 連続キンクの上限
 refine_mode: null # 精密化戦略: peak, minima, null (自動)
```

---

## TS 最適化セクション

### `hessian_dimer`

ヘシアン・ダイマー TS 最適化（tsopt --opt-mode light）。

```yaml
hessian_dimer:
 thresh_loose: gau_loose # 緩い収束プリセット
 thresh: baker # メイン収束プリセット
 update_interval_hessian: 500 # ヘシアン再構築間隔
 neg_freq_thresh_cm: 5.0 # 負振動数閾値 (cm^-1)
 flatten_amp_ang: 0.1 # フラット化振幅 (Å)
 flatten_max_iter: 50 # フラット化反復上限
 flatten_sep_cutoff: 0.0 # 代表原子間の最小距離
 flatten_k: 10 # モードあたりのサンプル代表原子数
 flatten_loop_bofill: false # フラット化変位に Bofill 更新
 mem: 100000 # ソルバーのメモリ上限
 device: auto # 固有値ソルバーのデバイス選択
 root: 0 # ターゲット TS ルートインデックス
 partial_hessian_flatten: true # 部分ヘシアンを虚モード検出に使用
```

---

### `rsirfo`

RS-I-RFO TS 最適化（tsopt --opt-mode heavy）。

```yaml
rsirfo:
 thresh: baker # RS-IRFO 収束プリセット
 max_cycles: 10000 # 反復上限
 hessian_update: bofill # ヘシアン更新スキーム
 hessian_init: calc # ヘシアン初期化
 hessian_recalc: 200 # ヘシアン再構築間隔
 trust_radius: 0.10 # 信頼領域半径
 trust_update: true # 信頼領域更新
 trust_min: 0.00 # 最小信頼半径
 trust_max: 0.30 # 最大信頼半径
out_dir:./result_tsopt/ # 出力ディレクトリ
```

---

### `sopt`

path-search で使用する単一構造最適化設定（HEI±1 ノードおよび kink ノード）。

```yaml
sopt:
 lbfgs:
 # 上記 lbfgs セクションと同じキー
 thresh: gau
 max_cycles: 10000
 out_dir:./result_path_search/
 dump: false
 # ...（詳細は lbfgs セクション参照）
 rfo:
 # 上記 rfo セクションと同じキー
 thresh: gau
 max_cycles: 10000
 out_dir:./result_path_search/
 dump: false
 # ...（詳細は rfo セクション参照）
```

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
 root: 0 # 基準振動モードのルートインデックス
 hessian_init: calc # ヘシアン初期化ソース
 hessian_update: bofill # ヘシアン更新スキーム
 hessian_recalc: null # ヘシアン再構築間隔
 rms_grad_thresh: 0.001 # RMS 勾配の収束閾値
 energy_thresh: 0.000001 # エネルギー変化閾値
 out_dir:./result_irc/ # 出力ディレクトリ
```

---

## 振動解析セクション

(ja-freq-section)=
### `freq` (section)

振動数解析設定。

```yaml
freq:
 amplitude_ang: 0.8 # モード変位振幅 (Å)
 n_frames: 20 # モードアニメーションのフレーム数
 max_write: 10 # 書き出すモードの最大数
 sort: value # ソート順: "value" または "abs"
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

## DFT セクション

(ja-dft-section)=
### `dft` (section)

DFT 計算設定。

```yaml
dft:
 func_basis: B3LYP/6-31G* # 汎関数/基底関数の組み合わせ文字列
 max_cycle: 100 # 最大 SCF 反復数
 conv_tol: 1.0e-09 # SCF 収束許容値 (Hartree)
 grid_level: 3 # PySCF グリッドレベル
```

---

## スキャン関連セクション

### `bias`

調和バイアス設定。

```yaml
bias:
 k: 100.0 # 調和バイアス強度 (eV/Å^2)
```

---

### `bond`

UMA ベースの結合変化検出。

```yaml
bond:
 device: cuda # UMA デバイス
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
 uma_model: uma-s-1p1
 ml_device: auto
 ml_hessian_mode: Analytical # VRAM に余裕がある場合に推奨
 hess_cutoff: 3.6 # Layer 2 の距離カットオフ (Å)
 movable_cutoff: 8.0 # Layer 3 の距離カットオフ (Å)

gs:
 max_nodes: 12
 climb: true
 climb_lanczos: true

opt:
 thresh: gau
 max_cycles: 300
 dump: false
 out_dir:./result_all/

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
 func_basis: B3LYP/6-31G*
 grid_level: 3
```

---

## 参照

- [all](all.md) - メインワークフロー
- [opt](opt.md) - 単一構造最適化
- [tsopt](tsopt.md) - 遷移状態最適化
- [path-search](path_search.md) - 再帰的 MEP 探索
- [freq](freq.md) - 振動解析
- [dft](dft.md) - DFT 計算
- [mlmm_calc](mlmm_calc.md) - ML/MM 計算機の詳細
