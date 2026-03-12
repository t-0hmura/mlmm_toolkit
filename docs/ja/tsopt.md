# `tsopt`

## 概要

> **概要:** Dimer（`--opt-mode grad`）または RS-I-RFO（`--opt-mode hess`、デフォルト）を使用して遷移状態*候補*を最適化します。マイクロイテレーション（`--microiter`、デフォルト有効）は `hess` モードで ML 1 ステップ RS-I-RFO と MM 緩和を交互に実行します。検証済み TS は**正確に 1 つ**の虚数振動数を示すべきです。必ず freq/IRC でモード/結合性を確認してください。

### `--opt-mode` の選択
- **`--opt-mode hess`（RS-I-RFO）** を使用: デフォルトの保守的なオプティマイザーで、ヘシアン計算のコストを許容できる場合。`--microiter`（デフォルト有効）により ML と MM 領域を交互に最適化。
- **`--opt-mode grad`（Dimer）** を使用: 軽量な探索が必要な場合、または複数の TS 推測構造から素早く反復する場合。`--ml-only-hessian-dimer` で ML 領域のみのヘシアンを Dimer 方向決定に使用（高速だが精度は低下）。

`mlmm tsopt` は ML/MM 計算機に特化した遷移状態最適化を実行します。`-b/--backend` で ML バックエンドを選択可能です（`uma`、`orb`、`mace`、`aimnet2`）。`--embedcharge` で xTB 点電荷埋め込み補正を有効化し、MM 環境から ML 領域への静電的影響を考慮できます。オプティマイザーは TS 推測構造から開始し、一次鞍点へ精密化します。

### 主な特徴
- **部分ヘシアンガイド付き Dimer:** ゆるい/最終 Dimer ループ中、hessian_ff 有限差分ヘシアンは無効化（`mm_fd=False`）されます。MLIP ヘシアンは MM 原子をゼロパディングして完全な 3N x 3N 空間に埋め込まれ、Dimer の方向更新をガイドする部分ヘシアンを提供します。
- **完全ヘシアンによるフラットニングループ:** 探索がフラットニングループに入ると、完全な ML/MM ヘシアン（MM 有限差分ブロックを含む）が正確に 1 回計算され、その後 Dimer セグメント間で Bofill ステップによりアクティブ部分空間で更新されます。
- **PHVA + TR 射影:** アクティブ自由度射影と質量加重並進/回転除去は `freq.py` をミラーリングし、一貫した虚振動数モード解析とモード書き出しを保証します。
- **出力変換:** `--convert-files`（デフォルト）により、PDB 入力は `.pdb` にミラーリングされ（`--dump` 時）、虚振動数モードは `_trj.xyz` とともに `.pdb` としてもエクスポートされます。

## 最小例

```bash
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_tsopt
```

## 出力の見方

- `result_tsopt/final_geometry.pdb`（または `final_geometry.xyz`）
- `result_tsopt/vib/final_imag_mode_*_trj.xyz`
- `result_tsopt/vib/final_imag_mode_*.pdb`

## よくある例

1. VRAM に余裕がある場合に light モード + 解析的ヘシアンで実行する。

```bash
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --opt-mode light --hessian-calc-mode Analytical --out-dir ./result_tsopt_light
```

2. 最適化軌跡を保存して確認する。

```bash
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --dump --out-dir ./result_tsopt_dump
```

3. heavy モードを YAML 上書きと併用する。

```bash
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --opt-mode heavy --config tsopt.yaml --out-dir ./result_tsopt_heavy
```

4. MACE バックエンドで TS 最適化を実行する。

```bash
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --backend mace --out-dir ./result_tsopt_mace
```

## ワークフロー

1. **入力処理** -- 酵素 PDB、Amber トポロジー、ML 領域定義を読み込みます。電荷/スピンを解決します。CLI と YAML の凍結原子がマージされます。
2. **ML/MM 計算機の構築** -- ML/MM 計算機（MLIP バックエンド + hessian_ff）を構築します。`-b/--backend` で ML バックエンドを選択し（デフォルト: `uma`）、`--hessian-calc-mode` は MLIP がヘシアンを解析的に評価するか有限差分で評価するかを制御します。`--embedcharge` で xTB 点電荷埋め込み補正を有効化できます。
3. **Light モード（Dimer）:**
   - ヘシアン Dimer ステージは、部分ヘシアン（アクティブ部分空間、TR 射影済み）を評価して Dimer 方向を定期的に更新します。
   - フラットニングループが有効な場合（`--flatten`）、保存されたアクティブヘシアンは変位と勾配差分を使用した Bofill 更新により更新されます。各ループで虚振動数モードを推定し、1 回フラットニングし、Dimer 方向を更新し、Dimer + LBFGS マイクロセグメントを実行します。
4. **Heavy モード（RS-I-RFO）:**
   - RS-I-RFO オプティマイザーを、`rsirfo` YAML セクションで定義されたオプションのヘシアン参照ファイルとマイクロサイクル制御とともに実行します。
   - `--flatten` が有効で収束後に 2 つ以上の虚振動数モードが残る場合、余分なモードをフラットニングし、1 つだけ残るかフラットニング反復上限に達するまで RS-I-RFO を再実行します。
5. **モードエクスポートと変換** -- 収束した虚振動数モードは常に `vib/final_imag_mode_*_trj.xyz` に書き出され、入力が PDB で変換が有効な場合は `.pdb` にもミラーリングされます。最適化軌跡と最終ジオメトリも `--dump` 時に入力テンプレート経由で PDB に変換されます。

## CLIオプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 開始ジオメトリ（PDB または XYZ）。XYZ の場合はトポロジーに `--ref-pdb` を使用。 | 必須 |
| `--ref-pdb FILE` | 入力が XYZ の場合の参照 PDB トポロジー。 | _None_ |
| `--parm PATH` | 全酵素の Amber parm7 トポロジー。 | 必須 |
| `--model-pdb PATH` | ML 領域原子を含む PDB。`--detect-layer` 有効時はオプション。 | _None_ |
| `--model-indices TEXT` | ML 領域のカンマ区切り原子インデックス（範囲指定可）。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` を 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer / --no-detect-layer` | 入力 PDB の B 因子から ML/MM レイヤーを検出。 | `True` |
| `-q, --charge INT` | ML 領域の総電荷。 | 必須 |
| `-m, --multiplicity INT` | ML 領域のスピン多重度 (2S+1)。 | `1` |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切りインデックス（YAML `geom.freeze_atoms` とマージ）。 | _None_ |
| `--hess-cutoff FLOAT` | ML 領域からの Hessian-MM 原子の距離カットオフ (Å)。可動 MM 原子に適用。`0.0` は ML のみの部分ヘシアン。エイリアス: `--radius-hessian`。 | `0.0` |
| `--movable-cutoff FLOAT` | 可動 MM 原子の距離カットオフ (Å)。 | _None_ |
| `--hessian-calc-mode CHOICE` | MLIP ヘシアンモード: `Analytical` または `FiniteDifference`。 | `FiniteDifference` |
| `--max-cycles INT` | 最大総オプティマイザーサイクル。 | `10000` |
| `--opt-mode CHOICE` | TS オプティマイザーモード: `grad`（Dimer）または `hess`（RS-I-RFO）。エイリアス `light`/`heavy` も使用可。 | `hess` |
| `--microiter/--no-microiter` | マイクロイテレーション: ML 1 ステップ（RS-I-RFO）+ MM 緩和（LBFGS）を交互に実行。`hess` モードでのみ有効。 | `True` |
| `--ml-only-hessian-dimer/--no-ml-only-hessian-dimer` | `grad` モードで Dimer 方向決定に ML 領域のみのヘシアンを使用。高速だが精度は低下。 | `False` |
| `--flatten/--no-flatten` | 余分な虚振動数モードフラットニングループの有効化/無効化。`--flatten` はデフォルト反復回数（50）を使用、`--no-flatten` は 0 に強制。light と heavy の両モードに適用。 | _None_（YAML/デフォルト依存; 実質的に 50 回で有効） |
| `--dump/--no-dump` | 連結軌跡 `optimization_all_trj.xyz` を書き出し。 | `False` |
| `--convert-files/--no-convert-files` | PDB 入力時の XYZ/TRJ から PDB コンパニオンの切り替え。 | `True` |
| `-o, --out-dir TEXT` | 出力ディレクトリ。 | `./result_tsopt/` |
| `--thresh TEXT` | 収束プリセット（`gau_loose\|gau\|gau_tight\|gau_vtight\|baker\|never`）。 | _None_ |
| `--partial-hessian-flatten / --full-hessian-flatten` | フラットニングループでの虚振動数モード検出に部分ヘシアン（ML のみ）を使用。 | `True`（部分） |
| `--active-dof-mode CHOICE` | 最終振動解析のアクティブ自由度: `all`、`ml-only`、`partial`、`unfrozen`。 | `partial` |
| `--config FILE` | 明示 CLI オプションより前に適用するベース YAML 設定ファイル。 | _None_ |
| `--show-config/--no-show-config` | 解決後の設定レイヤーを表示して実行を継続。 | `False` |
| `-b, --backend CHOICE` | ML 領域の MLIP バックエンド: `uma`（デフォルト）、`orb`、`mace`、`aimnet2`。 | `uma` |
| `--embedcharge/--no-embedcharge` | xTB 点電荷埋め込み補正の有効化。MM 環境から ML 領域への静電的影響を考慮。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに入力/設定を検証し、実行計画を表示。`--help-advanced` に表示。 | `False` |

## 出力

```
out_dir/ (デフォルト: ./result_tsopt/)
├── final_geometry.xyz             # 常に書き出し
├── final_geometry.pdb             # 入力が PDB の場合
├── optimization_all_trj.xyz       # 連結 Dimer セグメント（--dump 時）
├── optimization_all.pdb           # PDB コンパニオン（--dump かつ入力が PDB の場合）
├── vib/
│   ├── final_imag_mode_±XXXX.Xcm-1_trj.xyz  # 虚振動数モード軌跡
│   └── final_imag_mode_±XXXX.Xcm-1.pdb      # 虚振動数モード PDB コンパニオン
└── .dimer_mode.dat                # Dimer 方向シード（light モード）
```

## YAML設定

設定は **デフォルト < config < 明示CLI < override** の順で適用されます。
共有セクションは [YAML リファレンス](yaml_reference.md) を再利用します。ワークフローに合致している場合は以下のブロック全体をそのまま保持し、変更が必要な値のみ調整してください。

```yaml
geom:
 coord_type: cart                  # 座標タイプ: デカルト vs dlc 内部座標
 freeze_atoms: []                  # 0 始まり凍結原子（CLI/リンク検出とマージ）
calc:
 charge: 0                         # 総電荷（CLI 上書き）
 spin: 1                           # スピン多重度 2S+1
mlmm:
 real_parm7: real.parm7            # Amber parm7 トポロジー
 model_pdb: ml_region.pdb          # ML 領域定義
 backend: uma                      # ML バックエンド (uma/orb/mace/aimnet2)
 embedcharge: false                # xTB 点電荷埋め込み補正
 uma_model: uma-s-1p1              # uma-s-1p1 | uma-s-1p1 | uma-m-1p1
 uma_task_name: omol                # UMA タスク名 (backend=uma 時)
 ml_device: auto                   # ML デバイス選択
 ml_hessian_mode: Analytical        # ヘシアンモード選択
opt:
 thresh: baker                     # 収束プリセット（Gaussian/Baker 式）
 max_cycles: 10000                 # オプティマイザーサイクル上限
 print_every: 100                  # ログ出力間隔
 min_step_norm: 1.0e-08            # ステップ受け入れの最小ノルム
 assert_min_step: true             # ステップが閾値以下で停止
 rms_force: null                   # 明示的 RMS 力目標
 rms_force_only: false             # RMS 力収束のみに依存
 max_force_only: false             # 最大力収束のみに依存
 force_only: false                 # 変位チェックをスキップ
 converge_to_geom_rms_thresh: 0.05  # 参照への収束時の geom RMS 閾値
 overachieve_factor: 0.0           # 閾値を厳しくする係数
 check_eigval_structure: false     # ヘシアン固有値構造の検証
 line_search: true                 # ラインサーチを有効化
 dump: false                       # 軌跡/リスタートデータのダンプ
 dump_restart: false               # リスタートチェックポイントのダンプ
 prefix: ""                        # ファイル名プレフィックス
 out_dir: ./result_tsopt/          # 出力ディレクトリ
hessian_dimer:
 thresh_loose: gau_loose           # ゆるい収束プリセット
 thresh: baker                     # メイン収束プリセット
 update_interval_hessian: 500      # ヘシアン再構築間隔
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
  rotation_method: fourier         # 回転オプティマイザー手法
  rotation_thresh: 0.0001          # 回転収束閾値
  rotation_tol: 1                  # 回転許容係数
  rotation_max_element: 0.001      # 回転行列の最大要素
  rotation_interpolate: true       # 回転ステップの補間
  rotation_disable: false          # 回転を完全に無効化
  rotation_disable_pos_curv: true  # 正の曲率検出時に無効化
  rotation_remove_trans: true      # 並進成分を除去
  trans_force_f_perp: true         # 並進に垂直な力を射影
  bonds: null                      # 拘束用の結合リスト
  N_hessian: null                  # ヘシアンサイズの上書き
  bias_rotation: false             # 回転探索にバイアス
  bias_translation: false          # 並進探索にバイアス
  bias_gaussian_dot: 0.1           # ガウスバイアスの内積
  seed: null                       # 回転用の RNG シード
  write_orientations: true         # 回転方向を書き出し
  forward_hessian: true            # ヘシアンを前方伝播
 lbfgs:
  thresh: baker                    # LBFGS 収束プリセット
  max_cycles: 10000                # 反復上限
  print_every: 100                 # ログ出力間隔
  min_step_norm: 1.0e-08           # 受け入れ最小ステップノルム
  assert_min_step: true            # ステップ停滞時にアサート
  max_step: 0.3                    # 最大ステップ長
  control_step: true               # 適応的ステップ長制御
  double_damp: true                # ダブルダンピングセーフガード
  keep_last: 7                     # LBFGS バッファの履歴サイズ
  beta: 1.0                        # 初期ダンピングベータ
  mu_reg: null                     # 正則化強度
  max_mu_reg_adaptions: 10         # mu 適応の上限
rsirfo:
 thresh: baker                     # RS-IRFO 収束プリセット
 max_cycles: 10000                 # 反復上限
 print_every: 100                  # ログ出力間隔
 min_step_norm: 1.0e-08            # 受け入れ最小ステップノルム
 assert_min_step: true             # ステップ停滞時にアサート
 roots: [0]                        # 対象ルートインデックス
 hessian_ref: null                 # 参照ヘシアン
 hessian_update: bofill            # ヘシアン更新方式の上書き
 hessian_recalc_reset: true        # 正確なヘシアン後にリカルクカウンターをリセット
 max_micro_cycles: 50              # マクロサイクルごとのマイクロイテレーション
 augment_bonds: false              # 結合解析に基づく反応経路の拡張
 min_line_search: false            # 最小ラインサーチステップを強制
 max_line_search: false            # 最大ラインサーチステップを強制
 assert_neg_eigval: false          # 収束時に負の固有値を要求
```

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細なトラブルシューティングガイド

- [opt](opt.md) -- 単一構造の構造最適化
- [freq](freq.md) -- 検証済み TS の単一虚数振動数を確認
- [irc](irc.md) -- 最適化された TS からの反応経路追跡
- [all](all.md) -- 抽出 -> MEP -> tsopt -> IRC -> freq を連鎖させるエンドツーエンドワークフロー
- [YAML リファレンス](yaml_reference.md) -- `hessian_dimer`（ヘシアンガイド付き Dimer）と `rsirfo` の完全な設定オプション
- [用語集](glossary.md) -- TS、Dimer、RS-I-RFO、ヘシアンの定義
