# `tsopt`

## 概要

> **要約:** ML/MM 計算機（FAIR-Chem UMA + OpenMM）を使用した、部分ヘシアンガイド付き Dimer（`--opt-mode light`）または RS-I-RFO（`--opt-mode heavy`、デフォルト）による遷移状態最適化。

`mlmm tsopt` は ML/MM 計算機に特化した遷移状態最適化を実行します。オプティマイザーは TS 推測構造（例: `path-opt`/`path-search` からの最高エネルギーイメージ、またはユーザー独自の構造）から開始し、一次鞍点へ精密化します。

### 主な特徴

- **部分ヘシアンガイド付き Dimer:** ゆるい/最終 Dimer ループ中、OpenMM 有限差分ヘシアンは無効化（`mm_fd=False`）されます。UMA ヘシアンは MM 原子をゼロパディングして完全な 3N x 3N 空間に埋め込まれ、Dimer の方向更新をガイドする部分ヘシアンを提供します。
- **完全ヘシアンによるフラットンループ:** 探索がフラットンループに入ると、完全な ML/MM ヘシアン（MM 有限差分ブロックを含む）が正確に 1 回計算され、その後 Dimer セグメント間で Bofill ステップによりアクティブ部分空間で更新されます。
- **PHVA + TR 射影:** アクティブ自由度射影と質量加重並進/回転除去は `freq.py` をミラーリングし、一貫した虚数モード解析とモード書き出しを保証します。

### `--opt-mode` の選択

- **`heavy`（RS-I-RFO、デフォルト）:** 完全ヘシアンを使用する保守的なオプティマイザー。一般的により堅牢。
- **`light`（Dimer）:** 部分ヘシアンガイド付き Dimer を使用する軽量探索。ステップあたりのコストが低い場合が多い。

## 最小例

```bash
mlmm tsopt -i ts_guess.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --out-dir ./result_tsopt
```

## 出力の見方

- `result_tsopt/summary.md`
- `result_tsopt/key_ts.xyz`（または `key_ts.pdb`）
- `result_tsopt/key_imag_mode.trj`
- `result_tsopt/final_geometry.pdb`（または `final_geometry.xyz`）
- `result_tsopt/vib/final_imag_mode_*.trj`
- `result_tsopt/vib/final_imag_mode_*.pdb`

## よくある例

1. VRAM に余裕がある場合に light モード + 解析的ヘシアンで実行する。

```bash
mlmm tsopt -i ts_guess.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --opt-mode light --hessian-calc-mode Analytical --out-dir ./result_tsopt_light
```

2. 最適化軌跡を保存して確認する。

```bash
mlmm tsopt -i ts_guess.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --dump --out-dir ./result_tsopt_dump
```

3. heavy モードを YAML 上書きと併用する。

```bash
mlmm tsopt -i ts_guess.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --opt-mode heavy --config tsopt.yaml --out-dir ./result_tsopt_heavy
```

## 使用法

```bash
mlmm tsopt -i INPUT.pdb --real-parm7 real.parm7 --model-pdb model.pdb \
    -q CHARGE [-m SPIN] [--freeze-atoms "1,3,5"] [--max-cycles N] \
    [--dump/--no-dump] [--out-dir DIR] [--thresh PRESET] \
    [--opt-mode light|heavy] [--hessian-calc-mode Analytical|FiniteDifference] \
    [--config FILE] [--override-yaml FILE] [--show-config] [--dry-run] [--args-yaml FILE]
```

### 例

```bash
# デフォルト heavy モード（RS-I-RFO）
mlmm tsopt -i ts_guess.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --max-cycles 8000 --dump --out-dir ./result_tsopt/

# light モード（Dimer）、解析的ヘシアン
mlmm tsopt -i ts_guess.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --opt-mode light --hessian-calc-mode Analytical --out-dir ./result_tsopt/
```

## ワークフロー

1. **入力処理** -- 酵素 PDB、Amber トポロジー、ML 領域定義を読み込みます。電荷/スピンを解決します。CLI と YAML の凍結原子がマージされます。
2. **ML/MM 計算機の構築** -- ML/MM 計算機（FAIR-Chem UMA + OpenMM）を構築します。`--hessian-calc-mode` は UMA がヘシアンを解析的に評価するか有限差分で評価するかを制御します。
3. **Light モード（Dimer）:**
   - ヘシアン Dimer ステージは、正確なヘシアン（アクティブ部分空間、TR 射影済み）を評価して Dimer 方向を定期的に更新します。
   - フラットンループが有効な場合、保存されたアクティブヘシアンは変位と勾配差分を使用した Bofill で更新されます。各ループで虚数モードを推定し、一度フラットンし、Dimer 方向を更新し、Dimer + LBFGS マイクロセグメントを実行します。
4. **Heavy モード（RS-I-RFO）:**
   - RS-I-RFO オプティマイザーを、`rsirfo` YAML セクションで定義されたオプションのヘシアン参照ファイルとマイクロサイクル制御とともに実行します。
5. **モードエクスポート** -- 収束した虚数モードが `vib/` に `.trj`/`.pdb` ペアとして書き出されます。最終ジオメトリとオプションの軌跡も保存されます。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 開始ジオメトリ（PDB または XYZ）。XYZ の場合はトポロジーに `--ref-pdb` を使用。 | 必須 |
| `--ref-pdb FILE` | 入力が XYZ の場合の参照 PDB トポロジー。 | _None_ |
| `--real-parm7 PATH` | 全酵素の Amber parm7 トポロジー。 | 必須 |
| `--model-pdb PATH` | ML 領域原子を含む PDB。`--detect-layer` 有効時はオプション。 | _None_ |
| `--model-indices TEXT` | ML 領域のカンマ区切り原子インデックス（範囲指定可）。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` を 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer / --no-detect-layer` | 入力 PDB の B 因子から ML/MM レイヤーを検出。 | `True` |
| `-q, --charge INT` | ML 領域の総電荷。 | _None_ |
| `-m, --multiplicity INT` | ML 領域のスピン多重度 (2S+1)。 | _None_ |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切りインデックス（YAML `geom.freeze_atoms` とマージ）。 | _None_ |
| `--hess-cutoff FLOAT` | MM ヘシアン原子の距離カットオフ (A)。カットオフ指定時は `--detect-layer` が無効化。 | _None_ |
| `--movable-cutoff FLOAT` | 可動 MM 原子の距離カットオフ (A)。 | _None_ |
| `--hessian-calc-mode CHOICE` | UMA ヘシアンモード: `Analytical` または `FiniteDifference`。 | _None_ |
| `--max-cycles INT` | 最大総オプティマイザーサイクル。 | `10000` |
| `--dump/--no-dump` | 連結軌跡 `optimization_all.trj` を書き出し。 | `False` |
| `--out-dir TEXT` | 出力ディレクトリ。 | `./result_tsopt/` |
| `--thresh TEXT` | 収束プリセット（`gau_loose\|gau\|gau_tight\|gau_vtight\|baker\|never`）。 | _None_ |
| `--opt-mode CHOICE` | TS オプティマイザーモード: `light`（Dimer）または `heavy`（RS-I-RFO）。 | `heavy` |
| `--partial-hessian-flatten / --full-hessian-flatten` | フラットンループでの虚数モード検出に部分ヘシアン（ML のみ）を使用。 | `True`（部分） |
| `--active-dof-mode CHOICE` | 最終振動解析のアクティブ自由度: `all`、`ml-only`、`partial`、`unfrozen`。 | `partial` |
| `--config FILE` | 明示 CLI オプションより前に適用するベース YAML 設定ファイル。 | _None_ |
| `--override-yaml FILE` | 最後に適用する最優先 YAML 上書きファイル。 | _None_ |
| `--args-yaml FILE` | `--override-yaml` の legacy alias。 | _None_ |
| `--show-config/--no-show-config` | 解決後の設定レイヤーを表示して実行を継続。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに入力/設定を検証し、実行計画を表示。 | `False` |

## 出力

```
out_dir/  (デフォルト: ./result_tsopt/)
├── summary.md                   # 主要成果物のインデックス
├── key_ts.xyz                   # 最終TS構造へのショートカット（または key_ts.pdb）
├── key_imag_mode.trj            # 代表的な虚数モードへのショートカット
├── key_opt.trj                  # 最適化軌跡へのショートカット（存在する場合）
├── final_geometry.xyz            # 常に書き出し
├── final_geometry.pdb            # 入力が PDB の場合
├── optimization_all.trj          # 連結 Dimer セグメント（--dump 時）
├── optimization_all.pdb          # PDB コンパニオン（--dump かつ入力が PDB の場合）
├── vib/
│   ├── final_imag_mode_±XXXX.Xcm-1.trj   # 虚数モード軌跡
│   └── final_imag_mode_±XXXX.Xcm-1.pdb   # 虚数モード PDB コンパニオン
└── .dimer_mode.dat               # Dimer 方向シード（light モード）
```

## YAML 設定（`--config` / `--override-yaml` / `--args-yaml`）

`--config` をベース設定に、`--override-yaml` を最終上書きに使います（`--args-yaml` は `--override-yaml` の legacy alias）。マージ優先順位は次の通りです:

`defaults < config < 明示 CLI < override`。

```yaml
geom:
  coord_type: cart
  freeze_atoms: []
calc:
  charge: 0
  spin: 1
mlmm:
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
opt:
  thresh: baker
  max_cycles: 10000
  dump: false
  out_dir: ./result_tsopt/
hessian_dimer:
  thresh_loose: gau_loose
  thresh: baker
  update_interval_hessian: 500
  neg_freq_thresh_cm: 5.0
  flatten_amp_ang: 0.1
  flatten_max_iter: 50
  root: 0
  dimer:
    length: 0.0189
    rotation_max_cycles: 15
  lbfgs:
    thresh: baker
    max_cycles: 10000
    max_step: 0.3
rsirfo:
  thresh: baker
  max_cycles: 10000
  roots: [0]
  hessian_update: bofill
```

## 注意事項

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- 虚数モード検出はデフォルト閾値 ~5 cm^-1 を使用します（`hessian_dimer.neg_freq_thresh_cm` で設定可能）。選択された `root` がどの虚数モードをエクスポートするかを決定します。
- `--freeze-atoms` は 1 始まりインデックスを受け付け、YAML `geom.freeze_atoms` とマージされます。
- 収束プリセットは外側の管理（`opt`）と内側の LBFGS セグメント（`hessian_dimer.lbfgs`）の両方に伝播されます。
- PHVA 並進/回転射影は `freq` の実装をミラーリングし、アクティブ空間で正しい固有ベクトルを保持しながら GPU メモリ消費を削減します。

---

## 関連項目

- [opt](opt.md) -- 単一構造の構造最適化
- [freq](freq.md) -- 検証済み TS の単一虚数振動数を確認
- [irc](irc.md) -- 最適化された TS からの反応経路追跡
- [all](all.md) -- 抽出 -> MEP -> tsopt -> IRC -> freq を連鎖させるエンドツーエンドワークフロー
