# `irc`

## 概要

> **概要:** ML/MM 計算機を使用した EulerPC ベースの IRC（固有反応座標）積分により、遷移状態から反応物と生成物の方向へ追跡します。デフォルトでは正方向と逆方向の両方のブランチが計算されます。

### 概要
- **用途:** 最適化された TS があり、ML/MM で反応物・生成物方向への最小エネルギー経路を追跡したい場合。
- **手法:** 完全 ML/MM ヘシアン（FAIR-Chem UMA + hessian_ff）による EulerPC 予測子-補正子積分器。
- **出力:** `finished_irc_trj.xyz`、`forward_irc_trj.xyz`、PDB 入力時は `.pdb` コンパニオン。
- **デフォルト:** `--max-cycles 125`、`--step-size 0.10`、`--freeze-links` 有効、`--convert-files` 有効。
- **次のステップ:** IRC 端点で [freq](freq.md) を実行し、[opt](opt.md) で真の極小に精密化。

`mlmm irc` は EulerPC 積分器を使用して IRC 計算を実行します。CLI は意図的に狭く設計されており、コマンドラインに表面化されていないパラメータは YAML で提供し、実行を明示的かつ再現可能に保つべきです。入力は `pysisyphus.helpers.geom_loader` で読み取り可能な任意の構造（`.pdb`、`.xyz`、`_trj.xyz`、...）です。入力が `.pdb` の場合、生成される軌跡は追加で PDB に変換されます。

典型的なワークフローは `tsopt` -> `freq`（**1 つ**の虚数モードを確認）-> `irc` です。

## 最小の例

```bash
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc
```

## 出力チェックリスト

- `result_irc/finished_irc_trj.xyz`
- `result_irc/forward_irc_trj.xyz`

## 使用例

1. 正方向のみを実行する。

```bash
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --no-forward --out-dir ./result_irc_forward
```

2. ステップサイズを増やして解析的ヘシアンを使う。

```bash
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --step-size 0.20 \
 --hessian-calc-mode Analytical --out-dir ./result_irc_analytical
```

3. 両ブランチを保持してステップ上限を引き上げる。

```bash
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --max-cycles 150 \
 --out-dir ./result_irc_long
```

## ワークフロー

1. **入力準備** -- `geom_loader` でサポートされる任意の形式を受け付けます。参照 PDB が利用可能な場合（入力が `.pdb` または `--ref-pdb` 指定時）、EulerPC 軌跡はそのトポロジーを使用して PDB に変換されます。
2. **ML/MM 計算機の構築** -- `--real-parm7` と `--model-pdb` から ML/MM 計算機を構築します。`--hessian-calc-mode` は UMA ヘシアン評価を制御します。
3. **リンク凍結検出** -- `--freeze-links`（デフォルト）により、リンク水素の親原子が検出され `geom.freeze_atoms` にマージされます。
4. **IRC 積分** -- EulerPC 積分器が両方向に沿って IRC を伝播します（`--no-forward` で正方向を無効化可能）。ステップサイズとサイクル数で積分長を制御します。
5. **出力と変換** -- 軌跡は XYZ で書き出されます。PDB テンプレートが利用可能で `--convert-files` が有効な場合、PDB コンパニオンが生成されます。

## CLIオプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 構造ファイル（`.pdb`/`.xyz`/`_trj.xyz`/...）。 | 必須 |
| `--real-parm7 PATH` | 全酵素/MM 領域の Amber トポロジー。YAML の `calc.real_parm7` が無い場合は必須。 | _None_ |
| `--model-pdb PATH` | ML 領域を定義する PDB。`--no-detect-layer` かつ `--model-indices` 未指定時は必須。 | _None_ |
| `--model-indices TEXT` | ML 領域原子インデックス（カンマ区切り、範囲指定可: `1-10,15`）。`--model-pdb` 省略時に使用。 | _None_ |
| `--model-indices-one-based/--model-indices-zero-based` | `--model-indices` を 1 始まり/0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer/--no-detect-layer` | 入力 PDB の B 因子（`B=0/10/20`）から ML/MM レイヤーを検出。 | `True` |
| `-q, --charge INT` | 総電荷。YAML の `calc.charge` を上書き。 | 必須 |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。`calc.spin` を上書き。 | `1` |
| `--max-cycles INT` | IRC ステップの最大数。`irc.max_cycles` を上書き。 | `125` |
| `--step-size FLOAT` | 質量加重座標でのステップ長。`irc.step_length` を上書き。 | `0.10` |
| `--root INT` | 初期変位の虚数モードインデックス。`irc.root` を上書き。 | `0` |
| `--forward/--no-forward` | 正方向 IRC を実行。`irc.forward` を上書き。 | `True` |
| `--freeze-links/--no-freeze-links` | PDB 入力専用。リンク水素の親原子を凍結（`geom.freeze_atoms` にマージ）。 | `True` |
| `--out-dir PATH` | 出力ディレクトリ。`irc.out_dir` を上書き。 | `./result_irc/` |
| `--convert-files/--no-convert-files` | 参照 PDB 利用可能時の XYZ/TRJ から PDB コンパニオンの切り替え。 | `True` |
| `--hessian-calc-mode CHOICE` | UMA がヘシアンを構築する方法（`Analytical` または `FiniteDifference`）。`calc.hessian_calc_mode` を上書き。 | _デフォルト_ |
| `--config FILE` | 明示 CLI 適用前に読み込むベース YAML。 | _None_ |
| `--show-config/--no-show-config` | 解決済み YAML レイヤー/設定を表示して続行。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画のみ表示。 | `False` |

## 出力

```
out_dir/ (デフォルト: ./result_irc/)
├─ <prefix>irc_data.h5              # irc.dump_every ステップごとに書き出される HDF5 ダンプ
├─ <prefix>finished_irc_trj.xyz     # 完全 IRC 軌跡（XYZ/TRJ）
├─ <prefix>forward_irc_trj.xyz      # 正方向パスセグメント
├─ <prefix>finished_irc.pdb         # PDB 変換（入力が .pdb の場合のみ）
└─ <prefix>forward_irc.pdb          # PDB 変換（入力が .pdb の場合のみ）
```

## YAML設定

マージ順 **デフォルト < config < 明示CLI < override** でマッピングを提供します。
共有セクションはジオメトリ/計算機キーについて [YAML リファレンス](yaml_reference.md) を再利用します。`irc` では YAML/CLI マージ後に `geom.coord_type` が `cart` に、`calc.return_partial_hessian` が `false` に強制されます。

### CLI から YAML へのマッピング

| CLI オプション | YAML キー |
|------------|----------|
| `--charge` | `calc.charge` |
| `--multiplicity` | `calc.spin` |
| `--step-size` | `irc.step_length` |
| `--max-cycles` | `irc.max_cycles` |
| `--root` | `irc.root` |
| `--forward` | `irc.forward` |
| `--out-dir` | `irc.out_dir` |
| `--hessian-calc-mode` | `calc.hessian_calc_mode` |

### YAML 例

```yaml
geom:
 coord_type: cart                  # irc では cart に強制（YAML 値は無視）
 freeze_atoms: []                  # 0 始まり凍結原子（CLI/リンク検出とマージ）
calc:
 charge: 0                         # 総電荷（CLI 上書き）
 spin: 1                           # スピン多重度 2S+1
mlmm:
 real_parm7: real.parm7            # Amber parm7 トポロジー
 model_pdb: ml_region.pdb          # ML 領域定義
 uma_model: uma-s-1p1              # UMA モデルタグ
 uma_task_name: omol                # UMA タスク名
 ml_device: auto                   # UMA デバイス選択
 ml_hessian_mode: FiniteDifference  # ヘシアンモード選択
 return_partial_hessian: false     # irc では false に強制（完全ヘシアン）
irc:
 step_length: 0.1                  # 積分ステップ長
 max_cycles: 125                   # IRC に沿った最大ステップ数
 downhill: false                   # 下り方向のみに追従
 forward: true                     # 正方向に伝播
 root: 0                           # 基準振動ルートインデックス
 hessian_init: calc                # ヘシアン初期化ソース
 displ: energy                     # 変位構築方法
 displ_energy: 0.001               # エネルギーベースの変位スケーリング
 displ_length: 0.1                 # 長さベースの変位フォールバック
 rms_grad_thresh: 0.001            # RMS 勾配収束閾値
 hard_rms_grad_thresh: null        # ハード RMS 勾配停止
 energy_thresh: 0.000001           # エネルギー変化閾値
 imag_below: 0.0                   # 虚数振動数カットオフ
 force_inflection: true            # 変曲点検出を強制
 check_bonds: false                # 伝播中の結合チェック
 out_dir: ./result_irc/            # 出力ディレクトリ
 prefix: ""                        # ファイル名プレフィックス
 hessian_update: bofill            # ヘシアン更新方式
 hessian_recalc: null              # ヘシアン再構築間隔
 max_pred_steps: 500               # 予測子-補正子の最大ステップ数
 loose_cycles: 3                   # 厳密化前のゆるいサイクル数
 corr_func: mbs                    # 相関関数の選択
```

## 注意事項

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- 電荷/多重度の運用ルールは [CLI Conventions](cli_conventions.md) に集約しています。
- UMA オプションは mlmm 計算機に直接渡されます。`device: "auto"` の場合、計算機は GPU/CPU を自動選択します。
- VRAM に余裕がある場合は `--hessian-calc-mode` を `Analytical` に設定することを強く推奨します。
- `irc` は partial-first です。YAML で `calc.return_partial_hessian` を明示しない場合、初期ヘシアンは既定で部分ヘシアンになります。完全ヘシアンを使う場合は `calc.return_partial_hessian: false` を明示してください。
- `--freeze-links` は PDB 入力にのみ適用され、ヘシアン構築中にリンク水素の親原子を凍結します。
- `hessian_calc_mode: "FiniteDifference"` の場合、`geom.freeze_atoms` を使用して FD ヘシアン構築で凍結自由度をスキップできます。
- `--step-size` は質量加重座標です。`--root` は初期変位に使用する虚数振動数インデックスを選択します。
- 標準出力には進捗とタイミングが含まれます。終了コード: 成功時 `0`、`KeyboardInterrupt` 時 `130`、未処理例外時 `1`。

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細なトラブルシューティングガイド

- [tsopt](tsopt.md) -- IRC 実行前に TS を最適化
- [freq](freq.md) -- TS 候補が 1 つの虚数振動数を持つことを検証; IRC 端点を解析
- [opt](opt.md) -- IRC 端点を真の極小に最適化
- [all](all.md) -- tsopt の後に IRC を実行するエンドツーエンドワークフロー
- [YAML リファレンス](yaml_reference.md) -- `irc` の完全な設定オプション
- [用語集](glossary.md) -- IRC（固有反応座標）の定義
