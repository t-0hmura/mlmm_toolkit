# `irc`

## 概要

> **要約:** ML/MM 計算機を使用した EulerPC ベースの IRC（固有反応座標）積分により、遷移状態から反応物と生成物の方向へ追跡します。デフォルトでは正方向と逆方向の両方のブランチが計算されます。

- **用途:** 最適化された TS があり、ML/MM で反応物・生成物方向への最小エネルギー経路を追跡したい場合。
- **手法:** 完全 ML/MM ヘシアン（MLIP バックエンド（デフォルト: UMA）+ hessian_ff）による EulerPC 予測子-補正子積分器。
- **出力:** `finished_irc_trj.xyz`、`forward_irc_trj.xyz`、PDB 入力時は `.pdb` コンパニオン。
- **次のステップ:** IRC 端点で [freq](freq.md) を実行し、[opt](opt.md) で真の極小に精密化。

`mlmm irc` は EulerPC 積分器を使用して IRC 計算を実行します。CLI は意図的に狭く設計されており、コマンドラインに表面化されていないパラメータは YAML で提供し、実行を明示的かつ再現可能に保つべきです。入力は `pysisyphus.helpers.geom_loader` で読み取り可能な任意の構造（`.pdb`、`.xyz`、`_trj.xyz`、...）です。入力が `.pdb` の場合、生成される軌跡は追加で PDB に変換されます。

典型的なワークフローは `tsopt` -> `freq`（**1 つ**の虚振動数モードを確認）-> `irc` です。

## 最小例

```bash
mlmm irc -i ts.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc
```

## 出力の見方

- `result_irc/finished_irc_trj.xyz`
- `result_irc/forward_irc_trj.xyz`

## よくある例

1. 正方向のみを実行する。

```bash
mlmm irc -i ts.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --no-backward --out-dir ./result_irc_forward
```

2. ステップサイズを増やして解析的ヘシアンを使う。

```bash
mlmm irc -i ts.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --step-size 0.20 \
 --hessian-calc-mode Analytical --out-dir ./result_irc_analytical
```

3. 両ブランチを保持してステップ上限を引き上げる。

```bash
mlmm irc -i ts.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --max-cycles 150 \
 --out-dir ./result_irc_long
```

## ワークフロー

1. **入力準備** -- `geom_loader` でサポートされる任意の形式を受け付けます。参照 PDB が利用可能な場合（入力が `.pdb` または `--ref-pdb` 指定時）、EulerPC 軌跡はそのトポロジーを使用して PDB に変換されます。
2. **ML/MM 計算機の構築** -- `--parm` と `--model-pdb` から ML/MM 計算機を構築します。`-b/--backend` で ML バックエンドを選択し（デフォルト: `uma`）、`--hessian-calc-mode` は MLIP ヘシアン評価を制御します。`--embedcharge` で xTB 点電荷埋め込み補正を有効化できます。
3. **IRC 積分** -- EulerPC 積分器が両方向に沿って IRC を伝播します（`--no-forward` または `--no-backward` でブランチを無効化可能）。ステップサイズとサイクル数で積分長を制御します。
4. **出力と変換** -- 軌跡は XYZ で書き出されます。PDB テンプレートが利用可能で `--convert-files` が有効な場合、PDB コンパニオンが生成されます。

## CLIオプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-b, --backend CHOICE` | ML バックエンド: `uma`（デフォルト）、`orb`、`mace`、`aimnet2`。 | `uma` |
| `--embedcharge/--no-embedcharge` | xTB 点電荷埋め込み補正の有効化。MM 環境から ML 領域への静電的影響を考慮。 | `False` |
| `--embedcharge-cutoff FLOAT` | xTB 埋め込み用 MM 原子のカットオフ半径（Å）。 | `12.0` |
| `-i, --input PATH` | 構造ファイル（`.pdb`/`.xyz`/`_trj.xyz`/...）。 | 必須 |
| `--parm PATH` | 全酵素/MM 領域の Amber トポロジー。YAML の `calc.real_parm7` が無い場合は必須。 | _None_ |
| `--model-pdb PATH` | ML 領域を定義する PDB。`--no-detect-layer` かつ `--model-indices` 未指定時は必須。 | _None_ |
| `--model-indices TEXT` | ML 領域原子インデックス（カンマ区切り、範囲指定可: `1-10,15`）。`--model-pdb` 省略時に使用。 | _None_ |
| `--model-indices-one-based/--model-indices-zero-based` | `--model-indices` を 1 始まり/0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer/--no-detect-layer` | 入力 PDB の B 因子（`B=0/10/20`）から ML/MM レイヤーを検出。 | `True` |
| `-q, --charge INT` | 総電荷。YAML の `calc.charge` を上書き。 | _None_（`-l` 未指定時は必須） |
| `-l, --ligand-charge TEXT` | 残基ごとの電荷マッピング（例: `GPP:-3,SAM:1`）。`-q` 省略時に合計電荷を導出。 | _None_ |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。`calc.spin` を上書き。 | `1` |
| `--max-cycles INT` | IRC ステップの最大数。`irc.max_cycles` を上書き。 | `125` |
| `--step-size FLOAT` | 質量加重座標でのステップ長。`irc.step_length` を上書き。 | `0.10` |
| `--root INT` | 初期変位の虚振動数モードインデックス。`irc.root` を上書き。 | `0` |
| `--forward/--no-forward` | 正方向 IRC を実行。`irc.forward` を上書き。 | `True` |
| `--backward/--no-backward` | 逆方向 IRC を実行。`irc.backward` を上書き。 | `True` |
| `-o, --out-dir PATH` | 出力ディレクトリ。`irc.out_dir` を上書き。 | `./result_irc/` |
| `--ref-pdb FILE` | `--input` が XYZ の場合に使用する参照 PDB トポロジー（XYZ 座標を保持）。 | _None_ |
| `--convert-files/--no-convert-files` | 参照 PDB 利用可能時の XYZ/TRJ から PDB コンパニオンの切り替え。 | `True` |
| `--hessian-calc-mode CHOICE` | MLIP がヘシアンを構築する方法（`Analytical` または `FiniteDifference`）。`calc.hessian_calc_mode` を上書き。 | `FiniteDifference` |
| `--config FILE` | 明示 CLI 適用前に読み込むベース YAML。 | _None_ |
| `--show-config/--no-show-config` | 解決済み YAML レイヤー/設定を表示して続行。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画のみ表示。`--help-advanced` に表示。 | `False` |

## 出力

```
out_dir/ (デフォルト: ./result_irc/)
├─ <prefix>irc_data.h5              # irc.dump_every ステップごとに書き出される HDF5 ダンプ
├─ <prefix>finished_irc_trj.xyz     # 完全 IRC 軌跡（XYZ/TRJ）
├─ <prefix>forward_irc_trj.xyz      # 正方向パスセグメント
├─ <prefix>backward_irc_trj.xyz     # 逆方向パスセグメント
├─ <prefix>finished_irc.pdb         # PDB 変換（入力が .pdb の場合のみ）
├─ <prefix>forward_irc.pdb          # PDB 変換（入力が .pdb の場合のみ）
├─ <prefix>backward_irc.pdb         # PDB 変換（入力が .pdb の場合のみ）
├─ <prefix>forward_last.xyz         # 正方向 IRC 終点（XYZ、単一フレーム）
├─ <prefix>forward_last.pdb         # 正方向 IRC 終点（PDB、利用可能時）
├─ <prefix>backward_last.xyz        # 逆方向 IRC 終点（XYZ、単一フレーム）
└─ <prefix>backward_last.pdb        # 逆方向 IRC 終点（PDB、利用可能時）
```

## YAML設定

マージ順 **デフォルト < config < 明示CLI < override** でマッピングを提供します。
共有セクションはジオメトリ/計算機キーについて [YAML リファレンス](yaml_reference.md) を再利用します。`irc` では YAML/CLI マージ後に `geom.coord_type` が `cart` に強制されます。`calc.return_partial_hessian` は `true` に強制されます（partial Hessian、active-DOF 処理）。

### CLI から YAML へのマッピング

| CLI オプション | YAML キー |
|------------|----------|
| `--charge` | `calc.charge` |
| `--multiplicity` | `calc.spin` |
| `--step-size` | `irc.step_length` |
| `--max-cycles` | `irc.max_cycles` |
| `--root` | `irc.root` |
| `--forward` | `irc.forward` |
| `--backward` | `irc.backward` |
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
 backend: uma                      # ML バックエンド (uma/orb/mace/aimnet2)
 embedcharge: false                # xTB 点電荷埋め込み補正
 uma_model: uma-s-1p1              # uma-s-1p1 | uma-m-1p1
 uma_task_name: omol                # UMA タスク名 (backend=uma 時)
 ml_device: auto                   # ML デバイス選択
 ml_hessian_mode: Analytical         # ヘシアンモード選択
 return_partial_hessian: true      # irc では true に強制（partial Hessian、active-DOF 処理）
irc:
 step_length: 0.1                  # 積分ステップ長
 max_cycles: 125                   # IRC に沿った最大ステップ数
 downhill: false                   # 下り方向のみに追従
 forward: true                     # 正方向に伝播
 backward: true                    # 逆方向に伝播
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

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細なトラブルシューティングガイド

- [tsopt](tsopt.md) -- IRC 実行前に TS を最適化
- [freq](freq.md) -- TS 候補が 1 つの虚数振動数を持つことを検証; IRC 端点を解析
- [opt](opt.md) -- IRC 端点を真の極小に最適化
- [all](all.md) -- tsopt の後に IRC を実行するend-to-endワークフロー
- [YAML リファレンス](yaml_reference.md) -- `irc` の完全な設定オプション
- [用語集](glossary.md) -- IRC（固有反応座標）の定義
