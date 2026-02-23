# `irc`

## 概要

> **要約:** ML/MM 計算機を使用した EulerPC 予測子-補正子積分器により、遷移状態から反応物と生成物の方向へ固有反応座標（IRC）を追跡します。

`mlmm irc` は EulerPC 積分器を使用して IRC 計算を実行します。CLI は意図的に狭く設計されており、コマンドラインに表面化されていないパラメータは YAML で提供し、実行を明示的かつ再現可能に保つべきです。入力は `pysisyphus.helpers.geom_loader` で読み取り可能な任意の構造（`.pdb`、`.xyz`、`_trj.xyz`、...）です。入力が `.pdb` の場合、生成される軌跡は追加で PDB に変換されます。


## 最小例

```bash
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc
```

```bash
# YAML で最小設定（real_parm7/model_pdb を config に記述）
cat > irc_min.yaml << 'YAML'
calc:
 real_parm7: real.parm7
 model_pdb: ml_region.pdb
 use_bfactor_layers: false
YAML
mlmm irc -i ts.pdb -q 0 -m 1 --config irc_min.yaml \
 --max-cycles 50 --out-dir ./result_irc_yaml
```

## 出力の見方

- `result_irc/summary.md`
- `result_irc/key_irc_trj.xyz`
- `result_irc/key_irc_forward_trj.xyz`
- `result_irc/finished_irc_trj.xyz`
- `result_irc/forward_irc_trj.xyz`

## よくある例

1. 正方向のみを実行する。

```bash
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --out-dir ./result_irc_forward
```

2. ステップサイズと root を調整する。

```bash
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --step-size 0.20 --root 1 \
 --out-dir ./result_irc_step
```

3. 有限差分ヘシアンで確認する。

```bash
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --hessian-calc-mode FiniteDifference \
 --max-cycles 80 --out-dir ./result_irc_fd
```

## 使用法

```bash
mlmm irc -i INPUT.pdb --real-parm7 real.parm7
 [--model-pdb ml_region.pdb | --model-indices "1,2,3" | --detect-layer]
 [-q CHARGE] [-m MULT]
 [--max-cycles N] [--step-size Ds] [--root k]
 [--detect-layer/--no-detect-layer]
 [--out-dir DIR]
 [--hessian-calc-mode Analytical|FiniteDifference]
 [--show-config] [--dry-run]
```

### 例

```bash
# 有限差分ヘシアンとカスタムステップサイズによる正方向のみの IRC
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --step-size 0.2 --hessian-calc-mode FiniteDifference --out-dir ./irc_fd/

# PDB 入力を使用（軌跡も PDB としてエクスポート）
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc/
```

## ワークフロー

1. **入力準備** -- `geom_loader` でサポートされる任意の形式を受け付けます。入力が `.pdb` の場合、軌跡も PDB 形式に変換されます。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 構造ファイル（`.pdb`/`.xyz`/`_trj.xyz`/...）。 | 必須 |
| `--real-parm7 PATH` | 全酵素/MM 領域の Amber トポロジー。YAML の `calc.real_parm7` が無い場合は必須。 | _None_ |
| `--model-pdb PATH` | ML 領域を定義する PDB。`--no-detect-layer` かつ `--model-indices` 未指定時は必須。 | _None_ |
| `--model-indices TEXT` | ML 領域原子インデックス（カンマ区切り、範囲指定可: `1-10,15`）。`--model-pdb` 省略時に使用。 | _None_ |
| `--model-indices-one-based/--model-indices-zero-based` | `--model-indices` を 1 始まり/0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer/--no-detect-layer` | 入力 PDB の B 因子（`B=0/10/20`）から ML/MM レイヤーを検出。 | `True` |
| `-q, --charge INT` | 総電荷。YAML の `calc.charge` を上書き。 | 強く推奨 |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。`calc.spin` を上書き。 | `1` |
| `--max-cycles INT` | IRC ステップの最大数。`irc.max_cycles` を上書き。 | _デフォルト_ |
| `--step-size FLOAT` | 質量加重座標でのステップ長。`irc.step_length` を上書き。 | _デフォルト_ |
| `--root INT` | 初期変位の虚数モードインデックス。`irc.root` を上書き。 | _デフォルト_ |
| `--forward/--no-forward` | 正方向 IRC を実行。`irc.forward` を上書き。 | `True` |
| `--out-dir PATH` | 出力ディレクトリ。`irc.out_dir` を上書き。 | `./result_irc/` |
| `--hessian-calc-mode CHOICE` | UMA がヘシアンを構築する方法（`Analytical` または `FiniteDifference`）。`calc.hessian_calc_mode` を上書き。 | _デフォルト_ |
| `--config FILE` | 明示CLI適用前に読み込むベース YAML。 | _None_ |
| `--show-config/--no-show-config` | 解決済み YAML レイヤー/設定を表示して続行。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画のみ表示。 | `False` |

## 出力

```text
out_dir/ (デフォルト:./result_irc/)
 summary.md # 主要成果物のインデックス
 key_irc_trj.xyz # finished_irc_trj.xyz へのショートカット
 key_irc_forward_trj.xyz # forward_irc_trj.xyz へのショートカット
 key_irc.pdb # finished_irc.pdb へのショートカット（存在時）
 key_irc_data.h5 # irc_data.h5 へのショートカット（存在時）
 <prefix>irc_data.h5 # irc.dump_every ステップごとに書き出される HDF5 ダンプ
 <prefix>finished_irc_trj.xyz # 完全 IRC 軌跡（XYZ/TRJ）
 <prefix>forward_irc_trj.xyz # 正方向パスセグメント
 <prefix>finished_irc.pdb # PDB 変換（入力が.pdb の場合のみ）
 <prefix>forward_irc.pdb # PDB 変換（入力が.pdb の場合のみ）
```

## YAML 設定

セクション `geom`、`calc`、`irc` を含む YAML マッピングを提供します。マージ順は **デフォルト < config < 明示CLI < override** です（

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

## 注意事項

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- 電荷/多重度の運用ルールは [CLI Conventions](cli_conventions.md) に集約しています。
- UMA オプションは mlmm 計算機に直接渡されます。`device: "auto"` の場合、計算機は GPU/CPU を自動選択します。
- `hessian_calc_mode: "FiniteDifference"` の場合、`geom.freeze_atoms` を使用して FD ヘシアン構築で凍結自由度をスキップできます。
- `--step-size` は質量加重座標です。`--root` は初期変位に使用する虚数振動数インデックスを選択します。
- 標準出力には進捗とタイミングが含まれます。終了コード: 成功時 `0`、`KeyboardInterrupt` 時 `130`、未処理例外時 `1`。

---

## 関連項目

- [tsopt](tsopt.md) -- IRC 実行前に TS を最適化
- [freq](freq.md) -- TS 候補が 1 つの虚数振動数を持つことを検証
- [all](all.md) -- tsopt の後に IRC を実行するエンドツーエンドワークフロー
