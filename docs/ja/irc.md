# `irc`

## 概要

> **要約:** ML/MM 計算機を使用した EulerPC 予測子-補正子積分器により、遷移状態から反応物と生成物の方向へ固有反応座標（IRC）を追跡します。

`mlmm irc` は EulerPC 積分器を使用して IRC 計算を実行します。CLI は意図的に狭く設計されており、コマンドラインに表面化されていないパラメータは YAML で提供し、実行を明示的かつ再現可能に保つべきです。入力は `pysisyphus.helpers.geom_loader` で読み取り可能な任意の構造（`.pdb`、`.xyz`、`.trj`、...）です。入力が `.pdb` の場合、生成される軌跡は追加で PDB に変換されます。

設定の優先順位は: **組み込みデフォルト -> YAML -> CLI** です。

## 使用法

```bash
mlmm irc -i INPUT.pdb -q CHARGE [-m MULT]
    [--max-cycles N] [--step-size Ds] [--root k]
    [--forward True|False] [--backward True|False]
    [--out-dir DIR]
    [--hessian-calc-mode Analytical|FiniteDifference]
    [--args-yaml FILE]
```

### 例

```bash
# 有限差分ヘシアンとカスタムステップサイズによる正方向のみの IRC
mlmm irc -i ts.xyz -q -1 -s 2 --forward True --backward False \
    --step-size 0.2 --hessian-calc-mode FiniteDifference --out-dir ./irc_fd/

# PDB 入力を使用（軌跡も PDB としてエクスポート）
mlmm irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc/
```

## ワークフロー

1. **入力準備** -- `geom_loader` でサポートされる任意の形式を受け付けます。入力が `.pdb` の場合、軌跡も PDB 形式に変換されます。
2. **設定のマージ** -- デフォルト -> CLI -> YAML、セクション `geom`、`calc`、`irc` をマージします。`-q/--charge` と `-m/--multiplicity` の両方を明示的に指定することを強く推奨します。
3. **IRC 積分** -- EulerPC が `irc.forward/backward`、`irc.step_length`、`irc.root`、および UMA で設定されたヘシアンワークフローに従って正方向/逆方向の分岐を積分します。
4. **出力** -- 軌跡（`finished`、`forward`、`backward`）が `.trj` として書き出され、入力が `.pdb` の場合は追加で `.pdb` としても書き出されます。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 構造ファイル（`.pdb`/`.xyz`/`.trj`/...）。 | 必須 |
| `-q, --charge INT` | 総電荷。YAML の `calc.charge` を上書き。 | 強く推奨 |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。`calc.spin` を上書き。 | `1` |
| `--max-cycles INT` | IRC ステップの最大数。`irc.max_cycles` を上書き。 | _デフォルト_ |
| `--step-size FLOAT` | 質量加重座標でのステップ長。`irc.step_length` を上書き。 | _デフォルト_ |
| `--root INT` | 初期変位の虚数モードインデックス。`irc.root` を上書き。 | _デフォルト_ |
| `--forward {True\|False}` | 正方向 IRC を実行。`irc.forward` を上書き。 | `True` |
| `--backward {True\|False}` | 逆方向 IRC を実行。`irc.backward` を上書き。 | `True` |
| `--out-dir PATH` | 出力ディレクトリ。`irc.out_dir` を上書き。 | `./result_irc/` |
| `--hessian-calc-mode CHOICE` | UMA がヘシアンを構築する方法（`Analytical` または `FiniteDifference`）。`calc.hessian_calc_mode` を上書き。 | _デフォルト_ |
| `--args-yaml FILE` | セクション `geom`、`calc`、`irc` を含む YAML ファイル。 | _None_ |

## 出力

```text
out_dir/ (デフォルト: ./result_irc/)
  <prefix>irc_data.h5              # irc.dump_every ステップごとに書き出される HDF5 ダンプ
  <prefix>finished_irc.trj         # 完全 IRC 軌跡（XYZ/TRJ）
  <prefix>forward_irc.trj          # 正方向パスセグメント
  <prefix>backward_irc.trj         # 逆方向パスセグメント
  <prefix>finished_irc.pdb         # PDB 変換（入力が .pdb の場合のみ）
  <prefix>forward_irc.pdb          # PDB 変換（入力が .pdb の場合のみ）
  <prefix>backward_irc.pdb         # PDB 変換（入力が .pdb の場合のみ）
```

## YAML 設定

セクション `geom`、`calc`、`irc` を含む YAML マッピングを提供します。YAML 値が CLI を上書きします。

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

## 注意事項

- **強い推奨:** 意図しない条件での実行を避けるため、`-q/--charge` と `-m/--multiplicity` の両方を明示的に指定してください。
- UMA オプションは mlmm 計算機に直接渡されます。`device: "auto"` の場合、計算機は GPU/CPU を自動選択します。
- `hessian_calc_mode: "FiniteDifference"` の場合、`geom.freeze_atoms` を使用して FD ヘシアン構築で凍結自由度をスキップできます。
- `--step-size` は質量加重座標です。`--root` は初期変位に使用する虚数振動数インデックスを選択します。
- 標準出力には進捗とタイミングが含まれます。終了コード: 成功時 `0`、`KeyboardInterrupt` 時 `130`、未処理例外時 `1`。

---

## 関連項目

- [tsopt](tsopt.md) -- IRC 実行前に TS を最適化
- [freq](freq.md) -- TS 候補が 1 つの虚数振動数を持つことを検証
- [all](all.md) -- tsopt の後に IRC を実行するエンドツーエンドワークフロー
