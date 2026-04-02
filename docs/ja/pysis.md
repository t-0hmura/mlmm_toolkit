# mlmm pysis — pysisyphus YAMLワークフロー

> **概要:** mlmm-toolkit calculatorを事前登録した状態でpysisyphus YAMLワークフローファイルを実行。v0.1.x YAMLベースのワークフローと互換。

## 概要

`mlmm pysis`サブコマンドは、`mlmm` calculator typeが自動登録された状態でpysisyphus YAMLワークフローファイルを実行します。v0.1.xで`mlmm opt.yaml`が主要インターフェースだったYAMLベースのワークフローとの後方互換性を提供します。

```bash
mlmm pysis opt.yaml
mlmm pysis tsopt.yaml
mlmm pysis irc.yaml
```

## YAMLフォーマット

pysisyphus YAMLファイルは3つの主要セクション（`geom`、`calc`、ワークフローセクション（`opt`、`tsopt`、`irc`など））で構成されます。

### ジオメトリセクション

```yaml
geom:
  type: cart          # 座標系（cart, redund, dlc）
  fn: complex_layered.pdb
```

凍結原子を指定する場合は`freeze_atoms`（1始まりインデックス）を使用:

```yaml
geom:
  type: cart
  fn: complex_layered.pdb
  freeze_atoms: [101, 102, 103, 200, 201]  # 1始まり
```

### Calculatorセクション

```yaml
calc:
  type: mlmm
  input_pdb: complex_layered.pdb
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
  model_charge: 0
  model_mult: 1
  backend: uma
  embedcharge: false
```

`MLMMCore`の全パラメータをYAMLキーとして使用可能です。v0.1.xのパラメータ名（`real_pdb`など）も非推奨警告付きで受け付けます。

### ワークフローセクション

#### 構造最適化

```yaml
opt:
  type: lbfgs
  max_cycles: 300
  thresh: gau_loose    # gau, gau_tight, gau_vtight, gau_loose, baker
  dump: true
```

#### 遷移状態最適化

```yaml
tsopt:
  type: dimer
  max_cycles: 100
  thresh: gau_loose
  dump: true
```

#### IRC

```yaml
irc:
  type: eulerpc
  max_cycles: 75
  step_length: 0.15
  dump: true
```

## 完全なYAML例

### 構造最適化

```yaml
geom:
  type: cart
  fn: complex_layered.pdb
  freeze_atoms: [101, 102, 200, 201]

calc:
  type: mlmm
  input_pdb: complex_layered.pdb
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
  model_charge: 0
  backend: uma

opt:
  type: lbfgs
  max_cycles: 300
  thresh: gau_loose
  dump: true
```

実行: `mlmm pysis opt.yaml`

### 遷移状態最適化

```yaml
geom:
  type: cart
  fn: ts_candidate.pdb
  freeze_atoms: [101, 102, 200, 201]

calc:
  type: mlmm
  input_pdb: ts_candidate.pdb
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
  model_charge: 0

tsopt:
  type: dimer
  max_cycles: 100
  thresh: gau_loose
  dump: true
```

実行: `mlmm pysis tsopt.yaml`

### TSからのIRC

```yaml
geom:
  type: cart
  fn: ts_optimized.pdb

calc:
  type: mlmm
  input_pdb: ts_optimized.pdb
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
  model_charge: 0

irc:
  type: eulerpc
  max_cycles: 75
  step_length: 0.15
  dump: true
```

実行: `mlmm pysis irc.yaml`

## v0.1.xからの移行

| v0.1.x | v0.2.x |
|--------|--------|
| `mlmm opt.yaml` | `mlmm pysis opt.yaml` |
| `real_pdb: complex.pdb` | `input_pdb: complex.pdb`（`real_pdb`も警告付きで使用可） |
| `real_rst7: real.rst7` | 不要（警告付きで無視） |

既存のv0.1.x YAMLファイルは変更なしで動作します — 非推奨パラメータ名は警告付きで受け付けます。

## 関連ドキュメント

- [Python API](python-api.md) — mlmm-toolkitをPythonライブラリとして使用
- [ML/MM Calculator](mlmm-calc.md) — Calculatorアーキテクチャ
- [YAML Reference](yaml-reference.md) — YAML設定リファレンス
