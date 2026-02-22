# ML/MM 計算機

## 概要

> **要約:** PySisyphus 用の ONIOM 型 ML/MM 計算機。FAIR-Chem UMA（高レベル ML）と hessian_ff（低レベル MM）を結合し、酵素活性部位モデルのエネルギー、力、ヘシアンを計算します。

`mlmm_calc.mlmm` は、機械学習原子間ポテンシャル（FAIR-Chem UMA）と分子力学力場（Amber prmtop ベースの `hessian_ff`）を組み合わせた減算型 ONIOM スタイルの ML/MM 計算機を実装しています。`mlmm_toolkit` のすべての ML/MM 最適化、経路探索、スキャン、振動解析、IRC ワークフローのコア計算機として機能します。

この計算機はリンク原子なしで動作します。ML 領域はモデル PDB（`model.pdb`）で定義され、MM トポロジーは Amber prmtop（`real.parm7`）から取得され、座標は入力 PDB（`input.pdb`）から読み取られます。内部 `real.rst7` は ParmEd により `real.parm7` と `input.pdb` の座標を組み合わせて生成されます -- 外部の `real.rst7` や `real.pdb` は不要です。

## 3 層スキーム（エネルギー / 力 / ヘシアン）

計算機は ONIOM 減算法を使用して 3 つの評価を組み合わせます:

| レイヤー | システム | 手法 | 説明 |
| --- | --- | --- | --- |
| **REAL-low** | 全系 | MM (hessian_ff) | Amber prmtop ベースの MM で評価した全系 |
| **MODEL-low** | ML サブセット | MM (hessian_ff) | MM で評価した ML 領域 |
| **MODEL-high** | ML サブセット + リンク H | ML (UMA) | FAIR-Chem UMA で評価した ML 領域 |

結合エネルギーは:

```
E_ONIOM = E(REAL-low) - E(MODEL-low) + E(MODEL-high)
```

力とヘシアンも同じ減算パターンに従います。

## ヘシアン / 最適化用の層設定

実装では 3 層の B 因子と Hessian 対象 MM の別設定を使い、どの原子が MM ヘシアン計算に含まれるか、どの原子が凍結されるかを制御できます:

- **ML 領域**（B 因子 = 0.0）: UMA 機械学習ポテンシャルで処理
- **Movable-MM**（B 因子 = 10.0）: 最適化中に移動する MM 原子
- **Frozen**（B 因子 = 20.0）: 固定された MM 原子
- **Hessian 対象 MM**（専用 B 因子なし）: `hess_cutoff` / `hess_mm_atoms` で選択

レイヤー割り当ては `hess_cutoff`、`movable_cutoff`、`use_bfactor_layers`、および明示的な `*_mm_atoms` リストで制御されます。

## 機能

### リンク原子の再分配

リンク原子からの力とヘシアン寄与はヤコビアンを介して ML/MM 親原子に再分配されます。再分配は以下を追加します:
- セルフ項 `J^T H J`
- ジオメトリ依存の 2 次項 `sum (dJ^T/dx * f_L)` を親原子にインプレースで追加

### MM ヘシアン

MM バックエンドは `mm_backend` パラメータで選択できます：

- **`"hessian_ff"`** (デフォルト): `hessian_ff` による解析的ヘシアン（アクティブ原子のみ。任意で凍結行/列をゼロ埋めした完全デカルト形状に展開）。CPU のみ対応。
- **`"openmm"`**: OpenMM による有限差分 (FD) ヘシアン。CPU と CUDA の両プラットフォームに対応。`hessian_ff` が対応していない力場や、ワークフローに OpenMM を既に使用している場合に有用。

**YAML 設定例:**
```yaml
mlmm:
  mm_backend: openmm  # OpenMM を MM 計算に使用
  mm_device: cuda     # CUDA を使用 (または "cpu")
```

### UMA ヘシアンモード

- `"Analytical"`: 選択されたデバイスでの 2 次自動微分。
- `"FiniteDifference"`: 力の中心差分。

## 入力

| 入力 | 説明 |
| --- | --- |
| `input.pdb` | 入力構造（残基名/原子名がここから読み取られます） |
| `real.parm7` | Amber prmtop（完全 REAL 系のトポロジー） |
| `model.pdb` | ML 領域を定義する PDB（原子 ID の決定に使用） |

## 単位

| 量 | 内部単位 | PySisyphus インターフェース |
| --- | --- | --- |
| エネルギー | eV | Hartree |
| 力 | eV/A | Hartree/Bohr |
| ヘシアン | eV/A^2 | Hartree/Bohr^2 |

PySisyphus インターフェースは原子単位（Hartree/Bohr）に変換された値を返します。

## 実装ノート

- メモリ使用量を削減するため、ヘシアン組み立てはインプレースの `add_`/`sub_`/`mul_` を使用します。
- REAL/MODEL 間の原子 ID マッピングは `input.pdb` のみに基づきます。
- `real.rst7` は ParmEd を使用して `input.pdb` の座標から内部的に作成されます。

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [opt](opt.md) -- ML/MM 計算機を使用した単一構造の構造最適化
- [tsopt](tsopt.md) -- 遷移状態最適化
- [freq](freq.md) -- 振動解析
- [YAML リファレンス](yaml_reference.md) -- `calc`/`mlmm` 設定キー
