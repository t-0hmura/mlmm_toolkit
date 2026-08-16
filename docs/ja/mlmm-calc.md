# ML/MM calculator

## 概要

PySisyphus 用の ONIOM 型 ML/MM calculatorです。MLIP バックエンド（デフォルト: FAIR-Chem UMA、選択肢: `orb`、`mace`、`aimnet2`）と hessian_ff（低レベル MM）を結合し、酵素活性部位モデルのエネルギー、力、特に**解析的 Hessian**を計算します。

`mlmm_calc.mlmm` は、機械学習原子間ポテンシャル（MLIP）と分子力学力場（Amber prmtop ベースの `hessian_ff`）を組み合わせた減算型 ONIOM の ML/MM calculator を実装しています。`mlmm` のすべての ML/MM 最適化、経路探索、スキャン、振動解析、IRC ワークフローのコア計算機として機能します。

### マルチバックエンドアーキテクチャ

ML（高レベル）コンポーネントは、`-b/--backend` CLI オプションまたは `mlmm.backend` YAML キーで選択される複数の MLIP バックエンドのいずれかによって提供されます:

| バックエンド | 値 | パッケージ | インストール |
| --- | --- | --- | --- |
| FAIR-Chem UMA | `uma`（デフォルト） | `fairchem-core` | `pip install mlmm-toolkit` |
| ORB | `orb` | `orb-models` | `pip install "mlmm-toolkit[orb]"` |
| MACE | `mace` | `mace-torch` | 専用環境: `pip uninstall -y fairchem-core && pip install mace-torch` |
| AIMNet2 | `aimnet2` | `aimnet` | `pip install "mlmm-toolkit[aimnet]"` |

内部的に、すべてのバックエンドは `_MLBackend` 抽象クラスに準拠しており、エネルギー、力、Hessian の評価に対して統一的なインターフェースを提供します。ファクトリ関数が `backend` パラメータに基づいて適切なバックエンドを選択・インスタンス化します。

この計算機は共有結合的 ML/MM 境界にリンク水素原子を自動生成します。ML 領域はモデル PDB（`model.pdb`）で定義され、MM トポロジーは Amber prmtop（`real.parm7`）から取得され、座標は入力 PDB（`input.pdb`）から読み取られます。内部 `real.rst7` は ParmEd により `real.parm7` と `input.pdb` の座標を組み合わせて生成されます -- 外部の `real.rst7` や `real.pdb` は不要です。

## 3 層スキーム（エネルギー / 力 / Hessian）

計算機は ONIOM 減算法を使用して 3 つの評価を組み合わせます:

| レイヤー | システム | 手法 | 説明 |
| --- | --- | --- | --- |
| **REAL-low** | 全系 | MM (hessian_ff) | Amber prmtop ベースの MM で評価した全系 |
| **MODEL-low** | ML サブセット | MM (hessian_ff) | MM で評価した ML 領域 |
| **MODEL-high** | ML サブセット + リンク H | ML (MLIP) | 選択された MLIP バックエンド（デフォルト: UMA）で評価した ML 領域 |

結合エネルギーは:

```
E_ONIOM = E(REAL-low) - E(MODEL-low) + E(MODEL-high)
```

力と Hessian も同じ減算パターンに従います。

## Hessian / 最適化用の層設定

実装では 3 層の B 因子に加え、Hessian 計算対象の MM を指定する別系統の設定を用いて、どの原子が MM Hessian 計算に含まれるか、どの原子が凍結されるかを制御できます:

- **ML 領域**（B 因子 = 0.0）: 選択された MLIP バックエンド（デフォルト: UMA）で処理
- **Movable-MM**（B 因子 = 10.0）: 最適化中に移動する MM 原子
- **Frozen**（B 因子 = 20.0）: 固定された MM 原子
- **Hessian 対象 MM**（専用 B 因子なし）: `hess_cutoff` / `hess_mm_atoms` で選択

レイヤー割り当ては `hess_cutoff`、`movable_cutoff`、`use_bfactor_layers`、および明示的な `*_mm_atoms` リストで制御されます。

## 機能

### リンク原子の再分配

リンク原子からの力と Hessian 寄与はヤコビアンを介して ML/MM 親原子に再分配されます。再分配は以下を追加します:
- セルフ項 `J^T H J`
- ジオメトリ依存の 2 次項 `sum (dJ^T/dx * f_L)` を親原子にインプレースで追加

### MM Hessian

MM バックエンドは `mm_backend` パラメータで選択できます：

- **`"hessian_ff"`**（デフォルトの MM バックエンド）: 解析 Hessian 機能を持つ CPU 専用 MM エンジンです。実効デフォルトは有限差分（`mm_fd: true`）で、`mm_fd: false` を指定すると解析 MM Hessian を使います。アクティブブロックは、必要に応じて凍結行/列をゼロ埋めした完全デカルト形状へ展開できます。
  - CMAP トーション補正（parm7 に含まれる場合は保持）
- **`"openmm"`**: OpenMM による有限差分 (FD) Hessian。CPU と CUDA の両プラットフォームに対応。`hessian_ff` が対応していない力場や、ワークフローで OpenMM を既に使用している場合に有用。

**YAML 設定例:**
```yaml
mlmm:
 mm_backend: openmm # OpenMM を MM 計算に使用
 mm_device: cuda # CUDA を使用 (または "cpu")
```

### 2つの MM 層の CMAP

CMAP（クロスマップ骨格二面角補正）は ff19SB などで使われる 5 原子のトーション補正項です。差し引き式では、REAL と MODEL の MM 計算に同じ CMAP 方針を適用します。

| 領域 | E_MM(real) | E_MM(model) | ONIOM への正味の影響 |
|--------|-----------|------------|--------------------|
| `use_cmap: true`（デフォルト） | parm7 にあれば CMAP あり | parm7 にあれば CMAP あり | MODEL 内で完結する CMAP は相殺され、境界項は低レベル結合として残る |
| `use_cmap: false` | CMAP 除外 | CMAP 除外 | CMAP を除いた明示的な改変力場計算 |

ff19SB では CMAP が、対応するゼロ化済み骨格 cosine 項を置き換えます（[Tian et al., 2020](https://doi.org/10.1021/acs.jctc.9b00591)）。このため CMAP を保持するのが力場に忠実なデフォルトです。`use_cmap: false` は両方の MM 層から CMAP を除去します。

**YAML 設定例:**
```yaml
mlmm:
 use_cmap: false  # 両方の MM 層から CMAP を明示的に除去
```

### ML Hessian モード

- `"Analytical"`: UMA、ORB、MACE、AIMNet2 の自動微分またはネイティブ Hessian 経路。必要な API がない場合は暗黙に計算法を変えずエラーになります。
- `"FiniteDifference"`: 力の中心差分。全 MLIP バックエンドで使用可能です。

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
| 力 | eV/Å | Hartree/Bohr |
| Hessian | eV/Å² | Hartree/Bohr^2 |

PySisyphus インターフェースは原子単位（Hartree/Bohr）に変換された値を返します。

---

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [opt](opt.md) -- ML/MM calculatorを使用した単一構造の構造最適化
- [tsopt](tsopt.md) -- 遷移状態最適化
- [freq](freq.md) -- 振動解析
- [YAML リファレンス](yaml-reference.md) -- `calc`/`mlmm` 設定キー
