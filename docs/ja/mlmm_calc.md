# ML/MM 計算機

## 概要

> **要約:** PySisyphus 用の ONIOM 型 ML/MM 計算機。MLIP バックエンド（デフォルト: FAIR-Chem UMA、選択肢: `orb`、`mace`、`aimnet2`）と hessian_ff（低レベル MM）を結合し、酵素活性部位モデルのエネルギー、力、特に**解析的ヘシアン**を計算します。

`mlmm_calc.mlmm` は、機械学習原子間ポテンシャル（MLIP）と分子力学力場（Amber prmtop ベースの `hessian_ff`）を組み合わせた減算型 ONIOM スタイルの ML/MM 計算機を実装しています。`mlmm` のすべての ML/MM 最適化、経路探索、スキャン、振動解析、IRC ワークフローのコア計算機として機能します。

### マルチバックエンドアーキテクチャ

ML（高レベル）コンポーネントは、`-b/--backend` CLI オプションまたは `mlmm.backend` YAML キーで選択される複数の MLIP バックエンドのいずれかによって提供されます:

| バックエンド | 値 | パッケージ | インストール |
| --- | --- | --- | --- |
| FAIR-Chem UMA | `uma`（デフォルト） | `fairchem-core` | `pip install mlmm` |
| ORB | `orb` | `orb-models` | `pip install "mlmm-toolkit[orb]"` |
| MACE | `mace` | `mace-torch` | 別環境が必要（README 参照） |
| AIMNet2 | `aimnet2` | `aimnet2` | `pip install "mlmm-toolkit[aimnet2]"` |

内部的に、すべてのバックエンドは `_MLBackend` 抽象クラスに準拠しており、エネルギー、力、ヘシアンの評価に対して統一的なインターフェースを提供します。ファクトリ関数が `backend` パラメータに基づいて適切なバックエンドを選択・インスタンス化します。

### 点電荷埋め込み補正

`--embedcharge` を有効化すると（YAML では `mlmm.embedcharge: true`）、xTB 点電荷埋め込み補正が適用されます:

```
dE = E_xTB(ML + MM_charges) - E_xTB(ML_only)
```

この補正は、MM 環境からの静電的影響を点電荷を介して ML 領域にモデル化し、ML/MM 境界での分極効果の記述を改善します。対応する力とヘシアンの補正も同様に適用されます。デフォルトは `--no-embedcharge`（無効）です。`$PATH` 上に xTB 実行ファイルが必要です。

この計算機は共有結合的 ML/MM 境界にリンク水素原子を自動生成します。ML 領域はモデル PDB（`model.pdb`）で定義され、MM トポロジーは Amber prmtop（`real.parm7`）から取得され、座標は入力 PDB（`input.pdb`）から読み取られます。内部 `real.rst7` は ParmEd により `real.parm7` と `input.pdb` の座標を組み合わせて生成されます -- 外部の `real.rst7` や `real.pdb` は不要です。

## 3 層スキーム（エネルギー / 力 / ヘシアン）

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

力とヘシアンも同じ減算パターンに従います。

## ヘシアン / 最適化用の層設定

実装では 3 層の B 因子と Hessian 対象 MM の別設定を使い、どの原子が MM ヘシアン計算に含まれるか、どの原子が凍結されるかを制御できます:

- **ML 領域**（B 因子 = 0.0）: 選択された MLIP バックエンド（デフォルト: UMA）で処理
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
  - CMAP トーション補正（実装済みだが Gaussian と同様にデフォルトでは無効）
- **`"openmm"`**: OpenMM による有限差分 (FD) ヘシアン。CPU と CUDA の両プラットフォームに対応。`hessian_ff` が対応していない力場や、ワークフローに OpenMM を既に使用している場合に有用。

**YAML 設定例:**
```yaml
mlmm:
 mm_backend: openmm # OpenMM を MM 計算に使用
 mm_device: cuda # CUDA を使用 (または "cpu")
```

### model 系の CMAP

CMAP（クロスマップ骨格二面角補正）は、タンパク質 AMBER 力場で骨格コンフォメーションサンプリングを改善するために使用される 5 原子のトーション補正項です。ONIOM では、model 系の parm7 は実系トポロジーを ML 領域にスライスして生成されます。

デフォルト（`use_cmap: false`）では、CMAP 項を model parm7 から**除外**します:

| 領域 | E_MM(real) | E_MM(model) | ONIOM への正味の影響 |
|--------|-----------|------------|--------------------|
| `use_cmap: false`（デフォルト） | CMAP あり | CMAP **除外** | model 骨格 CMAP が E_total に残留 |
| `use_cmap: true` | CMAP あり | CMAP あり | model 骨格 CMAP が差し引きでキャンセル |

このデフォルト動作は Gaussian ONIOM と一致しており、Gaussian ONIOM も model MM パラメータから CMAP を除外します。典型的な活性部位モデル（リガンド + 機能的残基、ML 領域に骨格原子を含まない）では、どちらの設定でも CMAP in model はゼロになります。

**YAML 設定例:**
```yaml
mlmm:
 use_cmap: true  # model parm7 に CMAP を含める（Gaussian 非互換の挙動）
```

### ML ヘシアンモード

- `"Analytical"`: 選択されたデバイスでの 2 次自動微分。**UMA バックエンドでのみ使用可能**。
- `"FiniteDifference"`: 力の中心差分。**ORB、MACE、AIMNet2** バックエンドで使用（VRAM が限られている場合は UMA でも使用可能）。

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
| ヘシアン | eV/Å² | Hartree/Bohr^2 |

PySisyphus インターフェースは原子単位（Hartree/Bohr）に変換された値を返します。

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [opt](opt.md) -- ML/MM 計算機を使用した単一構造の構造最適化
- [tsopt](tsopt.md) -- 遷移状態最適化
- [freq](freq.md) -- 振動解析
- [YAML リファレンス](yaml_reference.md) -- `calc`/`mlmm` 設定キー
