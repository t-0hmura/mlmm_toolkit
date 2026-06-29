# Python API

mlmm-toolkit を Python ライブラリとして使用します — `MLMMCore`（基盤エンジン）、`MLMMASECalculator`（ASE インターフェース）、`mlmm`（pysisyphus Calculator）。

## クイックスタート

```python
from mlmm import MLMMCore, MLMMASECalculator, mlmm

# 基盤エンジン — energy (eV), forces (eV/Å), Hessian (eV/Å²)を返す
core = MLMMCore(
    input_pdb="complex_layered.pdb",
    real_parm7="real.parm7",
    model_pdb="ml_region.pdb",
    model_charge=0,
)

import numpy as np
coords = np.loadtxt(...)  # shape (N, 3), Å
result = core.compute(coords, return_forces=True, return_hessian=False)
print(result["energy"], result["forces"].shape)
```

## APIレベル

mlmm-toolkitは使用状況に応じて3つのAPIレベルを提供します:

| レベル | クラス | 入力単位 | 出力単位 | 用途 |
|--------|--------|----------|----------|------|
| 基盤エンジン | `MLMMCore` | Å | eV, eV/Å, eV/Å² | Pythonスクリプトから直接使用 |
| ASE | `MLMMASECalculator` | Å (`Atoms`経由) | eV, eV/Å | ASEベースのワークフロー（DMF, MD） |
| pysisyphus | `mlmm` (Calculator) | Bohr (`Geometry`経由) | Hartree, Hartree/Bohr | pysisyphusでの構造最適化、IRC、振動解析 |

## MLMMCore

コアML/MMエンジン。トポロジー、力場、MLIPバックエンドを初回に初期化し、以降の`compute()`呼び出しでは座標のみ更新します。

```python
from mlmm import MLMMCore

core = MLMMCore(
    input_pdb="complex_layered.pdb",
    real_parm7="real.parm7",
    model_pdb="ml_region.pdb",
    model_charge=0,
    model_mult=1,
    backend="uma",               # uma | orb | mace | aimnet2
    embedcharge=False,           # xTB点電荷埋め込み補正
    return_partial_hessian=True, # 部分Hessian（ML領域のみ）
)
```

### 主要パラメータ

| パラメータ | 型 | デフォルト | 説明 |
|------------|------|---------|------|
| `input_pdb` | `str` | *必須* | 入力PDB（B-factorレイヤー付き全系） |
| `real_parm7` | `str` | *必須* | 全系のAmber prmtop |
| `model_pdb` | `str` | *必須* | ML領域を定義するPDB |
| `model_charge` | `int` | `0` | ML領域の電荷 |
| `model_mult` | `int` | `1` | スピン多重度 |
| `backend` | `str` | `"uma"` | MLIPバックエンド |
| `embedcharge` | `bool` | `False` | xTB点電荷埋め込みの有効化 |
| `mm_backend` | `str` | `"hessian_ff"` | MMエンジン（`hessian_ff`または`openmm`） |
| `return_partial_hessian` | `bool` | `True` | 部分Hessian（ML + 境界）を返す |
| `link_mlmm` | `list` | `None` | リンク原子の手動指定 |

### compute()

```python
result = core.compute(
    coord_ang,                   # numpy (N, 3), Å
    return_forces=True,
    return_hessian=False,
)
# result["energy"]   : float (eV)
# result["forces"]   : numpy (N, 3) (eV/Å)
# result["hessian"]  : torch 4D (eV/Å²) — return_hessian=Trueの場合のみ。
#   return_partial_hessian=True（デフォルト）では形状 (n_active, 3, n_active, 3)、
#   False では展開された (N, 3, N, 3)。部分Hessianの場合は付随キー
#   "within_partial_hessian"（active_atoms / active_dofs / full_to_active のマッピングを含む dict）が併せて返される。
```

## MLMMASECalculator

`MLMMCore`をラップするASE `Calculator`。ASEの最適化器、MD、DMFと互換。

```python
from mlmm import MLMMCore, MLMMASECalculator
from ase.io import read

core = MLMMCore(
    input_pdb="complex_layered.pdb",
    real_parm7="real.parm7",
    model_pdb="ml_region.pdb",
)
calc = MLMMASECalculator(core)

atoms = read("complex_layered.pdb")
atoms.calc = calc
print(atoms.get_potential_energy())   # eV
print(atoms.get_forces().shape)       # (N, 3), eV/Å
```

## pysisyphus Calculator (`mlmm`)

pysisyphusの構造最適化、IRC、振動解析に使用。

```python
from mlmm import mlmm as MLMMCalc
from pysisyphus.helpers import geom_loader

calc = MLMMCalc(
    input_pdb="complex_layered.pdb",
    real_parm7="real.parm7",
    model_pdb="ml_region.pdb",
    model_charge=0,
)
geom = geom_loader("complex_layered.pdb")
geom.set_calculator(calc)
energy = geom.energy            # Hartree
forces = geom.forces            # Hartree/Bohr (flat)
```

## v0.1.x 互換性

### パラメータエイリアス

以下の v0.1.x パラメータ名は `DeprecationWarning` 付きで受け付けられます:

| v0.1.x 名 | v0.2.x 名 | 備考 |
|-------------|-------------|-------|
| `real_pdb` | `input_pdb` | レガシーエイリアス — `input_pdb` にマップ |
| `real_rst7` | *(削除)* | 警告付きで無視（内部で自動生成） |
| `vib_run` | *(削除)* | 警告付きで無視 |
| `vib_dir` | *(削除)* | 警告付きで無視 |

```python
# v0.1.x スタイル — 動作するが DeprecationWarning を発する
from mlmm import MLMMCore
core = MLMMCore(real_pdb="complex.pdb", real_parm7="real.parm7", model_pdb="ml.pdb")
```

### mlmm_ase() ファクトリ

v0.1.x の `mlmm_ase(real_pdb=..., ...)` ユーティリティ関数は後方互換のため保持されています:

```python
from mlmm import mlmm_ase

# v0.1.x スタイル — DeprecationWarning を発する
calc = mlmm_ase(real_pdb="complex.pdb", real_parm7="real.parm7", model_pdb="ml.pdb")
# 等価: MLMMASECalculator(MLMMCore(input_pdb="complex.pdb", ...))
```

## 関連ドキュメント

- [ML/MM Calculator](mlmm-calc.md) — アーキテクチャと内部詳細
- [YAML Reference](yaml-reference.md) — `--config` YAML の設定キー
