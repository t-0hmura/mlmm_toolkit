# Python API

> **概要:** mlmm-toolkitをPythonライブラリとして使用 — `MLMMCore`（基盤エンジン）、`MLMMASECalculator`（ASEインターフェース）、`mlmm`（pysisyphus Calculator）、`mlmm_ase()`互換ラッパー。

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
# result["hessian"]  : torch (3N, 3N) (eV/Å²)  — return_hessian=Trueの場合のみ
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

`mlmm` Calculatorはpysisyphusの`CALC_DICT`にも登録されているため、YAMLワークフローで使用可能です:

```yaml
calc:
  type: mlmm
  input_pdb: complex_layered.pdb
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
  model_charge: 0
```

YAMLワークフローの実行方法は[mlmm pysis](pysis.md)を参照してください。

## v0.1.x互換性

### パラメータエイリアス

以下のv0.1.xパラメータ名は`DeprecationWarning`付きで受け付けます:

| v0.1.x名 | v0.2.x名 | 備考 |
|-----------|-----------|------|
| `real_pdb` | `input_pdb` | エイリアス — `input_pdb`にマッピング |
| `real_rst7` | *(削除)* | 警告付きで無視（内部で自動生成） |
| `vib_run` | *(削除)* | 警告付きで無視 |
| `vib_dir` | *(削除)* | 警告付きで無視 |

```python
# v0.1.xスタイル — 動作しますがDeprecationWarningが出ます
from mlmm import MLMMCore
core = MLMMCore(real_pdb="complex.pdb", real_parm7="real.parm7", model_pdb="ml.pdb")
```

### mlmm_ase()ファクトリ

v0.1.xの`mlmm_ase(real_pdb=..., ...)`便利関数も保持しています:

```python
from mlmm import mlmm_ase

# v0.1.xスタイル — DeprecationWarningが出ます
calc = mlmm_ase(real_pdb="complex.pdb", real_parm7="real.parm7", model_pdb="ml.pdb")
# 等価: MLMMASECalculator(MLMMCore(input_pdb="complex.pdb", ...))
```

### importパス

v0.2.xはv0.1.xと同じimportパスを使用します:

```python
import mlmm
from mlmm import MLMMCore
```

## 関連ドキュメント

- [ML/MM Calculator](mlmm-calc.md) — アーキテクチャと内部詳細
- [mlmm pysis](pysis.md) — pysisyphus YAMLワークフローの実行
- [YAML Reference](yaml-reference.md) — YAMLモードの設定キー
