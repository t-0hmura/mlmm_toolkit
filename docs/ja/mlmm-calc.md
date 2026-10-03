# ML/MM 計算機

ML 領域を MLIP で、周りのタンパク質と溶媒を Amber の MM 力場で計算し、2 つを ONIOM の差し引きで合わせる計算機の仕組みと、Python から呼ぶ方法をまとめたページです。

ML/MM の構造最適化、経路探索、スキャン、振動解析、IRC のコマンドは、すべてこの計算機を使います。ML 領域は `-b/--backend` で選んだバックエンドで計算します。選べるのは `uma`（既定）、`orb`、`mace`、`aimnet2`、DFT/MM の `dft` です。バックエンドごとのインストール、モデル名、設定は [MLIP バックエンド](backends.md)にあります。

## ONIOM のエネルギー分解

計算機は、3 つの評価を ONIOM の差し引きで合わせます。

| 評価 | 系 | 手法 | 説明 |
| --- | --- | --- | --- |
| **REAL-low** | 全系 | MM（`hessian_ff` または OpenMM） | Amber parm7 の力場で評価した全系 |
| **MODEL-low** | ML 領域 | MM（同じエンジン） | MM で評価した ML 領域 |
| **MODEL-high** | ML 領域 + リンク H | MLIP または DFT | 選んだバックエンド（既定は UMA）で評価した ML 領域 |

合わせたエネルギーは次の式です。

```
E_ONIOM = E(REAL-low) - E(MODEL-low) + E(MODEL-high)
```

全系を MM で評価し、ML 領域を高レベルと MM の両方で評価して、ML 領域の MM エネルギーを差し引くことで二重に数えないようにします。力と Hessian も同じ差し引きで合わせます。

### 一般的な QM/MM との比較

| 側面 | 一般的な QM/MM | mlmm-toolkit の ML/MM |
| --- | --- | --- |
| 高レベル手法 | DFT、HF、post-HF | MLIP（UMA、ORB、MACE、AIMNet2）または DFT（`-b dft`） |
| 低レベル手法 | OpenMM / Amber | `hessian_ff`（既定）/ OpenMM |
| リンク原子 | 通常は必要 | parm7 の結合のうち ML/MM 境界をまたぐものすべてに自動で付く |
| 埋め込み | 静電埋め込みが一般的 | 既定は機械的埋め込み。`--embedcharge` で MM の点電荷を入れる（MLIP では xTB の補正、DFT では PySCF の Hamiltonian） |
| 速度 | 遅い（QM が律速） | MLIP なら速い（GPU で推論）。`-b dft` は DFT と同じコスト |

## Hessian / 最適化用の層設定

各原子は、入力 PDB の B-factor から読んだ 3 つの層（ML、可動 MM、凍結 MM）のどれかに属します（{ref}`ML 領域と層の組み方 › MM の層 <ja-mm-layers>`）。さらに 2 つの設定で、Hessian に入れる MM 原子と動かす MM 原子を決めます。

- **Hessian 対象 MM**（専用の B-factor はありません）：Hessian の行と列を計算する可動 MM 原子です。既定では可動 MM 原子をすべて含めます。`--hessian-cutoff` を指定すると ML 領域からその距離（Å）以内の可動 MM 原子だけを含めます。
- **距離で決める可動 MM**：`--movable-cutoff` を指定すると、ML 領域からその距離（Å）以内の MM 原子を可動にし、残りを凍結します。B-factor の層の代わりにこの距離を使います。

## 機能

### リンク原子の再分配

ML/MM 境界が共有結合を切るとき、MODEL-high の計算では、切れた結合の ML 側をリンク水素でキャップします。計算機は、parm7 の結合のうち片方の端だけが ML 領域にあるものすべてか、`link_mlmm` で指定した組にリンク水素を付けるので、`model.pdb` にはリンク水素を入れません。配置の方式は `--link-atom-method` で選べます。

| 方式 | 配置 | 推奨 |
| --- | --- | --- |
| **scaled**（g-factor、既定） | `r_L = r_QM + g·(r_MM − r_QM)`、`g = (CR_QM + CR_H)/(CR_QM + CR_MM)`（共有結合半径） | 推奨（滑らかな PES、一定のヤコビアン） |
| **fixed** | `r_L = r_QM + d·û`、`d` = 1.09 Å（親が C）/ 1.01 Å（親が N）、`û` は MM 原子へ向かう単位ベクトル | 推奨しない（座標に依存するヤコビアン） |

scaled は Gaussian ONIOM と同じ Morokuma–Dapprich の g-factor 法で、リンク水素は QM–MM 距離に比例して動きます。fixed は、リンク水素を結合軸に沿って一定の距離に置きます。

リンク水素にかかる力は、リンク位置のヤコビアン `J` を通して ML 側と MM 側の親原子に渡します。

```
F_QM += (1−g) · F_link    （scaled の場合）
F_MM += g · F_link
```

Hessian では、2 つの親原子のブロックに自己項 `Jᵀ H_link J` を足します。fixed では `J` が座標に依存するので、力で重み付けした 2 次微分の項 `Σ (∂Jᵀ/∂x) f_L` も足します。scaled は親原子の座標に対して線形なので、この項は要りません。

(ja-microiteration)=
### マイクロイテレーション

可動 MM 原子が多い系では、全座標を同時に最適化すると、MM 環境だけが緩和している間も毎ステップ MLIP の勾配を計算することになります。マイクロイテレーションは、Gaussian 16 と同じく 2 種類のステップを交互に行います。

```
最初に MM 環境を 1 回緩和し、収束するまで繰り返す:
    マクロステップ — ML 原子 + リンク原子の MM 側の親原子を 1 ステップ動かす（ONIOM の全力）
    マイクロステップ — 残りの可動 MM 原子を L-BFGS で緩和（MM の力のみ）
```

| | マクロステップ | マイクロステップ |
|---|---|---|
| **計算機** | ONIOM の全体（`E_MM(real) + E_ML(model) − E_MM(model)`） | MM 力場のみ（`E_MM(real)`） |
| **動かす座標** | ML 原子 + リンク原子の MM 側の親原子 | リンク原子の MM 側の親原子を除く可動 MM 原子 |
| **オプティマイザ** | `opt`：厳密な Hessian から始めて TS-BFGS で更新する RFO。`tsopt`：選んだ Hessian 型の TS オプティマイザ | L-BFGS（Hessian を使わず、マイクロステップごとに新しく始める） |
| **収束判定** | `--thresh`（既定は `opt` で `gau`、`tsopt` で `baker`） | `microiter.micro_thresh`（既定は `--thresh` と同じ） |

マイクロイテレーションは `--microiter/--no-microiter`（既定はオン）で切り替えます。使われるのは `opt --opt-mode hess` と、`tsopt` の Hessian 型のモード（`hess`、`rsirfo`、`rsprfo`、`trim`）で、`tsopt` の既定は `hess` です。マイクロステップの設定キーは、YAML リファレンスの [`microiter`](yaml-reference.md#microiter) にあります。

```{note}
**リンク原子の MM 側の親原子をマクロステップで動かす理由：**
scaled（g-factor）のリンク原子では、`r_L = (1−g)·r_QM + g·r_MM` によって、リンク原子の位置が QM 側と MM 側の**両方**の親原子に結び付いています。マイクロステップで MM 側の親原子が（ML の寄与のない MM の力だけで）動くと、サイクルの間でリンク原子の位置がずれ、マクロステップのエネルギーが振動します。MM 側の親原子をマクロステップで ML 原子と一緒に動かし、マイクロステップでは止めておくことで、このずれをなくします。
```

### MM Hessian

MM のエンジンは `--mm-backend`（YAML では `calc.mm_backend`）で選べます。

- **`hessian_ff`**（既定）：mlmm-toolkit に同梱された、Amber parm7 の力場用の CPU 専用 MM エンジンです。結合、角度、二面角、不正二面角、Lennard-Jones、静電、CMAP の項を計算し、MM の Hessian を解析的に計算できます。MM の Hessian は既定では有限差分（`calc.mm_fd: true`）で、`calc.mm_fd: false` にすると `hessian_ff` の解析 Hessian を使います。C++ のカーネルは初回の使用時に自動でビルドされ、C++20 のコンパイラが要ります（[インストール](installation.md)）。MM を CPU で計算するので、GPU のメモリは ML 領域に使えます。
- **`openmm`**：CPU または CUDA で動く OpenMM で、Hessian は有限差分です。`hessian_ff` が対応していない力場を使うときや、ワークフローですでに OpenMM を使っているときに選んでください。`mm_backend` と `mm_device` の YAML の例と VRAM の兼ね合いは、[デバイス設定と HPC](device-hpc.md) にあります。

一部の原子だけが動くときは、動く原子の Hessian のブロックを、凍結した原子の行と列を 0 で埋めた全デカルト座標の形に広げることもできます（`return_partial_hessian`）。

### ML Hessian

ML 領域の Hessian の作り方は、`--hessian-calc-mode`（YAML では `calc.hessian_calc_mode`）で選べます。

- `FiniteDifference`（既定）：力の中心差分です。すべてのバックエンドで使えます。
- `Analytical`：バックエンドの自動微分またはネイティブの Hessian です。UMA（ワーカー 1 つのとき）、ORB、MACE、AIMNet2 は入れた版が対応していれば使え、DFT のバックエンドも `--embedcharge` なしなら使えます。バックエンドごとの対応とメモリの兼ね合いは、[MLIP バックエンド › Hessian の計算方式](backends.md#hessian-の計算方式)にあります。

合わせた Hessian は、既定では float64 で組み立てます（`calc.H_double: true`）。

### 2 つの MM 層の CMAP

CMAP（クロスマップ骨格二面角補正）は、ff19SB などで使われる 5 原子のトーション補正項です。差し引きの式では、REAL と MODEL の MM 計算に同じ CMAP の扱いを適用します。

| 領域 | E_MM(real) | E_MM(model) | ONIOM への正味の影響 |
|--------|-----------|------------|--------------------|
| `use_cmap: true`（既定） | parm7 にあれば CMAP あり | parm7 にあれば CMAP あり | MODEL 内で完結する CMAP は相殺され、境界の項は低レベルの結合として残る |
| `use_cmap: false` | CMAP を除く | CMAP を除く | CMAP を除いた、明示的に改変した力場での計算 |

ff19SB では、CMAP が対応する骨格の cosine 項（0 にしてある項）を置き換えます（[Tian et al., 2020](https://doi.org/10.1021/acs.jctc.9b00591)）。このため、CMAP を保つのが力場に忠実な既定です。`use_cmap: false`（CLI では `--no-cmap`）は、両方の MM 層から CMAP を除きます。

**YAML の設定例：**
```yaml
calc:
 use_cmap: false  # 両方の MM 層から CMAP を明示的に除く
```

## 入力

| 入力 | CLI | 説明 |
| --- | --- | --- |
| `input.pdb` | `-i` | 入力構造。残基名、原子名、B-factor の層をここから読みます |
| `real.parm7` | `--parm7` | 全系（REAL）の Amber トポロジー |
| `model.pdb` | `--model-pdb` | ML 領域を定める PDB（ML 原子の特定に使います） |

`input.pdb` の原子の並びは `real.parm7` と一致させてください。計算機は `real.parm7` と `input.pdb` の座標から ParmEd で `real.rst7` を自分で書くので、別に `real.rst7` や `real.pdb` を用意する必要はありません。コマンドラインでは、ML 領域を `--model-indices` や B-factor の ML 層（`--detect-layer`）から決めることもできます。

## 単位

| 量 | 内部単位 | PySisyphus インターフェース |
| --- | --- | --- |
| エネルギー | eV | Hartree |
| 力 | eV/Å | Hartree/Bohr |
| Hessian | eV/Å² | Hartree/Bohr² |

## Python API

この計算機は、CLI を使わずに Python から呼ぶこともできます。`mlmm` パッケージは、`MLMMCore`（基盤エンジン）、`MLMMASECalculator`（ASE インターフェース）、`mlmm`（pysisyphus の Calculator）などのクラスを公開しています。コンストラクタの引数は CLI のオプションに対応します：`-i` → `input_pdb`、`--parm7` → `real_parm7`、`--model-pdb` → `model_pdb`、`-q`/`-m` → `model_charge`/`model_mult`、`-b` → `backend`、`--mm-backend` → `mm_backend`。ほかの引数の多くは、YAML の [`calc` セクション](yaml-reference.md#calc)のキーと同じ名前です。

### クイックスタート

```python
# cd examples/methyltransferase
from ase.io import read
from mlmm import MLMMCore

# 基盤エンジン：energy (eV)、forces (eV/Å)、Hessian (eV/Å²) を返す
core = MLMMCore(
    input_pdb="r_layered.pdb",
    real_parm7="complex.parm7",
    model_pdb="pocket.pdb",
    model_charge=-1,   # この例の ML 領域の電荷
)

coords = read("r_layered.pdb").get_positions()   # shape (N, 3), Å
result = core.compute(coords, return_forces=True, return_hessian=False)
print(result["energy"], result["forces"].shape)
```

### API レベル

mlmm-toolkit は 3 つの API レベルを提供します。

| レベル | クラス | 入力単位 | 出力単位 | 用途 |
|--------|--------|----------|----------|------|
| 基盤エンジン | `MLMMCore` | Å | eV, eV/Å, eV/Å² | Python スクリプトから直接使う |
| ASE | `MLMMASECalculator` | Å（`Atoms` 経由） | eV, eV/Å | ASE ベースのワークフロー（DMF、MD） |
| pysisyphus | `mlmm`（Calculator） | Bohr（`Geometry` 経由） | Hartree, Hartree/Bohr | pysisyphus での構造最適化、IRC、振動解析 |

### MLMMCore

ML/MM の基盤エンジンです。トポロジー、力場、MLIP バックエンドを最初に 1 回だけ用意し、以降の `compute()` では座標だけを更新します。

```python
from mlmm import MLMMCore

core = MLMMCore(
    input_pdb="r_layered.pdb",
    real_parm7="complex.parm7",
    model_pdb="pocket.pdb",
    model_charge=-1,
    model_mult=1,
    backend="uma",               # uma | orb | mace | aimnet2 | dft（dft には dft_settings も必要）
    return_partial_hessian=True, # Hessian 対象の原子の部分 Hessian
)
```

#### 主なパラメータ

| パラメータ | 型 | 既定 | 説明 |
|------------|------|---------|------|
| `input_pdb` | `str` | *必須* | 入力 PDB（B-factor の層を付けた全系） |
| `real_parm7` | `str` | *必須* | 全系の Amber prmtop |
| `model_pdb` | `str` | *必須* | ML 領域を定める PDB |
| `model_charge` | `int` | `0` | ML 領域の正味の電荷。コンストラクタが INFO レベルでログに出すので、電荷を持つ系では期待どおりか確かめてください。 |
| `model_mult` | `int` | `1` | ML 領域のスピン多重度 |
| `backend` | `str` | `"uma"` | MLIP バックエンド |
| `mm_backend` | `str` | `"hessian_ff"` | MM エンジン（`hessian_ff` または `openmm`） |
| `return_partial_hessian` | `bool` | `True` | `True` なら、`compute()` は 4 次元の `(n_active, 3, n_active, 3)` の部分 Hessian と、メタデータの dict `within_partial_hessian`（動く原子の番号、自由度の対応）を返します。`False` なら、全系に広げた 4 次元の `(N, 3, N, 3)` の Hessian を返します。 |
| `link_mlmm` | `list[tuple[str, str]]` | `None` | リンク原子の組を手で指定します。`[("RESN RESID ATOMNAME", "RESN RESID ATOMNAME"), ...]` の形で、1 つ目が ML 側、2 つ目が MM 側です（例：`[("SAM 359 CA", "SAM 359 N")]`）。`None` なら、parm7 のトポロジーで ML/MM 境界をまたぐ結合すべてにリンクを付けます（距離では決めません）。 |

#### compute()

```python
result = core.compute(
    coord_ang,                   # numpy (N, 3), Å
    return_forces=True,
    return_hessian=False,
)
# result["energy"]   : float (eV)
# result["forces"]   : numpy (N, 3) (eV/Å)
# result["hessian"]  : torch 4D (eV/Å²)。return_hessian=True のときだけ。
#   return_partial_hessian=True（既定）では形状 (n_active, 3, n_active, 3)、
#   False では全系に広げた (N, 3, N, 3)。部分 Hessian では、付随するキー
#   "within_partial_hessian"（active_atoms / active_dofs / full_to_active の対応を含む dict）も返す。
```

### MLMMASECalculator

`MLMMCore` を包む ASE の `Calculator` で、エネルギーと力を返します。ASE のオプティマイザ、MD、DMF で使えます。

```python
from mlmm import MLMMCore, MLMMASECalculator
from ase.io import read

core = MLMMCore(
    input_pdb="r_layered.pdb",
    real_parm7="complex.parm7",
    model_pdb="pocket.pdb",
    model_charge=-1,
)
calc = MLMMASECalculator(core)

atoms = read("r_layered.pdb")
atoms.calc = calc
print(atoms.get_potential_energy())   # eV
print(atoms.get_forces().shape)       # (N, 3), eV/Å
```

### pysisyphus Calculator (`mlmm`)

pysisyphus の構造最適化、IRC、振動解析で使います。引数は `MLMMCore` と同じです。

```python
from mlmm import mlmm as MLMMCalc
from pysisyphus.io import geom_from_pdb

calc = MLMMCalc(
    input_pdb="r_layered.pdb",
    real_parm7="complex.parm7",
    model_pdb="pocket.pdb",
    model_charge=-1,
)
geom = geom_from_pdb("r_layered.pdb")
geom.set_calculator(calc)
energy = geom.energy            # Hartree
forces = geom.forces            # Hartree/Bohr (flat)
```

## 使用上の注意点

- `--hessian-calc-mode Analytical` は、選んだバックエンド（またはその入れた版）に解析 Hessian がないとエラーで止まります。計算機が自分で有限差分に切り替えることはありません。
- `use_cmap: false`（`--no-cmap`）は ff19SB に沿った設定ではありません。両方の MM 層から CMAP を除きます。
- `--embedcharge` を付けると、マイクロイテレーションは使われません。MM だけのマイクロステップには埋め込みの力が入らないためで、このときオプティマイザは ML と MM の原子を一緒に動かします。
- fixed のリンク原子の配置は、ML 側の親原子が C か N のときだけ使えます。
- 既定の `return_partial_hessian=True` では、`compute()` は Hessian 対象の原子だけの Hessian を、4 次元の配列 `(n_active, 3, n_active, 3)` で返します。全系に戻すには `within_partial_hessian` を使ってください。
- `backend="dft"` には `dft_settings` も要ります。コマンドラインでは、`-b dft` が YAML の `calc.dft` ブロックと DFT のオプションからこれを作ります（[MLIP の TS を DFT で確かめる](dft-backend.md)）。

## 関連ドキュメント

- [トラブルシューティング](troubleshooting.md)：詳しい対処の手引き
- [opt](opt.md)：ML/MM 計算機を使う単一構造の構造最適化
- [tsopt](tsopt.md)：遷移状態の最適化
- [freq](freq.md)：振動解析
- [YAML リファレンス](yaml-reference.md)：`calc` と `microiter` の設定キー
- [MLIP バックエンド](backends.md)：バックエンドの選び方、インストール、精度、バックエンドの追加
- [デバイス設定と HPC](device-hpc.md)：ML/MM のデバイス設定と HPC での投入
