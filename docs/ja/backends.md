# MLIP バックエンド

ML 領域を計算するバックエンドの選び方と、バックエンドごとのインストール、モデル名、精度、再現性の設定、Hessian の計算方式をまとめたページです。デフォルトのバックエンドは **UMA**（Meta の Universal Models for Atoms）で、`-b/--backend` で **ORB**、**MACE**、**AIMNet2** も選べます。4 つとも機械学習原子間ポテンシャル（MLIP）です。どのバックエンドを選んでも、MM 領域は Amber のトポロジー（`--parm7`）で計算し、2 つを ONIOM で合わせます。

## バックエンドごとの特性

バックエンドは、ML/MM の計算を行うどのコマンドでも `-b/--backend` で選ぶか、YAML の `calc.backend` で設定します。

```bash
# UMA（デフォルト）
mlmm opt -i complex.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0

# ORB
mlmm opt -i complex.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -b orb

# MACE
mlmm opt -i complex.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -b mace

# AIMNet2
mlmm opt -i complex.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -b aimnet2
```

| バックエンド | インストール | モデル名 | 精度の設定 |
|---------|---------|------------------|------------------|
| `uma` | `pip install mlmm-toolkit`（`fairchem-core` は本体の依存）＋ [Hugging Face へのログイン](installation.md) | `uma-s-1p2`（デフォルト）/ `uma-m-1p1` | `uma_precision="fp32" \| "fp64"` |
| `orb` | `pip install "mlmm-toolkit[orb]"` | `orb_v3_conservative_omol` | `orb_precision="float32-high" \| "float32-highest" \| "float64"`（`fp32`・`float32` の名前でも受け付けます） |
| `mace` | 専用の環境で `pip uninstall -y fairchem-core && pip install mace-torch`（`mace-torch` が固定する `e3nn` の版が UMA とぶつかるため、この環境では UMA は動きません） | `MACE-OMOL-0` | `mace_dtype="float32" \| "float64"` |
| `aimnet2` | `pip install "mlmm-toolkit[aimnet]"` | `aimnet2` | なし |

`--backend-model NAME` は、選んだ `--backend` のモデルを替えます（例：`--backend uma --backend-model uma-m-1p1`）。`-b dft` を付けると ML 領域を [DFT](#dftmm-バックエンド) で計算でき、`--calc-file` で {ref}`任意の ASE calculator <ja-backends-custom-calculator>` を使えます。

実行時には、読み込むバックエンドとモデルが `[backend] Preparing MLIP model (UMA / UMA-S-1.2 (OMol))...` のように表示され、[JSON の出力](json-output.md#共通エンベロープ)の `mlip_backend`・`mlip_model`・`mlip_precision` に記録されます。

(ja-precision)=
### 精度（precision）

`--precision fp32|fp64` は MLIP の推論の浮動小数点精度を決めます。`--precision` を指定しないときは、バックエンドごとのデフォルト値を使います。

| バックエンド | `--precision` なし | `--precision fp64` |
|---------|------|------|
| `uma` | fp32 | 使えます |
| `orb` | fp64 | 使えます |
| `mace` | fp64 | 使えます |
| `aimnet2` | fp32（精度の設定がなく、`--precision fp32` を指定しても何も変わりません） | エラー |
| `custom`（`--calc-file`） | calculator 自身の設定 | エラー（`--precision fp32` もエラー） |

どの値を選ぶかは目的で決めます。

| 目的 | 推奨 | 理由 |
| --- | --- | --- |
| 通常の計算 | 指定しない | 上のデフォルト値（UMA・AIMNet2 は fp32、ORB・MACE は fp64）のままにします。 |
| 速さを優先するスクリーニング | 必要なときだけ `--precision fp32` | ORB・MACE の精度が下がります（[使用上の注意点](#使用上の注意点)）。 |
| 最終の TS と Hessian | 指定しない。UMA で n_imag ≥ 2 のときは `--precision fp64` と比べる（{ref}`tsopt <ja-wrong-imaginary-mode-count>`） | 精度によらず、`tsopt` の最後の Hessian で n_imag を確かめ、IRC と端点の最適化で TS が狙った R と P をつなぐことを確かめます。 |

fp64 は次のように指定します。

```bash
mlmm tsopt -i ts.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 --precision fp64
mlmm freq  -i opt.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 --precision fp64
mlmm irc   -i ts.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 --precision fp64
```

YAML では次のように書きます。

```yaml
calc:
 precision: fp64
```

## 決定論的実行と再現性

`--deterministic` を付けると、同じソフトウェアと GPU で、同じ入力から同じ結果が得られます。付けないと、同じ入力の 2 回の GPU の計算でも最後の桁が違うことがあります。

```bash
mlmm opt -i complex.pdb --parm7 enzyme.parm7 --model-pdb ml_region.pdb -q 0 --deterministic
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -q -1 --deterministic
```

| ML バックエンド | `--deterministic` |
|---|---|
| `uma` | 対応。入れたモデルの版で 2 回の計算が一致するか確かめてください |
| `orb` / `mace` | PyTorch の決定論的モードは有効になります。入れた版で 2 回の計算が一致するか確かめてください |
| `aimnet2` | **非対応**：エラーで止まります（[使用上の注意点](#使用上の注意点)） |
| `custom`（`--calc-file`） | **非対応**：エラーで止まります。渡された calculator は mlmm-toolkit の制御の外にあるためです |

## ワーカーと Hessian の計算方式

`--uma-workers N`（デフォルト 1）は UMA の予測器を N 個並列に動かし（`fairchem-core[extras]` が要ります）、`--uma-workers-per-node`（デフォルト 1）はそのうち 1 ノードで動かす数を決めます。どちらのフラグも `opt`、`tsopt`、`freq`、`irc`、`sp`、`all`、`path-opt`、`path-search`、`scan`、`scan2d`、`scan3d` にあります。ほかのバックエンドはこれらを警告を出して無視します。HPC のジョブのテンプレートは [デバイス設定 & HPC セットアップ](device-hpc.md) にあります。

### Hessian の計算方式

`--hessian-calc-mode` で、ML 領域の Hessian の計算方式を選びます。このフラグは `freq`、`irc`、`tsopt`、`sp`、`all` にあり、YAML では `calc.hessian_calc_mode` です。`FiniteDifference`（デフォルト）は力の中心差分を取り、`Analytical` はバックエンドの自動微分またはネイティブの Hessian を使います。UMA（ワーカー 1 つのとき）、ORB、MACE、AIMNet2 は入れた版が対応していれば解析 Hessian を計算でき、DFT のバックエンドも `--embedcharge` なしなら計算できます。自作の calculator は `FiniteDifference` だけに対応します。選んだバックエンドで解析 Hessian が使えないときは、エラーで止まります。[Hessian の MM の部分](mlmm-calc.md)は別に設定します。

UMA では次の 2 つのどちらかを選んでください。

```bash
--uma-workers 1 --hessian-calc-mode Analytical       # 解析 Hessian
--uma-workers 4 --hessian-calc-mode FiniteDifference # 並列の UMA 予測器 + 有限差分
```

モデルの精度と Hessian の精度は別の設定です。Hessian はデフォルト（`calc.H_double: true`）では float64 で組み立て、`H_double: false` にすると float32 で返します。`--precision fp64` のときは Hessian も常に float64 になり、設定ファイルの `H_double: false` は警告を出して上書きされます。

## xTB 静電補正

MLIP/MM の計算では、`--embedcharge/--no-embedcharge`（デフォルトはオフ）で `E_xTB(ML + MM point charges) - E_xTB(ML)` と、それに対応する力と Hessian の差を足し、ML 領域に MM の点電荷の影響を入れます。

`-b dft` では、`--embedcharge` は xTB の補正を使わず、MM の点電荷を PySCF のハミルトニアンに直接入れます。

## DFT/MM バックエンド

ML/MM の計算を行う 11 のコマンド（`all`、`opt`、`tsopt`、`irc`、`freq`、`scan`、`scan2d`、`scan3d`、`path-opt`、`path-search`、`sp`）は `-b dft --func-basis FUNCTIONAL/BASIS --dft-engine gpu|cpu`（デフォルトは `wb97m-v/def2-svp` と `gpu`）を受け付け、MM 領域は Amber の力場のまま、ML 領域を PySCF/GPU4PySCF で計算します。ポピュレーション解析付きの一点計算には、別の `mlmm dft` コマンドがあります。低メモリモードは [`dft`](dft.md#処理の仕組みと計算仕様) に、CPU のスレッド数とホスト RAM、SCF のチェックポイントは [`all` のオプションの一覧（英語のみ）](../reference/commands/all.md) にあります。

(ja-backends-custom-calculator)=
## カスタムバックエンド — 任意の ASE Calculator を使う（`--calc-file`）

組み込みの MLIP バックエンドのほかに、`--calc-file` で渡した任意の [ASE](https://wiki.fysik.dtu.dk/ase/) Calculator で **ML 領域**を計算できます。mlmm-toolkit 本体を変える必要はありません。ML/MM の ONIOM の ML 側に、GFN-xTB（`tblite` か `xtb-python` 経由）、DFTB+、ORCA、Psi4 など、ASE に対応した計算エンジンをつなげます。受け渡しは標準の ASE Calculator の形（エネルギーは eV、力は eV/Å）です。

ASE Calculator を返す `get_calculator` 関数を持つ Python ファイルを書きます。

```python
# my_calc.py（最小の例）
from ase.calculators.emt import EMT

def get_calculator(charge=0, spin=1, device="auto", **kwargs):
    return EMT()
```

`EMT()` を使いたいエンジンに替えてください（GFN-xTB なら `tblite.ase.TBLite(...)`、DFTB+ の ASE calculator、`ase.calculators.orca.ORCA(...)` など）。このファイルを各コマンドか `all` に渡すと `custom` の ML バックエンドが選ばれ、`--backend` の指定より優先されます。

```bash
mlmm sp    -i complex.pdb --parm7 system.parm7 --model-pdb ml_region.pdb --calc-file my_calc.py -q 0 -m 1
mlmm opt   -i complex.pdb --parm7 system.parm7 --model-pdb ml_region.pdb --calc-file my_calc.py -q 0 -m 1
mlmm freq  -i complex.pdb --parm7 system.parm7 --model-pdb ml_region.pdb --calc-file my_calc.py -q 0 -m 1
mlmm all   -i R.pdb P.pdb --parm7 system.parm7 --model-pdb ml_region.pdb --calc-file my_calc.py -q 0 -m 1
```

- 関数が `charge`・`spin`（多重度。`mult`・`multiplicity` の名前でも渡します）・`device` を引数に持つか、`**kwargs` を持てば、これらが渡されるので、全電荷が要るエンジン（xTB など）も設定できます。関数の名前は `--calc-file-func-name NAME` で変えられ、その名前に Calculator のインスタンスを置いてもかまいません。
- 自作の calculator が計算するのは **ML 領域だけ**です。MM 側はふつうどおり `hessian_ff` か OpenMM で計算し、ONIOM の結合も変わりません。Hessian は有限差分で求めるので、`freq` と `tsopt --opt-mode hess` はどのエンジンでも動きます。固定原子もふつうどおり効きます。
- `all` と、ML/MM の計算を行う各サブコマンドで使えます。`all` は、calculator を使うすべての段に同じ関数を渡します。独自の `--backend` 名を持つ、インストールできるバックエンドにするときは [開発者向け](#開発者向け) を見てください。

## Python API

Python では、`MLMMCore` か pysisyphus 用の calculator `mlmm` の `backend` 引数でバックエンドを選びます。クラスと引数、動かせる例は [ML/MM 計算機 › Python API](mlmm-calc.md#python-api) にあります。

## 開発者向け

### バックエンドディスパッチャのパターン

`MLMMCore` は、ML 領域の計算を選んだバックエンドのアダプタに渡し、MM の計算と ONIOM の結合は自分で行います。知らないバックエンド名は `ValueError` になります。mlmm-toolkit に `auto` のバックエンドは無く、各ワークフローはコマンドラインで選んだバックエンドを渡します。

### ファイルマップ

| ファイル | 役割 |
|------|------|
| `mlmm/backends/__init__.py` | `--precision`、`--backend-model`、`--calc-file`、`--uma-workers` を、選んだバックエンドの calculator の設定に変えます |
| `mlmm/backends/mlmm_calc.py` | `MLMMCore`（ML/MM の ONIOM の結合）、`MLMMASECalculator`（ASE）、`mlmm`（pysisyphus の Calculator）、バックエンドごとのアダプタ、有限差分 Hessian の組み立て、単位の変換 |
| `mlmm/backends/pyscf_dft.py` | ML 領域を計算する PySCF/GPU4PySCF の DFT バックエンド。静電埋め込みに対応し、計算のステップの間で SCF の状態を引き継ぎます |

組み込みのバックエンドを独自の `--backend` 名で足すときは、[CONTRIBUTING](https://github.com/t-0hmura/mlmm_toolkit/blob/main/CONTRIBUTING.md) のレシピ 3.2「Add an MLIP backend」に従ってください。

### ML/MM の段での GPU メモリ

ML/MM の段では、選んだ ML バックエンドと Hessian の途中の配列が GPU のメモリを使い、トポロジーの処理と解析的な MM の力場は CPU で動きます。単独の DFT は別の段です。`mlmm/backends/mlmm_calc.py` の有限差分 Hessian のループは、変位の方向を 1 つずつ評価して、同時に行う評価の数を抑えます。まとめて評価する実装に変えるときは、GPU の smoke テストを流し直して VRAM のピークを確かめてください。各段の実行のあとで calculator を解放するので、後の段が前の段のモデルをメモリに残すことはありません。

### ONIOM の結合と MLIP 単体の違い

`mlmm/backends/mlmm_calc.py` の MLIP のアダプタは **ML 領域だけ**を計算します。減算型の ONIOM のエネルギーの式（`# CHEMISTRY-RULE:1`）、リンク原子の Hessian の射影（`# CHEMISTRY-RULE:2`）、3 層の部分 Hessian の組み立て（`# CHEMISTRY-RULE:8`）は同じファイルにあります。新しい MLIP のバックエンドは ONIOM の結合を知らなくてよく、ML 領域のエネルギー、力、Hessian を正しい単位で返せば足ります。

## 使用上の注意点

- ORB と MACE の `--precision fp32` はスクリーニング専用です。結果を使う前に n_imag を確かめてください。
- ORB で `--precision fp32` を指定すると、精度を下げた `float32-high` のモードになります。
- AIMNet2 は `--precision fp64` にも `--deterministic` にも対応せず、どちらもエラーで止まります。AIMNet2 はモデルへの入力を float32 にし、力を PyTorch の決定論的モードの外にある独自の CUDA のコードで計算します。繰り返しの計算を一致させたいときは、UMA、ORB、MACE のいずれかで `--deterministic` を付け、同じ環境で 2 回実行して比べてください。
- `--calc-file` は `--precision`（`fp32` も `fp64` も）と `--deterministic` のどちらも受け付けず、エラーで止まります。精度は自作の calculator の中で設定してください。
- `--deterministic` は PyTorch の決定論的アルゴリズム（`torch.use_deterministic_algorithms`）を有効にし、GPU で決定論的に動く版の無い PyTorch の演算 1 つを置き換えます。
- `--deterministic` はプロセス全体に効きます。`all` に付ければ `all` が実行する段すべてに効くので、段ごとに付ける必要はありません。
- `--deterministic` を付けると遅くなることがあります。決定論的な GPU の演算は別の遅い実装を使うことがあるので、繰り返しの計算を一致させる必要があるときだけ使ってください。
- `--deterministic` では、実行する演算に PyTorch の決定論的な版が無いときは、再現しない結果を黙って出さずに、エラーで止まります。
- 環境変数 `MLMM_STRICT_DETERMINISTIC=1` でも、CI のジョブや Python API で同じモードになります。この変数があると、`--no-deterministic` を付けてもモードは切れません。
- `--deterministic` だけでは、別の計算機やソフトウェアの版でのビット単位の一致は保証されません。使う環境で 2 回実行して比べてください。
- UMA で `--uma-workers` を 2 以上にすると、`--hessian-calc-mode Analytical` とは併用できず、エラーで止まります。並列の予測器は autograd のモデルを持たないためです。解析 Hessian には `--uma-workers 1` を、複数のワーカーには `FiniteDifference` を使ってください。
- MACE は UMA と同じ環境には入れられません。専用の conda 環境に入れてください。
- パッケージを入れていないバックエンドを選ぶと、``orb-models is required for the ORB backend. Install with `pip install orb-models`.`` のようなエラーで止まります。
- `--embedcharge` はエネルギー、力、Hessian を求めるたびに、MM の点電荷あり・なしで xTB を実行します。ML 領域は 200〜300 原子くらいまでにし、実際の系で先に計算時間を測ってください。

## 関連ドキュメント

- [ML/MM 計算機](mlmm-calc.md)：ONIOM の結合、MM の Hessian、Python API（`MLMMCore`、`MLMMASECalculator`、`mlmm`）
- [アーキテクチャ](architecture.md)：ディレクトリの構成と依存の向き
- [デバイス設定 & HPC セットアップ](device-hpc.md)：GPU と CPU の割り当てとジョブのテンプレート
- [MLIP の TS を DFT で確かめる](dft-backend.md)：DFT の設定、メモリ、チェックポイント
- [トラブルシューティング](troubleshooting.md)：計算が失敗したとき
