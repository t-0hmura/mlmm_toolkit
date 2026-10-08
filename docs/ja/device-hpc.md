# デバイス設定 & HPC セットアップ

ML/MM 計算機の ML と MM の部分をどのデバイス（GPU か CPU）で動かすかの設定と、PBS・Slurm のジョブスクリプトをまとめたページです。デフォルトでは、ML の推論は CUDA が使えれば GPU で、MM の力場（`hessian_ff`）は CPU で動きます。

## デバイスパラメータ

デバイスは YAML ファイル（`--config`）の `calc` セクションで指定します。

| パラメータ | デフォルト | 説明 |
| --- | --- | --- |
| `ml_device` | `auto` | ML 推論のデバイス：`auto`、`cuda`、`cpu`。`auto` は CUDA が使えれば CUDA、なければ CPU を選びます。`--backend dft` では、代わりに `--dft-engine`（`calc.dft.engine`）がデバイスを決めます。 |
| `ml_cuda_idx` | `0` | CUDA で ML 推論を行うときの CUDA デバイス番号。 |
| `mm_backend` | `hessian_ff` | MM のエンジン：`hessian_ff`（CPU のみ）か `openmm`（CPU か CUDA）。 |
| `mm_device` | `cpu` | MM のエンジンのデバイス。`hessian_ff` は `cpu` か `auto` を受け、CPU で動きます。`openmm` は `cuda` も受け、`auto` は OpenMM に CUDA プラットフォームがあれば CUDA を選びます。 |
| `mm_cuda_idx` | `0` | OpenMM を CUDA で動かすときの CUDA デバイス番号。 |
| `mm_threads` | `16` | MM のエンジンの CPU スレッド数。 |

### YAML 設定例

```yaml
calc:
  ml_device: cuda
  ml_cuda_idx: 0
  mm_backend: hessian_ff
  mm_device: cpu
  mm_threads: 16
```

### OpenMM バックエンドを CUDA で使用

```yaml
calc:
  ml_device: cuda
  ml_cuda_idx: 0
  mm_backend: openmm
  mm_device: cuda
  mm_cuda_idx: 0
```

---

## VRAM 管理

### Hessian のデバイス（`--hess-device`）

`freq` と `irc` の `--hess-device` は `cuda`・`cpu`・`auto` を受け、デフォルトの `auto` は `ml_device` に従います。`freq` では、計算済みの Hessian を置いて対角化するデバイスを決めます。`irc` では、初期 Hessian を置くデバイスと IRC の演算を行うデバイスを決めます。

```bash
# デフォルト: ML のデバイス
mlmm freq -i r_complex_layered.pdb --parm7 real.parm7 -q -1

# 計算済み Hessian を CPU へ移して対角化
mlmm freq -i r_complex_layered.pdb --parm7 real.parm7 -q -1 --hess-device cpu
```

`--hess-device cpu` を使う場面：
- Hessian を GPU に置いて対角化すると、計算に要る VRAM を使ってしまう場合

### VRAM 節約のヒント

1. **ML 領域を小さくする:** `mlmm extract` で小さい `--radius` を使います。詳しくは [モデルを削る](model-setup.md#モデルを削る) を参照してください。
2. **hessian_ff（デフォルト）を使う:** hessian_ff は CPU で動くため、GPU に MM の分のメモリを取りません。
3. **VRAM を監視する:** `print_vram` はデフォルトで `true` で、Hessian の計算中に VRAM 使用量のピークを表示します。

---

## ジョブでの精度

精度はバックエンドと用途で選び、割り当てられた GPU で計算時間を測ってください。詳しくは {ref}`MLIP バックエンド › 精度 <ja-precision>` を参照してください。

---

## HPC ジョブ投入

### PBS 例

```bash
#!/bin/bash
#PBS -N mlmm_opt
#PBS -q default
#PBS -l nodes=1:ppn=32:gpus=1,mem=120GB,walltime=72:00:00
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID}
#PBS -e ${PBS_JOBNAME}.e${PBS_JOBID}

set -euo pipefail
hostname
cd "${PBS_O_WORKDIR}"

# hessian_ff は初回利用時に C++ カーネルを JIT ビルドします。システムの
# コンパイラがない、または C++20 モードでコンパイルできない場合は、サイトのコンパイラモジュールを読み込みます:
# module load <COMPILER_MODULE>

# conda 環境の有効化
source ~/miniconda3/etc/profile.d/conda.sh
conda activate <your-env>
command -v g++ >/dev/null || { echo "hessian_ff には g++ が必要です" >&2; exit 1; }
if ! g++ -std=c++20 -x c++ -fsyntax-only /dev/null; then
  echo "hessian_ff には PyTorch の C++20 JIT フラグに対応するコンパイラが必要です" >&2
  exit 1
fi
command -v ninja >/dev/null || { echo "hessian_ff には ninja が必要です" >&2; exit 1; }

# 最適化の実行
mlmm opt \
  -i r_complex_layered.pdb \
  --parm7 real.parm7 \
  -q -1 -m 1 \
  --opt-mode grad \
  --out-dir opt_result
```

### Slurm 例

```bash
#!/bin/bash
#SBATCH --job-name=mlmm_opt
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail
hostname

# hessian_ff は初回利用時に C++ カーネルを JIT ビルドします。必要な場合:
# module load <COMPILER_MODULE>
source ~/miniconda3/etc/profile.d/conda.sh
conda activate <your-env>
command -v g++ >/dev/null || { echo "hessian_ff には g++ が必要です" >&2; exit 1; }
if ! g++ -std=c++20 -x c++ -fsyntax-only /dev/null; then
  echo "hessian_ff には PyTorch の C++20 JIT フラグに対応するコンパイラが必要です" >&2
  exit 1
fi
command -v ninja >/dev/null || { echo "hessian_ff には ninja が必要です" >&2; exit 1; }

mlmm opt \
  -i r_complex_layered.pdb \
  --parm7 real.parm7 \
  -q -1 -m 1 \
  --opt-mode grad \
  --out-dir opt_result
```

### 重要なポイント

- **GPU:** ML 推論用に GPU を 1 基要求します（PBS は `gpus=1`、Slurm は `--gres=gpu:1`）。2 基目を要求するのは、OpenMM の MM バックエンドを別の CUDA デバイス（`mm_device: cuda`、`mm_cuda_idx: 1`）に置く場合だけです。
- **CPU スレッド:** MM バックエンドのスレッド数（`mm_threads`、デフォルト 16）に足りる CPU を要求します。上の例は余裕を見て 32（`ppn=32`、`--cpus-per-task=32`）を要求しています。
- **メモリ:** 代表的な試験計算と scheduler のピークメモリの記録から RAM を決めます。
- **CUDA ランタイム:** 公式 PyTorch wheel には CUDA のユーザー空間ライブラリが含まれるため、通常は対応する NVIDIA ドライバーだけで足ります。CUDA toolkit のモジュールを読み込むのは、それを必要とする拡張を使う場合だけです。

### GPU インデックスの指定

マルチ GPU ノードで特定の GPU を使う場合：

```bash
# 方法 A: 環境変数（全 CUDA プログラムに影響）
# scheduler の下では scheduler が設定した値をそのまま使い、自分で設定するのは scheduler の外だけ
export CUDA_VISIBLE_DEVICES=0

# 方法 B: YAML 設定（mlmm 固有）
# config.yaml に記述:
# calc:
#   ml_cuda_idx: 0
mlmm opt -i r_complex_layered.pdb --parm7 real.parm7 -q -1 --config config.yaml
```

---

## 使用上の注意点

* **hessian_ff は CPU でだけ動く**: デフォルトの `mm_backend: hessian_ff` では、`mm_device` は `cpu` か `auto` だけを受け、`mm_device: cuda` はエラーで止まります。MM を CUDA で動かすには `mm_backend: openmm` を使ってください。
* **ML 推論は GPU 1 基**: デフォルトの `--uma-workers 1` では、ML 推論は `ml_cuda_idx` の GPU 1 基で動きます。`--uma-workers` を 2 以上にする場合は [MLIP バックエンド › ワーカーと Hessian の計算方式](backends.md#ワーカーと-hessian-の計算方式) を参照してください。
* **ML と MM を同じ GPU に置くとメモリを共有する**: 両方を同じ CUDA デバイスで動かすときは、代表的な試験計算でピークメモリを測り、大きな系では `mm_device: cpu` を使ってください。
* **`--hess-device cpu` でもすべての out-of-memory は防げない**: Hessian の計算中にバックエンドの中で起きる out-of-memory は、Hessian を移す前に起きます。モデルを小さくするか、メモリの少ない Hessian やバックエンドの設定を選んでください。
* **各計算ノードに C++ コンパイラが要る**: デフォルトの `hessian_ff` MM バックエンドは、初回利用時に C++ カーネルを JIT ビルドします。これは CUDA とは関係なく行われます。各計算ノードに C++20 対応のコンパイラと Ninja が必要です（GCC 13.3 で検証済み）。システムの `g++` がない、または古い場合はコンパイラモジュールを読み込んでください。上のジョブスクリプトは両方を確かめます。

---

## 関連ドキュメント

- [インストール](installation.md) — インストール、CUDA、C++ コンパイラ
- [ML/MM 計算機](mlmm-calc.md) — 計算機の構成とパラメータ
- [MLIP バックエンド](backends.md) — 精度、ワーカー、Hessian の計算方式
- [YAML 設定の一覧](yaml-reference.md) — 設定の全体
- [freq](freq.md) · [irc](irc.md) — `--hess-device`
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法
