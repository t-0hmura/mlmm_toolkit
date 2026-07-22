# デバイス設定 & HPC セットアップ

## 概要

ML/MM calculatorの GPU/CPU デバイス設定と、HPC クラスタでのジョブ投入方法を説明します。

### 要点
- **ML バックエンド (UMA):** デフォルトで CUDA を使用（`ml_device: auto` → CUDA が利用可能なら CUDA）。
- **MM バックエンド (hessian_ff):** CPU のみ。OpenMM バックエンドは CUDA を使用可能。
- **Hessian 組み立て:** `--hess-device cpu` で CPU にオフロードし、VRAM を節約可能。
- **マルチ GPU:** ML 推論は単一 GPU（モデル並列は非対応）。OpenMM MM バックエンドは別の CUDA デバイス（`mm_device: cuda`, `mm_cuda_idx: 1`）に配置可能。

---

## デバイスパラメータ

ML/MM calculator（`mlmm_calc.mlmm`）は ML と MM で別々のデバイス設定を使用します:

| パラメータ | デフォルト | 説明 |
| --- | --- | --- |
| `ml_device` | `auto` | UMA 推論のデバイス。`auto` は CUDA が利用可能なら CUDA、なければ CPU。 |
| `ml_cuda_idx` | `0` | `ml_device=cuda` 時の CUDA デバイスインデックス。 |
| `mm_backend` | `hessian_ff` | MM 力場エンジン。`hessian_ff`（解析的、CPU のみ）または `openmm`（CUDA 対応）。 |
| `mm_device` | `cpu` | MM バックエンドのデバイス。hessian_ff は `cpu` 必須。openmm は `cuda` 使用可能。 |
| `mm_cuda_idx` | `0` | `mm_device=cuda` 時の CUDA デバイスインデックス（openmm のみ）。 |
| `mm_threads` | `16` | MM バックエンドの CPU スレッド数。 |

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

> **注意:** ML と MM の両方が CUDA を使用する場合、GPU メモリを共有します。大きな系では `mm_device: cpu` を使用して VRAM 消費を抑えることを推奨します。

---

## VRAM 管理

### Hessian デバイス（`--hess-device`）

`freq` コマンドは `--hess-device` で Hessian の組み立て・対角化のデバイスを制御できます：

```bash
# デフォルト: ml_device と同じ（通常 CUDA）
mlmm freq -i input.pdb --parm real.parm7 -q -1

# CPU で Hessian 組み立て（大きな系で VRAM を節約）
mlmm freq -i input.pdb --parm real.parm7 -q -1 --hess-device cpu
```

`--hess-device cpu` を使用する場面：
- 活性領域が大きい場合（非凍結原子 > 約 500）
- 振動数計算で CUDA out-of-memory エラーが発生する場合
- VRAM が限られている場合（< 16 GB）

### VRAM 節約のヒント

1. **ML 領域を小さくする:** `mlmm extract` で小さい `--radius` を使用、または `mlmm define-layer` で `--radius-freeze` を絞る。
2. **hessian_ff（デフォルト）を使用:** hessian_ff は CPU のみなので、VRAM はすべて UMA に使用可能。
3. **大きな系では OpenMM CUDA を避ける:** ML と MM の両方が CUDA を使うと VRAM 圧力が倍増する。
4. **VRAM を監視:** `print_vram` はデフォルトで true（Hessian 計算中に VRAM 使用量（ピーク）を表示）。抑制するには YAML で `print_vram: False` を設定。

---

## バックエンドごとの精度デフォルト値

`--precision` は `fp32` または `fp64`（大文字小文字無視）を選びます。
未指定時の有効なデフォルト値はバックエンドごとに異なります。

| backend | デフォルト | 理由 |
|---|---|---|
| UMA | fp32 | 上流 fairchem の baseline。 |
| ORB | fp64 | ORB fp32 は縮約 `float32-high`（TF32）matmul を使い、force noise が有限差分 Hessian に偽の虚振動を作る場合がある。 |
| MACE | fp64 | 上流の `default_dtype="float64"` と一致。 |
| AIMNet2 | fp32 | 精度切替を持たず、明示的 fp64 は拒否。 |

ORB/MACE で `--precision fp32` を明示するのは、Hessian の noise より
screening throughput を優先する場合に限ります。UMA fp64 は数値的に敏感な
TS/Hessian を安定化する場合がありますが、consumer GPU では遅くなります。
どの精度でも freq と IRC による独立検証が必要です。

```bash
# データセンター H200 — フル精度のベース推論
mlmm tsopt -i ts.pdb --parm enzyme.parm7 -l 'LIG:Q' -b uma --precision fp64 -o result_ts

# ORB の縮約精度を明示した screening
mlmm scan -i r.pdb --parm enzyme.parm7 -l 'LIG:Q' -b orb --precision fp32 --scan-lists '[(1,5,1.4)]' -o result_scan
```

`--precision` はすべての計算系サブコマンド（`sp`、`opt`、`tsopt`、`freq`、`irc`、`scan` / `scan2d` / `scan3d`、`path-opt`、`path-search`、`all`）で受け付けられ、バックエンドごとにルーティングされます（UMA precision、ORB precision、MACE `default_dtype`）。

```{note}
`-b aimnet2` では `fp32` は no-op、`fp64` は*拒否*されます — モデル入力が上流で float32 にキャストされるためです。fp64 が必要なら `uma`、`orb`、`mace` を使ってください。`--precision fp64` は GPU のリダクション順序ドリフトを*低減*しますが、実行をビット単位で同一には**しません**。ビット単位の厳密性は `--deterministic` のみが与えます — [再現性](reproducibility.md) を参照。
```

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
# コンパイラがない、または GCC 9 未満なら、サイトのコンパイラモジュールを読み込みます:
# module load <COMPILER_MODULE>

# conda 環境の有効化
source ~/miniconda3/etc/profile.d/conda.sh
conda activate <your-env>
command -v g++ >/dev/null || { echo "hessian_ff には g++ が必要です" >&2; exit 1; }
gxx_major=$(g++ -dumpversion | cut -d. -f1)
if [[ ! $gxx_major =~ ^[0-9]+$ ]] || (( gxx_major < 9 )); then
  echo "hessian_ff には GCC >= 9 が必要です（検出したメジャー: $gxx_major）" >&2
  exit 1
fi
command -v ninja >/dev/null || { echo "hessian_ff には ninja が必要です" >&2; exit 1; }

# 最適化の実行
mlmm opt \
  -i r_complex_layered.pdb \
  --parm p_complex.parm7 \
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
gxx_major=$(g++ -dumpversion | cut -d. -f1)
if [[ ! $gxx_major =~ ^[0-9]+$ ]] || (( gxx_major < 9 )); then
  echo "hessian_ff には GCC >= 9 が必要です（検出したメジャー: $gxx_major）" >&2
  exit 1
fi
command -v ninja >/dev/null || { echo "hessian_ff には ninja が必要です" >&2; exit 1; }

mlmm opt \
  -i r_complex_layered.pdb \
  --parm p_complex.parm7 \
  -q -1 -m 1 \
  --opt-mode grad \
  --out-dir opt_result
```

### 重要なポイント

- **GPU 1 基:** mlmm-toolkit はジョブあたり GPU 1 基を使用。PBS なら `gpus=1`、Slurm なら `--gres=gpu:1` を指定。
- **CPU スレッド:** MM バックエンド用に十分な CPU を確保（`mm_threads` デフォルト 16）。PBS なら `ppn=32`、Slurm なら `--cpus-per-task=32` を推奨。
- **メモリ:** 酵素活性部位モデルには通常 120 GB で十分。非常に大きな系では増量。
- **CUDA ランタイム:** 公式 PyTorch wheel には CUDA のユーザー空間ライブラリが含まれるため、通常は互換性のある NVIDIA ドライバーだけで十分です。必要な拡張をソースビルドする場合だけ CUDA toolkit module を読み込みます。
- **C++ コンパイラ:** デフォルトの `hessian_ff` MM バックエンドは初回利用時に C++ カーネルを JIT ビルドします。CUDA とは独立に、各計算ノードで GCC 9 以上と Ninja が必要です。システムの `g++` がない、または古い場合はコンパイラモジュールを読み込みます。

### GPU インデックスの指定

マルチ GPU ノードで特定の GPU を使用する場合：

```bash
# 方法 A: 環境変数（全 CUDA プログラムに影響）
export CUDA_VISIBLE_DEVICES=0

# 方法 B: YAML 設定（mlmm 固有）
# config.yaml に記述:
# calc:
#   ml_cuda_idx: 0
mlmm opt -i input.pdb --parm real.parm7 -q -1 --config config.yaml
```

---

## 制限事項

- **ML モデル並列は非対応:** ML 推論は単一 GPU で動作する。OpenMM MM バックエンドは別の CUDA デバイス（`mm_device: cuda`, `mm_cuda_idx`）を使用可能だが、デフォルトの hessian_ff MM バックエンドは CPU のみ。
- **分散計算非対応:** すべての計算は単一ノードの単一プロセス内で実行。
- **hessian_ff は CPU のみ:** デフォルトの MM バックエンドでは `mm_device` は `cpu`/`auto` のみ可。`mm_device: cuda` を指定すると ValueError を送出（暗黙の CPU フォールバックはしない）。

---

## 関連項目

- [はじめに](getting-started.md) -- インストールと CUDA セットアップ
- [ML/MM calculator](mlmm-calc.md) -- 計算機のアーキテクチャとパラメータ
- [YAML リファレンス](yaml-reference.md) -- 設定リファレンス
- [freq](freq.md) -- `--hess-device` オプションの詳細
- [トラブルシューティング](troubleshooting.md) -- よくあるエラーの修正
