# デバイス設定 & HPC セットアップ

## 概要

> **要約:** ML/MM カルキュレータの GPU/CPU デバイス設定と HPC クラスタでのジョブ投入方法。

### 要点
- **ML バックエンド (UMA):** デフォルトで CUDA を使用（`ml_device: auto` → CUDA が利用可能なら CUDA）。
- **MM バックエンド (hessian_ff):** CPU のみ。OpenMM バックエンドは CUDA を使用可能。
- **Hessian 組み立て:** `--hess-device cpu` で CPU にオフロードし、VRAM を節約可能。
- **マルチ GPU:** 非対応（ジョブあたり GPU 1 基）。

---

## デバイスパラメータ

ML/MM カルキュレータ（`mlmm_calc.mlmm`）は ML と MM で別々のデバイス設定を使用します：

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
mlmm:
  ml_device: cuda
  ml_cuda_idx: 0
  mm_backend: hessian_ff
  mm_device: cpu
  mm_threads: 16
```

### OpenMM バックエンドを CUDA で使用

```yaml
mlmm:
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
4. **VRAM を監視:** YAML で `print_vram: True` を設定すると、各 ML 評価後に VRAM 使用量を表示。

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

# 環境モジュールのロード
. /home/apps/Modules/init/profile.sh
module load cuda/12.9

# conda 環境の有効化
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mlmm

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

module load cuda/12.9

source ~/miniconda3/etc/profile.d/conda.sh
conda activate mlmm

mlmm opt \
  -i r_complex_layered.pdb \
  --parm p_complex.parm7 \
  -q -1 -m 1 \
  --opt-mode grad \
  --out-dir opt_result
```

### 重要なポイント

- **GPU 1 基:** mlmm_toolkit はジョブあたり GPU 1 基を使用。PBS なら `gpus=1`、Slurm なら `--gres=gpu:1` を指定。
- **CPU スレッド:** MM バックエンド用に十分な CPU を確保（`mm_threads` デフォルト 16）。PBS なら `ppn=32`、Slurm なら `--cpus-per-task=32` を推奨。
- **メモリ:** 酵素活性部位モデルには通常 120 GB で十分。非常に大きな系では増量。
- **CUDA モジュール:** PyTorch が正しい CUDA ランタイムを検出するよう、conda 有効化**前に** CUDA をロード。

### GPU インデックスの指定

マルチ GPU ノードで特定の GPU を使用する場合：

```bash
# 方法 A: 環境変数（全 CUDA プログラムに影響）
export CUDA_VISIBLE_DEVICES=0

# 方法 B: YAML 設定（mlmm 固有）
# config.yaml に記述:
# mlmm:
#   ml_cuda_idx: 0
mlmm opt -i input.pdb --parm real.parm7 -q -1 --config config.yaml
```

---

## 制限事項

- **マルチ GPU 非対応:** pdb2reaction（Ray 経由のマルチワーカー推論をサポート）とは異なり、mlmm_toolkit は GPU 1 基で動作。hessian_ff バックエンドは並列プレディクタに非対応。
- **分散計算非対応:** すべての計算は単一ノードの単一プロセス内で実行。
- **hessian_ff は CPU のみ:** デフォルトの MM バックエンドは `mm_device` の設定に関係なく常に CPU で実行。

---

## 関連項目

- [Getting Started](getting_started.md) -- インストールと CUDA セットアップ
- [ML/MM Calculator](mlmm_calc.md) -- カルキュレータのアーキテクチャとパラメータ
- [YAML Reference](yaml_reference.md) -- 設定リファレンス
- [freq](freq.md) -- `--hess-device` オプションの詳細
- [Troubleshooting](troubleshooting.md) -- よくあるエラーの修正
