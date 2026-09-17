# mlmm-toolkit ドキュメント

*バージョン: v{{ release }}*

**mlmm-toolkit** は、機械学習原子間ポテンシャル（Machine Learning Interatomic Potential）と分子力学（Molecular Mechanics）を ONIOM 的に結合した **ML/MM 法** を用いて、PDB 構造から酵素反応経路を自動モデリングする Python 製 CLI ツールキットです。

<img src="../mlmm_toolkit_overview.png" alt="mlmm-toolkit workflow overview" width="90%">

```{toctree}
:maxdepth: 2
:caption: ガイド
:hidden:

getting-started
cif
concepts
quickstart-all
quickstart-scan-spec
quickstart-tsopt-freq
recipes-common-errors
troubleshooting
cli-conventions
reproducibility
```

```{toctree}
:maxdepth: 2
:caption: コマンド
:hidden:

all
extract
add-elem-info
mm-parm
define-layer
opt
tsopt
path-opt
path-search
scan
scan2d
scan3d
freq
irc
dft
sp
trj2fig
oniom-export
oniom-import
fix-altloc
energy-diagram
bond-summary
oniom-gaussian
oniom-orca
```

```{toctree}
:maxdepth: 2
:caption: リファレンス
:hidden:

yaml-reference
json-output
mlmm-calc
python-api
backends
device-hpc
architecture
output-layout
mcp_server
glossary
```

---

## 目的別クイックスタート

インストールと入力の準備は[はじめに](getting-started.md)を参照してください。

| 目的 | ガイド |
|---|---|
| 一気通貫の初回実行 | [クイックスタート: all](quickstart-all.md) |
| 単一構造と結合スキャンから開始 | [クイックスタート: scan](quickstart-scan-spec.md) |
| 手元のTS候補を検証 | [クイックスタート: tsopt](quickstart-tsopt-freq.md) |
| GPU GUIで対話的に実行 | [Colabを開く](https://colab.research.google.com/github/t-0hmura/mlmm_toolkit/blob/main/examples/mlmm_colab.ipynb) |
| 実行失敗・エラーを調べる | [典型エラー別レシピ](recipes-common-errors.md) |

## CLI サブコマンド

### メインワークフロー

| サブコマンド | 説明 |
|---|---|
| [`all`](all.md) | ML/MMモデル構築とMEP探索。TS・IRC・熱化学・DFTは任意 |

### 構造準備

| サブコマンド | 説明 |
|---|---|
| [`extract`](extract.md) | タンパク質–リガンド複合体からML領域を定義 |
| [`add-elem-info`](add-elem-info.md) | PDBの元素列（77–78）を補完 |
| [`mm-parm`](mm-parm.md) | Amberトポロジー・座標（parm7/rst7）を構築 |
| [`define-layer`](define-layer.md) | ML・可動MM・凍結MM層をB-factorで指定 |

### 構造最適化

| サブコマンド | 説明 |
|---|---|
| [`opt`](opt.md) | L-BFGSまたはRFOで構造最適化 |
| [`tsopt`](tsopt.md) | RS-P-RFO・DimerなどでTS候補を最適化 |

### 経路探索・最適化

| サブコマンド | 説明 |
|---|---|
| [`path-opt`](path-opt.md) | 2端点間のMEPをGSMまたはDMFで最適化 |
| [`path-search`](path-search.md) | MEP探索と再帰的な精密化 |

### スキャン

| サブコマンド | 説明 |
|---|---|
| [`scan`](scan.md) | 拘束付き距離スキャン。複数距離の協奏変化と多段階に対応 |
| [`scan2d`](scan2d.md) | 2次元エネルギー面 |
| [`scan3d`](scan3d.md) | 3次元エネルギー面 |

### 解析・後処理

| サブコマンド | 説明 |
|---|---|
| [`irc`](irc.md) | 固有反応座標を追跡 |
| [`freq`](freq.md) | 振動解析と熱化学 |
| [`dft`](dft.md) | GPU4PySCFまたはPySCFによるDFT一点計算 |
| [`sp`](sp.md) | ML/MM ONIOMエネルギー・力。任意でHessian |
| [`trj2fig`](trj2fig.md) | XYZ軌跡のエネルギープロファイルを描画 |
| [`energy-diagram`](energy-diagram.md) | 数値から状態エネルギー図を描画 |
| [`bond-summary`](bond-summary.md) | 構造間の共有結合変化を記録 |

### ユーティリティ

| サブコマンド | 説明 |
|---|---|
| [`fix-altloc`](fix-altloc.md) | PDBの代替コンフォメーションを解決 |

### エクスポート・インポート

| サブコマンド | 説明 |
|---|---|
| [`oniom-export`](oniom-export.md) | Gaussian ONIOMまたはORCA QM/MM入力を生成 |
| [`oniom-import`](oniom-import.md) | ONIOM入力からXYZ・層付きPDBを再構築 |

## 設定・リファレンス

| トピック | ページ |
|---|---|
| CLI規約と入力形式 | [CLI規約](cli-conventions.md) · [mmCIF](cif.md) |
| 概念・用語 | [概念](concepts.md) · [用語集](glossary.md) |
| YAML設定 | [YAMLリファレンス](yaml-reference.md) |
| 出力ファイル・JSON | [出力構造](output-layout.md) · [JSONスキーマ](json-output.md) |
| バックエンド・再現性 | [バックエンド](backends.md) · [再現性](reproducibility.md) |
| デバイス・HPC | [デバイスとHPC](device-hpc.md) |
| Python API・構成 | [Python API](python-api.md) · [ML/MM計算機](mlmm-calc.md) · [アーキテクチャ](architecture.md) |
| MCPサーバー | [MCPサーバー](mcp_server.md) |
| トラブルシューティング | [トラブルシューティング](troubleshooting.md) |
| 自動生成CLIリファレンス（英語） | [コマンドリファレンス](../reference/commands/index.md) |
| スターター設定（英語） | [YAML抜粋](../reference/yaml.md) |

## システム要件

インストールとバックエンドの互換性は[はじめに](getting-started.md#インストール)を参照してください。
GPUとドライバーは選択したバックエンドの要件を満たす必要があります。
VRAM・RAM・実行時間は、代表的な計算から見積もってください。
`mm-parm` にはAmberToolsが必要です。

## 重要な概念

- **層:** B=0はML、B=10は可動MM、B=20は凍結MMです。凍結原子もMMの非結合相互作用に寄与します。Hessianに含めるMM原子は `hess_cutoff` / `hess_mm_atoms` で別に選択します。
- **電荷・スピン:** `--ligand-charge` で残基電荷（例: `'SAM:1,GPP:-3'`）、`-q/--charge` でML領域の正味電荷、`-m/--multiplicity` で多重度（既定1）を指定します。
- **ブール値:** `--flag` / `--no-flag` で指定します。例: `--tsopt --thermo --no-dft`。
- **設定:** [YAMLリファレンス](yaml-reference.md)を参照してください。最適化せずに実効設定を確認する例:

```bash
mlmm opt -i layered.pdb --parm system.parm7 -q 0 --show-config --dry-run
```

## 出力構造

MEPモードの `all` は `summary.log`・`summary.json`、MEP（`mep.pdb` / `mep_trj.xyz`、bridge入力では `mep.cif` も）、
`energy_diagram_MEP.png` を出力します。
再利用できる準備ファイルは `ml_region.pdb`、`mm_parm/`、`layered/` です。
`segments/seg_NN/` にR/TS/P構造と指定したTS・IRC・freq・DFT結果、`_work/` に準備・scan・pathの中間出力を保存します。
TS-onlyモードではE1/TS/E2と表記し、MEPは出力しません。
全体のツリーは [all](all.md#出力)、ファイルの規約は[出力構造](output-layout.md)を参照してください。

## 引用

Ohmura, T., Inoue, S., Terada, T. (2025). *ML/MM toolkit — Toward Accelerated Mechanistic Investigation of Enzymatic Reactions.* [ChemRxiv](https://doi.org/10.26434/chemrxiv-2025-jft1k)。

## ライセンス

GNU General Public License version 3 or later (GPL-3.0-or-later)。

## ヘルプ

```bash
mlmm --help
mlmm <command> --help
```
