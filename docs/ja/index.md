---
orphan: true
---

# [mlmm-toolkit]{.p2r-wordmark} ドキュメント

:::{container} p2r-hero-meta
[バージョン: v{{ release }}]{.p2r-pill} [GitHub](https://github.com/t-0hmura/mlmm_toolkit){.p2r-meta-gh} [ChemRxiv 論文](https://doi.org/10.26434/chemrxiv-2025-jft1k){.p2r-meta-paper}
:::

:::{container} p2r-hero
<img src="../overview.jpg" alt="mlmm-toolkit ワークフロー概要" class="p2r-hero-figure">

{.p2r-tagline}
**mlmm-toolkit** は、酵素複合体を対象として、機械学習原子間ポテンシャル（MLIP）と分子力学（MM）を ONIOM で結合した ML/MM 法を用いて PDB/mmCIF の構造から反応機構を解析するための Python 製 CLI ツールキットです。

{.p2r-lead}
初めての方は [はじめに](getting-started.md) からお読みください。

{.p2r-cta}
[はじめに](getting-started.md){.p2r-btn .p2r-btn-primary} [インストール](installation.md){.p2r-btn .p2r-btn-install} [Google Colabで実行](https://colab.research.google.com/github/t-0hmura/mlmm_toolkit/blob/main/examples/mlmm_colab.ipynb){.p2r-btn .p2r-btn-colab}
:::

## クイックスタート

::::{container} p2r-cards
:::{container} p2r-card p2r-card-endpoint
**反応の前後の構造から反応機構解析を一気通貫で行う**

<!-- p2r-mode-stages endpoint -->

[クイックスタート: all の Endpoint モード](quickstart-all.md)
:::

:::{container} p2r-card p2r-card-scan
**1 つの構造から一気通貫で反応機構解析を行う**

<!-- p2r-mode-stages scan -->

[クイックスタート: all の Scan-list モード](quickstart-scan.md)
:::

:::{container} p2r-card p2r-card-tsonly
**TS 構造から一気通貫で反応機構解析を行う**

<!-- p2r-mode-stages tsonly -->

[クイックスタート: TS-only モード](quickstart-tsopt.md)
:::
::::

| 目的 | ページ |
|------|------|
| **ML 領域と層を決める・計算を軽くする** | [ML 領域と層の組み方](model-setup.md) |
| **反応機構解析・遷移状態探索を行う** | [反応機構を調べるコツ](mechanism-tips.md) |
| **MLIP で求めた TS 構造を DFT で構造最適化する** | [MLIP での TS 構造を DFT で検証する](dft-backend.md) |
| **計算の失敗原因を調査する** | [トラブルシューティング](troubleshooting.md) |

## サブコマンド

<!-- p2r-stage-strip -->

| サブコマンド | 説明 |
|---------|------|
| [`all`](all.md) | ML/MM モデルの構築と MEP 探索。TS・IRC・熱化学・DFT は任意 |
| [`mm-parm`](mm-parm.md) | Amber のトポロジー・座標（parm7/rst7）を構築 |
| [`extract`](extract.md) | タンパク質–リガンド複合体から ML 領域を定義 |
| [`define-layer`](define-layer.md) | ML・可動 MM・固定 MM の層を B-factor で指定 |
| [`opt`](opt.md) | 単一構造の構造最適化（L-BFGS または RFO が利用可能） |
| [`scan`](scan.md) | 距離・角度・二面角を変数としてエネルギー変化をスキャン（協奏スキャン・多段階スキャンに対応） |
| [`scan2d`](scan2d.md) | 2 次元のエネルギー地形の探索・2D PES のマッピング |
| [`scan3d`](scan3d.md) | 3 次元のエネルギー地形の探索・3D PES のマッピング |
| [`path-opt`](path-opt.md) | Growing String 法（GSM）または Direct Max Flux 法による 1 段階の最小エネルギー経路（MEP）最適化。反応物と生成物の 2 構造から反応に伴う構造変化のトラジェクトリを構築する。 |
| [`path-search`](path-search.md) | 複数のステップからなる反応を自動検出して多段階の最小エネルギー経路を構築する再帰的 MEP 探索。反応物・中間体・生成物構造を入力できる。 |
| [`tsopt`](tsopt.md) | 単一構造の遷移状態（TS）構造最適化（Dimer または RS-P-RFO）。虚振動が一つの構造を得る。 |
| [`irc`](irc.md) | 固有反応座標（IRC: Intrinsic Reaction Coordinate）計算 |
| [`freq`](freq.md) | 振動数解析と熱化学解析 |
| [`dft`](dft.md) | DFT 一点計算（GPU4PySCF / PySCF） |
| [`fix-altloc`](fix-altloc.md) | PDB の代替位置指示子を解決 |
| [`add-elem-info`](add-elem-info.md) | PDB の元素列（77–78）を修復 |
| [`sp`](sp.md) | ML/MM ONIOM のエネルギー・力。任意で Hessian |
| [`bond-summary`](bond-summary.md) | 複数の状態の構造間の共有結合変化を検出・レポート |
| [`trj2fig`](trj2fig.md) | mlmm-toolkit が生成する XYZ トラジェクトリからエネルギープロファイルをプロット |
| [`energy-diagram`](energy-diagram.md) | 数値入力からエネルギーダイアグラムを作成 |
| [`oniom-export`](oniom-export.md) | Gaussian ONIOM または ORCA QM/MM の入力を生成 |
| [`oniom-import`](oniom-import.md) | ONIOM の入力から XYZ・層付き PDB を再構築 |

## 設定・参照資料

| トピック | ページ |
|-------|------|
| **共通オプションと入力要件** | [共通オプションと残基・原子の指定](cli-conventions.md) |
| **原子の固定と距離の拘束（`--freeze-atoms`・`--distance-restraint`）** | {ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>` |
| **よくあるエラーと対処** | [トラブルシューティング](troubleshooting.md) |
| **CLI コマンドの一覧** | [コマンドの一覧](../reference/commands/index.md) |
| **YAML 設定オプション** | [YAML 設定の一覧](yaml-reference.md) · [YAML の抜粋（英語のみ）](../reference/yaml.md) |
| **MLIP バックエンド設定** | [MLIP バックエンド](backends.md) |
| **各コマンドが書き出すファイル・ディレクトリの詳細** | [出力ディレクトリのレイアウト](output-layout.md) |
| **結果を報告する JSON の項目詳細** | [JSON 出力の一覧](json-output.md) |
| **GPU・CPU の割り当てと HPC** | [デバイス設定 & HPC セットアップ](device-hpc.md) |
| **Python から ML/MM 計算機を使う** | [ML/MM 計算機](mlmm-calc.md) |
| **AI エージェントに mlmm-toolkit を使わせる** | [エージェントスキル](#エージェントスキル) |
| **コードの構成（開発者向け）** | [アーキテクチャ](architecture.md) |
| **用語** | [用語集](glossary.md) |

## システム要件

- **OS**: Linux（Windows では WSL2、Mac では Docker を推奨）
- **Python**: 3.12（推奨）または 3.13。ORB は 3.12 が必要です
- **PyTorch**: CUDA 版（推奨）または CPU 版
- **GPU**: VRAM 8 GB 以上。CPU のみでも実行可能ですが低速です。対応する NVIDIA ドライバーが必要で、UMA は実行時に GCC と CUDA Toolkit も必要です
- **RAM / VRAM**: デフォルトの UMA での一点計算は、約 3,400 原子（ML 領域は約 650 原子）と約 9,200 原子（ML 領域は約 600 原子）で RAM 8 GB・VRAM 8 GB、約 14,000 原子（ML 領域は約 370 原子）で RAM 16 GB・VRAM 8 GB が目安で、Hessian もほぼ同じです。RAM は全系、VRAM は ML 領域の大きさで増えます
- AmberTools（`mm-parm` が Amber のトポロジーを作るのに使います）

セットアップは [インストール](installation.md) を参照してください。

## 重要な概念

3 つの層（ML・可動 MM・固定 MM）と ONIOM での組み合わせ方は [はじめに](getting-started.md)、ML 領域と層の決め方は [ML 領域と層の組み方](model-setup.md) を参照してください。

## エージェントスキル

`mlmm-toolkit` は、ツールの使い方や反応機構解析の手順をまとめた AI エージェント向けの手順書を `skills/` に同梱しています。

次の内容をカバーしています：CLI サブコマンド・構造 I/O・バックエンドインストール・ワークフロー・出力解析・HPC 運用

導入するときは、AI エージェントに次のように指示してください。

> `https://github.com/t-0hmura/mlmm_toolkit/tree/main/skills` をスキルとして取り込み、`mlmm-install` の手順に従って mlmm-toolkit をインストールしてください。

GitHub のリポジトリを clone 済みなら、URL の代わりに手元の `skills/` のパスを渡しても構いません。導入後は、たとえば次のように指示できます。

> 〈論文〉を読んで、〈PDB ID〉の構造からモデルを作成し、〈反応段階〉の経路について、mlmm-toolkit のスキルを用いて反応機構解析を行ってください。

詳しくは [AI エージェント向けスキル](agent-skills.md) を参照してください。

## 引用

`mlmm-toolkit` を研究で利用する場合は、以下の論文を引用してください:

```bibtex
@article{ohmura2025mlmm,
  author       = {Ohmura, Takuto and Inoue, Sei and Terada, Tohru},
  title        = {ML/MM toolkit -- Towards Accelerated Mechanistic Investigation of Enzymatic Reactions},
  journal      = {ChemRxiv},
  year         = {2025},
  doi          = {10.26434/chemrxiv-2025-jft1k}
}
```

ソフトウェアまたは特定のリリースを引用する場合は、Zenodo レコードを使用してください:

```bibtex
@software{ohmura2026mlmm_software,
  author       = {Ohmura, Takuto},
  title        = {mlmm-toolkit},
  year         = {2026},
  version      = {0.6.0},
  url          = {https://github.com/t-0hmura/mlmm_toolkit},
  license      = {GPL-3.0-or-later},
  doi          = {10.5281/zenodo.19197863}
}
```

## ライセンス

`mlmm-toolkit` は **GNU General Public License version 3 or later (GPL-3.0-or-later)** の下で配布されています。

## ヘルプ

```bash
# 一般的なヘルプ
mlmm --help

# コマンドのヘルプ
mlmm <subcommand> --help

# 詳細オプション（内部チューニング用）
mlmm <subcommand> --help-advanced
```

問題や機能リクエストは [GitHub Issues](https://github.com/t-0hmura/mlmm_toolkit/issues) に報告してください。
