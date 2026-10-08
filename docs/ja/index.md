---
orphan: true
---

# [mlmm-toolkit]{.p2r-wordmark} ドキュメント

:::{container} p2r-hero-meta
[バージョン: v{{ release }}]{.p2r-pill} [GitHub](https://github.com/t-0hmura/mlmm_toolkit){.p2r-meta-gh} [ChemRxiv 論文](https://doi.org/10.26434/chemrxiv-2025-jft1k){.p2r-meta-paper}
:::

:::{container} p2r-hero
<img src="../mlmm_toolkit_overview.png" alt="mlmm-toolkit ワークフロー概要" class="p2r-hero-figure">

{.p2r-tagline}
**mlmm-toolkit** は、機械学習原子間ポテンシャルと分子力学を ONIOM 的に結合した **ML/MM 法** を用いて、酵素複合体などの PDB 構造から反応機構解析を行うための Python 製 CLI ツールキットです。

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
| **反応機構を調べる・TS が取れない** | [反応機構を調べるコツ](mechanism-tips.md) |
| **求めた TS 構造を DFT で構造最適化する** | [求めた TS 構造を DFT で構造最適化する](dft-backend.md) |
| **計算が失敗した** | [トラブルシューティング](troubleshooting.md) |

## サブコマンド

<!-- p2r-stage-strip -->

| サブコマンド | 説明 |
|---------|------|
| [`all`](all.md) | ML/MM モデルの構築と MEP 探索。TS・IRC・熱化学・DFT は任意 |
| [`fix-altloc`](fix-altloc.md) | PDB の代替コンフォメーションを解決 |
| [`add-elem-info`](add-elem-info.md) | PDB の元素列（77–78）を補完 |
| [`mm-parm`](mm-parm.md) | Amber のトポロジー・座標（parm7/rst7）を構築 |
| [`extract`](extract.md) | タンパク質–リガンド複合体から ML 領域を定義 |
| [`define-layer`](define-layer.md) | ML・可動 MM・凍結 MM の層を B-factor で指定 |
| [`opt`](opt.md) | L-BFGS または RFO で構造最適化 |
| [`scan`](scan.md) | 拘束付き距離スキャン。複数距離の協奏変化と多段階に対応 |
| [`scan2d`](scan2d.md) | 2 次元エネルギー面 |
| [`scan3d`](scan3d.md) | 3 次元エネルギー面 |
| [`path-opt`](path-opt.md) | 2 端点間の MEP を GSM または DMF で最適化 |
| [`path-search`](path-search.md) | MEP 探索と再帰的な精密化 |
| [`tsopt`](tsopt.md) | RS-P-RFO・Dimer などで TS 候補を最適化 |
| [`irc`](irc.md) | 固有反応座標を追跡 |
| [`freq`](freq.md) | 振動解析と熱化学 |
| [`dft`](dft.md) | GPU4PySCF または PySCF による DFT 一点計算 |
| [`sp`](sp.md) | ML/MM ONIOM のエネルギー・力。任意で Hessian |
| [`bond-summary`](bond-summary.md) | 構造間の共有結合の変化を記録 |
| [`trj2fig`](trj2fig.md) | XYZ 軌跡のエネルギープロファイルを描画 |
| [`energy-diagram`](energy-diagram.md) | 数値から状態エネルギー図を描画 |
| [`oniom-export`](oniom-export.md) | Gaussian ONIOM または ORCA QM/MM の入力を生成 |
| [`oniom-import`](oniom-import.md) | ONIOM の入力から XYZ・層付き PDB を再構築 |

## 設定・参照資料

| トピック | ページ |
|-------|------|
| **共通オプションと入力要件** | [共通オプションと残基・原子の指定](cli-conventions.md) |
| **原子の固定と距離の拘束（`--freeze-atoms`・`--distance-restraint`）** | {ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>` |
| **よくあるエラーと対処** | [トラブルシューティング](troubleshooting.md) |
| **CLI コマンドの一覧（英語のみ、自動生成）** | [コマンドの一覧（英語のみ）](../reference/commands/index.md) |
| **YAML 設定オプション** | [YAML 設定の一覧](yaml-reference.md) · [YAML の抜粋（英語のみ）](../reference/yaml.md) |
| **MLIP バックエンド設定** | [MLIP バックエンド](backends.md) |
| **各コマンドが書き出すファイル** | [出力ディレクトリのレイアウト](output-layout.md) |
| **`result.json` と `summary.json` の欄** | [JSON 出力の一覧](json-output.md) |
| **GPU・CPU の割り当てと HPC** | [デバイス設定 & HPC セットアップ](device-hpc.md) |
| **Python から ML/MM 計算機を使う** | [ML/MM 計算機](mlmm-calc.md) |
| **AI エージェントから呼ぶ（MCP）** | [MCP サーバー](mcp_server.md) |
| **コードの構成（開発者向け）** | [アーキテクチャ](architecture.md) |
| **用語** | [用語集](glossary.md) |

## システム要件

### ハードウェア
- **OS**: Linux（Windows では WSL2 上の Linux に導入してください）
- **GPU**: 使用するバックエンドと PyTorch wheel に対応する NVIDIA ドライバー。CPU のみでも実行可能ですが低速です
- **VRAM / RAM**: モデル、系の大きさ、Hessian の計算方式で変わります。代表的な計算で最大使用量を測ってください

### ソフトウェア
- Python >= 3.11
- CPU 版または CUDA 対応の PyTorch。ビルド済みの wheel は CUDA ランタイムを含むので、手元の CUDA toolkit は通常いりません（ソースからビルドするときだけ必要です）
- AmberTools（`mm-parm` が Amber のトポロジーを作るのに使います）

セットアップは [インストール](installation.md) を参照してください。

## 重要な概念

3 つの層（ML・可動 MM・凍結 MM）と ONIOM での組み合わせ方は [はじめに](getting-started.md)、ML 領域と層の決め方は [ML 領域と層の組み方](model-setup.md) を参照してください。

## エージェントスキル

`skills/` に、AI エージェント向けの手順書（CLI コマンド・構造 I/O・バックエンド・ワークフローと出力・HPC 運用）を同梱しています。導入するときは、AI エージェントに次のように指示してください。

> `https://github.com/t-0hmura/mlmm_toolkit/tree/main/skills` をスキルとして取り込み、`mlmm-install-backends` の手順に従って mlmm-toolkit をインストールして

GitHub のリポジトリを clone 済みなら、URL の代わりに手元の `skills/` の path を渡しても構いません。導入した後は、たとえば次のように頼めます。

> 〈論文〉を読んで、〈PDB ID〉の構造からモデルを作成し、〈反応段階〉の経路について、mlmm-toolkit のスキルを用いて反応機構解析を行ってください。

## 引用

Ohmura, T., Inoue, S., Terada, T. (2025). *ML/MM toolkit — Toward Accelerated Mechanistic Investigation of Enzymatic Reactions.* [ChemRxiv](https://doi.org/10.26434/chemrxiv-2025-jft1k)。

## ライセンス

GNU General Public License version 3 or later (GPL-3.0-or-later)。

## ヘルプ

```bash
mlmm --help
mlmm <command> --help
mlmm <command> --help-advanced
```

問題や機能リクエストは [GitHub Issues](https://github.com/t-0hmura/mlmm_toolkit/issues) に報告してください。
