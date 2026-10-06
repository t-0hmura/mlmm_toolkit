---
orphan: true
---

# mlmm-toolkit ドキュメント

[GitHub](https://github.com/t-0hmura/mlmm_toolkit) · [ChemRxiv 論文](https://doi.org/10.26434/chemrxiv-2025-jft1k) · [Google Colab で実行](https://colab.research.google.com/github/t-0hmura/mlmm_toolkit/blob/main/examples/mlmm_colab.ipynb)

*バージョン: v{{ release }}*

---

<img src="../mlmm_toolkit_overview.png" alt="mlmm-toolkit workflow overview" width="90%">

**mlmm-toolkit** は、機械学習原子間ポテンシャルと分子力学を ONIOM 的に結合した **ML/MM 法** を用いて、PDB 構造から酵素反応経路を自動モデリングする Python 製 CLI ツールキットです。

初めての方は [はじめに](getting-started.md) からお読みください。

## 目的別クイックスタート

| 目的 | ガイド |
|---|---|
| 反応の前後の構造から反応機構解析を一気通貫で行う | [クイックスタート: all の Endpoint モード](quickstart-all.md) |
| 1 つの構造から一気通貫で反応機構解析を行う | [クイックスタート: all の Scan-list モード](quickstart-scan.md) |
| TS 構造から一気通貫で反応機構解析を行う | [クイックスタート: TS-only モード](quickstart-tsopt.md) |
| ML 領域と層を決める・計算を軽くする | [ML 領域と層の組み方](model-setup.md) |
| 反応機構を調べる・TS が取れない | [反応機構を調べるコツ](mechanism-tips.md) |
| 求めた TS 構造を DFT で構造最適化する | [求めた TS 構造を DFT で構造最適化する](dft-backend.md) |
| 計算が失敗した | [トラブルシューティング](troubleshooting.md) |

前提条件は [インストール](installation.md) を参照してください。

## CLI サブコマンド

| サブコマンド | 説明 |
|---|---|
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

## 設定・リファレンス

| トピック | ページ |
|---|---|
| CLI 規約と入力形式 | [共通オプションと残基・原子の指定](cli-conventions.md) |
| 原子の固定と距離の拘束（`--freeze-atoms`・`--distance-restraint`） | {ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>` |
| 用語 | [用語集](glossary.md) |
| YAML 設定 | [YAML リファレンス](yaml-reference.md) |
| 出力ファイル・JSON | [出力構造](output-layout.md) · [JSON スキーマ](json-output.md) |
| バックエンド・再現性 | [バックエンド](backends.md) |
| デバイス・HPC | [デバイスと HPC](device-hpc.md) |
| Python API・構成 | [ML/MM 計算機](mlmm-calc.md) · [アーキテクチャ](architecture.md) |
| MCP サーバー | [MCP サーバー](mcp_server.md) |
| トラブルシューティング | [トラブルシューティング](troubleshooting.md) |
| 自動生成 CLI リファレンス（英語） | [コマンドリファレンス](../reference/commands/index.md) |
| スターター設定（英語） | [YAML 抜粋](../reference/yaml.md) |

## システム要件

インストールとバックエンドの互換性は [インストール](installation.md) を参照してください。
GPU とドライバーは選択したバックエンドの要件を満たす必要があります。
VRAM・RAM・実行時間は、代表的な計算から見積もってください。
`mm-parm` には AmberTools が必要です。

## 重要な概念

3 つの層（ML・可動 MM・凍結 MM）と ONIOM での組み合わせ方は [はじめに](getting-started.md#概要)、ML 領域と層の決め方は [ML 領域と層の組み方](model-setup.md) を参照してください。

## エージェントスキル

`skills/` に、AI エージェント向けの手順書（CLI コマンド・構造 I/O・バックエンド・ワークフローと出力・HPC 運用）を同梱しています。導入するときは、AI エージェントに次のように指示してください。

> `https://github.com/t-0hmura/mlmm_toolkit/tree/main/skills` をスキルとして取り込んで

clone 済みなら、URL の代わりに手元の `skills/` の path を渡しても構いません。

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
