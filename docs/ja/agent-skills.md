# AI エージェント向けスキル

`mlmm-toolkit` には、AI エージェント（Claude Code、Codex、Cursor など）向けの手順書（スキル）が `skills/` ディレクトリに同梱されています。`skills/` をエージェントに読み込ませることで、自然言語による指示で計算の実行と解析を行えます。

## 最初の例

導入するときは、AI エージェントに次のように指示してください。

> `https://github.com/t-0hmura/mlmm_toolkit/tree/main/skills` をスキルとして取り込み、`mlmm-install` の手順に従って mlmm-toolkit をインストールして

GitHub のリポジトリを clone 済みなら、URL の代わりに手元の `skills/` の path を渡しても構いません。導入後は、たとえば次のように指示できます。

> 〈論文〉を読んで、〈PDB ID〉の構造からモデルを作成し、〈反応段階〉の経路について、mlmm-toolkit のスキルを用いて反応機構解析を行ってください。

## スキルの一覧

`skills/` 内の各フォルダが 1 つのスキル（`<名前>/SKILL.md`）に対応します。

| スキル | 内容 |
| --- | --- |
| `mlmm-overview`（最初に参照） | mlmm-toolkit の概要、構造に応じた `all` の 3 つの入力モードの選択、TS 探索の方針、出力の解釈 |
| `mlmm-cli` | 各サブコマンドの用途、成否の判定、注意点と対処 |
| `mlmm-model-setup` | 構造ファイルと Amber parm7・rst7 の取り扱い、電荷とスピン多重度の決定、ML 領域と層の設定、状態間・変異体間での原子の対応 |
| `mlmm-install` | mlmm-toolkit 本体・MLIP バックエンド・DFT・xTB・AmberTools の導入、計算環境（スケジューラ・GPU・CUDA・conda）の確認 |
| `mlmm-hpc` | PBS・SLURM のジョブスクリプト、walltime の目安、ジョブの監視、動的なジョブ割り当て |
| `colab-local-gpu-runtime` | WSL2 と Docker Desktop を用いた、ローカルの NVIDIA GPU 上での Colab インターフェースの実行（Windows） |

## 手動で導入する場合

スキルのフォルダを、使用するエージェントのスキル用ディレクトリにコピーします。

* Claude Code：プロジェクト単位では `.claude/skills/`、全プロジェクト共通では `~/.claude/skills/`
* Codex：リポジトリ単位では `.agents/skills/`、全リポジトリ共通では `~/.agents/skills/`

例として、リポジトリのルートで次のコマンドを実行します。

```bash
mkdir -p ~/.claude/skills
cp -r skills/mlmm-* skills/colab-local-gpu-runtime ~/.claude/skills/
```

オプションの詳細と既定値は、インストールした CLI の `mlmm <subcommand> --help-advanced` と[コマンドリファレンス](../reference/commands/index.md)で確認してください。

## 関連ドキュメント

* [はじめに](getting-started.md) — 最短の実行と、次に読むページ
* [インストール](installation.md) — 環境の準備
* [`skills/README.md`](https://github.com/t-0hmura/mlmm_toolkit/blob/main/skills/README.md) — GitHub 上のスキルの一覧
