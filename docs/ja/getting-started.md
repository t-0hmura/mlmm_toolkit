# はじめに

## 概要

<img src="../mlmm_toolkit_overview.png" alt="mlmm-toolkit workflow overview" width="90%">

`mlmm-toolkit` は、機械学習原子間ポテンシャル（MLIP）と内蔵 MM 力場エンジンを ONIOM 的に結合した **ML/MM 法** を用いて、**PDB 構造** から **酵素反応経路** を自動的に構築する Python 製の CLI ツールキットです。デフォルトの MLIP バックエンドは **UMA**（Meta の FAIR-Chem）で、`--backend` オプションにより **ORB**、**MACE**、**AIMNet2** も選択できます。

多くのワークフローで、**1 コマンド**で反応経路の**初期推定**を得られます。
```bash
mlmm -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
```

---
さらに `--tsopt --thermo --dft` を追加すると、**ML/MM モデル構築 → MEP 探索 → TS 最適化 → IRC → 熱化学補正 → DFT 一点計算** までまとめて実行できます。
```bash
mlmm -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt --thermo --dft
```
---

入力として、(i) 反応順に並べたタンパク質-リガンド複合体の PDB を 2 つ以上（R →... → P）、(ii) `--scan-lists` を指定した 1 つの PDB、または (iii) TS 候補 1 構造 + `--tsopt` を与えると、`mlmm-toolkit` が次の処理を自動化します。

- ユーザーが指定した基質の周辺から **活性部位ポケット** を抽出し、**ML 領域** を定義
- AmberTools を用いて **Amber トポロジー（parm7/rst7）** を自動生成し、**hessian_ff** の MM エンジンに渡す
- ML 領域を MLIP バックエンド（デフォルト: UMA）で、MM 領域を hessian_ff で扱う **ONIOM 的 ML/MM** のエネルギー・力・Hessian を構築
- Growing String Method (GSM) や Direct Max Flux (DMF) などの経路最適化手法で **最小エネルギー経路 (MEP)** を探索
- 必要に応じて **遷移状態** を最適化し、**振動解析**・**IRC 計算**・**DFT 一点計算** を実行

```{important}
TSOPT 終端の振動解析で虚振動がちょうど 1 つあることを確認し、IRC と端点最適化で目的の反応物・生成物につながるか検証してください。追加の `freq` は、全振動モードや熱化学量が必要な場合に実行します。
```

MM 領域の計算には hessian_ff（内蔵の C++ ネイティブ MM 力場エンジン）を用います。全エネルギーは ONIOM 的な減算分解に従います:

```
E_total = E_REAL_low + E_MODEL_high - E_MODEL_low
```

ここで REAL は全系、MODEL は ML 領域、"high" は MLIP バックエンド、"low" は hessian_ff です。

一連の処理は CLI から呼び出せるように統一されており、手作業を最小化して **多段階の酵素反応メカニズム** を組み立てられるように設計されています。同じワークフローは小分子系にも適用可能です。`.xyz` 入力を個別計算で使う場合は、対応する全系トポロジーを `--parm`、構造テンプレートを `--ref-pdb`、ML 領域を `--model-pdb`、`--model-indices`、または有効な B-factor layer で指定します。

```{important}
- 入力 PDB ファイルには**水素原子**が含まれている必要があります。
- 複数の PDB を提供する場合、**同じ原子が同じ順序**で含まれている必要があります（座標のみ異なる可能性があります）。そうでない場合はエラーが発生します。
- 個別の ML/MM 計算には **`--parm`**（全系の Amber トポロジー）と、`--model-pdb`、`--model-indices`、または有効な B-factor layer のいずれかによる ML 領域指定が必要です。`all` ワークフローではトポロジーと ML 領域を自動生成できます。
- `mlmm all` と個別コマンドのどちらでも、`-q/--charge` は全系ではなく ML 領域（ONIOM モデル系）の正味電荷です。
- MD スナップショットには、MD 計算で用いた全系の `.parm7` を再利用してください。
```

```{tip}
初めて使う場合は、まず [概念とワークフロー](concepts.md) を参照してください。
症状から切り分ける場合は、まず [典型エラー別レシピ](recipes-common-errors.md) を参照してください。
セットアップや実行中にエラーが発生した場合は [トラブルシューティング](troubleshooting.md) を参照してください。
```

### 対話型 Colab GUI

[mlmm Colab ノートブック](https://colab.research.google.com/github/t-0hmura/mlmm_toolkit/blob/main/examples/mlmm_colab.ipynb)では、PDB/mmCIF 構造と対応する全系 `parm7` のアップロード、3D での ML 領域選択、生成コマンドの検証と実行、現在の呼び出しで生成された結果だけの確認ができます。各ユーザーは専用の GPU ランタイムで実行します。MACE と ORB はモデル利用のログインが不要ですが、UMA には Hugging Face のアクセス許可が必要です。互換性のないバックエンドへ切り替える場合は、ランタイムを再起動してください。DFT の操作項目は、Setup で DFT の追加依存関係を選択した場合だけ表示されます。Setup は指定バージョンの PyPI wheel をインストールし、対応する Git tag からサンプルを取得します。このため、本番ノートブックを実行できるのは対象 wheel の公開後です。

### CLI の慣習

| 慣習 | 例 | 備考 |
|-----|-----|------|
| **残基セレクタ** | `'SAM,GPP'`, `'A:123,B:456'` | 複数値はシェル展開防止のためクォート |
| **電荷マッピング** | `-l 'SAM:1,GPP:-3'` | `all` / `extract` などではコロン（`:`）で名前と電荷を区切る。`mm-parm` は互換用に `=` も受理 |
| **原子セレクタ** | `'TYR,285,CA'` または `'TYR 285 CA'` | 区切り文字: 空白、カンマ、スラッシュ、バッククォート、バックスラッシュ |

詳細は [CLI 規約](cli-conventions.md) を参照してください。

### 水素原子付与の推奨ツール

PDB に水素原子がない場合は、mlmm を実行する前に次のいずれかを使ってください。

| ツール | コマンド例 | 備考 |
|--------|------------|------|
| **reduce** (Richardson Lab) | `reduce input.pdb > output.pdb` | 高速、結晶構造に広く使用 |
| **pdb2pqr** | `pdb2pqr --ff=AMBER input.pdb output.pqr` | 水素を追加し部分電荷を割り当て |
| **Open Babel** | `obabel input.pdb -O output.pdb -h` | 汎用ケモインフォマティクスツールキット |
| **mm-parm --add-h** | `mlmm mm-parm -i input.pdb --add-h` | PDBFixer が必要（`pip install "mlmm-toolkit[pdbfixer]"` または `conda install -c conda-forge pdbfixer`） |

複数の PDB 入力で同一の原子順序を確保するには、すべての構造に同じ水素付与ツールを一貫した設定で適用してください。

```{warning}
このソフトウェアはまだ開発中です。自己責任でご使用ください。
```

---

## インストール

Linux の CPU/GPU 環境で利用できます。GPU 実行には対応する NVIDIA ドライバーが必要です。公式 PyTorch wheel には CUDA ランタイムが含まれるため、通常は CUDA toolkit を別途インストールする必要はありません。

以下は PyTorch 2.13 の `cu130` wheel を使う例です。CPU 実行や別の GPU 環境では、対応する PyTorch wheel を選んでください。MM 計算には C++20 対応コンパイラー、トポロジー生成には AmberTools が必要です。

```bash
conda create -n mlmm-toolkit python=3.12 -y
conda activate mlmm-toolkit
conda install -c conda-forge ambertools=24.8 "numpy>=2,<2.5" pdbfixer -y
pip install torch==2.13.0 --index-url https://download.pytorch.org/whl/cu130
pip install mlmm-toolkit

# UMA の利用許諾を取得した後、Hugging Face にログイン
hf auth login
mlmm --version
```

UMA を使う場合は、[モデルページ](https://huggingface.co/facebook/UMA)で FAIR Chemistry License v1 に同意してください。ログインは環境ごとに一度行います。

### 追加コンポーネント

| 用途 | インストール・設定 |
| --- | --- |
| ORB / AIMNet2 | ORB は Python 3.11／3.12 が必要です（3.12 推奨）。`pip install --only-binary=dm-tree "mlmm-toolkit[orb]"` / `pip install "mlmm-toolkit[aimnet]"` |
| MACE | UMA と `e3nn` の依存バージョンが競合するため、専用環境で使用します。 |
| DMF 経路探索 | `conda install -c conda-forge cyipopt -y` と `pip install 'pydmf>=1.2'` |
| Plotly の PNG 出力 | `plotly_get_chrome -y` |
| hessian_ff の手動ビルド | 初回使用時に JIT コンパイルされます。ネイティブ拡張を利用できない場合は、下記を実行してください。 |

```bash
conda install -c conda-forge ninja -y
cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make
```

`hessian_ff` はビルド済みキャッシュを確認し、必要な場合は自動ビルドします。失敗時の確認事項と手動再ビルドは、[トラブルシューティング](troubleshooting.md)を参照してください。C/CUDA 拡張をソースからビルドする場合の toolkit/compiler 設定やジョブスクリプトは、[デバイスと HPC](device-hpc.md)を参照してください。

---

## マルチバックエンドの使用例

デフォルトの MLIP バックエンドは UMA です。`-b/--backend` で代替バックエンドに切り替えます:

この例では、`ml_region.pdb` は `real.parm7` に対応する全系の構造で、`ml.pdb` が ML 領域を指定します。

```bash
# ORB バックエンドを使用
mlmm opt -i ml_region.pdb --parm real.parm7 --model-pdb ml.pdb -q 0 -b orb

# MACE バックエンドを使用
mlmm opt -i ml_region.pdb --parm real.parm7 --model-pdb ml.pdb -q 0 -b mace

```

---

## 推奨クイックスタート

- [クイックスタート: `mlmm all`](quickstart-all.md)
- [クイックスタート: `mlmm scan`](quickstart-scan-spec.md)
- [クイックスタート: `mlmm tsopt`](quickstart-tsopt-freq.md)

---

## 典型的な手動ワークフロー

再利用可能なトポロジーと PDB を個別サブコマンドで準備する場合は、まず次を実行します。

```bash
mlmm mm-parm -i input.pdb -l 'LIG:0' --out-prefix system
mlmm extract -i system.pdb -c LIG -l 'LIG:0' -o model.pdb
mlmm define-layer -i system.pdb --model-pdb model.pdb -o system_layered.pdb
```

```text
1. mm-parm - parm7/rst7 と LEaP のトポロジー対応 PDB を生成
2. extract - その生成 PDB から活性部位ポケットを抽出
3. define-layer - 同じ生成 PDB に 3 層 ML/MM 分割を付与（B-factor エンコード）
4. all の MEP stage - 単一パス path-opt がデフォルト。`mlmm all --refine-path` で再帰 path-search に切替
5. tsopt - 遷移状態最適化
6. freq - 振動解析と熱化学
7. dft - DFT 一点計算
```

LEaP が水素を変更する場合があるため、2 以降では `mm-parm` が出力した
PDB を使用します。明示的な `--out-prefix` でこの PDB を出力でき、空の元素記号列は
原子レコードと順序を保ったまま補完されます。`all` は同等の準備を内部管理し、内部では
`extract → mm-parm → define-layer` の順に処理します。この内部順序は、単独ファイルを
手動で再利用するための手順ではありません。各ステップは単独でも実行できます。

---

## コマンドラインの基本

`mlmm` のデフォルトのサブコマンドは `all` です。

```bash
mlmm [OPTIONS]...
# は以下と同等
mlmm all [OPTIONS]...
```

`all` ワークフローは、ML 領域抽出、MM パラメータ化、レイヤー定義、MEP 探索、TS 最適化、振動解析、DFT 一点計算（任意）を 1 つのコマンドで連続実行する**統合コマンド**です。

ML 領域抽出を使用する場合、すべての上位ワークフローで共通する重要なオプションが 2 つあります:

- `-i/--input`: 1 つ以上の**完全系構造**（反応物、中間体、生成物）。
- `-c/--center`: **基質/抽出中心**の定義方法（例: 残基名や残基 ID）。

`--center/-c` を省略すると、ML 領域抽出はスキップされ、**入力構造全体**がそのまま使用されます。

---

## メインワークフローモード

### 複数構造からの MEP 探索

反応順に並べた、同じ原子・原子順序の全系構造を 2 つ以上指定します。

```bash
mlmm -i R.pdb I1.pdb I2.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
     --out-dir ./result_all --tsopt --thermo --dft
```

デフォルトは隣接ペアごとの単一パス `path-opt` です。`--refine-path` を指定すると再帰的な `path-search` に切り替わります。どちらも GSM/DMF を選択できます。

### 単一構造とスキャン定義

変化させる原子間距離が分かっている場合は、1 構造に `--scan-lists` を併用します。

```bash
mlmm -i R.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
     --scan-lists '[("TYR 285 CA","MMT 309 C10",2.20),("TYR 285 CB","MMT 309 C11",1.80)]' \
                  '[("TYR 285 CB","MMT 309 C11",1.20)]'
```

各タプル `(i, j, target_Å)` には PDB 原子セレクタまたは 1 始まりの原子番号を指定します。1 リテラル内の距離は同時に変化させ、複数のリテラルは順に実行します。複数リテラルは、1 つの `--scan-lists` の後に続けてください。

### TS 候補からの最適化と IRC

TS 候補を 1 つ指定して `--tsopt` を有効にすると、MEP 探索を省略します。

```bash
mlmm -i TS_CANDIDATE.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt --thermo
```

IRC 後の端点 E1/E2 は未割当です。最適化した構造を確認してから、反応物・生成物を割り当ててください。各モードの処理と出力は [all](all.md) を参照してください。

```{important}
単一入力には `--scan-lists` または `--tsopt` が必要です。
```

---

## 重要な CLI オプションと動作

| オプション | 説明 |
|----------|------|
| `-i, --input PATH...` | 入力構造。**2 つ以上の PDB** → MEP 探索; **1 つの PDB + `--scan-lists`** → 段階的スキャン; **1 つの PDB + `--tsopt`** → TSOPT のみ |
| `-c, --center TEXT` | 基質/抽出中心を定義。残基名（`'SAM,GPP'`）、残基ID（`A:123,B:456`）、または PDB パスをサポート |
| `-l, --ligand-charge TEXT` | 電荷情報: マッピング（`'SAM:1,GPP:-3'`）または単一整数 |
| `-q, --charge INT` | ML 領域の総電荷の強制上書き |
| `-m, --multiplicity INT` | スピン多重度（例: 一重項は `1`） |
| `-s, --scan-lists TEXT...` | `all`の単一入力経路ではインライン`(i,j,target)`リテラル。YAML/JSONと双方向4-tupleはstandalone `scan`で使用 |
| `--parm PATH` | 全系の Amber parm7 トポロジー（`all` では自動生成） |
| `--model-pdb PATH` | ML 領域を定義する PDB ファイル。個別計算では `--model-indices` または有効な B-factor layer も選択可能（`all` では自動生成可） |
| `--tsopt/--no-tsopt` | TS 最適化と IRC を有効化 |
| `--thermo/--no-thermo` | 振動解析と熱化学を実行 |
| `--dft/--no-dft` | DFT 一点計算を実行 |
| `--refine-path/--no-refine-path` | `mlmm all` で単一パス `path-opt`（デフォルト）または再帰 `path-search` を選択 |
| `--mep-mode gsm\|dmf` | どちらの経路探索にも用いる MEP 最適化法（デフォルト: `gsm`） |
| `--dmf-backend gpu\|cpu` | DMF 実装。GPU メモリ不足時は `cpu` を選択 |
| `-o, --out-dir PATH` | トップレベル出力ディレクトリ |
| `-b, --backend uma\|orb\|mace\|aimnet2` | MLIP バックエンド選択（デフォルト: `uma`） |
| `--opt-mode grad\|hess` | TSOPT と IRC 後の端点最適化の fallback。`--opt-mode-post` が優先されます。 |
| `--hessian-calc-mode Analytical\|FiniteDifference` | ML Hessian 計算モード。全 MLIP バックエンドで `Analytical` を利用可能。`--workers > 1` とは併用不可。 |

`mlmm all --mep-mode dmf` は、デフォルトの単一パス `path-opt` と
`--refine-path` で選択する再帰的 `path-search` のどちらにも Direct Max Flux
を適用します。デフォルトは GSM です。

すべてのオプションと YAML スキーマについては [all](all.md) および [YAML リファレンス](yaml-reference.md) を参照してください。

---

## 実行サマリー

出力ディレクトリの `summary.log` と `summary.json` に、実行コマンド、セグメントごとの障壁高、MEP 統計、後処理結果がまとまります。入力検証で早期に終了した場合は、作られないことがあります。

`segments/seg_NN/` には各段階の計算結果が置かれます。段階別 JSON の出力条件は、[出力ディレクトリ構成](output-layout.md)を参照してください。

---

## ヘルプ

`--help` は主要オプション、`--help-advanced` は全オプションを表示します。

```bash
mlmm all --help
mlmm all --help-advanced
```

個別計算については、[コマンド一覧](index.md#cli-サブコマンド)から各ページを参照してください。
