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
単一コマンドの TS 結果は「候補」として扱ってください。酵素反応では、endpoint 品質、ポケット定義、拘束、scan ターゲットの調整を伴う反復が一般的です。最終解釈の前に、`freq` と `irc` の両方で TS を必ず検証してください。
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

`path-search` の命名に関する注意: CLI サブコマンドとドキュメントは `path-search`（ハイフン）、内部ワークフローモジュールは `path_search`（アンダースコア）です。

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

`mlmm-toolkit` は、CUDA 対応 GPU を備えた Linux 環境（ローカルワークステーションまたは HPC クラスター）向けに設計されています。特に **PyTorch**、**fairchem-core (UMA)**、**gpu4pyscf-cuda12x** などの依存関係は、動作する CUDA インストールを前提としています。

### 前提条件

mlmm-toolkit は以下のコンポーネントを使用します:

- **MLIP バックエンド**: ML 領域のエネルギー・力・Hessian 計算。デフォルトは UMA（fairchem-core）。ORB（`pip install "mlmm-toolkit[orb]"`）、AIMNet2（`pip install "mlmm-toolkit[aimnet]"`）も利用可能。MACE も利用可能ですが、`e3nn` バージョン競合のため `fairchem-core` を先にアンインストールする必要があります（`pip uninstall fairchem-core && pip install mace-torch`）。
- **hessian_ff**: MM 領域の Amber 力場計算（C++ 拡張のビルドが必要）
- **AmberTools**: `mm-parm` サブコマンドによる parm7/rst7 の自動生成（tleap、antechamber、parmchk2）

詳細は上流プロジェクトを参照してください:
- fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
- Hugging Face トークンとセキュリティ: <https://huggingface.co/docs/hub/security-tokens>

### クイックスタート

以下は多くの CUDA 12.9 クラスターで動作する最小限のセットアップ例です。この例はトポロジーを自動生成するデフォルトの `all` ルートと GSM MEP モード（DMF なし）を想定しています。先に AmberTools をインストールしてください。DMF を使用する場合は `cyipopt` と `pydmf>=1.2` も必要です。

```bash
# 1) AmberTools と CUDA 対応の PyTorch ビルドをインストール
# 2) mlmm-toolkit をインストール
# 3) hessian_ff の C++ 拡張をビルド
# 4) Plotly 図表エクスポート用のヘッドレス Chrome をインストール

conda install -c conda-forge ambertools -y
pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cu129
pip install mlmm-toolkit

# オプション: 代替 MLIP バックエンドのインストール
pip install "mlmm-toolkit[orb]"       # ORB バックエンド
pip install "mlmm-toolkit[aimnet]"   # AIMNet2 バックエンド
# MACE バックエンド (UMA と競合 — 先に fairchem-core をアンインストール)
# pip uninstall fairchem-core && pip install mace-torch

cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make
plotly_get_chrome -y
```

> **Note:** 環境が変わる場合（別ノード/別コンテナ/別 Python・PyTorch）は、その環境で `hessian_ff` を再ビルドしてください。  
> 多くのクラスターでは、先に Ninja を入れてから再ビルドすると確実です:
>
> ```bash
> conda install -c conda-forge ninja -y
> cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make clean && make
> ```

最後に、UMA バックエンドを使用する場合は、モデルをダウンロードできるように **Hugging Face Hub** にログインします（UMA バックエンド使用時のみ必要）:

```bash
# Hugging Face CLI
hf auth login --token '<YOUR_ACCESS_TOKEN>' --add-to-git-credential
```

または

```bash
# クラシック CLI
huggingface-cli login
```

これはマシン/環境ごとに 1 回だけ行う必要があります。

- MEP 探索で Direct Max Flux (DMF) 法を使用する場合は、mlmm のインストール前に conda 環境を作成して `cyipopt` と `pydmf>=1.2` をインストールしてください。
  ```bash
  # 専用の conda 環境を作成してアクティブ化
  conda create -n mlmm python=3.11 -y
  conda activate mlmm

  # cyipopt と pydmf をインストール（MEP 探索の DMF 法に必要）
  conda install -c conda-forge cyipopt -y
  pip install 'pydmf>=1.2'
  ```

- 公式 PyTorch wheel には CUDA のユーザー空間ライブラリが含まれます。通常は互換性のある NVIDIA ドライバーと割り当て済み GPU だけで動作し、ローカル CUDA toolkit や CUDA モジュールは不要です。C/CUDA 拡張をソースからビルドする場合だけ、サイトが指定する toolkit/compiler モジュールをビルド時と実行時の両方で使用してください。

### ステップバイステップインストール

環境を段階的に構築する場合:

1. **NVIDIA ドライバーと GPU 割り当てを確認**

    ```bash
    nvidia-smi
    ```

2. **conda 環境を作成してアクティブ化**

    ```bash
    conda create -n mlmm python=3.11 -y
    conda activate mlmm
    ```

3. **AmberTools をインストール**

    ```bash
    conda install -c conda-forge ambertools -y
    ```

4. **cyipopt と pydmf をインストール（オプション: DMF 法に必要）**

    ```bash
    conda install -c conda-forge cyipopt -y
    pip install 'pydmf>=1.2'
    ```

5. **適切な CUDA ビルドの PyTorch をインストール**

    ```bash
    pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cu129
    ```

6. **mlmm 本体をインストール**

    ```bash
    pip install mlmm-toolkit
    ```

7. **hessian_ff の C++ 拡張をビルド**

    多くの環境では初回使用時に JIT コンパイルされます。ネイティブ拡張が利用できない旨の警告が表示された場合は、手動でビルドしてください:

    ```bash
    cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make
    ```

    > **Note:** 環境が変わった場合は、その環境で Ninja を入れて再ビルドしてください:
    >
    > ```bash
    > conda install -c conda-forge ninja -y
    > cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make clean && make
    > ```

8. **Plotly 可視化用 Chrome をインストール**

    ```bash
    plotly_get_chrome -y
    ```

9. **Hugging Face Hub にログイン（UMA バックエンド使用時のみ必要）**

    ```bash
    huggingface-cli login
    ```

10. **（任意）代替 MLIP バックエンドのインストール**

    ```bash
    pip install "mlmm-toolkit[orb]"      # ORB バックエンド
    pip install "mlmm-toolkit[aimnet]"  # AIMNet2 バックエンド
    # MACE バックエンド (UMA と競合 — 先に fairchem-core をアンインストール)
    # pip uninstall fairchem-core && pip install mace-torch
    ```

11. **インストールの確認**

    ```bash
    mlmm --version
    ```

    インストールされたバージョンが表示されます（例: `0.x.y`; 正確な出力は git タグによって異なります）。

---

## マルチバックエンドの使用例

デフォルトの MLIP バックエンドは UMA です。`-b/--backend` で代替バックエンドに切り替えます:

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

メインのエントリーポイントは `pip` でインストールされる `mlmm` コマンドです。内部的には **Click** ライブラリを使用しており、デフォルトのサブコマンドは `all` です。

つまり:

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

### 複数構造 MEP ワークフロー（反応物 → 生成物）

推定反応座標に沿った複数の完全な PDB 構造（例: R → I1 → I2 → P）がすでにある場合に使用します。

**最小例**

```bash
mlmm -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
```

**詳細例**

```bash
mlmm -i R.pdb I1.pdb I2.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --out-dir ./result_all --tsopt --thermo --dft
```

動作:

- 反応順序で 2 つ以上の**完全系**を受け取る
- 各構造の ML 領域を抽出・定義
- Amber parm7/rst7 トポロジーを生成し、3 層 ML/MM 分割を付与
- デフォルトで単一パス `path-opt` による **MEP 探索**を実行（隣接ペアごとに GSM、出力は `_work/path_opt/` 以下）
- `--refine-path` を指定すると、再帰的 `path-search`（自動精密化）に切り替え（出力は `_work/path_search/` 以下）
- PDB テンプレートが利用可能な場合、ML 領域 MEP を**完全系**にマージ
- オプションで各セグメントに対して TS 最適化、振動解析、DFT 一点計算を実行

このモードは、適度に間隔を空けた中間体（例: ドッキング、MD、手動モデリングから）を生成できる場合に推奨されます。

```{important}
`mlmm-toolkit` は複数の入力 PDB が**同じ原子を同じ順序で**含むことを前提とします（座標のみ異なる）。座標以外のフィールドが入力間で異なる場合はエラーが発生します。入力 PDB ファイルには**水素原子**も含まれている必要があります。
```

---

### 単一構造 + 段階的スキャン（MEP 精密化に供給）

**1 つの PDB 構造**しかないが、反応に沿ってどの原子間距離が変化するかが分かっている場合に使用します。

`-i` に 1 つの構造を指定し、`--scan-lists` を併用します:

**最小例**

```bash
mlmm -i R.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --scan-lists '[("TYR 285 CA","MMT 309 C10",2.20),("TYR 285 CB","MMT 309 C11",1.80)]' '[("TYR 285 CB","MMT 309 C11",1.20)]'
```

**詳細例**

```bash
mlmm -i SINGLE.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --scan-lists '[("TYR 285 CA","MMT 309 C10",2.20),("TYR 285 CB","MMT 309 C11",1.80)]' '[("TYR 285 CB","MMT 309 C11",1.20)]' --multiplicity 1 --out-dir ./result_scan_all --tsopt --thermo --dft
```

要点:

- `--scan-lists` は抽出された ML 領域上での**段階的距離スキャン**を定義します。
- 各タプル `(i, j, target_A)` は:
 - `'TYR,285,CA'` のような PDB 原子セレクタ文字列（**区切り文字: 空白/カンマ/スラッシュ/バッククォート/バックスラッシュ**）**または** 1-based の原子インデックス
 - ML 領域のインデックスに自動的にリマッピングされます。
- 1 つの `--scan-lists` リテラルで単一スキャンステージ、複数リテラルで逐次ステージを実行。複数リテラルは 1 つのフラグの後に続けて指定します（フラグの繰り返しは不可）。
- 各ステージは `stage_XX/result.pdb` を出力し、中間体または生成物の候補として扱われます。
- デフォルトの `all` ワークフローは連結されたステージに対して単一パス `path-opt` GSM チェーンを実行します。
- `--refine-path` を使用すると、再帰的 `path-search`（自動精密化）に切り替わります。

このモードは、単一構造から反応経路を構築する場合に有用です。

---

### 単一構造 TSOPT のみモード

すでに**遷移状態候補**があり、それを最適化して IRC 計算を行いたい場合に使用します。

PDB を 1 つだけ指定し、`--tsopt` を有効にします:

**最小例**

```bash
mlmm -i TS_CANDIDATE.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt
```

**詳細例**

```bash
mlmm -i TS_CANDIDATE.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt --thermo --dft --out-dir ./result_tsopt_only
```

動作:

- MEP/経路探索を完全にスキップ
- ML 領域の **TS** を TS 最適化で最適化
- 両方向に **IRC** を実行し、未割当 endpoint E1/E2 を極小化
- その後 `freq` と `dft` を E1/TS/E2 に対して実行可能。構造を確認してから反応物/生成物を割り当てる
- MLIP、Gibbs、DFT//MLIP/MM エネルギー図を生成

```{important}
単一入力での実行には、**`--scan-lists`**（段階的スキャン → GSM）**または** **`--tsopt`**（TSOPT のみ）のいずれかが必要です。単一の `-i` のみでこれらを指定しないと、完全なワークフローはトリガーされません。
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
| `-s, --scan-lists TEXT...` | 単一入力実行時の段階的距離スキャン（YAML/JSON ファイルまたはインラインリテラル） |
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

`mlmm all` 実行後、トップレベル出力には次が保存されます。

- `summary.log` — 人が読むための実行要約
- `summary.json` — 機械処理向けの要約

通常は、実行コマンド、セグメントごとの障壁高、MEP 統計、後処理（thermo/DFT）結果がまとまります。セグメント別の記録はルートの `summary.json` に集約されます。`segments/seg_NN/` には正規の reactant/TS/product 構造と、実行された各段階のディレクトリが置かれます。各段階の `result.json`/`summary.json` は、その段階の writer が JSON を出力した場合だけ存在します。詳細は [出力ディレクトリ構成](output-layout.md) を参照してください。

---

## CLI サブコマンド

ほとんどのユーザーは主に `mlmm all` を使用します。CLI は個別のサブコマンドも公開しており、各サブコマンドは `-h/--help` に対応しています。
`mlmm all --help` は主要オプションのみを表示します。`mlmm all --help-advanced` で全オプションを表示できます。
`scan` / `scan2d` / `scan3d` と計算系サブコマンド（`opt` / `path-opt` / `path-search` / `tsopt` / `freq` / `irc` / `dft`）に加え、ユーティリティ系（`mm-parm` / `define-layer` / `add-elem-info` / `trj2fig` / `energy-diagram` / `oniom-export`）も同様に `--help` は主要オプションのみ、`--help-advanced` で全オプションを表示します。`extract` と `fix-altloc` も段階的 help に対応し、`--help-advanced` で parser の全オプションを表示します。

| サブコマンド | 役割 | ドキュメント |
|------------|------|------------|
| `all` | 一気通貫ワークフロー | [all](all.md) |
| `extract` | 活性部位ポケット抽出 | [extract](extract.md) |
| `mm-parm` | Amber parm7/rst7 構築 | [mm-parm](mm-parm.md) |
| `define-layer` | 3 層 ML/MM 領域定義 | [define-layer](define-layer.md) |
| `opt` | 構造最適化 | [opt](opt.md) |
| `tsopt` | 遷移状態最適化 | [tsopt](tsopt.md) |
| `path-opt` | MEP 最適化 (GSM/DMF) | [path-opt](path-opt.md) |
| `path-search` | 再帰的 MEP 探索 | [path-search](path-search.md) |
| `scan` | 1D 結合長スキャン | [scan](scan.md) |
| `scan2d` | 2D 距離スキャン | [scan2d](scan2d.md) |
| `scan3d` | 3D 距離スキャン | [scan3d](scan3d.md) |
| `irc` | IRC 計算 | [irc](irc.md) |
| `freq` | 振動解析 | [freq](freq.md) |
| `dft` | DFT 一点計算 | [dft](dft.md) |
| `oniom-export` | Gaussian ONIOM / ORCA QM/MM 入力生成（`--mode g16\|orca`） | [oniom-export](oniom-export.md) |
| `oniom-import` | Gaussian/ORCA ONIOM 入力から XYZ + 層付き PDB を再構築 | [oniom-import](oniom-import.md) |
| `trj2fig` | エネルギープロファイルプロット | [trj2fig](trj2fig.md) |
| `energy-diagram` | 数値系列から状態エネルギー図を描画 | [energy-diagram](energy-diagram.md) |
| `add-elem-info` | PDB の元素列（77-78）を修復 | [add-elem-info](add-elem-info.md) |
| `fix-altloc` | PDB の代替位置標識（altLoc）を除去 | [fix-altloc](fix-altloc.md) |

```{tip}
`all`、`tsopt`、`freq`、`irc` では、VRAM に余裕がある場合 **`--hessian-calc-mode Analytical`**（ML 領域用）を使用できます。UMA、ORB、MACE、AIMNet2 が対応しますが、`--workers > 1` と同時に指定するとエラーになります。
```

---

## クイックリファレンス

よく使う実行パターン:

```bash
# 2 構造以上で基本 MEP 探索
mlmm -i R.pdb P.pdb -c 'SUBSTRATE' -l 'SUB:-1'

# TS/熱化学/DFT まで実行
mlmm -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt --thermo --dft

# 1 構造 + staged scan
mlmm -i SINGLE.pdb -c 'LIG' -l 'LIG:-1' --scan-lists '[("RES1,100,CA","LIG,200,C1",2.0)]'

# TS 候補の単独最適化
mlmm -i TS.pdb -c 'LIG' -l 'LIG:-1' --tsopt --thermo

# 個別サブコマンド（extract + mm-parm + define-layer 実行後）
mlmm path-search -i R.pdb P.pdb --parm real.parm7 --model-pdb model.pdb -q 0 -m 1
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb model.pdb -q 0 -m 1
```

主要オプション:

| オプション | 用途 |
|----------|------|
| `-i` | 入力構造（単数または複数） |
| `-c` | 抽出中心（基質）指定 |
| `-l, --ligand-charge` | 基質電荷指定（例: `'SAM:1,GPP:-3'`） |
| `--parm` | Amber parm7（個別サブコマンドで必要） |
| `--model-pdb` | ML 領域定義 PDB（個別サブコマンドでは `--model-indices` または有効な B-factor layer も選択可能） |
| `-b, --backend` | MLIP バックエンド選択（`uma`, `orb`, `mace`, `aimnet2`） |
| `--tsopt` | TS 最適化 + IRC |
| `--thermo` | 振動解析/熱化学 |
| `--dft` | DFT 一点計算 |
| `-o, --out-dir` | 出力ディレクトリ |

---

## ヘルプ

任意のサブコマンドについて:

```bash
mlmm <subcommand> --help
mlmm <subcommand> --help-advanced
mlmm all --help-advanced
```

`all` では `--help` は短縮版です。全オプションを確認するときは `--help-advanced` を使用してください。
