# はじめに

## 概要

`mlmm_toolkit` は、機械学習原子間ポテンシャル（MLIP）と分子力学（hessian_ff）を ONIOM 的に結合した **ML/MM 法** を用いて、**PDB 構造** から **酵素反応経路** を自動的に構築する Python 製の CLI ツールキットです。デフォルトの MLIP バックエンドは **UMA**（Meta の FAIR-Chem）で、`--backend` オプションにより **ORB**、**MACE**、**AIMNet2** も選択できます。

多くのワークフローにおいて、**1 コマンド** で反応経路の**初期推定（first-pass）**を得ることができます。
```bash
mlmm -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

---
さらに `--tsopt --thermo --dft` を追加すると、**MEP 探索 → TS 最適化 → IRC → 熱化学解析 → DFT 一点計算** までまとめて実行できます。
```bash
mlmm -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --tsopt --thermo --dft
```
---

入力として、(i) 反応順に並べたタンパク質-リガンド複合体の PDB を 2 つ以上（R →... → P）、(ii) `--scan-lists` を指定した 1 つの PDB、または (iii) TS 候補 1 構造 + `--tsopt` を与えると、`mlmm` が次を自動化します。

- ユーザーが指定した基質の周辺から **活性部位ポケット** を抽出し、計算用の **クラスターモデル** を構築
- AmberTools を用いて **Amber トポロジ（parm7/rst7）** を自動生成し、**hessian_ff** の MM エンジンに渡す
- ML 領域を MLIP バックエンド（デフォルト: UMA）で、MM 領域を hessian_ff で計算する **ONIOM 的 ML/MM** エネルギー・力・ヘシアンを構築
- Growing String Method (GSM) や Direct Max Flux (DMF) などの経路最適化手法で **最小エネルギー経路 (MEP)** を探索
- 必要に応じて **遷移状態** を最適化し、**振動解析**・**IRC 計算**・**DFT 一点計算** を実行

```{important}
単一コマンドの TS 結果は「候補」として扱ってください。酵素反応では、endpoint 品質、ポケット定義、拘束、scan ターゲットの調整を伴う反復が一般的です。最終解釈の前に、`freq` と `irc` の両方で TS を必ず検証してください。
```

ML 領域の計算には機械学習原子間ポテンシャル（MLIP）を使用します。デフォルトのバックエンドは UMA（Meta の FAIR-Chem）で、ORB、MACE、AIMNet2 も `--backend` で選択可能です。MM 領域の計算には hessian_ff（Amber 力場ベースの C++ ネイティブ拡張）を用います。全エネルギーは ONIOM 的な減算分解に従います:

```
E_total = E_REAL_low + E_MODEL_high - E_MODEL_low
```

ここで REAL は全系、MODEL は ML 領域、"high" は MLIP バックエンド、"low" は hessian_ff です。

一連の処理は CLI から呼び出せるように統一されており、手作業を最小化して **多段階の酵素反応メカニズム** を組み立てられるように設計しています。

```{important}
- 入力 PDB ファイルには**水素原子**が含まれている必要があります。
- 複数の PDB を提供する場合、**同じ原子が同じ順序**で含まれている必要があります（座標のみ異なる可能性があります）。そうでない場合はエラーが発生します。
- ML/MM 計算には **`--parm`**（全系の Amber トポロジ）と **`--model-pdb`**（ML 領域を定義する PDB）が必要です。`all` ワークフローではこれらが自動生成されます。
```

```{tip}
初めて使う場合は、まず [概念とワークフロー](concepts.md) を参照してください。
症状から切り分ける場合は、まず [典型エラー別レシピ](recipes_common_errors.md) を参照してください。
セットアップや実行中にエラーが発生した場合は [トラブルシューティング](troubleshooting.md) を参照してください。
```

### CLI の慣習

| 慣習 | 例 | 備考 |
|-----|-----|------|
| **残基セレクタ** | `'SAM,GPP'`, `'A:123,B:456'` | 複数値はシェル展開防止のためクォート |
| **電荷マッピング** | `--ligand-charge 'SAM:1,GPP:-3'` | コロン（`:`）またはイコール（`=`）で名前と電荷を区切り、カンマでエントリを区切る |
| **原子セレクタ** | `'TYR,285,CA'` または `'TYR 285 CA'` | 区切り文字: 空白、カンマ、スラッシュ、バッククォート、バックスラッシュ |

### 水素原子付与の推奨ツール

PDB に水素原子がない場合は、mlmm を実行する前に次のいずれかを使ってください。

| ツール | コマンド例 | 備考 |
|--------|------------|------|
| **reduce** (Richardson Lab) | `reduce input.pdb > output.pdb` | 高速、結晶構造に広く使用 |
| **pdb2pqr** | `pdb2pqr --ff=AMBER input.pdb output.pqr` | 水素を追加し部分電荷を割り当て |
| **Open Babel** | `obabel input.pdb -O output.pdb -h` | 汎用ケモインフォマティクスツールキット |
| **mm-parm --add-h** | `mlmm mm-parm -i input.pdb --add-h` | PDBFixer による水素付加（AmberTools 経由） |

複数の PDB 入力で同一の原子順序を確保するには、すべての構造に同じ水素付与ツールを一貫した設定で適用してください。

```{warning}
このソフトウェアはまだ開発中です。自己責任でご使用ください。
```

---

## インストール

`mlmm_toolkit` は、CUDA 対応 GPU を備えた Linux 環境（ローカルワークステーションまたは HPC クラスター）向けに設計されています。特に **PyTorch**、**fairchem-core (UMA)**、**gpu4pyscf-cuda12x** などの依存関係は、動作する CUDA インストールを前提としています。

### 前提条件

mlmm_toolkit は以下のコンポーネントを使用します:

- **MLIP バックエンド**: ML 領域のエネルギー・力・ヘシアン計算。デフォルトは UMA（fairchem-core）。ORB（`pip install mlmm[orb]`）、AIMNet2（`pip install mlmm[aimnet2]`）も利用可能。MACE は e3nn 競合のため別環境が必要
- **hessian_ff**: MM 領域の Amber 力場計算（C++ 拡張のビルドが必要）
- **AmberTools**: `mm-parm` サブコマンドによる parm7/rst7 の自動生成（tleap、antechamber、parmchk2）

詳細は上流プロジェクトを参照してください:
- fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
- Hugging Faceトークンとセキュリティ: <https://huggingface.co/docs/hub/security-tokens>

### クイックスタート

以下は多くの CUDA 12.9 クラスターで動作する最小限のセットアップ例です。

```bash
# 1) CUDA 対応の PyTorch ビルドをインストール
pip install torch --index-url https://download.pytorch.org/whl/cu129

# 2) GitHub から mlmm_toolkit をインストール
pip install git+https://github.com/t-0hmura/mlmm_toolkit.git

# 3) hessian_ff の C++ 拡張をビルド
cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make

# 4) Plotly 図表エクスポート用のヘッドレス Chrome をインストール
plotly_get_chrome -y
```

環境が変わる場合（別ノード/別コンテナ/別 Python・PyTorch）には、その環境で `hessian_ff` を再ビルドしてください。  
多くのクラスターでは、先に Ninja を入れてから再ビルドすると確実です:

```bash
conda install -c conda-forge ninja -y
cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make clean && make
```

最後に、UMA バックエンドを使用する場合は、モデルをダウンロードできるように **Hugging Face Hub** にログインします:

```bash
huggingface-cli login
```

これはマシン/環境ごとに 1 回だけ行う必要があります。ORB/MACE/AIMNet2 バックエンドを使用する場合は不要です。

### AmberTools のインストール

`mm-parm` サブコマンド（parm7/rst7 の自動生成）には AmberTools が必要です。conda での導入が最も簡単です:

```bash
conda install -c conda-forge ambertools -y
```

AmberTools がインストールされていなくても、`--parm` を手動で用意すれば他のサブコマンドは動作します。

### ステップバイステップインストール

環境を段階的に構築する場合:

1. **CUDA をロード（HPC で環境モジュールを使用する場合）**

 ```bash
 module load cuda/12.9
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

4. **cyipopt をインストール（オプション: DMF 法に必要）**

 ```bash
 conda install -c conda-forge cyipopt -y
 ```

5. **適切な CUDA ビルドの PyTorch をインストール**

 ```bash
 pip install torch --index-url https://download.pytorch.org/whl/cu129
 ```

6. **mlmm_toolkit 本体をインストール**

 ```bash
 pip install git+https://github.com/t-0hmura/mlmm_toolkit.git
 ```

7. **hessian_ff の C++ 拡張をビルド**

 ```bash
 cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make
 ```

 環境が変わった場合は、その環境で Ninja を入れて再ビルドしてください:

 ```bash
 conda install -c conda-forge ninja -y
 cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make clean && make
 ```

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
 pip install mlmm[orb]      # ORB バックエンド
 pip install mlmm[aimnet2]  # AIMNet2 バックエンド
 # MACE: pip uninstall fairchem-core && pip install mace-torch（別環境が必要）
 ```

11. **インストールの確認**

 ```bash
 mlmm --version
 ```

 インストールされたバージョンが表示されます（例: `0.x.y`; 正確な出力は git タグによって異なります）。

---

## マルチバックエンドの使用例

デフォルトの MLIP バックエンドは UMA です。`--backend` で代替バックエンドに切り替え、`--embedcharge` で xTB ポイントチャージ埋め込みを有効化できます:

```bash
# ORB バックエンドを使用
mlmm opt -i pocket.pdb --parm real.parm7 --model-pdb ml.pdb -q 0 --backend orb

# MACE バックエンドを使用
mlmm opt -i pocket.pdb --parm real.parm7 --model-pdb ml.pdb -q 0 --backend mace

# xTB ポイントチャージ埋め込みを有効化
mlmm opt -i pocket.pdb --parm real.parm7 --model-pdb ml.pdb -q 0 --embedcharge
```

---

## 推奨クイックスタート導線

- [クイックスタート: `mlmm all`](quickstart_all.md)
- [クイックスタート: `mlmm scan` + `--spec`](quickstart_scan_spec.md)
- [クイックスタート: `mlmm tsopt` -> `mlmm freq`](quickstart_tsopt_freq.md)

---

## コマンドラインの基本

メインのエントリーポイントは `pip` でインストールされる `mlmm` コマンドです。内部的には **Click** ライブラリを使用しており、デフォルトのサブコマンドは `all` です。

つまり:

```bash
mlmm [OPTIONS]...
# は以下と同等
mlmm all [OPTIONS]...
```

`all` ワークフローは、クラスター抽出、MM パラメータ化、MEP 探索、TS 最適化、振動解析、オプションの DFT 一点計算を 1 つのコマンドで連続実行する**オーケストレーター**です。

---

## 典型ワークフロー

`mlmm all` を個別サブコマンドへ分解すると、典型的には次の順で進みます。

```text
1. extract - 完全系 PDB から活性部位ポケットを抽出
2. mm-parm - Amber parm7/rst7 を生成
3. define-layer - 3 層 ML/MM 分割を付与（B-factor エンコード）
4. path-search - 再帰的 MEP 探索（GSM/DMF）
5. tsopt - 遷移状態最適化
6. freq - 振動解析と熱化学
7. dft - DFT 一点計算
```

`all` は上記 1-7 を自動実行します。必要に応じて各ステップを単独で実行してデバッグできます。

---

## メインワークフローモード

### 複数構造 MEP ワークフロー（反応物 → 生成物）

推定反応座標に沿った複数の完全な PDB 構造（例: R → I1 → I2 → P）がすでにある場合に使用します。

**最小例**

```bash
mlmm -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

**詳細例**

```bash
mlmm -i R.pdb I1.pdb I2.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --out-dir ./result_all --tsopt --thermo --dft
```

動作:

- 反応順序で 2 つ以上の**完全系**を受け取る
- 各構造の触媒クラスターモデルを抽出
- AmberTools で parm7/rst7 を自動生成
- ML/MM 計算機を用いて**再帰的 MEP 探索**を実行
- PDB テンプレートが利用可能な場合、クラスターモデル MEP を**完全系**にマージ
- オプションで各セグメントに対して TS 最適化、振動解析、DFT 一点計算を実行

---

### 単一構造 + 段階的スキャン

**1 つの PDB 構造**しかないが、反応に沿ってどの原子間距離が変化するかが分かっている場合に使用します。

```bash
mlmm -i R.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
 --scan-lists '[("TYR 285 CA","MMT 309 C10",2.20)]'
```

---

### 単一構造 TSOPT のみモード

すでに**遷移状態候補**があり、それを最適化して IRC 計算を行いたい場合に使用します。

```bash
mlmm -i TS_CANDIDATE.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --tsopt
```

---

## 重要な CLI オプションと動作

| オプション | 説明 |
|----------|------|
| `-i, --input PATH...` | 入力構造。**2 つ以上の PDB** → MEP 探索; **1 つの PDB + `--scan-lists`** → 段階的スキャン; **1 つの PDB + `--tsopt`** → TSOPT のみ |
| `-c, --center TEXT` | 基質/抽出中心を定義。残基名（`'SAM,GPP'`）、残基ID（`A:123,B:456`）、または PDB パスをサポート |
| `--ligand-charge TEXT` | 電荷情報: マッピング（`'SAM:1,GPP:-3'`）または単一整数 |
| `-q, --charge INT` | ML 領域の総電荷の強制上書き |
| `--parm PATH` | 全系の Amber parm7 トポロジ（`all` では自動生成） |
| `--model-pdb PATH` | ML 領域を定義する PDB ファイル（`all` では自動生成） |
| `--tsopt/--no-tsopt` | TS 最適化と IRC を有効化 |
| `--thermo/--no-thermo` | 振動解析と熱化学を実行 |
| `--dft/--no-dft` | DFT 一点計算を実行 |
| `--out-dir PATH` | トップレベル出力ディレクトリ |
| `--backend uma\|orb\|mace\|aimnet2` | MLIP バックエンド選択（デフォルト: `uma`） |
| `--embedcharge/--no-embedcharge` | xTB ポイントチャージ埋め込み補正（デフォルト: 無効） |
| `--opt-mode grad\|hess` | `all` のワークフロープリセット: `grad`（LBFGS/Dimer）または `hess`（RFO/RS-I-RFO、デフォルト） |
| `--mep-mode gsm\|dmf` | MEP 手法: Growing String Method または Direct Max Flux |

すべてのオプションと YAML スキーマについては [all](all.md) および [YAML リファレンス](yaml_reference.md) を参照してください。

---

## 実行サマリー

`mlmm all` 実行後、トップレベル出力には次が保存されます。

- `summary.log` -- 人が読むための実行要約
- `summary.yaml` -- 機械処理向けの要約

代表的には、実行コマンド、セグメントごとの障壁高、MEP 統計、後処理（thermo/DFT）結果がまとまります。`path_search/` 以下の各セグメントにも個別 summary が出力されます。

---

## CLI サブコマンド

`mlmm all --help` は主要オプションのみを表示します。`mlmm all --help-advanced` で全オプションを表示できます。
`scan` / `scan2d` / `scan3d` と計算系サブコマンド（`opt` / `path-opt` / `path-search` / `tsopt` / `freq` / `irc` / `dft`）に加え、ユーティリティ系（`mm-parm` / `define-layer` / `add-elem-info` / `trj2fig` / `energy-diagram` / `oniom-export`）も同様に `--help` は主要オプションのみ、`--help-advanced` で全オプションを表示します。`extract` と `fix-altloc` も段階的 help に対応し、`--help-advanced` で parser の全オプションを表示します。

| サブコマンド | 役割 | ドキュメント |
|------------|------|------------|
| `all` | エンドツーエンドワークフロー | [all](all.md) |
| `init` | `mlmm all` 用 YAML テンプレート生成 | [init](init.md) |
| `extract` | 活性部位ポケット抽出 | [extract](extract.md) |
| `mm-parm` | Amber parm7/rst7 構築 | [mm_parm](mm_parm.md) |
| `define-layer` | 3 層 ML/MM 領域定義 | [define_layer](define_layer.md) |
| `opt` | 構造最適化 | [opt](opt.md) |
| `tsopt` | 遷移状態最適化 | [tsopt](tsopt.md) |
| `path-opt` | MEP 最適化 (GSM/DMF) | [path_opt](path_opt.md) |
| `path-search` | 再帰的 MEP 探索 | [path_search](path_search.md) |
| `scan` | 1D 結合長スキャン | [scan](scan.md) |
| `scan2d` | 2D 距離スキャン | [scan2d](scan2d.md) |
| `scan3d` | 3D 距離スキャン | [scan3d](scan3d.md) |
| `irc` | IRC 計算 | [irc](irc.md) |
| `freq` | 振動解析 | [freq](freq.md) |
| `dft` | DFT 一点計算 | [dft](dft.md) |
| `trj2fig` | エネルギープロファイルプロット | [trj2fig](trj2fig.md) |
| `energy-diagram` | 数値系列から状態エネルギー図を描画 | [energy-diagram](energy_diagram.md) |
| `oniom-export` | Gaussian ONIOM / ORCA QM/MM 入力生成（`--mode g16|orca`） | [oniom_export](oniom_export.md) |

---

## クイックリファレンス

よく使う実行パターン:

```bash
# 2 構造以上で基本 MEP 探索
mlmm -i R.pdb P.pdb -c 'SUBSTRATE' --ligand-charge 'SUB:-1'

# TS/熱化学/DFT まで実行
mlmm -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' --tsopt --thermo --dft

# 1 構造 + staged scan
mlmm -i SINGLE.pdb -c 'LIG' --scan-lists '[("RES1,100,CA","LIG,200,C1",2.0)]'

# TS 候補の単独最適化
mlmm -i TS.pdb -c 'LIG' --tsopt --thermo
```

主要オプション:

| オプション | 用途 |
|----------|------|
| `-i` | 入力構造（単数または複数） |
| `-c` | 抽出中心（基質）指定 |
| `--ligand-charge` | 基質電荷指定（例: `'SAM:1,GPP:-3'`） |
| `--parm` | Amber parm7（個別サブコマンドで必要） |
| `--model-pdb` | ML 領域定義 PDB（個別サブコマンドで必要） |
| `--backend` | MLIP バックエンド選択（`uma`, `orb`, `mace`, `aimnet2`） |
| `--embedcharge` | xTB ポイントチャージ埋め込み補正を有効化 |
| `--tsopt` | TS 最適化 + IRC |
| `--thermo` | 振動解析/熱化学 |
| `--dft` | DFT 一点計算 |
| `--out-dir` | 出力ディレクトリ |

---

## ヘルプ

任意のサブコマンドについて:

```bash
mlmm <subcommand> --help
mlmm <subcommand> --help-advanced
mlmm all --help-advanced
```

`all` では `--help` は短縮版です。全オプションを確認するときは `--help-advanced` を使用してください。
