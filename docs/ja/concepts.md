# 概念とワークフロー

このページでは、`mlmm_toolkit` を使ううえでの **全体像** を説明します。
ML/MM 3 層システム、ONIOM 分解、「ポケット」「テンプレート」「セグメント」「画像（image）」が何を指すのか、そしてトップレベルの `all` が各サブコマンドをどう組み合わせるのかを把握するためのページです。

---

## ワークフローの全体像

一般的なワークフローは以下の通りです。

```text
全系入力 (PDB/XYZ/GJF)
 │
 ├─ (任意) ポケット抽出 [extract] ← --center/-c を使う場合は PDB が必要
 │ ↓
 │ ポケット/クラスターモデル (PDB)
 │ │
 │ ├─ MM トポロジー構築 [mm-parm] ← parm7/rst7 の自動生成
 │ │ ↓
 │ │ real.parm7 + real.rst7
 │ │
 │ ├─ (任意) 3層定義 [define-layer] ← B-factor によるレイヤーエンコード
 │ │
 │ ├─ (任意) 段階的スキャン [scan] ← 単一構造ワークフロー
 │ │ ↓
 │ │ 順序付けられた中間体
 │ │
 │ └─ MEP 探索 [path-search] または [path-opt]
 │ ↓
 │ MEP 軌跡 (mep_trj.xyz) + エネルギーダイアグラム  ← ML/MM (MLIP + hessian_ff) 計算
 │
 └─ (任意) TS 最適化 + IRC [tsopt] → [irc]
 └─ (任意) 熱化学 [freq]
 └─ (任意) DFT 一点計算 [dft]
```

各ステージはサブコマンドとして単独実行できます。また `mlmm all` を使うと、複数ステージをまとめて実行できます。

```{important}
遷移状態: HEI や `tsopt` の出力は **TS 候補** として扱い、`freq`（虚振動数モードが 1 本）と `irc`（両端が意図した極小へ落ちる）で検証してから解釈してください。
```

---

## ML/MM 3層システム

`mlmm_toolkit` の中核は、ONIOM 的な ML/MM 結合スキームです。系を以下の 3 層に分割し、各層に異なるレベルの理論を適用します。

| 層 | B-factor | 計算レベル | 説明 |
|----|----------|-----------|------|
| **Layer 1: ML 領域** | 0.0 | MLIP（UMA, ORB, MACE, AIMNet2） | 活性部位。エネルギー・力・ヘシアンすべてを MLIP バックエンドで計算 |
| **Layer 2: Movable-MM** | 10.0 | hessian_ff（MM） | 最適化時に移動可能な MM 原子 |
| **Layer 3: Frozen** | 20.0 | なし | 座標固定。計算不参加 |

B-factor 値は PDB ファイルの温度因子カラム（列 61-66）にエンコードされます。`define-layer` サブコマンドで自動設定できます。Hessian 対象 MM 原子は `hess_cutoff` や `hess_mm_atoms` で別途制御します。

```{tip}
B-factor エンコーディングにより、B-factor カラーリングに対応した分子ビューアで層割り当てを視覚的に確認できます。
```

---

## ONIOM 的エネルギー分解

全系の ML/MM エネルギーは次のように計算されます:

```
E(ML/MM) = E_MM(real) + E_ML(model) - E_MM(model)
```

ここで:
- **E_MM(real)**: 全系（real system）の MM エネルギー（hessian_ff）
- **E_ML(model)**: ML 領域（model system）の MLIP エネルギー（デフォルトは UMA、ORB/MACE/AIMNet2 も利用可能）
- **E_MM(model)**: ML 領域の MM エネルギー（重複分を差し引くため）

力とヘシアンも同様の ONIOM 分解で結合されます。リンク水素の寄与はヤコビアンを用いて ML 原子と MM 原子に再分配されます。

MLIP バックエンドは `-b/--backend`（デフォルト: `uma`）で選択します。代替バックエンド（`orb`、`mace`、`aimnet2`）はオプション依存としてインストールします（例: `pip install "mlmm-toolkit[orb]"`）。

`--embedcharge` を有効にすると、MM 環境が ML 領域に及ぼす静電影響を考慮するための xTB ポイントチャージ埋め込み補正が適用されます。

### 従来の QM/MM との比較

| 側面 | 従来の QM/MM | mlmm_toolkit ML/MM |
|------|-------------|---------------------|
| 高レベル手法 | DFT、HF、post-HF | MLIP（UMA, ORB, MACE, AIMNet2） |
| 低レベル手法 | OpenMM / Amber | hessian_ff（C++ ネイティブ拡張） |
| リンク原子 | 通常必要 | 共有結合境界で自動生成 |
| 静電埋め込み | 一般的 | デフォルトでは不使用（ONIOM 減算による機械的埋め込み）。`--embedcharge` で xTB ポイントチャージ埋め込みを有効化可能 |
| 速度 | 低速（QM がボトルネック） | 高速（GPU 上の ML 推論） |

---

## hessian_ff: MM エンジン

`hessian_ff` は Amber 力場パラメータ（parm7）をベースとした C++ ネイティブ拡張の力場計算エンジンです。主な特徴:

- **解析ヘシアン**: 正確な二階微分
- **CPU 実行**: GPU メモリを MLIP 推論に専有させるため、MM 計算は CPU で実行
- **Amber 互換**: ff19SB/ff14SB、GAFF2 などの Amber 力場に対応

OpenMM とは異なり、`hessian_ff` は ONIOM 結合と振動解析に必要な **MM ヘシアン** を提供するために特化しています。

---

## Amber parm7/rst7 トポロジー

- **parm7（prmtop）**: Amber のトポロジーファイル。原子タイプ、結合、角度、二面角、VDW パラメータ、部分電荷などを含む
- **rst7（inpcrd）**: Amber の座標/速度ファイル。原子座標を含む

`mm-parm` サブコマンドが **AmberTools**（tleap、antechamber、parmchk2）を使用して PDB から parm7/rst7 を自動生成します。具体的には:

- 非標準残基（基質、補因子）の自動識別
- **GAFF2**（General Amber Force Field 2）によるパラメータ化
- AM1-BCC 部分電荷の割り当て
- タンパク質残基には ff19SB で完全なトポロジーを構築

---

## 主要オブジェクトと用語

### 全系とポケット（クラスターモデル）

- **全系**: 元の入力構造（酵素ならタンパク質-基質複合体など）。
- **ポケット / クラスターモデル**: 基質周辺を切り出した小系。MEP/TS 探索を軽量化します。

ポケット抽出は主に以下で制御します。
- `-c/--center`: 基質の指定（残基ID、残基名、または基質のみのPDB）
- `-r/--radius`, `--radius-het2het`, `--include-H2O`, `--exclude-backbone`, `--add-linkH`, `--selected-resn`

### Real system と Model system（ONIOM 用語）

- **Real system（実系）**: すべての原子（3 層すべて）。MM（低レベル）で評価されます。
- **Model system（モデル系）**: ML 領域（Layer 1 のみ）。MLIP（高レベル）と MM（低レベル）の両方で評価されます。

### 画像（image）とセグメント

- **画像（image）**: 経路上の 1 つの構造（1 ノード）。
- **セグメント**: 2つの端点を結ぶ MEP。多構造入力は、隣接する端点ペアごとのセグメントに分解されます。

### テンプレートとファイル変換（`--convert-files`）

`mlmm` は軌跡（例: `mep_trj.xyz`, `irc_trj.xyz`）を出力します。トポロジー対応の入力（PDB テンプレートまたは Gaussian 入力）がある場合、コンパニオンファイルも出力できます:
- PDB テンプレートがあれば `.pdb` コンパニオン
- Gaussian テンプレートがあれば `.gjf` コンパニオン

この動作は `--convert-files {True|False}`（デフォルト: `True`）でグローバルに制御します。

---

## 代表的な3つの使い方

### 1) 複数構造の MEP 探索（R →... → P）

すでに反応座標に沿った **2つ以上** の構造がある場合。

例:

```bash
mlmm -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
```

### 2) 単一構造の段階的スキャン → MEP

入力が **1つ** しかないが、スキャンにより端点系列を生成できる場合。

例:

```bash
mlmm -i holo.pdb -c '308,309' \
 --scan-lists '[("TYR,285,CA","MMT,309,C10",2.20)]'
```

### 3) TSOPT のみ（ポケット TS 最適化）

TS 候補がすでにある、あるいは 1 構造で TS 最適化だけ試したい場合。

例:

```bash
mlmm -i ts_guess.pdb -c 'SAM,GPP' --tsopt
```

---

## `all` と個別サブコマンドの使い分け

### `mlmm all` を推奨するケース
- extract → mm-parm → MEP → TSOPT/IRC → freq/DFT まで **まとめて** 実行したい
- 出力ディレクトリやログ管理を 1 コマンドに寄せたい

### 個別サブコマンドを推奨するケース
- あるステージだけを **切り出してデバッグ** したい（例: `mm-parm` のみ）
- `--parm` と `--model-pdb` を手動で用意して、カスタムワークフローを組みたい
- `oniom-export --mode g16|orca` で Gaussian/ORCA ONIOM 入力を生成したい

---

## CLI に関する重要な注意点

```{important}
- ブール値オプションは `--flag` / `--no-flag` と `--flag True/False`（`yes/no`, `1/0` 含む）の両方を受理します。新規スクリプトでは toggle 形式を推奨します。
- 複数 PDB を与える場合、各ファイルは **同じ原子が同じ順序** で並んでいることが重要です（座標だけが異なる）。
- 酵素の反応機構解析では、水素を含んだ入力 PDB を用意することを強く推奨します。
- ML/MM 計算には parm7 トポロジーが必須です。`all` ワークフローでは自動生成されますが、個別サブコマンドでは `--parm` を明示的に指定する必要があります。
```

---

## 次に読むページ

### 入門
- [はじめに](getting_started.md) -- インストールと初回実行
- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- よくあるエラーと対処法

### 主要サブコマンド
| サブコマンド | 用途 | ドキュメント |
|------------|------|-------------|
| `all` | エンドツーエンドワークフロー | [all.md](all.md) |
| `extract` | ポケット抽出 | [extract.md](extract.md) |
| `mm-parm` | Amber トポロジー構築 | [mm_parm.md](mm_parm.md) |
| `define-layer` | 3 層 ML/MM 領域定義 | [define_layer.md](define_layer.md) |
| `path-search` | 再帰的 MEP 探索 | [path_search.md](path_search.md) |
| `tsopt` | TS 最適化 | [tsopt.md](tsopt.md) |
| `freq` | 振動解析 | [freq.md](freq.md) |
| `dft` | DFT 一点計算 | [dft.md](dft.md) |
| `oniom-export` | Gaussian ONIOM / ORCA QM/MM 入力生成（`--mode g16\|orca`） | [oniom_export.md](oniom_export.md) |

### リファレンス
- [YAML リファレンス](yaml_reference.md) -- 全オプションの YAML 設定
- [用語集](glossary.md) -- 用語リファレンス
