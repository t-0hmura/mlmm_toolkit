# 概念とワークフロー

このページでは、`mlmm_toolkit` を使ううえでの **全体像** を説明します。
ML/MM 3 層システム、ONIOM 分解、「ポケット」「テンプレート」「セグメント」「画像（image）」が何を指すのか、そしてトップレベルの `all` が各サブコマンドをどう組み合わせるのかを把握するためのページです。

---

## ワークフローの全体像

一般的なワークフローは以下の通りです。

```text
全系入力 (PDB)
   │
   ├─ (任意) ポケット抽出              [extract]     ← --center/-c を使う場合は PDB が必要
   │        ↓
   │   ポケット/クラスターモデル (PDB)
   │        │
   │        ├─ MM トポロジ構築          [mm-parm]     ← parm7/rst7 の自動生成
   │        │        ↓
   │        │   real.parm7 + real.rst7
   │        │
   │        ├─ (任意) 3層定義           [define-layer] ← B-factor によるレイヤーエンコード
   │        │
   │        ├─ (任意) 段階的スキャン    [scan]         ← 単一構造ワークフロー
   │        │        ↓
   │        │   順序付けられた中間体
   │        │
   │        └─ MEP 探索                 [path-search] または [path-opt]
   │                 ↓                    ← ML/MM (UMA + hessian_ff) 計算
   │            MEP 軌跡 (mep.trj) + エネルギーダイアグラム
   │
   └─ (任意) TS 最適化 + IRC            [tsopt] → [irc]
             └─ (任意) 熱化学           [freq]
             └─ (任意) DFT 一点計算     [dft]
```

各ステージはサブコマンドとして単独実行できます。また `mlmm all` を使うと、複数ステージをまとめて実行できます。

```{important}
遷移状態: HEI や `tsopt` の出力は **TS 候補** として扱い、`freq`（虚数モードが 1 本）と `irc`（両端が意図した極小へ落ちる）で検証してから解釈してください。
```

---

## 重要な概念

### ML/MM 3層システム

`mlmm_toolkit` の中核は、ONIOM 的な ML/MM 結合スキームです。系を以下の 3 層に分割し、各層に異なるレベルの理論を適用します。

| 層 | B-factor | 計算レベル | 説明 |
|----|----------|-----------|------|
| **Layer 1: ML 領域** | 0.0 | UMA（MLIP） | 活性部位。エネルギー・力・ヘシアンすべてを UMA で計算 |
| **Layer 2: Movable-MM** | 10.0 | hessian_ff（MM） | 最適化時に移動可能な MM 原子 |
| **Layer 3: Frozen** | 20.0 | なし | 座標固定。計算不参加 |

B-factor 値は PDB ファイルの温度因子カラム（列 61-66）にエンコードされます。`define-layer` サブコマンドで自動設定できます。Hessian 対象 MM 原子は `hess_cutoff` や `hess_mm_atoms` で別途制御します。

### ONIOM エネルギー分解

全系の ML/MM エネルギーは次のように計算されます:

```
E(ML/MM) = E_MM(real) + E_ML(model) - E_MM(model)
```

ここで:
- **E_MM(real)**: 全系（real system）の MM エネルギー（hessian_ff）
- **E_ML(model)**: ML 領域（model system）の UMA エネルギー
- **E_MM(model)**: ML 領域の MM エネルギー（重複分を差し引くため）

力とヘシアンも同様の ONIOM 分解で結合されます。リンク水素の寄与はヤコビアンを用いて ML 原子と MM 原子に再分配されます。

### hessian_ff（MM エンジン）

`hessian_ff` は Amber 力場パラメータ（parm7）をベースとした C++ ネイティブ拡張の力場計算エンジンです。主な特徴:

- **解析ヘシアン**: 正確な二階微分
- **CPU 実行**: GPU メモリを UMA 推論に専有させるため、MM 計算は CPU で実行
- **Amber 互換**: ff19SB/ff14SB、GAFF2 などの Amber 力場に対応

### parm7 と rst7

- **parm7（prmtop）**: Amber のトポロジファイル。原子タイプ、結合、角度、二面角、VDW パラメータ、部分電荷などを含む
- **rst7（inpcrd）**: Amber の座標/速度ファイル。原子座標を含む

`mm-parm` サブコマンドが PDB から parm7/rst7 を自動生成します。非標準リガンドには GAFF2 + AM1-BCC で自動パラメータ化します。

### 全系とポケット（クラスターモデル）

- **全系**: 元の入力構造（酵素ならタンパク質-基質複合体など）。
- **ポケット / クラスターモデル**: 基質周辺を切り出した小系。MEP/TS 探索を軽量化します。

ポケット抽出は主に以下で制御します。
- `-c/--center`: 基質の指定（残基ID、残基名、または基質のみのPDB）
- `-r/--radius`, `--radius-het2het`, `--include-H2O`, `--exclude-backbone`, `--add-linkH`, `--selected-resn`

### 画像（image）とセグメント

- **画像（image）**: 経路上の1つの幾何（1ノード）。
- **セグメント**: 2つの端点を結ぶ MEP。多構造入力は、隣接する端点ペアごとのセグメントに分解されます。

### テンプレートとファイル変換（`--convert-files`）

`mlmm` は軌跡（例: `mep.trj`, `irc.trj`）を出力します。
PDB テンプレートがある場合、必要に応じて companion ファイルも出力できます。

---

## 代表的な3つの使い方

### 1) 複数構造の MEP 探索（R → ... → P）

すでに反応座標に沿った **2つ以上** の構造がある場合。

例:

```bash
mlmm -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
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
mlmm -i ts_guess.pdb -c 'SAM,GPP' --tsopt True
```

---

## `all` と個別サブコマンドの使い分け

### `mlmm all` を推奨するケース
- extract → mm-parm → MEP → TSOPT/IRC → freq/DFT まで **まとめて** 実行したい
- 出力ディレクトリやログ管理を 1 コマンドに寄せたい

### 個別サブコマンドを推奨するケース
- あるステージだけを **切り出してデバッグ** したい（例: `mm-parm` のみ）
- `--real-parm7` と `--model-pdb` を手動で用意して、カスタムワークフローを組みたい
- `oniom-gaussian` や `oniom-orca` で Gaussian/ORCA ONIOM 入力を生成したい

---

## CLI に関する重要な注意点

```{important}
- ブール値オプションは `True`/`False` を明示して渡します（例: `--tsopt True`）。
- 複数 PDB を与える場合、各ファイルは **同じ原子が同じ順序** で並んでいることが重要です（座標だけが異なる）。
- 酵素の反応機構解析では、水素を含んだ入力 PDB を用意することを強く推奨します。
- ML/MM 計算には parm7 トポロジが必須です。`all` ワークフローでは自動生成されますが、個別サブコマンドでは `--real-parm7` を明示的に指定する必要があります。
```

---

## 次に読むページ

### 入門
- [はじめに](getting-started.md) -- インストールと初回実行
- [トラブルシューティング](troubleshooting.md) -- よくあるエラーと対処法

### 主要サブコマンド
| サブコマンド | 用途 | ドキュメント |
|------------|------|-------------|
| `all` | エンドツーエンドワークフロー | [all.md](all.md) |
| `extract` | ポケット抽出 | [extract.md](extract.md) |
| `mm-parm` | Amber トポロジ構築 | [mm_parm.md](mm_parm.md) |
| `define-layer` | 3 層 ML/MM 領域定義 | [define_layer.md](define_layer.md) |
| `path-search` | 再帰的 MEP 探索 | [path_search.md](path_search.md) |
| `tsopt` | TS 最適化 | [tsopt.md](tsopt.md) |
| `oniom-gaussian` | Gaussian ONIOM 入力生成 | [oniom_export.md](oniom_export.md) |
| `oniom-orca` | ORCA QM/MM 入力生成 | [oniom_export.md](oniom_export.md) |

### リファレンス
- [YAML リファレンス](yaml-reference.md) -- 全オプションの YAML 設定
- [用語集](glossary.md) -- 用語リファレンス
