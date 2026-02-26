# mlmm_toolkit ドキュメント

**バージョン: {{ version }}**

**mlmm_toolkit** は、機械学習原子間ポテンシャル（UMA）と分子力学（hessian_ff）を ONIOM 的に結合した **ML/MM 法** を用いて、PDB 構造から酵素反応経路を自動モデリングする Python 製 CLI ツールキットです。

```{toctree}
:maxdepth: 2
:caption: ガイド
:hidden:

getting_started
quickstart_all
quickstart_scan_spec
quickstart_tsopt_freq
concepts
recipes_common_errors
troubleshooting
cli_conventions
```

```{toctree}
:maxdepth: 2
:caption: コマンド
:hidden:

all
init
extract
add_elem_info
mm_parm
define_layer
opt
tsopt
path_opt
path_search
scan
scan2d
scan3d
freq
irc
dft
trj2fig
oniom_export
oniom_import
fix_altloc
energy_diagram
```

```{toctree}
:maxdepth: 2
:caption: リファレンス
:hidden:

yaml_reference
mlmm_calc
glossary
```


---

## 目的別クイックスタート

| ユースケース | 推奨コマンド | ガイド |
|--------------|--------------|--------|
| 最初の 1 回を実行（エンドツーエンド） | `mlmm all` | [クイックスタート: all](quickstart_all.md) |
| 単一構造スキャン（`--spec`） | `mlmm scan` | [クイックスタート: scan + spec](quickstart_scan_spec.md) |
| TS 検証（`tsopt` -> `freq`） | `mlmm tsopt`, `mlmm freq` | [クイックスタート: tsopt -> freq](quickstart_tsopt_freq.md) |
| PDB から反応経路探索を一通り実行 | `mlmm all` | [all.md](all.md) |
| `all` 用 YAML テンプレートを生成 | `mlmm init` | [init.md](init.md) |
| タンパク質-リガンド複合体からQM領域を抽出 | `mlmm extract` | [extract.md](extract.md) |
| MM トポロジ（parm7/rst7）を構築 | `mlmm mm-parm` | [mm_parm.md](mm_parm.md) |
| ML/MM 3層領域を定義 | `mlmm define-layer` | [define_layer.md](define_layer.md) |
| 単一構造を最適化 | `mlmm opt` | [opt.md](opt.md) |
| 遷移状態を探索・最適化 | `mlmm tsopt` | [tsopt.md](tsopt.md) |
| 最小エネルギー経路を探索 | `mlmm path-search` | [path_search.md](path_search.md) |
| 遷移状態からIRCを実行 | `mlmm irc` | [irc.md](irc.md) |
| エネルギープロファイルを可視化 | `mlmm trj2fig` | [trj2fig.md](trj2fig.md) |
| Gaussian ONIOM / ORCA QM/MM 入力を生成 | `mlmm oniom-export --mode <g16|orca>` | [oniom_export.md](oniom_export.md) |
| ONIOM 入力から XYZ/層付き PDB を再構築 | `mlmm oniom-import` | [oniom_import.md](oniom_import.md) |
| 数値から状態エネルギーダイアグラムを描画 | `mlmm energy-diagram` | [energy_diagram.md](energy_diagram.md) |
| チュートリアルに従う | -- | [はじめに](getting_started.md) |
| 症状からエラー対処を探す | -- | [典型エラー別レシピ](recipes_common_errors.md) |
| 全体像（概念・用語）の把握 | -- | [概念とワークフロー](concepts.md) |
| よくあるエラーの解決 | -- | [トラブルシューティング](troubleshooting.md) |
| 略語や用語を調べる | -- | [用語集](glossary.md) |

---

## ドキュメントガイド

| トピック | ページ |
|---------|--------|
| **インストールと初回実行** | [はじめに](getting_started.md) |
| **主要概念とワークフロー概要** | [概念とワークフロー](concepts.md) |
| **症状起点の切り分け** | [典型エラー別レシピ](recipes_common_errors.md) |
| **よくあるエラーと対処** | [トラブルシューティング](troubleshooting.md) |
| **CLI 規約と入力要件** | [CLI 規約](cli_conventions.md) |

---

## CLI サブコマンド

### メインワークフロー
| サブコマンド | 説明 |
|---------|------|
| [`all`](all.md) | エンドツーエンドワークフロー: 抽出 -> MM parm -> MEP -> TS 最適化 -> IRC -> freq -> DFT |
| [`init`](init.md) | `mlmm all` 用 YAML テンプレートを生成 |

### 構造準備
| サブコマンド | 説明 |
|---------|------|
| [`extract`](extract.md) | タンパク質-リガンド複合体から活性部位ポケット（クラスターモデル）を抽出 |
| [`add-elem-info`](add_elem_info.md) | PDB の元素カラム（77-78）を修復 |
| [`mm-parm`](mm_parm.md) | AmberTools (tleap + GAFF2) を使用して Amber トポロジ（parm7/rst7）を構築 |
| [`define-layer`](define_layer.md) | ML 領域からの距離に基づき 3 層 ML/MM 領域を定義し、B-factor でエンコード |

### 構造最適化
| サブコマンド | 説明 |
|---------|------|
| [`opt`](opt.md) | 単一構造の構造最適化（L-BFGS / RFO） |
| [`tsopt`](tsopt.md) | 遷移状態最適化（Dimer / RS-I-RFO） |

### 経路探索・最適化
| サブコマンド | 説明 |
|---------|------|
| [`path-opt`](path_opt.md) | GSM または DMF による MEP 最適化 |
| [`path-search`](path_search.md) | 自動精密化を伴う再帰的 MEP 探索 |

### スキャン
| サブコマンド | 説明 |
|---------|------|
| [`scan`](scan.md) | 拘束条件付き 1D 結合長スキャン |
| [`scan2d`](scan2d.md) | 2D 距離グリッドスキャン |
| [`scan3d`](scan3d.md) | 3D 距離グリッドスキャン |

### 解析・後処理
| サブコマンド | 説明 |
|---------|------|
| [`irc`](irc.md) | 固有反応座標（IRC）計算 |
| [`freq`](freq.md) | 振動数解析と熱化学 |
| [`dft`](dft.md) | DFT 一点計算（GPU4PySCF / PySCF） |
| [`trj2fig`](trj2fig.md) | XYZ 軌跡からエネルギープロファイルをプロット |
| [`energy-diagram`](energy_diagram.md) | 数値入力からエネルギーダイアグラムを作成 |

### エクスポート
| サブコマンド | 説明 |
|---------|------|
| [`oniom-export`](oniom_export.md) | Amber parm7 から Gaussian ONIOM / ORCA QM/MM 入力を生成（`--mode g16|orca`） |
| [`oniom-import`](oniom_import.md) | Gaussian/ORCA ONIOM 入力を読み込み、XYZ と層付き PDB を再構築 |

---

## 設定・リファレンス

| トピック | ページ |
|---------|--------|
| **YAML 設定オプション** | [YAML リファレンス](yaml_reference.md) |
| **ML/MM 計算機アーキテクチャ** | [ML/MM 計算機](mlmm_calc.md) |
| **用語集** | [用語集](glossary.md) |

---

## システム要件

### ハードウェア
- **OS**: Linux（Ubuntu 20.04+、CentOS 8+ で動作確認）
- **GPU**: CUDA 12.x 互換
- **VRAM**: 最小 8 GB（1000 原子以上には 16 GB 以上推奨）
- **RAM**: 16 GB 以上推奨

### ソフトウェア
- Python 3.11
- CUDA サポート付き PyTorch
- CUDA 12.x ツールキット
- AmberTools（`mm-parm` サブコマンドに必要）

---

## クイック例

### 基本的な ML/MM MEP 探索
```bash
mlmm -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

### TS 最適化を含む完全ワークフロー
```bash
mlmm -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft
```

### 単一構造スキャンモード
```bash
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --spec scan.yaml --print-parsed
```

### TS 最適化のみ
```bash
mlmm -i TS_candidate.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
 --tsopt
```

---

## 重要な概念

### ML/MM 3 層システム
mlmm_toolkit は PDB の B-factor を用いた 3 層分割スキームを使用します:
- **ML 領域**（B=0.0）: UMA 機械学習ポテンシャルで計算
- **Movable-MM**（B=10.0）: 最適化時に移動可能な MM 原子
- **Frozen**（B=20.0）: 座標固定の MM 原子

Hessian 計算に含める MM 原子は、B-factor 専用層ではなく `hess_cutoff` や `hess_mm_atoms` で制御します。

### 電荷とスピン
- 未知残基の電荷を指定するには `--ligand-charge` を使用: `'SAM:1,GPP:-3'`
- ML 領域の総電荷を上書きするには `-q/--charge` を使用
- スピン多重度は `-m/--multiplicity`（デフォルト: 1）で設定

### ブール値オプション
ブール値 CLI オプションはトグル形式（`--flag` / `--no-flag`）を使用します:
```bash
--tsopt --thermo --no-dft
```

### YAML 設定

すべてのオプションについては [YAML リファレンス](yaml_reference.md) を参照してください。

---

## 出力構造

典型的な `mlmm all` 実行の出力:
```
result_all/
├── ml_region.pdb # ML 領域定義
├── summary.log # 人間が読めるサマリー
├── summary.yaml # 機械可読サマリー
├── pockets/ # 抽出されたクラスターモデル
├── mm_parm/ # AMBER トポロジファイル
├── scan/ # （オプション）スキャン結果
├── path_search/ # MEP 軌跡とダイアグラム
│ ├── mep_trj.xyz # MEP 軌跡
│ ├── mep.pdb # PDB 形式の MEP
│ ├── mep_w_ref.pdb # 全系とマージされた MEP
│ └── seg_*/ # セグメントごとの詳細
└── path_search/post_seg_*/ # 後処理出力
 ├── tsopt/ # TS 最適化結果
 ├── irc/ # IRC 軌跡
 ├── freq/ # 振動モード
 └── dft/ # DFT 結果
```

---

## 引用

本ソフトウェアを研究に使用した場合は、以下を引用してください:

[1] Ohmura, T., Inoue, S., Terada, T. (2025). ML/MM toolkit -- Towards Accelerated Mechanistic Investigation of Enzymatic Reactions. ChemRxiv. https://doi.org/10.26434/chemrxiv-2025-jft1k

## ライセンス

`mlmm_toolkit` は Pysisyphus から派生した **GNU General Public License version 3 (GPL-3.0)** の下で配布されています。

---

## 参考文献

1. Wood, B. M. et al. (2025). UMA: A Family of Universal Models for Atoms. [arXiv:2506.23971](http://arxiv.org/abs/2506.23971)
2. Steinmetzer, J., Kupfer, S., & Grafe, S. (2021). pysisyphus: Exploring potential energy surfaces in ground and excited states. *Int. J. Quantum Chem.*, 121(3). [DOI:10.1002/qua.26390](https://doi.org/10.1002/qua.26390)

---

## ヘルプ

```bash
# 一般的なヘルプ
mlmm --help

# コマンドのヘルプ
mlmm <command> --help
```
