# 用語集

このページでは、mlmm_toolkit ドキュメント内で使われる略語・専門用語を簡潔に説明します。

---

## ML/MM・ONIOM

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **ML/MM** | Machine Learning / Molecular Mechanics | 機械学習ポテンシャルと分子力学を組み合わせたマルチスケール手法。mlmm_toolkit の中核概念。 |
| **ONIOM** | Our own N-layered Integrated molecular Orbital and molecular Mechanics | 異なる計算レベルを多層的に組み合わせる手法。mlmm_toolkit は ONIOM 的な 2 体系（real/model）エネルギー分解を使用。 |
| **QM/MM** | Quantum Mechanics / Molecular Mechanics | 量子化学と分子力学の結合手法。ML/MM は QM 部分を機械学習ポテンシャルで置き換えた変種。 |
| **real system** | -- | ONIOM 分解における全系。parm7 トポロジで記述され、hessian_ff で MM エネルギーを計算。 |
| **model system** | -- | ONIOM 分解における ML 領域（活性部位）。UMA で ML エネルギーを計算。 |
| **リンク水素** | Link Hydrogen | ML 領域と MM 領域の境界で切断された結合をキャップする水素原子。ヤコビアンで力を再分配。 |

---

## 3層システム

| 用語 | B-factor | 説明 |
|------|----------|------|
| **ML 領域（Layer 1）** | 0.0 | UMA で完全計算（エネルギー・力・ヘシアン）される活性部位原子。 |
| **Movable-MM（Layer 2）** | 10.0 | 最適化時に移動可能な MM 原子。 |
| **Frozen（Layer 3）** | 20.0 | 座標固定の MM 原子。計算不参加。 |

---

## 力場・Amber

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **parm7（prmtop）** | Amber Parameter/Topology file | Amber のトポロジファイル。原子タイプ、結合、角度、二面角、VDW パラメータ、部分電荷を含む。 |
| **rst7（inpcrd）** | Amber Restart / Initial Coordinates | Amber の座標ファイル。原子の 3D 座標を格納。 |
| **hessian_ff** | -- | mlmm_toolkit に同梱される C++ ネイティブ拡張の Amber 力場計算エンジン。解析ヘシアンをサポート。 |
| **GAFF2** | General Amber Force Field 2 | 有機小分子向けの汎用 Amber 力場。リガンドのパラメータ化に使用。 |
| **AM1-BCC** | AM1 Bond Charge Corrections | 半経験的 AM1 法に基づく部分電荷割り当て方法。antechamber で使用。 |
| **AmberTools** | -- | Amber のオープンソースツール群。tleap（トポロジ構築）、antechamber（リガンドパラメータ化）、parmchk2（不足パラメータ補完）を含む。 |
| **tleap** | -- | AmberTools のトポロジ構築プログラム。PDB からparm7/rst7 を生成。 |
| **ff19SB** | -- | タンパク質向け Amber 力場（2019 年版）。mlmm_toolkit のデフォルト。 |
| **ff14SB** | -- | タンパク質向け Amber 力場（2014 年版）。`--ff-set ff14SB` で選択可能。 |

---

## 反応経路・最適化

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **MEP** | Minimum Energy Path | 反応物から生成物へ至る最小エネルギー経路（ポテンシャルエネルギー面上の最も低い経路）。 |
| **TS** | Transition State | 反応座標に沿ったエネルギー極大に対応する一次の鞍点。 |
| **IRC** | Intrinsic Reaction Coordinate | TS から反応物側・生成物側へ向かう、質量重み付き最急降下経路。TS の接続検証によく使われます。 |
| **GSM** | Growing String Method | 端点からストリング（画像列）を伸長・最適化して MEP を近似する手法。 |
| **DMF** | Direct Max Flux | 反応座標方向のフラックスを最大化することで MEP を最適化する chain-of-states 手法。 |
| **NEB** | Nudged Elastic Band | 画像間にばね力を導入し、画像間隔を保ちながら経路を最適化する chain-of-states 手法。 |
| **HEI** | Highest-Energy Image | MEP 上でエネルギーが最大の画像。TS の初期推定としてよく使われます。 |
| **画像（Image）** | -- | 経路上の1つの幾何（1ノード）。 |
| **セグメント** | -- | 2つの隣接する端点を結ぶ MEP（例: R → I1, I1 → I2,...）。 |

---

## 最適化アルゴリズム

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **L-BFGS** | Limited-memory BFGS | 勾配履歴からヘシアンを近似する準ニュートン法。`--opt-mode light` で使用。 |
| **RFO** | Rational Function Optimization | 明示的なヘシアン情報を使用する信頼領域最適化法。`--opt-mode heavy` で使用。 |
| **RS-I-RFO** | Restricted-Step Image-RFO | 1つの負固有値方向に沿う、鞍点（TS）最適化用の RFO 変種。 |
| **Dimer** | Dimer Method | 完全なヘシアンを計算せずに最低曲率モードを推定する TS 最適化法。`--opt-mode light` の TSOPT で使用。 |
| **PHVA** | Partial Hessian Vibrational Analysis | アクティブ（非凍結）原子のヘシアンブロックのみを使用した振動数計算。`freq` のデフォルト。 |

---

## 機械学習・計算機

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **MLIP** | Machine Learning Interatomic Potential | 量子化学データから学習し、構造からエネルギー・力を予測する原子間ポテンシャル。 |
| **UMA** | Universal Machine-learning potential for Atoms | Meta が公開している事前学習 MLIP 群。mlmm_toolkit の ML 領域計算バックエンドです。 |
| **解析ヘシアン** | Analytical Hessian | エネルギーの正確な二階微分を計算。高速だが VRAM を多く消費。 |
| **有限差分** | Finite Difference | 微小変位による微分近似。低速だがメモリ効率が良い。 |

---

## 量子化学

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **QM** | Quantum Mechanics | DFT、HF、post-HF などの第一原理電子状態計算。 |
| **DFT** | Density Functional Theory | 電子密度汎関数に基づく電子状態計算法。 |
| **Hessian** | -- | エネルギーの二階微分行列。振動解析や TS 最適化に使用します。 |
| **SP** | Single Point | 固定構造での計算（最適化なし）。高精度エネルギー補正によく使用。 |
| **スピン多重度** | Spin Multiplicity | 2S+1（S は全スピン）。シングレット = 1、ダブレット = 2、トリプレット = 3 など。 |

---

## 構造生物学・ポケット抽出

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **PDB** | Protein Data Bank | タンパク質などの三次元構造を表す標準フォーマット（およびデータベース）。 |
| **XYZ** | -- | 元素記号と直交座標を並べたシンプルなテキスト形式。 |
| **GJF** | Gaussian Job File | Gaussian の入力形式。電荷/多重度と座標の読み取りに利用可能。 |
| **ポケット** | Active-site Pocket | 基質周辺を切り出した部分構造。MEP/TS 探索の計算量削減に使用。「クラスターモデル」とも。 |
| **クラスターモデル** | Cluster Model | ポケットの別名。酵素-基質複合体から計算可能なサイズに切り出した部分系。 |
| **リンク水素** | Link Hydrogen | ポケット抽出時に切断された結合をキャップするために付加する水素原子。 |
| **主鎖** | Backbone | タンパク質の主骨格（N-Ca-C-O 原子）。`--exclude-backbone` で除外可能。 |
| **B-factor** | Temperature Factor | PDB の温度因子カラム。mlmm_toolkit では 3 層の層割り当てをエンコードするために使用。 |

---

## 熱化学

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **ZPE** | Zero-Point Energy | 0 K での振動エネルギー。電子エネルギーへの量子補正。 |
| **ギブズエネルギー** | Free Energy (G) | G = H - TS。熱・エントロピー寄与を含む。 |
| **エンタルピー** | (H) | H = E + PV。定圧での全熱含量。 |
| **エントロピー** | (S) | 無秩序さの尺度。ギブズエネルギーに -TS として寄与。 |

---

## 単位・定数

| 用語 | 説明 |
|------|------|
| **Hartree** | 原子単位系のエネルギー。1 Hartree = 627.5 kcal/mol = 27.21 eV。 |
| **kcal/mol** | 反応エネルギー表現でよく使われる単位。 |
| **kJ/mol** | キロジュール/モル。1 kcal/mol = 4.184 kJ/mol。 |
| **eV** | 電子ボルト。1 eV = 23.06 kcal/mol。 |
| **Bohr** | 原子単位系の長さ。1 Bohr = 0.529 Å。 |
| **Å（オングストローム）** | 10^-10 m。原子間距離の表現でよく使われる長さ単位。 |

---

## CLI 規則

| 用語 | 説明 |
|------|------|
| **真偽値オプション** | `True` または `False`（大文字始まり）を取る CLI フラグ。例: `--tsopt`。 |
| **残基セレクタ** | `'SAM,GPP'`（名前）や `'A:123,B:456'`（チェーン:ID）のような指定方法。 |
| **原子セレクタ** | `'TYR,285,CA'` のように残基名・番号・原子名で特定の原子を指定する方法。 |
| **B-factor 層エンコーディング** | PDB の B-factor カラムを使用して 3 層の層割り当て（0.0, 10.0, 20.0）をエンコードする方式。Hessian 対象 MM 原子は別設定で制御。 |

---

## 関連項目

- [はじめに](getting_started.md) -- インストールと初回実行
- [概念とワークフロー](concepts.md) -- ポケット抽出、ML/MM レイヤー、MEP 探索、後処理の全体像
- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- よくあるエラーと対処法
- [YAML リファレンス](yaml_reference.md) -- 設定ファイルの仕様
- [ML/MM 計算機](mlmm_calc.md) -- 機械学習ポテンシャルの詳細
