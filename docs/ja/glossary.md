# 用語集

ドキュメントに出てくる略語・手法名・単位を、分野ごとに 1 行で引くページです。オプション、出力の欄、ステータスの値（`stalled` など）は、各コマンドのページと [JSON 出力リファレンス](json-output.md) にあります。

## ML/MM・ONIOM

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **ML/MM** | Machine Learning / Molecular Mechanics | 機械学習ポテンシャルと分子力学を組み合わせたマルチスケール手法。mlmm-toolkit の中核概念。 |
| **ONIOM** | Our own N-layered Integrated molecular Orbital and molecular Mechanics | 異なる計算レベルを多層的に組み合わせる手法。mlmm-toolkit は ONIOM 的な引き算方式を使用: E_total = E_REAL_low + E_MODEL_high - E_MODEL_low。 |
| **QM/MM** | Quantum Mechanics / Molecular Mechanics | 量子化学と分子力学の結合手法。ML/MM は QM 部分を機械学習ポテンシャルで置き換えた変種。 |
| **real system** | — | ONIOM 分解における全系（3 層すべて）。parm7 トポロジーで記述され、MM バックエンドで MM エネルギーを計算。 |
| **model system** | — | ONIOM 分解における ML 領域（Layer 1）。MLIP バックエンド（デフォルト: UMA）と MM の両方で評価。 |
| **リンク水素** | Link Hydrogen | ML/MM 境界を横切る parm7 の結合ごとに置く水素原子です。その結合に沿って配置し、力はヤコビアンで再分配します。`extract --add-linkh` が付ける水素はポケット確認用のキャップで、境界の定義には使いません。 |
| **リンク原子** | Link atom | **リンク水素** を参照。mlmm-toolkit では、切断した ML/MM 境界に置くリンク原子は水素です。 |
| **hessian_ff** | — | mlmm-toolkit に同梱される C++ ネイティブ拡張の Amber 力場計算エンジン。解析 Hessian をサポート。 |
| **3 層システム** | 3-layer system | mlmm-toolkit の B-factor による層分割方式: ML（B=0.0）、Movable-MM（B=10.0）、Frozen-MM（B=20.0）。 |
| **B-factor エンコーディング** | B-factor encoding | PDB の B-factor（温度因子）列に層の所属を格納する方式: 0.0 = ML、10.0 = Movable-MM、20.0 = Frozen-MM。Hessian 対象 MM 原子はカットオフ/明示的インデックスで制御。{ref}`MM の層 <ja-mm-layers>` を参照。 |

## 力場・Amber

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **parm7（prmtop）** | Amber Parameter/Topology file | Amber のトポロジーファイル。原子タイプ、結合、角度、二面角、VDW パラメータ、部分電荷を含む。 |
| **rst7（inpcrd）** | Amber Restart / Initial Coordinates | Amber の座標ファイル。原子の 3D 座標を格納。 |
| **AmberTools** | — | Amber のオープンソースツール群。tleap（トポロジー構築）、antechamber（リガンドパラメータ化）、parmchk2（不足パラメータ補完）を含む。`mlmm mm-parm` に必要。 |
| **tleap** | — | AmberTools のトポロジー構築プログラム。PDB から parm7/rst7 を生成。 |
| **antechamber** | — | AmberTools のプログラム。小分子に GAFF2 原子型と AM1-BCC 部分電荷を割り当てる。 |
| **parmchk2** | — | AmberTools のプログラム。GAFF2 型付けで不足する力場パラメータをチェック・補完する。 |
| **GAFF2** | General Amber Force Field 2 | 有機小分子向けの汎用 Amber 力場。リガンド（基質、補因子）のパラメータ化に使用。 |
| **ff19SB** | — | タンパク質向け Amber 力場（2019 年版）。mlmm のデフォルト。 |
| **ff14SB** | — | タンパク質向け Amber 力場（2014 年版）。`--ff-set ff14SB` で選択可能。 |
| **AM1-BCC** | AM1 Bond Charge Corrections | 半経験的 AM1 法に基づく部分電荷割り当て方法。antechamber で使用。HF/6-31G* RESP 電荷を近似。 |

## 反応経路・最適化

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **MEP** | Minimum Energy Path | 反応物から生成物へ至る最小エネルギー経路（ポテンシャルエネルギー面上の最も低い経路）。 |
| **TS** | Transition State | ポテンシャルエネルギー面上の一次鞍点（first-order saddle point）。反応座標方向にのみ負の曲率（虚振動数）を 1 つ持つ停留点。 |
| **n_imag** | Number of imaginary modes（虚振動の数） | 虚振動の分類基準（既定では ν < −5.00 cm⁻¹）より下の振動モードの数。TS では n_imag = 1 で、`result.json` には `tsopt` では `n_imaginary_modes`、`freq` では `n_imaginary` として記録されます。 |
| **IRC** | Intrinsic Reaction Coordinate | TS から反応物側・生成物側へ向かう、質量重み付き最急降下経路。TS の接続検証によく使われます。 |
| **GSM** | Growing String Method | 端点からストリング（画像列）を伸長・最適化して MEP を近似する手法。 |
| **DMF** | Direct Max Flux | 反応座標方向のフラックスを最大化することで MEP を最適化する chain-of-states 手法。`--mep-mode dmf` で選択。 |
| **HEI** | Highest-Energy Image | MEP 上でエネルギーが最大の画像。TS の初期推定としてよく使われます。 |
| **画像（Image）** | — | 経路上の 1 つの構造（1 ノード）。 |
| **セグメント** | — | 2 つの隣接する端点を結ぶ MEP（例: R → I1, I1 → I2, …）。 |
| **ねじれ** | Kink | 配座だけが変わる経路の区間（セグメント）。`path-search` が HEI の両側で最適化した 2 つの構造（End1 と End2。[path-search の処理の仕組み](path-search.md#処理の仕組みと計算仕様) の 2）の間で、共有結合が変わらない区間を指します。`path-search` は新しい GSM・DMF の経路の代わりに、線形補間のノードを数個（`search.kink_max_nodes`、既定 3）入れて 1 つずつ最適化します。 |

## 最適化アルゴリズム

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **L-BFGS** | Limited-memory BFGS | 勾配履歴から Hessian を近似する準ニュートン法。`--opt-mode grad` で使用。 |
| **RFO** | Rational Function Optimization | 明示的な Hessian 情報を使用する信頼領域最適化法。`--opt-mode hess` で使用。 |
| **RS-I-RFO** | Restricted-Step Image-RFO | 1 つの負固有値方向に沿う、鞍点（TS）最適化用の RFO 変種。 |
| **Dimer** | Dimer Method | 低曲率方向を追跡する TS 最適化法。mlmm-toolkit の Hessian-guided Dimer は初期および定期的な活性部分空間 Hessian を使うため、活性自由度が多い系ではランダムな初期方向より頑健です。`--opt-mode grad` の TSOPT で使用。 |
| **PHVA** | Partial Hessian Vibrational Analysis | アクティブ（非凍結）原子の Hessian ブロックのみを使用した振動解析。`freq` のデフォルト。 |

## 機械学習・計算機

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **MLIP** | Machine Learning Interatomic Potential | 量子化学データから学習し、構造からエネルギー・力を予測する原子間ポテンシャル。 |
| **UMA** | Universal Models for Atoms | Meta が公開している事前学習 MLIP 群。mlmm のデフォルト MLIP バックエンド。`--backend uma` で選択（デフォルト）。 |
| **ORB** | ORB Models | Orbital Materials が提供する MLIP バックエンド。`--backend orb` で選択。`pip install "mlmm-toolkit[orb]"` で追加インストール。 |
| **MACE** | MACE (Message-passing Atomic Cluster Expansion) | 等変メッセージパッシングに基づく MLIP バックエンド。`--backend mace` で選択。専用の conda 環境で `pip uninstall fairchem-core`（UMA の pin が `e3nn` で衝突するため）を実行してから `pip install mace-torch` でインストールします。 |
| **AIMNet2** | Atoms In Molecules Network 2 | ニューラルネットワークベースの MLIP バックエンド。`--backend aimnet2` で選択。`pip install "mlmm-toolkit[aimnet]"` で追加インストール。 |
| **解析 Hessian** | Analytical Hessian | バックエンドの微分可能な、またはネイティブの Hessian 計算で二階微分を求めます。計算時間とメモリはバックエンドと系によって変わります。UMA、ORB、MACE、AIMNet2 で使えます。 |
| **有限差分** | Finite Difference | 変位させた構造の力から二階微分を近似します。計算時間とメモリはバックエンドと系によって変わり、すべての MLIP バックエンドで使えます。 |

## 量子化学

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **QM** | Quantum Mechanics | DFT、HF、post-HF などの第一原理電子状態計算。 |
| **DFT** | Density Functional Theory | 電子密度汎関数に基づく電子状態計算法。 |
| **Hessian** | — | エネルギーの二階微分行列。振動解析や TS 最適化に使用します。 |
| **SP** | Single Point | 固定構造での計算（最適化なし）。高精度エネルギー補正によく使用。 |
| **スピン多重度** | Spin Multiplicity | 2S+1（S は全スピン）。一重項 = 1、二重項 = 2、三重項 = 3 など。 |

## 構造生物学・ポケット抽出

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **PDB** | Protein Data Bank | タンパク質などの三次元構造を表す標準フォーマット（およびデータベース）。 |
| **XYZ** | — | 元素記号と直交座標を並べたシンプルなテキスト形式。計算のコマンドは、XYZ を `--ref-pdb` と一緒に受け付けます。 |
| **GJF** | Gaussian Job File | Gaussian の入力形式。`oniom-export` が `g16` モードで書き出し、`oniom-import` と `bond-summary` が読み込みます。 |
| **ポケット** | Active-site Pocket | `extract` サブコマンドで基質周辺から切り出した部分構造。ML/MM ワークフローでは ML 領域と周辺 MM 環境を定義する。 |
| **抽出用リンク水素** | Extractor-only Link Hydrogen | `extract --add-linkh` がポケット確認用に付けるキャップ水素です。ML/MM 計算では使わず、リンク水素を置く結合は ML/MM 境界にある parm7 の結合から決まります。 |
| **主鎖** | Backbone | タンパク質の主骨格（N–Cα–C–O 原子）。`--exclude-backbone` で除外可能。 |
| **B-factor** | Temperature Factor | PDB の温度因子列。mlmm では 3 層への割り当てをエンコードするために使用（0.0, 10.0, 20.0）。 |

## 熱化学

| 用語 | 正式名称 | 説明 |
|------|----------|------|
| **ZPE** | Zero-Point Energy（零点エネルギー） | 0 K での振動エネルギー。電子エネルギーへの量子補正。 |
| **ギブズ自由エネルギー** | Gibbs Free Energy (G) | G = H - TS。熱・エントロピー寄与を含む。 |
| **エンタルピー** | (H) | H = E + PV。定圧での全熱含量。 |
| **エントロピー** | (S) | 無秩序さの尺度。ギブズ自由エネルギーに −TS として寄与。 |
| **QRRHO** | Quasi-Rigid-Rotor Harmonic Oscillator | Grimme の低振動数補正を含む熱化学近似。`freq` で自動適用。 |

## 単位・定数

| 用語 | 説明 |
|------|------|
| **Hartree** | 原子単位系のエネルギー。1 Hartree ≈ 627.5 kcal/mol ≈ 27.21 eV。 |
| **kcal/mol** | 反応エネルギー表現でよく使われる単位。 |
| **kJ/mol** | キロジュール/モル。1 kcal/mol ≈ 4.184 kJ/mol。 |
| **eV** | 電子ボルト。1 eV ≈ 23.06 kcal/mol。 |
| **Bohr** | 原子単位系の長さ。1 Bohr ≈ 0.529 Å。 |
| **Å（オングストローム）** | 10⁻¹⁰ m。原子間距離の標準単位。 |
| **cm⁻¹** | 波数（逆センチメートル）。振動数の標準単位。虚振動数は負の値で表されます。 |
| **虚振動数** | Hessian 行列の負の固有値に対応する振動数。TS では 1 本のみ存在（一次鞍点）。負の cm⁻¹ 値で報告されます。 |

(ja-frequency-thresholds)=
### 虚振動の分類基準と QRRHO のローター閾値

虚振動の分類基準と QRRHO のローター閾値は用途が異なります。

| 閾値 | 役割 | 定義場所 |
|------|------|----------|
| **ν < −5.00 cm⁻¹** | 既定の虚振動分類基準。 | 設定で変えられます（`freq.zero_cutoff_cm`） |
| **100 cm⁻¹** | *QRRHO のローター閾値*（Grimme）。`freq` の熱化学計算では、これ未満の **正の** 低振動モードのエントロピーを、調和振動子の値から自由回転子の値へ滑らかに切り替えます。変わるのはエントロピーとギブズ自由エネルギーだけです | 固定（mlmm-toolkit の設定では変えられません） |

## CLI 規則

ブール値オプション、残基セレクタ、原子セレクタの書き方は [共通オプションと残基・原子の指定](cli-conventions.md) にまとめています。

## 使用上の注意点

* **虚振動と負の符号**: n_imag に数えるのは、虚振動の分類基準より下のモードだけです。振動数はすべて符号付きで出力され、負の値の総数（`result.json` の `n_negative_modes`）は収束の判定を変えない別の診断値です。

## 関連ドキュメント

- [はじめに](getting-started.md) — 最短の実行と次に読むページ
- [インストール](installation.md) — セットアップと依存関係
- [all](all.md) — ポケット抽出、MEP 探索、後処理の全体像
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法
- [YAML 設定リファレンス](yaml-reference.md) — 設定ファイルの仕様
- [MLIP バックエンド](backends.md) — MLIP バックエンドの詳細
- [ML/MM 計算機](mlmm-calc.md) — ONIOM の結合、リンク原子、MM Hessian
