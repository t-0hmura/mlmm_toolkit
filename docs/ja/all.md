# `all`

`mlmm all` は全系レイヤード PDB 上で ML/MM 機構パイプライン全体を 1 コマンドで実行します。手作業で `extract` → `mm-parm` → `define-layer` → `scan` / `path-search` → `tsopt` → `irc` / `freq` / `dft` を連結する代わりに、`all` 1 つで完結します。活性部位抽出、MM トポロジー準備、ML/MM レイヤー割り当て、任意の段階的スキャン、MEP 探索（デフォルトは再帰的 `path-search`）、任意の後処理（TS 最適化、EulerPC IRC、熱化学、DFT 一点計算、DFT//MLIP ダイアグラム）を順に実行します。ML 領域のデフォルト MLIP バックエンドは UMA で、`-b/--backend` で他のバックエンドを選択できます。

`all` は渡す入力に応じて 3 つのモードのいずれかで動作します:

- **マルチ構造アンサンブル** -- 反応順に 2 つ以上の完全 PDB を提供し、複数構造にまたがる GSM MEP 探索を駆動する。
- **単一構造 + 段階的スキャン** -- 1 つの PDB と `--scan-lists` を提供する。各リテラルがスキャンステージとなり、緩和済みの端点が MEP の端点となる。
- **TSOPT のみ** -- 1 つの PDB を提供し `--tsopt` を設定（`--scan-lists` なし）して、MEP 探索なしで TS 最適化を直接実行する。

```{important}
`--tsopt` は **TS 候補**を生成します。`all` は検証のために IRC と freq を自動実行しますが、機構解釈の前に必ず結果（虚振動数モード + 端点の結合性）を確認してください。
```

## 実行例

コマンド形式:

```bash
mlmm all -i INPUT1 [INPUT2...] -c SUBSTRATE [options]
```

コアオプションは `mlmm all --help`、全オプション一覧は `mlmm all --help-advanced` で確認できます。

完全な後処理付きのマルチ構造 MEP:

```bash
mlmm all -i R.pdb P.pdb -c "SAM,GPP" -l "SAM:1,GPP:-3" \
 --tsopt --thermo --dft --out-dir ./result_all
```

単一構造 + 段階的スキャン（2 ステージ）:

```bash
mlmm all -i A.pdb -c "308,309" --scan-lists "[(12,45,1.35)]" "[(10,55,2.20)]" \
 --multiplicity 1 --out-dir ./result_scan_all
# 1 つのリテラルで複数の結合を同時に駆動可能: '[(10,55,2.20),(23,34,1.80)]'
```

TSOPT のみの検証（単一入力、MEP 探索なし）:

```bash
mlmm all -i A.pdb -c "GPP,MMT" -l "GPP:-3,MMT:-1" \
 --tsopt --thermo --dft --out-dir result_tsopt_only
```

xTB 点電荷埋め込み付きの ORB バックエンド:

```bash
mlmm all -i R.pdb P.pdb -c "SAM,GPP" -l "SAM:1,GPP:-3" \
 --backend orb --embedcharge --out-dir ./result_all_orb
```

PDB コンパニオンはテンプレートが利用可能な場合に生成され、`--convert-files/--no-convert-files`（デフォルト有効）で制御されます。

## 処理の流れ

1. **活性部位抽出と ML 領域定義**（複数入力時はマルチ構造の和集合）
   - 基質を定義します（`-c/--center`、PDB、残基 ID、または残基名で指定）。
   - 任意で `--ligand-charge` を総数値（分配）またはマッピング（例: `GPP:-3,MMT:-1`）として提供します。
   - 抽出器は入力ごとのポケット PDB を `<out-dir>/_work/pockets/` に書き出します。最初のポケットが `<out-dir>/ml_region.pdb`（`--model-pdb` として再利用可能な成果物）としてコピーされ、後続の全 ML/MM 計算の ML 領域を定義します。
   - 抽出器の**最初のモデルの ML 領域の総電荷**が後続ステップの総電荷として使用され、丸め処理が発生した場合はコンソールに通知されます。
   - `-c/--center` を省略した場合は抽出をスキップし、完全入力構造をそのまま使用します。

2. **ML/MM 準備（parm7 + レイヤー割り当て）**
   - 最初の完全入力 PDB に対して `mm_parm` を一度実行し、`<out-dir>/mm_parm/<input_basename>.parm7` / `.rst7`（`--parm` として再利用可能な成果物）を構築します。これは自動的に `--parm` として渡されます。
   - 各完全系 PDB に対して `define-layer` を実行し、ML 領域定義に基づく 3 層 B 因子（ML=0.0、MovableMM=10.0、FrozenMM=20.0）を付与します。レイヤード全系 PDB は `<out-dir>/layered/` に書き出されます。

3. **任意の段階的スキャン（単一入力のみ）**
   - 完全入力 PDB が 1 つのみで `--scan-lists` が指定された場合、レイヤード全系 PDB に対して ML/MM 計算機を使用した段階的な結合距離駆動スキャンを実行します。
   - 各ステージの最終緩和構造（`stage_XX/result.pdb`）が中間体/生成物候補として収集されます。
   - 経路探索の入力系列は `[初期レイヤード PDB, stage_01/result.pdb, stage_02/result.pdb,...]` となります。

4. **全系レイヤード PDB での MEP 探索**
   - すべての MEP 計算は全系レイヤード PDB（`--parm` + `--detect-layer`）上で実行されます（ポケット上ではありません）。
   - **`path-search`（デフォルト）:** 自動精密化を含む再帰的 `path_search` を実行し、多段階反応を自動検出して各素反応の詳細な MEP を構築します。複雑な多段階反応では手動での試行錯誤が必要な場合があります。両モードとも Stage 5 後処理に対応。
   - **`--no-refine-path`:** 隣接ペアごとに単一パス `path-opt` GSM を実行後、軌跡を結合、セグメントごとの HEI 抽出、結合変化検出、`summary.json` 書き出しまで行い、Stage 5 後処理（TSOPT、thermo、DFT）が利用可能。
   - マルチ入力実行では、元の完全 PDB がマージ参照として自動的に供給されます。スキャン由来の系列（単一構造の場合）では、元の完全 PDB 1 つがすべての入力の参照テンプレートとして再利用されます。

5. **サマリーと任意の後処理**
   - MEP エンジン生出力（セグメントごとの軌跡、全 MEP 軌跡、エンジンの `summary.json`）は `<out-dir>/_work/path_search/`（`--no-refine-path` 使用時は `<out-dir>/_work/path_opt/`）に書き出され、マージ済み成果物（`mep.pdb`・`mep_trj.xyz`・`mep_plot.png`・`energy_diagram_MEP.png`）は `<out-dir>/` へ移動され、`summary.{json,log}` はコピーされます。
   - `--tsopt`: 各 HEI で TS を最適化し、EulerPC IRC を実行し、セグメントエネルギーダイアグラムを描画します。
   - `--thermo`: (R, TS, P) で ML/MM 熱化学を計算し、Gibbs ダイアグラムを追加します。
   - `--dft`: (R, TS, P) で DFT 一点計算を実行し、DFT ダイアグラムを追加します。`--thermo` と組み合わせると、DFT//MLIP Gibbs ダイアグラムも生成されます。
   - VRAM に余裕がある場合は `--hessian-calc-mode` を `Analytical` に設定することを強く推奨します（デフォルトの FiniteDifference より優先）。

6. **TSOPT のみモード**（単一入力、`--tsopt`、`--scan-lists` なし）
   - ステップ (4)-(5) をスキップし、レイヤード全系 PDB で `tsopt` を実行し、EulerPC IRC と両端の極小化を行い、R-TS-P の ML/MM エネルギーダイアグラムを構築し、任意で Gibbs、DFT、DFT//MLIP ダイアグラムを追加します。
   - このモードでのみ、**より高いエネルギー**の IRC 端点が反応物 (R) として採用されます。

## 出力

ツリーは 3 つのゾーンで構成されます: **ルート直下の成果物**、**`segments/seg_NN/` 配下のセグメント別成果物**、**`_work/` 配下のパイプライン作業領域**（結果を取り出したあとは削除して構いません）。最初に確認する 3 つは `summary.log`、`summary.json`、`mep.pdb`（連結した反応経路。ルートへ移動。MEP エンジン生出力は `_work/path_search/` に残る）です。

```text
<out-dir>/
 summary.json                          # トップレベルサマリーのミラー（path_search 実行時）
 summary.log
 mep.pdb                               # 連結 MEP 経路（ルートにコピー）
 mep_trj.xyz
 mep_plot.png                          # MEP 生エネルギープロファイル
 energy_diagram_MEP.png                # 全セグメント MEP 障壁
 energy_diagram_UMA_all.png            # 集約後処理ダイアグラム（有効時）
 energy_diagram_G_UMA_all.png
 energy_diagram_DFT_all.png
 energy_diagram_G_DFT_plus_UMA_all.png
 irc_plot_all.png
 ml_region.pdb                         # ML 領域定義（--model-pdb として再利用可能）
 mm_parm/<input1>.parm7,.rst7          # 最初の完全酵素入力 PDB から生成した MM トポロジー（--parm として再利用可能）
 layered/                              # レイヤード全系 PDB（B 因子アノテーション付き、再利用可能な入力）
 segments/                             # 反応セグメント別の成果物
  seg_NN/                              # 2 桁インデックス (1 始まり)、例: seg_01, seg_02
   reactant.pdb · ts.pdb · product.pdb  # 正規 R/TS/P
   ts/...                              # TS 最適化 + EulerPC IRC（--tsopt）
   irc/...
   freq/...                            # --thermo の場合
   dft/...                             # --dft の場合
   structures/{reactant,ts,product}.pdb  # 入れ子コピー + 生 IRC 端点
   energy_diagram_{UMA,G_UMA,DFT,G_DFT_plus_UMA}.png
 _work/                               # パイプライン作業領域（削除可）
  pockets/                             # 入力ごとのポケット PDB（複数構造は統合）
  scan/                                # 単一構造+スキャンモードの場合のみ（stage_01/result.pdb …）
  path_search/                         # MEP エンジン生出力（--no-refine-path 時は path_opt/）
   summary.{json,log} · seg_NN_mep/    # 生の per-segment GSM 軌跡（マージ済み MEP 成果物は root へ移動）
```

**TSOPT のみモード**（単一入力 + `--tsopt`、`--scan-lists` なし）では MEP ステージが無く、最適化済み R/TS/P と `ts/`・`irc/`・`freq/`・`dft/` は `segments/seg_01/` 配下に生成され、`_work/path_search/` は存在しません。

`-v 2` ではコンソールに抽出、MM 準備、スキャンステージ、MEP の進捗（GSM）、ステージごとのタイミングが要約されます。{ref}`ja-verbosity-levels` を参照してください。

### `summary.log` の読み方
ログは番号付きセクションで構成されています:
- **[1] グローバル MEP 概要** -- イメージ/セグメント数、MEP 軌跡プロットパス、集約 MEP エネルギーダイアグラム。
- **[2] セグメントレベル MEP サマリー（MLIP 経路）** -- セグメントごとの障壁、反応エネルギー、結合変化サマリー。
- **[3] セグメントごとの後処理（TSOPT / Thermo / DFT）** -- セグメントごとの TS 虚振動数チェック、IRC 出力、エネルギーテーブル。
- **[4] エネルギーダイアグラム（概要）** -- MEP/UMA/Gibbs/DFT シリーズのダイアグラムテーブルと任意のクロスメソッドサマリーテーブル。
- **[5] 出力ディレクトリ構造** -- インラインアノテーション付きの生成ファイルのコンパクトツリー。

### `summary.json` の読み方
summary.json はコンパクトな機械可読サマリーです。主なトップレベルキー:
- `out_dir`、`n_images`、`n_segments` -- 実行メタデータと総数。
- `segments` -- `index`、`tag`、`kind`、`barrier_kcal`、`delta_kcal`、`bond_changes` を持つセグメントごとのエントリリスト。
- `energy_diagrams`（任意）-- `labels`、`energies_kcal`、`energies_au`、`ylabel`、`image` パスを持つダイアグラムペイロード。

## CLI オプション

> **注意:** 表示されるデフォルト値はオプション未指定時に使用されます。完全なフラグ一覧は生成された command reference（`reference/commands/`）にあり、以下の表は説明が必要なオプションを扱います。

### 入出力オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順の 2 つ以上の完全 PDB（`--scan-lists`（段階的スキャン）または `--tsopt`（TSOPT のみ）の場合のみ単一入力可）。 | 必須 |
| `-c, --center TEXT` | 基質指定（PDB パス、残基 ID（`308,309`）、または残基名（`SAM,GPP`））。省略時は抽出をスキップし完全構造をそのまま使用。 | _None_ |
| `-l, --ligand-charge TEXT` | 非標準残基の総電荷または残基別マッピング（例: `GPP:-3,MMT:-1`）。 | _None_ |
| `-q, --charge INT` | 総電荷を強制指定（最優先の上書き）。 | _None_ |
| `-o, --out-dir PATH` | トップレベル出力ディレクトリ。 | `./result_all/` |
| `--parm FILE` | 全系の AMBER parm7 トポロジーファイル。省略時は `mm_parm` で自動生成。 | _None_ |
| `--model-pdb FILE` | 構築済み ML 領域 PDB。指定時は ML 領域決定をスキップし、このファイルで ML 領域を直接定義。 | _None_ |
| `--ref-pdb FILE` | XYZ 入力用の参照 PDB。入力が XYZ の場合に PDB メタデータ（残基、鎖、B 因子）を復元するために必要。 | _None_ |
| `--convert-files/--no-convert-files` | テンプレート利用可能時の XYZ/TRJ から PDB コンパニオンのグローバルトグル。 | `True` |
| `--dump/--no-dump` | オプティマイザーダンプを保存。常に `path-search`/`path-opt` に転送。`scan`/`tsopt` にはここで明示設定時のみ転送。`freq` は明示的に `--no-dump` を指定しない限りデフォルトで dump=True。 | `False` |
| `--config FILE` | 先に適用するベース YAML。 | _None_ |
| `--show-config/--no-show-config` | 実行前に解決済み設定を表示。 | `False` |
| `--dry-run/--no-dry-run` | 実行せず検証と計画表示のみ行う。`--help-advanced` に表示。 | `False` |

### 抽出オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-r, --radius FLOAT` | ポケット包含カットオフ (Å)。 | `2.6` |
| `--radius-het2het FLOAT` | 独立したヘテロ-ヘテロカットオフ (Å)。 | `0.0` |
| `--include-H2O/--no-include-H2O` | 水分子（HOH/WAT/TIP3/SOL）を含める。 | `True` |
| `--exclude-backbone/--no-exclude-backbone` | 非基質アミノ酸の主鎖原子を除去。 | `False` |
| `--add-linkH/--no-add-linkH` | 切断結合にリンク水素を付加。 | `False` |
| `--selected-resn TEXT` | 強制包含する残基。 | `""` |
| `--modified-residue TEXT` | 修飾アミノ酸残基名をカンマ区切りで指定（任意で電荷付き）。主鎖切断と電荷計算にアミノ酸として扱う。例: `HD1,HD2,HD3` または `HD1:0,SEP:-2`。 | `""` |

### MM 準備オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--auto-mm-ff-set {ff19SB\|ff14SB}` | `mm_parm` 用の力場セット（ff19SB は OPC3、ff14SB は TIP3P を使用）。 | `ff19SB` |
| `--auto-mm-add-ter/--auto-mm-no-add-ter` | リガンド/水/イオンブロック周囲の TER 挿入を制御。 | `True` |
| `--auto-mm-keep-temp` | `mm_parm` の一時作業ディレクトリを保持（デバッグ用）。 | `False` |
| `--auto-mm-ligand-mult TEXT` | `mm_parm` に転送するスピン多重度マッピング（例: `GPP:2,SAM:1`）。省略時は全リガンドに 1 を使用。 | _None_ |

### MEP 探索オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `-b, --backend CHOICE` | ML バックエンド: `uma`（デフォルト）、`orb`、`mace`、`aimnet2`。全計算サブコマンドに転送。 | `uma` |
| `--embedcharge/--no-embedcharge` | xTB 点電荷埋め込み補正の有効化。MM 環境から ML 領域への静電的影響を考慮。 | `False` |
| `--embedcharge-cutoff FLOAT` | xTB 埋め込み用 MM 原子のカットオフ半径（Å）。 | `12.0` |
| `--cmap/--no-cmap` | model parm7 に CMAP（骨格クロスマップ二面角補正）を含めるかどうか。デフォルト: 無効（Gaussian ONIOM と同一）。 | `--no-cmap` |
| `--max-nodes INT` | セグメント GSM の内部ノード数。 | `20` |
| `--max-cycles INT` | GSM マクロサイクルの最大数。 | `300` |
| `--climb/--no-climb` | セグメント GSM の TS 精密化を有効化。 | `True` |
| `--opt-mode [grad\|hess]` | スキャン/path-opt/path-search と単一構造最適化のプリセット（`grad` → LBFGS/Dimer、`hess` → RFO/RSIRFO）。 | `grad` |
| `--opt-mode-post [grad\|hess]` | TSOPT/IRC 後端点最適化向けのプリセット上書き（`grad` → Dimer/LBFGS、`hess` → RS-I-RFO/RFO）。 | `hess` |
| `--thresh TEXT` | 収束プリセット（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`、`never`）。実効デフォルト: path-opt は `gau_loose`、scan は `gau`。 | _None_ |
| `--thresh-post TEXT` | IRC 後端点最適化の収束プリセット。 | `baker` |
| `--preopt/--no-preopt` | セグメント化前に端点を事前最適化。 | `True` |
| `--refine-path/--no-refine-path` | True（デフォルト）の場合は再帰的 `path-search`、False の場合は `path-opt` セグメントチェーン（単一パス GSM + 軌跡結合 + HEI 抽出 + 結合変化検出 + summary.json）。両モードとも Stage 5（TSOPT/thermo/DFT）対応。 | `True` |
| `--hessian-calc-mode CHOICE` | ML/MM ヘシアンモード（`Analytical` または `FiniteDifference`）。 | `FiniteDifference` |
| `--detect-layer/--no-detect-layer` | 入力 PDB の B 因子（B=0/10/20）から ML/MM レイヤーを検出。無効時は下流ツールで `--model-pdb` または `--model-indices` が必要。 | `True` |

TSOPT の最適化モード選択順: `--opt-mode-post`（設定時）-> `--opt-mode`（明示指定時のみ）-> TSOPT デフォルト（`hess`）。

### スキャンオプション（単一入力実行）

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-s, --scan-lists TEXT...` | 段階的スキャン: `(i,j,target_A)` タプル。 | _None_ |
| `--scan-out-dir PATH` | スキャン出力ディレクトリの上書き。 | _None_ |
| `--scan-one-based/--no-scan-one-based` | スキャンインデックスの解釈を上書き（True = 1 始まり、False = 0 始まり）。 | _None_ |
| `--scan-max-step-size FLOAT` | 最大ステップサイズ (Å)。 | _デフォルト_ |
| `--scan-bias-k FLOAT` | 調和バイアス強度 (eV/Å^2)。 | _デフォルト_ |
| `--scan-relax-max-cycles INT` | ステップごとの緩和最大サイクル。 | _デフォルト_ |
| `--scan-preopt/--no-scan-preopt` | スキャン事前最適化トグルの上書き。 | _None_ |
| `--scan-endopt/--no-scan-endopt` | スキャンステージ終端最適化の上書き。 | _None_ |

### 後処理 + freq / DFT 上書きオプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--tsopt/--no-tsopt` | 反応セグメントごとに TS 最適化 + EulerPC IRC を実行。 | `False` |
| `--thermo/--no-thermo` | R/TS/P で振動解析 (`freq`) を実行。 | `False` |
| `--dft/--no-dft` | R/TS/P で DFT 一点計算を実行。 | `False` |
| `--flatten/--no-flatten` | `tsopt` での余分な虚振動数モードフラットニングを有効化。 | `False` |
| `--tsopt-max-cycles INT` | `tsopt --max-cycles` の上書き。 | _デフォルト_ |
| `--tsopt-out-dir PATH` | tsopt サブディレクトリのカスタマイズ。 | _None_ |
| `--freq-out-dir PATH` | freq 出力ディレクトリの上書き。 | _None_ |
| `--freq-max-write INT` | 出力する最大モード数。 | _デフォルト_ |
| `--freq-amplitude-ang FLOAT` | モードアニメーション振幅 (Å)。 | _デフォルト_ |
| `--freq-n-frames INT` | モードアニメーションのフレーム数。 | _デフォルト_ |
| `--freq-sort TEXT` | モードソート方法。 | _デフォルト_ |
| `--freq-temperature FLOAT` | 熱化学温度 (K)。 | _デフォルト_ |
| `--freq-pressure FLOAT` | 熱化学圧力 (atm)。 | _デフォルト_ |
| `--dft-out-dir PATH` | DFT 出力ディレクトリの上書き。 | _None_ |
| `--dft-func-basis TEXT` | 汎関数/基底関数ペア。 | _デフォルト_ |
| `--dft-max-cycle INT` | 最大 SCF 反復数。 | _デフォルト_ |
| `--dft-conv-tol FLOAT` | SCF 収束閾値。 | _デフォルト_ |
| `--dft-grid-level INT` | PySCF グリッドレベル。 | _デフォルト_ |
| `--dft-engine [gpu\|cpu]` | DFTエンジン（GPU or CPU PySCF）。 | _None_ |

## YAML 設定

`all` は階層化された YAML をサポートします:

- `--config FILE`: ベース設定。

`defaults < config < CLI < override-yaml`

解決後の YAML が下流サブコマンドへ転送されます。各ツールは独自ドキュメントに記載されたセクションを読み取ります:

| サブコマンド | YAML セクション |
|------------|---------------|
| [`path-search`](path-search.md) | `geom`, `calc`/`mlmm`, `gs`, `opt`, `lbfgs`, `bond`, `search` |
| [`scan`](scan.md) | `geom`, `calc`/`mlmm`, `opt`, `lbfgs` |
| [`tsopt`](tsopt.md) | `geom`, `calc`/`mlmm`, `opt`, `hessian_dimer`, `rsirfo` |
| [`freq`](freq.md) | `geom`, `calc`/`mlmm`, `freq`, `thermo` |
| [`dft`](dft.md) | `dft` |

> **注意:** CLI 値の後に適用されます。

**最小の YAML 例:**
```yaml
calc:
 charge: 0
 spin: 1
mlmm:
 real_parm7: real.parm7
 model_pdb: ml_region.pdb
 backend: uma                    # ML バックエンド (uma/orb/mace/aimnet2)
 embedcharge: false              # xTB 点電荷埋め込み補正
 uma_model: uma-s-1p1            # uma-s-1p1 | uma-m-1p1
 hessian_calc_mode: Analytical     # VRAM に余裕がある場合に推奨
gs:
 max_nodes: 12
 climb: true
dft:
 grid_level: 6
```

すべての YAML オプションの完全なリファレンスは **[YAML 設定リファレンス](yaml-reference.md)** を参照してください。

## 注記

入力形式は抽出の有無に依存します:

- 抽出有効時（`-c/--center`）: 入力は残基を特定するため **PDB** ファイルが必要。
- 抽出スキップ時: 入力は **PDB/XYZ** が使用可能。
- マルチ構造実行には 2 つ以上の構造が必要。

電荷は優先度の高い順に解決されます -- `-q/--charge`（明示的な CLI 上書き）-> ポケット抽出（`-c` 指定時、アミノ酸 + イオン + `--ligand-charge` の合計）-> `-l, --ligand-charge` フォールバック（抽出スキップ時）-> デフォルト（未解決の電荷はエラー）。スピンの解決: `--multiplicity`（CLI）-> デフォルト（1）。正しい電荷伝播のため、非標準基質には常に `--ligand-charge` を指定してください。最初のモデルの ML 領域の総電荷は最も近い整数に丸められ、丸め処理が発生した場合はコンソールに通知されます。

## 関連項目

- [extract](extract.md) -- 単独の ML 領域決定（`all` が内部で呼び出し）
- [mm_parm](mm-parm.md) -- AMBER トポロジー構築（`all` が内部で呼び出し）
- [path-search](path-search.md) -- 単独の再帰的 MEP 探索
- [tsopt](tsopt.md) -- 単独の TS 最適化
- [freq](freq.md) -- 振動解析と熱化学
- [dft](dft.md) -- DFT 一点計算
- [trj2fig](trj2fig.md) -- 軌跡からエネルギープロファイルをプロット
- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- よくあるエラーと対処法
- [YAML リファレンス](yaml-reference.md) -- 完全な YAML 設定オプション
- [用語集](glossary.md) -- MEP、TS、IRC、GSM の定義
