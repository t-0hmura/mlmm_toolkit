# `all`

## 概要

> **要約:** エンドツーエンドの酵素反応ワークフロー -- 活性部位抽出、ML/MM レイヤー割り当て、MM トポロジー準備、任意の段階的スキャン、全系レイヤード PDB での MEP 探索（GSM）、任意の TS 最適化、擬似 IRC、熱化学、DFT、DFT//UMA ダイアグラム。

`mlmm all` は全系レイヤード PDB 上で ML/MM を用いるワンショットパイプラインを実行します。以下の 3 つのモードをサポートしています:

- **マルチ構造アンサンブル** -- 反応順に 2 つ以上の完全 PDB を提供します。活性部位領域を抽出（ML 領域定義用）、MM トポロジー構築、ML/MM レイヤー割り当て後、レイヤード全系 PDB で GSM MEP 探索を実行し、任意でセグメントごとの後処理（TSOPT/freq/DFT）を実行します。
- **単一構造 + 段階的スキャン** -- 1 つの PDB と `--scan-lists` を提供します。スキャンで中間体/生成物候補を生成し、MEP の端点として使用します。
  - 1 つの `--scan-lists` リテラルで単一のスキャンステージを実行します。
  - 複数ステージは 1 つの `--scan-lists` フラグの後に複数の値として渡します（フラグ自体は繰り返せません）。
- **TSOPT のみ** -- 1 つの PDB を提供し、`--scan-lists` を省略して `--tsopt` を設定します。レイヤード全系 PDB で TS 最適化を実行し、擬似 IRC、両端の極小化、エネルギーダイアグラムの構築を行います。

```{important}
`--tsopt` は **TS 候補**を生成します。`all` は検証のために IRC と freq を自動実行しますが、機構解釈の前に必ず結果（虚振動数モード + 端点の結合性）を確認してください。
```

### 概要
- **用途:** 活性部位抽出から MEP 探索、任意の TS/freq/DFT 後処理まで、ML/MM による完全なエンドツーエンドの酵素反応研究を行いたい場合。
- **手法:** 活性部位抽出 + ML/MM レイヤー割り当て + MM トポロジー（AmberTools）+ 全系レイヤード PDB 上の GSM MEP 探索（MLIP バックエンド + hessian_ff）+ 任意の TSOPT/IRC/freq/DFT。`--backend` で ML バックエンドを選択可能。
- **出力:** `summary.log`、`summary.yaml`、MEP 軌跡、エネルギーダイアグラム、セグメントごとの後処理。
- **デフォルト:** `--opt-mode grad`、`--thresh gau`、`--max-cycles 300`、`--preopt` 有効、`--climb` 有効。
- **次のステップ:** `summary.log` とエネルギーダイアグラムを確認し、精密化が必要な場合は [tsopt](tsopt.md)/[freq](freq.md)/[dft](dft.md) を単独実行。

## 最小例

```bash
mlmm all -i R.pdb P.pdb -c "SAM,GPP" --ligand-charge "SAM:1,GPP:-3" --out-dir ./result_all
```

## 出力の見方

- `result_all/summary.log`
- `result_all/summary.yaml`
- `result_all/path_search/mep.pdb`（または `result_all/path_search/seg_*/`）

## よくある例

1. TS 最適化・熱化学・DFT まで一括実行する。

```bash
mlmm all -i R.pdb P.pdb -c "SAM,GPP" --ligand-charge "SAM:1,GPP:-3" \
 --tsopt --thermo --dft --out-dir ./result_all
```

2. 単一構造 + 段階的スキャンを実行する。

```bash
mlmm all -i A.pdb -c "308,309" --scan-lists "[(12,45,1.35)]" "[(10,55,2.20)]" \
 --multiplicity 1 --out-dir ./result_scan_all
```

3. 重い処理を流す前に計画だけ確認する。

```bash
mlmm all -i R.pdb P.pdb -c "SAM,GPP" --ligand-charge "SAM:1,GPP:-3" --dry-run
```

4. ORB バックエンドと xTB 点電荷埋め込みを使用する。

```bash
mlmm all -i R.pdb P.pdb -c "SAM,GPP" --ligand-charge "SAM:1,GPP:-3" \
 --backend orb --embedcharge --out-dir ./result_all_orb
```

PDB コンパニオンはテンプレートが利用可能な場合に生成され、`--convert-files/--no-convert-files`（デフォルト有効）で制御されます。

## ワークフロー

1. **活性部位抽出と ML 領域定義**（複数入力時はマルチ構造の和集合）
   - 基質を定義します（`-c/--center`、PDB、残基 ID、または残基名で指定）。
   - 任意で `--ligand-charge` を総数値（分配）またはマッピング（例: `GPP:-3,MMT:-1`）として提供します。
   - 抽出器は入力ごとのポケット PDB を `<out-dir>/pockets/` に書き出します。最初のポケットが `<out-dir>/ml_region.pdb` としてコピーされ、後続の全 ML/MM 計算の ML 領域を定義します。
   - 抽出器の**最初のモデルのポケット総電荷**が後続ステップの総電荷として使用され、丸め処理が発生した場合はコンソールに通知されます。
   - 追加の抽出トグル: `--radius`、`--radius-het2het`、`--include-H2O/--no-include-H2O`、`--exclude-backbone/--no-exclude-backbone`、`--add-linkH/--no-add-linkH`、`--selected-resn`、`--verbose/--no-verbose`。
   - `-c/--center` を省略した場合は抽出をスキップし、完全入力構造をそのまま使用します。

2. **ML/MM 準備（parm7 + レイヤー割り当て）**
   - 最初の完全入力 PDB に対して `mm_parm` を一度実行し、`<out-dir>/mm_parm/<input_basename>.parm7` / `.rst7` を構築します。これは自動的に `--parm` として渡されます。
   - 各完全系 PDB に対して `define-layer` を実行し、ML 領域定義に基づく 3 層 B 因子（ML=0.0、MovableMM=10.0、FrozenMM=20.0）を付与します。レイヤード全系 PDB は `<out-dir>/layered/` に書き出されます。
   - このステージは `--auto-mm-ff-set`、`--auto-mm-add-ter`、`--auto-mm-keep-temp` で調整できます。

3. **任意の段階的スキャン（単一入力のみ）**
   - 完全入力 PDB が 1 つのみで `--scan-lists` が指定された場合、レイヤード全系 PDB に対して ML/MM 計算機を使用した段階的な結合距離駆動スキャンを実行します。
   - 各ステージの最終緩和構造（`stage_XX/result.pdb`）が中間体/生成物候補として収集されます。
   - 経路探索の入力系列は `[初期レイヤード PDB, stage_01/result.pdb, stage_02/result.pdb,...]` となります。

4. **全系レイヤード PDB での MEP 探索**
   - すべての MEP 計算は全系レイヤード PDB（`--parm` + `--detect-layer`）上で実行されます（ポケット上ではありません）。
   - **`--refine-path`（デフォルト）:** キンク検出と精密化を含む再帰的 `path_search` を実行。
   - **`--no-refine-path`:** 隣接ペアごとに単一パス `path-opt` GSM を実行後、軌跡を結合、セグメントごとの HEI 抽出、結合変化検出、`summary.yaml` 書き出しまで行い、Stage 4 後処理（TSOPT、thermo、DFT）も両モードで利用可能。
   - マルチ入力実行では、元の完全 PDB がマージ参照として自動的に供給されます。スキャン由来の系列（単一構造の場合）では、単一の元の完全 PDB がすべての入力の参照テンプレートとして再利用（繰り返し）されます。

5. **サマリーと任意の後処理**
   - セグメントごとの軌跡、全 MEP 軌跡、`summary.yaml` が `<out-dir>/path_search/` に書き出されます。
   - `--tsopt`: 各 HEI で TS を最適化し、EulerPC IRC を実行し、セグメントエネルギーダイアグラムを描画します。
   - `--thermo`: (R, TS, P) で ML/MM 熱化学を計算し、Gibbs ダイアグラムを追加します。
   - `--dft`: (R, TS, P) で DFT 一点計算を実行し、DFT ダイアグラムを追加します。`--thermo` と組み合わせると、DFT//UMA Gibbs ダイアグラムも生成されます。
   - 共有の上書き: `--opt-mode`、`--opt-mode-post`（TSOPT と IRC 後 endpoint-opt のモード上書き）、`--flatten/--no-flatten`、`--hessian-calc-mode`、`--tsopt-max-cycles`、`--tsopt-out-dir`、`--freq-*`、`--dft-*`。
   - VRAM に余裕がある場合は `--hessian-calc-mode` を `Analytical` に設定することを強く推奨します。

6. **TSOPT のみモード**（単一入力、`--tsopt`、`--scan-lists` なし）
   - ステップ (4)-(5) をスキップし、レイヤード全系 PDB で `tsopt` を実行し、擬似 IRC と両端の極小化を行い、R-TS-P の ML/MM エネルギーダイアグラムを構築し、任意で Gibbs、DFT、DFT//UMA ダイアグラムを追加します。
   - このモードでのみ、**より高いエネルギー**の IRC 端点が反応物 (R) として採用されます。

### 電荷とスピンの優先順位

**電荷の解決（優先度の高い順）:**

| 優先度 | ソース | 使用タイミング |
|--------|--------|----------------|
| 1 | `-q/--charge` | 明示的な CLI 上書き |
| 2 | ポケット抽出 | `-c` 指定時（アミノ酸、イオン、`--ligand-charge` の合計） |
| 3 | `--ligand-charge` | 抽出スキップ時のフォールバック |
| 4 | デフォルト | なし（未解決なら中断） |

**スピンの解決:** `--multiplicity`（CLI）-> デフォルト（1）

> **ヒント:** 正しい電荷伝播のため、非標準基質には常に `--ligand-charge` を指定してください。

### 入力の前提条件
- 抽出有効時（`-c/--center`）: 入力は残基を特定するため **PDB** ファイルが必要。
- 抽出スキップ時: 入力は **PDB/XYZ** が使用可能。
- マルチ構造実行には 2 つ以上の構造が必要。

## CLI オプション

> **注意:** 表示されるデフォルト値はオプション未指定時に使用されます。

### 入出力オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順の 2 つ以上の完全 PDB（`--scan-lists` または `--tsopt` の場合のみ単一入力可）。 | 必須 |
| `-c, --center TEXT` | 基質指定（PDB パス、残基 ID、または残基名）。省略時は抽出をスキップ。 | _None_ |
| `--ligand-charge TEXT` | 総電荷または残基別マッピング（例: `GPP:-3,MMT:-1`）。 | _None_ |
| `-q, --charge INT` | 総電荷を強制指定（最優先の上書き）。 | _None_ |
| `--out-dir PATH` | トップレベル出力ディレクトリ。 | `./result_all/` |
| `--convert-files/--no-convert-files` | テンプレート利用可能時の XYZ/TRJ から PDB コンパニオンのグローバルトグル。 | `True` |
| `--dump/--no-dump` | オプティマイザーダンプを保存。常に `path-search`/`path-opt` に転送。`scan`/`tsopt` にはここで明示設定時のみ転送。`freq` は明示的に `--no-dump` を指定しない限りデフォルトで dump=True。 | `False` |
| `--config FILE` | 先に適用するベース YAML。 | _None_ |
| `--show-config/--no-show-config` | 実行前に解決済み設定を表示。 | `False` |
| `--dry-run/--no-dry-run` | 実行せず検証と計画表示のみ行う。 | `False` |

### 抽出オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--radius FLOAT` | ポケット包含カットオフ (Å)。 | 抽出器デフォルト |
| `--radius-het2het FLOAT` | 独立したヘテロ-ヘテロカットオフ。 | 抽出器デフォルト |
| `--include-H2O/--no-include-H2O` | 水分子を含める。 | 抽出器デフォルト |
| `--exclude-backbone/--no-exclude-backbone` | 非基質アミノ酸の主鎖原子を除去。 | 抽出器デフォルト |
| `--add-linkH/--no-add-linkH` | 切断結合にリンク水素を付加。 | 抽出器デフォルト |
| `--selected-resn TEXT` | 強制包含する残基。 | `""` |
| `--verbose/--no-verbose` | INFO レベルの抽出器ログを有効化。 | 抽出器デフォルト |

### MM 準備オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--auto-mm-ff-set TEXT` | `mm_parm` 用の力場セット。 | `mm_parm` デフォルト |
| `--auto-mm-add-ter {True\|False}` | TER レコードを追加。 | `mm_parm` デフォルト |
| `--auto-mm-keep-temp {True\|False}` | 一時ファイルを保持。 | `mm_parm` デフォルト |

### MEP 探索オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `--backend CHOICE` | ML バックエンド: `uma`（デフォルト）、`orb`、`mace`、`aimnet2`。全計算サブコマンドに転送。 | `uma` |
| `--embedcharge/--no-embedcharge` | xTB 点電荷埋め込み補正の有効化。MM 環境から ML 領域への静電的影響を考慮。 | `False` |
| `--max-nodes INT` | セグメント GSM の内部ノード数。 | `10` |
| `--max-cycles INT` | GSM マクロサイクルの最大数。 | `300` |
| `--climb/--no-climb` | セグメント GSM の TS 精密化を有効化。 | `True` |
| `--opt-mode [grad\|hess]` | スキャン/path-search と単一構造最適化のプリセット（`grad` → LBFGS/Dimer、`hess` → RFO/RSIRFO）。 | `grad` |
| `--opt-mode-post [grad\|hess]` | TSOPT/IRC 後端点最適化向けのプリセット上書き（`grad` → Dimer/LBFGS、`hess` → RS-I-RFO/RFO）。 | `hess` |
| `--thresh TEXT` | 収束プリセット（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`、`never`）。 | `gau` |
| `--thresh-post TEXT` | IRC 後端点最適化の収束プリセット。 | `baker` |
| `--preopt/--no-preopt` | セグメント化前に端点を事前最適化。 | `True` |
| `--refine-path/--no-refine-path` | True の場合は再帰的 `path-search`、False の場合は `path-opt` セグメントチェーン（単一パス GSM + 軌跡結合 + HEI 抽出 + 結合変化検出 + summary.yaml）。両モードとも Stage 4（TSOPT/thermo/DFT）対応。 | `True` |
| `--hessian-calc-mode CHOICE` | ML/MM ヘシアンモード（`Analytical` または `FiniteDifference`）。 | _デフォルト_ |

TSOPT の最適化モード選択順: `--opt-mode-post`（設定時）-> `--opt-mode`（明示指定時のみ）-> TSOPT デフォルト（`hess`）。

### スキャンオプション（単一入力実行）

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--scan-lists TEXT...` | 段階的スキャン: `(i,j,target_A)` タプル。 | _None_ |
| `--scan-out-dir PATH` | スキャン出力ディレクトリの上書き。 | _None_ |
| `--scan-one-based/--no-scan-one-based` | 原子セレクターに 1 始まりインデックスを強制。 | `True` |
| `--scan-max-step-size FLOAT` | 最大ステップサイズ (Å)。 | _デフォルト_ |
| `--scan-bias-k FLOAT` | 調和バイアス強度 (eV/Å^2)。 | _デフォルト_ |
| `--scan-relax-max-cycles INT` | ステップごとの緩和最大サイクル。 | _デフォルト_ |
| `--scan-preopt/--no-scan-preopt` | スキャン事前最適化トグルの上書き。 | _デフォルト_ |
| `--scan-endopt/--no-scan-endopt` | スキャンステージ終端最適化の上書き。 | _デフォルト_ |

### 後処理オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--tsopt/--no-tsopt` | 反応セグメントごとに TS 最適化 + 擬似 IRC を実行。 | `False` |
| `--thermo/--no-thermo` | R/TS/P で振動解析 (`freq`) を実行。 | `False` |
| `--dft/--no-dft` | R/TS/P で DFT 一点計算を実行。 | `False` |
| `--flatten/--no-flatten` | `tsopt` での余分な虚振動数モードフラットニングを有効化。 | `False` |
| `--tsopt-max-cycles INT` | `tsopt --max-cycles` の上書き。 | _デフォルト_ |
| `--tsopt-out-dir PATH` | tsopt サブディレクトリのカスタマイズ。 | _None_ |

### Freq 上書き

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--freq-out-dir PATH` | freq 出力ディレクトリの上書き。 | _None_ |
| `--freq-max-write INT` | 出力する最大モード数。 | _デフォルト_ |
| `--freq-amplitude-ang FLOAT` | モードアニメーション振幅 (Å)。 | _デフォルト_ |
| `--freq-n-frames INT` | モードアニメーションのフレーム数。 | _デフォルト_ |
| `--freq-sort TEXT` | モードソート方法。 | _デフォルト_ |
| `--freq-temperature FLOAT` | 熱化学温度 (K)。 | _デフォルト_ |
| `--freq-pressure FLOAT` | 熱化学圧力 (atm)。 | _デフォルト_ |

### DFT 上書き

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--dft-out-dir PATH` | DFT 出力ディレクトリの上書き。 | _None_ |
| `--dft-func-basis TEXT` | 汎関数/基底関数ペア。 | _デフォルト_ |
| `--dft-max-cycle INT` | 最大 SCF 反復数。 | _デフォルト_ |
| `--dft-conv-tol FLOAT` | SCF 収束閾値。 | _デフォルト_ |
| `--dft-grid-level INT` | PySCF グリッドレベル。 | _デフォルト_ |

## 出力

```text
<out-dir>/
 ml_region.pdb                         # ML 領域定義（最初のポケットのコピー）
 pockets/
  pocket_<input1_basename>.pdb
  pocket_<input2_basename>.pdb
  ...
 mm_parm/
  <input1_basename>.parm7              # 最初の完全酵素入力 PDB から生成
  <input1_basename>.rst7
 scan/                                 # 単一構造+スキャンモードの場合のみ
  stage_01/result.pdb
  stage_02/result.pdb
  ...
 summary.yaml                          # トップレベルサマリーのミラー（path_search 実行時）
 summary.log
 mep_plot.png
 energy_diagram_MEP.png
 energy_diagram_UMA_all.png            # 集約後処理ダイアグラム（有効時）
 energy_diagram_G_UMA_all.png
 energy_diagram_DFT_all.png
 energy_diagram_G_DFT_plus_UMA_all.png
 irc_plot_all.png
 path_search/                          # path_search 実行時
  mep_trj.xyz
  mep.pdb
  mep_w_ref.pdb
  mep_w_ref_seg_XX.pdb
  summary.yaml
  summary.log
  mep_plot.png
  energy_diagram_MEP.png
  post_seg_XX/                         # 後処理有効時
   ts/...
   irc/...
   freq/...                            # --thermo の場合
   dft/...                             # --dft の場合
   energy_diagram_UMA.png
   energy_diagram_G_UMA.png
   energy_diagram_DFT.png
   energy_diagram_G_DFT_plus_UMA.png
 tsopt_single/                         # 単一構造 TSOPT のみモードの場合
  ts/...
  irc/...
  structures/
   reactant.pdb
   ts.pdb
   product.pdb
  freq/...                             # --thermo の場合
  dft/...                              # --dft の場合
  energy_diagram_UMA.png
  energy_diagram_G_UMA.png
  energy_diagram_DFT.png
  energy_diagram_G_DFT_plus_UMA.png
```

### `summary.log` の読み方
ログは番号付きセクションで構成されています:
- **[1] グローバル MEP 概要** -- イメージ/セグメント数、MEP 軌跡プロットパス、集約 MEP エネルギーダイアグラム。
- **[2] セグメントレベル MEP サマリー（UMA 経路）** -- セグメントごとの障壁、反応エネルギー、結合変化サマリー。
- **[3] セグメントごとの後処理（TSOPT / Thermo / DFT）** -- セグメントごとの TS 虚数振動数チェック、IRC 出力、エネルギーテーブル。
- **[4] エネルギーダイアグラム（概要）** -- MEP/UMA/Gibbs/DFT シリーズのダイアグラムテーブルと任意のクロスメソッドサマリーテーブル。
- **[5] 出力ディレクトリ構造** -- インラインアノテーション付きの生成ファイルのコンパクトツリー。

### `summary.yaml` の読み方
YAML はコンパクトな機械可読サマリーです。主なトップレベルキー:
- `out_dir`、`n_images`、`n_segments` -- 実行メタデータと総数。
- `segments` -- `index`、`tag`、`kind`、`barrier_kcal`、`delta_kcal`、`bond_changes` を持つセグメントごとのエントリリスト。
- `energy_diagrams`（任意）-- `labels`、`energies_kcal`、`energies_au`、`ylabel`、`image` パスを持つダイアグラムペイロード。

## YAML 設定

`all` は階層化された YAML をサポートします:

- `--config FILE`: ベース設定。

`defaults < config < CLI < override-yaml`

解決後の YAML が下流サブコマンドへ転送されます。各ツールは独自ドキュメントに記載されたセクションを読み取ります:

| サブコマンド | YAML セクション |
|------------|---------------|
| [`path-search`](path_search.md) | `geom`, `calc`/`mlmm`, `gs`, `opt`, `lbfgs`, `bond`, `search` |
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
 uma_model: uma-s-1p1            # UMA モデルタグ (backend=uma 時)
 ml_hessian_mode: Analytical     # VRAM に余裕がある場合に推奨
gs:
 max_nodes: 12
 climb: true
dft:
 grid_level: 6
```

すべての YAML オプションの完全なリファレンスは **[YAML 設定リファレンス](yaml_reference.md)** を参照してください。

## 注意事項

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- **Python >= 3.11** が必要です。
- **基質 (`-c/--center`) と（必要な場合は）`--ligand-charge` は実質的に必須です。**
- 形式電荷を推測できない場合は常に `--ligand-charge`（数値または残基別マッピング）を指定し、正しい総電荷がスキャン/MEP/TSOPT/DFT に伝播するようにしてください。
- 単一構造モードでは `--scan-lists` または `--tsopt` が必要です。それ以外の場合は少なくとも 2 つの構造が必要です。
- preflight チェックで、繰り返し指定した `-i/--input` が実ファイルであること、および AmberTools コマンド（`tleap`, `antechamber`, `parmchk2`）の存在を実行前に検証します。
- マージ用の参照 PDB テンプレートは元の入力から自動的に導出されます。
- 収束プリセット: `--thresh` のデフォルトは `gau`、`--thresh-post` のデフォルトは `baker`。
- 抽出半径: `--radius` または `--radius-het2het` に `0` を渡すと、抽出器により内部的に `0.001 Å` にクランプされます。
- ダイアグラムのエネルギーは最初の状態（反応物）を基準とした kcal/mol（Hartree から変換）でプロットされます。
- 電荷処理: 抽出器の最初のモデルのポケット総電荷が、経路/スキャン/TSOPT の総電荷として使用されます（整数に丸め）。
- `-c/--center` を省略すると抽出をスキップし、入力構造全体を MEP/tsopt/freq/DFT ステージに直接渡します。単一構造実行では `--scan-lists` または `--tsopt` が依然として必要です。

---

## 関連項目

- [extract](extract.md) -- 単独のポケット抽出（`all` が内部で呼び出し）
- [mm_parm](mm_parm.md) -- AMBER トポロジー構築（`all` が内部で呼び出し）
- [path-search](path_search.md) -- 単独の再帰的 MEP 探索
- [tsopt](tsopt.md) -- 単独の TS 最適化
- [freq](freq.md) -- 振動解析と熱化学
- [dft](dft.md) -- DFT 一点計算
- [trj2fig](trj2fig.md) -- 軌跡からエネルギープロファイルをプロット
- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- よくあるエラーと対処法
- [YAML リファレンス](yaml_reference.md) -- 完全な YAML 設定オプション
- [用語集](glossary.md) -- MEP、TS、IRC、GSM の定義
