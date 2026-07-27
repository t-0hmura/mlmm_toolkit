# `all`

`mlmm all` は全系レイヤード PDB 上で ML/MM 機構パイプライン全体を 1 コマンドで実行します。内部では、活性部位抽出、MM トポロジー準備、ML/MM レイヤー割り当て、任意の段階的スキャン、MEP 探索（デフォルトは単一パス `path-opt`、`--refine-path` で再帰 `path-search`）、任意の後処理（TS 最適化、EulerPC IRC、熱化学、DFT 一点計算、DFT//MLIP/MM ダイアグラム）を順に管理します。ML 領域のデフォルト MLIP バックエンドは UMA で、`-b/--backend` で他のバックエンドを選択できます。

この順序は `all` が内部管理するステージを示します。単独で再利用するファイルを
作る場合は、`mm-parm` に明示的な出力接頭辞を与え、その出力 PDB に対して
`extract` と `define-layer` を実行します。この PDB は生成された parm7 と原子の
同一性・順序が一致し、空の元素記号列は `mm-parm` が補完します。

```bash
mlmm mm-parm -i input.pdb -l 'LIG:0' --out-prefix system
mlmm extract -i system.pdb -c LIG -l 'LIG:0' -o model.pdb
mlmm define-layer -i system.pdb --model-pdb model.pdb -o system_layered.pdb
```

`all` は渡す入力に応じて 3 つのモードのいずれかで動作します:

- **マルチ構造アンサンブル** -- 反応順に 2 つ以上の完全 PDB を提供し、複数構造にまたがる GSM（デフォルト）または DMF の MEP 探索を実行する。
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

ORB バックエンド:

```bash
mlmm all -i R.pdb P.pdb -c "SAM,GPP" -l "SAM:1,GPP:-3" \
 --backend orb --out-dir ./result_all_orb
```

CPU 実装の DMF（DMF バックエンドのデフォルトは GPU）:

```bash
mlmm all -i R.pdb P.pdb -c "SAM,GPP" -l "SAM:1,GPP:-3" \
 --mep-mode dmf --dmf-backend cpu --out-dir ./result_all_dmf
```

trajectory/structureを変換したPDB companionはテンプレートが利用可能な場合に生成され、
`--convert-files/--no-convert-files`（デフォルト有効）で制御されます。後述する確認用の
ML領域PDB pairは別の成果物であり、PDB入力では常に出力されます。

## 処理の流れ

1. **活性部位抽出と ML 領域定義**（複数入力時はマルチ構造の和集合）
   - 基質を定義します（`-c/--center`、PDB、残基 ID、または残基名で指定）。
   - 任意で `--ligand-charge` を総数値（分配）またはマッピング（例: `GPP:-3,MMT:-1`）として提供します。
   - 抽出器は入力ごとのポケット PDB を `<out-dir>/_work/pockets/` に書き出します。最初のポケットが `<out-dir>/ml_region.pdb`（`--model-pdb` として再利用可能な成果物）としてコピーされ、後続の全 ML/MM 計算の ML 領域を定義します。
   - `<out-dir>/ml_region_without_linkH.xyz` と `ml_region_with_linkH.xyz` に、リンク H 挿入前後のモデル系を出力します。PDB 入力では対応する `.pdb` companion も出力します。自動リンクペアは ML/MM 選択を横切る parm7 結合から決まり、距離による結合認識は行いません。
   - 抽出器の**最初のモデルの ML 領域の総電荷**が後続ステップの総電荷として使用され、丸め処理が発生した場合はコンソールに通知されます。
   - `-c/--center` を省略した場合は抽出をスキップし、完全入力構造をそのまま使用します。

2. **ML/MM 準備（parm7 + レイヤー割り当て）**
   - 最初の完全入力 PDB に対して `mm_parm` を一度実行し、`<out-dir>/mm_parm/<input_basename>.parm7` / `.rst7`（`--parm` として再利用可能な成果物）を構築します。これは自動的に `--parm` として渡されます。
   - 各完全系 PDB に対して `define-layer` を実行し、ML 領域定義に基づく 3 層 B 因子（ML=0.0、MovableMM=10.0、FrozenMM=20.0）を付与します。レイヤード全系 PDB は `<out-dir>/layered/` に書き出されます。

3. **任意の段階的スキャン（単一入力のみ）**
   - 完全入力 PDB が 1 つのみで `--scan-lists` が指定された場合、レイヤード全系 PDB に対して ML/MM calculatorを使用した段階的な結合距離駆動スキャンを実行します。
   - 各ステージの最終緩和構造（`stage_XX/result.pdb`）が中間体/生成物候補として収集されます。
   - 経路探索の入力系列は `[初期レイヤード PDB, stage_01/result.pdb, stage_02/result.pdb,...]` となります。

4. **全系レイヤード PDB での MEP 探索**
   - すべての MEP 計算は全系レイヤード PDB（`--parm` + `--detect-layer`）上で実行されます（ポケット上ではありません）。
   - **`--refine-path`:** 自動精密化を含む再帰的 `path_search` を実行し、多段階反応を自動検出して各素反応の詳細な MEP を構築します。複雑な多段階反応では手動での試行錯誤が必要な場合があります。両モードとも Stage 5 後処理に対応。
   - `--mep-mode` で GSM（デフォルト）または DMF を選択します。`--dmf-backend gpu` は `dmf.torch`、`--dmf-backend cpu` は NumPy 実装を使用します。GPU メモリ不足時は CPU を選択してください。
   - **`--no-refine-path`（デフォルト）:** 隣接ペアごとに選択した最適化法で単一パス `path-opt` を実行後、軌跡を結合、セグメントごとの HEI 抽出、結合変化検出、`summary.json` 書き出しまで行い、Stage 5 後処理（TSOPT、thermo、DFT）が利用可能。
   - マルチ入力実行では、元の完全 PDB がマージ参照として自動的に供給されます。スキャン由来の系列（単一構造の場合）では、元の完全 PDB 1 つがすべての入力の参照テンプレートとして再利用されます。

5. **サマリーと任意の後処理**
   - MEP エンジン生出力（セグメントごとの軌跡、全 MEP 軌跡、エンジンの `summary.json`）は `<out-dir>/_work/path_opt/`（`--refine-path` 使用時は `<out-dir>/_work/path_search/`）に書き出され、マージ済み成果物（`mep.pdb`、bridge 入力時の `mep.cif`、`mep_trj.xyz`、`mep_plot.png`、`energy_diagram_MEP.png`）は `<out-dir>/` へ移動され、`summary.{json,log}` はコピーされます。
   - `--tsopt`: 各 HEI で TS を最適化し、EulerPC IRC を実行し、セグメントエネルギーダイアグラムを描画します。
   - `--thermo`: (R, TS, P) で ML/MM 熱化学を計算し、Gibbs ダイアグラムを追加します。
   - `--dft`: (R, TS, P) のモデル領域で DFT 一点計算を実行し、モデル DFT 電子エネルギーダイアグラムを追加します。`--thermo` と組み合わせると、subtractive DFT//MLIP/MM 全エネルギーに ML/MM 熱補正を加えた DFT//MLIP/MM Gibbs ダイアグラムも生成されます。
   - `--tr-projection` は TS 最適化、IRC、振動解析、flatten PHVA に転送されます。デフォルトの `constrained` は凍結 anchor を動かさない全系剛体運動だけを除去し、実用的な ML/MM 境界では有効 rank は通常 0 です。
   - VRAM に余裕がある場合は `--hessian-calc-mode` を `Analytical` に設定することを強く推奨します（デフォルトの FiniteDifference より優先）。

6. **TSOPT のみモード**（単一入力、`--tsopt`、`--scan-lists` なし）
   - ステップ (4)-(5) をスキップし、レイヤード全系 PDB で `tsopt` を実行し、EulerPC IRC と両端の極小化を行い、R-TS-P の ML/MM エネルギーダイアグラムを構築し、任意で Gibbs、DFT、DFT//MLIP/MM ダイアグラムを追加します。
   - このモードでのみ、**より高いエネルギー**の IRC 端点が反応物 (R) として採用されます。

## 出力

ツリーは 3 つのゾーンで構成されます: **ルート直下の成果物**、**`segments/seg_NN/` 配下のセグメント別成果物**、**`_work/` 配下のパイプライン作業領域**（結果を取り出したあとは削除して構いません）。最初に確認する 3 つは `summary.log`、`summary.json`、`mep.pdb`（連結した反応経路。ルートへ移動）です。CIF/mmCIF bridge 入力では、元の識別子を復元した `mep.cif` もルートへ移動します。

```text
<out-dir>/
 summary.json                          # トップレベルサマリーのミラー（MEP ステージ実行時）
 summary.log
 mep.pdb · mep.cif                     # CIF は bridge 入力で元の ID を復元
 mep_trj.xyz
 mep_plot.png                          # MEP 生エネルギープロファイル
 energy_diagram_MEP.png                # 全セグメント MEP 障壁
 energy_diagram_MLIP_all.png           # 集約後処理ダイアグラム（有効時）
 energy_diagram_G_MLIP_all.png
 energy_diagram_DFT_all.png
 energy_diagram_G_DFT_plus_MLIP_all.png
 irc_plot_all.png
 ml_region.pdb                         # ML 領域定義（--model-pdb として再利用可能）
 ml_region_without_linkH.xyz           # リンク H 挿入前の ML モデル
 ml_region_with_linkH.xyz              # parm7 結合由来リンク H 挿入後の ML モデル
 ml_region_without_linkH.pdb           # PDB 入力時の topology 付き companion
 ml_region_with_linkH.pdb              # 生成した HL/LKH を含む PDB companion
 mm_parm/<input1>.parm7,.rst7          # 最初の完全酵素入力 PDB から生成した MM トポロジー（--parm として再利用可能）
 layered/                              # レイヤード全系 PDB（B 因子アノテーション付き、再利用可能な入力）
 segments/                             # 反応セグメント別の成果物
  seg_NN/                              # 2 桁インデックス (1 始まり)、例: seg_01, seg_02
   reactant.{pdb,cif} · ts.{pdb,cif} · product.{pdb,cif} # CIF は bridge 入力時
   ts/...                              # TS 最適化 + EulerPC IRC（--tsopt）
   irc/...
   freq/...                            # --thermo の場合
   dft/...                             # --dft の場合
   structures/{reactant,ts,product}.pdb  # 入れ子コピー + 生 IRC 端点
   energy_diagram_{MLIP,G_MLIP,DFT,G_DFT_plus_MLIP}.png
 _work/                               # パイプライン作業領域（削除可）
  pockets/                             # 入力ごとのポケット PDB（複数構造は統合）
  scan/                                # 単一構造+スキャンモードの場合のみ（stage_01/result.pdb …）
  path_opt/                            # MEP エンジン生出力（--refine-path 時は path_search/）
   summary.{json,log} · seg_NN_mep/    # セグメント別の生 MEP 軌跡（マージ済み成果物はルートへ移動）
```

**TSOPT のみモード**（単一入力 + `--tsopt`、`--scan-lists` なし）では MEP ステージが無く、最適化済み R/TS/P と `ts/`・`irc/`・`freq/`・`dft/` は `segments/seg_01/` 配下に生成され、`_work/path_opt/` は存在しません。

`-v 2` ではコンソールに抽出、MM 準備、スキャンステージ、MEP の進捗、ステージごとの所要時間が要約されます。{ref}`ja-verbosity-levels` を参照してください。

### `summary.log` の読み方
ログは番号付きセクションで構成されています:
- **[1] グローバル MEP 概要** -- イメージ/セグメント数、MEP 軌跡プロットパス、集約 MEP エネルギーダイアグラム。
- **[2] セグメントレベル MEP サマリー（MLIP 経路）** -- セグメントごとの障壁、反応エネルギー、結合変化サマリー。
- **[3] セグメントごとの後処理（TSOPT / Thermo / DFT）** -- セグメントごとの TS 虚振動数チェック、IRC 出力、エネルギーテーブル。
- **[4] エネルギーダイアグラム（概要）** -- MEP/MLIP/Gibbs/DFT シリーズのダイアグラムテーブルと任意のクロスメソッドサマリーテーブル。
- **[5] 出力ディレクトリ構造** -- インラインアノテーション付きの生成ファイルのコンパクトツリー。

### `summary.json` の読み方
summary.json はコンパクトな機械可読サマリーです。主なトップレベルキー:
- `out_dir`、`n_images`、`n_segments` -- 実行メタデータと総数。
- `segments` -- `index`、`tag`、`kind`、`barrier_kcal`、`delta_kcal`、`bond_changes` を持つセグメントごとのエントリリスト。
- `energy_diagrams`（任意）-- `labels`、`energies_kcal`、`energies_au`、`ylabel`、`image` パスを持つダイアグラムペイロード。

stage の `result.json` または `thermoanalysis.yaml` が書き出される場合、
`rigid_projection` ブロックに treatment、有効 rank、Hessian source、Hessian shape が
記録されます。全原子凍結はアクティブ自由度が残らないためエラーになります。

## CLI オプション

> **注意:** 表示されるデフォルト値はオプション未指定時に使用されます。完全なフラグ一覧は生成された command reference（`reference/commands/`）にあり、以下の表は説明が必要なオプションを扱います。

### 入出力オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順の 2 つ以上の完全 PDB（`--scan-lists`（段階的スキャン）または `--tsopt`（TSOPT のみ）の場合のみ単一入力可）。 | 必須 |
| `-c, --center TEXT` | 基質指定（PDB パス、残基 ID（`308,309`）、または残基名（`SAM,GPP`））。省略時は抽出をスキップし完全構造をそのまま使用。 | _None_ |
| `-l, --ligand-charge TEXT` | 非標準残基の総電荷または残基別マッピング（例: `GPP:-3,MMT:-1`）。 | _None_ |
| `-q, --charge INT` | ML 領域/model system の正味電荷を強制指定（最優先の上書き）。 | _None_ |
| `-o, --out-dir PATH` | トップレベル出力ディレクトリ。 | `./result_all/` |
| `--parm FILE` | 全系の AMBER parm7 トポロジーファイル。省略時は `mm_parm` で自動生成。 | _None_ |
| `--model-pdb FILE` | 構築済み ML 領域 PDB。指定時は ML 領域決定をスキップし、このファイルで ML 領域を直接定義。 | _None_ |
| `--ref-pdb FILE` | XYZ 入力用の参照 PDB。入力が XYZ の場合に PDB メタデータ（残基、鎖、B 因子）を復元するために必要。 | _None_ |
| `--convert-files/--no-convert-files` | テンプレート利用可能時に XYZ/TRJ から対応する PDB の生成を切り替えるグローバルトグル。 | `True` |
| `--dump/--no-dump` | 任意のオプティマイザ軌跡・リスタートを保存。常に `path-search`/`path-opt` に転送し、`scan`/`tsopt` にはここで明示設定時のみ転送します。`--thermo` 時は Gibbs 集約を欠損させないため、必須の子 `thermoanalysis.yaml` handoff を `--no-dump` でも保持します。 | `False` |
| `--config FILE` | 先に適用するベース YAML。 | _None_ |
| `--show-config/--no-show-config` | 実行前に解決済み設定を表示。 | `False` |
| `--dry-run/--no-dry-run` | 一時ディレクトリで抽出/setup と電荷・parity 検証を実行し、計画を表示して計算 stage は省略。`--help-advanced` に表示。 | `False` |

### 抽出オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-r, --radius FLOAT` | ポケット包含カットオフ (Å)。 | `2.6` |
| `--radius-het2het FLOAT` | 独立したヘテロ-ヘテロカットオフ (Å)。 | `0.0` |
| `--include-h2o/--no-include-h2o` | 水分子（HOH/WAT/H2O/DOD/TIP/TIP3/SOL）を含める。 | `True` |
| `--exclude-backbone/--no-exclude-backbone` | 非基質アミノ酸の主鎖原子を除去。 | `False` |
| `--add-linkh/--no-add-linkh` | 切断結合にリンク水素を付加。 | `False` |
| `--selected-resn TEXT` | 強制包含する残基。 | `""` |
| `--modified-residue TEXT` | 修飾アミノ酸残基名をカンマ区切りで指定（任意で電荷付き）。主鎖切断と電荷計算にアミノ酸として扱う。例: `HD1,HD2,HD3` または `HD1:0,SEP:-2`。 | `""` |

### MM 準備オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--auto-mm-ff-set {ff19SB\|ff14SB}` | `mm_parm` 用の力場セット（ff19SB は OPC3、ff14SB は TIP3P を使用）。 | `ff19SB` |
| `--auto-mm-add-ter/--auto-mm-no-add-ter` | リガンド/水/イオンブロック周囲の TER 挿入を制御。 | `True` |
| `--auto-mm-disulfide/--auto-mm-no-disulfide` | mm_parm に転送：CYS/CYM/CYX にわたり SG-SG 幾何からジスルフィドを検出して結合（結合された CYS は CYX にリネーム）。`--auto-mm-no-disulfide` では既に CYX の残基のみを結合。 | `True` |
| `--auto-mm-keep-temp` | `mm_parm` の一時作業ディレクトリを保持（デバッグ用）。 | `False` |
| `--auto-mm-ligand-mult TEXT` | `mm_parm` に転送するスピン多重度マッピング（例: `GPP:2,SAM:1`）。省略時は全リガンドに 1 を使用。 | _None_ |

### MEP 探索オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `-b, --backend CHOICE` | ML バックエンド: `uma`（デフォルト）、`orb`、`mace`、`aimnet2`。全計算サブコマンドに転送。 | `uma` |
| `--embedcharge/--no-embedcharge` | v0.3.3 では使用不可。旧コマンドを明示的に拒否するためにのみ残されています。 | `False` |
| `--embedcharge-cutoff FLOAT` | 廃止した電子埋め込み経路とともに使用不可。 | — |
| `--cmap/--no-cmap` | model parm7 に CMAP（骨格クロスマップ二面角補正）を含めるかどうか。デフォルト: 無効（Gaussian ONIOM と同一）。 | `--no-cmap` |
| `--mep-mode [gsm\|dmf]` | `path-opt` と再帰的 `path-search` の両方へ転送する MEP 最適化法。 | `gsm` |
| `--dmf-backend [gpu\|cpu]` | DMF 実装。明示指定時だけ子コマンドへ転送するため、省略時は子コマンドの YAML 設定 `dmf.backend` が有効。 | `gpu` |
| `--max-nodes INT` | GSM/DMF セグメントの内部ノード数。 | `20` |
| `--max-cycles INT` | MEP 最適化サイクルの最大数。 | `300` |
| `--climb/--no-climb` | 選択した最適化法が対応する場合に climbing-image TS 精密化を有効化。 | `True` |
| `--opt-mode [grad\|hess]` | スキャン/path-search と単一構造最適化のプリセット（`grad` → L-BFGS/Dimer、`hess` → RFO/RSIRFO）。 | `grad` |
| `--opt-mode-post [grad\|hess]` | TSOPT/IRC 後端点最適化向けのプリセット上書き（`grad` → Dimer/L-BFGS、`hess` → RS-I-RFO/RFO）。 | `hess` |
| `--thresh TEXT` | 収束プリセット（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`、`never`）。実効デフォルト: path-opt は `gau_loose`、scan は `gau`。 | _None_ |
| `--thresh-post TEXT` | IRC 後端点最適化の収束プリセット。 | `baker` |
| `--preopt/--no-preopt` | セグメント化前に端点を事前最適化。 | `True` |
| `--refine-path/--no-refine-path` | `--no-refine-path`（デフォルト）= 単一パス `path-opt`（軌跡結合 + HEI 抽出 + 結合変化検出 + `summary.json`）、`--refine-path` = 再帰的 `path-search`。どちらも `--mep-mode` の選択と Stage 5（TSOPT/thermo/DFT）に対応。 | `False` |
| `--hessian-calc-mode CHOICE` | ML/MM Hessian モード（`Analytical` または `FiniteDifference`）。 | `FiniteDifference` |
| `--precision [fp32\|fp64]` | バックエンド精度。省略時は UMA/AIMNet2 fp32、ORB/MACE fp64。AIMNet2 は fp64 を拒否。 | バックエンド依存 |
| `--workers INT` | UMA predictor worker 数。2 以上は `fairchem-core[extras]` が必要で、解析 Hessian と併用不可。 | `1` |
| `--workers-per-node INT` | UMA 並列 predictor のノード当たり worker 数。 | _None_ |
| `--detect-layer/--no-detect-layer` | 入力 PDB の B 因子（B=0/10/20）から ML/MM レイヤーを検出。無効時は下流ツールで `--model-pdb` または `--model-indices` が必要。 | `True` |

TSOPT の最適化モード選択順: `--opt-mode-post`（設定時）-> `--opt-mode`（明示指定時のみ）-> TSOPT デフォルト（`hess` → RS-I-RFO）。

### スキャンオプション（単一入力実行）

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-s, --scan-lists TEXT...` | 段階的スキャン: `(i,j,target_A)` タプル。 | _None_ |
| `--scan-out-dir PATH` | スキャン出力ディレクトリの上書き。 | _None_ |
| `--scan-one-based/--scan-zero-based` | スキャン原子インデックスを 1 始まりまたは 0 始まりとして解釈。 | _None_ |
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
| `--reject-uphill/--no-reject-uphill` | IRC 後の**エンドポイント再最適化のみ**で RFO の上り坂ステップを拒否（opt 子へ転送。低エネルギー形状へロールバックし trust radius を縮小）。TS 最適化では拒否を常に無効化し、経路探索には影響しない。emergency floor 到達時は、保持したエンドポイントを通常の収束条件で最終確認。 | `True` |
| `--tr-projection [constrained\|legacy-active]` | 凍結境界 TR 処理を `tsopt`、`irc`、`freq`、flatten PHVA へ転送。`legacy-active` は非推奨の比較専用で、pass/HOSP 遷移状態認定には使用不可。 | `constrained` |
| `--irc-step-size FLOAT` | TS 後の各 IRC に EulerPC 最大ステップ（Bohr）を転送。数フレームで停止する場合は `0.05` など小さい値で再試行。 | IRC デフォルト `0.10` |
| `--irc-never-stop/--no-irc-never-stop` | エネルギー上昇/plateau 停止だけを無視して IRC を継続。収束、非有限値、サイクル上限では停止。 | `False` |
| `--tsopt-max-cycles INT` | `tsopt --max-cycles` の上書き。 | _デフォルト_ |
| `--tsopt-out-dir PATH` | tsopt サブディレクトリのカスタマイズ。 | _None_ |
| `--freq-out-dir PATH` | freq 出力ディレクトリの上書き。 | _None_ |
| `--freq-max-write INT` | 出力する最大モード数。 | _デフォルト_ |
| `--freq-amplitude-ang FLOAT` | モードアニメーション振幅 (Å)。 | _デフォルト_ |
| `--freq-n-frames INT` | モードアニメーションのフレーム数。 | _デフォルト_ |
| `--freq-sort TEXT` | モードソート方法。 | _デフォルト_ |
| `--freq-temperature FLOAT` | 熱化学温度 (K)。 | _デフォルト_ |
| `--freq-pressure FLOAT` | 熱化学圧力 (atm)。 | _デフォルト_ |
| `--freq-symmetry-number INT` | R/TS/P の全 freq 計算に共通の回転対称数。省略時は各子計算の YAML/デフォルトに従う。 | _None_ |
| `--dft-out-dir PATH` | DFT 出力ディレクトリの上書き。 | _None_ |
| `--dft-func-basis TEXT` | 汎関数/基底関数ペア。 | _デフォルト_ |
| `--dft-max-cycle INT` | 最大 SCF 反復数。 | _デフォルト_ |
| `--dft-conv-tol FLOAT` | SCF 収束閾値。 | _デフォルト_ |
| `--dft-grid-level INT` | PySCF グリッドレベル。 | _デフォルト_ |
| `--dft-engine [gpu\|cpu]` | DFT エンジン（GPU or CPU PySCF）。 | _None_ |

## YAML 設定

`all` は YAML 設定をサポートします:

- `--config FILE`: ベース設定。

`defaults < config < 明示指定 CLI`

解決後の YAML が下流サブコマンドへ転送されます。各ツールは独自ドキュメントに記載されたセクションを読み取ります:

| サブコマンド | YAML セクション |
|------------|---------------|
| [`path-search`](path-search.md) | `geom`, `calc`/`mlmm`, `gs`, `opt`, `lbfgs`, `bond`, `search` |
| [`scan`](scan.md) | `geom`, `calc`/`mlmm`, `opt`, `lbfgs` |
| [`tsopt`](tsopt.md) | `geom`, `calc`/`mlmm`, `opt`, `hessian_dimer`, `rsirfo` |
| [`freq`](freq.md) | `geom`, `calc`/`mlmm`, `freq`, `thermo` |
| [`dft`](dft.md) | `dft` |

明示指定した CLI 値だけが `--config` の値を上書きします。

**最小の YAML 例:**
```yaml
geom:
 tr_projection: constrained      # legacy-active は非推奨・比較専用
calc:
 charge: 0
 spin: 1
 real_parm7: real.parm7
 model_pdb: ml_region.pdb
 backend: uma                    # ML バックエンド (uma/orb/mace/aimnet2)
 embedcharge: false              # 互換性用。true は拒否される
 uma_model: uma-s-1p2            # uma-s-1p2 | uma-m-1p1
 hessian_calc_mode: Analytical     # VRAM に余裕がある場合に推奨
gs:
 max_nodes: 20
 climb: true
dft:
 grid_level: 6
```

すべての YAML オプションの完全なリファレンスは **[YAML 設定リファレンス](yaml-reference.md)** を参照してください。

`--tr-projection` と `tsopt --ref-mode` は別の機能です。前者は凍結境界の剛体モード、
後者は鞍点回復で使う内部的な MEP 接線 handoff を制御します。

## 注記

入力形式は抽出の有無に依存します:

- 抽出有効時（`-c/--center`）: 入力は残基を特定するため **PDB** ファイルが必要。
- 抽出スキップ時: 入力は **PDB/XYZ** が使用可能。
- マルチ構造実行には 2 つ以上の構造が必要。

電荷は優先度の高い順に解決されます -- `-q/--charge`（明示的な CLI 上書き）-> ポケット抽出（`-c` 指定時、アミノ酸 + イオン + `--ligand-charge` の合計）-> `-l, --ligand-charge` フォールバック（抽出スキップ時）-> デフォルト（未解決の電荷はエラー）。スピンの解決: `--multiplicity`（CLI）-> デフォルト（1）。正しい電荷伝播のため、非標準基質には常に `--ligand-charge` を指定してください。最初のモデルの ML 領域の総電荷は最も近い整数に丸められ、丸め処理が発生した場合はコンソールに通知されます。

## 関連項目

- [extract](extract.md) -- 単独の ML 領域決定（`all` が内部で呼び出し）
- [mm-parm](mm-parm.md) -- AMBER トポロジー構築（`all` が内部で呼び出し）
- [path-search](path-search.md) -- 単独の再帰的 MEP 探索
- [tsopt](tsopt.md) -- 単独の TS 最適化
- [freq](freq.md) -- 振動解析と熱化学
- [dft](dft.md) -- DFT 一点計算
- [trj2fig](trj2fig.md) -- 軌跡からエネルギープロファイルをプロット
- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- よくあるエラーと対処法
- [YAML リファレンス](yaml-reference.md) -- 完全な YAML 設定オプション
- [用語集](glossary.md) -- MEP、TS、IRC、GSM の定義
