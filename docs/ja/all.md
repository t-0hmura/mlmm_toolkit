# `all`

## 概要

> **要約:** エンドツーエンドの酵素反応ワークフロー -- ポケット抽出、任意の段階的スキャン、再帰的 MEP 探索（GSM）、全系へのマージ、任意の TS 最適化、擬似 IRC、熱化学、DFT、DFT//UMA ダイアグラム。

`mlmm all` はポケットモデルを中心としたワンショットパイプラインを実行します。以下の 3 つのモードをサポートしています:

- **マルチ構造アンサンブル** -- 反応順に 2 つ以上の完全 PDB を提供します。ポケット抽出、MM トポロジー構築、再帰的 GSM MEP 探索、全系へのマージを行い、任意でセグメントごとの後処理（TSOPT/freq/DFT）を実行します。
- **単一構造 + 段階的スキャン** -- 1 つの PDB と `--scan-lists` を提供します。スキャンで中間体/生成物候補を生成し、MEP の端点として使用します。
- **TSOPT のみ** -- 1 つの PDB を提供し、`--scan-lists` を省略して `--tsopt` を設定します。ポケットで TS 最適化を実行し、擬似 IRC、両端の極小化、エネルギーダイアグラムの構築を行います。

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

## 使用法

```bash
# 標準（反応順のマルチ構造アンサンブル）
mlmm all -i R.pdb [I1.pdb...] P.pdb -c <substrate-spec> [--ligand-charge <map-or-number>]
 [--multiplicity <2S+1>] [--max-nodes N] [--max-cycles N]
 [--climb/--no-climb] [--sopt-mode lbfgs|rfo|light|heavy]
 [--opt-mode light|lbfgs|heavy|rfo]
 [--pre-opt/--no-pre-opt] [--hessian-calc-mode Analytical|FiniteDifference] [--out-dir DIR]
 [--tsopt/--no-tsopt] [--thermo/--no-thermo] [--dft/--no-dft]
 [--tsopt-max-cycles N] [--freq-* overrides] [--dft-* overrides]

# 単一構造 + 段階的スキャン
mlmm all -i A.pdb -c "308,309" --scan-lists "[(12,45,1.35)]" "[(10,55,2.20)]" \
 --multiplicity 1 --sopt-mode lbfgs --pre-opt \
 --out-dir result_all --tsopt --thermo --dft

# 単一構造 TSOPT のみモード（経路探索なし）
mlmm all -i single.pdb -c "GPP,MMT" --ligand-charge "GPP:-3,MMT:-1" \
 --tsopt --thermo --dft --out-dir result_tsopt_single
```

### 例

```bash
# 明示的な基質とリガンド電荷によるミニマルなエンドツーエンド実行（マルチ構造）
mlmm all -i reactant.pdb product.pdb -c "GPP,MMT" --ligand-charge "GPP:-3,MMT:-1"

# 中間体付きアンサンブル、残基 ID 基質指定、完全な後処理
mlmm all -i A.pdb B.pdb C.pdb -c "308,309" --ligand-charge "-1" \
 --multiplicity 1 --max-nodes 10 --max-cycles 100 --climb \
 --sopt-mode lbfgs --no-dump --config params.yaml --pre-opt \
 --out-dir result_all --tsopt --thermo --dft

# 単一構造 + スキャンで順序付き系列を構築
mlmm all -i A.pdb -c "308,309" --scan-lists "[(10,55,2.20),(23,34,1.80)]" \
 --multiplicity 1 --out-dir result_scan_all --tsopt --thermo --dft

# 単一構造 TSOPT のみモード（path_search なし）
mlmm all -i A.pdb -c "GPP,MMT" --ligand-charge "GPP:-3,MMT:-1" \
 --tsopt --thermo --dft --out-dir result_tsopt_only
```

## ワークフロー

1. **活性部位ポケット抽出**（複数入力時はマルチ構造の和集合）
 - 基質を定義します（`-c/--center`、PDB、残基 ID、または残基名で指定）。
 - 任意で `--ligand-charge` を総数値（分配）またはマッピング（例: `GPP:-3,MMT:-1`）として提供します。
 - 抽出器は入力ごとのポケット PDB を `<out-dir>/pockets/` に書き出します。
 - 抽出器の**最初のモデルのポケット総電荷**が後続ステップの総電荷として使用され、丸め処理が発生した場合はコンソールに通知されます。
 - 追加の抽出トグル: `--radius`、`--radius-het2het`、`--include-H2O/--no-include-H2O`、`--exclude-backbone/--no-exclude-backbone`、`--add-linkH/--no-add-linkH`、`--selected-resn`、`--verbose/--no-verbose`。
 - `-c/--center` を省略した場合は抽出をスキップし、完全入力構造をそのまま使用します。

### 電荷の解決順序
1. `-q/--charge`（明示的な CLI 上書き）— 最優先。
2. ポケット抽出の電荷サマリー（アミノ酸、イオン、`--ligand-charge` の合計）。
3. 抽出スキップ時の `--ligand-charge` フォールバック。
4. デフォルト: なし（未解決なら中断）。

2. **ML/MM 準備**
 - 最初のポケットを `<out-dir>/ml_region.pdb` として `--model-pdb` に使用します。リンク水素をこの定義から除外したい場合は、抽出時に `--no-add-linkH`（デフォルト）を保持します。
 - 最初の完全入力 PDB に対して `mm_parm` を一度実行し、`<out-dir>/mm_parm/<input_basename>.parm7` / `.rst7` を構築します。これは自動的に `--real-parm7` として渡されます。
 - このステージは `--auto-mm-ff-set`、`--auto-mm-add-ter`、`--auto-mm-keep-temp` で調整できます。

3. **任意の段階的スキャン（単一入力のみ）**
 - 完全入力 PDB が 1 つのみで `--scan-lists` が指定された場合、抽出されたポケット PDB に対して UMA 計算機を使用した段階的な結合距離駆動スキャンを実行します。
 - 各ステージの最終緩和構造（`stage_XX/result.pdb`）が中間体/生成物候補として収集されます。
 - 経路探索の入力系列は `[初期ポケット, stage_01/result.pdb, stage_02/result.pdb,...]` となります。

4. **ポケット入力での MEP 探索（再帰的 GSM）**
 - このコマンドから転送されたオプションで `path_search` を実行します。
 - マルチ入力実行では、元の完全 PDB がマージ参照として自動的に供給されます。スキャン由来の系列（単一構造の場合）では、単一の元の完全 PDB がすべてのポケット入力の参照テンプレートとして再利用（繰り返し）されます。

5. **全系へのマージと任意の後処理**
 - ポケット MEP は `<out-dir>/path_search/` 内で元の全系テンプレートにマージされます。
 - ポケットのみと全系の軌跡、セグメントごとのマージ PDB、サマリーが書き出されます。
 - `--tsopt`: HEI ポケットで TS を最適化し、虚数モードに沿った変位による擬似 IRC を実行し、結合状態マッチングにより正方向/逆方向の対応を割り当て、セグメントダイアグラムを描画します。
 - `--thermo`: (R, TS, P) で UMA 熱化学を計算し、Gibbs ダイアグラムを追加します。
 - `--dft`: (R, TS, P) で DFT 一点計算を実行し、DFT ダイアグラムを追加します。`--thermo` と組み合わせると、DFT//UMA Gibbs ダイアグラムも生成されます。

6. **TSOPT のみモード**（単一入力、`--tsopt`、`--scan-lists` なし）
 - ステップ (4)-(5) をスキップし、ポケットで `tsopt` を実行し、擬似 IRC と両端の極小化を行い、R-TS-P の UMA エネルギーダイアグラムを構築し、任意で UMA Gibbs、DFT、DFT//UMA ダイアグラムを追加します。
 - このモードでのみ、**より高いエネルギー**の IRC 端点が反応物 (R) として採用されます。

## CLI オプション

### 入出力オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH...` | 反応順の 2 つ以上の完全 PDB（`--scan-lists` または `--tsopt` の場合のみ単一入力可）。 | 必須 |
| `-c, --center TEXT` | 基質指定（PDB パス、残基 ID、または残基名）。省略時は抽出をスキップ。 | _None_ |
| `--ligand-charge TEXT` | 総電荷または残基別マッピング（例: `GPP:-3,MMT:-1`）。 | _None_ |
| `-q, --charge INT` | 総電荷を強制指定（最優先の上書き）。 | _None_ |
| `--out-dir PATH` | トップレベル出力ディレクトリ。 | `./result_all/` |
| `--dump/--no-dump` | オプティマイザーダンプを保存。 | `False` |
| `--config FILE` | 先に適用するベース YAML。 | _None_ |
| `--show-config/--no-show-config` | 実行前に解決済み設定を表示。 | `False` |
| `--dry-run/--no-dry-run` | 実行せず検証と計画表示のみ行う。 | `False` |

### 抽出オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--radius FLOAT` | ポケット包含カットオフ (Angstrom)。 | 抽出器デフォルト |
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
| `--max-nodes INT` | セグメント GSM の内部ノード数。 | `10` |
| `--max-cycles INT` | GSM マクロサイクルの最大数。 | `100` |
| `--climb/--no-climb` | セグメント GSM の TS 精密化を有効化。 | `True` |
| `--sopt-mode TEXT` | 単一構造オプティマイザープリセット（`lbfgs`、`rfo`、`light`、`heavy`）。 | _デフォルト_ |
| `--opt-mode TEXT` | スキャンおよび tsopt 用のオプティマイザープリセット（`light`、`lbfgs`、`heavy`、`rfo`）。 | _デフォルト_ |
| `--pre-opt/--no-pre-opt` | セグメント化前に端点を事前最適化。 | `True` |
| `--hessian-calc-mode CHOICE` | UMA ヘシアンモード（`Analytical` または `FiniteDifference`）。 | _デフォルト_ |

### スキャンオプション（単一入力実行）

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--scan-lists TEXT...` | 段階的スキャン: `(i,j,target_A)` タプル。 | _None_ |
| `--scan-out-dir PATH` | スキャン出力ディレクトリの上書き。 | _None_ |
| `--scan-one-based/--no-scan-one-based` | 原子セレクターに 1 始まりインデックスを強制。 | `True` |
| `--scan-max-step-size FLOAT` | 最大ステップサイズ (Angstrom)。 | _デフォルト_ |
| `--scan-bias-k FLOAT` | 調和バイアス強度 (eV/Angstrom^2)。 | _デフォルト_ |
| `--scan-relax-max-cycles INT` | ステップごとの緩和最大サイクル。 | _デフォルト_ |
| `--scan-preopt/--no-scan-preopt` | スキャン事前最適化トグルの上書き。 | _デフォルト_ |
| `--scan-endopt/--no-scan-endopt` | スキャンステージ終端最適化の上書き。 | _デフォルト_ |

### 後処理オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--tsopt/--no-tsopt` | 反応セグメントごとに TS 最適化 + 擬似 IRC を実行。 | `False` |
| `--thermo/--no-thermo` | R/TS/P で振動解析 (`freq`) を実行。 | `False` |
| `--dft/--no-dft` | R/TS/P で DFT 一点計算を実行。 | `False` |
| `--tsopt-max-cycles INT` | `tsopt --max-cycles` の上書き。 | _デフォルト_ |
| `--tsopt-out-dir PATH` | tsopt サブディレクトリのカスタマイズ。 | _None_ |

### Freq 上書き

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--freq-out-dir PATH` | freq 出力ディレクトリの上書き。 | _None_ |
| `--freq-max-write INT` | 出力する最大モード数。 | _デフォルト_ |
| `--freq-amplitude-ang FLOAT` | モードアニメーション振幅 (Angstrom)。 | _デフォルト_ |
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
 ml_region.pdb # ML 領域定義（最初のポケットのコピー）
 pockets/
 pocket_<input1_basename>.pdb
 pocket_<input2_basename>.pdb
...
 mm_parm/
 <input1_basename>.parm7 # 最初の完全酵素入力 PDB から生成
 <input1_basename>.rst7
 scan/ # 単一構造+スキャンモードの場合のみ
 stage_01/result.pdb
 stage_02/result.pdb
...
 summary.yaml # トップレベルサマリーのミラー（path_search 実行時）
 summary.log
 mep_plot.png
 energy_diagram_MEP.png
 energy_diagram_UMA_all.png # 集約後処理ダイアグラム（有効時）
 energy_diagram_G_UMA_all.png
 energy_diagram_DFT_all.png
 energy_diagram_G_DFT_plus_UMA_all.png
 irc_plot_all.png
 path_search/ # path_search 実行時
 mep_trj.xyz
 mep.pdb
 mep_w_ref.pdb
 mep_w_ref_seg_XX.pdb
 summary.yaml
 summary.log
 mep_plot.png
 energy_diagram_MEP.png
 post_seg_XX/ # 後処理有効時
 ts/...
 irc/...
 freq/... # --thermo の場合
 dft/... # --dft の場合
 energy_diagram_UMA.png
 energy_diagram_G_UMA.png
 energy_diagram_DFT.png
 energy_diagram_G_DFT_plus_UMA.png
 tsopt_single/ # 単一構造 TSOPT のみモードの場合
 ts/...
 irc/...
 structures/
 reactant.pdb
 ts.pdb
 product.pdb
 freq/... # --thermo の場合
 dft/... # --dft の場合
 energy_diagram_UMA.png
 energy_diagram_G_UMA.png
 energy_diagram_DFT.png
 energy_diagram_G_DFT_plus_UMA.png
```

## YAML 設定

`all` は 2 層の YAML を受け付けます:

- `--config FILE`: ベース設定。

適用順序:

`defaults < config < CLI < override-yaml`

解決後の YAML が下流サブコマンドへ転送されます。各ツールは独自ドキュメントに記載されたセクションを読み取ります:

| サブコマンド | YAML セクション |
|------------|---------------|
| [`path-search`](path_search.md) | `geom`, `calc`/`mlmm`, `gs`, `opt`, `lbfgs`, `bond`, `search` |
| [`scan`](scan.md) | `geom`, `calc`/`mlmm`, `opt`, `lbfgs` |
| [`tsopt`](tsopt.md) | `geom`, `calc`/`mlmm`, `opt` |
| [`freq`](freq.md) | `geom`, `calc`/`mlmm`, `freq` |
| [`dft`](dft.md) | `dft` |

## 注意事項

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- **Python >= 3.11** が必要です。
- **基質 (`-c/--center`) と（必要な場合は）`--ligand-charge` は実質的に必須です。**
- 単一構造モードでは `--scan-lists` または `--tsopt` が必要です。それ以外の場合は少なくとも 2 つの構造が必要です。
- preflight チェックで、繰り返し指定した `-i/--input` が実ファイルであること、および AmberTools コマンド（`tleap`, `antechamber`, `parmchk2`）の存在を実行前に検証します。
- ダイアグラムのエネルギーは最初の状態を基準とした kcal/mol（Hartree から変換）でプロットされます。
- 電荷処理: 抽出器の最初のモデルのポケット総電荷が、経路/スキャン/TSOPT の総電荷として使用されます（整数に丸め）。

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
