# 出力ディレクトリのレイアウト

このページでは、各 `mlmm` サブコマンドが出力ディレクトリに書き込むファイルと、エージェントや後段スクリプトが従うべき規約を説明します。

## ファイル名の規約

| ファイル名 | 書き込み元 | 用途 |
|---|---|---|
| `summary.json` | 集約結果の書き込み処理まで到達した `all` / `path-search` | 集約ワークフローの正規 JSON エンベロープ（[JSON 出力リファレンス](json-output.md)）。早期の CLI 引数または入力の検証では作られない場合があります。 |
| `summary.json` | 正常終了した段階別・レポート系コマンドでは `--out-json` 指定時（デフォルトは `--no-out-json`）。捕捉した実行時エラーでは、フラグなしでも可能な範囲でエラーエンベロープを書く場合があります | 個別結果の `result.json` と互換性のあるミラー。書き込み処理が正常終了した場合は同一のバイト列です。`fix-altloc`、`add-elem-info`、`bond-summary` などは出力しません。 |
| `result.json` | 段階別 `summary.json` と同条件（`opt`、`tsopt`、`freq`、`irc`、`sp`、scan 系、`path-opt`、`dft`、`extract`、`trj2fig`、`energy-diagram`） | 個別結果・レポートの正規エンベロープ。互換ミラーより後に公開されるため、中断した世代を判定するときはこちらを読みます。 |
| `summary.log` | `path-search`、`all` | 人が読むための実行ログ（セグメント / ステージごとに 1 行）。 |
| `final_geometry.xyz` | `opt`、`tsopt` | 最適化された構造（XYZ、フル精度）。 |
| `mep.pdb` / `mep.cif` / `mep_trj.xyz` | `path-search`、`all` | 反応経路のフレーム。bridge 入力では `mep.cif` が元の ID を復元します。単独実行の `path-opt` は代わりに `final_geometries_trj.xyz` / `final_geometries.pdb` を書き込みます。 |
| `ml_region_without_linkH.{xyz,pdb}` / `ml_region_with_linkH.{xyz,pdb}` | `all`、`dft` | parm7 の ML/MM 境界結合からリンク H を生成する前後の ML モデル。PDB companion は PDB 入力時に出力します。 |
| `mep_plot.png` | `path-search`、`all` | 生の MEP エネルギープロファイル（PNG）。`all` はエンジン出力からルートにコピーします。 |
| `forward_irc_trj.xyz` / `backward_irc_trj.xyz`（および `finished_irc_trj.xyz`） | `irc` | IRC 軌跡（XYZ）。対応する `*_irc.pdb` ファイルが同じフレームを PDB 形式で保持します。 |
| `frequencies_cm-1.txt` | `freq` | 振動数の一覧（cm⁻¹）。 |
| `*.gjf` | 各種（`--convert-files` 指定時） | Gaussian 形式の構造ファイル。 |

## デフォルトの `--out-dir`

| サブコマンド | デフォルトの `--out-dir` |
|---|---|
| `all` | `./result_all/` |
| `opt` | `./result_opt/` |
| `tsopt` | `./result_tsopt/` |
| `freq` | `./result_freq/` |
| `irc` | `./result_irc/` |
| `dft` | `./result_dft/` |
| `scan` / `scan2d` / `scan3d` | `./result_scan*/` |
| `path-opt` / `path-search` | `./result_path_*/` |
| `sp` | `./result_sp/` |
| `extract` | `./`（作業ディレクトリに `pocket.pdb` を書き込み。入力が複数の場合は `pocket_<input>.pdb`） |
| `mm-parm` | `./`（`<prefix>.parm7` / `<prefix>.rst7` を書き込み） |
| `define-layer` | `./`（`<input>_layered.pdb` を書き込み） |

`--out-dir <path>`（または `-o`）で上書きできます。明示的に指定したパスは、ステージ別デフォルトと YAML の両方に優先します。

## 単独実行と `all`

サブコマンドを単独で実行すると、**フラットな**結果ディレクトリが書き込まれます。同じライターでも `all` によってオーケストレーションされると、構造化されたツリーにネストされます。

- **単独サブコマンド** → 上記のファイルを含むフラットな `result_<subcmd>/`。`segments/` も `_work/` もありません。これらは `all` が 1 回の実行で複数のライターを協調させるときのみ現れます。
- **`all` の内部では、リーフライターはそのままネストされます。** `segments/seg_NN/<subcmd>/` のセグメント別リーフ出力は、単独の `result_<subcmd>/` と構造的に同一です。`all` はライターの出力先を別のディレクトリに向けているだけです。
- **`path-search` / `path-opt` はエンジン側の例外です。** 単独実行では `path-search` 自体が成果物となります（`result_path_search/` に独自の `summary.log`、`mep.pdb`、bridge 入力時の `mep.cif`、`mep_trj.xyz`、`mep_plot.png`、`energy_diagram_MEP.png` を持ちます）。`all` の内部では、その生の出力は `_work/path_opt/` 以下のエンジン用スクラッチであり（`--refine-path` 指定時のみ `_work/path_search/`）、マージされた成果物（`mep.pdb`、bridge 入力時の `mep.cif`、`mep_trj.xyz`、`mep_plot.png`、`energy_diagram_MEP.png`）はパイプラインのルートに移動され、`summary.{json,log}` がそこにコピーされます。この非対称性は意図的なものです。

したがって `all` のツリーには 3 つのゾーンがあります。

```text
result_all/
├─ summary.log · summary.json                 # ルートへコピー
├─ mep.pdb · mep.cif · mep_trj.xyz · mep_plot.png · energy_diagram_MEP.png
├─ energy_diagram_*_all.png · irc_plot_all.png
├─ ml_region.pdb                              # ML-region definition (reusable as --model-pdb)
├─ ml_region_without_linkH.{xyz,pdb} · ml_region_with_linkH.{xyz,pdb}
├─ mm_parm/                                   # MM topology <input>.parm7 / .rst7 (reusable as --parm)
├─ layered/                                   # layered full-system PDBs (B-factor annotated; reusable inputs)
├─ segments/
│  └─ seg_NN/                                  # 反応セグメント別の成果物（2桁番号）
│     ├─ reactant.{pdb,cif} · ts.{pdb,cif} · product.{pdb,cif} # CIF は bridge 入力時
│     └─ ts/ · irc/ · freq/ · dft/ · structures/    # 段階別の作業ファイル（--tsopt / --thermo / --dft）
└─ _work/                                      # パイプラインのスクラッチ（削除可）
   ├─ pockets/ · scan/
   └─ path_opt/                                # MEP エンジンの生出力（--refine-path 時は path_search/）
```

TSOPT のみのモードでは MEP ステージがないため、`_work/path_opt/` は存在せず、成果物は `segments/seg_01/` 以下に置かれます。モードごとの完全な内訳は [all](all.md) を参照してください。

## エージェント向けレシピ

```python
# 実行したコマンドに対応する正規のファイル名を選ぶ。
import json
from pathlib import Path

subcommand = "opt"  # 実行したコマンドに置き換える
primary = "summary.json" if subcommand in {"all", "path-search"} else "result.json"
summary = json.loads((Path(out_dir) / primary).read_text())

if summary["status"] == "error":
    chain = summary.get("error_class_chain", [])
    if "OptimizationError" in chain:
        # retry with looser convergence threshold
        ...
    else:
        raise RuntimeError(summary["error"])
```

`all` / `path-search` は、集約結果の書き込み処理まで到達すると正規の
`summary.json` を書きます。段階別・レポート系コマンドは、`--out-json` を指定して
正常終了した場合に `result.json` と互換ミラーを書きます。捕捉した実行時例外では、
フラグなしでも可能な範囲でエラーエンベロープを書くことがあります。使用法の検証で
終了した場合や出力ディレクトリの確定前は、JSON が存在するとは限りません。
