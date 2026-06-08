# 出力ディレクトリのレイアウト

このページでは、各 `mlmm` サブコマンドが出力ディレクトリに書き込むファイルと、エージェントや後段スクリプトが従うべき規約を説明します。

## ファイル名の規約

| ファイル名 | 書き込み元 | 用途 |
|---|---|---|
| `summary.json` | `all`、`path-search`、および `write_result_json` を実行するすべてのステージ別サブコマンド | 正規の JSON エンベロープ（[JSON 出力リファレンス](json-output.md) を参照）。まずこれを読んでください。`write_result_json` を呼ばない純粋なユーティリティ系サブコマンド（例: `fix-altloc`、`add-elem-info`、`bond-summary`）は出力しません。 |
| `result.json` | `write_result_json` を呼ぶステージ別サブコマンド（`opt`、`tsopt`、`freq`、`irc`、`sp`、`scan` / `scan2d` / `scan3d`、`path-opt`、`dft`、`extract`） | 別名ファイル — ペイロードは `summary.json` と同一です。単一ファイル名の規約に従い、`summary.json` を読んでください。`result.json` は同じ内容を持つため、`summary.json` のみを利用する場合は削除できます。 |
| `summary.log` | `path-search`、`all` | 人が読むための実行ログ（セグメント / ステージごとに 1 行）。 |
| `final_geometry.xyz` | `opt`、`tsopt` | 最適化された構造（XYZ、フル精度）。 |
| `mep.pdb` / `mep_trj.xyz` | `path-search`、`all` | 反応経路のフレーム（PDB / XYZ）。単独実行の `path-opt` は代わりに `final_geometries_trj.xyz` / `final_geometries.pdb` を書き込みます。 |
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
| `extract` | `./`（作業ディレクトリに `<input>_<cluster>.pdb` を書き込み） |
| `mm-parm` | `./`（`<base>.parm7` / `<base>.rst7` を書き込み） |
| `define-layer` | `./`（ラベル付き PDB をインラインで書き込み） |

`--out-dir <path>`（または `-o`）で上書きできます。明示的に指定したパスは、ステージ別デフォルトと YAML の両方に優先します。

## 単独実行と `all`

サブコマンドを単独で実行すると、**フラットな**結果ディレクトリが書き込まれます。同じライターでも `all` によってオーケストレーションされると、構造化されたツリーにネストされます。

- **単独サブコマンド** → 上記のファイルを含むフラットな `result_<subcmd>/`。`segments/` も `_work/` もありません。これらは `all` が 1 回の実行で複数のライターを協調させるときのみ現れます。
- **`all` の内部では、リーフライターはそのままネストされます。** `segments/seg_NN/<subcmd>/` のセグメント別リーフ出力は、単独の `result_<subcmd>/` と構造的に同一です。`all` はライターの出力先を別のディレクトリに向けているだけです。
- **`path-search` / `path-opt` はエンジン側の例外です。** 単独実行では `path-search` 自体が成果物です（`result_path_search/` に独自の `summary.log`、`mep.pdb`、`mep_trj.xyz`、`mep_plot.png`、`energy_diagram_MEP.png` を持ちます）。`all` の内部では、その生の出力は `_work/path_search/` 以下のエンジン用スクラッチであり、マージされた成果物（`mep.pdb`、`mep_trj.xyz`、`mep_plot.png`、`energy_diagram_MEP.png`）はパイプラインのルートに移動され、`summary.{json,log}` がそこにコピーされます。この非対称性は意図的なものです。

したがって `all` のツリーには 3 つのゾーンがあります。

```text
result_all/
├─ summary.log · summary.json                 # copied to the root
├─ mep.pdb · mep_trj.xyz · mep_plot.png · energy_diagram_MEP.png   # MEP products moved from the engine
├─ energy_diagram_*_all.png · irc_plot_all.png
├─ ml_region.pdb                              # ML-region definition (reusable as --model-pdb)
├─ mm_parm/                                   # MM topology <input>.parm7 / .rst7 (reusable as --parm)
├─ layered/                                   # layered full-system PDBs (B-factor annotated; reusable inputs)
├─ segments/
│  └─ seg_NN/                                  # 2-digit per-reactive-segment deliverables
│     ├─ reactant.pdb · ts.pdb · product.pdb         # canonical R/TS/P
│     └─ ts/ · irc/ · freq/ · dft/ · structures/    # per-stage working files (--tsopt / --thermo / --dft)
└─ _work/                                      # pipeline scratch (safe to remove)
   ├─ pockets/ · scan/
   └─ path_search/                             # raw MEP-engine output (path_opt/ with --no-refine-path)
```

TSOPT のみのモードでは MEP ステージがないため、`_work/path_search/` は存在せず、成果物は `segments/seg_01/` 以下に置かれます。モードごとの完全な内訳は [all](all.md) を参照してください。

## エージェント向けレシピ

```python
# Read whichever subcommand's output, single filename across the board.
import json
from pathlib import Path

summary = json.loads((Path(out_dir) / "summary.json").read_text())

if summary["status"] == "error":
    chain = summary.get("error_class_chain", [])
    if "OptimizationError" in chain:
        # retry with looser convergence threshold
        ...
    else:
        raise RuntimeError(summary["error"])
```

`summary.json` は、`write_result_json` を呼んだすべてのサブコマンドで存在することが保証されています（失敗パスを含みます。失敗エンベロープもスキーマバージョンとエラークラスチェーンを保持します）。
