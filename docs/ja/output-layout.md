# 出力ディレクトリのレイアウト

各コマンドが出力ディレクトリに書く主なファイルの名前と、デフォルトの出力ディレクトリ、`all` の中での置き場所を引くページです。各コマンドのファイルの全体は、そのコマンドのページの「主な出力ファイル」にあります。

## ファイル名の規約

| ファイル名 | 書き出し元 | 用途 |
|---|---|---|
| `summary.json` | `all` / `path-search` | 集約ワークフローの JSON 結果（[JSON 出力の一覧](json-output.md)）。 |
| `summary.json` | `--out-json`（デフォルト: `--no-out-json`）を指定した個別計算・レポートのコマンド | そのコマンドの `result.json` の写し。書き込みが正常に終われば同一内容です。`mm-parm`、`define-layer`、`fix-altloc`、`add-elem-info`、`bond-summary`、`oniom-export`、`oniom-import` は書き出しません。 |
| `result.json` | `--out-json` を指定した `opt`、`tsopt`、`freq`、`irc`、`sp`、`scan` / `scan2d` / `scan3d`、`path-opt`、`dft`、`extract`、`trj2fig`、`energy-diagram` | 個別計算・レポートの JSON 結果。収束せずに終わった場合も書き出します。`extract`、`trj2fig`、`energy-diagram` は、最初の出力ファイルと同じ場所に書きます。 |
| `run.log` | 出力ディレクトリが作られた後の CLI / Colab 実行 | シェルでそのまま使える形の実行コマンドと、コマンド実行中の標準出力・標準エラー。実行がどう終わったかを示す行も入り、その行は各コマンドのページにあります（例：[opt の収束の判定](opt.md#収束の判定)）。ヘルプ、バージョン、dry-run の呼び出しと、出力ディレクトリを持たないコマンド（`extract`、`mm-parm`、`define-layer`、`add-elem-info`、`fix-altloc`、`bond-summary`、`trj2fig`、`energy-diagram`、`oniom-export`、`oniom-import`）では生成しません。 |
| `summary.log` | `path-search`、`all` | テキストの要約。ヘッダーに `Scientific status` があり、番号付きの節にセグメントごとの障壁と結合の変化、出力のツリーが入ります（[all の実行結果の判定](all.md#実行結果の判定)）。 |
| `final_geometry.xyz` / `final_geometry.pdb` | `opt`、`tsopt` | 最適化された構造で、次のコマンドに渡す構造です。最適化が収束しなかった場合も書きます。`.xyz` は常に書き、`.pdb` は PDB/mmCIF 入力のときに書き、B-factor に層を持ちます。 |
| `mep_trj.pdb` / `mep_trj.cif` / `mep_trj.xyz` | `path-search` | 反応経路のフレーム。mmCIF・大きな PDB の入力でファイル変換が有効なときは `.cif` も書き出します。 |
| `final_geometries_trj.xyz` / `hei.xyz` | `path-opt` | 反応経路のフレーム（全イメージ）と最高エネルギーイメージ（HEI）。ファイル変換が有効なときは `.pdb` の写し（`final_geometries.pdb`、`hei.pdb`）と、mmCIF・大きな PDB の入力では `.cif` の写しも書き出します。 |
| `mep_plot.png` / `energy_diagram_MEP.png` | `path-search` | 最小エネルギー経路（MEP）のエネルギー。`mep_plot.png` は経路に沿った図、`energy_diagram_MEP.png` は状態エネルギー図です。`all` はこの 2 つをルートに置きます。 |
| `finished_irc_trj.xyz` / `forward_irc_trj.xyz` / `backward_irc_trj.xyz` | `irc` | IRC（固有反応座標）の軌跡（経路全体と各方向）。参照トポロジーがあれば `.pdb` の写し、mmCIF・大きな PDB の入力では `.cif` の写しも書き出します。 |
| `forward_first.xyz` / `backward_last.xyz` | `irc` | 2 つの分岐の端で、[`opt`](opt.md) で最適化する端点の候補です。 |
| `frequencies_cm-1.txt` | `freq` | 振動モードの一覧。 |
| `ml_region_without_linkH.xyz` / `ml_region_with_linkH.xyz` | `dft` | DFT 計算が見る ML 領域の、リンク水素を付ける前と後の構造。PDB 入力では `.pdb` の写しも書きます。`all` の中では各 `segments/seg_NN/dft/{R,TS,P}/` にあります。 |
| `*.pdb` / `*.cif` | `--convert-files`（デフォルト。`--no-convert-files` で無効）を持つコマンドと `extract` | 入力の形式に応じた出力の写しで、元の出力の横に書き出します。PDB/mmCIF 入力では PDB、mmCIF・大きな PDB の入力ではさらに元の chain ID・残基番号・挿入コードを保つ `.cif` です。複数フレームの軌跡は複数モデルの CIF になります。`extract` には変換の切り替えがなく、mmCIF・大きな PDB の入力では常に `.cif` を書き出します。 |

## デフォルトの `--out-dir`

| サブコマンド | デフォルトの `--out-dir` |
|---|---|
| `all` | `./result_all/` |
| `opt` | `./result_opt/` |
| `tsopt` | `./result_tsopt/` |
| `freq` | `./result_freq/` |
| `irc` | `./result_irc/` |
| `dft` | `./result_dft/` |
| `scan` | `./result_scan/` |
| `scan2d` | `./result_scan2d/` |
| `scan3d` | `./result_scan3d/` |
| `path-opt` | `./result_path_opt/` |
| `path-search` | `./result_path_search/` |
| `sp` | `./result_sp/` |
| `extract` | `./`（`pocket.pdb` を書き出し。入力が複数の場合は `pocket_<input>.pdb`） |
| `mm-parm` | `./`（`<prefix>.parm7` と `<prefix>.rst7` を書き出し。接頭辞のデフォルトは入力ファイル名から `.pdb` を除いたもの） |
| `define-layer` | `./`（`<input>_layered.pdb` を書き出し） |

ほかのディレクトリに書くには `-o/--out-dir <path>` を指定します。準備のコマンドはファイルパスを取ります：`extract` は `-o/--output <file>` に 1 つ以上のパス、`mm-parm` は `-o/--out-prefix <prefix>`、`define-layer` は `-o/--output <file>` です。

## スタンドアロン と `all` の比較

単独実行では `result_<subcmd>/` にファイルが並び、`segments/` や `_work/` はありません。`all` では、後処理の各段階が同じファイル構成で `segments/seg_NN/` の下の `ts/`・`irc/`・`freq/`・`dft/` に配置されます。

- **`path-search` / `path-opt` は配置が異なります。** 単独で実行すると、上の表のファイルを `result_path_search/` か `result_path_opt/` に書きます。`all` の中の MEP 探索は、デフォルトでは `path-opt`、`--refine-path` を付けると再帰的な `path-search` で行います。生の出力は `_work/path_opt/` か `_work/path_search/` に残り、MEP のファイル（`mep_trj.*`、`mep_plot.png`、`energy_diagram_MEP.png`）だけがルートに移されます。

したがって `all` のツリーには 3 つのゾーンがあります。下のツリーは主な項目だけで、`seg_NN/` の中を含むすべてのファイルとエネルギー図の名前は [all の主な出力ファイル](all.md#主な出力ファイル) にあります。

```text
result_all/
├─ summary.log · summary.json                        # 実行の要約
├─ mep_trj.pdb · mep_trj.cif · mep_trj.xyz           # MEP 座標
├─ mep_plot.png · energy_diagram_MEP.png · energy_diagram_*_all.png · irc_plot_all.png
├─ ml_region.pdb                                     # ML 領域（--model-pdb で再利用できる）
├─ mm_parm/                                          # Amber トポロジー <input>.parm7 / .rst7（--parm7 で再利用できる）
├─ layered/                                          # B-factor に層を書いた全系の構造
├─ segments/
│  └─ seg_NN/                                        # 反応の 1 段: seg_01, seg_02, ...
│     ├─ reactant.{pdb,cif} · ts.* · product.*       # 最適化した R・TS・P（--tsopt）
│     └─ ts/ · irc/ · endpoint_opt/ · freq/{R,TS,P}/ · dft/{R,TS,P}/  # ステージ別の作業ファイル（--tsopt / --thermo / --dft）
└─ _work/                                            # 途中のファイル（TS 候補の HEI を含む）
   ├─ pockets/ · scan/
   └─ path_opt/                                      # MEP 探索と hei_seg_NN.*（--refine-path のときは path_search/）
```

TS-only モード（遷移状態（TS）の候補 1 つに `--tsopt` を付け、`-s/--scan-lists` を付けずに与える実行）では MEP ステージがないため、`_work/path_opt/` は存在せず、成果物は `segments/seg_01/` 下に置かれます。

## 使用上の注意点

* **早い段階で止まった実行**: 引数や入力の検査で止まった実行は、出力ディレクトリができる前に終わるため、上のファイルを 1 つも書かないことがあります。
* **`--out-json` なしの `summary.json` / `result.json`**: 正常終了した個別計算のコマンドがこれらを書くのは、`--out-json` 指定時だけです。出力ディレクトリを用意した後に例外で実行が止まった場合は、フラグが無くても `"execution_status": "failed"` と `"error_type"` を含む両方を書き出します。

## 関連ドキュメント

- [all](all.md#主な出力ファイル) — `result_all/` のツリー全体とエネルギー図
- [JSON 出力の一覧](json-output.md) — `summary.json` と `result.json` のキー、Python と jq の例
- [共通オプションと残基・原子の指定](cli-conventions.md) — `--out-dir`、`--convert-files`、終了コード
- [トラブルシューティング](troubleshooting.md) — よくあるエラーと対処法
