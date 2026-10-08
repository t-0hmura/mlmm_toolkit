# `trj2fig`（軌跡のエネルギープロファイル）

`trj2fig` サブコマンドは、`opt`・`scan`・`path-opt`・`path-search`・`irc` などが書き出した **XYZ 軌跡に沿ったエネルギーを図にします**。各フレームのエネルギーはコメント行から読みます。図と CSV の表を書き出せます。

## 主な用途

* **経路のエネルギープロファイル**: 最小エネルギー経路（MEP）・スキャン・固有反応座標（IRC）の軌跡の ΔE の図
* **最適化の確認**: 最適化のサイクルごとにエネルギーがどう下がったか
* **自分の作図用のデータ**: フレームごとのエネルギーの CSV

---

## 基本的な実行例

### 1. デフォルトの PNG

最初のフレームを基準にした ΔE を描き、`energy.png` に書き出します。

```bash
mlmm trj2fig -i traj.xyz --out-json
```

端末に `[trj2fig] Saved figure -> energy.png` が出て、`result.json` の `n_frames` に読み込んだフレーム数が入れば成功です。

### 2. CSV と SVG、フレーム 5 を基準に Hartree で

0 から数えたフレーム番号 5（6 番目のフレーム）を基準に、表と図を Hartree で書き出します。

```bash
mlmm trj2fig -i traj.xyz -o energy.csv energy.svg -r 5 --unit hartree
```

### 3. 複数の形式、x 軸を反転

最後のフレームを左端に置きます。このとき ΔE の基準も最後のフレームになります。

```bash
mlmm trj2fig -i traj.xyz --reverse-x -o energy.png energy.html energy.pdf
```

### 4. MLIP でエネルギーを計算し直す

コメント行を読む代わりに、UMA のモデルと精度を指定して（[MLIP バックエンド](backends.md)）、中性の一重項として全フレームを計算し直し、使ったバックエンドを `result.json` に記録します。

```bash
mlmm trj2fig -i traj.xyz -q 0 -m 1 -b uma --backend-model uma-s-1p2 \
    --precision fp32 -o energy.png energy.csv --out-json
```

---

## 処理の仕組みと計算仕様

1. **エネルギーの読み込み**:
各フレームのコメント行からエネルギーを読みます。mlmm が書き出す軌跡は、そのまま読めます（例: `optimization_trj.xyz`（[`opt`](opt.md) の `--dump`）、`scan_trj.xyz`（[`scan`](scan.md)）、`mep_trj.xyz`、`finished_irc_trj.xyz`。後の 2 つは [出力ディレクトリのレイアウト](output-layout.md)）。ほかのファイルでは `E=<値>` の形で書き、単位（`Ha`・`Eh`・`hartree`・`eV`・`kcal/mol`）を付けられます。単位が無ければ Hartree ですが、拡張 XYZ（コメントに `Properties=` か `Lattice=` がある）の `energy=` だけは eV とします。`-q` か `-m` を指定したときは、代わりに MLIP（機械学習原子間ポテンシャル）バックエンドで全フレームを計算し直します。
2. **基準の選択**:
`-r init`（デフォルト）は左端のフレームで、最初のフレーム、`--reverse-x` では最後のフレームです。整数は 0 始まりのフレーム番号、`none` は絶対エネルギーです。
3. **単位の変換**:
エネルギーを kcal/mol（デフォルト）か Hartree に変換し、基準の値を引いて ΔE にします。y 軸は `ΔE (kcal/mol)`、絶対エネルギーでは `E (…)` になります。
4. **書き出し**:
出力ごとに拡張子で形式を決めます。`.png`・`.jpg`・`.jpeg`・`.svg`・`.pdf`・`.html` は図、`.csv` は表です。PNG は図の縦横 2 倍の画素数で書き出します。

---

## 主な出力ファイル

```text
energy.png      # 図（出力の指定が無いときのデフォルト）
energy.csv      # エネルギーの表（.csv の出力を指定したとき）
result.json     # 要約（--out-json 指定時）
summary.json    # result.json の写し。result.json を読む（--out-json 指定時）
```

* **CSV の列**: `frame`、`energy_hartree`、図に描いた値（`--unit` の単位）の 3 列です。3 列目の名前は、基準があるときは `delta_kcal` か `delta_hartree`、`-r none` では `energy_kcal` か `energy_hartree` です。
* **`result.json`** は最初の出力と同じディレクトリに書き出します。`n_frames`、`min_energy_hartree`、`max_energy_hartree`、`energy_source`（`trajectory_comment` か `mlip_recomputed`）と、書き出したファイルを順に並べた `output_files` を持ちます。`mlip_backend`・`mlip_model`・`mlip_model_label`・`mlip_task`・`mlip_precision` はエネルギーを計算し直したときだけ値が入り、それ以外では null です（[JSON 出力リファレンス](json-output.md#trj2fig)）。

---

## 主な CLI オプション

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | XYZ 軌跡 |
| `-o, --output` | パス | `energy.png` | 出力ファイル（`.png`, `.jpg`, `.jpeg`, `.html`, `.svg`, `.pdf`, `.csv`）。`-o` を繰り返すか、その後にファイル名を続けて並べる（`-o energy.csv energy.svg`） |
| `--unit` | `kcal` / `hartree` | `kcal` | 図と表の値の単位 |
| `-r, --reference` | 文字列 | `init` | 基準: `init`、`none`、または 0 始まりのフレーム番号 |
| `-q, --charge` | 整数 | `None` | フレーム全体の総電荷。指定すると MLIP でエネルギーを計算し直す |
| `-m, --multiplicity` | 整数 | `None` | スピン多重度（2S+1）。指定すると MLIP でエネルギーを計算し直す |
| `--reverse-x/--no-reverse-x` | フラグ | `False` | 最後のフレームを左端に置く |
| `-b, --backend` | 文字列 | `uma` | 計算し直しに使う MLIP（`uma`, `orb`, `mace`, `aimnet2`。[MLIP バックエンド](backends.md)） |
| `--backend-model` | 文字列 | `None` | 選んだバックエンドのモデル（例: `uma-s-1p2`）。省略時はそのバックエンドのデフォルトのモデル |
| `--precision` | `fp32` / `fp64` | バックエンドごと | 計算し直しの精度（UMA は `fp32`、ORB と MACE は `fp64`）。AIMNet2 は `fp32` だけ |
| `--out-json/--no-out-json` | フラグ | `False` | `result.json` と `summary.json` を出力 |

全オプションの一覧は [自動生成 CLI リファレンス](../reference/commands/trj2fig.md) を参照してください。

---

## 使用上の注意点

* **計算し直し**: `-q` と `-m` の片方だけを指定すると、もう片方は電荷 0 または多重度 1 になります。
* **計算し直したエネルギーは MLIP のエネルギー**: MLIP はフレームの全原子を ML/MM に分けずに計算するので、結果は ML/MM ONIOM のエネルギーではありません。mlmm の軌跡の ONIOM エネルギーを描くときは、`-q`・`-m` を付けずにコメント行から読んでください。
* **エネルギーを読めないコメント行**は、そのフレームを示すエラーで止まります。
* **対応しない拡張子**はエラーで止まります。
* **終了コード**: {ref}`終了コード <ja-exit-codes>`を参照してください。

---

## 関連ドキュメント

* [path-search](path-search.md) — 図にする MEP の軌跡
* [irc](irc.md) — 図にする IRC の軌跡
* [energy-diagram](energy-diagram.md) — 与えた数値から状態エネルギー図を描く
* [all](all.md) — 全工程のワークフロー
* [トラブルシューティング](troubleshooting.md) — 実行に失敗したときの対処。図の書き出しの失敗は {ref}`図のエクスポートが失敗する <ja-plot-export-fails-chrome-missing>`
