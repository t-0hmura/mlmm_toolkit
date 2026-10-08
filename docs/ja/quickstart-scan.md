# クイックスタート: `mlmm all --scan-lists`

`mlmm all --scan-lists`（`-s`）は、1 つの構造から反応経路を作ります。指定した距離を ML/MM の全系の上で調和拘束のもとで目標値まで動かし（スキャン）、スキャンの端点から最小エネルギー経路（MEP）を探索します。`--tsopt` を付けると、遷移状態（TS）の最適化と固有反応座標（IRC）の計算まで進みます。[クイックスタート: `mlmm all`](quickstart-all.md) と同じく、同じ実行の最初に ML 領域を切り出し、Amber トポロジーと層を組みます。反応座標を自分で決めるので、ほかの箇所の構造変化が起きにくく、目当ての TS を得やすくなります。出発の構造は、反応物でも中間体や生成物でもかまいません。

以下のコマンドは同梱の [`examples/beza/`](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples/beza) にある `1.R.pdb` を使うので、そのディレクトリで実行します。

---

## 主な用途

* **生成物の構造が無い**: できる結合と切れる結合を動かして、反応物（R）から生成物（P）を作る
* **2 段の反応**: メチル基転移の後のプロトン移動のように、1 段ずつ順に動かす
* **スキャンから MEP と TS へ**: 同じ実行のまま、スキャンした経路から MEP と TS へ進む

## スキャンコマンドの選び方

| 目的 | コマンド | トポロジーと ML 領域 |
| --- | --- | --- |
| 拘束した構造とスキャンの軌跡だけを得る | `mlmm scan` | `--parm7` と ML 領域の指定（`--model-pdb`、`--model-indices`、入力の B-factor の層のどれか）を渡す |
| スキャンから MEP へ、必要なら TS と IRC まで進む | `mlmm all -s ...` | `all` が自動で作る |
| 2 つか 3 つの座標でエネルギーの 2D・3D マップを作る | `mlmm scan2d` / `scan3d` | `--parm7` と ML 領域の指定を渡す |

## 基本的な実行例

### 1. まず入力を確かめる

`--dry-run` は、入力・層・ML 領域の切り出し・`-s` の書き方・電荷とスピン・トポロジーの原子数と並び・AmberTools のコマンドを確かめ、計算をせずに止まります。

```bash
mlmm all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
  -s '[(4360, 4419, 1.50)]' --dry-run
```

端末に `[all] Dry-run validation passed: ...` と `[all] Planned stages: extract -> mm_parm -> scan -> path_opt.` が出て、最後が `[Dry run] --dry-run completed. Input command is valid.` なら確認は通っています。通らないときはエラーメッセージが出て止まるので、そのメッセージを [トラブルシューティング](troubleshooting.md) で探してください。

### 2. 1 段のスキャン

SAM の CS1（原子 4360）と GPP の C7（原子 4419。IUPAC の番号では GPP の C6）の距離を 1.50 Å まで動かします。原子は番号でも名前でも指定でき、次の 2 つは同じスキャンです。

```bash
# 原子の番号（既定は 1 始まり）
mlmm all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 -s '[(4360, 4419, 1.50)]' -o ./result_scan

# 原子の名前
mlmm all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 -s '[("SAM,320,CS1", "GPP,321,C7", 1.50)]' -o ./result_scan
```

端末の最後のほうの `====== Pipeline summary ======` の下に `Scientific status: success` と出れば成功で、`summary.json` の `scientific_status` にも同じ値が入ります。

### 3. 2 段のスキャン

`-s` の後に並べる角括弧のリスト 1 つが 1 段で、段は順に実行されます。このリストは文字で書き下したものなので、リテラルと呼びます。段 1 では SAM のメチル炭素（CS1）を GPP の C7 に近づけ（1.50 Å）、SAM の SD から離します（3.30 Å）。段 2 では GPP の H11 を C7 から離し（2.90 Å）、Glu186 の OE2 に移します（1.00 Å）。

```bash
mlmm all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -s \
  '[("SAM,320,CS1","GPP,321,C7",1.50),("SAM,320,CS1","SAM,320,SD",3.30)]' \
  '[("GPP,321,C7","GPP,321,H11",2.90),("GLU,186,OE2","GPP,321,H11",1.00)]' \
  -o ./result_scan
```

どちらの段にも、できる結合と切れる結合の両方を入れています。自分の反応の段の分け方は {ref}`反応の分け方を決める <ja-mechanism-split>` を参照してください。

## `--scan-lists` の書き方

* **タプル**: 各リテラルはタプルのリストです。距離は `(atom1, atom2, target_Å)`、角度は `(atom1, atom2, atom3, target_deg)`、二面角は `(atom1, atom2, atom3, atom4, target_deg)` と書きます。距離の単位は Å、角度と二面角の単位は度です。
* **原子**: 全系の入力構造での順番を整数で書くか、原子の名前を二重引用符で囲んで書きます。番号は 1 から数え、`--scan-zero-based` では 0 から数えます。スキャンは全系の上で行うので、`-c` を付けても番号は入力構造の番号のままです。
* **段**: リテラル 1 つが 1 段で、同じリテラルの中のタプルは同時に動きます。リテラルを並べると段が順に実行され、各段は前の段の結果から始まります。`-s` は 1 回だけ書き、その後にリテラルを並べてください。
* **引用符**: 括弧や空白をシェルに解釈させないよう、各リテラルを一重引用符で囲みます。
* **同梱の PDB**: chain 欄が空なので、chain ID を含む書き方ではなく、`"SAM,320,CS1"` のような名前か番号の書き方を使います。

原子の名前と引用符の書き方の詳細は {ref}`スキャンリスト仕様 <ja-scan-list-spec>` を参照してください。

## 主な出力ファイル

上のコマンドは次のファイルを書き出します。

```text
result_scan/
├── summary.log
├── summary.json                 # 結果（scientific_status を含む）
├── mep_trj.pdb                  # 全セグメントの MEP
├── energy_diagram_MEP.png       # MEP のエネルギープロファイル
├── ml_region.pdb                # ML 領域（--model-pdb で再利用できる）
├── mm_parm/                     # Amber トポロジー 1.R.parm7 と 1.R.rst7（--parm7 で再利用できる）
├── layered/                     # 3 つの層を B-factor に書いた 1.R_layered.pdb
└── _work/                       # 途中のファイル（TS 候補の HEI を含む。実行後も残る）
    ├── scan/
    │   ├── preopt/              # 最適化した出発構造
    │   ├── stage_01/            # スキャンの段 1
    │   │   ├── result.{xyz,pdb} # 拘束した端点（--scan-endopt のときだけ拘束なしで最適化）
    │   │   ├── scan_trj.xyz     # スキャンの軌跡
    │   │   └── scan.pdb
    │   ├── stage_02/            # スキャンの段 2（2 段の実行）
    │   └── result.json          # 各段のスキャンの結果
    └── path_opt/                # MEP 探索（MEP を再帰的に詰める --refine-path のときは path_search/）
        └── hei_seg_01.{xyz,pdb} # セグメント 1 の最高エネルギーのイメージ
```

これらのコマンドは MEP 探索で終わるため、`segments/` は作られません。`--tsopt` を付けると、反応セグメントごとに `segments/seg_NN/` に R/TS/P の構造と IRC が入り、`--thermo` を付けると `freq/` も加わります。

## 結果の確認

1. **完了状況**: `scientific_status` には、求めた段がすべて収束すると `success`、そうでなければ `partial` か `failed` が入り、[理由](json-output.md#実行と要求段階の完了状況)は `scientific_status_reasons` に出ます。`--tsopt` のとき、虚振動のモードができる結合と切れる結合を動かすかと、端点が狙った R と P かの 2 つは自分で確かめてください。
2. **スキャン**: `_work/scan/stage_01/scan_trj.xyz` をビューアで開き、距離が狙いどおりに変わるかを確かめます。各段の終わりに端末に `[stage 1] Covalent-bond changes (start vs final): Yes` か `No` が出て、`_work/scan/result.json` の `stages[].bond_changes` に記録されます。
3. **MEP**: `mep_trj.pdb` と、最高エネルギーのイメージ（HEI、TS の候補）`_work/path_opt/hei_seg_01.pdb` を開き、`energy_diagram_MEP.png` にはっきりした障壁があるかを確かめます。
4. **TS（`--tsopt` のとき）**: TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます。このとき端末に `[microiter] Converged!` が出て、続いて `[Imaginary modes] n=1 (...)` が出ます。`segments/seg_01/ts/vib/imag_*_trj.xyz` をビューアで開き、できる結合と切れる結合に沿って原子が動くかを確認してください。
5. **端点（`--tsopt` のとき）**: `segments/seg_01/irc/finished_irc_trj.xyz` と、最適化した端点の `segments/seg_01/reactant.pdb`・`product.pdb` を開き、狙った R と P かを確かめます。IRC が収束しなくても、端点の最適化で狙った R と P に着けば、その結果は使えます。

## 使用上の注意点

* **入力**: 全系の PDB か mmCIF の構造 1 つです。`-c` を付けると、指定した残基のまわりから ML 領域を切り出します。省くと、入力の B-factor の層か `--model-pdb` から ML 領域を取ります。
* **`all` と `scan` の既定値**: 2 つのコマンドはスキャンの処理を共有しますが、オプションの名前と、スキャンの前後の最適化の既定値が違います。

  | コマンド | 刻み幅 / 拘束 | 緩和の上限 | スキャンの前 / 後の最適化 |
  | --- | --- | --- | --- |
  | `mlmm all` | `--scan-max-step-size 0.20` Å、`--scan-restraint-k 300` eV/Å² | `--scan-relax-max-cycles 100000` | 前: on（`--preopt/--no-preopt` に従う。`--scan-preopt/--no-scan-preopt` で別に指定できる）、後: off（`--scan-endopt/--no-scan-endopt`） |
  | `mlmm scan` | `--max-step-size 0.20` Å、`--restraint-k 300` eV/Å² | `--relax-max-cycles 100000` | 前: off（`--preopt/--no-preopt`）、後: off（`--endopt/--no-endopt`） |

  拘束の強さは YAML の [`bias.k`](yaml-reference.md#bias) でも指定できます。
* **結合変化と `--refine-path`**: スキャンの端点は、その段で結合変化が出たかどうかに関係なく、すべて MEP 探索に渡されます。MEP を再帰的に詰める（`path_search/`）かどうかは `--refine-path` だけで決まります。
* **スキャンの端点**: スキャンが終わって得られるのは拘束した構造です。拘束なしの最適化、または TS 最適化と IRC で確かめるまでは、極小点とも遷移状態とも言えません。
* **`scan` を単独で使うとき**: 単独の [`scan`](scan.md) は YAML・JSON のスペックファイルも受け付けます。`all -s` が受け付けるのはインラインのリテラルだけで、スペックファイルを渡すとエラーで止まります。`scan` の既定の `--no-preopt --no-endopt` のままでは、入力構造をそのまま出発点にし、各段の終わりは拘束したままの構造になります。
* **`scan --dry-run`**: 入力、電荷とスピン、`--scan-lists` の解釈を確かめ、予定のスキャンを表示して止まります。`all` の切り出しとトポロジーの作成は `all --dry-run` でしか確かめられません。

## 次のステップ

- [クイックスタート: TS-only モード](quickstart-tsopt.md): スキャンの最高点などの TS 候補を最適化して確かめる
- [反応機構を調べるコツ](mechanism-tips.md): 反応の段の分け方（「反応の分け方を決める」）
- {ref}`原子の固定と距離の拘束 <ja-freeze-atoms-and-restraints>`: 拘束の強さの意味
- [`scan`](scan.md): スキャンを単独で実行する
- [`all`](all.md): オプションの説明。`mlmm all --help-advanced` でも見られます
- [用語集](glossary.md): MEP・TS・IRC などの用語
- [トラブルシューティング](troubleshooting.md): エラーメッセージや症状から対処を探す
