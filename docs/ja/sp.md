# `sp`（一点計算）

`sp` サブコマンドは、1 つの構造の **ML/MM ONIOM エネルギーと原子に働く力**を計算し、`--hess` を付けると動ける原子の **Hessian** も計算します。ML 領域は選んだバックエンドで、酵素の残りは `--parm7` の Amber 力場で計算します。構造最適化は行わず、入力の構造のまま評価します。

---

## 主な用途

* **最適化の前の確認**: ML 領域・電荷・多重度が受け付けられ、バックエンドが有限のエネルギーと力を返すかを確かめる
* **バックエンドの比較**: 同じ構造と ML 領域を UMA・ORB・MACE・AIMNet2・DFT（`-b dft`）で評価する
* **参照値の作成**: 力と Hessian を `.npy` ファイルとして、エネルギーを端末か `result.json` から得て、自分の解析に使う

---

## 基本的な実行例

### 1. エネルギーと力

デフォルトのバックエンド（UMA）で、中性の一重項の ML 領域を評価します。`enzyme.pdb` は全系、`real.parm7` はその Amber トポロジー、`ml_region.pdb` は ML 領域の原子です。`-q` と `-m` は ML 領域の電荷と多重度です。

```bash
mlmm sp -i enzyme.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 --out-json
```

端末に `[sp] energy = … a.u.  |force|_max = … a.u./bohr` が出て、`result_sp/` に `forces.npy` と、`energy_au` を持つ `result.json` があれば成功です。

### 2. Hessian も計算する

`--hess` を付けると、動ける原子の Hessian も計算します。

```bash
mlmm sp -i enzyme.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 --hess
```

---

## 処理の仕組みと計算仕様

1. **ML/MM の系を組む**:
`-i` から全系を、`--parm7` から Amber トポロジーを、`--model-pdb`・`--model-indices`・入力の B-factor のどれかから ML 領域を読みます。電荷は `-q`、または PDB/mmCIF 入力での `-l` から決まります。固定 MM 層と `--freeze-atoms` で指定した原子を固定します。
2. **エネルギーと力**:
入力の構造で、ML 領域をバックエンドで、MM 原子を力場で 1 回ずつ計算し、ONIOM のエネルギーと力に組み合わせます。エネルギーと力の最大成分を端末に表示し、力を `forces.npy` に保存します。固定原子に働く力は 0 です。
3. **Hessian（`--hess` 指定時）**:
Hessian に入るのは ML 領域と可動 MM 原子で、固定原子は入りません。`--hessian-calc-mode FiniteDifference`（デフォルト）は力を数値微分し、`Analytical` は ML 領域に UMA・ORB・MACE・AIMNet2 の解析 Hessian を使います。`Analytical` は `--uma-workers`（MLIP の並列ワーカー数）を 2 以上にすると使えません。MM の部分はデフォルトで有限差分で、YAML の `calc.mm_fd: false` で `hessian_ff` の解析 Hessian になります。

---

## 主な出力ファイル

`--out-dir` に以下のファイルを書き出します。

| ファイル | 内容 | 書き出す条件 |
| --- | --- | --- |
| `forces.npy` | 全系の全原子についての ONIOM の力の `(N, 3)` 配列（Hartree/bohr） | 常に |
| `hessian.npy` | 質量重み付けなしの ONIOM Hessian（Hartree/bohr²）。Hessian に入る M 原子（入力の順）の `(3M, 3M)` | `--hess` 指定時 |
| `result.json` | エネルギー（`energy_au`）、バックエンド、モデル、電荷、多重度、ML 領域（指定元と原子数）、`.npy` ファイルのパス、経過時間 | `--out-json` 指定時 |
| `summary.json` | `result.json` と同じ内容 | `--out-json` 指定時 |

---

## 主な CLI オプション

ML/MM の計算コマンドに共通のオプションは {ref}`ML/MM の共通オプション <ja-mlmm-options>` に 1 か所でまとめてあります。下の表は `sp` に固有のものだけです。

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `-i, --input` | パス | （必須） | 全系の入力構造ファイル（`.pdb`, `.cif`, `.mmcif`、または `--ref-pdb` と組み合わせた `.xyz`） |
| `-q, --charge` | 整数 | `None` | ML 領域の電荷。`-l` を使う場合のほかは必須 |
| `-m, --multiplicity` | 整数 | `1` | ML 領域のスピン多重度（2S+1） |
| `-l, --ligand-charge` | 文字列 | `None` | 残基ごとの形式電荷（例: `'SAM:1,GPP:-3'`）またはリガンドの総電荷。`-q` を省いたときに ML 領域の電荷を求めるのに使用（PDB/mmCIF 入力または `--ref-pdb`） |
| `-b, --backend` | 文字列 | `uma` | ML 領域のバックエンド（`uma`, `orb`, `mace`, `aimnet2`, `dft`）。`-b dft` の設定は [MLIP の TS を DFT で確かめる](dft-backend.md) を参照 |
| `--hess/--no-hess` | フラグ | `False` | Hessian も計算して `hessian.npy` に書き出す |
| `--hessian-calc-mode` | `FiniteDifference` / `Analytical` | `FiniteDifference` | Hessian の計算法（有限差分 / 解析的）。`--hess` と併用 |
| `--hessian-cutoff` | 浮動小数点数 | `None` | ML 領域からこの距離（Å）以内の可動 MM 原子だけを Hessian に入れる。デフォルトでは可動 MM 原子すべて |
| `--freeze-atoms` | 文字列 | `None` | 固定する原子インデックス（1 始まり、カンマ区切り: 例 `'1,3,5'`） |
| `--embedcharge/--no-embedcharge` | フラグ | `False` | MM の点電荷による静電埋め込み（MLIP では xTB の補正、`-b dft` では PySCF の点電荷） |
| `-o, --out-dir` | パス | `./result_sp/` | 出力先ディレクトリ |
| `--out-json/--no-out-json` | フラグ | `False` | `result.json` と `summary.json` を出力 |

全オプションは [自動生成のオプションの一覧（英語のみ）](../reference/commands/sp.md) を参照してください。

> **補足:** YAML（`--config`）では、`calc` でバックエンドを設定し、`geom.freeze_atoms`（1 始まり）で固定原子を追加できます。`geom.freeze_atoms` は `--freeze-atoms` と合わせて使われます。

---

## 使用上の注意点

* **失敗したとき**: `ML region electron count inconsistent` のような 1 行の `Error: …` か、トレースバック付きの `Unhandled error during single-point:` が出て、0 以外の終了コードで終わります。
* **エネルギーがおかしいとき**: 有限の値でもおかしいときは、ML 領域とその {ref}`電荷・多重度 <ja-charge--spin>` を見直してください。
* **固定原子**: インデックスは 1 始まりで、固定原子に働く力は 0 になります。固定 MM 層も固定されます。
* **原子電荷**: `sp -b dft` が出すのは ML(DFT)/MM のエネルギーと力だけです。ML 領域の Mulliken・meta-Löwdin・IAO の電荷が必要なときは [`dft`](dft.md) を使ってください。
* **終了コード**: {ref}`終了コード <ja-exit-codes>`を参照してください。

---

## 関連ドキュメント

* [opt](opt.md) — 構造最適化
* [tsopt](tsopt.md) — 遷移状態（TS）候補の構造最適化
* [freq](freq.md) — 振動解析と熱化学
* [dft](dft.md) — 原子電荷も出す ML 領域の DFT 一点計算
* [MLIP の TS を DFT で確かめる](dft-backend.md) — `-b dft` の設定（`--func-basis`・`--dft-engine`）と GPU メモリ
* [MLIP バックエンド](backends.md) — バックエンド・精度・ワーカーの選び方
* [トラブルシューティング](troubleshooting.md) — 実行に失敗したときの対処
