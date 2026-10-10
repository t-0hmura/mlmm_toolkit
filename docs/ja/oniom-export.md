# `oniom-export`（Gaussian ONIOM・ORCA QM/MM の入力の書き出し）

`oniom-export` サブコマンドは、mlmm の ML/MM 系を **Gaussian ONIOM（`--mode g16`）または ORCA QM/MM（`--mode orca`）の入力ファイルとして書き出します**。Amber トポロジー（`--parm7`）と、層を B-factor に持つ PDB を読みます。ML 領域を QM 領域として、座標・QM 原子・可動原子・MM パラメータを 1 つの入力ファイルにまとめます。トポロジーは CMAP 項を含まないものが必要で、作り方は [mm-parm](mm-parm.md#oniom-export-用の-cmap-を含まないトポロジー) にあります。

---

## 主な用途

* **Gaussian ONIOM**: mlmm で得た構造（TS 候補など）を、高層を DFT にした Gaussian の ONIOM 計算にかける
* **ORCA QM/MM**: 同じ系を ORCA の QM/MM 機能で計算する
* **往復**: 書き出した入力を mlmm の外で直し、[`oniom-import`](oniom-import.md) の `--ref-pdb` で原子名・残基名つきで戻す

---

## 基本的な実行例

### 1. Gaussian ONIOM（--mode g16）

`mlmm tsopt` で得た TS 候補を Gaussian ONIOM の入力に書き出します。`result_tsopt/final_geometry.pdb` は層を B-factor に持つ全系の構造、`real.parm7` は同じ系のトポロジー、`ml_region.pdb` は QM 原子を選ぶ PDB です。

```bash
mlmm oniom-export --mode g16 --parm7 real.parm7 -i result_tsopt/final_geometry.pdb \
    --model-pdb ml_region.pdb -o ts_refine.com -q 0 -m 1
```

端末に `[oniom-gaussian] Wrote 'ts_refine.com'` と、QM 原子・可動原子・リンク境界の数が出ます。

原子の並びがすでに確かなときは、`--no-element-check` でトポロジーとの元素の 1 原子ずつの照合を省けます。

```bash
mlmm oniom-export --mode g16 --parm7 real.parm7 -i result_tsopt/final_geometry.pdb \
    --model-pdb ml_region.pdb -o ts_refine.gjf -q 0 -m 1 --no-element-check
```

### 2. ORCA QM/MM（--mode orca）

同じ構造を ORCA QM/MM の入力に書き出します。拡張子 `.inp` で ORCA モードに決まるので、`--mode orca` は省けます。

```bash
mlmm oniom-export --parm7 real.parm7 -i result_tsopt/final_geometry.pdb \
    --model-pdb ml_region.pdb -o ts_refine.inp -q 0 -m 1
```

端末に `[oniom-orca] Wrote 'ts_refine.inp'` が出て、ORCA 用の力場ファイルがそろうと `[oniom-orca] ORCAFF.prms: <path>` が続きます。

QM+MM 系全体の電荷と多重度（`Charge_Total`、`Mult_Total`）を自分で指定します。

```bash
mlmm oniom-export --mode orca --parm7 real.parm7 -i result_tsopt/final_geometry.pdb \
    --model-pdb ml_region.pdb -o ts_refine.inp -q 0 -m 1 --total-charge -1 --total-mult 1
```

手元にある `ORCAFF.prms` を使い、変換の手順を省きます。

```bash
mlmm oniom-export --mode orca --parm7 real.parm7 -i result_tsopt/final_geometry.pdb \
    --model-pdb ml_region.pdb -o ts_refine.inp -q 0 -m 1 \
    --orcaff ./ORCAFF.prms --no-convert-orcaff
```

### 3. 手法・プロセッサ数・メモリ

QM の手法と、Gaussian の入力に書くリソース（`%nprocshared`、`%mem`）を変えます。

```bash
mlmm oniom-export --mode g16 --parm7 real.parm7 -i result_tsopt/final_geometry.pdb \
    --model-pdb ml_region.pdb -o ts_refine.com -q 0 -m 1 \
    --method 'wb97xd/def2-svp' --nproc 16 --mem 32GB
```

---

## 処理の仕組みと計算仕様

1. **トポロジーと層**:
parm7 から原子・結合・電荷・Amber パラメータを、`-i` の PDB から座標と B-factor を読みます。PDB は parm7 と同じ原子を同じ順に並べたものが必要で、`--element-check` で元素を 1 原子ずつ照合します。B-factor の 0・10・20（±1.0 以内）が ML・可動 MM・固定 MM を表します。
2. **QM 領域**:
`--model-pdb` を指定すると、その原子が QM 領域です。原子名・残基名・chain・残基番号・挿入コードで `-i` の原子に対応づけます。指定しなければ、B-factor が 0 の原子が QM 領域です。固定 MM 以外のすべての原子が可動で、QM 原子は常に可動です。
3. **QM/MM 境界**:
Gaussian では、切った QM–MM 結合ごとに MM 原子をリンク水素に置き換えます。`--link-atom-method scaled`（デフォルト）は Morokuma/Dapprich の g-factor で位置を決め、`fixed` は QM 原子から 1.09 Å（QM 原子が炭素）または 1.01 Å（窒素）の位置に置きます。ORCA は `QMAtoms` と `ORCAFF.prms` からキャップを自分で作り、入力には推定したキャップの位置をコメントとしてだけ書きます。
4. **入力の書き出し**:
Gaussian の入力には、route の `#p oniom(<method>:amber=softonly)`、可動フラグ（`0` は可動、`-1` は固定）と層（`H` か `L`）つきの座標、結合の情報、Amber パラメータを書きます。電荷と多重度の行は 3 組で、系全体（トポロジーの総電荷と `-m`）、続いて QM 領域を 2 回（`-q` と `-m`）です。ORCA の入力には、`! <method>` と `! QMMM`、`ORCAFFFilename`・`QMAtoms`・`ActiveAtoms`・`Charge_Total`・`Mult_Total` を含む `%qmmm` ブロック、QM 領域の電荷と多重度（`-q`、`-m`）つきの `* xyz` の座標を書きます。ORCA モードでは `ORCAFF.prms` を探すか作ります。`%qmmm` のキーワードは [ORCA 6.0 マニュアル（QM/MM）](https://www.faccts.de/docs/orca/6.0/manual/contents/typical/qmmm.html) にあります。

---

## 主な出力ファイル

* **Gaussian の入力**（`--mode g16`）: `-o` のファイル（`.com` か `.gjf`）。端末に `[oniom-gaussian] Wrote '<file>'`、`QM atoms: N, Movable atoms: M`、`Link boundaries: K` が出ます。
* **ORCA の入力**（`--mode orca`）: `-o` のファイル（`.inp`）。端末に `[oniom-orca] Wrote '<file>'`、`QM atoms: N, Active atoms: M`、`Link boundaries (auto-capped by ORCA): K` が出ます。
* **`ORCAFF.prms`**（ORCA）: `--orcaff` のファイル、または `-o` と同じディレクトリの `<parm7 stem>.ORCAFF.prms`。あれば使い、無ければ `--convert-orcaff` が有効で `orca_mm` が `PATH` にあるときに `orca_mm -convff -AMBER <parm7>` で作ります。それでも無いときは端末に `[oniom-orca] NOTE: ORCAFF.prms not found at '<path>'. Run manually: cd <dir> && orca_mm -convff -AMBER <parm7>` が出ます。このコマンドを実行するまで `.inp` は完成しません。

---

## 主な CLI オプション

| オプション | 引数の型 | デフォルト | 説明 |
| --- | --- | --- | --- |
| `--parm7` | パス | （必須） | 全系の Amber トポロジー。CMAP 項を含まないもの |
| `-i, --input` | パス | （必須） | parm7 と同じ原子順の全系の PDB。B-factor に層（0、10、20）を持つもの |
| `--model-pdb` | パス | `None` | QM 原子の PDB。省略時は `-i` の B-factor が 0 の原子 |
| `-o, --output` | パス | （必須） | 書き出す入力ファイル。`.com` / `.gjf`（g16）または `.inp`（ORCA） |
| `--mode` | `g16` / `orca` | `-o` の拡張子から | 書き出す先のプログラム |
| `--method` | 文字列 | `wB97XD/def2-TZVPD`（g16）、`B3LYP D3BJ def2-SVP`（ORCA） | QM の手法と基底関数。`oniom(<method>:amber=softonly)`（g16）または `!` の行（ORCA）に書く |
| `-q, --charge` | 整数 | （必須） | QM 領域の電荷 |
| `-m, --multiplicity` | 整数 | `1` | QM 領域のスピン多重度 |
| `--nproc` | 整数 | `8` | プロセッサ数（g16 は `%nprocshared`、ORCA は `%pal nprocs`） |
| `--mem` | 文字列 | `16GB` | g16: メモリ（`%mem`） |
| `--total-charge`, `--total-mult` | 整数 | トポロジーの総電荷、`-m` | ORCA: QM+MM 系全体の電荷と多重度（`Charge_Total`、`Mult_Total`） |
| `--orcaff` | パス | `-o` と同じディレクトリの `<parm7 stem>.ORCAFF.prms` | ORCA: 使う `ORCAFF.prms`（既存のファイル） |
| `--convert-orcaff/--no-convert-orcaff` | フラグ | `True` | ORCA: `--orcaff` を指定せず、デフォルトのファイルが無いときに `orca_mm -convff -AMBER` で作る |
| `--element-check/--no-element-check` | フラグ | `True` | `-i` の元素をトポロジーと 1 原子ずつ照合する |
| `--link-atom-method` | `scaled` / `fixed` | `scaled` | g16: リンク水素を g-factor（`scaled`）または固定の結合長（`fixed`）で置く |

全オプションは [自動生成のオプションの一覧（英語のみ）](../reference/commands/oniom_export.md) を参照してください。

---

## 使用上の注意点

* **CMAP**: Gaussian ONIOM は parm7 の CMAP 項を表せず、ORCA の MM エンジンは CMAP 項を使わないため、CMAP を含むトポロジーではファイルを書く前に止まります。書き出しには CMAP を含まないトポロジーを作ってください。mlmm の ML/MM 計算は両方の MM 層に CMAP を適用するので、ふだんの計算では CMAP を残せます。
* **モードの決まり方**: `--mode` が拡張子より優先です。`--mode` が無いと、`.gjf`・`.com`・`.inp` 以外の拡張子はエラーです。
* **原子の並び**: `-i` は parm7 と原子数が同じ PDB（`.pdb` か `.ent`）が必要です。原子数が違うと `--no-element-check` でも止まり、照合が有効なら最初に元素が食い違った原子で `Element sequence mismatch at atom index …`（0 から数えた番号）を出して止まります。
* **系全体の電荷**: Gaussian の real 系の電荷と、デフォルトの ORCA `Charge_Total` は、parm7 の部分電荷の和を整数に丸めた値です。和が整数から 0.05 を超えて離れていると止まります。ORCA モードでは `--total-charge` で指定してください。
* **Gaussian の境界**: 切る QM–MM 結合ごとに別の MM 原子が必要です。2 つの QM 原子が同じ MM 原子に結合していると、Gaussian の書き出しは止まります。
* **ジョブの種類**: 書き出した入力には、Gaussian の route の行にも ORCA の `!` の行にもジョブのキーワード（`opt`、`freq` など）がないため、そのままでは 1 点計算です。実行する前に、目的のジョブのキーワードを足してください。
* **ORCA の前に `ORCAFF.prms` を確かめる**: `.inp` は `ORCAFF.prms` を絶対パスで参照します。`.inp` を実行する前に、ほかのマシンへ移した場合も含めて、このファイルがあることを確かめてください。
* **原子順のマーカー**: 書き出したファイルには `MLMM_REF_PDB_ORDER_V1_SHA256=<digest>` が入ります。`-i` の各原子の名前・番号などの識別欄のハッシュで、座標・occupancy・B-factor は含みません。[`oniom-import`](oniom-import.md) の `--ref-pdb` は、名前を書き戻す前にこのマーカーを照合します。
* **多重度**: `-m` に 1 未満の値を指定すると、コマンドラインの段階で拒否されます。
* **必要なもの**: Gaussian と ORCA は mlmm-toolkit に含まれません。別にインストールし、ライセンスを取得してください。
* **終了コード**: {ref}`終了コード <ja-exit-codes>`を参照してください。

---

## 関連ドキュメント

* [oniom-import](oniom-import.md) — 直した ONIOM 入力を XYZ と層付き PDB に戻す
* [mm-parm](mm-parm.md) — CMAP を含まないものも含めた Amber トポロジーの作成
* [define-layer](define-layer.md) — 全系の PDB に層の B-factor を書き込む
* [トラブルシューティング](troubleshooting.md) — 実行に失敗したときの対処
