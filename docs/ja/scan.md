# `scan`

ML/MM 計算機を用い、調和拘束による結合距離スキャンで階層化された酵素 PDB の反応座標を駆動します。単一の出発構造から 1 つ以上の原子間距離を目標値まで駆動して反応軌跡の粗い候補を生成し、下流の MEP 精密化用の中間体/生成物候補を得たいときに使用します。`mlmm scan` は ML/MM 計算機（`mlmm.backends.mlmm_calc.mlmm`）による調和拘束付きの段階的な結合距離駆動スキャンを実行します。各ステップで一時的なターゲットを更新し、拘束ウェルを適用して LBFGS で構造を緩和します。ML/MM 計算機は MLIP バックエンド（デフォルト: UMA、`-b/--backend` で選択）と hessian_ff を結合します。`-s/--scan-lists` で YAML/JSON スペックファイル（推奨）またはインライン Python リテラルとしてターゲット距離を定義します。

## 実行例

コマンド形式:

```bash
mlmm scan -i INPUT.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q CHARGE [-m MULT] \
 [-s scan.yaml | -s "[(I,J,TARGET_ANG)]"] [options]
```

スペックファイルによるスキャン（`--print-parsed` を追加すると、解釈したスキャンスペックを検証し GPU 計算を実行せずに終了します）:

```bash
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan.yaml -o ./result_scan
```

インライン Python リテラル:

```bash
# インライン Python リテラル
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s "[(12,45,2.20)]"
```

ステージごとの軌跡を保存して確認する:

```bash
# ステージごとの軌跡を保存して確認する
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan.yaml --dump -o ./result_scan_dump
```

## 処理の流れ

1. `geom_loader` で構造を読み込み、CLI またはデフォルトから電荷/スピンを解決します。ML/MM 計算機に `--parm`、`--model-pdb`、`-q/--charge`、任意で `-m/--multiplicity` を提供します。
2. 任意でバイアスなし事前最適化（`--preopt`）を実行し、開始点を緩和します。
3. `-s/--scan-lists`（YAML/JSON スペックファイルまたはインラインリテラル）からステージターゲットを解析し、`(i, j)` インデックスを正規化します（デフォルトは 1 始まり）。入力が PDB の場合、各エントリは整数インデックスまたは `'TYR,285,CA'` のような原子セレクター文字列のいずれかで指定可能です。セレクターフィールドはスペース、カンマ、スラッシュ、バッククォート、バックスラッシュで区切ることができ、順序は任意です。
4. 結合ごとの変位を計算してステップに分割します:
 - スキャンタプル `[(i, j, target_A)]` に対し、`delta = target - current_distance_A` を計算。
 - `--max-step-size = h` の場合、ステージは `N = ceil(max(|delta|) / h)` 回のバイアス付き緩和を実行。
 - 各ペアの増分変化は `step_k = delta / N` (Å)。ステップ `s` での一時ターゲットは `r_k(s) = r_k(0) + s * step_k`。
5. すべてのステップを進み、調和ウェル `E_bias = sum 1/2 * k * (|r_i - r_j| - target_k)^2` を適用して LBFGS で極小化。`k` は `--bias-k`（eV/Å²）から取得され、Hartree/Bohr^2 に一度変換されます。座標は PySisyphus 用に Bohr で保存され、レポート時に内部変換されます。
6. 各ステージの最後のステップ後、任意でバイアスなし緩和（`--endopt`）を実行してから共有結合変化を報告し `result.*` ファイルを書き出します。
7. すべてのステージで繰り返します。任意の軌跡は `--dump` が `True` の場合のみダンプされます。

## 出力

各ステージは最終ジオメトリとバイアスステップ軌跡を `stage_XX/` 配下に書き出し、ルートに連結軌跡を生成します。最初に確認するファイルはステージ別の `result.pdb`（または `result.xyz`）と、常に生成される `scan_trj.xyz` / `scan.pdb` です。

```
out_dir/ (デフォルト:./result_scan/)
├─ scan_trj.xyz              # 全ステージ結合軌跡（常に書き出し）
├─ scan.pdb                  # 結合 PDB コンパニオン（PDB 入力のみ、常に書き出し）
├─ preopt/                   # --preopt が True の場合
│  ├─ result.xyz
│  └─ result.pdb             # PDB 入力のみ
└─ stage_XX/                 # ステージごとに 1 フォルダ（k = 01..K）
   ├─ result.xyz             # 最終（endopt 済みの可能性あり）ジオメトリ
   ├─ result.pdb             # 入力が PDB の場合
   ├─ scan_trj.xyz           # ステージ別バイアスステップフレーム（常に書き出し）
   └─ scan.pdb               # scan_trj.xyz の PDB 版（PDB 入力のみ、常に書き出し）
```

## CLI オプション

完全なフラグ一覧は生成された [コマンドリファレンス](../reference/commands/index.md) にあります。下表は説明が必要なオプションを扱います。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力 PDB（またはトポロジー用に `--ref-pdb` 付きの XYZ）。 | 必須 |
| `--parm PATH` | 完全 REAL 系の Amber prmtop。 | 必須 |
| `--model-pdb PATH` | ML 領域を定義する PDB（原子 ID）。`--detect-layer` 有効時または `--model-indices` 指定時は省略可能。 | _None_ |
| `--model-indices TEXT` | ML 領域原子インデックス（カンマ区切り、範囲指定可）。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` を 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer / --no-detect-layer` | 入力 PDB の B 因子から ML/MM レイヤーを検出。 | `True` |
| `-q, --charge INT` | ML 領域の総電荷。 | _None_（`-l` 未指定時は必須） |
| `-l, --ligand-charge TEXT` | 残基ごとの電荷マッピング（例: `GPP:-3,SAM:1`）。`-q` 省略時に合計電荷を導出。 | _None_ |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切り原子インデックス（YAML `geom.freeze_atoms` とマージ）。 | _None_ |
| `--hess-cutoff FLOAT` | Hessian 計算に含める MM 原子の ML 領域からの距離カットオフ (Å)。`--detect-layer` と併用可能。 | _None_ |
| `--movable-cutoff FLOAT` | 可動 MM 距離カットオフ (Å)。指定すると `--detect-layer` を無効化。 | _None_ |
| `-s, --scan-lists TEXT` | スキャンターゲット: YAML/JSON スペックファイルパス（自動検出）または `(i, j, target_A)` 三つ組もしくは `(i, j, start, end)` 四つ組（双方向スキャン）を含むインライン Python リテラル。各リテラルが 1 ステージ。単一フラグの後に複数リテラルを供給可能。`i`/`j` は整数インデックスまたは `"TYR,285,CA"` のような PDB 原子セレクターが使用可能。 | 必須 |
| `--one-based/--zero-based` | 原子インデックスを 1 始まり（既定）または 0 始まりとして解釈。 | `True`（1 始まり） |
| `--print-parsed/--no-print-parsed` | `-s/--scan-lists` 解釈後のステージ情報を表示。 | `False` |
| `--max-step-size FLOAT` | ステップごとのスキャン結合の最大変化量 (Å)。積分ステップ数を制御。 | `0.20` |
| `--bias-k FLOAT` | 調和バイアス強度 `k`（eV/Å²）。 | `300` |
| `--opt-mode {grad,hess,lbfgs,rfo,light,heavy}` | `mlmm all` からの転送互換オプション。現状の `scan` 緩和は mode に関わらず LBFGS を使用。 | _None_ |
| `--max-cycles INT` | 各バイアスステップおよび pre/end 最適化ステージの最大 LBFGS サイクル。 | `10000` |
| `--relax-max-cycles INT` | `--max-cycles` の互換エイリアス（指定時は上書き）。 | _None_ |
| `--preopt/--no-preopt` | スキャン前にバイアスなし最適化を実行。 | `False` |
| `--endopt/--no-endopt` | 各ステージ後にバイアスなし最適化を実行。 | `False` |
| `--dump/--no-dump` | ステップごとのオプティマイザー軌跡ファイルをダンプ。注: `scan_trj.xyz`/`scan.pdb` はこのフラグに関係なく常に書き出されます。 | `False` |
| `-o, --out-dir TEXT` | 出力ディレクトリルート。 | `./result_scan/` |
| `--thresh TEXT` | 収束プリセット（`gau_loose\|gau\|gau_tight\|gau_vtight\|baker\|never`）。 | _None_（`gau` を継承） |
| `--config FILE` | ベース YAML 設定ファイル（最初に適用）。 | _None_ |
| `--ref-pdb FILE` | `--input` が XYZ の場合の参照 PDB トポロジー。 | _None_ |
| `-b, --backend CHOICE` | ML 領域の MLIP バックエンド: `uma`、`orb`、`mace`、`aimnet2`。 | `uma` |
| `--embedcharge/--no-embedcharge` | xTB 点電荷埋め込み補正の有効化。MM 環境から ML 領域への静電的影響を考慮。 | `False` |
| `--embedcharge-cutoff FLOAT` | xTB 埋め込み用 MM 原子のカットオフ半径（Å）。 | `12.0` |
| `--cmap/--no-cmap` | model parm7 に CMAP（骨格クロスマップ二面角補正）を含めるかどうか。デフォルト: 無効（Gaussian ONIOM と同一）。 | `--no-cmap` |
| `--mm-backend [hessian_ff\|openmm]` | MM バックエンド（解析的 Hessian vs OpenMM 有限差分）。 | `hessian_ff` |
| `--link-atom-method [scaled\|fixed]` | リンク原子の配置: scaled（$g$ 因子）または固定 1.09/1.01 Å。 | `scaled` |
| `--out-json/--no-out-json` | 機械可読な `result.json` を `out_dir` に書き出す。 | `False` |
| `--dry-run/--no-dry-run` | オプションの検証と実行計画の表示のみ行い、スキャンは実行しない。`--help-advanced` に表示。 | `False` |
| `--convert-files/--no-convert-files` | PDB テンプレート利用可能時の XYZ/TRJ から PDB コンパニオン生成の切り替え。 | `True` |

## スキャン対象の構文

**YAML/JSON スペックフォーマット（推奨）**

`-s/--scan-lists` は YAML/JSON ファイルを自動検出します。ファイルパスを渡すとスペックモードになります:

```yaml
one_based: true # 任意; デフォルトは CLI の --one-based/--zero-based
stages:
 - [[12, 45, 2.20]]
 - [[10, 55, 1.35], [23, 34, 1.80]]
```

- `stages` は必須です。
- 各ステージは `(i, j, target_A)` の三つ組のリストです。
- インデックスは整数または PDB セレクター（PDB 入力時）が使用可能で、インラインリテラルと同じです。

**インラインリテラルフォーマット**

`-s/--scan-lists` がファイルパスでない値を受け取ると、**Python リテラル**文字列として評価されます。シェルクォートに注意してください。

各リテラルは三つ組 `(atom1, atom2, target_A)` の Python リストです:

```
-s '[(atom1, atom2, target_A),...]'
```

- シェルが括弧やスペースを解釈しないよう、リテラル全体を**シングルクォート**で囲んでください。
- 各三つ組は `atom1`--`atom2` 間の距離を `target_A` に向けて駆動します。
- 1 つのリテラル = 1 つの**ステージ**です。複数ステージの場合、**単一の** `-s/--scan-lists` フラグの後に複数リテラルを渡します（フラグを繰り返さないでください）。

原子は**整数インデックス**または **PDB セレクター文字列**で指定できます:

| 方法 | 例 | 備考 |
| --- | --- | --- |
| 整数インデックス | `(1, 5, 2.0)` | デフォルトは 1 始まり（`--one-based`） |
| PDB セレクター | `("TYR,285,CA", "MMT,309,C10", 2.0)` | 残基名、残基番号、原子名 |

PDB セレクターのトークンは、カンマ `,`、スペース、スラッシュ `/`、バッククォート `` ` ``、バックスラッシュ `\` のいずれかで区切れます。トークンの順序は自由です。

```bash
# 以下はすべて同じ原子を指定:
"TYR,285,CA"
"TYR 285 CA"
"TYR/285/CA"
"285,TYR,CA" # 順序は自由
```

クォート規則:

```bash
# 正しい: リスト全体をシングルクォート、内側のセレクター文字列をダブルクォート
-s '[("TYR,285,CA","MMT,309,C10",1.35)]'

# 正しい: 整数インデックスは内側のクォート不要
-s '[(1, 5, 2.0)]'

# 非推奨: 外側をダブルクォートにすると内側のクォートをエスケープする必要あり
-s "[(\"TYR,285,CA\",\"MMT,309,C10\",1.35)]"
```

単一の `-s/--scan-lists` フラグの後に複数リテラルを渡します。各リテラルが 1 ステージになります:

```bash
# ステージ 1: 1 つの結合を 1.35 Å に駆動
# ステージ 2: 2 つの結合を同時に駆動
-s \
 '[("TYR,285,CA","MMT,309,C10",1.35)]' \
 '[("TYR,285,CA","MMT,309,C10",2.20),("TYR,285,CB","MMT,309,C11",1.80)]'
```

ステージは順次実行され、各ステージは前のステージの緩和結果から開始します。**`-s/--scan-lists` フラグを繰り返さないでください** -- 単一のフラグの後にすべてのステージリテラルを供給してください。

**協奏的スキャン vs 段階的スキャン**

渡すリテラルの数で、座標を同時に駆動するか（協奏的）順次駆動するか（段階的）が決まります:

| 形式 | 指定方法 | 意味 | 機構を事前に要する? |
| --- | --- | --- | --- |
| 協奏的 | **単一**のリテラルに複数の `(i, j, target)` タプルを含める | すべての座標を 1 ステージ内で同時に駆動 | 不要 |
| 段階的 | **複数**のリテラル（ステージごとに 1 つ） | 各ステージが独自の拘束付き緩和で、`stage_NN/` に書き出される | 必要 — ステージごとに機構を定義 |

```bash
# 協奏的: 1 ステージ、2 つの距離を同時に駆動
mlmm scan -i r.pdb --parm enzyme.parm7 -l 'LIG:Q' \
    -s '[(1,5,1.40),(7,9,1.60)]' -o result_concerted

# 段階的: 2 つの順次ステージ
mlmm scan -i r.pdb --parm enzyme.parm7 -l 'LIG:Q' \
    -s '[(1,5,1.40)]' \
       '[(7,9,0.95)]' -o result_staged
```

協奏的スキャンは機構の分解を必要としません — [`path-search`](path-search.md) が多段の自動セグメント化を行ってくれます。段階的スキャンは事前に機構の定義が必要ですが、**機構が分かっている場合は段階的スキャンの方がステップごとの制御がクリーンで、一般に推奨されます。**（4-tuple は双方向スキャン用に 2 ステージへ展開されます。）

**双方向スキャン（4-tuple）**

3-tuple `(i, j, target)` の代わりに **4-tuple** `(i, j, start, end)` を指定すると、現在の構造から両方向にスキャンします。CLI は各 4-tuple を自動的に 2 ステージに展開します:

1. **パス 1:** `i`--`j` の距離を現在の値から `start` に向けて駆動。
2. **パス 2:** 初期構造を復元し、`i`--`j` の距離を `end` に向けて駆動。

結合軌跡は `start → 初期構造 → end` の順に連結され、出発構造を通る連続的な経路が得られます。

```bash
# 双方向スキャン: 結合 12--45 を現在の構造から
# 1.35 Å（パス 1）と 2.50 Å（パス 2）に向けて駆動
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s '[(12, 45, 1.35, 2.50)]'
```

これは 2 つの手動ステージの間にジオメトリリセットを行うのと同等ですが、スクリプトを書く必要がありません。同じリテラル内で 3-tuple と 4-tuple を混在させることもできます。

## バリアの方向を読む

読み取るバリアは、スキャンがどちらの端点から始まったかに依存します。スキャン（またはそれが種とする経路）が**生成物から始まる**場合、生のバリア値は**逆方向**です。

| 量 | 式 |
| --- | --- |
| 順方向バリア | `E(TS) − E(reactant)` |
| 逆方向バリア（生成物始点の生の値） | `E(TS) − E(product)` |

これは*読み取り時*の解釈であり、CLI フラグではありません。どちらの端点が反応物か生成物かは、スキャン方向を信頼せず、IRC の `segments/seg_NN/{reactant,product}.pdb` を読んで必ず確認してください。生成物始点のキャンペーンでは、欲しい順方向バリアは `E(TS) − E(reactant)` であり、生成物始点に対して印字される値ではありません。

## YAML 設定

スキャンは共有の `geom`（`coord_type`、`freeze_atoms`）、`calc` / `mlmm`（ML/MM 計算機設定）、`opt` / `lbfgs`（オプティマイザー）の各セクションに加え、`bias`（`k`、調和強度（eV/Å²））と MLIP ベースの結合変化検出用 `bond` セクションを読み込みます。

- `coord_type`: 座標タイプ（デカルト vs dlc 内部座標）。
- `freeze_atoms`: CLI `--freeze-atoms` とマージされる 1 始まり凍結原子。

### セクション `calc` / `mlmm`
- ML/MM 計算機の設定: `charge`、`spin`、`backend`、`embedcharge`、MLIP モデル設定、`device`、近傍半径、ヘシアンオプション等。

### セクション `opt` / `lbfgs`
- オプティマイザー設定: `thresh`、`max_cycles`、`print_every`、ステップ制御、ラインサーチ、ダンプフラグ。

### セクション `bias`
- `k`（`300`）: 調和強度（eV/Å²）。

### セクション `bond`
- MLIP ベースの結合変化検出:
 - `device`（`"auto"`）: グラフ分析用の MLIP デバイス。
 - `bond_factor`（`1.20`）: カットオフ用の共有結合半径スケーリング。
 - `margin_fraction`（`0.05`）: 比較用の許容分数。
 - `delta_fraction`（`0.05`）: 結合形成/切断をフラグする最小相対変化。

全スキーマ（すべてのキーとデフォルト）: [YAML リファレンス](yaml-reference.md)。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド
- [scan2d](scan2d.md) -- 2D 距離グリッドスキャン
- [scan3d](scan3d.md) -- 3D 距離グリッドスキャン
- [opt](opt.md) -- 単一構造の構造最適化
- [all](all.md) -- 単一構造入力の `--scan-lists` 付き一気通貫ワークフロー
- [path-search](path-search.md) -- スキャン端点を中間体として使用する MEP 探索
