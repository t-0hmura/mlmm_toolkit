# `scan`

## 概要

> **要約:** ML/MM 計算機を使用して、調和拘束による結合距離スキャンで反応座標を駆動します。`--scan-lists` でターゲット距離を定義します。複数ステージは順次実行され、各ステージは前ステージの緩和構造から開始します。

### 概要
- **用途:** 単一構造があり、特定の距離を変化させて反応経路を探索したい場合に使用（多くの場合 `path-search`/`path-opt` の前段階）。
- **入力:** 1 つの構造 + 1 つ以上の `--scan-lists` リテラル（各リテラル = 1 ステージ）。
- **デフォルト:** LBFGS オプティマイザー、`--preopt True`、`--endopt True`、`--max-step-size 0.20` A。
- **出力:** ステージごとの `result.xyz`（+ 任意で `.pdb`）、`--dump True` 時は連結軌跡。
- **注記:** `--scan-lists` は **Python リテラル**として解析されます。クオートとエスケープに注意してください（例を参照）。

`mlmm scan` は ML/MM 計算機（`mlmm_toolkit.mlmm_calc.mlmm`）を使用して調和拘束による段階的な結合距離駆動スキャンを実行します。各ステップで一時的なターゲットが更新され、拘束ウェルが適用され、LBFGS で構造が緩和されます。ML/MM 計算機は FAIR-Chem UMA と OpenMM をリンク原子なしで結合します。

単一フラグの後に複数の `--scan-lists` リテラルを提供すると、ステージは順次実行され、各ステージは前ステージの緩和構造から開始します。バイアス付きウォークの後、任意のバイアスなし事前/事後最適化（`--preopt`、`--endopt`）でジオメトリを整えてから `result.*` をディスクに書き出します。設定の優先順位は **CLI > YAML > 内部デフォルト** です。YAML セクションは `geom`、`calc`/`mlmm`、`opt`、`lbfgs`、`bias`、`bond` を含みます。

## 使用法

```bash
mlmm scan -i INPUT.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
          -q CHARGE [-m MULT] --scan-lists "[(I,J,TARGET_ANG)]" [options]
```

### 例

```bash
# 1 つの結合を 1.6 から 2.2 A にプッシュする単一ステージ
mlmm scan -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 --scan-lists "[(12,45,2.20)]"

# ダンプ付き 2 ステージ、凍結原子、YAML 上書き
mlmm scan -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
    -q -1 -m 1 --freeze-atoms "1,3,5" --scan-lists "[(12,45,2.20)]" \
    "[(10,55,1.35),(23,34,1.80)]" --max-step-size 0.20 --dump True \
    --args-yaml params.yaml --out-dir ./result_scan/
```

## ワークフロー

1. `geom_loader` で構造を読み込み、CLI またはデフォルトから電荷/スピンを解決します。ML/MM 計算機に `--real-parm7`、`--model-pdb`、`-q/--charge`、任意で `-m/--multiplicity` を提供します。
2. 任意でバイアスなし事前最適化（`--preopt True`）を実行し、開始点を緩和します。
3. `--scan-lists` で提供された各ステージリテラルについて、`(i, j)` インデックスを解析して正規化します（デフォルトは 1 始まり）。入力が PDB の場合、各エントリは整数インデックスまたは `'TYR,285,CA'` のような原子セレクター文字列のいずれかで指定可能です。セレクターフィールドはスペース、カンマ、スラッシュ、バッククォート、バックスラッシュで区切ることができ、順序は任意です。
4. 結合ごとの変位を計算してステップに分割します:
   - スキャンタプル `[(i, j, target_A)]` に対し、`delta = target - current_distance_A` を計算。
   - `--max-step-size = h` の場合、ステージは `N = ceil(max(|delta|) / h)` 回のバイアス付き緩和を実行。
   - 各ペアの増分変化は `delta_k = delta_k / N` (A)。ステップ `s` での一時ターゲットは `r_k(s) = r_k(0) + s * delta_k`。
5. すべてのステップを進み、調和ウェル `E_bias = sum 1/2 * k * (|r_i - r_j| - target_k)^2` を適用して LBFGS で極小化。`k` は `--bias-k`（eV/A^2）から取得され、Hartree/Bohr^2 に一度変換されます。座標は PySisyphus 用に Bohr で保存され、レポート時に内部変換されます。
6. 各ステージの最後のステップ後、任意でバイアスなし緩和（`--endopt True`）を実行してから共有結合変化を報告し `result.*` ファイルを書き出します。
7. すべてのステージで繰り返します。任意の軌跡は `--dump` が `True` の場合のみダンプされます。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力 PDB（またはトポロジー用に `--ref-pdb` 付きの XYZ）。 | 必須 |
| `--real-parm7 PATH` | 完全 REAL 系の Amber prmtop。 | 必須 |
| `--model-pdb PATH` | ML 領域を定義する PDB（原子 ID）。 | 必須 |
| `-q, --charge INT` | ML 領域の総電荷。 | 必須 |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `--scan-lists TEXT` | `(i, j, target_A)` タプルを含む Python リテラル。各リテラルが 1 ステージ。単一フラグの後に複数リテラルを供給可能。`i`/`j` は整数インデックスまたは `"TYR,285,CA"` のような PDB 原子セレクターが使用可能。 | 必須 |
| `--zero-based` | 原子インデックスを 1 始まりではなく 0 始まりとして解釈。 | `False`（1 始まり） |
| `--max-step-size FLOAT` | ステップごとのスキャン結合の最大変化量 (A)。積分ステップ数を制御。 | `0.20` |
| `--bias-k FLOAT` | 調和バイアス強度 `k`（eV/A^2）。 | `100` |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切り原子インデックス（YAML `geom.freeze_atoms` とマージ）。 | _None_ |
| `--preopt {True\|False}` | スキャン前にバイアスなし最適化を実行。 | `True` |
| `--endopt {True\|False}` | 各ステージ後にバイアスなし最適化を実行。 | `True` |
| `--dump {True\|False}` | 連結バイアス軌跡（`scan.trj`/`scan.pdb`）をダンプ。 | `False` |
| `--out-dir TEXT` | 出力ディレクトリルート。 | `./result_scan/` |
| `--args-yaml FILE` | `geom`、`calc`/`mlmm`、`opt`、`lbfgs`、`bias`、`bond` の YAML 上書き。 | _None_ |

## 出力

```
out_dir/ (デフォルト: ./result_scan/)
├─ preopt/                   # --preopt が True の場合
│  ├─ result.xyz
│  └─ result.pdb             # PDB 入力のみ
└─ stage_XX/                 # ステージごとに 1 フォルダ（k = 01..K）
    ├─ result.xyz            # 最終（endopt 済みの可能性あり）ジオメトリ
    ├─ result.pdb            # 入力が PDB の場合
    ├─ scan.trj              # --dump True 時のバイアスステップフレーム連結
    └─ scan.pdb              # scan.trj の PDB 版（PDB 入力のみ）
```

## YAML 設定（`--args-yaml`）

YAML ルートはマッピングでなければなりません。設定の優先順位は **CLI > YAML > 内部デフォルト** です。

### セクション `geom`

- `coord_type`: 座標タイプ（デカルト vs dlc 内部座標）。
- `freeze_atoms`: CLI `--freeze-atoms` とマージされる 0 始まり凍結原子。

### セクション `calc` / `mlmm`

- ML/MM 計算機の設定: `charge`、`spin`、UMA `model`、`task_name`、`device`、近傍半径、ヘシアンオプション等。

### セクション `opt` / `lbfgs`

- オプティマイザー設定: `thresh`、`max_cycles`、`print_every`、ステップ制御、ラインサーチ、ダンプフラグ。

### セクション `bias`

- `k`（`100`）: 調和強度（eV/A^2）。

### セクション `bond`

- UMA ベースの結合変化検出:
  - `device`（`"cuda"`）: グラフ分析用の UMA デバイス。
  - `bond_factor`（`1.20`）: カットオフ用の共有結合半径スケーリング。
  - `margin_fraction`（`0.05`）: 比較用の許容分数。
  - `delta_fraction`（`0.05`）: 結合形成/切断をフラグする最小相対変化。

---

## 関連項目

- [scan2d](scan2d.md) -- 2D 距離グリッドスキャン
- [scan3d](scan3d.md) -- 3D 距離グリッドスキャン
- [opt](opt.md) -- 単一構造の構造最適化
- [all](all.md) -- 単一構造入力の `--scan-lists` 付きエンドツーエンドワークフロー
- [path-search](path_search.md) -- スキャン端点を中間体として使用する MEP 探索
