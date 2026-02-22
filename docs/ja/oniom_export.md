# `oniom-gaussian` / `oniom-orca`

## 概要

> **要約:** ML/MM システムを Gaussian または ORCA ONIOM 入力形式にエクスポートします。力場パラメータは Amber parm7 トポロジーファイルから直接抽出され、任意の Amber 互換力場（ff14SB、ff19SB、GAFF2 等）をサポートします。

### 概要
- **2 つのサブコマンド:** `mlmm oniom-gaussian`（Gaussian ONIOM `.com`/`.gjf`）と `mlmm oniom-orca`（ORCA QM/MM `.inp`）。
- **入力:** Amber parm7 トポロジー、任意の座標ファイル（PDB/XYZ）、任意の ML 領域 PDB。
- **元素検証:** parm7 と座標ファイル間のベストエフォートの元素配列チェック（原子順序は不変と仮定）。
- **用途:** mlmm_toolkit で準備したシステムを使用して Gaussian または ORCA で QM/MM 計算を実行する場合。

## 使用法

### oniom-gaussian

```bash
mlmm oniom-gaussian --parm7 real.parm7 [-i coords.pdb] [--model-pdb ml_region.pdb] \
    -o output.com [-q CHARGE] [-m MULT] [--method "B3LYP/6-31G(d,p)"] \
    [--near FLOAT] [--nproc INT] [--mem TEXT] \
    [--element-check|--no-element-check]
```

### oniom-orca

```bash
mlmm oniom-orca --parm7 real.parm7 [-i coords.pdb] [--model-pdb ml_region.pdb] \
    -o output.inp [-q CHARGE] [-m MULT] [--method "B3LYP D3BJ def2-SVP"] \
    [--total-charge INT] [--total-mult INT] \
    [--nproc INT] [--near FLOAT] [--orcaff PATH] \
    [--convert-orcaff|--no-convert-orcaff] [--element-check|--no-element-check]
```

### 例

```bash
# Gaussian ONIOM 入力を生成
mlmm oniom-gaussian --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
    -o system.com -q 0 -m 1

# ORCA QM/MM 入力を生成
mlmm oniom-orca --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
    -o system.inp -q 0 -m 1
```

## 説明

両方のサブコマンドは Amber parm7 トポロジーファイルを読み込み、外部プログラム用の QM/MM 入力ファイルを生成します。ワークフローは以下の通りです:

1. **トポロジーの読み込み**: ParmEd で parm7 ファイルを読み込み、原子タイプ、電荷、結合、角度、二面角、van der Waals パラメータを抽出。
2. **座標の読み込み**（任意）: `-i/--input` が提供された場合、PDB または XYZ 座標ファイルを読み込みます。`--element-check` が有効（デフォルト）な場合、元素配列が parm7 トポロジーと一致することを検証。
3. **QM 領域の同定**: `--model-pdb` が提供された場合、ML 領域 PDB の原子を全系とマッチングして QM（高）レイヤーを定義。
4. **入力の生成**: 適切な入力ファイルを書き出します:
   - **Gaussian**: `%nproc`、`%mem`、メソッド指定、レイヤー割り当て付き原子座標、コネクティビティを含む ONIOM 入力。QM 領域から `--near` Angstrom 以内の可動原子を同定。
   - **ORCA**: `%pal nprocs`、メソッド指定、QM 原子選択、MM パラメータ用 parm7 ファイル参照を含む QM/MM 入力。

## CLI オプション

### oniom-gaussian

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--parm7 PATH` | Amber parm7 トポロジーファイル。 | 必須 |
| `-i, --input PATH` | 座標ファイル（PDB/XYZ）; 原子順序は parm7 と一致必須。 | _None_ |
| `--element-check / --no-element-check` | 入力と parm7 間の元素配列を検証。 | `True` |
| `--model-pdb PATH` | QM 領域原子を定義する PDB ファイル。 | _None_ |
| `-o, --output PATH` | 出力 Gaussian 入力ファイル（`.com` または `.gjf`）。 | 必須 |
| `--method TEXT` | QM メソッドと基底関数。 | `B3LYP/6-31G(d,p)` |
| `-q, --charge INT` | QM 領域の電荷。 | `0` |
| `-m, --mult INT` | QM 領域の多重度。 | `1` |
| `--near FLOAT` | 可動原子の距離カットオフ (Angstrom)。 | `6.0` |
| `--nproc INT` | プロセッサ数。 | `8` |
| `--mem TEXT` | メモリ割り当て。 | `16GB` |

### oniom-orca

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--parm7 PATH` | Amber parm7 トポロジーファイル。 | 必須 |
| `-i, --input PATH` | 座標ファイル（PDB/XYZ）; 原子順序は parm7 と一致必須。 | _None_ |
| `--element-check / --no-element-check` | 入力と parm7 間の元素配列を検証。 | `True` |
| `--model-pdb PATH` | QM 領域原子を定義する PDB ファイル。 | _None_ |
| `-o, --output PATH` | 出力 ORCA 入力ファイル（`.inp`）。 | 必須 |
| `--method TEXT` | QM メソッドと基底関数。 | `B3LYP D3BJ def2-SVP` |
| `-q, --charge INT` | QM 領域の電荷。 | `0` |
| `-m, --mult INT` | QM 領域の多重度。 | `1` |
| `--total-charge INT` | 全 QM+MM 系の総電荷（`Charge_Total`）。 | トポロジ由来 |
| `--total-mult INT` | 全 QM+MM 系の総多重度（`Mult_Total`）。 | `--mult` と同じ |
| `--nproc INT` | プロセッサ数。 | `8` |
| `--near FLOAT` | B-factor レイヤーが無い場合に ActiveAtoms を決める距離カットオフ（Å）。 | `6.0` |
| `--orcaff PATH` | ORCAFF.prms へのパス。省略時は `<parm7_stem>.ORCAFF.prms` を使用/生成。 | 自動 |
| `--convert-orcaff/--no-convert-orcaff` | ORCAFF.prms が無い場合に `orca_mm -convff -AMBER` 自動実行を試行。 | `True` |

## 出力

```
<output>.com / <output>.gjf   # Gaussian ONIOM 入力（oniom-gaussian）
<output>.inp                  # ORCA QM/MM 入力（oniom-orca）
```

- コンソールに QM 原子数のサマリー、（ORCA の場合）入力と一緒に parm7 ファイルをコピーするリマインダーが出力されます。

## 注意事項

- トポロジー解析に `ParmEd`（`parmed`）が必要です。未インストール時は `ImportError` が発生します。
- 座標ファイルと parm7 間の元素検証はベストエフォートです。原子順序は不変と仮定されます。
- Gaussian の場合、`--near` カットオフは ONIOM 最適化中にどの MM 原子が移動可能かを決定します。

## 失敗時の最短復旧レシピ
1. エクスポート開始時に `ImportError: parmed` が出る。

```bash
python -c "import parmed; print(parmed.__version__)"
pip install parmed
```

2. 元素/原子順検証（`--element-check`）で失敗する。

```bash
mlmm oniom-gaussian --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
  -o system.com --element-check
```

3. QM/MM 境界やレイヤー割り当てが不自然。

```bash
mlmm define-layer -i pocket.pdb --real-parm7 real.parm7 -o layered.pdb
mlmm oniom-orca --parm7 real.parm7 -i layered.pdb --model-pdb ml_region.pdb \
  -o system.inp
```

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細な対処ガイド

- [mm_parm](mm_parm.md) -- AMBER トポロジー（parm7/rst7）の構築
- [define_layer](define_layer.md) -- ML/MM レイヤーの定義
- [opt](opt.md) -- ML/MM 構造最適化
