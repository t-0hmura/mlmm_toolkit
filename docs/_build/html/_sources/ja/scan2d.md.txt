# `scan2d`

## 概要

> **要約:** 調和拘束と ML/MM 緩和による 2 距離（d1, d2）グリッドスキャンを実行します。`--spec`（YAML/JSON、推奨）または `--scan-lists` を使用します。

`mlmm scan2d` は `--max-step-size` を使用して 2 つの結合距離の線形グリッドを構築し、適切な拘束を適用して各グリッド点を緩和し、バイアスなしの ML/MM エネルギーを可視化用に記録します。スキャンはまず d1 を反復し d1 拘束のみで構造を緩和し、次に各 d1 値について d2 を反復し両方の拘束を適用します。

各グリッド点のエネルギーはバイアスなしで再評価され、PES グリッドとコンタープロットが作成されます。出力にはグリッド点ごとの XYZ スナップショット、PES をまとめた `surface.csv`、2D コンターマップ（`scan2d_map.png`）、底面投影付き 3D ランドスケープ（`scan2d_landscape.html`）が含まれます。


## 最小例

```bash
mlmm scan2d -i input.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --spec scan2d.yaml --print-parsed --out-dir ./result_scan2d/
```

## 出力の見方

- `result_scan2d/surface.csv`
- `result_scan2d/grid/point_i000_j000.xyz`
- `result_scan2d/scan2d_map.png` と `result_scan2d/scan2d_landscape.html`

## よくある例

1. YAML spec の解釈結果を先に確認する。
2. `--scan-lists` を使う。
3. `--dump` を有効にして d1 ごとの内側軌跡を保存する。

## 使用法

```bash
mlmm scan2d -i INPUT.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q CHARGE [-m SPIN] \
 [--spec scan2d.yaml | --scan-lists "[(I1,J1,LOW1,HIGH1),(I2,J2,LOW2,HIGH2)]"] \
 [--one-based|--zero-based] [--max-step-size FLOAT] [--bias-k FLOAT] \
 [--freeze-atoms "1,3,5"] [--relax-max-cycles INT] [--thresh PRESET] \
 [--dump/--no-dump] [--out-dir DIR] \
 [--preopt/--no-preopt] [--baseline {min|first}] [--zmin FLOAT] [--zmax FLOAT]
```

### 例

```bash
# ミニマル例（d1 と d2 の 2 つの範囲）
mlmm scan2d -i input.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --scan-lists "[(12,45,1.30,3.10),(10,55,1.20,3.20)]"

# TRJ ダンプ付き LBFGS スキャン、コンタープロットの固定カラースケール
mlmm scan2d -i input.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --scan-lists "[(12,45,1.30,3.10),(10,55,1.20,3.20)]" \
 --max-step-size 0.20 --dump --out-dir ./result_scan2d/ --preopt --baseline min \
 --zmin 0.0 --zmax 40.0
```

## ワークフロー

1. **入力と事前最適化** -- 酵素 PDB を読み込み、電荷/スピンを解決し、ML/MM 計算機（FAIR-Chem UMA + OpenMM）を構築し、`--preopt` の場合は任意でバイアスなし事前最適化を実行。
2. **グリッド構築** -- `--scan-lists` を 2 つの四つ組に解析し、インデックスを正規化（デフォルト 1 始まりまたは `"TYR,285,CA"` のような PDB 原子セレクター）。`ceil(|high - low| / h) + 1` 点の線形グリッドを構築（`h = --max-step-size`）。
3. **外側ループ（d1）** -- 各 d1 値について、**d1 拘束のみ**で系を緩和。
4. **内側ループ（d2）** -- 現在の d1 での各 d2 値について、最も近い収束済み構造から開始し**両方の拘束**で緩和。
5. **エネルギー評価** -- 各 (i, j) ペアで ML/MM エネルギーをバイアスなしで評価し `surface.csv` に記録。
6. **可視化** -- `scan2d_map.png`（2D コンター）と `scan2d_landscape.html`（3D サーフェス）を書き出し。`--zmin/--zmax` でカラースケールをクランプ。ベースライン: `--baseline min` は最小エネルギーをゼロに; `--baseline first` は (i=0, j=0) グリッド点をゼロに。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 入力酵素複合体 PDB（必須）。 | 必須 |
| `--real-parm7 PATH` | 酵素の Amber parm7 トポロジー（必須）。 | 必須 |
| `--model-pdb PATH` | ML 領域を定義する PDB。`--detect-layer` 有効時はオプション。 | _None_ |
| `--model-indices TEXT` | ML 領域のカンマ区切り原子インデックス（範囲指定可）。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` を 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer / --no-detect-layer` | 入力 PDB の B 因子から ML/MM レイヤーを検出。 | `True` |
| `-q, --charge INT` | ML 領域の総電荷。 | _None_ |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | _None_（デフォルト `1`） |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切りインデックス。 | _None_ |
| `--hess-cutoff FLOAT` | MM ヘシアン原子の距離カットオフ (A)。カットオフ指定で `--detect-layer` が無効化。 | _None_ |
| `--movable-cutoff FLOAT` | 可動 MM 原子の距離カットオフ (A)。 | _None_ |
| `--spec FILE` | `pairs`（2 四つ組）を持つ YAML/JSON 仕様。`one_based` を任意指定可能。 | 推奨 |
| `--scan-lists TEXT` | : 2 つの四つ組を含む Python リテラル: `"[(i1,j1,low1,high1),(i2,j2,low2,high2)]"`。インデックスは整数または PDB 原子セレクター。 | `--spec` の代替 |
| `--one-based / --zero-based` | `--scan-lists` の `(i,j)` インデックスを 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--print-parsed/--no-print-parsed` | `--spec`/`--scan-lists` 解釈後のペア情報を表示。 | `False` |
| `--max-step-size FLOAT` | ステップごとの最大距離増分 (A)。グリッド密度を決定。 | `0.20` |
| `--bias-k FLOAT` | 調和ウェル強度 k (eV/A^2)。 | `100.0` |
| `--relax-max-cycles INT` | バイアス緩和ごとの最大 LBFGS サイクル。 | `10000` |
| `--dump/--no-dump` | d1 スライスごとの内側 d2 スキャン TRJ を書き出し。 | `False` |
| `--out-dir TEXT` | 基本出力ディレクトリ。 | `./result_scan2d/` |
| `--thresh TEXT` | 収束プリセット（`gau_loose\|gau\|gau_tight\|gau_vtight\|baker\|never`）。 | _None_ |
| `--config FILE` | ベース YAML 設定ファイル（最初に適用）。 | _None_ |
| `--ref-pdb FILE` | `--input` が XYZ の場合の参照 PDB トポロジー。 | _None_ |
| `--preopt/--no-preopt` | スキャン前にバイアスなし事前最適化を実行。 | `True` |
| `--baseline {min,first}` | 相対エネルギーの基準（kcal/mol）。 | `min` |
| `--zmin FLOAT` | コンターカラースケールの下限（kcal/mol）。 | 自動スケール |
| `--zmax FLOAT` | コンターカラースケールの上限（kcal/mol）。 | 自動スケール |

## 出力

```
out_dir/ (デフォルト:./result_scan2d/)
├── surface.csv # PES グリッド: i, j, d1_A, d2_A, energy_hartree, energy_kcal, bias_converged
├── scan2d_map.png # 2D コンターマップ
├── scan2d_landscape.html # 3D サーフェス可視化（Plotly）
├── grid/
│ ├── point_i###_j###.xyz # 各 (i, j) ペアの緩和ジオメトリ
│ ├── point_i###_j###.pdb # PDB コンパニオン（入力が PDB の場合）
│ ├── preopt_i###_j###.xyz # 事前最適化構造（--preopt 時）
│ └── inner_path_d1_###_trj.xyz # d1 スライスごとの内側 d2 軌跡（--dump 時）
└── (stdout) # 進捗とエネルギーサマリー
```



ミニマル例（[opt](opt.md) と同じキーで拡張）:

```yaml
geom:
 coord_type: cart
 freeze_atoms: []
calc:
 charge: 0
 spin: 1
mlmm:
 real_parm7: real.parm7
 model_pdb: ml_region.pdb
opt:
 thresh: baker
 max_cycles: 10000
 dump: false
 out_dir:./result_scan2d/
lbfgs:
 max_step: 0.3
 out_dir:./result_scan2d/
bias:
 k: 100.0
```

## 注意事項

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- ML/MM 計算機（`mlmm_toolkit.mlmm_calc.mlmm`）は酵素複合体全体を保持します。ML 領域は `--model-pdb` から、Amber パラメータは `--real-parm7` から読み取られます。
- バイアスは最終エネルギー記録前に常に除去されるため、`surface.csv` はグリッド点間で直接比較可能です。
- 入力が PDB の場合、各グリッド点 XYZ と（存在する場合）内部パス TRJ は B 因子アノテーション付きで PDB ファイルにも変換されます: ML 領域原子 = 100.00、凍結原子 = 50.00、両方 = 150.00。
- `--scan-lists` の `i`/`j` エントリは整数インデックス（デフォルト 1 始まり）または `"TYR,285,CA"` のような PDB 原子セレクターが使用可能です。

---

## 関連項目

- [scan](scan.md) -- 1D 結合距離駆動スキャン
- [scan3d](scan3d.md) -- 3D 距離グリッドスキャン
- [opt](opt.md) -- 単一構造の構造最適化
- [all](all.md) -- エンドツーエンドワークフロー
