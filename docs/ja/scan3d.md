# `scan3d`

## 概要

> **要約:** 調和拘束と ML/MM 緩和による 3 距離（d1, d2, d3）グリッドスキャンを実行します。`--spec`（YAML/JSON、推奨）または legacy の `--scan-lists` を使用します。

### 概要
- **入力:** 1 つの完全酵素 PDB + `--spec scan3d.yaml`（推奨）または **1 つ**の legacy `--scan-lists` リテラル（3 つの四つ組）。
- **グリッド順序:** d1 が最初にスキャンされ、次に各 d1 値について d2（両方の拘束が有効）、次に各 (d1, d2) について d3（3 つの拘束すべてが有効）。
- **エネルギー:** 記録されるエネルギーは**バイアスなし**で評価されるため、グリッド点間で直接比較可能。
- **出力:** `surface.csv`、`grid/` 下のグリッド点ごとのジオメトリ、HTML アイソサーフェスプロット（`scan3d_density.html`）。
- **注意:** 3D グリッドは急速に大きくなります。まず粗い `--max-step-size` または小さい範囲を検討してください。

`mlmm scan3d` は d1、d2、d3 のネストループを実行し、ML/MM 計算機（`mlmm_toolkit.mlmm_calc.mlmm`）を使用して適切な拘束で各点を緩和します。ML 領域は `--model-pdb` から、Amber パラメータは `--real-parm7` から読み取られ、オプティマイザーは PySisyphus LBFGS です。

設定セクション（`geom`、`calc`/`mlmm`、`opt`、`lbfgs`、`bias`）は `--config`（ベース）と `--override-yaml`（最終上書き）で 2 層指定できます。`--args-yaml` は `--override-yaml` の legacy エイリアスです。優先順位は **内部デフォルト < `--config` < `--override-yaml` < 明示 CLI** です。YAML と `--freeze-atoms` で提供される凍結原子はオプティマイザーと ML/MM 計算機の間で共有されます。

## 最小例

```bash
mlmm scan3d -i input.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 --spec scan3d.yaml --print-parsed --out-dir ./result_scan3d/
```

## 出力の見方

- `result_scan3d/surface.csv`
- `result_scan3d/grid/point_i000_j000_k000.xyz`
- `result_scan3d/scan3d_density.html`

## よくある例

1. YAML spec の `pairs` を先に検証する。
2. 後方互換のため legacy `--scan-lists` を使う。
3. `--dump` を有効にして `(d1,d2)` スライスごとの d3 軌跡を保存する。

## 使用法

```bash
mlmm scan3d -i INPUT.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
    -q CHARGE [-m MULT] \
    [--spec scan3d.yaml | --scan-lists "[(I1,J1,LOW1,HIGH1),(I2,J2,LOW2,HIGH2),(I3,J3,LOW3,HIGH3)]"] \
    [--one-based|--zero-based] [--max-step-size FLOAT] [--bias-k FLOAT] \
    [--freeze-atoms "1,3,5"] [--relax-max-cycles INT] [--thresh PRESET] \
    [--dump/--no-dump] [--out-dir DIR] \
    [--config FILE] [--override-yaml FILE | --args-yaml FILE] \
    [--preopt/--no-preopt] [--baseline {min|first}] [--zmin FLOAT] [--zmax FLOAT]
```

### 例

```bash
# ミニマルな 3 距離スキャン
mlmm scan3d -i input.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 --scan-lists "[(12,45,1.30,3.10),(10,55,1.20,3.20),(15,60,1.10,3.00)]"

# 事前最適化とカスタム出力ディレクトリ付き
mlmm scan3d -i input.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 --scan-lists "[(12,45,1.30,3.10),(10,55,1.20,3.20),(15,60,1.10,3.00)]" \
    --max-step-size 0.20 --dump --out-dir ./result_scan3d/ \
    --preopt --baseline min
```

## ワークフロー

1. `geom_loader` で構造を読み込み、CLI から電荷/スピンを解決し、`--preopt` の場合は任意でバイアスなし事前最適化を実行。
2. 単一の `--scan-lists` リテラルを解析して 3 つの四つ組にします（デフォルト 1 始まりインデックス、`--zero-based` 指定時は 0 始まり）。PDB 入力の場合、各原子エントリは整数インデックスまたは `"TYR,285,CA"` のようなセレクター文字列が使用可能。区切り文字はスペース、カンマ、スラッシュ、バッククォート、バックスラッシュ。
3. 外側ループ `d1[i]`: d1 拘束のみで緩和。d1 値が最も近い以前のスキャン済みジオメトリから開始。
4. 中間ループ `d2[j]`: d1 と d2 の拘束で緩和。最も近い (d1, d2) ジオメトリから開始。
5. 内側ループ `d3[k]`: 3 つの拘束すべてで緩和。バイアスなしエネルギーを測定（評価時にバイアス除去）し、拘束ジオメトリと収束フラグを書き出し。
6. スキャン完了後、`surface.csv` を組み立て、kcal/mol ベースラインシフト（`--baseline {min|first}`）を適用し、3D RBF 補間アイソサーフェスプロット（`scan3d_density.html`）を生成（`--zmin/--zmax` を尊重）。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 完全酵素 PDB（リンク原子なし）。 | 必須 |
| `--real-parm7 PATH` | 完全酵素の Amber parm7 トポロジー。 | 必須 |
| `--model-pdb PATH` | ML 領域を定義する PDB。 | _None_ |
| `--model-indices TEXT` | 明示的な ML 領域原子インデックス（`--model-pdb` の代替）。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` のインデックス規約。 | `True`（1 始まり） |
| `--detect-layer / --no-detect-layer` | B 因子から ML/MM レイヤーを自動検出。 | `True` |
| `-q, --charge INT` | ML 領域の総電荷。 | 必須 |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `--freeze-atoms TEXT` | 1 始まりカンマ区切りの凍結原子インデックス。 | _None_ |
| `--hess-cutoff FLOAT` | Hessian-MM レイヤーのカットオフ距離。 | _None_ |
| `--movable-cutoff FLOAT` | Movable-MM レイヤーのカットオフ距離。 | _None_ |
| `--spec FILE` | `pairs`（3 四つ組）を持つ YAML/JSON 仕様。`one_based` を任意指定可能。 | 推奨 |
| `--scan-lists TEXT` | legacy: 3 つの四つ組 `(i,j,low,high)` を含む単一 Python リテラル。`i`/`j` は整数インデックスまたは PDB 原子セレクター。 | `--spec` の代替 |
| `--one-based / --zero-based` | `(i, j)` インデックスを 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--print-parsed/--no-print-parsed` | `--spec`/`--scan-lists` 解釈後のペア情報を表示。 | `False` |
| `--max-step-size FLOAT` | ステップごとの最大距離増分 (A)。グリッド密度を制御。 | `0.20` |
| `--bias-k FLOAT` | 調和ウェル強度 k (eV/A^2)。 | `100.0` |
| `--relax-max-cycles INT` | バイアス緩和ごとの最大オプティマイザーサイクル。 | `10000` |
| `--dump/--no-dump` | (d1, d2) スライスごとの内側 d3 スキャン TRJ を書き出し。 | `False` |
| `--out-dir TEXT` | グリッドとプロットの出力ディレクトリルート。 | `./result_scan3d/` |
| `--thresh TEXT` | 収束プリセット上書き（`gau_loose`、`gau`、`gau_tight`、`gau_vtight`、`baker`、`never`）。 | _None_ |
| `--config FILE` | ベース YAML 設定ファイル（最初に適用）。 | _None_ |
| `--override-yaml FILE` | 最終 YAML 上書きファイル（YAML レイヤーの最優先）。 | _None_ |
| `--args-yaml FILE` | `--override-yaml` の legacy エイリアス。 | _None_ |
| `--ref-pdb FILE` | 非 PDB 入力用の参照 PDB トポロジー。 | _None_ |
| `--preopt/--no-preopt` | スキャン前にバイアスなし最適化を実行。 | `True` |
| `--baseline {min,first}` | kcal/mol エネルギーをグローバル最小値または `(i,j,k)=(0,0,0)` がゼロになるようシフト。 | `min` |
| `--zmin FLOAT` | アイソサーフェスカラーバンドの手動下限（kcal/mol）。 | 自動スケール |
| `--zmax FLOAT` | アイソサーフェスカラーバンドの手動上限（kcal/mol）。 | 自動スケール |

## 出力

```
out_dir/ (デフォルト: ./result_scan3d/)
  surface.csv                          # グリッドメタデータ（d1, d2, d3, energy, convergence）
  scan3d_density.html                  # 3D エネルギーアイソサーフェス可視化
  grid/point_i###_j###_k###.xyz        # 各グリッド点の緩和ジオメトリ
  grid/point_i###_j###_k###.pdb        # PDB コンパニオン（B 因子: ML=100, frozen=50, both=150）
  grid/inner_path_d1_###_d2_###.trj    # --dump が True の場合のみ
```

## 注意事項

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes_common_errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- ML/MM 計算機（`mlmm_toolkit.mlmm_calc.mlmm`）は 1D/2D スキャンと同じ `HarmonicBiasCalculator` を再利用します。
- `--baseline` のデフォルトはグローバル最小値です。`--baseline first` は `(i,j,k)=(0,0,0)` グリッド点が存在する場合にそれを基準にします。
- 3D 可視化は 50x50x50 グリッドの RBF 補間を使用し、半透明のステップカラーアイソサーフェスを表示します。
- 入力が PDB の場合、各グリッド点 XYZ と（存在する場合）内部パス TRJ は入力 PDB をテンプレートとして PDB ファイルにも変換されます。B 因子は `opt` ツールと一貫してアノテーションされます: ML 領域原子 = 100.00、凍結原子 = 50.00、両方 = 150.00。
- プロットのカラースケールは `--zmin/--zmax` でクランプして、スキャン間で一貫した比較が可能です。

---

## 関連項目

- [scan](scan.md) -- 1D 結合距離駆動スキャン
- [scan2d](scan2d.md) -- 2D 距離グリッドスキャン
- [opt](opt.md) -- 構造最適化（スキャン前に実行する場合が多い）
- [all](all.md) -- エンドツーエンドワークフロー
- [YAML リファレンス](yaml_reference.md) -- スキャンの完全な設定オプション
