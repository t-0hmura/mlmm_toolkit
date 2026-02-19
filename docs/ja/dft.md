# `dft`

## 概要

> **要約:** PySCF/GPU4PySCF を使用して ML 領域の DFT 一点計算を実行し、MM エネルギーと再結合して ML(dft)/MM 総エネルギーを取得します。

`mlmm dft` は完全酵素 PDB から ML 領域を抽出し、リンク水素を付加して PySCF（または GPU4PySCF）による一点計算を実行します。DFT 評価後、PySCF 高レベルエネルギーと全系の MM 評価（REAL-low）および ML サブセットの MM 評価（MODEL-low）を組み合わせて **ML(dft)/MM 総エネルギー** を再計算します:

```
E_total = E_REAL_low + E_ML(DFT) - E_MODEL_low
```

GPU4PySCF バックエンドは利用可能な場合に自動的に有効化されます。それ以外は PySCF CPU が使用されます。デフォルトの汎関数/基底関数は `wb97m-v/6-31g**` です。

## 使用法

```bash
mlmm dft -i INPUT.pdb --real-parm7 real.parm7 --model-pdb model.pdb \
    -q CHARGE [-m SPIN] [--freeze-atoms "1,3,5"] [--func-basis "FUNC/BASIS"] \
    [--max-cycle N] [--conv-tol Eh] [--grid-level L] [--out-dir DIR] [--args-yaml FILE]
```

### 例

```bash
# 基本的な DFT 一点計算
mlmm dft -i enzyme.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1

# カスタム汎関数/基底関数、凍結原子、厳密な SCF
mlmm dft -i enzyme.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q -1 -m 2 \
    --func-basis "wb97m-v/def2-tzvpd" --freeze-atoms "1,3,5" --max-cycle 150 --conv-tol 1e-9
```

## ワークフロー

1. **入力処理** -- 完全酵素 PDB（`-i`）、Amber トポロジー（`--real-parm7`）、ML 領域定義（`--model-pdb` または `--model-indices` または `--detect-layer` による B 因子検出）を読み込みます。リンク水素は自動付加されます（C/N 親原子が 1.7 A 以内）。YAML で明示的な `link_mlmm` ペアが提供されない限り有効です。
2. **設定のマージ** -- デフォルト -> CLI -> YAML（`geom`、`calc`/`mlmm`、`dft` ブロック）。YAML 上書きは **CLI > YAML > デフォルト** の優先順位に従います。
3. **SCF 構築** -- `--func-basis` が汎関数と基底関数に解析されます。密度適合と非局所補正は PySCF/GPU4PySCF のデフォルトに従います。
4. **ML(dft)/MM 再結合** -- DFT が収束した後、全系（REAL-low）と ML サブセット（MODEL-low）の MM 評価が計算されます。結合エネルギーは Hartree と kcal/mol で報告されます。
5. **集団解析と出力** -- Mulliken、meta-Lowdin、IAO 電荷とスピン密度（UKS のみ）が結合エネルギーブロックとともに書き出されます。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 完全酵素 PDB ファイル（`.pdb` 必須）。 | 必須 |
| `--real-parm7 PATH` | 全系の Amber parm7 トポロジー。 | 必須 |
| `--model-pdb PATH` | ML 領域を定義する PDB（原子 ID が酵素 PDB と一致必須）。`--detect-layer` 有効時はオプション。 | _None_ |
| `--model-indices TEXT` | ML 領域のカンマ区切り原子インデックス（範囲指定可、例: `1-5`）。`--model-pdb` 省略時に使用。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` を 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer / --no-detect-layer` | 入力 PDB の B 因子（B=0/10/20）から ML/MM レイヤーを検出。 | `True` |
| `-q, --charge INT` | ML 領域の電荷。 | 必須 |
| `-m, --multiplicity INT` | ML 領域のスピン多重度 (2S+1)。 | `1` |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切りインデックス（例: `"1,3,5"`）。YAML `geom.freeze_atoms` とマージ。 | _None_ |
| `--func-basis TEXT` | 汎関数/基底関数ペア（`"FUNC/BASIS"`）。 | `wb97m-v/6-31g**` |
| `--max-cycle INT` | 最大 SCF 反復数。 | `100` |
| `--conv-tol FLOAT` | SCF 収束閾値 (Hartree)。 | `1e-9` |
| `--grid-level INT` | PySCF 数値積分グリッドレベル。 | `3` |
| `--out-dir DIR` | 出力ディレクトリ。 | `./result_dft/` |
| `--args-yaml FILE` | `geom`、`calc`/`mlmm`、`dft` の YAML 上書き。 | _None_ |

## 出力

```
out_dir/  (デフォルト: ./result_dft/)
├── ml_region_with_linkH.xyz    # DFT に使用された ML 領域座標（リンク水素付き）
├── result.yaml                 # DFT + ML(dft)/MM エネルギーサマリー、電荷、スピン密度
└── (stdout)                    # 整形された設定ブロックとエネルギーの出力
```

## YAML 設定（`--args-yaml`）

マッピングルートを受け付けます。`dft` セクション（およびオプションの `geom`、`calc`/`mlmm`）が存在する場合に適用されます。YAML 値が CLI 値を上書きします。

`dft` キー（括弧内はデフォルト）:
- `func_basis`（`"wb97m-v/6-31g**"`）: 結合 `FUNC/BASIS` 文字列。
- `conv_tol`（`1e-9`）: SCF 収束閾値 (Hartree)。
- `max_cycle`（`100`）: 最大 SCF 反復数。
- `grid_level`（`3`）: PySCF `grids.level`。
- `verbose`（`4`）: PySCF 冗長度 (0-9)。
- `out_dir`（`"./result_dft/"`）: 出力ディレクトリルート。

```yaml
geom:
  coord_type: cart
calc:
  charge: 0
  spin: 1
mlmm:
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
dft:
  func_basis: wb97m-v/6-31g**
  conv_tol: 1.0e-09
  max_cycle: 100
  grid_level: 3
  verbose: 4
  out_dir: ./result_dft/
```

## 注意事項

- リンク水素は自動検出されます（C/N 親原子が 1.7 A 以内）。YAML で明示的な `link_mlmm` ペアが提供されない限り有効です。サポートされていない親元素ではエラーが発生します。
- DFT オプション（汎関数/基底関数、SCF 制御）は `dft` キーの下で YAML 上書き可能です。
- 終了コード: SCF 収束時 `0`、未収束時 `3`、キーボード割り込み時 `130`、その他のエラー時 `1`。

---

## 関連項目

- [freq](freq.md) -- 振動解析（DFT 精密化の前に実行する場合が多い）
- [opt](opt.md) -- 単一構造の構造最適化
- [all](all.md) -- `--dft True` 付きエンドツーエンドワークフロー
