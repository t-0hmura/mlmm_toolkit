# `dft`

## 概要

> **要約:** PySCF/GPU4PySCF を使用して ML 領域の DFT 一点計算を実行し、MM エネルギーと再結合して ML(dft)/MM 総エネルギーを取得します。結果にはエネルギーと集団解析（Mulliken、meta-Lowdin、IAO 電荷）が含まれます。

`mlmm dft` は完全酵素 PDB から ML 領域を抽出し、リンク水素を付加して PySCF（または GPU4PySCF）による一点計算を実行します。DFT 評価後、PySCF 高レベルエネルギーと全系の MM 評価（REAL-low）および ML サブセットの MM 評価（MODEL-low）を組み合わせて **ML(dft)/MM 総エネルギー** を再計算します:

```
E_total = E_REAL_low + E_ML(DFT) - E_MODEL_low
```

GPU4PySCF バックエンドは利用可能な場合に自動的に有効化されます。それ以外は PySCF CPU が使用されます。デフォルトの汎関数/基底関数は `wb97m-v/def2-tzvpd` です。

## 最小例

```bash
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_dft
```

## 出力の見方

- `result_dft/ml_region_with_linkH.xyz`
- `result_dft/result.yaml`
- 標準出力の ML(dft)/MM 合成エネルギー表示

## よくある例

1. 汎関数/基底関数を変更して一点計算する。

```bash
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --func-basis "wb97m-v/def2-tzvpd" --out-dir ./result_dft_tz
```

2. ML/MM 側で凍結原子を指定して DFT を実行する。

```bash
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q -1 -m 2 --freeze-atoms "1,3,5" --out-dir ./result_dft_freeze
```

3. SCF 収束を厳しくして反復回数を増やす。

```bash
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --conv-tol 1e-10 --max-cycle 200 --out-dir ./result_dft_tight
```

## ワークフロー

1. **入力処理** -- 完全酵素 PDB（`-i`）、Amber トポロジー（`--parm`）、ML 領域定義（`--model-pdb` または `--model-indices` または `--detect-layer` による B 因子検出）を読み込みます。リンク水素は自動付加されます（C/N 親原子が 1.7 Å 以内）。YAML で明示的な `link_mlmm` ペアが提供されない限り有効です。
2. **SCF 構築** -- `--func-basis` が汎関数と基底関数に解析されます。密度適合は PySCF デフォルトで自動有効化されます。GPU4PySCF バックエンドは利用可能な場合に使用され、それ以外は CPU PySCF が使用されます。`--embedcharge` が有効な場合、Amber トポロジーの MM 点電荷が `pyscf.qmmm.mm_charge()` を介して QM ハミルトニアンに埋め込まれ、DFT 波動関数が MM 環境で自己無撞着に分極します。
3. **ML(dft)/MM 再結合** -- DFT が収束した後、全系（REAL-low）と ML サブセット（MODEL-low）の MM 評価が計算されます。結合エネルギーは Hartree と kcal/mol で報告されます。
4. **集団解析と出力** -- Mulliken、meta-Lowdin、IAO 電荷とスピン密度（UKS のみ）が結合エネルギーブロックとともに `result.yaml` に書き出されます。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-b, --backend CHOICE` | ONIOM 低レベル再結合用 MLIP バックエンド: `uma`（デフォルト）、`orb`、`mace`、`aimnet2`。 | `uma` |
| `--embedcharge/--no-embedcharge` | 静電埋め込みの有効化: Amber トポロジーの MM 点電荷を PySCF QM ハミルトニアンに追加し、DFT 波動関数が MM 環境により分極。 | `False` |
| `--embedcharge-cutoff FLOAT` | xTB 埋め込み用 MM 原子のカットオフ半径（Å）。 | `12.0` |
| `--cmap/--no-cmap` | model parm7 に CMAP（骨格クロスマップ二面角補正）を含めるかどうか。デフォルト: 無効（Gaussian ONIOM と同一）。 | `--no-cmap` |
| `-i, --input PATH` | 完全酵素構造ファイル（PDB または XYZ）。XYZ の場合は `--ref-pdb` でトポロジーを指定。 | 必須 |
| `--parm PATH` | 全系の Amber parm7 トポロジー。 | 必須 |
| `--model-pdb PATH` | ML 領域を定義する PDB（原子 ID が酵素 PDB と一致必須）。`--detect-layer` 有効時はオプション。 | _None_ |
| `--model-indices TEXT` | ML 領域のカンマ区切り原子インデックス（範囲指定可、例: `1-5`）。`--model-pdb` 省略時に使用。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` を 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer / --no-detect-layer` | 入力 PDB の B 因子（B=0/10/20）から ML/MM レイヤーを検出。 | `True` |
| `-q, --charge INT` | ML 領域の電荷。 | 必須 |
| `-m, --multiplicity INT` | ML 領域のスピン多重度 (2S+1)。 | `1` |
| `--freeze-atoms TEXT` | 凍結する 1 始まりカンマ区切りインデックス（例: `"1,3,5"`）。YAML `geom.freeze_atoms` とマージ。 | _None_ |
| `--func-basis TEXT` | 汎関数/基底関数ペア（`"FUNC/BASIS"`）。 | `wb97m-v/def2-tzvpd` |
| `--max-cycle INT` | 最大 SCF 反復数。 | `100` |
| `--conv-tol FLOAT` | SCF 収束閾値 (Hartree)。 | `1e-9` |
| `--grid-level INT` | DFT 積分グリッドレベル (0=粗, 3=デフォルト, 5=fine, 9=very fine)。 | `3` |
| `-o, --out-dir DIR` | 出力ディレクトリ。 | `./result_dft/` |
| `--config FILE` | 明示的な CLI オプション適用前に読み込むベース YAML。 | _None_ |
| `--show-config/--no-show-config` | 解決済み設定を表示して実行を継続。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに設定検証と実行計画表示のみ行う。`--help-advanced` に表示。 | `False` |
| `--ref-pdb FILE` | XYZ/GJF 入力時の参照 PDB（原子順序と残基マッピングのテンプレート）。 | _None_ |
| `--convert-files/--no-convert-files` | PDB テンプレートがあれば XYZ/TRJ → PDB コンパニオンファイルを生成。 | `True` |

## 出力

```
out_dir/ (デフォルト: ./result_dft/)
├── ml_region_with_linkH.xyz    # DFT に使用された ML 領域座標（リンク水素付き）
├── result.yaml                 # DFT + ML(dft)/MM エネルギーサマリー、電荷、スピン密度
└── (stdout)                    # 整形された設定ブロックとエネルギーの出力
```

- `result.yaml` の内容:
  - `energy`: Hartree/kcal/mol 値、収束フラグ、実行時間、バックエンド情報（gpu4pyscf vs pyscf(cpu)）。
  - `charges`: Mulliken、meta-Lowdin、IAO 原子電荷（手法失敗時は `null`）。
  - `spin_densities`: Mulliken、meta-Lowdin、IAO スピン密度（UKS のみ）。
- 電荷、多重度、スピン (2S)、汎関数、基底関数、収束パラメータ、解決済み出力ディレクトリも要約されます。

## YAML 設定

マッピングルートを受け付けます。`dft` セクション（およびオプションの `geom`、`calc`/`mlmm`）が存在する場合に適用されます。マージ順:
- デフォルト
- `--config`
- 明示的に指定した CLI オプション

`dft` キー（括弧内はデフォルト）:
- `func_basis`（`"wb97m-v/def2-tzvpd"`）: 結合 `FUNC/BASIS` 文字列。
- `conv_tol`（`1e-9`）: SCF 収束閾値 (Hartree)。
- `max_cycle`（`100`）: 最大 SCF 反復数。
- `grid_level`（`3`）: PySCF `grids.level`。
- `verbose`（`4`）: PySCF 冗長度 (0-9)。
- `out_dir`（`"./result_dft/"`）: 出力ディレクトリルート。

```yaml
geom:
 coord_type: cart                  # オプションの geom_loader 設定
calc:
 charge: 0                         # ML 領域の電荷
 spin: 1                           # スピン多重度 2S+1
mlmm:
 real_parm7: real.parm7            # Amber parm7 トポロジー
 model_pdb: ml_region.pdb          # ML 領域定義
dft:
 func_basis: wb97m-v/def2-tzvpd      # 交換相関汎関数 / 基底関数セット
 conv_tol: 1.0e-09                # SCF 収束閾値 (Hartree)
 max_cycle: 100                    # 最大 SCF 反復数
 grid_level: 3                     # PySCF グリッドレベル
 verbose: 4                        # PySCF 冗長度 (0-9)
 out_dir: ./result_dft/            # 出力ディレクトリルート
```

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細なトラブルシューティングガイド

- [freq](freq.md) -- 振動解析（DFT 精密化の前に実行する場合が多い）
- [opt](opt.md) -- 単一構造の構造最適化
- [all](all.md) -- `--dft` 付き一気通貫ワークフロー
- [YAML リファレンス](yaml_reference.md) -- `dft` の完全な設定オプション
- [用語集](glossary.md) -- DFT、SP（一点計算）の定義
