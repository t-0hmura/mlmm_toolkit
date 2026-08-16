# `dft`

GPU4PySCF（または CPU PySCF）を使用して ML 領域で DFT 一点エネルギー計算を実行し、QM 領域（ML 領域）の DFT エネルギーを MM エネルギーと合成して ML(dft)/MM 総エネルギーを取得します。DFT 勾配と力は要求しません。`mlmm dft` は酵素全体の PDB から ML 領域を抽出し、リンク水素を付加したうえで PySCF（または GPU4PySCF）で計算します。MLIP 経路探索後の停留点（R / TS / P / IM）に対する DFT 一点エネルギー評価や、MLIP 障壁の基準汎関数/基底による sanity check に使用します。デフォルトの汎関数/基底関数は `wb97m-v/def2-tzvpd` です。結果にはエネルギーと集団解析（Mulliken、meta-Lowdin、IAO 電荷）が含まれます。

```
E_total = E_REAL_low + E_ML(DFT) - E_MODEL_low
```

## 実行例

ML 領域に対する最小構成の DFT 一点計算:

```bash
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_dft
```

```bash
# 汎関数/基底関数を変更して一点計算する
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --func-basis "wb97m-v/def2-tzvpd" --out-dir ./result_dft_tz
```

```bash
# SCF 収束を厳しくして反復回数を増やす
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --conv-tol 1e-10 --max-cycle 200 --out-dir ./result_dft_tight
```

## 処理の流れ

1. **入力処理** -- ML/MM 中核の `MLMMCore` が酵素全体の PDB（`-i`）、Amber トポロジー（`--parm`）、ML 領域定義（`--model-pdb` または `--model-indices` または `--detect-layer` による B 因子検出）を読み込みます。YAML で明示的な `link_mlmm` ペアが指定されない限り、ML/MM 選択を横切る parm7 結合にリンク水素を自動付加します。結合の認識に距離は使用しません。
2. **SCF 構築** -- `--func-basis` でスラッシュ区切りで汎関数/基底関数を定義します。GPU4PySCF バックエンドは利用可能な場合に使用され、closed-shell の GPU 経路では `--lowmem`（デフォルト）が有効なら低メモリ実装 `gpu4pyscf.dft.rks_lowmem.RKS` を使用します。CPU モードを強制するには `--engine cpu` を使用してください。`mlmm dft` は SCF オブジェクトに対して `density_fit()` を呼びません。標準 GPU/CPU 経路ではバックエンドのデフォルト JK 実装を使用し、lowmem 経路では `rks_lowmem.RKS` のメモリ効率の良い直接 JK が使用されます。
3. **ML(dft)/MM 再結合** -- DFT は `MLMMCore` の高レベル MODEL エネルギーだけを置き換えます。`MLMMCore` は選択した MM バックエンドで REAL-low と MODEL-low を評価し、差し引き式を適用します。別個のトポロジー構築、MM calculator 経路、DFT 力計算はありません。
4. **集団解析と出力** -- Mulliken、meta-Lowdin、IAO 電荷とスピン密度（UKS のみ）が結合エネルギーブロックとともに `result.yaml` に書き出されます。

## 出力

```
out_dir/ (デフォルト: ./result_dft/)
├── ml_region_without_linkH.xyz # リンク水素生成前の正確な ML 選択
├── ml_region_with_linkH.xyz    # リンク水素生成後の PySCF 入力 snapshot
├── ml_region_without_linkH.pdb # PDB 入力かつ --convert-files 時。topology 付き companion
├── ml_region_with_linkH.pdb    # PDB 入力かつ --convert-files 時。生成リンク水素は HL/LKH
├── result.yaml                 # DFT + ML(dft)/MM エネルギーサマリー、電荷、スピン密度
├── result.json                 # --out-json 指定時のみ
└── (stdout)                    # 整形された設定ブロックとエネルギーの出力
```

- `result.yaml` の内容:
  - `energy`: Hartree/kcal/mol 値、収束フラグ、実行時間、バックエンド情報（`engine`: `gpu4pyscf(rks_lowmem)` / `gpu4pyscf` / `pyscf(cpu)`、`used_gpu`、`used_lowmem`）。
  - `mlmm_energy`: REAL-low / MODEL-low の MM 評価値と再結合エネルギー `E_total = E_REAL_low + E_ML(DFT) - E_MODEL_low`（Hartree と kcal/mol）。
  - `charges [index, element, mulliken, lowdin, iao]`: Mulliken / meta-Lowdin / IAO 原子電荷（計算に失敗した場合は `null`）。
  - `spin_densities [index, element, mulliken, lowdin, iao]`: 同形式のスピン密度（UKS のみ）。
- 電荷、多重度、スピン (2S)、汎関数、基底関数、収束パラメータ、解決済み出力ディレクトリも要約されます。

## CLI オプション

`mlmm dft --help` でコアオプション、`mlmm dft --help-advanced` で全オプションを表示します。全フラグの一覧は生成された[コマンドリファレンス](../reference/commands/index.md)にあります。以下の表では説明が必要なオプションを取り上げます。

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `--cmap/--no-cmap` | REAL と MODEL の両 MM 層で CMAP を保持します。 | `--cmap` |
| `-i, --input PATH` | 酵素全体の構造ファイル（PDB/mmCIF、または `--ref-pdb` でトポロジーを指定した XYZ）。 | 必須 |
| `--parm PATH` | 全系の Amber parm7 トポロジー。 | 必須 |
| `--model-pdb PATH` | ML 領域を定義する PDB（原子 ID が酵素 PDB と一致必須）。`--detect-layer` 有効時はオプション。 | _None_ |
| `--model-indices TEXT` | ML 領域のカンマ区切り原子インデックス（範囲指定可、例: `1-5`）。`--model-pdb` 省略時に使用。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` を 1 始まりまたは 0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer` | 入力 PDB の B 因子（B=0/10/20）から ML/MM レイヤーを自動検出。 | 有効 |
| `-q, --charge INT` | ML 領域の電荷。`-l/--ligand-charge` 指定時は不要（PDB 入力または `--ref-pdb` 付き XYZ）。 | `-l` 指定時を除き必須 |
| `-l, --ligand-charge TEXT` | 全体電荷、または残基名ごとのマッピング（例: `SAM:1,GPP:-3`）。`-q` 省略時に ML 領域の電荷を導出するために使用（PDB 入力または `--ref-pdb` が必要）。 | _None_ |
| `-m, --multiplicity INT` | ML 領域のスピン多重度 (2S+1)。 | `1` |
| `--func-basis TEXT` | 汎関数/基底関数ペア（`"FUNC/BASIS"`）。 | `wb97m-v/def2-tzvpd` |
| `--max-cycle INT` | 最大 SCF 反復数。 | `100` |
| `--conv-tol FLOAT` | SCF 収束閾値 (Hartree)。 | `1e-9` |
| `--grid-level INT` | DFT 積分グリッドレベル (0=粗, 3=デフォルト, 5=細かい, 9=非常に細かい)。 | `3` |
| `--engine {gpu,cpu}` | GPU4PySCF（`gpu`）または CPU PySCF（`cpu`）を強制。 | `gpu` |
| `--lowmem/--no-lowmem` | closed-shell の GPU 経路で `gpu4pyscf.dft.rks_lowmem.RKS` を使用（メモリ効率の良い直接 JK；`mlmm dft` はどちらの経路でも `density_fit()` を呼ばない）。open-shell、CPU、`rks_lowmem` 非搭載の旧 `gpu4pyscf` では標準 RKS/UKS に自動フォールバック。 | `True` |
| `-o, --out-dir DIR` | 出力ディレクトリ。 | `./result_dft/` |
| `--link-atom-method {scaled,fixed}` | リンク原子位置モード: `scaled`（g-factor、Gaussian ONIOM 標準）または `fixed`（旧式 1.09/1.01 Å 固定）。 | `scaled` |
| `--mm-backend {hessian_ff,openmm}` | ONIOM 低レベル評価用 MM バックエンド。Hessian 構築法は `calc.mm_fd` が別に制御します（デフォルト `true`: 有限差分）。 | `hessian_ff` |
| `--out-json/--no-out-json` | 機械可読な `result.json` を `out_dir` に出力。 | `False` |
| `--config FILE` | 明示的な CLI オプション適用前に読み込むベース YAML。 | _None_ |
| `--show-config/--no-show-config` | 解決済み設定を表示して実行を継続。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに設定検証と実行計画表示のみ行う。`--help-advanced` に表示。 | `False` |
| `--ref-pdb FILE` | XYZ 入力時の参照 PDB トポロジー（原子順序と残基マッピングのテンプレート）。 | _None_ |
| `--convert-files/--no-convert-files` | PDB テンプレートがあれば XYZ/TRJ → 対応する PDB ファイルを生成。 | `True` |

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
- `verbose`（`0`）: PySCF verbose レベル (0-9)。デフォルトは quiet。CLI `-v 2/3` では実行時に PySCF verbose レベルが最低 `4` へ上がります。
- `out_dir`（`"./result_dft/"`）: 出力ディレクトリルート。

```yaml
geom:
 coord_type: cart                  # オプションの geom_loader 設定
calc:
 model_charge: 0                   # ML 領域の電荷
 model_mult: 1                     # スピン多重度 2S+1
 real_parm7: real.parm7            # Amber parm7 トポロジー
 model_pdb: ml_region.pdb          # ML 領域定義
dft:
 func_basis: wb97m-v/def2-tzvpd      # 交換相関汎関数 / 基底関数セット
 conv_tol: 1.0e-09                # SCF 収束閾値 (Hartree)
 max_cycle: 100                    # 最大 SCF 反復数
 grid_level: 3                     # PySCF グリッドレベル
 verbose: 0                        # PySCF verbose レベル (0-9); CLI -v 2/3 では実行時 PySCF verbose レベルが >=4
 out_dir: ./result_dft/            # 出力ディレクトリルート
```

## 注記

- 基底関数名が `def2` で始まる場合、対応する def2 有効内殻ポテンシャル（ECP）が自動的に付加されます（元素の有無はチェックしません）。
- **Blackwell アーキテクチャ GPU**（RTX 50xx）: インストールした
  GPU4PySCF/CuPy が対象デバイスをサポートすることを確認してください。
  GPU 経路が失敗する場合は `--engine cpu` または外部 DFT
  プログラムを使用してください。
- **def2-TZVPD でメモリ不足になる場合**: 必要メモリは元素、基底関数、
  汎関数、グリッド、ソフトウェア構成に依存します。対象系で試行し、
  必要なら目的の物理量への影響を検証した上で小さい基底関数を選択してください。
- GPU4PySCF のコンパイル済みホイールは非 x86 環境では動作しない場合があります。ソースからビルドしてください（参照: https://github.com/pyscf/gpu4pyscf）。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細なトラブルシューティングガイド
- [freq](freq.md) -- 振動解析（DFT 一点エネルギー評価の前に実行する場合が多い）
- [opt](opt.md) -- 単一構造の構造最適化
- [all](all.md) -- `--dft` 付き一気通貫ワークフロー
- [YAML リファレンス](yaml-reference.md) -- `dft` の完全な設定オプション
- [用語集](glossary.md) -- DFT、SP（一点計算）の定義
