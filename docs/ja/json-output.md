# JSON 出力リファレンス

mlmm は、AI エージェント・スクリプト・下流ツールがプログラムで利用するための機械可読 JSON 出力を提供します。

## `--out-json` フラグ

主要な MLIP 系サブコマンド（`opt`, `sp`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `dft`, `extract`）が `--out-json / --no-out-json`（デフォルト: off）に対応しています。
有効にすると、出力ディレクトリに `result.json` が生成されます。

```bash
mlmm opt -i r_complex_layered.pdb --max-cycles 5 --out-json --out-dir result_opt
cat result_opt/result.json | python -m json.tool
```

`all` / `path-search` は常に `summary.json` を出力します（`--out-json` 不要）。

### `summary.json` ミラー

`write_result_json` は各ステージの `result.json` ペイロードを同じディレクトリの `summary.json` にミラーします。MCP クライアントやエージェントスクリプトは全サブコマンドで単一のファイル名（`summary.json`）を読めば済みます。隣に書かれる `result.json` は同一のペイロードを保持します。

## 共通エンベロープ

すべての `result.json`（およびミラーされた `summary.json`）に自動付与されるフィールド:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `schema_version` | string | エンベロープのスキーマバージョン。現在値は `mlmm.core.utils.RESULT_JSON_SCHEMA_VERSION` に由来する（この文書のリテラルではなく定数を参照すること）。値の更新は構造変更を示す。 |
| `command` | string | サブコマンド名（例: `"opt"`） |
| `mlmm_version` | string | パッケージバージョン |
| `status` | string | 各サブコマンドが設定（自動付与ではない）。all/path-search は success/partial/failed、エラー時は error、ステージ別は converged/not_converged/completed 等 |
| `elapsed_seconds` | float | 実行時間（秒） |
| `environment` | object | ハードウェア情報（下表参照） |

### エラーエンベロープ（`status == "error"` のとき）

| フィールド | 型 | 説明 |
|-----------|------|------|
| `error` | string | 元の例外の `str(exc)` |
| `error_type` | string | 例外クラス名（例: `"OptimizationError"`） |
| `error_class_chain` | list[string] | 完全な MRO クラス名（例: `["OptimizationError", "RuntimeError", "Exception", "BaseException"]`）。テキスト解析なしで階層をマッチできる。 |
| `error_module` | string | 例外クラスが定義されたモジュール |
| `error_label` | string | 高レベルの CLI ステージラベル（例: `"opt"`、`"tsopt-stage"`） |

**`environment`**:

| フィールド | 型 | 例 |
|-----------|------|------|
| `device` | string | `"cuda"` または `"cpu"` |
| `gpu_name` | string | `"NVIDIA GeForce RTX 5080"` |
| `gpu_vram_gb` | float | `16.6` |
| `cuda_version` | string | `"12.9"` |
| `cpu` | string | `"AMD Ryzen 9 7950X 16-Core Processor"` |
| `n_cpus` | int | `32` |
| `ram_gb` | float | `133.7` |

## サブコマンド別スキーマ

### `opt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"converged"` / `"not_converged"` |
| `energy_hartree` | float | 最終 ONIOM エネルギー (Hartree) |
| `n_opt_cycles` | int | 最適化サイクル数 |
| `opt_mode` | string | `"grad"`, `"hess"`, `"light"`, `"heavy"`, `"lbfgs"`, `"rfo"` のいずれか（`light`/`lbfgs` は `grad`、`heavy`/`rfo` は `hess` の別名） |
| `backend` | string | ML バックエンド |
| `charge` | int | モデル領域電荷 |
| `spin` | int | モデル領域スピン多重度 |
| `n_atoms` | int | 全原子数（全レイヤー） |
| `n_freeze_atoms` | int | 凍結原子数 |
| `thresh` | string | 収束閾値プリセット名 |
| `max_cycles` | int | 最大サイクル数 |
| `input_file` | string | 入力ファイル名 |
| `final_max_force` | float | 最終 max gradient (Hartree/Bohr) |
| `final_rms_force` | float | 最終 RMS gradient |
| `final_max_step` | float | 最終 max 変位 (Bohr) |
| `final_rms_step` | float | 最終 RMS 変位 |
| `convergence_thresholds` | object | 収束閾値の数値 |
| `files` | object | 出力ファイルマップ |

### `tsopt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"completed"`（hess モード）/ `"converged"` / `"not_converged"`（grad モードは収束判定を返す） |
| `energy_hartree` | float | TS エネルギー (Hartree) |
| `n_imaginary_modes` | int | 虚振動数 |
| `imaginary_frequencies_cm` | float[] | 虚振動数 (cm$^{-1}$, 負の値) |
| `opt_mode` | string | `"grad"`, `"hess"`, `"light"`, `"heavy"`, `"dimer"`, `"rsirfo"`, `"trim"`, `"rsprfo"` のいずれか（`light`/`dimer` は `grad` (PHG-Dimer)、`heavy`/`rsirfo` は `hess` (RS-I-RFO)、`trim` は TRIM、`rsprfo` は RS-P-RFO の別名） |
| `n_atoms` | int | 全原子数 |
| `n_opt_cycles` | int | 最適化サイクル数 |
| `backend` | string | ML バックエンド |
| `files` | object | 最終構造 + vib モードファイル |

### `freq`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"completed"` |
| `n_modes` | int | 全基準振動数 |
| `n_imaginary` | int | 虚振動数 |
| `frequencies_cm` | float[] | 全振動数 (cm$^{-1}$) |
| `imaginary_frequencies_cm` | float[] | 負の振動数のみ |
| `thermochemistry` | object\|null | 熱化学データ |
| `backend` | string | ML バックエンド |
| `n_atoms` | int | 原子数 |
| `files` | object | `{"frequencies_txt": "frequencies_cm-1.txt"}` |

**`thermochemistry`** (thermoanalysis 利用不可時は null):

`temperature_K`, `pressure_atm`, `zpe_ha`, `thermal_correction_energy_ha`, `thermal_correction_enthalpy_ha`, `thermal_correction_free_energy_ha`, `sum_EE_and_ZPE_ha`, `sum_EE_and_thermal_energy_ha`, `sum_EE_and_thermal_free_energy_ha`, `E_thermal_cal_per_mol`, `Cv_cal_per_mol_K`, `S_cal_per_mol_K`

### `irc`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_frames_forward` / `backward` / `total` | int | IRC フレーム数 |
| `energy_reactant_hartree` | float | 反応物エネルギー |
| `energy_ts_hartree` | float | TS エネルギー |
| `energy_product_hartree` | float | 生成物エネルギー |
| `forward_converged` / `backward_converged` | bool | IRC 収束判定 |
| `backend` | string | ML バックエンド |
| `bond_changes` | object | `{formed: [...], broken: [...]}` |

### `scan` / `scan2d` / `scan3d`

scan は `stages[]` 配列にステージごとのデータと `n_stages` を含みます。scan2d/scan3d は `n_grid_points` と `pair1`/`pair2`(/`pair3`)（各 `{i, j, low, high}`）に加えて、表面の最小エネルギー `min_energy_hartree` を含みます。いずれも `backend`・`charge`・`spin` を備えます。

### `path-opt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `converged` | bool | 収束判定 |
| `mep_mode` | string | `"dmf"` / `"gsm"` |
| `image_energies_hartree` | float[] | 全イメージエネルギー |
| `barrier_kcal` | float | 前方障壁 (kcal/mol) |
| `delta_kcal` | float | 反応エネルギー (kcal/mol) |
| `backend` | string | ML バックエンド |

### `dft`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `converged` | bool | SCF 収束? |
| `energy_hartree` | float | DFT エネルギー |
| `xc_functional` | string | 汎関数 |
| `basis_set` | string | 基底関数 |
| `used_gpu` | bool | GPU 使用? |
| `backend` | string | ONIOM 高レベル領域の ML バックエンド |
| `charges` / `spin_densities` | object | `{mulliken, lowdin, iao}` 原子電荷/スピン密度 |

### `extract`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_atoms_extracted` | int | 抽出後の原子数 |
| `total_charge` | float | 合計電荷 |
| `protein_charge` | float | タンパク質電荷 |
| `ligand_total_charge` | float | リガンド電荷合計 |
| `unknown_residue_charges` | object | `{残基名: 電荷}` |
| `center` | string | 基質指定（生の `-c` 値）: PDB パス、残基ID リスト（例 `'A:123,B:456'`）、または残基名リスト（例 `'GPP,MMT'`） |
| `radius` | float | 抽出半径 (angstrom) |

## `summary.json` (`path-search` / `all`)

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"success"` / `"partial"` / `"failed"`（all。path-search は success/partial） |
| `segments` | object[] | セグメントごとの障壁、反応エネルギー、結合変化 |
| `energy_diagrams` | object[] | エネルギーダイアグラム |
| `mlip_backend` | string | モデル名 |
| `environment` | object | ハードウェア情報 |

`all` はさらに `n_segments_reactive`（bridge 以外の反応セグメント数）, `rate_limiting_step`, `overall_reaction_energy_kcal`, `post_segments` を含みます。

## 使用例

### Python

```python
import json

with open("result_opt/result.json") as f:
    result = json.load(f)

if result["status"] == "converged":
    print(f"Energy: {result['energy_hartree']:.6f} Hartree")
else:
    print(f"Not converged after {result['n_opt_cycles']} cycles")
```

### jq

```bash
jq '.status' result.json                    # 収束確認
jq '.barrier_kcal' result.json               # 障壁エネルギー
jq '.imaginary_frequencies_cm' result.json   # 虚振動数
jq '.thermochemistry.sum_EE_and_thermal_free_energy_ha' result.json  # 自由エネルギー
```
