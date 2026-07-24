# JSON 出力リファレンス

mlmm は、AI エージェント・スクリプト・下流ツールがプログラムから利用するための機械可読 JSON 出力を提供します。

## `--out-json` フラグ

主要な MLIP 系・レポート系サブコマンド（`opt`, `sp`, `tsopt`, `freq`,
`irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `dft`, `extract`, `trj2fig`,
`energy-diagram`）が `--out-json / --no-out-json`（デフォルト: off）に
対応しています。有効にすると、正規の `result.json` と、同一内容の互換ミラー
`summary.json` が通常の出力と同じ場所に生成されます。

```bash
mlmm opt -i r_complex_layered.pdb --parm real.parm7 -q 0 -m 1 \
  --max-cycles 5 --out-json --out-dir result_opt
cat result_opt/result.json | python -m json.tool
```

`all` / `path-search` は、集約結果を書き込む段階まで到達すると、`--out-json` なしで `summary.json` を出力します。早期の CLI 引数または入力の検証で失敗した場合は、ファイルが作られないことがあります。

### `summary.json` ミラー

`write_result_json` は両方の名前に同じバイト列を準備し、互換ミラー
`summary.json` を先に、正規の `result.json` を最後に公開します。
正常終了時は両ファイルのバイト列が同一です。公開が中断された場合は
`result.json` を正規とし、処理と書き込みが正常終了したことを確認してください。
実行管理側が `run_id` を割り当てた場合は、その一致も検証します。

## 共通エンベロープ

共通の書き込み処理と集約結果の生成側が供給するフィールドを示します。
任意と記したフィールドは、生成側が対応するデータを渡した場合にだけ含まれます。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `schema_version` | string | エンベロープのスキーマバージョン。現在値は `mlmm.core.utils.RESULT_JSON_SCHEMA_VERSION` に由来する（この文書のリテラルではなく定数を参照すること）。値の更新は構造変更を示す。 |
| `command` | string | leaf envelope はサブコマンド名（例: `"opt"`）、aggregate `all` / `path-search` summary は完全な invocation string。 |
| `mlmm_version` / `mlmm_toolkit_version` | string | パッケージバージョン（leaf は `mlmm_version`、aggregate summary は `mlmm_toolkit_version`）。 |
| `status` | string | コマンド固有。`all` は success/partial/failed、`path-search` は success/partial、opt は converged/not_converged/stalled、tsopt はさらに unverified、完了した解析/積分 stage は completed、例外 envelope は error。 |
| `elapsed_seconds` | float | 任意の実行時間（秒）。shared writer に時間を渡さない producer では省略。 |
| `environment` | object | ハードウェア情報（下表参照） |
| `run_id` | string | 任意。MCP などの orchestrator が現在の呼び出し identity を割り当てた場合に含まれる。矛盾する caller 値は拒否される。 |

MLIP/ML/MM calculator stageでは、さらに以下を記録します:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `mlip_backend` | string \| null | backend識別子（`uma`, `orb`, `mace`, `aimnet2`, `custom`）。plot-only commandがcalculatorを評価していない場合はnull |
| `mlip_model` | string \| null | 正確なmodel/checkpoint。`--calc-file`では`filename:factory` |
| `mlip_precision` | string \| null | 実効精度（`fp32` / `fp64`）。custom calculatorではnull |
| `mm_backend` | string \| null | MM energy/Hessian backend（`hessian_ff` / `openmm`）。plot-onlyではnull |
| `link_atom_method` | string \| null | link atom配置（`scaled` / `fixed`）。plot-onlyではnull |
| `use_cmap` | bool \| null | CMAP項を有効にしたか。plot-onlyではnull |

### 実行結果と科学的妥当性

複数段階のワークフローと scan の出力処理は、構成要素を評価できる場合に以下のフィールドを追加します。出力されるフィールドはコマンドによって異なり、各コマンド固有の `status` も互換性のため維持されます。科学的に利用できるかを判断する際は、`scientific_status` と各 outcome を確認してください。収束を確認できない個別結果は安全側に倒して扱われ、`usable` にはなりません。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `execution_status` | string | 通常は `completed` または `failed`。必須の構成コマンドが実行されたかを示します。 |
| `scientific_status` | string | `success`、`partial`、`failed`。得られた科学的結果が完全かつ利用可能かを示します。 |
| `scientific_status_reasons` | string[] | 利用できない、または欠落した個別結果の理由。正常終了時は省略されます。集約ワークフローの従来の `status_reasons` とは別です。 |
| `expected_item_ids` / `observed_item_ids` | string[] | 集約結果の欠落を検出するための、期待された項目と観測された項目の ID。 |
| `stage_outcomes` | object[] | `stage`、`item_id`、`required`、`executed`、`converged`、`usable`、`reason`、`artifacts` を持つ段階別 outcome。 |
| `point_outcomes` | object[] | `point_id`、`executed`、`converged`、`energy_valid`、`artifact_written`、`seed_eligible`、`reason` を持つ scan 点別 outcome。 |

`run_id` が存在する場合は、現在の呼び出しを識別します。`all` の集約結果では
`current_output_paths` と `key_output_files` をその呼び出しの manifest から
再構築するため、再利用した出力ディレクトリに残る既存ファイルは除外されます。

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
| `gpu_name` | string | `"<gpu model>"` |
| `gpu_vram_gb` | float | `<vram in GB>` |
| `cuda_version` | string | `"<cuda version>"` |
| `cpu` | string | `"<cpu model>"` |
| `n_cpus` | int | `<int>` |
| `ram_gb` | float | `<ram in GB>` |

オプティマイザは `"status": "stalled"` を返すこともあります。これは、設定した force/step の収束基準を満たさないまま、設定ウィンドウにわたってエネルギーが減少しなくなった状態（エネルギープラトー）です。stalled は converged とは別の非収束アウトカムであり、`converged` として報告されることは決してありません。停滞した最適化を繰り返さないよう、以降の flatten/再試行も停止します。存在する場合は `stop_reason` にエネルギー範囲・ウィンドウ・満たせなかった基準が記録されます。stalled は（例えば摂動した構造やより厳しいステップ制御で）再試行し得るものであり、`max_cycles` 枯渇や一般的な失敗のエイリアスではありません。microiteration では、macro ステップの stall と直近の micro（MM）緩和の stall はいずれも真実に報告され、macro 収束として偽装されることはありません。

## サブコマンド別スキーマ

### `sp`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` / `stage` | string / string | `"ok"` / `"sp"` |
| `input` | string | 入力構造path |
| `real_parm7` | string | 全系Amber topology path |
| `charge` / `spin` | int / int | model領域の電荷とspin多重度 |
| `energy_au` | float | ONIOM一点energy (Hartree) |
| `forces_path` | string | `forces.npy`のpath |
| `hessian_path` | string \| null | `hessian.npy`のpath。`--hess`無指定時はnull |
| `elapsed` | string | 人間可読の経過時間 |

### `opt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"converged"` / `"not_converged"` / `"stalled"`（エネルギープラトー、上記参照） |
| `stop_reason` | string | 非収束停止（stalled/stopped）時のみ出力。エネルギープラトーの範囲・ウィンドウと満たせなかった基準を記録 |
| `energy_hartree` | float | 最終 ONIOM エネルギー (Hartree) |
| `n_opt_cycles` | int | 最適化サイクル数 |
| `opt_mode` | string | `"grad"`, `"hess"`, `"light"`, `"heavy"`, `"lbfgs"`, `"rfo"` のいずれか（`light`/`lbfgs` は `grad`、`heavy`/`rfo` は `hess` の別名） |
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
| `rigid_projection` | object\|null | `--flatten` が PHVA を実行した場合の凍結境界 TR provenance |
| `files` | object | 出力ファイルマップ |

### `tsopt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | optimizer 収束かつ `n_imaginary_modes == 1` の場合だけ `"converged"`。それ以外は `"not_converged"`、`--skip-final-freq` 時は `"unverified"`。エネルギープラトーによる `"stalled"`（上記参照）はこれらすべてに優先し、`converged` として報告されることはありません（dimer (grad) モードも `stalled` を返します）。 |
| `energy_hartree` | float | TS エネルギー (Hartree) |
| `n_imaginary_modes` | int | 虚振動数 |
| `imaginary_frequencies_cm` | float[] | 虚振動数 (cm$^{-1}$, 負の値) |
| `opt_mode` | string | `"grad"`, `"hess"`, `"light"`, `"heavy"`, `"dimer"`, `"rsirfo"`, `"trim"`, `"rsprfo"` のいずれか（`light`/`dimer` は `grad` (PHG-Dimer)、`heavy`/`rsirfo` は `hess` (RS-I-RFO)、`trim` は TRIM、`rsprfo` は RS-P-RFO の別名） |
| `n_atoms` | int | 全原子数 |
| `n_opt_cycles` | int | 最適化サイクル数 |
| `rigid_projection` | object | Dimer/flatten/最終鞍点解析の凍結境界 TR provenance |
| `reference_mode_file` | string\|null | `--ref-mode`で渡した高度なpath由来mode |
| `safeguards` | object | heavy modeのtrial拒否/recovery、exact saddle、target-mode診断 |
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
| `n_atoms` | int | 原子数 |
| `rigid_projection` | object | 振動解析と熱化学で使った凍結境界 TR provenance |
| `files` | object | 出力map。`--dump-hess`時は`hessian_npz`を含む |

**`thermochemistry`** (thermoanalysis 利用不可時は null):

`temperature_K`, `pressure_atm`, `symmetry_number`, `symmetry_number_source`, `zpe_ha`, `thermal_correction_energy_ha`, `thermal_correction_enthalpy_ha`, `thermal_correction_free_energy_ha`, `sum_EE_and_ZPE_ha`, `sum_EE_and_thermal_energy_ha`, `sum_EE_and_thermal_free_energy_ha`, `E_thermal_cal_per_mol`, `Cv_cal_per_mol_K`, `S_cal_per_mol_K`

### `irc`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_frames_forward` / `n_frames_backward` / `n_frames_total` | int | IRC フレーム数 |
| `energy_first_hartree` | float | 連結経路の最初の端点。単独 IRC は反応物/生成物の化学的な同一性を割り当てない |
| `energy_ts_hartree` | float | TS エネルギー |
| `energy_last_hartree` | float | 連結経路の最後の端点。単独 IRC は反応物/生成物の化学的な同一性を割り当てない |
| `endpoint_energy_orientation` | string | `"finished_first_to_finished_last"` |
| `energy_reactant_hartree` / `energy_product_hartree` | float | 最初/最後の端点を表す互換エイリアス。キー名から反応物/生成物の同一性を推定しない |
| `forward_converged` / `backward_converged` | bool\|null | 各方向の収束フラグ |
| `never_stop` | bool | 任意指定のエネルギー上昇・平坦化回避モードを有効にしたか |
| `never_stop_energy_bypasses` | int | 実際に回避したエネルギー上昇・平坦化停止イベントの数 |
| `rigid_projection` | object | 初期/更新Hessianの凍結境界 TR provenance |
| `rigid_projection.electronic_state_verified` | bool | ファイルから初期化したHessianの model charge・多重度を identity 検証できたか |
| `bond_changes` | object | 最初→最後の方向の `{formed: [...], broken: [...]}`。比較できない場合は省略 |
| `bond_changes_direction` | string | 結合変化がある場合は `"finished_first_to_finished_last"` |
| `files` | object | 軌跡と端点ファイル（XYZと、利用可能なPDB/CIF companion） |

**`rigid_projection` provenance:** 選択した処理は `treatment`、有効 rank は
`effective_rank` として、アクティブ/凍結原子数・インデックス、および各 workflow が
使ったHessian source/shape とともに記録します。デフォルトは `constrained`。
`legacy-active` は非推奨の比較専用で、pass/HOSP 遷移状態認定には使用できません。`freq --dump` は同じ object を
`thermoanalysis.yaml` にも書き出します。最後の 2 値のキー名は生成 workflow により
`hessian_source` / `hessian_shape` または `source` / `raw_hessian_shape` です。

### `scan` / `scan2d` / `scan3d`

scan は `stages[]` 配列にステージごとのデータと `n_stages` を含みます。各 stage には（追加フィールド）`optimizer_status`（`converged`/`not_converged`/`stalled`）と、そのステージの最後の optimizer が非収束停止した場合の `stop_reason` を含みます。scan2d/scan3d は `n_grid_points` と `pair1`/`pair2`(/`pair3`)（各 `{i, j, low, high}`）に加えて、表面の最小エネルギー `min_energy_hartree` を含みます。fresh runでは共通calculator provenanceと`charge`・`spin`を記録します。plot-onlyの`scan3d --csv`も同じキーを保持しますが、importしたenergy gridからcalculatorを特定できないため、`mlip_backend`、`mlip_model`、`mlip_precision`、`mm_backend`、`link_atom_method`、`use_cmap`、`charge`、`spin`はnullです。

### `path-opt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `converged` | bool | 収束判定 |
| `mep_mode` | string | `"dmf"` / `"gsm"` |
| `image_energies_hartree` | float[] | 全イメージエネルギー |
| `barrier_kcal` | float | 前方障壁 (kcal/mol) |
| `delta_kcal` | float | 反応エネルギー (kcal/mol) |

### `dft`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `converged` | bool | SCF 収束? |
| `status` | string | `"converged"` または `"not_converged"`。後者は exit code 3 より前に書き込まれる。 |
| `energy_hartree` | float | DFT エネルギー |
| `xc_functional` | string | 汎関数 |
| `basis_set` | string | 基底関数 |
| `used_gpu` | bool | GPU 使用? |
| `grid_level` | int | YAML/CLI 解決後の DFT grid level |
| `conv_tol` | float | YAML/CLI 解決後の SCF 収束閾値 |
| `max_cycle` | int | YAML/CLI 解決後の SCF 最大反復数 |
| `engine` | string | 実際に使用した runtime engine label |
| `charges` / `spin_densities` | object | `{mulliken, lowdin, iao}` 原子電荷/スピン密度 |

### `trj2fig`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"ok"` |
| `n_frames` | int | 軌跡のフレーム数。 |
| `min_energy_hartree` / `max_energy_hartree` | float | フレームエネルギーの最小値と最大値。 |
| `energy_source` | string | `"trajectory_comment"` または `"mlip_recomputed"`。 |
| `mlip_backend` / `mlip_model` / `mlip_precision` | string \| null | 再計算で確定した来歴。コメントモードではすべて null。 |
| `charge` / `multiplicity` | int \| null | 再計算で解決された電荷とスピン多重度。コメントモードでは null。再計算時に省略した値は 0 と 1。 |
| `output_files` | string[] | すべての出力パスを順序どおりに保持する正規フィールド。別ディレクトリに同名ファイルがあっても保持される。 |
| `files` | object | 後方互換用のベース名からパスへの対応表。同じベース名が重複すると一方だけが残る。 |

`-q/--charge` または `-m/--multiplicity` のいずれかを指定すると、選択した MLIP で全フレームを再計算します。これは MLIP による各フレームの直接再評価であり、トポロジーやモデル領域の入力を受け取らず、ONIOM エネルギーは計算しません。

### `extract`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_atoms_extracted` | int | 抽出後の原子数 |
| `total_charge` | float | 合計電荷 |
| `protein_charge` | float | タンパク質電荷 |
| `ligand_total_charge` | float | リガンド電荷合計 |
| `unknown_residue_charges` | object | `{残基名: 電荷}` |
| `center` | string | 基質指定（生の `-c` 値）: PDB パス、残基ID リスト（例 `'A:123,B:456'`）、または残基名リスト（例 `'GPP,MMT'`） |
| `radius` | float | 抽出半径 (Å) |
| `status` | string | `"ok"` |
| `ion_total_charge` | float | イオン電荷合計 |
| `input_files` | string[] | 入力 PDB パス |
| `n_atoms_raw` | int | 抽出前の生入力の原子数 |
| `n_link_hydrogens` | int | 切断結合に付加されたリンク H 原子数 |
| `files` | object | 出力ファイル名のマップ（入力ごとのポケット PDB 等） |
| `exclude_backbone` | bool | 実行時の `--exclude-backbone` の値 |
| `include_h2o` | bool | 実行時の `--include-H2O` の値 |
| `ligand_charge_input` | string | 生の `-l/--ligand-charge` 引数 |
| `ion_charges` | array | イオン残基の `[残基名, 電荷]` ペアのリスト |

### `energy-diagram`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"ok"` |
| `n_points` | int | エネルギーデータ点の数 |
| `files` | object | 出力ダイアグラムのファイル名からパスへの対応表 |

## `summary.json` (`path-search` / `all`)

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"success"` / `"partial"` / `"failed"`（all。path-search は success/partial） |
| `execution_status` / `scientific_status` | string / string | 実行の完了度と科学的な利用可能性。従来の `status` とは分けて評価します。 |
| `scientific_status_reasons` | string[] | 不完全または利用できない科学的結果の理由。正常終了時は省略されます。 |
| `expected_item_ids` / `observed_item_ids` | string[] | 期待された集約項目と観測された集約項目。 |
| `n_segments` | int | セグメント数 |
| `segments` | object[] | セグメントごとの障壁、反応エネルギー、結合変化 |
| `energy_diagrams` | object[] | エネルギーダイアグラム |
| `mlip_backend` | string | バックエンド名（`uma`, `orb`, `mace`, `aimnet2`, `custom`） |
| `mlip_model` | string \| null | バックエンドと分離して記録するモデル/checkpoint名 |
| `mlip_precision` | string \| null | 実効`fp32` / `fp64`。custom calculatorではnull |
| `charge` | int | モデル領域の電荷 |
| `spin` | int | モデル領域のスピン多重度 |
| `environment` | object | ハードウェア情報 |

`all` はさらに以下を含みます。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_segments_reactive` | int | bridge 以外の反応セグメント数。 |
| `rate_limiting_step` | object | 互換性のため維持するキー。各段階の始状態を基準にした局所障壁が最大のセグメントと method。microkinetics に基づく律速段階の判定ではない。 |
| `overall_reaction_energy_kcal` | float | 全体の反応エネルギー。 |
| `post_segments` | list | セグメントごとの TS/IRC/freq/DFT 結果。 |
| `post_segments[].thermo_symmetry` | object | 子 freq が報告した状態別の回転対称 provenance。有効な `symmetry_number` と `symmetry_number_source` の両方を持つ R/TS/P 状態だけを含み、欠けた状態は省略する。どの状態にも有効な provenance が無い場合だけフィールド全体を省略する。 |
| `key_output_files` | object | 現在の呼び出しの出力索引。ルートファイルはファイル名 → 説明、各 `seg_NN` は `{description, files}` で、`files` はそのセグメントディレクトリからの相対パス。 |
| `current_output_paths` | string[] | `--out-dir` からの相対パスを並べたリスト。現在の呼び出しが記録した成果物だけを含みます。 |

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
