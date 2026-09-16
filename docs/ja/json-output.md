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
| `status` | string | コマンド固有。`all` は success/partial/failed、`path-search` は success/partial、`opt` と `tsopt` は数値outcomeの converged/not_converged/stalled、完了した解析/積分 stage は completed、例外 envelope は error。TSの鞍点次数は `saddle_validation` / `hessian_status` に分離して記録します。 |
| `elapsed_seconds` | float | 任意の実行時間（秒）。shared writer に時間を渡さない producer では省略。 |
| `environment` | object | ハードウェア情報（下表参照） |
| `run_id` | string | 任意。MCP などの orchestrator が現在の呼び出し identity を割り当てた場合に含まれる。矛盾する caller 値は拒否される。 |

MLIP/ML/MM calculator stageでは、さらに以下を記録します:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `mlip_backend` | string \| null | backend識別子（`uma`, `orb`, `mace`, `aimnet2`, `dft`, `custom`）。DFT leaf は `dft`、plot-only commandがcalculatorを評価していない場合はnull |
| `mlip_model` | string \| null | 正確なmodel/checkpoint。`--calc-file`では`filename:factory` |
| `mlip_precision` | string \| null | 実効精度（`fp32` / `fp64`）。custom calculatorではnull |
| `mm_backend` | string \| null | MM energy/Hessian backend（`hessian_ff` / `openmm`）。plot-onlyではnull |
| `link_atom_method` | string \| null | link atom配置（`scaled` / `fixed`）。plot-onlyではnull |
| `use_cmap` | bool \| null | CMAP項を有効にしたか。plot-onlyではnull |

### 実行と要求段階の完了状況

複数段階のワークフローと scan の出力処理は、構成要素を評価できる場合に以下のフィールドを追加します。出力されるフィールドはコマンドによって異なり、各コマンド固有の `status` も互換性のため維持されます。要求段階の完了状況は `scientific_status` と各 outcome で確認できます。必須の最適化・計算結果が欠ける場合は未完了として記録します。IRC の停止理由・端点 stationary 判定は診断情報です。IRC 独立の `scientific_status` は出力せず、all は TSOPT と両端点 OPT の数値収束を集約します。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `execution_status` | string | 通常は `completed` または `failed`。必須の構成コマンドが実行されたかを示します。 |
| `scientific_status` | string | `success`、`partial`、`failed`。要求した計算段階と最適化・SCF 結果の完了度。振動・結合対応の解釈は別に確認します。 |
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
| `opt_mode` | string | `"grad"`, `"hess"`, `"lbfgs"`, `"rfo"` のいずれか |
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
| `status` | string | 全体の outcome。`optimization_status` と一致するが、収束したのに最終 TS エネルギーを評価できなかった場合は `"energy_missing"` に降格する。数値収束と鞍点次数は `optimization_status` / `saddle_validation` を個別に参照 |
| `optimization_status` | string | 数値 optimizer の結果: `"converged"` / `"not_converged"` / `"stalled"`。鞍点次数とは独立 |
| `saddle_validation` | string | 終端 exact PHVA による `"first_order"` / `"higher_order"` / `"no_imaginary"` / `"unavailable"` |
| `saddle_order_verified` | bool | `saddle_validation: "first_order"` の場合だけ `true` |
| `hessian_status` | string | `"completed"` / `"failed"` / `"skipped"` / `"unavailable"`。失敗理由は `hessian_error` |
| `reaction_mode_index` | int\|null | downstream IRC に使う負の exact-PHVA root。root 0 fallback は明示され、反応 identity を保証しない |
| `reaction_mode_frequency_cm` | float\|null | 選択した負 root の振動数 |
| `reaction_mode_source` | string\|null | 参照方向整合または明示 fallback による root 選択元 |
| `energy_hartree` | float \| null | TS エネルギー (Hartree)。最終エネルギー評価に失敗した場合は `null`（writer が非有限 float をすべて `null` に置換する）で、そのとき `status` は `"energy_missing"` |
| `n_imaginary_modes` | int\|null | 虚振動モードの数。PHVA を実行しなかった場合は `null` |
| `imaginary_frequencies_cm` | float[]\|null | 虚振動数 (cm$^{-1}$, 負の値)。PHVA 未実行時は `null` |
| `opt_mode` | string | `"grad"`, `"hess"`, `"dimer"`, `"rsprfo"`, `"rsirfo"`, `"trim"` のいずれか。`hess` は RS-P-RFO を選択。 |
| `opt_mode_requested` | string | CLI で要求した preset |
| `optimizer` | string | 実際に使用した optimizer algorithm |
| `n_atoms` | int | 全原子数 |
| `n_opt_cycles` | int | 最適化サイクル数 |
| `charge` / `spin` | int / int | model 領域の電荷/多重度 |
| `rigid_projection` | object | Dimer/flatten/最終鞍点解析の凍結境界 TR provenance |
| `reference_mode_file` | string\|null | `--ref-mode` で渡した高度な path 由来 mode。Hessian family のみ |
| `safeguards` | object | Hessian family の trial 拒否/recovery、exact saddle、target-mode 診断 |
| `files` | object | 最終構造 + vib モードファイル |

終端exact PHVAは数値収束後だけ実行します。非収束または`stalled`なら終端構造を
保持してPHVAをskipします。PHVA失敗時は構造を破棄したり振動数を捏造したりせず、`hessian_status: "failed"`
と理由を記録します。数値 status と鞍点次数は独立で、数値収束済み高次停留点は
`optimization_status: "converged"`、`saddle_validation: "higher_order"` のまま
保持され、一次 TS 認定にはなりません。`all` は有効な負 root がある場合だけ警告付き
診断 IRC に進むことがあります。数値非収束、虚振動 0 本、PHVA 失敗/skip、または
有効な負 root なしでは、TS 成果物登録後に IRC 前で停止します。明示的な
`--skip-final-freq` は最終構造を保持し、`n_imaginary_modes: null`、
`imaginary_frequencies_cm: []` を記録します。

### `freq`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"completed"` |
| `n_modes` | int | 基準振動モードの総数 |
| `n_imaginary` | int | 虚振動モードの数 |
| `frequencies_cm` | float[] | 全振動数 (cm$^{-1}$) |
| `imaginary_frequencies_cm` | float[] | 負の振動数のみ |
| `thermochemistry` | object\|null | 熱化学データ |
| `charge` / `spin` | int / int | model 領域の電荷/多重度 |
| `n_atoms` | int | 原子数 |
| `n_freeze_atoms` | int | 凍結原子数 |
| `rigid_projection` | object | 振動解析と熱化学で使った凍結境界 TR provenance |
| `files` | object | 出力map。`--dump-hess`時は`hessian_npz`を含む |

**`thermochemistry`** (thermoanalysis 利用不可時は null):

`temperature_K`, `pressure_atm`, `point_group`, `point_group_source`, `symmetry_number`, `symmetry_number_source`, `electronic_energy_ha`（報告される `E + G_corr = G` の `E`）, `zpe_ha`, `thermal_correction_energy_ha`, `thermal_correction_enthalpy_ha`, `thermal_correction_free_energy_ha`, `sum_EE_and_ZPE_ha`, `sum_EE_and_thermal_energy_ha`, `sum_EE_and_thermal_free_energy_ha`, `E_thermal_cal_per_mol`, `Cv_cal_per_mol_K`, `S_cal_per_mol_K`。`point_group_source` は `auto` または保守的なフォールバックを示す `auto-fallback`、`symmetry_number_source` はこれらに加えて YAML 上書きの `config` / `override` を取ります。

### `irc`

`status: "completed"` は実行が戻ったことを示します。IRC 独自の `scientific_status`、`stage_outcomes`、`forward_status` / `backward_status` は出力しません。方向ごとの停止理由と軌跡を保持し、端点最適化の結果は `all` の `endpoint_opt` に記録します。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `status` | string | `"completed"` |
| `n_frames_forward` / `n_frames_backward` / `n_frames_total` | int | IRC フレーム数 |
| `forward_short_branch` / `backward_short_branch` | bool | cycle 上限前に3フレーム以内で停止した分岐。診断用のみ |
| `energy_first_hartree` | float | 連結経路の最初の端点。単独 IRC は反応物/生成物の化学的な同一性を割り当てない |
| `energy_ts_hartree` | float | TS エネルギー |
| `energy_last_hartree` | float | 連結経路の最後の端点。単独 IRC は反応物/生成物の化学的な同一性を割り当てない |
| `endpoint_energy_orientation` | string | `"finished_first_to_finished_last"` |
| `energy_reactant_hartree` / `energy_product_hartree` | float | 最初/最後の端点を表す互換エイリアス。キー名から反応物/生成物の同一性を推定しない |
| `forward_requested` / `backward_requested` | bool | 各方向を要求したか |
| `forward_integration_converged` / `backward_integration_converged` | bool\|null | RMS 勾配の停留判定が発火して停止したか。診断専用で、`--never-stop` はこの判定を迂回するため常に `false`。削除した `*_converged` が表していた条件は、`*_downhill_departure_valid` との連言で再構成できる |
| `forward_downhill_departure_valid` / `backward_downhill_departure_valid` | bool\|null | TS から downhill に離れたことを確認できたか |
| `forward_integration_stop_reason` / `backward_integration_stop_reason` | string\|null | 数値伝播が失敗した場合だけ非空になる理由 |
| `never_stop` | bool | 任意指定の物理的端点停止回避モードを有効にしたか |
| `never_stop_energy_bypasses` | int | 実際に回避したenergy上昇・1 step energy変化量停止event数 |
| `rigid_projection` | object | 初期/更新 Hessian の凍結境界 TR provenance |
| `rigid_projection.electronic_state_verified` | bool | ファイルから初期化した Hessian の model charge・多重度を identity 検証できたか |
| `bond_changes` | object | 最初→最後の方向の `{formed: [...], broken: [...]}`。比較できない場合は省略 |
| `bond_changes_direction` | string | 結合変化がある場合は `"finished_first_to_finished_last"` |
| `files` | object | 軌跡と端点ファイル（XYZと、利用可能なPDB/CIF companion） |

**`rigid_projection` provenance:** 選択した処理は `treatment`、有効 rank は
`effective_rank` として、アクティブ/凍結原子数・インデックス、および各 workflow が
使った Hessian source/shape とともに記録します。処理は常に `constrained` で、
古い非constrained設定は明示的に拒否されます。`freq --dump` は同じ object を
`thermoanalysis.yaml` にも書き出します。最後の 2 値のキー名は生成 workflow により
`hessian_source` / `hessian_shape` または `source` / `raw_hessian_shape` です。

### `scan` / `scan2d` / `scan3d`

scan は固定の L-BFGS 経路を `scan_opt_mode: "grad"` / `scan_optimizer: "lbfgs"` として記録し、`stages[]` 配列にステージごとのデータと `n_stages` を含みます。各 stage には（追加フィールド）`optimizer_status`（`converged`/`not_converged`/`stalled`）と、そのステージの最後の optimizer が非収束停止した場合の `stop_reason` を含みます。scan2d/scan3d は `n_grid_points` と `pair1`/`pair2`(/`pair3`)（各 `{i, j, low, high}`）に加えて、表面の最小エネルギー `min_energy_hartree` を含みます。fresh run は事前最適化行を除く試行数 `n_points_attempted` と、明示的に収束し有限値・構造 artifact を持つ `n_points_usable`、共通 calculator provenance、`charge`・`spin`を記録します。plot-only の `scan3d --csv` は `n_points_attempted` を出力せず、収束・artifact provenance が完全な CSV の場合だけ `n_points_usable` を出力します。また import した energy grid から calculator を特定できないため、`mlip_backend`、`mlip_model`、`mlip_precision`、`mm_backend`、`link_atom_method`、`use_cmap`、`charge`、`spin` は null です。

### `path-opt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `converged` | bool \| null | 収束判定: エンジン自身の収束シグナルによる `true` / `false`。読み取れない場合は `null`（`status` は `"completed"` となり、収束を主張しない） |
| `mep_mode` | string | `"dmf"` / `"gsm"` |
| `image_energies_hartree` | float[] | 全イメージエネルギー |
| `n_images` | int | イメージ数 |
| `hei_index` | int | 最高エネルギーイメージの index |
| `barrier_kcal` | float | 前方障壁 (kcal/mol) |
| `delta_kcal` | float | 反応エネルギー (kcal/mol) |
| `files` | object | 軌跡と HEI のファイル map |

### `dft`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `converged` | bool | SCF 収束? |
| `status` | string | `"converged"` または `"not_converged"`。後者は exit code 3 より前に書き込まれる。 |
| `energy_hartree` | float | DFT エネルギー |
| `xc_functional` | string | 汎関数 |
| `basis_set` | string | 基底関数 |
| `used_gpu` | bool | GPU 使用? |
| `n_atoms` | int | QM 領域の原子数 |
| `grid_level` | int | YAML/CLI 解決後の DFT grid level |
| `conv_tol` | float | YAML/CLI 解決後の SCF 収束閾値 |
| `max_cycle` | int | YAML/CLI 解決後の SCF 最大反復数 |
| `engine` | string | 実際に使用した runtime engine label |
| `charges` / `spin_densities` | object | `{mulliken, lowdin, iao}` 原子電荷/スピン密度 |
| `files` | object | `{"result_yaml": "result.yaml"}` |

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
| `include_h2o` | bool | 実行時の `--include-h2o` の値 |
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
| `execution_status` / `scientific_status` | string / string | 実行の完了度と、要求した数値最適化・計算段階の完了度。 |
| `scientific_status_reasons` | string[] | 要求した結果の欠損・未収束などの理由。正常終了時は省略されます。 |
| `pipeline_stop` | object \| 不在 | 早期停止時のみ存在。`stage` は `post`（`reason` は `no_segments` / `no_reactive_segment`）、`before_irc`（TSOPT の理由と `segment`・`tsopt_result`）、または `endpoint_opt`（`endpoint_execution_failed` と端点別 `failures`）。`summary.log` では `Pipeline stop` |
| `expected_item_ids` / `observed_item_ids` | string[] | 期待された集約項目と観測された集約項目。 |
| `config` | object | 実効設定。`mep_mode` は GSM/DMF、`ts_opt_mode` / `endpoint_opt_mode` は設定済み後処理 preset を示す。generic `opt_mode*` は解決済み CLI 入力を保持する。`path_opt_mode` は端点 preoptimization に使う単一構造 optimizer であり（`preopt` を参照）、MEP path algorithm ではない。 |
| `n_segments` | int | セグメント数 |
| `search_max_depth` | int | 実効の再帰分割階層上限。`0` は分割無効 |
| `path_optimizers` | string[] | 経路の準備・精密化で実際に使用した単一構造オプティマイザ（`lbfgs`, `rfo`）。`all` ではスキャン・アライメントの実行も含む。`path-opt` の `result.json` と `all` の `summary.json` にも記録 |
| `preopt_requested` / `preopt_converged` | bool / bool \| null | 端点事前最適化を実行したか、および全端点が収束したか。読み取れない端点があれば `null`。`all` はこの事前収束情報を使います。ただし、要求した最終 TS 最適化と両端点最適化がすべての反応区間で収束した場合は、最終結果で判定します。元のフィールドは保持します |
| `segments` | object[] | セグメントごとの障壁、反応エネルギー、結合変化 |
| `energy_diagrams` | object[] | エネルギーダイアグラム |
| `mlip_backend` | string | バックエンド名（`uma`, `orb`, `mace`, `aimnet2`, `custom`） |
| `mlip_model` | string \| null | バックエンドと分離して記録するモデル/checkpoint名 |
| `mlip_precision` | string \| null | 実効`fp32` / `fp64`。custom calculatorではnull |
| `charge` | int | モデル領域の電荷 |
| `spin` | int | モデル領域のスピン多重度 |
| `environment` | object | ハードウェア情報 |
| `references` | object[] | 解決済みworkflowで実際に使った手法の `{method, citation, doi}` record。同じreference setを `summary.log` と最終標準出力の末尾（elapsed time直前）にまとめて出力します。 |

`all` はさらに以下を含みます。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_segments_reactive` | int | bridge 以外の反応セグメント数。 |
| `rate_limiting_step` | object | 互換性のため維持するキー。各段階の始状態を基準にした局所障壁が最大のセグメントと method。microkinetics に基づく律速段階の判定ではない。 |
| `overall_reaction_energy_kcal` | float | 全体の反応エネルギー。 |
| `post_segments` | list | セグメントごとの TS/IRC/freq/DFT 結果。 |
| `post_segments[].irc` / `.endpoint_assignment` / `.endpoint_opt` | object | 順に IRC 停止診断、端点の向き付け、端点 OPT の収束記録。IRC 停止・結合対応は独立した成功条件にしない。端点の connectivity 情報は機構解釈用に保持する。 |
| `post_segments[].thermo_symmetry` | object | 子 freq が報告した状態別の点群・回転対称 provenance。MEP 実行では R/TS/P、TS-only 実行では E1/TS/E2 を対象とし、有効な対称数 provenance を持つ状態だけを含む。欠けた状態は省略し、どの状態にも有効な provenance が無い場合だけフィールド全体を省略する。 |
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
