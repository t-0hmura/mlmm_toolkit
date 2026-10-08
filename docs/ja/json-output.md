# JSON 出力の一覧

このページでは、`--out-json` で書き出す `result.json` と `summary.json` の欄（key）を、全コマンドに共通の欄とコマンドごとの欄に分けて示します。

## `--out-json` フラグ

ML/MM 計算機を使用する主なサブコマンドとレポート系サブコマンドは `--out-json / --no-out-json`（デフォルト: 無効）に対応しています。有効にすると、通常の出力と同じ場所に `result.json` と `summary.json` を書き出します。2 つのファイルの中身は同じなので、`result.json` を読みます。

```bash
mlmm opt -i r_complex_layered.pdb --parm7 real.parm7 -q 0 -m 1 \
  --max-cycles 5 --out-json --out-dir result_opt
cat result_opt/result.json | python -m json.tool
```

`result_opt/result.json` を開き、まず `execution_status` と `scientific_status` を読みます。`opt`、`tsopt`、`path-opt` は `optimization_status` も記録します。

## 共通エンベロープ

どの結果ファイルにも下の欄があります。「任意」と書いた欄は、そのコマンドが対応するデータを持つときだけ出ます。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `schema_version` | string | ファイルのスキーマの版。版が上がると構造が変わったことを示します。 |
| `command` | string | 単体のコマンドはサブコマンド名（例: `"opt"`）、`all` / `path-search` の要約はコマンドライン全体を記録します |
| `mlmm_version` / `mlmm_toolkit_version` | string | パッケージバージョン（単体のコマンドの結果は `mlmm_version`、`all` / `path-search` の要約は `mlmm_toolkit_version`） |
| `execution_status` | string | 実行の完了状況: `completed` / `failed`。 |
| `scientific_status` | string | 結果の利用可否: `success` / `partial` / `failed`。 |
| `run_id` | string | 任意。現在の呼び出しの UUID。MCP サーバーからコマンドを起動したときと、`all` の実行（各段を含む）で書かれます。 |
| `elapsed_seconds` | float | 任意。実行時間（秒）。時間を記録しないコマンドでは省きます |
| `environment` | object | ハードウェア情報（下表参照） |

ML/MM 計算機を評価するコマンドは、さらに以下を記録します。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `mlip_backend` | string \| null | バックエンドの識別子（`uma`、`orb`、`mace`、`aimnet2`、`dft`、`custom`）。`dft` の結果は `dft`、描画だけのコマンドが計算機を評価しなかった場合は null |
| `mlip_model` | string \| null | 正確なモデル/チェックポイント名。`--calc-file` では `filename:factory` |
| `mlip_model_label` | string \| null | 正確な識別子から導出した論文表記用のモデル名 |
| `mlip_task` | string \| null | 複数ドメインのモデルで使ったバックエンドのタスク（UMA では `calc.uma_task_name` で変えない限り `omol`）。正確な識別子は `mlip_model` に保持 |
| `mlip_precision` | string \| null | 実効精度の共通表記（`fp32` / `fp64`）。dtype をユーザーのコードが管理する自作の計算機では null |
| `mm_backend` | string \| null | MM のエネルギー/Hessian のバックエンド（`hessian_ff` / `openmm`）。描画だけのコマンドが計算機を評価しなかった場合は null |
| `link_atom_method` | string \| null | リンク原子の置き方（`scaled` / `fixed`）。描画だけの出力では null |
| `use_cmap` | bool \| null | CMAP 項を有効にしたか。描画だけの出力では null |

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

### 実行と要求段階の完了状況

すべての結果に `execution_status` と `scientific_status` を出します。複数段階の計算と scan は、各段の結果も下の欄に残します。必要な最適化や計算が欠けていれば、結果は未完了です。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `execution_status` | string | `completed` または `failed`。`all` では、IRC の後の端点の最適化が収束しなくても `completed` のままで、エラーで止まると `failed` になります。 |
| `scientific_status` | string | 要求した段がすべて収束すれば `success`、そうでなければ `partial` か `failed`。`all` では、n_imag ≥ 2 の TS は `partial` になり、n_imag = 0 の TS は IRC の前で止まるので、`success` は n_imag = 1 を意味します。有効な TS と、片方の端点の最適化の失敗の組み合わせも `partial` です。単体の `tsopt` は収束だけを見て n_imag を見ないので、`n_imaginary_modes` を読んでください。 |
| `scientific_status_reasons` | string[] | 利用できない、または欠落した段の理由。正常終了時は省略されます。 |
| `expected_item_ids` / `observed_item_ids` | string[] | 欠けた作業を検出するための、期待された段と観測された段の ID。 |
| `stage_outcomes` | object[] | `stage`、`item_id`、`required`、`executed`、`converged`、`usable`、`reason`、`artifacts` を持つ段階別の結果。 |
| `point_outcomes` | object[] | `point_id`、`executed`、`converged`、`energy_valid`、`artifact_written`、`seed_eligible`、`reason` を持つ scan の点ごとの結果。 |

### エラーエンベロープ（`execution_status == "failed"` のとき）

| フィールド | 型 | 説明 |
|-----------|------|------|
| `error` | string | エラーメッセージ |
| `error_type` | string | 例外クラス名（例: `"OptimizationError"`） |
| `error_class_chain` | list[string] | 例外のクラスとその親クラスの名前を、具体的なものから順に並べたもの（例: `["OptimizationError", "RuntimeError", "Exception", "BaseException"]`）。エージェントはテキストを解析せずに階層を照合できます |
| `error_module` | string | 例外クラスが定義されたモジュール |
| `error_label` | string | 高レベルの CLI ステージラベル（例: `opt` は `"optimization"`、`tsopt` は `"TS optimization"`） |

## エラー処理

`opt`、`tsopt`、`freq`、`irc`、`sp`、`dft`、`scan`、`scan2d`、`scan3d`、`path-opt`、`path-search` が出力ディレクトリを用意した後に例外で止まった場合は、`--out-json` が無くても `"execution_status": "failed"` と `"error_type"` を含む `result.json` と `summary.json` を書き出します。それより前の失敗は[使用上の注意点](#使用上の注意点)を参照してください。

最適化が収束せずに終わった場合、`result.json` には `"optimization_status": "not_converged"` が記録されます。再実行の前に何を変えるかは、{ref}`トラブルシューティングの計算 / 収束 <ja-calculation--convergence>` を参照してください。

オプティマイザは `"optimization_status": "stalled"` を返すこともあります。これは、力/ステップの収束基準を満たさないまま、設定ウィンドウにわたってエネルギーが減少しなくなった状態（エネルギープラトー）です。停滞は未収束の一種で、`converged` にはなりません。`stop_reason` にはエネルギーの範囲、ウィンドウ、満たせなかった基準が記録されます。`--flatten` を指定した `opt` と `tsopt` は、`--max-cycles` の残りがあれば、停滞の後も flatten ループを実行します。マイクロイテレーションでは、マクロステップと直近の MM 緩和の両方が収束したときだけ収束とみなします。MM 緩和が停滞しても、その時点で力が閾値を下回っていれば収束として扱います。

## サブコマンド別スキーマ

### `sp`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `stage` | string | `"sp"` |
| `input` | string | 入力構造のパス |
| `real_parm7` | string | 全系の Amber トポロジーのパス |
| `charge` / `spin` | int / int | ML 領域の電荷とスピン多重度 |
| `energy_au` | float | ONIOM の一点エネルギー (Hartree) |
| `forces_path` | string | `forces.npy` のパス |
| `hessian_path` | string \| null | `hessian.npy` のパス。`--hess` 無指定時は null |
| `elapsed` | string | 人間が読める形の経過時間 |

### `opt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `optimization_status` | string | `"converged"` / `"not_converged"` / `"stalled"`（エネルギープラトー。[エラー処理](#エラー処理)を参照） |
| `stop_reason` | string | オプティマイザが収束せずに止まったとき（`stalled` か `not_converged`）だけ出力。エネルギープラトーの範囲・ウィンドウや満たせなかった基準など、止まった理由を記録 |
| `energy_hartree` | float | 最終 ONIOM エネルギー (Hartree) |
| `n_opt_cycles` | int | 最適化サイクル数 |
| `opt_mode` | string | `"grad"`、`"hess"`、`"lbfgs"`、`"rfo"` のいずれか |
| `charge` | int | ML 領域の電荷 |
| `spin` | int | ML 領域のスピン多重度 |
| `n_atoms` | int | 全原子数（全層） |
| `n_freeze_atoms` | int | 凍結原子数 |
| `thresh` | string | 収束閾値プリセット名 |
| `max_cycles` | int | 最大サイクル数 |
| `input_file` | string | 入力ファイル名 |
| `final_max_force` | float | 最終 max gradient (Hartree/Bohr) |
| `final_rms_force` | float | 最終 RMS gradient |
| `final_max_step` | float | 最終 max 変位 (Bohr) |
| `final_rms_step` | float | 最終 RMS 変位 |
| `convergence_thresholds` | object | プリセットの収束閾値の数値 |
| `files` | object | 出力ファイルマップ |
| `rigid_projection` | object \| null | 任意。`--flatten` を実行したときに出ます。[剛体モードの射影の記録](#剛体モードの射影の記録)を参照 |

### `tsopt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `flatten_requested` / `flatten_enabled` | bool | flatten 反復を設定したかどうか |
| `flatten_skip_reason` | string \| null | それ以上 flatten のステップを行わなかった理由（該当時） |
| `optimization_status` | string | 数値オプティマイザの結果: `"converged"` / `"not_converged"` / `"stalled"`。鞍点の次数とは独立 |
| `saddle_validation` | string | 最後の厳密な PHVA（部分 Hessian 振動解析）による `"first_order"` / `"higher_order"` / `"no_imaginary"` / `"unavailable"` |
| `saddle_order_verified` | bool | `saddle_validation: "first_order"` の場合だけ `true` |
| `hessian_status` | string | `"completed"` / `"failed"` / `"skipped"` / `"unavailable"`。失敗理由は `hessian_error` |
| `reaction_mode_index` | int\|null | `all` が IRC でたどる、厳密な PHVA の負の固有値のモード。参照方向にそろったモードが無いときは 0 番を使ってその旨を記録し、反応のモードであることは確かめていません |
| `reaction_mode_frequency_cm` | float\|null | 選んだ負のモードの振動数 |
| `reaction_mode_source` | string\|null | モードの選び方（参照方向にそろえたか、0 番の代用か） |
| `energy_hartree` | float \| null | TS エネルギー (Hartree)。最終エネルギーの評価に失敗した場合は `null`（有限でない数はすべて `null` で書かれます）で、そのとき実行と結果の状態は `failed` |
| `n_imaginary_modes` | int\|null | 虚振動モードの数。PHVA を実行しなかった場合は `null` |
| `n_negative_modes` | int\|null | 大きさを問わない負の振動数の数（閾値以内も含む）。`n_imaginary_modes` と並べて見る診断用の値で、PHVA を実行しなかった場合は `null` |
| `imaginary_frequencies_cm` | float[]\|null | 虚振動数 (cm⁻¹, 負の値)。PHVA を実行しなかった場合、`--skip-final-freq` では `[]`、それ以外は `null` |
| `frequency_zero_cutoff_cm` / `imaginary_mode_criterion` / `imaginary_frequency_threshold_cm` | float / string / float | 既定値は `5.0`、`"frequency_cutoff_cm"`、`-5.0` で、ν < −5.00 cm⁻¹ だけを虚振動として数えます。 |
| `opt_mode` | string | `"grad"`、`"hess"`、`"dimer"`、`"rsprfo"`、`"rsirfo"`、`"trim"` のいずれか。`hess` は RS-P-RFO を選択 |
| `opt_mode_requested` | string | CLI で要求したプリセット |
| `optimizer` | string | 実際に使用したオプティマイザのアルゴリズム |
| `n_atoms` | int | 全原子数 |
| `n_opt_cycles` | int | 最適化サイクル数 |
| `charge` | int | ML 領域の電荷 |
| `spin` | int | ML 領域のスピン多重度 |
| `reference_mode_file` | string\|null | `--ref-mode` で渡した、経路から得たモードのファイル。Hessian を使うオプティマイザのみ |
| `safeguards` | object | Hessian を使う TS 最適化の診断: 却下したステップと回復、厳密な鞍点の確認、目標モードの同一性 |
| `rigid_projection` | object | Dimer、flatten、最後の鞍点解析の剛体モードと Hessian の記録。[剛体モードの射影の記録](#剛体モードの射影の記録)を参照 |
| `files` | object | final geometry と振動モードのファイル。`--dump-hess` で書いたときは `hessian_npy`（絶対パス）を含む |

TS 最適化が成功すると、反応モードの虚振動が 1 つ出ます。このとき `saddle_validation: "first_order"`、`n_imaginary_modes: 1` です。最後の PHVA は、オプティマイザが収束したときかエネルギープラトーで止まったとき（`stalled`）に実行し、それ以外では実行しなかったこと（skipped）を記録します。PHVA が失敗した場合は、final geometry を残したまま `hessian_status: "failed"` を記録します。`optimization_status` と `saddle_validation` は独立しているので、収束した実行が `saddle_validation: "higher_order"` で終わることもあり、それは一次の TS ではありません。`all` が IRC に進む条件は [tsopt の TS の判定](tsopt.md#ts-の判定) を参照してください。

### `freq`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_modes` | int | 基準振動モードの総数 |
| `n_imaginary` | int | 虚振動モードの数 |
| `n_negative_modes` | int | 大きさを問わない負の振動数の数（閾値以内も含む） |
| `frequencies_cm` | float[] | 全振動数 (cm⁻¹) |
| `imaginary_frequencies_cm` | float[] | 負の振動数のみ |
| `thermochemistry` | object\|null | 熱化学データ（下表参照） |
| `charge` | int | ML 領域の電荷 |
| `spin` | int | ML 領域のスピン多重度 |
| `n_atoms` | int | 全原子数 |
| `n_freeze_atoms` | int | 凍結原子数 |
| `files` | object | 出力マップ。`--dump-hess` で書いたときは `hessian_npy`（絶対パス）を含む |
| `rigid_projection` | object | 振動解析と熱化学で使った剛体モードと Hessian の記録。`--dump` では `thermoanalysis.yaml` にも書きます |

**`thermochemistry`**:

| フィールド | 型 | 単位 |
|-----------|------|------|
| `temperature_K` | float | K |
| `pressure_atm` | float | atm |
| `point_group` | string | 自動検出した分子の点群 |
| `point_group_source` | string | `"auto"` または保守的なフォールバックを示す `"auto-fallback"` |
| `symmetry_number` | int | 外部回転の対称数 |
| `symmetry_number_source` | string | `"auto"`、`"auto-fallback"`、`"config"`、`"override"` |
| `electronic_energy_ha` | float | Hartree（報告される `E + G_corr = G` の `E`） |
| `zpe_ha` | float | Hartree |
| `thermal_correction_energy_ha` | float | Hartree |
| `thermal_correction_enthalpy_ha` | float | Hartree |
| `thermal_correction_free_energy_ha` | float | Hartree |
| `sum_EE_and_ZPE_ha` | float | Hartree |
| `sum_EE_and_thermal_energy_ha` | float | Hartree |
| `sum_EE_and_thermal_free_energy_ha` | float | Hartree |
| `E_thermal_cal_per_mol` | float | cal/mol |
| `Cv_cal_per_mol_K` | float | cal/(mol K) |
| `S_cal_per_mol_K` | float | cal/(mol K) |

### `irc`

IRC は `execution_status` と `scientific_status` を出し、方向ごとの停止理由と軌跡を残します。`all` は、その後の端点の最適化を `endpoint_opt` に記録します。連結した経路は最初のフレームから TS を通って最後のフレームまで続き、単体の `irc` はどちらの端が反応物でどちらが生成物かを決めません。IRC がどう止まったかは `scientific_status` に入りません。単体の `irc` は、積分が通常どおり終われば `completed` / `success` で、`all` は TS と端点の最適化から判定します。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_frames_forward` / `n_frames_backward` / `n_frames_total` | int | IRC フレーム数 |
| `forward_short_branch` / `backward_short_branch` | bool | サイクル上限に達する前に 3 フレーム以内で止まった分岐。診断用のみ |
| `energy_first_hartree` | float | 連結経路の最初のフレームのエネルギー |
| `energy_ts_hartree` | float | TS エネルギー |
| `energy_last_hartree` | float | 連結経路の最後のフレームのエネルギー |
| `endpoint_energy_orientation` | string | `"finished_first_to_finished_last"` |
| `forward_requested` / `backward_requested` | bool | 各方向を要求したか |
| `forward_integration_converged` / `backward_integration_converged` | bool\|null | RMS 勾配の停留判定が働いて止まったか。診断専用で、`--never-stop` はこの判定を迂回するため常に `false`。分岐が TS から下り方向に離れ、かつこの判定を満たしたかは、`*_downhill_departure_valid` と合わせて確かめます |
| `forward_downhill_departure_valid` / `backward_downhill_departure_valid` | bool\|null | TS から下り方向に離れたことを確認できたか |
| `forward_integration_stop_reason` / `backward_integration_stop_reason` | string\|null | 数値的な積分が失敗した場合だけ空でない理由 |
| `never_stop` | bool | `--never-stop`（物理的な端点での停止を迂回する）を有効にしたか |
| `never_stop_energy_bypasses` | int | 実際に迂回した、エネルギー上昇または 1 ステップのエネルギー変化による停止の回数 |
| `rigid_projection` | object | 初期/更新 Hessian の剛体モードと Hessian の記録。[剛体モードの射影の記録](#剛体モードの射影の記録)を参照 |
| `rigid_projection.hessian_source` | string | 初期 Hessian の出どころ。`"file"`（`--read-hess`）、`"cache"`（同じ実行の前の段）、`"fresh"`（新規計算。`irc.hessian_init` が `calc` 以外の場合に EulerPC が作る Hessian も含む） |
| `bond_changes` | object | 最初→最後の方向の `{formed: [...], broken: [...]}`。比較できない場合は省略 |
| `bond_changes_direction` | string | `bond_changes` がある場合は `"finished_first_to_finished_last"` |
| `files` | object | 軌跡と端点のファイル（XYZ と、利用できる場合は PDB/CIF 版） |

### `scan`

| フィールド | 型 | 説明 |
|---|---|---|
| `scan_opt_mode` | string | 指定した `grad` または `hess` のプリセット |
| `scan_optimizer` | string | 実際のオプティマイザ: `lbfgs` または `rfo` |
| `n_stages` | int | スキャンステージ数 |
| `stages` | object[] | ステージごとの結果 |
| `charge` | int | ML 領域の電荷 |
| `spin` | int | ML 領域の多重度 |
| `files` | object | 出力ファイル |

`stages[]` は `n_steps`、`converged`、`pairs_1based`、`energies_hartree`、`final_energy_hartree`、`bond_changes`、`optimizer_status`（`converged` / `not_converged` / `stalled`）を含みます。その段の最後の最適化が収束せずに止まった場合は `stop_reason` も記録します。

### `scan2d` / `scan3d`

| フィールド | 型 | 説明 |
|---|---|---|
| `n_grid_points` | int | 格子点数 |
| `n_points_attempted` | int | 新規計算で試行した格子点数（事前最適化を除く） |
| `n_points_usable` | int | 明示的に収束し、有限のエネルギー・座標と構造ファイルを持つ点の数 |
| `grid_points` | object[] | 格子インデックス、距離、エネルギー、収束結果、`geometry_file` の対応 |
| `current_output_paths` | string[] | 今回の実行が書いた CSV・HTML・PNG と格子点の構造。前の実行で残ったファイルは並びません |
| `pair1`、`pair2`（`pair3`） | object | `{i, j, low, high}` |
| `min_energy_hartree` | float | 表面上の最小エネルギー |
| `charge` | int \| null | ML 領域の電荷。作図のみの `scan3d --csv` では null |
| `spin` | int \| null | ML 領域の多重度。作図のみの `scan3d --csv` では null |
| `files` | object | CSV・プロットファイル |

新規の `scan2d` / `scan3d` には共通の計算機の欄も記録します。作図のみの `scan3d --csv` では同じ欄を null で書きます。作図のみの結果は `n_points_attempted` を出さず、`n_points_usable` は CSV にすべての点の収束と構造ファイルの記録がある場合だけ出します。

### `path-opt`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `converged` | bool \| null | 収束判定: エンジン自身の収束シグナルによる `true` / `false`。読み取れない場合は `null`（`optimization_status` は `"completed"` となり、収束を主張しない） |
| `mep_mode` | string | `"dmf"` / `"gsm"` |
| `image_energies_hartree` | float[] | 全イメージのエネルギー |
| `n_images` | int | イメージ数 |
| `hei_index` | int | 最高エネルギーイメージのインデックス |
| `barrier_kcal` | float | 前方障壁 (kcal/mol) |
| `delta_kcal` | float | 反応エネルギー (kcal/mol) |
| `files` | object | 軌跡と HEI のファイル |

### `path-search`

`path-search` には `--out-json` フラグがありません。`summary.json` を書き出し、その欄は [`summary.json` (`path-search` / `all`)](#ja-summary-json-path-search-all) にあります。

### `dft`

> **注:** `--out-json` 指定時、`dft` は SCF 収束・非収束の両方で `result.json` と `summary.json` を書き、非収束時は `scientific_status: "failed"`、`converged: false` を記録します。SCF が収束しなかった実行は終了コード 1 で終わります。未処理の例外では標準のエラーエンベロープを書きます。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `converged` | bool | SCF が収束したか |
| `energy_hartree` / `energy_kcal_per_mol` | float | ML 領域の DFT エネルギー。`model_dft_energy_*` と同じ値 |
| `model_dft_energy_hartree` / `model_dft_energy_kcal_per_mol` | float | ML 領域の DFT エネルギー |
| `total_dft_mm_energy_hartree` / `total_dft_mm_energy_kcal_per_mol` | float | 組み合わせた DFT/MM エネルギー |
| `xc_functional` | string | 汎関数 |
| `basis_set` | string | 基底関数 |
| `used_gpu` | bool | GPU を使ったか |
| `used_lowmem` / `lowmem_requested` | bool | 実際の低メモリ状態と、要求した低メモリ状態 |
| `dft_settings` / `dft_resources` | object | 正規化した計算設定と、実際に使ったホストの資源 |
| `effective_ecp` | string/object \| null | PySCF へ渡した ECP |
| `embedding` | object | 点電荷の埋め込み: 有効かどうか、cutoff、電荷の数、電荷を識別するダイジェスト |
| `charges` | object | `{mulliken, lowdin, iao}` 原子ごとの配列 |
| `spin_densities` | object | `{mulliken, lowdin, iao}` 原子ごとの配列 |
| `n_atoms` | int | ML 領域の原子数 |
| `grid_level` | int | DFT のグリッドレベル |
| `conv_tol` | float | SCF の収束閾値 |
| `max_cycle` | int | YAML と CLI を適用した後の SCF 最大反復数 |
| `engine` | string | 実行に使ったエンジンの表記（`pyscf(cpu)`、`gpu4pyscf`、または低メモリの GPU 版） |
| `charge` | int | ML 領域の電荷 |
| `spin` | int | ML 領域のスピン多重度 |
| `input_file` | string | 入力構造のパス |
| `files` | object | `{"result_yaml": "result.yaml"}` |

### `extract`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_atoms_extracted` | int | 抽出後の原子数 |
| `total_charge` | float | 合計電荷 |
| `protein_charge` | float | タンパク質電荷 |
| `ligand_total_charge` | float | リガンド電荷合計 |
| `ion_total_charge` | float | イオン電荷合計 |
| `unknown_residue_charges` | object | `{残基名: 電荷}` |
| `center` | string | 基質指定（`-c` の値）: PDB パス、残基 ID リスト（例 `'A:123,B:456'`）、または残基名リスト（例 `'GPP,MMT'`） |
| `radius` | float | 抽出半径 (Å) |
| `input_files` | string[] | 入力 PDB パス |
| `n_atoms_raw` | int | 主鎖の除外や切り詰めの前の、選んだ残基の原子数（入力構造全体ではない） |
| `n_link_hydrogens` | int | 切断結合に付加したリンク H 原子数 |
| `files` | object | 書き出したファイルのマップ（入力ごとのポケット PDB 等） |
| `exclude_backbone` | bool | 実行時の `--exclude-backbone` の値 |
| `include_h2o` | bool | 実行時の `--include-h2o` の値 |
| `ligand_charge_input` | string | 与えた `-l/--ligand-charge` の引数 |
| `ion_charges` | array | イオン残基の `[残基名, 電荷]` ペアのリスト |

### `trj2fig`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_frames` | int | 軌跡のフレーム数 |
| `min_energy_hartree` / `max_energy_hartree` | float | フレームエネルギーの最小値と最大値 |
| `energy_source` | string | `"trajectory_comment"` または `"mlip_recomputed"` |
| `mlip_backend` / `mlip_model` / `mlip_precision` | string \| null | エネルギーの再計算に使った計算機。コメントモードではすべて null |
| `charge` / `multiplicity` | int \| null | 再計算に使った電荷とスピン多重度。コメントモードでは null。片方だけを指定した場合、もう片方は 0 または 1 |
| `output_files` | string[] | すべての出力のパスを順序どおりに並べたもの。別ディレクトリに同名ファイルがあっても保持される |
| `files` | object | ベース名からパスへの対応表。同じベース名の出力が 2 つあると一方だけが残るので、`output_files` を使ってください |

`-q/--charge` または `-m/--multiplicity` のいずれかを指定すると、選択した MLIP で全フレームを再計算します。この再評価は MLIP だけで行い、トポロジーや ML 領域の入力を受け取らず、ONIOM エネルギーは計算しません。

### `energy-diagram`

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_points` | int | エネルギーデータ点の数 |
| `files` | object | 出力ダイアグラムのファイル名からパスへの対応表 |

### `bond-summary`

`--json` を付けると、`bond-summary` は JSON を**標準出力**に出し、`result.json` は書きません。残すには標準出力をリダイレクトします。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `execution_status` / `scientific_status` | string / string | すべての組を比較できれば `completed` / `success`。比較できない組があると `execution_status` は `failed`、`scientific_status` は `partial` か `failed` で、コマンドは終了コード 1 で終わります |
| `comparisons` | object[] | 組ごとの比較。`structure_a`、`structure_b`、`bonds_formed` / `bonds_broken`（数）、`formed` / `broken`（結合ごとの記録: `atom_i`、`atom_j`、`element_i`、`element_j`、`distance_a_angstrom`、`distance_b_angstrom`）。比較できなかった組は代わりに `error` を持ちます |

### 剛体モードの射影の記録

`freq`、`irc`、`tsopt` の結果は `rigid_projection` object を含み、`opt` では `--flatten` 実行時に含みます。`freq --dump` は同じ object を `thermoanalysis.yaml` にも書きます。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `treatment` | string | 固定の剛体モード処理: `"constrained"` |
| `algorithm` | string | 射影の方法の名前 |
| `effective_rank` | int | 動ける原子の Hessian から除いた剛体方向の数 |
| `full_rigid_rank` | int | 凍結原子を考える前の、系全体の剛体運動のランク |
| `frozen_constraint_rank` | int | 凍結原子を動かさない条件で除かれたランク |
| `svd_rtol` | float | ランクの判定に用いる相対 SVD 許容値 |
| `active_atom_count` / `frozen_atom_count` | int | 動ける原子と凍結原子の数 |
| `active_atoms` / `frozen_atoms` | int[] | 動ける原子と凍結原子の 0 始まりのインデックス |
| `hessian_space` | string | 入力 Hessian 空間: `"full"` / `"active"` |
| `hessian_source` / `source` | string | Hessian の出どころ。`freq`/`irc` は `hessian_source` で、`"file"`（`--read-hess`）、`"cache"`（同じ実行の前の段）、`"fresh"`（新規計算）のいずれか。`opt`/`tsopt` は `source` を使用 |
| `hessian_shape` / `raw_hessian_shape` | int[2] | 入力 Hessian の形状。`freq`/`irc` は `hessian_shape`、`opt`/`tsopt` は `raw_hessian_shape` を使用（`freq` は両方を記録） |
| `near_zero_mode_count` / `near_zero_frequencies_cm` | int / float[] | ±`frequency_zero_cutoff_cm`（既定 5.00 cm⁻¹）以内のモードの数と値。これらのモードは全振動数の一覧にも入ります |

`constrained` は、凍結原子を動かさない系全体の剛体運動だけを除きます。詳しくは [freq](freq.md#凍結境界での剛体モード) を参照してください。

(ja-summary-json-path-search-all)=
## `summary.json` (`path-search` / `all`)

`all` と `path-search` は、より多くの欄を持つ `summary.json` を書き出します。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `execution_status` / `scientific_status` | string / string | 実行の完了度と、要求した数値最適化・計算段階の完了度。 |
| `scientific_status_reasons` | string[] | 要求した結果の欠損・未収束などの理由。正常終了時は省略されます。 |
| `pipeline_stop` | object \| 不在 | 早期停止時のみ存在。`stage` は `post`（`reason` は `no_segments` / `no_reactive_segment`）、`before_irc`（TSOPT の理由と `segment`・`tsopt_result`）、または `endpoint_opt`（`endpoint_execution_failed` と端点別 `failures`）。`summary.log` では `Pipeline stop` |
| `expected_item_ids` / `observed_item_ids` | string[] | 期待された段と観測された段の ID。 |
| `config` | object | 実効設定。`mep_mode` は GSM/DMF、`ts_opt_mode` / `endpoint_opt_mode` は設定した後処理のプリセットを示す。一般の `opt_mode*` は実際に使った CLI の値を記録する。`path_opt_mode` は端点の事前最適化に使う単一構造オプティマイザであり（`preopt` を参照）、MEP の経路アルゴリズムではない。 |
| `n_segments` | int | セグメント数 |
| `search_max_depth` | int | 実効の再帰分割階層上限。`0` は分割無効 |
| `path_optimizers` | string[] | 経路の準備・精密化で実際に使用した単一構造オプティマイザ（`lbfgs`, `rfo`）。`all` ではスキャン・アライメントの実行も含む。`path-opt` の `result.json` にも記録 |
| `preopt_requested` / `preopt_converged` | bool / bool \| null | 端点の事前最適化を実行したか、および全端点が収束したか。読み取れない端点があれば `null`。`all` では `preopt_converged` も `scientific_status` に数えます。ただし、要求した最終 TS 最適化と両端点の最適化がすべての反応区間で収束した場合は数えません。欄そのものは残します |
| `segments` | object[] | セグメントごとの `index`（1 始まり）、`tag`（`seg_001` のような名前。共有結合が変わらないねじれ（kink）のセグメントは名前に `kink` を含みます）、`kind`（反応セグメントは `"seg"`、[ブリッジセグメント](path-search.md#処理の仕組みと計算仕様)は `"bridge"`、TS-only モードは `"tsopt"`）、`converged`（そのセグメントを作った最適化がすべて収束したか）、`barrier_kcal`、`delta_kcal`、`bond_changes`（ブリッジセグメントは `""`）。`barrier_kcal` は TS 最適化の前の MEP 上の障壁です。TS-only モードでは TS − R で、R は IRC の両端のうちエネルギーが高いほうです（向きの名前で、化学的な反応の向きではありません）。 |
| `energy_diagrams` | object[] | ラベルと kcal/mol の値を持つエネルギーダイアグラム |
| `mlip_backend` | string | バックエンド名（`uma`、`orb`、`mace`、`aimnet2`、`dft`、`custom`） |
| `mlip_model` | string \| null | 正確なモデル/チェックポイント名。`dft` では `FUNCTIONAL/BASIS` |
| `mlip_model_label` | string \| null | 論文表記用のモデル名 |
| `mlip_task` | string \| null | 複数ドメインのモデルで使ったバックエンドのタスク |
| `mlip_precision` | string \| null | 実効の `fp32` / `fp64`。DFT と自作の計算機では null（DFT のエンジンは別の欄に記録） |
| `charge` | int | ML 領域の電荷 |
| `spin` | int | ML 領域のスピン多重度 |
| `environment` | object | ハードウェア情報 |
| `references` | object[] | 実行で実際に使った手法の `{method, citation, doi}` の記録。同じ文献の組を `summary.log` と最終標準出力の末尾（経過時間の直前）にまとめて出力します。 |

`all` はさらに以下を含みます。

| フィールド | 型 | 説明 |
|-----------|------|------|
| `n_segments_reactive` | int | bridge 以外の反応セグメント数 |
| `rate_limiting_step` | object | 反応セグメントのうち局所障壁が最大のもの。`{segment, barrier_kcal, method}` と、`mep_barrier_kcal`（MEP 上の障壁。TS-only モードでは無し）を持ちます。障壁は、すべてのセグメントで得られる最も高い水準の手法（`DFT//MLIP/MM_Gibbs` > `DFT` > `MLIP_Gibbs` > `MLIP` > `MEP`。`--backend dft` では `MLIP_Gibbs` と `MLIP` の代わりに `DFT/MM_Gibbs` と `DFT/MM`）で取ります。microkinetics に基づく律速段階の判定ではありません。 |
| `overall_reaction_energy_kcal` | float | 全体の反応エネルギー |
| `overall_reaction_energy_method` | string | 全体の反応エネルギーの手法。`rate_limiting_step.method` と同じ一覧から、それより高くない水準を使います |
| `post_segments` | list | セグメントごとの TS/IRC/freq/DFT 結果 |
| `post_segments[].tsopt.n_imaginary_modes` / `.imaginary_frequencies_cm` | int / float[] | 最適化した TS の n_imag と虚振動数 (cm⁻¹, 負の値) |
| `post_segments[].mlip` / `.gibbs_mlip` / `.dft` / `.gibbs_dft_mlip` | object | 1 つのレベルでの R・TS・P のエネルギーと、`barrier_kcal`・`delta_kcal`・`energies_kcal`。順に ML/MM の電子エネルギー（`--tsopt`）、ML/MM の Gibbs エネルギー（`--thermo`）、DFT のエネルギー（`--dft`）、DFT//ML/MM の Gibbs エネルギー（`--thermo` と `--dft`）です |
| `post_segments[].tsopt.energy_valid` / `.structure_valid` | bool | `energy_valid`: 最終 TS のエネルギーが有限の数。`structure_valid`: 最終 TS の構造ファイルがあり、座標が有限。これらの確認のために追加の Hessian や最適化は実行しません |
| `post_segments[].tsopt.n_opt_cycles` / `.max_cycles` | int / int\|null | TS 最適化で実行したサイクル数と設定上限。通常の非収束時にも記録します。 |
| `post_segments[].irc` / `.endpoint_assignment` / `.endpoint_opt` | object | 順に IRC の停止の診断、端点の向き付け、端点の最適化の収束記録。TS-only モードでは `endpoint_assignment.policy` は `higher_energy_endpoint_as_reactant`、`chemical_direction_known` は false。`endpoint_opt.reactant` と `.product` に `optimization_status`、`n_opt_cycles`、`max_cycles`、`stop_reason`（ある場合）を記録します。IRC がどう止まったかと結合変化が合うかは `scientific_status` に入りません。最適化した端点が意図した R と P かは、利用者が確かめます。 |
| `post_segments[].thermo_symmetry` | object | 各段の freq が報告した状態別の点群と回転の対称数の記録。TS-only モードを含む R/TS/P を対象とし、有効な対称数の記録を持つ状態だけを含みます。欠けた状態は省き、どの状態にも有効な記録が無い場合だけ欄全体を省きます。 |
| `key_output_files` | object | 今回の実行の出力の索引。ルートのファイルはファイル名 → 説明、各 `seg_NN` は `{description, files}` で、`files` はそのセグメントディレクトリからの相対パス。 |
| `current_output_paths` | string[] | `--out-dir` からの相対パスを並べたリスト。今回の実行が書いたファイルだけを含みます。 |

## 使用例

### Python

`opt --out-json` が書き出す `result_opt/result.json` を読む例です。

```python
import json

with open("result_opt/result.json") as f:
    result = json.load(f)

status = result.get("optimization_status")
if result["execution_status"] == "failed":
    raise RuntimeError(f"{result['error_type']}: {result['error']}")
elif status == "converged":
    print(f"Energy: {result['energy_hartree']:.6f} Hartree")
elif status in {"not_converged", "stalled"}:
    print(f"Not converged after {result['n_opt_cycles']} cycles")
    print(f"Max force: {result['final_max_force']:.6f}")
else:
    print(f"Status: {status}")
```

### jq

```bash
# 収束確認
jq '{execution_status, scientific_status}' result.json

# path-opt の障壁
jq '.barrier_kcal' result.json

# tsopt の虚振動数
jq '.imaginary_frequencies_cm' result.json

# freq の自由エネルギー
jq '.thermochemistry.sum_EE_and_thermal_free_energy_ha' result.json

# all の各セグメントの障壁（--tsopt の後）
jq '.post_segments[] | {index, barrier_kcal: .mlip.barrier_kcal}' result_all/summary.json
```

## 使用上の注意点

- CLI のオプションや入力を確かめている段階（出力ディレクトリを用意する前）で失敗すると、JSON を書かずに止まることがあります。0 以外の終了コードは失敗として扱い、標準エラー出力かジョブログでメッセージを確認してください。
- `all` と `path-search` は要約の段に着いてから `summary.json` を書くため、早い段階の入力エラーではファイルが残りません。

## 関連ドキュメント

- {ref}`終了コード <ja-exit-codes>` — 終了コードの意味
- [出力ディレクトリのレイアウト](output-layout.md) — `result.json` と `summary.json` を書き出す場所
- [トラブルシューティング](troubleshooting.md) — 失敗・未収束の後に何を変えるか
- [YAML 設定の一覧](yaml-reference.md) — これらのスキーマに現れる設定入力
- [all](all.md), [path-search](path-search.md) — `--out-json` なしで `summary.json` を書き出すサブコマンド
- [opt](opt.md), [sp](sp.md), [tsopt](tsopt.md), [freq](freq.md), [irc](irc.md), [scan](scan.md), [scan2d](scan2d.md), [scan3d](scan3d.md), [path-opt](path-opt.md), [dft](dft.md), [extract](extract.md), [trj2fig](trj2fig.md), [energy-diagram](energy-diagram.md), [bond-summary](bond-summary.md) — `--out-json` を付けたときだけ JSON を出すサブコマンド（`bond-summary` は `--json`）
