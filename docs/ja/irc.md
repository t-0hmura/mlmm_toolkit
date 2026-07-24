# `irc`

`mlmm irc` は ML/MM calculatorを用いた EulerPC ベースの IRC（固有反応座標）積分により、遷移状態から反応物・生成物の方向へ経路を追跡します。最適化された TS が期待どおり反応物と生成物を接続するかを検証したいとき、あるいは下流の熱化学計算 / DFT 単点計算用の反応物 / 生成物構造を生成したいときに使用します。典型的には `tsopt` -> `freq`（**1 つ**の虚振動数モードを確認）-> `irc` というワークフローで実行します。デフォルトでは正方向と逆方向の両方のブランチが計算されます。共通input bridgeはPDB/mmCIFと`geom_loader`対応形式を受け入れます。直接入力または`--ref-pdb`でPDB/mmCIF topologyがあり、変換が有効ならPDB companionを生成し、mmCIF/oversized-PDB bridge入力では元IDを復元したCIF companionも生成します。

## 実行例

最小構成で TS の PDB から実行:

```bash
mlmm irc -i ts.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc
```

正方向のみ実行:

```bash
mlmm irc -i ts.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --no-backward --out-dir ./result_irc_forward
```

ステップサイズを小さくして解析 Hessian を使用:

```bash
mlmm irc -i ts.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --step-size 0.05 \
 --hessian-calc-mode Analytical --out-dir ./result_irc_analytical
```

IRC がほぼ直ちに停止する場合は、まず `--step-size` を小さくします（例:
0.10 から 0.05 Bohr）。検証済みの小さな shoulder でエネルギー上昇/plateau
停止だけが残る場合は、`--never-stop` を明示的に指定できます:

```bash
mlmm irc -i ts.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0 \
 --step-size 0.05 --never-stop --max-cycles 250 -o result_irc_continue
```

このモードでも integrator 収束、非有限値、サイクル上限では停止します。
両方向の軌跡と終点接続を確認してから採用してください。

両ブランチを保持してステップ上限を引き上げ:

```bash
mlmm irc -i ts.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --max-cycles 150 \
 --out-dir ./result_irc_long
```

コマンド形式:

```bash
mlmm irc -i TS_STRUCTURE --parm PARM7 --model-pdb ML_REGION [options]
```

`mlmm irc --help` でコアオプションを、`mlmm irc --help-advanced` で全オプション一覧を表示します。

## 処理の流れ

1. **入力準備** -- TS 構造、Amber トポロジー（`--parm`）、ML 領域定義（`--model-pdb` / `--model-indices`）を読み込み、電荷とスピンを確定します。直接PDB/mmCIF入力または`--ref-pdb`がcompanion出力用topologyを提供します。
2. **ML/MM calculatorの構築** -- `--parm` と `--model-pdb` から ML/MM calculatorを構築します。`-b/--backend` で ML バックエンドを選択し（デフォルト: `uma`）、`--hessian-calc-mode` は MLIP Hessian評価を制御します。v0.3.3 は機械的埋め込みを使用し、電子埋め込みの要求は calculator 構築前に拒否します。
3. **凍結境界の TR 処理** -- `--tr-projection constrained` は、凍結 anchor をすべて動かさない全系剛体運動だけを除去します。一般的な有効 rank は anchor が 0/1/2/非共線の 3 個以上のとき 6/3/1/0 で、実用的な ML/MM 境界では通常 0 です。`legacy-active` は非推奨の比較専用処理で、pass/HOSP 遷移状態認定には使用できません。
4. **IRC 積分** -- EulerPC 積分器が両方向に沿って IRC を伝播します（`--no-forward` または `--no-backward` でブランチを無効化可能）。ステップサイズとサイクル数で積分長を制御します。
5. **出力と変換** -- 軌跡はXYZで書き出されます。PDB/mmCIF topologyが利用可能で`--convert-files`が有効ならPDB companionを生成し、bridge入力では元ID付きCIF companionも生成します。

## 出力

```
out_dir/ (デフォルト: ./result_irc/)
├─ result.json                      # --out-json 時。rigid_projection provenance を含む
├─ <prefix>irc_data.h5              # 低レベルの周期 checkpoint。YAML でのみ opt-in
├─ <prefix>finished_irc_trj.xyz     # 完全 IRC 軌跡（XYZ/TRJ）
├─ <prefix>forward_irc_trj.xyz      # 正方向パスセグメント
├─ <prefix>backward_irc_trj.xyz     # 逆方向パスセグメント
├─ <prefix>finished_irc.pdb         # PDB 変換（入力が .pdb または --ref-pdb 指定時）
├─ <prefix>finished_irc.cif         # bridge入力。元IDを復元
├─ <prefix>forward_irc.pdb          # PDB 変換（入力が .pdb または --ref-pdb 指定時）
├─ <prefix>forward_irc.cif          # bridge入力の順方向CIF
├─ <prefix>backward_irc.pdb         # PDB 変換（入力が .pdb または --ref-pdb 指定時）
├─ <prefix>backward_irc.cif         # bridge入力の逆方向CIF
├─ <prefix>forward_first.xyz        # 正方向 IRC 終点（XYZ、単一フレーム）
├─ <prefix>forward_first.pdb/.cif   # 正方向IRC終点companion（利用可能時）
├─ <prefix>backward_last.xyz        # 逆方向 IRC 終点（XYZ、単一フレーム）
└─ <prefix>backward_last.pdb/.cif   # 逆方向IRC終点companion（利用可能時）
```

`irc.prefix`が空でない場合、EulerPCはファイル名との間に`_`を1つ補います。たとえば
`prefix: trial`は`trial_finished_irc_trj.xyz`を生成し、`result.json.files`にも
正規化後の名前を記録します。

`irc.dump_every` のデフォルトは `null` なので、HDF5 checkpoint は作成されません。
YAML で正の値を指定した場合のみ、現在方向の座標・エネルギー・勾配で周期的に
上書きされます。最終的な双方向 IRC 成果物ではなく、Hessian は含まず、
`result.json.files` にも登録しません。

standalone IRCはstitched pathの`first` / `last`端点と、その方向のbond changesを
記録します。化学的なreactant/product identityは割り当てないため、R/Pの命名前に
端点構造を確認または参照構造と対応付けてください。

主に確認するファイル:

- `result_irc/finished_irc_trj.xyz`
- `result_irc/forward_irc_trj.xyz`

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-b, --backend CHOICE` | ML バックエンド: `uma`（デフォルト）、`orb`、`mace`、`aimnet2`。 | `uma` |
| `--embedcharge/--no-embedcharge` | v0.3.3 では使用不可。旧コマンドを明示的に拒否するためにのみ残されています。 | `False` |
| `--embedcharge-cutoff FLOAT` | 廃止した電子埋め込み経路とともに使用不可。 | — |
| `--cmap/--no-cmap` | model parm7 に CMAP（骨格クロスマップ二面角補正）を含めるかどうか。デフォルト: 無効（Gaussian ONIOM と同一）。 | `--no-cmap` |
| `--hess-device CHOICE` | 初期Hessianの格納・IRC演算のデバイス: `auto`、`cuda`、`cpu`。大規模非凍結系では `cpu` を推奨。 | `auto` |
| `--read-hess PATH` | `mlmm freq --dump-hess`のidentified `.npz`を読み込む。geometry、原子順序、layer選択、active-DOF basisが一致する必要があり、cache／新規計算より優先。 | _None_ |
| `-i, --input PATH` | 構造ファイル（`.pdb`/`.xyz`/`_trj.xyz`/...）。`geom_loader` で読み取り可能な任意の形式。 | 必須 |
| `--parm PATH` | 全酵素/MM 領域の Amber トポロジー。YAML の `calc.real_parm7` が無い場合は必須。 | _None_ |
| `--model-pdb PATH` | ML 領域を定義する PDB。`--no-detect-layer` かつ `--model-indices` 未指定時は必須。 | _None_ |
| `--model-indices TEXT` | ML 領域原子インデックス（カンマ区切り、範囲指定可: `1-10,15`）。`--model-pdb` 省略時に使用。 | _None_ |
| `--model-indices-one-based/--model-indices-zero-based` | `--model-indices` を 1 始まり/0 始まりとして解釈。 | `True`（1 始まり） |
| `--detect-layer/--no-detect-layer` | 入力 PDB の B 因子（`B=0/10/20`）から ML/MM レイヤーを検出。 | `True` |
| `--freeze-atoms TEXT` | 1 始まりの凍結原子インデックスをカンマ区切りで指定。 | _None_ |
| `--tr-projection [constrained\|legacy-active]` | 凍結/部分Hessianの剛体モード処理。`legacy-active` は非推奨の比較専用で、pass/HOSP 遷移状態認定には使用不可。 | `constrained` |
| `-q, --charge INT` | ML 領域/model system の正味電荷。YAML の `calc.model_charge` を上書き。 | _None_（`-l` 未指定時は必須） |
| `-l, --ligand-charge TEXT` | 未知リガンド残基の合計電荷または残基別マッピング（例: `GPP:-3,SAM:1`）。`-q` 省略時に ML 領域の正味電荷を導出。 | _None_ |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。`calc.spin` を上書き。 | `1` |
| `--max-cycles INT` | IRC ステップの最大数。`irc.max_cycles` を上書き。 | `125` |
| `--step-size FLOAT` | ステップ長（Bohr、非質量加重デカルト座標）。`irc.step_length` を上書き。 | `0.10` |
| `--root INT` | 初期変位の虚振動数モードインデックス。`irc.root` を上書き。 | `0` |
| `--forward/--no-forward` | 正方向 IRC を実行。`irc.forward` を上書き。 | `True` |
| `--backward/--no-backward` | 逆方向 IRC を実行。`irc.backward` を上書き。 | `True` |
| `--never-stop/--no-never-stop` | エネルギー上昇/plateau 停止だけを無視。収束、非有限値、最大サイクルでは停止。 | `False` |
| `-o, --out-dir PATH` | 出力ディレクトリ。`irc.out_dir` を上書き。 | `./result_irc/` |
| `--ref-pdb FILE` | `--input`がXYZの場合に使用する参照PDB/mmCIF topology（XYZ座標を保持）。 | _None_ |
| `--convert-files/--no-convert-files` | 参照topologyがある場合のXYZ/TRJ→PDB/CIF companionを切り替え。 | `True` |
| `--hessian-calc-mode CHOICE` | MLIP がHessianを構築する方法（`Analytical` または `FiniteDifference`）。`calc.hessian_calc_mode` を上書き。 | `FiniteDifference` |
| `--workers INT` | UMA predictor worker 数。2 以上は `fairchem-core[extras]` が必要で、解析 Hessian と併用不可。 | `1` |
| `--workers-per-node INT` | UMA 並列 predictor のノード当たり worker 数。 | _None_ |
| `--config FILE` | 明示 CLI 適用前に読み込むベース YAML。 | _None_ |
| `--show-config/--no-show-config` | 解決済み YAML レイヤー/設定を表示して続行。 | `False` |
| `--mm-backend [hessian_ff\|openmm]` | MM バックエンド。Hessian 構築法は `calc.mm_fd` が別に制御します（既定 `true`: 有限差分）。 | `hessian_ff` |
| `--link-atom-method [scaled\|fixed]` | リンク原子配置: scaled（$g$ 係数）または fixed（1.09/1.01 Å）。 | `scaled` |
| `--out-json/--no-out-json` | 機械可読な `result.json` を `out_dir` に書き出し。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画のみ表示。`--help-advanced` に表示。 | `False` |
| `--allow-unverified-hess-state/--no-allow-unverified-hess-state` | charge/多重度を検証できない schema 1 Hessian を許可。`--read-hess` と独立した状態確認が必要。 | `False` |

NPZ の geometry、原子順序、active basis、model charge、多重度は現在の実行と
一致する必要があります。電子状態 identity を持たない schema 1 は
`--allow-unverified-hess-state` の明示的 opt-in が必要で、schema 2 の状態不一致は
常に致命的です。
`result.json["rigid_projection"]["electronic_state_verified"]` が検証結果を記録します。

## YAML 設定

マージ順 **デフォルト < config < 明示CLI < override** でマッピングを提供します。
共有セクションはジオメトリ/計算機キーについて [YAML リファレンス](yaml-reference.md) を再利用します。`irc` では YAML/CLI マージ後に `geom.coord_type` が `cart` に強制されます。`calc.return_partial_hessian` は明示的な YAML 指定が無い場合に `true` がデフォルト適用されます（active-DOF 処理を伴う partial Hessian）。

### CLI から YAML へのマッピング

| CLI オプション | YAML キー |
|------------|----------|
| `--charge` | `calc.charge` |
| `--multiplicity` | `calc.spin` |
| `--tr-projection` | `geom.tr_projection` |
| `--step-size` | `irc.step_length` |
| `--max-cycles` | `irc.max_cycles` |
| `--root` | `irc.root` |
| `--forward` | `irc.forward` |
| `--backward` | `irc.backward` |
| `--never-stop` | `irc.never_stop` |
| `--out-dir` | `irc.out_dir` |
| `--hessian-calc-mode` | `calc.hessian_calc_mode` |

### YAML 例

```yaml
geom:
 coord_type: cart                  # irc では cart に強制（YAML 値は無視）
 freeze_atoms: []                  # 1 始まり凍結原子（CLI/リンク検出とマージ）
 tr_projection: constrained        # legacy-active は非推奨・比較専用
calc:
 model_charge: 0                   # ML 領域/model system の正味電荷
 spin: 1                           # スピン多重度 2S+1
mlmm:
 real_parm7: real.parm7            # Amber parm7 トポロジー
 model_pdb: ml_region.pdb          # ML 領域定義
 backend: uma                      # ML バックエンド (uma/orb/mace/aimnet2)
 embedcharge: false                # 互換性用。true は拒否される
 uma_model: uma-s-1p2              # uma-s-1p2 | uma-m-1p1
 uma_task_name: omol                # UMA タスク名 (backend=uma 時)
 ml_device: auto                   # ML デバイス選択
 hessian_calc_mode: Analytical         # Hessianモード選択
 return_partial_hessian: true      # irc では true に強制（partial Hessian、active-DOF 処理）
irc:
 step_length: 0.1                  # 積分ステップ長
 max_cycles: 125                   # IRC に沿った最大ステップ数
 downhill: false                   # 下り方向のみに追従
 forward: true                     # 正方向に伝播
 backward: true                    # 逆方向に伝播
 never_stop: false                 # energy-rise/plateau 停止のみ無視
 root: 0                           # 基準振動ルートインデックス
 hessian_init: calc                # Hessian初期化ソース
 displ: energy                     # 変位構築方法
 displ_energy: 0.001               # エネルギーベースの変位スケーリング
 displ_length: 0.1                 # 長さベースの変位フォールバック
 rms_grad_thresh: 0.001            # RMS 勾配収束閾値
 hard_rms_grad_thresh: null        # ハード RMS 勾配停止
 energy_thresh: 0.000001           # エネルギー変化閾値
 imag_below: 0.0                   # 虚振動数カットオフ
 force_inflection: true            # 変曲点検出を強制
 check_bonds: false                # 伝播中の結合チェック
 out_dir: ./result_irc/            # 出力ディレクトリ
 prefix: ""                        # ファイル名プレフィックス
 hessian_update: bofill            # Hessian更新方式
 hessian_recalc: null              # Hessian再構築間隔
 max_pred_steps: 500               # 予測子-補正子の最大ステップ数
 loose_cycles: 3                   # 厳密化前のゆるいサイクル数
 corr_func: mbs                    # 相関関数の選択
```

完全なスキーマ（すべての `irc` キーとデフォルト）: [YAML リファレンス](yaml-reference.md#irc-section)。

## 注記

- デフォルトでは両方のブランチを実行します。片方向のみが必要な場合は `--no-forward` または `--no-backward` で一方を無効化します。
- 早期停止時はまず `--step-size` を小さくし、`--never-stop` は経路を確認した上で opt-in してください。
- 全原子凍結では IRC 方向が無いため、明示的なエラーになります。
- `--out-json` 時は `result.json.rigid_projection` に treatment、有効 rank、
  初期Hessian source、Hessian shape を記録します。
- `legacy-active` は非推奨の比較専用処理で、pass/HOSP 遷移状態認定には
  使用できません。現行の共通射影 kernel で rank 退化構造も処理しますが、
  bitwise 一致は保証しません。

## 関連項目

- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細なトラブルシューティングガイド
- [tsopt](tsopt.md) -- IRC 実行前に TS を最適化
- [freq](freq.md) -- TS 候補が 1 つの虚振動数を持つことを検証; IRC 端点を解析
- [opt](opt.md) -- IRC 端点を真の極小に最適化
- [all](all.md) -- tsopt の後に IRC を実行する一気通貫ワークフロー
- [YAML リファレンス](yaml-reference.md) -- `irc` の完全な設定オプション
- [用語集](glossary.md) -- IRC（固有反応座標）の定義
