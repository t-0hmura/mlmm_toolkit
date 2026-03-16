# `freq`

## 概要

> **概要:** PHVA 対応の ML/MM 振動解析と熱化学（ZPE、Gibbs エネルギー等）を計算します。VRAM に余裕がある場合は `--hessian-calc-mode Analytical` でヘシアン評価を高速化できます。虚数振動数は負の値で表示されます。

`mlmm freq` は ML/MM 計算機（`mlmm.mlmm_calc.mlmm`）による振動解析を実行し、PHVA による凍結原子に対応します。基準振動アニメーションを `_trj.xyz` と `.pdb`（酵素の原子順序にマップバック）としてエクスポートし、オプションの `thermoanalysis` パッケージがインストールされている場合は Gaussian スタイルの熱化学サマリーを出力します。

## 最小例

```bash
mlmm freq -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_freq
```

## 出力の見方

- `result_freq/frequencies_cm-1.txt`
- `result_freq/mode_*_trj.xyz`
- `result_freq/mode_*.pdb`（PDB 入力の場合）

## よくある例

1. まずは出力モード数を絞って確認する。

```bash
mlmm freq -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --max-write 6 --out-dir ./result_freq_quick
```

2. 凍結原子を指定した PHVA と熱化学ダンプを実行する。

```bash
mlmm freq -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --freeze-atoms "1,3,5,7" --dump --out-dir ./result_freq_phva
```

3. VRAM に余裕があるノードで解析的ヘシアンを使う。

```bash
mlmm freq -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --hessian-calc-mode Analytical --out-dir ./result_freq_analytical
```

## ワークフロー

1. **ML/MM 計算機の構築** -- ML 領域は `--model-pdb` で提供され、Amber パラメータは `--parm` から読み取られます。`--hessian-calc-mode` は解析的または有限差分のヘシアンを選択します。計算機は完全な 3N x 3N ヘシアンまたはアクティブ自由度のサブブロックを返す場合があります。
   - VRAM に余裕がある場合は `--hessian-calc-mode` を `Analytical` に設定することを強く推奨します。
2. **PHVA と TR 射影** -- 凍結原子がある場合、固有解析はアクティブ部分空間内で行われ、並進/回転モードがそこに射影されます。3N x 3N とアクティブブロックの両方のヘシアンが受け付けられ、振動数は cm^-1 で報告されます（負の値 = 虚数）。
3. **アクティブ自由度モード** -- `--active-dof-mode` は振動解析に含まれる原子を制御します: `all`（全原子）、`ml-only`（ML 層, B=0）、`partial`（ML + Hessian 対象 MM、デフォルト）、`unfrozen`（非凍結層、通常 B=0/10）。
4. **モードエクスポート** -- `--max-write` はアニメーション化するモード数を制限します。モードは値（または `--sort abs` で絶対値）でソートされます。エクスポートされた各モードは `_trj.xyz`（XYZ ライク軌跡）と `.pdb` ファイル（酵素の原子順序にマップバックされた PDB アニメーション）を書き出します。正弦波アニメーション振幅（`--amplitude-ang`）とフレーム数（`--n-frames`）は YAML デフォルトに合致します。
5. **熱化学** -- `thermoanalysis` がインストールされている場合、PHVA 振動数を使用した QRRHO ライクなサマリー（EE、ZPE、E/H/G 補正、熱容量、エントロピー）が出力されます。CLI の圧力（atm）は内部で Pa に変換されます。`--dump` の場合、`thermoanalysis.yaml` スナップショットも書き出されます。
6. **デバイス選択** -- `ml_device="auto"` は CUDA が利用可能な場合にトリガーし、それ以外は CPU。内部の TR 射影/モード組み立ては転送を抑えるため同じデバイスで実行されます。
7. **終了動作** -- キーボード割り込みはコード 130 で終了。その他の失敗はトレースバックを出力してコード 1 で終了。

## CLIオプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-b, --backend CHOICE` | ML バックエンド: `uma`（デフォルト）、`orb`、`mace`、`aimnet2`。 | `uma` |
| `--embedcharge/--no-embedcharge` | xTB 点電荷埋め込み補正の有効化。MM 環境から ML 領域への静電的影響を考慮。 | `False` |
| `--embedcharge-cutoff FLOAT` | xTB 埋め込み用 MM 原子のカットオフ半径（Å）。 | `12.0` |
| `-i, --input PATH` | 完全酵素 PDB（リンク原子なし）。 | 必須 |
| `--parm PATH` | 完全酵素の Amber parm7 トポロジー。 | 必須 |
| `--model-pdb PATH` | ML 領域を定義する PDB。`--detect-layer` 有効時はオプション。 | _None_ |
| `--model-indices TEXT` | 明示的な ML 領域原子インデックス（`--model-pdb` の代替）。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` のインデックス規約。 | `True`（1 始まり） |
| `--detect-layer / --no-detect-layer` | 入力 PDB の B 因子から ML/MM レイヤーを自動検出。 | `True` |
| `-q, --charge INT` | ML 領域の電荷。 | _None_（`-l` 未指定時は必須） |
| `-l, --ligand-charge TEXT` | 残基ごとの電荷マッピング（例: `GPP:-3,SAM:1`）。`-q` 省略時に合計電荷を導出。 | _None_ |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `--freeze-atoms TEXT` | 1 始まりカンマ区切りの凍結原子インデックス。 | _None_ |
| `--hess-cutoff FLOAT` | Hessian 対象 MM 原子のカットオフ距離。 | _None_ |
| `--movable-cutoff FLOAT` | Movable-MM レイヤーのカットオフ距離。 | _None_ |
| `--hessian-calc-mode CHOICE` | ヘシアンモード（`Analytical` または `FiniteDifference`）。 | `FiniteDifference` |
| `--max-write INT` | エクスポートするモード数。 | `10` |
| `--amplitude-ang FLOAT` | モードアニメーション振幅 (Å)。 | `0.8` |
| `--n-frames INT` | モードアニメーションのフレーム数。 | `20` |
| `--sort CHOICE` | モードのソート方法: `value`（cm^-1）または `abs`。 | `value` |
| `--temperature FLOAT` | 熱化学温度 (K)。 | `298.15` |
| `--pressure FLOAT` | 熱化学圧力 (atm)。 | `1.0` |
| `--dump/--no-dump` | `thermoanalysis.yaml` を書き出し。 | `False` |
| `--convert-files/--no-convert-files` | PDB テンプレートが利用可能な場合の XYZ/TRJ から PDB コンパニオンの切り替え。 | `True` |
| `-o, --out-dir TEXT` | 出力ディレクトリ。 | `./result_freq/` |
| `--active-dof-mode CHOICE` | アクティブ自由度選択: `all`、`ml-only`、`partial`、`unfrozen`。 | `partial` |
| `--hess-device CHOICE` | ヘシアン組み立て/対角化のデバイス: `auto`、`cuda`、`cpu`。大規模系で VRAM 不足を回避するには `cpu` を使用。 | `auto` |
| `--ref-pdb FILE` | 非 PDB 入力用の参照 PDB トポロジー。 | _None_ |
| `--config FILE` | 明示 CLI 適用前に読み込むベース YAML。 | _None_ |
| `--show-config/--no-show-config` | 解決済み YAML レイヤー/設定を表示して続行。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画のみ表示。`--help-advanced` に表示。 | `False` |

## 出力

```
out_dir/ (デフォルト: ./result_freq/)
├─ mode_XXXX_±freqcm-1_trj.xyz   # モードごとの正弦波アニメーション（XYZ ライク軌跡）
├─ mode_XXXX_±freqcm-1.pdb       # 酵素原子順序にマップバックされた PDB アニメーション
├─ frequencies_cm-1.txt           # 選択されたソート順での全振動数リスト
└─ thermoanalysis.yaml            # thermoanalysis がインポート可能で --dump が True の場合
```

- コンソールには解決済みの `geom`、`calc`、`freq`、熱化学設定をまとめたブロックが出力されます。

## YAML設定

マージ順 **デフォルト < config < 明示CLI < override** でマッピングを提供します。
共有セクションは [YAML リファレンス](yaml_reference.md) を再利用します。
熱化学制御用の追加 `thermo` セクションがサポートされます。

```yaml
geom:
 coord_type: cart                  # 座標タイプ: デカルト vs dlc 内部座標
 freeze_atoms: []                  # 0 始まり凍結原子（CLI/リンク検出とマージ）
calc:
 charge: 0                         # 総電荷（CLI 上書き）
 spin: 1                           # スピン多重度 2S+1
mlmm:
 real_parm7: real.parm7            # Amber parm7 トポロジー
 model_pdb: ml_region.pdb          # ML 領域定義
 backend: uma                      # ML バックエンド (uma/orb/mace/aimnet2)
 embedcharge: false                # xTB 点電荷埋め込み補正
 uma_model: uma-s-1p1              # uma-s-1p1 | uma-m-1p1
 uma_task_name: omol                # UMA タスク名 (backend=uma 時)
 ml_device: auto                   # ML デバイス選択
 ml_hessian_mode: Analytical         # ヘシアンモード選択
 out_hess_torch: true              # torch 形式ヘシアンを要求
 mm_fd: true                       # MM 有限差分トグル
 return_partial_hessian: true      # 部分ヘシアンを許可（PHVA デフォルト）
freq:
 amplitude_ang: 0.8                # モードの変位振幅 (Å)
 n_frames: 20                      # モードごとのフレーム数
 max_write: 10                     # 書き出す最大モード数
 sort: value                       # ソート順: value vs abs
thermo:
 temperature: 298.15               # 熱化学温度 (K)
 pressure_atm: 1.0                 # 熱化学圧力 (atm)
 dump: false                       # true の場合 thermoanalysis.yaml を書き出し
```

---

## 関連項目

- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- 詳細なトラブルシューティングガイド

- [tsopt](tsopt.md) -- TS 候補の最適化（freq/IRC で検証; 期待: 1 つの虚数振動数）
- [opt](opt.md) -- 構造最適化（多くの場合 freq の前に実行）
- [dft](dft.md) -- より高レベルのエネルギー精密化のための DFT 一点計算
- [all](all.md) -- `--thermo` 付きエンドツーエンドワークフロー
- [YAML リファレンス](yaml_reference.md) -- `freq` と `thermo` の完全な設定オプション
- [用語集](glossary.md) -- ZPE、Gibbs エネルギー、エンタルピー、エントロピーの定義
