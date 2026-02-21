# `freq`

## 概要

> **要約:** PHVA 対応の ML/MM 振動解析と熱化学（ZPE、Gibbs エネルギー等）を計算します。虚数振動数は負の値で表示されます。適切に収束した遷移状態は正確に 1 つの虚数振動数を持つべきです。

### 概要
- **用途:** ML/MM による極小/TS 候補の検証および熱力学補正の計算。
- **凍結原子:** PHVA（部分ヘシアン振動解析）でサポート。
- **出力:** `frequencies_cm-1.txt`、モードごとの `.trj` と `.pdb` アニメーション、有効時は `thermoanalysis.yaml`。
- **TS チェック:** 適切に収束した TS は**正確に 1 つ**の虚数振動数（負の cm^-1）を持つことが期待されます。

`mlmm freq` は ML/MM 計算機（`mlmm_toolkit.mlmm_calc.mlmm`）による振動解析を実行し、PHVA による凍結原子に対応します。基準振動アニメーションを `.trj` と `.pdb`（酵素の原子順序にマップバック）としてエクスポートし、オプションの `thermoanalysis` パッケージがインストールされている場合は Gaussian スタイルの熱化学サマリーを出力します。

設定は **デフォルト < `--config` < 明示CLI < `--override-yaml`** に従います（`geom`、`calc`/`mlmm`、`freq`、`thermo`）。`--args-yaml` は `--override-yaml` の legacy alias として引き続き利用できます。ML 領域は `--model-pdb` で、完全酵素トポロジーは `--real-parm7` で提供され、入力座標は完全酵素 PDB（リンク原子なし）でなければなりません。

## 最小例

```bash
mlmm freq -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --out-dir ./result_freq
```

## 出力の見方

- `result_freq/summary.md`
- `result_freq/key_frequencies.txt`
- `result_freq/key_mode_1.trj`
- `result_freq/frequencies_cm-1.txt`
- `result_freq/mode_*.trj`
- `result_freq/mode_*.pdb`（PDB 入力の場合）

## よくある例

1. まずは出力モード数を絞って確認する。

```bash
mlmm freq -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --max-write 6 --out-dir ./result_freq_quick
```

2. 凍結原子を指定した PHVA と熱化学ダンプを実行する。

```bash
mlmm freq -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --freeze-atoms "1,3,5,7" --dump --out-dir ./result_freq_phva
```

3. VRAM に余裕があるノードで解析的ヘシアンを使う。

```bash
mlmm freq -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --hessian-calc-mode Analytical --out-dir ./result_freq_analytical
```

## 使用法

```bash
mlmm freq -i INPUT.pdb --real-parm7 real.parm7 --model-pdb model.pdb \
    -q CHARGE [-m MULT] [--freeze-atoms "1,3,5"] \
    [--max-write N] [--amplitude-ang FLOAT] [--n-frames N] [--sort {value|abs}] \
    [--temperature K] [--pressure atm] [--dump/--no-dump] \
    [--hessian-calc-mode {Analytical|FiniteDifference}] \
    [--active-dof-mode {all|ml-only|partial|unfrozen}] \
    [--out-dir DIR] [--config FILE] [--override-yaml FILE|--args-yaml FILE] \
    [--show-config] [--dry-run]
```

### 例

```bash
# ミニマル実行
mlmm freq -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0

# カスタムオプション付き PHVA
mlmm freq -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 \
    --freeze-atoms "1,3,5,7" --max-write 10 --sort abs --dump --config ./freq.yaml
```

## ワークフロー

- **ジオメトリ読み込みと凍結処理**: 構造は `pysisyphus.helpers.geom_loader` で読み込まれます。`--freeze-atoms "1,3,5"` は 1 始まりインデックスを受け付け、YAML `geom.freeze_atoms` とマージされます。マージされたリストはジオメトリエコーと ML/MM 計算機の両方に渡され、PHVA が有効になります。
- **ML/MM 計算機**: ML 領域は `--model-pdb` で提供され、Amber パラメータは `--real-parm7` から読み取られます。`--hessian-calc-mode` は解析的または有限差分のヘシアンを選択します。計算機は完全な 3N x 3N ヘシアンまたはアクティブ自由度のサブブロックを返す場合があります。
- **PHVA と TR 射影**: 凍結原子がある場合、固有解析はアクティブ部分空間内で行われ、並進/回転モードがそこに射影されます。3N x 3N とアクティブブロックの両方のヘシアンが受け付けられ、振動数は cm^-1 で報告されます（負の値 = 虚数）。
- **アクティブ自由度モード**: `--active-dof-mode` は振動解析に含まれる原子を制御します: `all`（全原子）、`ml-only`（ML 層, B=0）、`partial`（ML + Hessian 対象 MM、デフォルト）、`unfrozen`（非凍結層、通常 B=0/10）。
- **モードエクスポート**: `--max-write` はアニメーション化するモード数を制限します。モードは値（または `--sort abs` で絶対値）でソートされます。エクスポートされた各モードは `.trj`（XYZ ライク軌跡）と `.pdb` ファイル（酵素の原子順序にマップバックされた PDB アニメーション）を書き出します。正弦波アニメーション振幅（`--amplitude-ang`）とフレーム数（`--n-frames`）は YAML デフォルトに合致します。
- **熱化学**: `thermoanalysis` がインストールされている場合、PHVA 振動数を使用した QRRHO ライクなサマリー（EE、ZPE、E/H/G 補正、熱容量、エントロピー）が出力されます。CLI の圧力（atm）は内部で Pa に変換されます。`--dump` の場合、`thermoanalysis.yaml` スナップショットも書き出されます。
- **デバイス選択**: `ml_device="auto"` は CUDA が利用可能な場合にトリガーし、それ以外は CPU。内部の TR 射影/モード組み立ては転送を最小化するため同じデバイスで実行されます。
- **終了動作**: キーボード割り込みはコード 130 で終了。その他の失敗はトレースバックを出力してコード 1 で終了。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input PATH` | 完全酵素 PDB（リンク原子なし）。 | 必須 |
| `--real-parm7 PATH` | 完全酵素の Amber parm7 トポロジー。 | 必須 |
| `--model-pdb PATH` | ML 領域を定義する PDB。 | 必須 |
| `--model-indices TEXT` | 明示的な ML 領域原子インデックス（`--model-pdb` の代替）。 | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | `--model-indices` のインデックス規約。 | `True`（1 始まり） |
| `--detect-layer / --no-detect-layer` | 入力 PDB の B 因子から ML/MM レイヤーを自動検出。 | `True` |
| `-q, --charge INT` | ML 領域の電荷。 | 必須 |
| `-m, --multiplicity INT` | スピン多重度 (2S+1)。 | `1` |
| `--freeze-atoms TEXT` | 1 始まりカンマ区切りの凍結原子インデックス。 | _None_ |
| `--hess-cutoff FLOAT` | Hessian 対象 MM 原子のカットオフ距離。 | _None_ |
| `--movable-cutoff FLOAT` | Movable-MM レイヤーのカットオフ距離。 | _None_ |
| `--hessian-calc-mode CHOICE` | ヘシアンモード（`Analytical` または `FiniteDifference`）。 | _None_ |
| `--max-write INT` | エクスポートするモード数。 | `20` |
| `--amplitude-ang FLOAT` | モードアニメーション振幅 (A)。 | `0.8` |
| `--n-frames INT` | モードアニメーションのフレーム数。 | `20` |
| `--sort CHOICE` | モードのソート方法: `value`（cm^-1）または `abs`。 | `value` |
| `--temperature FLOAT` | 熱化学温度 (K)。 | `298.15` |
| `--pressure FLOAT` | 熱化学圧力 (atm)。 | `1.0` |
| `--dump/--no-dump` | `thermoanalysis.yaml` を書き出し。 | `False` |
| `--out-dir TEXT` | 出力ディレクトリ。 | `./result_freq/` |
| `--active-dof-mode CHOICE` | アクティブ自由度選択: `all`、`ml-only`、`partial`、`unfrozen`。 | `partial` |
| `--config FILE` | 明示CLI適用前に読み込むベース YAML。 | _None_ |
| `--override-yaml FILE` | 最終 YAML 上書き（YAML レイヤー最優先）。 | _None_ |
| `--args-yaml FILE` | `--override-yaml` の legacy alias。 | _None_ |
| `--show-config/--no-show-config` | 解決済み YAML レイヤー/設定を表示して続行。 | `False` |
| `--dry-run/--no-dry-run` | 実行せずに検証と実行計画のみ表示。 | `False` |
| `--ref-pdb FILE` | 非 PDB 入力用の参照 PDB トポロジー。 | _None_ |

## 出力

```
out_dir/ (デフォルト: ./result_freq/)
  summary.md                      # 主要成果物のインデックス
  key_frequencies.txt             # frequencies_cm-1.txt へのショートカット
  key_mode_1.trj                  # 代表モード軌跡へのショートカット
  key_mode_1.pdb                  # 代表モードPDBへのショートカット（存在時）
  key_thermo.yaml                 # thermoanalysis.yaml へのショートカット（存在時）
  mode_XXXX_{+/-freq}cm-1.trj    # XYZ ライク軌跡、モードごとの正弦波アニメーション
  mode_XXXX_{+/-freq}cm-1.pdb    # 酵素原子順序にマップバックされた PDB アニメーション
  frequencies_cm-1.txt            # 選択されたキーでソートされた全計算振動数（cm^-1）
  thermoanalysis.yaml             # --dump 時のオプション熱化学ペイロード
```

- コンソールには解決済みの `geom`、`calc`、`freq`、熱化学設定をまとめたブロックが出力されます。

## 注意事項

- 症状起点で切り分ける場合は [典型エラー別レシピ](recipes-common-errors.md) を先に参照し、詳細は [トラブルシューティング](troubleshooting.md) を確認してください。

- ML/MM 計算機は Hartree/Bohr^2 でヘシアンを返します。cm^-1 への変換はツールキットの他の箇所で使用される PySisyphus/ASE 規約に従います。
- 熱化学はオプションの `thermoanalysis` パッケージに依存します。未インストール時は警告のみ出力され、実行は継続します。
- `--hessian-calc-mode` は標準的な優先順位（デフォルト < config < 明示CLI < override）に従います。

---

## 関連項目

- [tsopt](tsopt.md) -- TS 候補の最適化（freq/IRC で検証; 期待: 1 つの虚数振動数）
- [opt](opt.md) -- 構造最適化（多くの場合 freq の前に実行）
- [dft](dft.md) -- より高レベルのエネルギー精密化のための DFT 一点計算
- [all](all.md) -- `--thermo` 付きエンドツーエンドワークフロー
- [YAML リファレンス](yaml-reference.md) -- `freq` と `thermo` の完全な設定オプション
