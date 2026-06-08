# `sp`

`mlmm sp` は、単一構造における ML/MM ONIOM エネルギーと原子に働く力（任意で ONIOM ヘシアン全体）を評価します。最適化を実行する前に層構造を高速に確認する、同一の ONIOM 分割上でバックエンドどうしを直接比較する、オプティマイザのループ外で参照用ヘシアンを生成する、といった用途に使います。

## 実行例

層構造 PDB 上のエネルギーと力（B-factor が ML / movable-MM / frozen-MM をエンコード）:

```bash
# energy + forces on a layered PDB (B-factor encodes ML / movable-MM / frozen-MM)
mlmm sp -i layered.pdb --parm real.parm7 -q 0 -m 1
```

ONIOM ヘシアン全体も計算する（`--backend uma` のとき Analytical）:

```bash
# also compute the full ONIOM Hessian (Analytical when --backend uma)
mlmm sp -i layered.pdb --parm real.parm7 -q 0 -m 1 --hess
```

## 出力

`sp` はデフォルトで `result_sp/` 以下に出力を書き込みます。ONIOM エネルギーは stdout にも出力されます。JSON ファイル（同一のペイロードを両方の名前にミラー）は `--out-json` を指定したときのみ出力されます。

| ファイル | 内容 | 出力 |
|---|---|---|
| `forces.npy` | 原子単位（Hartree / Bohr）の ONIOM 力の `(N, 3)` 配列 | 常時 |
| `hessian.npy` | 質量で重み付けしていない `(3N, 3N)` ONIOM ヘシアン（Hartree / Bohr²） | `--hess` 指定時のみ |
| `result.json` / `summary.json` | ONIOM エネルギー（a.u.）、バックエンド、電荷/スピン、npy 出力へのパス、経過時間 | `--out-json` 指定時のみ |

`sp` は `summary.log` を書き込みません。

## CLI オプション

コマンド形式:

```bash
mlmm sp -i INPUT --parm PARM7 -q CHARGE [options]
```

| 入力 | 必須 | 備考 |
|---|---|---|
| `-i, --input FILE` | はい | ML / movable-MM / frozen-MM 分割を定義する層構造 PDB（または XYZ） |
| `--parm FILE` | はい | 全系の Amber `parm7` トポロジー（`--real-parm7` をエイリアスとして保持） |
| `-q, --charge INT` | はい | ML 領域の総電荷 |
| `-l, --ligand-charge TEXT` | いいえ | リガンドごとの電荷マッピング（例: `SAM:1,GPP:-3`） |
| `-m, --multiplicity INT` | いいえ | ML 領域のスピン多重度、2S+1（デフォルト `1`） |

### ML 領域の選択

分割を入力 PDB の B-factor に埋め込む（ML=0.0、movable-MM=10.0、frozen-MM=20.0）方法を `--detect-layer`（デフォルト）で使うか、明示的に渡します:

| フラグ | 意味 |
|---|---|
| `--detect-layer / --no-detect-layer` | B-factor エンコードを使用（デフォルト `on`） |
| `--model-pdb FILE` | ML 原子を定義する代替 PDB |
| `--model-indices TEXT` | カンマ区切りの 1-based 原子インデックス（例: `1-50,75,100-110`） |

### ヘシアンバックエンド

`--hess` を指定すると、バックエンドの選択がヘシアン計算戦略を決めます:

- `--backend uma`（デフォルト）→ UMA の torch autograd 経路による ML 領域の `Analytical` ヘシアン。MM 領域は `hessian_ff` の解析ヘシアンを使用
- `--backend orb` / `mace` / `aimnet2` → ML 領域は `FiniteDifference` にフォールバック

`--hessian-calc-mode` で呼び出しごとに上書きできます。

### その他のオプション

フラグの完全な一覧は自動生成された[コマンドリファレンス](../reference/commands/index.md)にあります。以下の表は説明が必要なオプションを扱います。

| フラグ | デフォルト | 意味 |
|---|---|---|
| `-b, --backend [uma\|orb\|mace\|aimnet2]` | `uma` | ML 領域の MLIP バックエンド |
| `--hess / --no-hess` | `--no-hess` | `hessian.npy` も計算して書き込む |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | auto | 特定のヘシアンモードを強制（`--hess` 指定時のみ） |
| `--embedcharge / --no-embedcharge` | off | MM→ML カップリングのための xTB 点電荷埋め込み補正 |
| `--link-atom-method [scaled\|fixed]` | `scaled` | リンク原子の配置 |
| `--mm-backend [hessian_ff\|openmm]` | `hessian_ff` | MM バックエンド（解析ヘシアン vs 有限差分ヘシアン） |
| `-o, --out-dir PATH` | `./result_sp/` | 出力ディレクトリ |
| `--precision [fp32\|fp64]` | `fp32` | バックエンドに渡す数値精度 |
| `--config PATH` | — | `calc.*`、`geom.*` のデフォルトを与える YAML 設定 |
| `--show-config / --dry-run` | off | 有効なマージ済み設定を表示 / 実行せずに検証 |

ヘシアンの cutoff 上書き、MCP 形式の result.json などを含む完全な一覧は `mlmm sp --help-advanced` を実行してください。

## 関連項目

- [`opt`](opt.md) — 層構造を最適化（マイクロイテレーション）
- [`tsopt`](tsopt.md) — TS 候補を精密化（ML/MM ONIOM）
- [`freq`](freq.md) — ONIOM 振動解析 + QRRHO 熱化学
- [`dft`](dft.md) — ML 領域上の DFT 一点計算に相当する処理
