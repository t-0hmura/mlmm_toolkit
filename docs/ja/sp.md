# `sp`

`mlmm sp` は、単一構造における ML/MM ONIOM エネルギーと原子に働く力（任意で ONIOM Hessian 全体）を評価します。次のような用途に使います。

- 最適化を実行する前に層構造を高速に確認する
- 同一の ONIOM 分割上でバックエンドどうしを直接比較する
- オプティマイザのループ外で参照用 Hessian を生成する

## 実行例

層構造 PDB 上のエネルギーと力（B-factor が ML / movable-MM / frozen をエンコード）:

```bash
mlmm sp -i layered.pdb --parm real.parm7 -q 0 -m 1
```

ONIOM Hessian 全体も計算する（デフォルトは FiniteDifference。バックエンドのネイティブ Hessian を使うには `--hessian-calc-mode Analytical` を指定）:

```bash
mlmm sp -i layered.pdb --parm real.parm7 -q 0 -m 1 --hess
```

## 出力

`sp` はデフォルトで `result_sp/` 以下に出力を書き込みます。ONIOM エネルギーは stdout にも出力されます。JSON ファイル（同一内容を両方のファイル名（result.json / summary.json）に出力）は `--out-json` を指定したときのみ出力されます。

| ファイル | 内容 | 出力 |
|---|---|---|
| `forces.npy` | 原子単位（Hartree / Bohr）の ONIOM 力の `(N, 3)` 配列 | 常時 |
| `hessian.npy` | 質量で重み付けしていない `(3N, 3N)` ONIOM Hessian（Hartree / Bohr²） | `--hess` 指定時のみ |
| `result.json` / `summary.json` | ONIOM エネルギー（a.u.）、バックエンド、電荷/スピン、npy 出力へのパス、経過時間 | `--out-json` 指定時のみ |

`sp` は `summary.log` を書き込みません。

## CLI オプション

コマンド形式:

```bash
mlmm sp -i INPUT --parm PARM7 -q CHARGE [options]
```

| 入力 | 必須 | 備考 |
|---|---|---|
| `-i, --input FILE` | はい | 層構造 PDB/mmCIF、または `--ref-pdb` を伴う XYZ 座標 |
| `--ref-pdb FILE` | XYZ の場合 | 原子順序が一致する全系 PDB/mmCIF（トポロジーと層情報を供給） |
| `--parm FILE` | はい | 全系の Amber `parm7` トポロジー（`--real-parm7` をエイリアスとして保持） |
| `-q, --charge INT` | はい（`-l` を指定する場合は不要） | ML 領域の総電荷 |
| `-l, --ligand-charge TEXT` | いいえ | リガンドごとの電荷マッピング（例: `SAM:1,GPP:-3`）。`-q` を省略した場合に正味電荷を導出 |
| `-m, --multiplicity INT` | いいえ | ML 領域のスピン多重度、2S+1（デフォルト `1`） |

### ML 領域の選択

分割を入力 PDB の B-factor に埋め込む（ML=0.0、movable-MM=10.0、frozen=20.0）方法を `--detect-layer`（デフォルト）で使うか、明示的に渡します:

| フラグ | 意味 |
|---|---|
| `--detect-layer / --no-detect-layer` | B-factor エンコードを使用（デフォルト `on`） |
| `--model-pdb FILE` | ML 原子を定義する代替 PDB |
| `--model-indices TEXT` | カンマ区切りの 1-based 原子インデックス（例: `1-50,75,100-110`） |

### Hessian バックエンド

`--hess` と `--hessian-calc-mode Analytical` を指定すると、選択した
バックエンド（UMA、ORB、MACE、AIMNet2）の解析/native Hessian 経路を使います。
`FiniteDifference` は全バックエンドで利用できます。MM バックエンドは
`hessian_ff` がデフォルトですが、MM Hessian はデフォルトでは有限差分です。
`calc.mm_fd: false` で `hessian_ff` の解析 MM Hessian を選択できます。
要求した backend API が無い場合はエラーになります。

### その他のオプション

フラグの完全な一覧は自動生成された[コマンドリファレンス](../reference/commands/index.md)にあります。以下の表は説明が必要なオプションを扱います。

| フラグ | デフォルト | 意味 |
|---|---|---|
| `-b, --backend [uma\|orb\|mace\|aimnet2]` | `uma` | ML 領域の MLIP バックエンド |
| `--hess / --no-hess` | `--no-hess` | `hessian.npy` も計算して書き込む |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | `FiniteDifference` | `--hess` 指定時の Hessian モード。`Analytical` はバックエンドのネイティブ経路を使用 |
| `--embedcharge / --no-embedcharge` | off | v0.3.3 では使用不可。旧コマンドの明示的拒否用 |
| `--link-atom-method [scaled\|fixed]` | `scaled` | リンク原子の配置 |
| `--mm-backend [hessian_ff\|openmm]` | `hessian_ff` | MM バックエンド。Hessian 法は `calc.mm_fd` で別に選択 |
| `-o, --out-dir PATH` | `./result_sp/` | 出力ディレクトリ |
| `--precision [fp32\|fp64]` | バックエンド依存 | バックエンドに渡す数値精度（未指定: UMA/AIMNet2 は fp32、ORB/MACE は fp64） |
| `--config PATH` | — | `calc.*`、`geom.*` のデフォルトを与える YAML 設定 |
| `--show-config / --dry-run` | off | 有効なマージ済み設定を表示 / 実行せずに検証 |

Hessian の cutoff 上書き、MCP 形式の result.json などを含む完全な一覧は `mlmm sp --help-advanced` を実行してください。

## 関連項目

- [`opt`](opt.md) — 層構造を最適化（マイクロイテレーション）
- [`tsopt`](tsopt.md) — TS 候補を精密化（ML/MM ONIOM）
- [`freq`](freq.md) — ONIOM 振動解析 + QRRHO 熱化学
- [`dft`](dft.md) — ML 領域上の DFT 一点計算に相当する処理
