# CLI 規約

このページでは、`mlmm` の全コマンドで使用される規約を説明します。これらの規約を理解することで、正しいコマンドを記述し、よくあるエラーを回避できます。

---

## ブール値オプション

ブール値オプションは root CLI で正規化されます。
基本は toggle 形式（`--flag/--no-flag`）を使ってください。

```bash
# 推奨
--tsopt --thermo --no-dft
--dump                 # 旧記法の --dump True と同義
--no-dump              # 旧記法の --dump False と同義

# 後方互換として受理（非推奨）
--tsopt True --thermo yes --dft 0
```

旧来の値指定形式（`--flag True/False`）は当面受理されますが、deprecation warning が出ます。
旧記法の値としては `True/False`, `1/0`, `yes/no`, `y/n`, `t/f`, `on/off`（大文字小文字非区別）を受け入れます。

よく使うブール値オプション:
- `--tsopt`, `--thermo`, `--dft` -- 後処理ステージの有効化
- `--dump` -- 軌跡ファイルの出力
- `--preopt`, `--endopt` -- 前処理/後処理最適化の切り替え
- `--climb` -- MEP 探索でクライミングイメージを有効化

---

## Progressive Help (`all`)

`mlmm all` は 2 段階ヘルプです:

```bash
mlmm all --help            # 主要オプションのみ
mlmm all --help-advanced   # 全オプション
```

`scan` / `scan2d` / `scan3d` と計算系サブコマンド（`opt` / `path-opt` / `path-search` / `tsopt` / `freq` / `irc` / `dft`）に加え、ユーティリティ系（`mm-parm` / `define-layer` / `add-elem-info` / `trj2fig` / `energy-diagram` / `oniom-gaussian` / `oniom-orca`）も同様に `--help` は主要オプションのみ、`--help-advanced` で全オプションを表示します。`extract` と `fix-altloc` も段階的 help に対応し、`--help-advanced` で parser の全オプションを表示します。

---

## init テンプレート

`all` 用の YAML テンプレートを生成して検証できます:

```bash
mlmm init --out mlmm_all.config.yaml
mlmm all --config mlmm_all.config.yaml --dry-run
```

---

## 残基セレクタ

残基セレクタは、基質や抽出中心として使用する残基を指定します。

### 残基名による指定
```bash
-c 'SAM,GPP'          # SAM または GPP という名前の残基をすべて選択
-c 'LIG'              # LIG という名前の残基をすべて選択
```

### 残基IDによる指定
```bash
-c '123,456'          # 残基 123 と 456
-c 'A:123,B:456'      # チェーン A の残基 123、チェーン B の残基 456
-c '123A'             # 挿入コード A を持つ残基 123
-c 'A:123A'           # チェーン A、残基 123、挿入コード A
```

### PDB ファイルによる指定
```bash
-c substrate.pdb      # 別の PDB から座標を使用して基質を特定
```

```{note}
残基名で選択する場合、同名の残基が複数あれば**すべて**が含まれ、警告がログに出力されます。
```

---

## 電荷の指定

### 残基別マッピング（推奨）

`--ligand-charge` はコロン（`:`）またはイコール（`=`）の両方をサポートします:

```bash
--ligand-charge 'SAM:1,GPP:-3'    # SAM は +1、GPP は -3
--ligand-charge 'SAM=1,GPP=-3'    # 同じ意味（= 区切り）
--ligand-charge 'LIG:-2'          # LIG は -2
```

### 総電荷の明示的上書き
```bash
-q 0                              # ML 領域の総電荷を 0 に強制
-q -1                             # ML 領域の総電荷を -1 に強制
```

### 電荷の解決順序
1. `-q/--charge`（明示的な CLI 上書き）-- 最優先
2. ポケット抽出（アミノ酸、イオン、`--ligand-charge` の合計）
3. フォールバックとしての `--ligand-charge`（抽出がスキップされた場合）
4. `.gjf` テンプレートのメタデータ
5. デフォルト: なし（未解決なら中断）

```{tip}
非標準の残基（基質、補因子、特殊なリガンド）には必ず `--ligand-charge` を指定し、正しい電荷伝播を確保してください。
```

---

## スピン多重度

```bash
-m 1      # 一重項（デフォルト）
-m 2      # 二重項
-m 3      # 三重項
```

```{note}
`all` と計算系サブコマンドでは `-m/--multiplicity` を使います。`mm-parm` の `--ligand-mult` は残基メタデータ用の別オプションです。
```

---

## B-factor 層エンコーディング

`mlmm_toolkit` は PDB の B-factor カラム（列 61-66）を使用して、原子の層割り当てをエンコードします:

| B-factor | 層 | 説明 |
|----------|-----|------|
| **0.0** | Layer 1: ML | UMA による完全計算（エネルギー・力・ヘシアン） |
| **10.0** | Layer 2: Movable-MM | 最適化時に移動可能な MM 原子 |
| **20.0** | Layer 3: Frozen | 座標固定 |

### 層の定義方法

1. **`define-layer` サブコマンド**（推奨）:
   ```bash
   mlmm define-layer -i system.pdb --model-pdb ml_region.pdb -o labeled.pdb
   ```

2. **距離カットオフ**（YAML/CLI）:
   ```yaml
   calc:
     hess_cutoff: 3.6      # Hessian 対象 MM の距離カットオフ
     movable_cutoff: 8.0   # Movable-MM の距離カットオフ（それ以外は Frozen）
   ```

3. **B-factor からの読み取り**:
   ```yaml
   calc:
     use_bfactor_layers: true   # 入力 PDB の B-factor から層を読み取り
   ```

4. **明示的インデックス指定**（YAML）:
   ```yaml
   calc:
     hess_mm_atoms: [100, 101, 102, ...]
     movable_mm_atoms: [200, 201, 202, ...]
     frozen_mm_atoms: [300, 301, 302, ...]
   ```

---

## --real-parm7 と --model-pdb

ML/MM 計算を行う大半のサブコマンド（`opt`, `tsopt`, `path-opt`, `path-search`, `scan`, `freq`, `irc` など）では、以下の 2 つのオプションが必要です:

| オプション | 説明 |
|----------|------|
| `--real-parm7` | 全系（real system）の Amber parm7 トポロジファイル |
| `--model-pdb` | ML 領域（model system）を定義する PDB ファイル |

`all` ワークフローでは、`extract` と `mm-parm` により自動生成されます。個別サブコマンドを使う場合は、手動で指定してください。

```bash
# 個別サブコマンドの例
mlmm opt -i input.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0
```

---

## 原子セレクタ

原子セレクタは、スキャンや拘束に使用する特定の原子を指定します。指定方法は以下の通りです:

### 整数インデックス（デフォルトは1始まり）
```bash
--scan-lists '[(1, 5, 2.0)]'      # 原子 1 と 5、ターゲット距離 2.0 Å
```

### PDB形式のセレクタ文字列
```bash
--scan-lists '[("TYR,285,CA", "MMT,309,C10", 2.20)]'
```

セレクタのフィールドは以下で区切れます:
- 空白: `'TYR 285 CA'`
- カンマ: `'TYR,285,CA'`
- スラッシュ: `'TYR/285/CA'`
- バッククォート: `` 'TYR`285`CA' ``
- バックスラッシュ: `'TYR\285\CA'`

---

## 入力ファイル要件

### PDB ファイル
- **水素原子**を含む必要があります（`reduce`、`pdb2pqr`、`mm-parm --add-h` 等で追加）
- 列 77-78 に**元素記号**が必要（欠けている場合は `mlmm add-elem-info` を使用）
- 複数の PDB は**同じ原子を同じ順序**で持つ必要があります（座標のみ異なる）

### Amber トポロジ
- **parm7**: 全系の力場トポロジ（`mm-parm` で自動生成可能）
- **rst7**: 対応する座標ファイル

---

## YAML 設定

詳細設定は多層 YAML で渡せます:

```bash
mlmm all -i R.pdb P.pdb -c 'LIG' --config config.yaml --override-yaml override.yaml
```

適用順序:
```
デフォルト < config < CLI オプション < override-yaml
```

利用可能なすべてのオプションは [YAML リファレンス](yaml-reference.md) を参照してください。

---

## 出力ディレクトリ

`--out-dir` で結果の保存先を指定します:

```bash
--out-dir ./my_results/    # カスタム出力ディレクトリ
```

デフォルトの出力ディレクトリ:
- `all`: `./result_all/`
- `opt`: `./result_opt/`
- `tsopt`: `./result_tsopt/`
- `path-opt`: `./result_path_opt/`
- `path-search`: `./result_path_search/`
- `scan`: `./result_scan/`
- `freq`: `./result_freq/`
- `irc`: `./result_irc/`
- `dft`: `./result_dft/`

---

## 関連項目

- [はじめに](getting-started.md) -- インストールと初回実行
- [概念とワークフロー](concepts.md) -- ML/MM 3層システム、ONIOM 分解の全体像
- [典型エラー別レシピ](recipes-common-errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- よくあるエラーと対処法
- [YAML リファレンス](yaml-reference.md) -- 全設定オプション
