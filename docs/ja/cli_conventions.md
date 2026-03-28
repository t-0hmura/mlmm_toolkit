# CLI 規約

このページでは、`mlmm-toolkit` の全コマンドで使用される規約を説明します。これらの規約を理解することで、正しいコマンドを記述し、よくあるエラーを回避できます。

---

## ブール値オプション

ブール値オプションは root CLI で正規化されます。
次の2記法を受け付けます。

```bash
# 推奨
--tsopt --thermo --no-dft

# 互換記法として受理
--tsopt True --thermo yes --dft 0
```

`--flag` 単独で定義されているオプションでも、互換のため `--no-flag` と `--flag False` を受理します。
`extract` と `fix-altloc` を含むすべてのサブコマンドが Click を CLI バックエンドとして使用します。

よく使うブール値オプション:
- `--tsopt`, `--thermo`, `--dft` -- 後処理ステージの有効化
- `--dump` -- 軌跡ファイルの出力
- `--preopt`, `--endopt` -- 前処理/後処理最適化の切り替え
- `--climb` -- MEP 探索でクライミングイメージを有効化

---

## Progressive Help (`all`)

`mlmm all` は 2 段階ヘルプです:

```bash
mlmm all --help           # 主要オプションのみ
mlmm all --help-advanced  # 全オプション
```

`scan` / `scan2d` / `scan3d` と計算系サブコマンド（`opt` / `path-opt` / `path-search` / `tsopt` / `freq` / `irc` / `dft`）に加え、ユーティリティ系（`mm-parm` / `define-layer` / `add-elem-info` / `trj2fig` / `energy-diagram` / `oniom-export`）も同様に `--help` は主要オプションのみ、`--help-advanced` で全オプションを表示します。`extract` と `fix-altloc` も段階的 help に対応し、`--help-advanced` で parser の全オプションを表示します。

---

## 設定の確認

現在の設定を確認できます（YAML オーバーライドの検証に便利）：

```bash
mlmm opt -i input.pdb --parm real.parm7 -q -1 --show-config --dry-run
```

---

## ML/MM 必須オプション

ML/MM 計算を行う大半のサブコマンド（`all`、`extract`、`mm-parm`、`define-layer` を除く）では、以下の 2 つのオプションが必要です:

```bash
--parm real.parm7      # 全系（real system）の Amber parm7 トポロジーファイル
--model-pdb model.pdb  # ML 領域（model system）を定義する PDB ファイル
```

`all` ワークフローでは、`mm-parm` と `define-layer` により自動生成されます。個別サブコマンドを使う場合は、手動で指定してください。

```bash
# 個別サブコマンドの例
mlmm path-search -i R.pdb P.pdb --parm real.parm7 --model-pdb model.pdb -q 0 -m 1
```

---

## B-factor 層エンコーディング

`mlmm-toolkit` は PDB の B-factor カラム（列 61-66）を使用して、3 層 ML/MM 分割をエンコードします:

| 層 | B-factor | 説明 |
|-----|----------|------|
| ML | 0.0 | MLIP によるエネルギー・力・ヘシアン計算（デフォルトバックエンド: UMA） |
| Movable-MM | 10.0 | 最適化時に移動可能な MM 原子 |
| Frozen | 20.0 | 座標固定 |

`define-layer` サブコマンドがこれらの B-factor を PDB に書き込みます。B-factor カラーリングに対応した分子ビューアで層割り当てを確認できます。

B-factor の読み取り時には許容差 1.0 が使用され、0/10/20 に近い値はそれぞれ ML/Movable/Frozen にマッピングされます。

### 層の定義方法

1. **`define-layer` サブコマンド**（推奨）:
   ```bash
   mlmm define-layer -i system.pdb --model-pdb ml_region.pdb -o labeled.pdb
   ```

2. **距離カットオフ**（YAML/CLI）:
   ```yaml
   calc:
     hess_cutoff: 3.6       # Hessian 対象 MM の距離カットオフ
     movable_cutoff: 8.0    # Movable-MM の距離カットオフ（それ以外は Frozen）
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

## 残基セレクタ

残基セレクタは、基質や抽出中心として使用する残基を指定します。

### 残基名による指定
```bash
-c 'SAM,GPP'   # SAM または GPP という名前の残基をすべて選択
-c 'LIG'       # LIG という名前の残基をすべて選択
```

### 残基IDによる指定
```bash
-c '123,456'       # 残基 123 と 456
-c 'A:123,B:456'   # チェーン A の残基 123、チェーン B の残基 456
-c '123A'          # 挿入コード A を持つ残基 123
-c 'A:123A'        # チェーン A、残基 123、挿入コード A
```

### PDB ファイルによる指定
```bash
-c substrate.pdb   # 別の PDB から座標を使用して基質を特定
```

```{note}
残基名で選択する場合、同名の残基が複数あれば**すべて**が含まれ、警告がログに出力されます。
```

---

## 電荷の指定

PDB 入力の場合、`--ligand-charge` で非標準残基（基質、補因子、金属イオン）の電荷のみを指定します。総系電荷は、標準アミノ酸の電荷、イオン電荷、指定したリガンド電荷を合計して**自動算出**されるため、複合体全体の原子を手動で数える必要がありません。これは、総電荷が一目でわからない大規模な酵素-基質系で特に有用です。

### 残基別マッピング（推奨）

`--ligand-charge` はコロン（`:`）またはイコール（`=`）の両方をサポートします:

```bash
-l 'SAM:1,GPP:-3'   # SAM は +1、GPP は -3
-l 'SAM=1,GPP=-3'   # 同じ意味（= 区切り）
-l 'LIG:-2'         # LIG は -2
```

### 総電荷の明示的上書き
```bash
-q 0    # 総系電荷を 0 に強制
-q -1   # 総系電荷を -1 に強制
```

### 電荷の解決順序
1. `-q/--charge`（明示的な CLI 上書き）-- 最優先。
2. ML 領域決定の電荷サマリー（アミノ酸、イオン、`--ligand-charge` の合計）。
3. 抽出をスキップした場合の `--ligand-charge` フォールバック。
4. デフォルト: なし（未解決なら中断）。

計算系サブコマンド（`scan` / `scan2d` / `scan3d` / `opt` / `path-opt` / `path-search` / `tsopt` / `freq` / `irc` / `dft` / `oniom-export`）では、引き続き `-q/--charge` の明示指定が必要です。

```{tip}
非標準の残基（基質、補因子、特殊なリガンド）には必ず `--ligand-charge` を指定し、正しい電荷伝播を確保してください。
```

---

## リガンド電荷フォーマット

`--ligand-charge` オプションは 2 つのフォーマットをサポートします:

### マッピング形式（推奨）
```bash
-l 'SAM:1,GPP:-3'   # 残基名ごとのマッピング
-l 'SAM=1,GPP=-3'   # 同じ意味（= 区切り）
-l 'LIG:-2'         # 単一残基のマッピング
```

コロン（`:`）とイコール（`=`）の両方をセパレータとして使用できます。

### 整数形式
```bash
-l -3   # 全未知残基の合計電荷
```

マッピング形式では、残基名の大文字小文字は区別しません。マッピングされていない非標準残基はデフォルトで電荷 0 となります。

---

## スピン多重度

```bash
-m 1   # 一重項（デフォルト）
-m 2   # 二重項
-m 3   # 三重項
```

```{note}
`all` と計算系サブコマンドでは `-m/--multiplicity` を使います。`mm-parm` の `--ligand-mult` は残基メタデータ用の別オプションです。
```

---

## 原子セレクタ

原子セレクタは、スキャンや拘束に使用する特定の原子を指定します。指定方法は以下の通りです:

### 整数インデックス（デフォルトは1始まり）
```bash
--scan-lists '[(1, 5, 2.0)]'   # 原子 1 と 5、ターゲット距離 2.0 Å
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

3 つのトークン（残基名、残基番号、原子名）は任意の順序で指定できます。パーサーは非標準の順序でもフォールバックヒューリスティックで解釈します。

---

## 入力ファイル要件

### PDB ファイル
- **水素原子**を含む必要があります（`reduce`、`pdb2pqr`、`mm-parm --add-h` 等で追加）
- 列 77-78 に**元素記号**が必要（欠けている場合は `mlmm add-elem-info` を使用）
- 複数の PDB は**同じ原子を同じ順序**で持つ必要があります（座標のみ異なる）

### XYZ ファイル
- ML 領域決定をスキップする場合（`-c/--center` を省略）に使用可能

### Amber トポロジー
- `--parm`: 全系の力場トポロジー（`mm-parm` で自動生成可能）
- parm7 は入力 PDB の原子順序と正確に一致する必要があります

---

## バックエンド選択

すべての計算系サブコマンド（`opt`、`tsopt`、`freq`、`irc`、`dft`、`scan`、`scan2d`、`scan3d`、`path-opt`、`path-search`、`all`）で以下を指定できます:

| オプション | 説明 | デフォルト |
|----------|------|----------|
| `-b, --backend` | ML 領域の MLIP バックエンド: `uma`、`orb`、`mace`、`aimnet2` | `uma` |
| `--embedcharge/--no-embedcharge` | xTB 点電荷埋め込み補正の有効化 | `--no-embedcharge` |

代替バックエンドはオプション依存グループでインストールします:

```bash
pip install "mlmm-toolkit[orb]"       # ORB バックエンド
pip install "mlmm-toolkit[aimnet2]"   # AIMNet2 バックエンド
# MACE: pip uninstall fairchem-core && pip install mace-torch (別環境が必要)
```

---

## YAML 設定

詳細設定は多層 YAML で渡せます。適用順序:
```
デフォルト < config < CLI オプション < override-yaml
```

利用可能なすべてのオプションは [YAML リファレンス](yaml_reference.md) を参照してください。

---

## 出力ディレクトリ

`--out-dir` で結果の保存先を指定します:

```bash
--out-dir ./my_results/   # カスタム出力ディレクトリ
```

デフォルトの出力ディレクトリ:
- `all`: `./result_all/`
- `extract`: カレントディレクトリまたは `-o` で指定
- `mm-parm`: カレントディレクトリまたは `--out-prefix` で指定
- `define-layer`: カレントディレクトリまたは `-o` で指定
- `opt`: `./result_opt/`
- `tsopt`: `./result_tsopt/`
- `path-opt`: `./result_path_opt/`
- `path-search`: `./result_path_search/`
- `scan`: `./result_scan/`
- `scan2d`: `./result_scan2d/`
- `scan3d`: `./result_scan3d/`
- `freq`: `./result_freq/`
- `irc`: `./result_irc/`
- `dft`: `./result_dft/`

---

## 関連項目

- [はじめに](getting_started.md) -- インストールと初回実行
- [概念とワークフロー](concepts.md) -- ML/MM 3層システム、ONIOM 分解の全体像
- [典型エラー別レシピ](recipes_common_errors.md) -- 症状起点の切り分け
- [トラブルシューティング](troubleshooting.md) -- よくあるエラーと対処法
- [YAML リファレンス](yaml_reference.md) -- 全設定オプション
