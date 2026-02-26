# ONIOM インポート（`oniom-import`）

## 概要

> **要約:** Gaussian/ORCA の ONIOM 入力ファイルを読み込み、XYZ と B-factor 層付き PDB を再構築します。

### 概要
- **用途:** 外部 ONIOM/QM/MM 入力を mlmm_toolkit の後続処理へ戻したい場合。
- **手法:** Gaussian（`.gjf`/`.com`）または ORCA（`.inp`）入力を解析し、座標と層情報を復元して標準出力を生成。
- **出力:** `<out_prefix>.xyz` と `<out_prefix>_layered.pdb`。
- **層エンコード:** PDB の B-factor で ML(QM)=0.00、Movable MM=10.00、Frozen MM=20.00。

## 最小例

```bash
mlmm oniom-import -i ts_guess.inp -o ts_guess_imported
```

## 出力の見方
- `<out_prefix>.xyz` が生成され、原子数が元 ONIOM 入力と一致する。
- `<out_prefix>_layered.pdb` が生成され、B-factor が `0/10/20` で層を表す。
- ログに mode、原子数、QM/Movable/Frozen の件数が表示される。

## よくある例

```bash
# 拡張子からモード自動判定（.gjf/.com -> g16, .inp -> orca）
mlmm oniom-import -i model.gjf -o model_imported

# モードを明示
mlmm oniom-import -i model.inp --mode orca -o model_imported

# 参照 PDB の命名/残基情報を保持（原子数一致が必要）
mlmm oniom-import -i model.inp --ref-pdb complex_layered.pdb -o model_imported
```

## 使用法

```bash
mlmm oniom-import -i INPUT.[gjf|com|inp] [--mode g16|orca] \
  [-o OUT_PREFIX] [--ref-pdb REF.pdb]
```

## ワークフロー
1. `--mode` または入力拡張子からインポートモードを決定します。
2. 座標と層情報を解析します。
 - Gaussian モード: ONIOM 座標行（`H`/`L` レイヤー）を使用。
 - ORCA モード: `%qmmm`（`QMAtoms`/`ActiveAtoms`）と `* xyz` ブロックを使用。
3. `<out_prefix>.xyz` を出力します。
4. `<out_prefix>_layered.pdb` を出力します。
 - `--ref-pdb` 未指定: 汎用的な原子/残基名で生成。
 - `--ref-pdb` 指定: 命名/残基メタデータを保持しつつ座標/B-factor を更新。

## CLI オプション

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input FILE` | 入力 ONIOM ファイル（g16: `.gjf`/`.com`, ORCA: `.inp`）。 | 必須 |
| `--mode [g16|orca]` | 入力モード。未指定時は拡張子から推定。 | 自動判定 |
| `-o, --out-prefix PATH` | 出力プレフィックス。 | カレントディレクトリ上の入力 stem |
| `--ref-pdb FILE` | 原子名/残基メタデータ保持用の参照 PDB（原子数一致必須）。 | _None_ |

## 出力

```text
<out_prefix>.xyz
<out_prefix>_layered.pdb
```

## 注意事項
- 本コマンドは `mlmm oniom-export` 由来の ONIOM 入力形式を前提とします。
- 拡張子からモード推定できない場合は `--mode g16|orca` を指定してください。
- `--ref-pdb` 使用時は原子数が完全一致している必要があります。

---

## 関連項目

- [oniom_export](oniom_export.md) -- Gaussian ONIOM / ORCA QM/MM へのエクスポート
- [define_layer](define_layer.md) -- 層定義ルール
- [troubleshooting](troubleshooting.md) -- 詳細なトラブルシューティング
