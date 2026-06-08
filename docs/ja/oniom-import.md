# `oniom-import`

Gaussian/ORCA の ONIOM 入力ファイルを読み込み、XYZ と B-factor 層付き PDB を再構築し、外部で用意した ONIOM 入力を XYZ + 層付き PDB のペアとして ML/MM ツールチェーンに戻します。XYZ のコメント行に QM 領域の電荷と多重度（`q=<charge> m=<multiplicity>`）を記録し、層付き PDB は ML/Movable-MM/Frozen-MM 層を B-factor 列にエンコードします。インポートモードは `--mode` または入力拡張子から決定します。`--ref-pdb` を渡すと、層付き構造を再構築する際に参照 PDB から原子/残基メタデータを復元します。

## 実行例

```bash
mlmm oniom-import -i ts_guess.inp -o ts_guess_imported
```

拡張子からモード自動判定（`.gjf`/`.com` -> g16, `.inp` -> orca）:

```bash
mlmm oniom-import -i model.gjf -o model_imported
```

モードを明示:

```bash
mlmm oniom-import -i model.inp --mode orca -o model_imported
```

参照 PDB の命名/残基情報を保持（原子数一致が必要）:

```bash
mlmm oniom-import -i model.inp --ref-pdb complex_layered.pdb -o model_imported
```

## 処理の流れ
1. `--mode` または入力拡張子からインポートモードを決定します。
2. 座標と層情報を解析します。
 - Gaussian モード: ONIOM 座標行（`H`/`L` レイヤー）を使用。
 - ORCA モード: `%qmmm`（`QMAtoms`/`ActiveAtoms`）と `* xyz` ブロックを使用。
3. `<out_prefix>.xyz` を出力します。
4. `<out_prefix>_layered.pdb` を出力します。
 - `--ref-pdb` 未指定: 汎用的な原子/残基名で生成。
 - `--ref-pdb` 指定: 命名/残基メタデータを保持しつつ座標/B-factor を更新。

## 出力

```text
<out_prefix>.xyz
<out_prefix>_layered.pdb
```

- `<out_prefix>.xyz` が生成され、原子数が元 ONIOM 入力と一致する。
- `<out_prefix>_layered.pdb` が生成され、B-factor が `0/10/20` で層を表す。
- ログに mode、原子数、QM/Movable/Frozen の件数が表示される。

## CLI オプション

コマンド形式:

```bash
mlmm oniom-import -i INPUT.[gjf|com|inp] [--mode g16|orca] \
  [-o OUT_PREFIX] [--ref-pdb REF.pdb]
```

| オプション | 説明 | デフォルト |
| --- | --- | --- |
| `-i, --input FILE` | 入力 ONIOM ファイル（g16: `.gjf`/`.com`, ORCA: `.inp`）。 | 必須 |
| `--mode [g16|orca]` | 入力モード。未指定時は拡張子から推定。 | 自動判定 |
| `-o, --out-prefix PATH` | 出力プレフィックス。 | カレントディレクトリ上の入力 stem |
| `--ref-pdb FILE` | 原子名/残基メタデータ保持用の参照 PDB（原子数一致必須）。 | _None_ |

高度なオプションを含む全フラグ一覧は、生成される[コマンドリファレンス](../reference/commands/index.md)を参照してください。

## 関連項目

- [oniom_export](oniom-export.md) -- Gaussian ONIOM / ORCA QM/MM へのエクスポート
- [define_layer](define-layer.md) -- 層定義ルール
- [troubleshooting](troubleshooting.md) -- 詳細なトラブルシューティング
